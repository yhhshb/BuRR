#ifndef BURR_HPP
#define BURR_HPP

#include <cstdint>
#include <algorithm>
#include <iterator>
#include <vector>
#include <string>
#include "constants.hpp"

namespace ribbon {

class base
{

};

#define FTMPL_CLASS_HEADER template <typename Hasher, typename Iterator, typename RibbonType>
#define FTMPL_CLASS_HEADER_HEADER template <typename Hasher, typename Iterator, typename RibbonType = uint64_t>
#define FTMPL_METHOD_HEADER filter<Hasher, Iterator, RibbonType>

FTMPL_CLASS_HEADER
class filter : base
{
    public:
        using key_type = typename Iterator::value_type::first_type;
        using mapped_type = typename Iterator::value_type::second_type;
        using value_type = std::pair<key_type, mapped_type>;
        using size_type = std::size_t; // Unsigned integer type (usually std::size_t)
        using difference_type = std::ptrdiff_t; // Signed integer type (usually std::ptrdiff_t)
        // using key_compare = Compare;
        // using allocator_type = Allocator;
        using reference = value_type&;
        using const_reference = const value_type&;

        class option_t
        {
            public:
                option_t() 
                    : seed(42), 
                      load_factor(0.95),
                      normal_threshold(true),
                      check(false),
                      verbose(false)
                {}
                std::size_t seed;
                double load_factor; // remember: bumping deals with collisions so load factor can be < 1
                std::size_t max_ram,
                std::string tmp_dir,
                bool normal_threshold,
                bool check;
                bool verbose;
        }

        filter(Iterator start, Iterator stop, option_t const& options);

        std::size_t count() const noexcept; // returns the number of elements matching specific key

        template <typename K>
        const_iterator find(K const& key) const; // finds element with specific key

        template <typename K>
        bool contains(K const& key) const; // checks if the container contains element with specific key

    private:
        option_t opts;
        std::size_t bucket_size;
        void build(Iterator start, Iterator stop, option_t const& options);
        void build_layer(external_memory_vector<value_type>& pairs);

};

FTMPL_CLASS_HEADER 
FTMPL_METHOD_HEADER::filter(Iterator start, Iterator stop, option_t const& options)
    : opts(options)
{
    auto a = tlx::integer_log2_ceil(coeff_bits);
    auto b = ((mode == opts.normal_threshold ? 2 : 4) * a);
    auto c = tlx::integer_log2_floor(coeff_bits * coeff_bits / b); // w * w / (factor * log(w))
    bucket_size = std::size_t(1) << c; // round down to next power of two
    build(start, stop);
}

FTMPL_CLASS_HEADER 
FTMPL_METHOD_HEADER::build(Iterator start, Iterator stop)
    : opts(options)
{
    std::string external_key_storage = join_path(opts.tmp_dir, random_name());
    external_memory_vector<value_type> hashed_keys(opts.max_ram, external_key_storage);
    std::transform(start, stop, std::back_inserter(hashed_keys), [](auto &v) {return std::make_pair(Hasher::hash(v.first, opts.seed), v.second);})

    do {
        build_layer(hashed_keys);
    } while();
}

FTMPL_CLASS_HEADER
void 
FTMPL_METHOD_HEADER::build_layer(external_memory_vector<value_type>& pairs)
{
    // banding add range (recursion)
    
    if (opts.verbose) {
        std::cerr << "Bumped " << num_bumped << 
        " out of " << input_size << 
        " (" << (num_bumped * 100.0 / input_size) << "%) with " << slots_per_item_ << " slots per item:\n"
        "\t" << empty_slots << "empty slots (" << empty_slots * 100.0 / storage_.GetNumSlots() << "%)\n";
    }
}

} // namespace ribbon

#endif // BURR_HPP

template <typename BandingStorage, typename Hasher, typename Iterator,
          typename BumpStorage = std::vector<typename std::iterator_traits<Iterator>::value_type>>
bool BandingAddRange(BandingStorage *bs, Hasher &hasher, Iterator begin,
                     Iterator end, BumpStorage *bump_vec) {
    using CoeffRow = typename BandingStorage::CoeffRow;
    using Index = typename BandingStorage::Index;
    using ResultRow = typename BandingStorage::ResultRow;
    using Hash = typename Hasher::Hash;
    constexpr bool kFCA1 = Hasher::kFirstCoeffAlwaysOne;
    constexpr bool oneBitThresh = Hasher::kThreshMode == ThreshMode::onebit;

    constexpr bool debug = false;
    constexpr bool log = Hasher::log;

    if (begin == end)
        return true;

    rocksdb::StopWatchNano timer(true);
    const Index num_starts = bs->GetNumStarts();
    const Index num_buckets = bs->GetNumBuckets();
    sLOG << "Constructing ribbon with" << num_buckets
         << "buckets,  num_starts = " << num_starts;

    const auto num_items = end - begin; // std::distance(begin, end);
#ifdef RIBBON_PASS_HASH
    constexpr bool sparse = Hasher::kSparseCoeffs && Hasher::kCoeffBits < 128;
    auto input = std::make_unique<
        std::tuple<Index, Index, std::conditional_t<sparse, uint32_t, Hash>>[]>(
        num_items);
#else
    auto input = std::make_unique<std::pair<Index, Index>[]>(num_items);
#endif

    {
        sLOG << "Processing" << num_items << "items";

        for (Index i = 0; i < static_cast<Index>(num_items); i++) {
            const Hash h = hasher.GetHash(*(begin + i));
            const Index start = hasher.GetStart(h, num_starts);
            const Index sortpos = Hasher::StartToSort(start);
#ifdef RIBBON_PASS_HASH
            if constexpr (sparse) {
                uint32_t compact_hash = hasher.GetCompactHash(h);
                input[i] = std::make_tuple(sortpos, i, compact_hash);
            } else {
                input[i] = std::make_tuple(sortpos, i, h);
            }
#else
            input[i] = std::make_pair(sortpos, i);
#endif
        }
    }
    LOGC(log) << "\tInput transformation took "
              << timer.ElapsedNanos(true) / 1e6 << "ms";
    my_sort(input.get(), input.get() + num_items);
    LOGC(log) << "\tSorting took " << timer.ElapsedNanos(true) / 1e6 << "ms";

    const auto do_bump = [&](auto &vec) {
        sLOG << "Bumping" << vec.size() << "items";
        for (auto [row, idx] : vec) {
            sLOG << "\tBumping row" << row << "item"
                 << tlx::wrap_unprintable(*(begin + idx));
            bs->SetCoeffs(row, 0);
            bs->SetResult(row, 0);
            bump_vec->push_back(*(begin + idx));
        }
        vec.clear();
    };

    Index last_bucket = 0;
    bool all_good = true;
    Index thresh = Hasher::NoBumpThresh();
    // Bump cache (row, input item) pairs that may have to be bumped retroactively
    Index last_cval = -1;
    std::vector<std::pair<Index, Index>> bump_cache;
    // For 1-bit thresholds, we also need an uncompressed bump cache for undoing
    // all insertions with the same uncompressed value if we end up in the
    // "plus" case with a separately stored threshold
    [[maybe_unused]] Index last_val = -1;
    [[maybe_unused]] std::conditional_t<oneBitThresh, decltype(bump_cache), int> unc_bump_cache;

#ifndef RIBBON_PASS_HASH
    auto next = *(begin + input[0].second);
#endif

    for (Index i = 0; i < static_cast<Index>(num_items); ++i) {
#ifdef RIBBON_PASS_HASH
        const auto [sortpos, idx, hash] = input[i];
#else
        const auto [sortpos, idx] = input[i];
#endif
        const Index start = Hasher::SortToStart(sortpos),
                    bucket = Hasher::GetBucket(sortpos),
                    val = Hasher::GetIntraBucket(sortpos),
                    cval = hasher.Compress(val);
        assert(bucket >= last_bucket);
        assert(oneBitThresh || cval < Hasher::NoBumpThresh());

#ifndef RIBBON_PASS_HASH
        const Hash hash = hasher.GetHash(next);
        if (i + 1 < num_items)
            next = *(begin + input[i + 1].second);

        // prefetch the cache miss far in advance, assuming the iterator
        // is to contiguous storage
        if (TLX_LIKELY(i + 32 < num_items))
            __builtin_prefetch(&*begin + input[i + 32].second, 0, 1);
#endif

        if (bucket != last_bucket) {
            // moving to next bucket
            sLOG << "Moving to bucket" << bucket << "was" << last_bucket;
            if constexpr (oneBitThresh) {
                unc_bump_cache.clear();
                last_val = val;
            }
            if (thresh == Hasher::NoBumpThresh()) {
                sLOG << "Bucket" << last_bucket << "has no bumped items";
                bs->SetMeta(last_bucket, thresh);
            }
            all_good = true;
            last_bucket = bucket;
            thresh = Hasher::NoBumpThresh(); // maximum == "no bumpage"
            last_cval = cval;
            bump_cache.clear();
        } else if (!all_good) {
            // direct hard bump
            sLOG << "Directly bumping" << tlx::wrap_unprintable(*(begin + idx))
                 << "from bucket" << bucket << "val" << val << cval << "start"
                 << start << "sort" << sortpos << "hash" << std::hex << hash
                 << "data"
                 << (uint64_t)(Hasher::kIsFilter
                                   ? hasher.GetResultRowFromHash(hash)
                                   : hasher.GetResultRowFromInput(*(begin + idx)))
                 << std::dec;
            bump_vec->push_back(*(begin + idx));
            continue;
        } else if (cval != last_cval) {
            // clear bump cache
            sLOG << "Bucket" << bucket << "cval" << cval << "!=" << last_cval;
            bump_cache.clear();
            last_cval = cval;
        }
        if constexpr (oneBitThresh) {
            // split into constexpr and normal if because unc_bump_cache isn't a
            // vector if !oneBitThresh
            if (val != last_val) {
                unc_bump_cache.clear();
                last_val = val;
            }
        }


        const CoeffRow cr = hasher.GetCoeffs(hash);
        const ResultRow rr = Hasher::kIsFilter
                                 ? hasher.GetResultRowFromHash(hash)
                                 : hasher.GetResultRowFromInput(*(begin + idx));

        auto [success, row] = BandingAdd<kFCA1>(bs, start, cr, rr);
        if (!success) {
            assert(all_good);
            if (bump_vec == nullptr)
                // bumping disabled, abort!
                return false;
            // if we got this far, this is the first failure in this bucket,
            // and we need to undo insertions with the same cval
            sLOG << "First failure in bucket" << bucket << "val" << val
                 << "start" << start << "sort" << sortpos << "hash" << std::hex
                 << hash << "data" << (uint64_t)rr << std::dec << "for item"
                 << tlx::wrap_unprintable(*(begin + idx)) << "-> threshold"
                 << cval << "clash in row" << row;
            thresh = cval;
            if constexpr (oneBitThresh) {
                if (cval == 2) {
                    sLOG << "First failure in bucket" << bucket << "val" << val
                         << "is a 'plus' case (below threshold)";
                    // "plus" case: store uncompressed threshold in hash table
                    hasher.Set(bucket, val);
                    // Set meta to 0 (some bumpage) but keep thresh at 2 so that
                    // we don't retroactively bump everything when moving to the
                    // next bucket
                    bs->SetMeta(bucket, 0);
                    all_good = false;

                    // bump all items with the same uncompressed value
                    do_bump(unc_bump_cache);
                    sLOG << "Also bumping"
                         << tlx::wrap_unprintable(*(begin + idx));
                    bump_vec->push_back(*(begin + idx));
                    // proceed to next item, don't do regular bumping (but only
                    // if cval == 2, not generally!)
                    continue;
                }
            }
            bs->SetMeta(bucket, thresh);
            all_good = false;

            do_bump(bump_cache);
            bump_vec->push_back(*(begin + idx));
        } else {
            sLOG << "Insertion succeeded of item"
                 << tlx::wrap_unprintable(*(begin + idx)) << "in pos" << row
                 << "bucket" << bucket << "val" << val << cval << "start"
                 << start << "sort" << sortpos << "hash" << std::hex << hash
                 << "data" << (uint64_t)rr << std::dec;
            bump_cache.emplace_back(row, idx);
            if constexpr (oneBitThresh) {
                // also record in bump cache for uncompressed values
                unc_bump_cache.emplace_back(row, idx);
            }
        }
    }
    // set final threshold
    if (thresh == Hasher::NoBumpThresh()) {
        bs->SetMeta(last_bucket, thresh);
    }

    // migrate thresholds to hash table
    if constexpr (oneBitThresh) {
        hasher.Finalise(num_buckets);
    }

    LOGC(log) << "\tActual insertion took " << timer.ElapsedNanos(true) / 1e6
              << "ms";
    return true;
}