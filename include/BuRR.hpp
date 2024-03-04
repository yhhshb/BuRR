#ifndef BURR_HPP
#define BURR_HPP

#include <cmath>
#include <chrono>
#include <tuple>
#include <string>
#include <vector>
#include <cassert>
#include "../bundled/biolib/include/external_memory_vector.hpp"
#include "../bundled/biolib/include/bit_operations.hpp"
#include "../bundled/biolib/include/size_iterator.hpp"

#include "constants.hpp"

namespace ribbon {

class option_t 
{
    public:
        double epsilon;
        std::size_t layers;
        uint64_t seed;
        std::string tmp_dir;
        std::size_t max_ram;
        std::size_t nthreads;
        std::size_t verbose;
        bool check;
};

#define CLASS_HEADER 
#define METHOD_HEADER BuRR // BuRR<Iterator, Hasher>
#define BUILD_FUNCTION_HEADER template <typename Iterator, typename Hasher>

class BuRR // Standard IR data structure, 2-bit variant
{
    public:
        using size_type = std::size_t; // Unsigned integer type (usually std::size_t)
        using difference_type = std::ptrdiff_t; // Signed integer type (usually std::ptrdiff_t)
        // using key_compare = Compare;
        // using allocator_type = Allocator;

        enum class ThreshMode : std::size_t { normal = 0, onebit, twobit };

        class const_iterator
        {
            public:
                const_iterator(const BuRR& parent, std::size_t idx);

            private:
                friend bool operator==(const const_iterator& itr1, const const_iterator& itr2) {return true;}
                friend bool operator!=(const const_iterator& itr1, const const_iterator& itr2) {return not(itr1 == itr2);}
        };

        BuRR(const option_t& build_options);

        BUILD_FUNCTION_HEADER
        void build(Iterator start, Iterator stop);

        BUILD_FUNCTION_HEADER
        void build(Iterator start, std::size_t n);

        std::size_t count() const noexcept; // returns the number of elements matching specific key

        template <typename K>
        const_iterator find(K const& key) const; // finds element with specific key

        template <typename K>
        bool contains(K const& key) const; // checks if the container contains element with specific key

    private:
        std::tuple<std::size_t, std::size_t> get_thresholds(std::size_t ribbon_width, double eps, std::size_t bucket_size);

        // void build_layer(emem::external_memory_vector<value_type>& pairs);
        option_t option_bundle;
        std::size_t bucket_size;
        std::size_t lower;
        std::size_t upper;
        std::size_t m;
};


CLASS_HEADER
BUILD_FUNCTION_HEADER
void 
METHOD_HEADER::build(Iterator start, Iterator stop)
{
    using key_type = typename Iterator::value_type::first_type;
    using mapped_type = typename Iterator::value_type::second_type;
    using value_type = std::pair<key_type, mapped_type>;
    // using reference = value_type&;
    // using const_reference = const value_type&;

    const std::size_t run_id = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    auto ribbon_width = ::bit::size<Hasher::hash_t>();
    auto a = std::ceil(std::log2(ribbon_width));
    assert(a >= 0);
    auto b = 4 * a; // 2bit mode
    auto c = std::floor(std::log2(((ribbon_width * ribbon_width) / b))); // w * w / (factor * log(w))
    assert(c >=0);
    bucket_size = std::size_t(1) << c; // round-down to next power of two
    [lower, upper] = get_thresholds(ribbon_width, option_bundle.epsilon, bucket_size);

    // std::string external_key_storage = join_path(opts.tmp_dir, random_name());
    emem::external_memory_vector<value_type> hashed_keys(option_bundle.max_ram, option_bundle.tmp_dir, util::get_name("hashes", run_id));
    std::transform(start, stop, std::back_inserter(hashed_keys), [this](auto &v) {return std::make_pair(Hasher::hash(v.first, option_bundle.seed), v.second);});
    m = (1 + option_bundle.epsilon) * hashed_keys.size(); // m = (1+e)*n
    m = ((m + ribbon_width - 1) / ribbon_width) * ribbon_width; // round up to next multiple of ribbon_width for interleaved storage

    // TODO
    // The final result is a (interleaved) matrix Z + bumping information and the last bumped elements
    // Z must depend on concrete types, independent of the template paramters
    // Z is generated from two local vectors using back substitution

    do {
        build_layer(hashed_keys);
    } while(false);
}

CLASS_HEADER
BUILD_FUNCTION_HEADER
void 
METHOD_HEADER::build(Iterator start, std::size_t n)
{
    build(iterators::size_iterator(start, 0), iterators::size_iterator(start, n));
}

// CLASS_HEADER
// void 
// METHOD_HEADER::build_layer(emem::external_memory_vector<value_type>& pairs)
// {
//     // banding add range (recursion)
    
//     if (option_bundle.verbose) {
//         // std::cerr 
//         // << "Bumped " << num_bumped << 
//         // " out of " << input_size << " (" << (num_bumped * 100.0 / input_size) << "%) with " << slots_per_item_ << " slots per item:\n"
//         // "\t" << empty_slots << "empty slots (" << empty_slots * 100.0 / storage_.GetNumSlots() << "%)\n";
//     }
// }

/*uncomment when everything works fine*/
// #undef CLASS_HEADER
// #undef METHOD_HEADER
// #undef BUILD_FUNCTION_HEADER

} // namespace ribbon

#endif // BURR_HPP


// template <
//     typename BandingStorage, 
//     typename Hasher, 
//     typename Iterator,
//     typename BumpStorage = std::vector<typename std::iterator_traits<Iterator>::value_type>
// >
// bool BandingAddRange(
//     BandingStorage *bs, 
//     Hasher &hasher, 
//     Iterator begin, 
//     Iterator end, 
//     BumpStorage *bump_vec) 
// {
//     using CoeffRow = typename BandingStorage::CoeffRow;
//     using Index = typename BandingStorage::Index;
//     using ResultRow = typename BandingStorage::ResultRow;
//     using Hash = typename Hasher::Hash;
//     constexpr bool kFCA1 = Hasher::kFirstCoeffAlwaysOne;
//     constexpr bool oneBitThresh = Hasher::kThreshMode == ThreshMode::onebit;

//     constexpr bool debug = false;
//     constexpr bool log = Hasher::log;

//     if (begin == end)
//         return true;

//     rocksdb::StopWatchNano timer(true);
//     const Index num_starts = bs->GetNumStarts();
//     const Index num_buckets = bs->GetNumBuckets();
//     sLOG << "Constructing ribbon with" << num_buckets
//          << "buckets,  num_starts = " << num_starts;

//     const auto num_items = end - begin; // std::distance(begin, end);
// #ifdef RIBBON_PASS_HASH
//     constexpr bool sparse = Hasher::kSparseCoeffs && Hasher::kCoeffBits < 128;
//     auto input = std::make_unique<
//         std::tuple<Index, Index, std::conditional_t<sparse, uint32_t, Hash>>[]>(
//         num_items);
// #else
//     auto input = std::make_unique<std::pair<Index, Index>[]>(num_items);
// #endif

//     {
//         sLOG << "Processing" << num_items << "items";

//         for (Index i = 0; i < static_cast<Index>(num_items); i++) {
//             const Hash h = hasher.GetHash(*(begin + i));
//             const Index start = hasher.GetStart(h, num_starts);
//             const Index sortpos = Hasher::StartToSort(start);
// #ifdef RIBBON_PASS_HASH
//             if constexpr (sparse) {
//                 uint32_t compact_hash = hasher.GetCompactHash(h);
//                 input[i] = std::make_tuple(sortpos, i, compact_hash);
//             } else {
//                 input[i] = std::make_tuple(sortpos, i, h);
//             }
// #else
//             input[i] = std::make_pair(sortpos, i);
// #endif
//         }
//     }
//     LOGC(log) << "\tInput transformation took "
//               << timer.ElapsedNanos(true) / 1e6 << "ms";
//     my_sort(input.get(), input.get() + num_items);
//     LOGC(log) << "\tSorting took " << timer.ElapsedNanos(true) / 1e6 << "ms";

//     const auto do_bump = [&](auto &vec) {
//         sLOG << "Bumping" << vec.size() << "items";
//         for (auto [row, idx] : vec) {
//             sLOG << "\tBumping row" << row << "item"
//                  << tlx::wrap_unprintable(*(begin + idx));
//             bs->SetCoeffs(row, 0);
//             bs->SetResult(row, 0);
//             bump_vec->push_back(*(begin + idx));
//         }
//         vec.clear();
//     };

//     Index last_bucket = 0;
//     bool all_good = true;
//     Index thresh = Hasher::NoBumpThresh();
//     // Bump cache (row, input item) pairs that may have to be bumped retroactively
//     Index last_cval = -1;
//     std::vector<std::pair<Index, Index>> bump_cache;
//     // For 1-bit thresholds, we also need an uncompressed bump cache for undoing
//     // all insertions with the same uncompressed value if we end up in the
//     // "plus" case with a separately stored threshold
//     [[maybe_unused]] Index last_val = -1;
//     [[maybe_unused]] std::conditional_t<oneBitThresh, decltype(bump_cache), int> unc_bump_cache;

// #ifndef RIBBON_PASS_HASH
//     auto next = *(begin + input[0].second);
// #endif

//     for (Index i = 0; i < static_cast<Index>(num_items); ++i) {
// #ifdef RIBBON_PASS_HASH
//         const auto [sortpos, idx, hash] = input[i];
// #else
//         const auto [sortpos, idx] = input[i];
// #endif
//         const Index start = Hasher::SortToStart(sortpos),
//                     bucket = Hasher::GetBucket(sortpos),
//                     val = Hasher::GetIntraBucket(sortpos),
//                     cval = hasher.Compress(val);
//         assert(bucket >= last_bucket);
//         assert(oneBitThresh || cval < Hasher::NoBumpThresh());

// #ifndef RIBBON_PASS_HASH
//         const Hash hash = hasher.GetHash(next);
//         if (i + 1 < num_items)
//             next = *(begin + input[i + 1].second);

//         // prefetch the cache miss far in advance, assuming the iterator
//         // is to contiguous storage
//         if (TLX_LIKELY(i + 32 < num_items))
//             __builtin_prefetch(&*begin + input[i + 32].second, 0, 1);
// #endif

//         if (bucket != last_bucket) {
//             // moving to next bucket
//             sLOG << "Moving to bucket" << bucket << "was" << last_bucket;
//             if constexpr (oneBitThresh) {
//                 unc_bump_cache.clear();
//                 last_val = val;
//             }
//             if (thresh == Hasher::NoBumpThresh()) {
//                 sLOG << "Bucket" << last_bucket << "has no bumped items";
//                 bs->SetMeta(last_bucket, thresh);
//             }
//             all_good = true;
//             last_bucket = bucket;
//             thresh = Hasher::NoBumpThresh(); // maximum == "no bumpage"
//             last_cval = cval;
//             bump_cache.clear();
//         } else if (!all_good) {
//             // direct hard bump
//             sLOG << "Directly bumping" << tlx::wrap_unprintable(*(begin + idx))
//                  << "from bucket" << bucket << "val" << val << cval << "start"
//                  << start << "sort" << sortpos << "hash" << std::hex << hash
//                  << "data"
//                  << (uint64_t)(Hasher::kIsFilter
//                                    ? hasher.GetResultRowFromHash(hash)
//                                    : hasher.GetResultRowFromInput(*(begin + idx)))
//                  << std::dec;
//             bump_vec->push_back(*(begin + idx));
//             continue;
//         } else if (cval != last_cval) {
//             // clear bump cache
//             sLOG << "Bucket" << bucket << "cval" << cval << "!=" << last_cval;
//             bump_cache.clear();
//             last_cval = cval;
//         }
//         if constexpr (oneBitThresh) {
//             // split into constexpr and normal if because unc_bump_cache isn't a
//             // vector if !oneBitThresh
//             if (val != last_val) {
//                 unc_bump_cache.clear();
//                 last_val = val;
//             }
//         }


//         const CoeffRow cr = hasher.GetCoeffs(hash);
//         const ResultRow rr = Hasher::kIsFilter
//                                  ? hasher.GetResultRowFromHash(hash)
//                                  : hasher.GetResultRowFromInput(*(begin + idx));

//         auto [success, row] = BandingAdd<kFCA1>(bs, start, cr, rr);
//         if (!success) {
//             assert(all_good);
//             if (bump_vec == nullptr)
//                 // bumping disabled, abort!
//                 return false;
//             // if we got this far, this is the first failure in this bucket,
//             // and we need to undo insertions with the same cval
//             sLOG << "First failure in bucket" << bucket << "val" << val
//                  << "start" << start << "sort" << sortpos << "hash" << std::hex
//                  << hash << "data" << (uint64_t)rr << std::dec << "for item"
//                  << tlx::wrap_unprintable(*(begin + idx)) << "-> threshold"
//                  << cval << "clash in row" << row;
//             thresh = cval;
//             if constexpr (oneBitThresh) {
//                 if (cval == 2) {
//                     sLOG << "First failure in bucket" << bucket << "val" << val
//                          << "is a 'plus' case (below threshold)";
//                     // "plus" case: store uncompressed threshold in hash table
//                     hasher.Set(bucket, val);
//                     // Set meta to 0 (some bumpage) but keep thresh at 2 so that
//                     // we don't retroactively bump everything when moving to the
//                     // next bucket
//                     bs->SetMeta(bucket, 0);
//                     all_good = false;

//                     // bump all items with the same uncompressed value
//                     do_bump(unc_bump_cache);
//                     sLOG << "Also bumping"
//                          << tlx::wrap_unprintable(*(begin + idx));
//                     bump_vec->push_back(*(begin + idx));
//                     // proceed to next item, don't do regular bumping (but only
//                     // if cval == 2, not generally!)
//                     continue;
//                 }
//             }
//             bs->SetMeta(bucket, thresh);
//             all_good = false;

//             do_bump(bump_cache);
//             bump_vec->push_back(*(begin + idx));
//         } else {
//             sLOG << "Insertion succeeded of item"
//                  << tlx::wrap_unprintable(*(begin + idx)) << "in pos" << row
//                  << "bucket" << bucket << "val" << val << cval << "start"
//                  << start << "sort" << sortpos << "hash" << std::hex << hash
//                  << "data" << (uint64_t)rr << std::dec;
//             bump_cache.emplace_back(row, idx);
//             if constexpr (oneBitThresh) {
//                 // also record in bump cache for uncompressed values
//                 unc_bump_cache.emplace_back(row, idx);
//             }
//         }
//     }
//     // set final threshold
//     if (thresh == Hasher::NoBumpThresh()) {
//         bs->SetMeta(last_bucket, thresh);
//     }

//     // migrate thresholds to hash table
//     if constexpr (oneBitThresh) {
//         hasher.Finalise(num_buckets);
//     }

//     LOGC(log) << "\tActual insertion took " << timer.ElapsedNanos(true) / 1e6
//               << "ms";
//     return true;
// }