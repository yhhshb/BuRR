#ifndef BURR_HPP
#define BURR_HPP

#include <cmath>
#include <chrono>
#include <tuple>
#include <string>
#include <vector>
#include <cassert>
#include "../bundled/biolib/include/external_memory_vector.hpp"
#include "../bundled/biolib/include/packed_vector.hpp"
#include "../bundled/biolib/include/bit_operations.hpp"
#include "../bundled/biolib/include/size_iterator.hpp"
#include "../bundled/biolib/include/toolbox.hpp"

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

#define CLASS_HEADER template <class SolutionStorage, class Hasher>
#define METHOD_HEADER BuRR<SolutionStorage, Hasher>
#define BUILD_FUNCTION_HEADER template <class Iterator>

CLASS_HEADER
class BuRR // Standard IR data structure, 2-bit variant
{
    public:
        using size_type = std::size_t; // Unsigned integer type (usually std::size_t)
        using difference_type = std::ptrdiff_t; // Signed integer type (usually std::ptrdiff_t)
        // using key_compare = Compare;
        // using allocator_type = Allocator;

        enum class bumping : uint8_t {ALL = 0, MOST, LITTLE, NONE};

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

        template <typename KeyType, typename ValueType>
        ValueType at (const KeyType& x) const noexcept;

        template <typename KeyType>
        std::size_t query(const KeyType key);

    private:
        std::tuple<std::size_t, std::size_t> 
        get_thresholds(std::size_t ribbon_width, double eps, std::size_t bucket_size);

        template <typename KeyType, typename ValueType>
        emem::external_memory_vector<std::pair<KeyType, ValueType>> 
        build_layer(
            const std::size_t run_id, 
            const emem::external_memory_vector<std::pair<KeyType, ValueType>>& pairs, 
            std::vector<KeyType>& coefficients, 
            std::vector<ValueType>& results
        );

        std::pair<std::size_t, bumping> 
        hash_to_bucket(typename Hasher::hash_t hkey, std::size_t m) const noexcept;

        option_t option_bundle;
        std::size_t bucket_size;
        std::size_t lower;
        std::size_t upper;
        // std::size_t m;
        bit::packed::vector<std::size_t> bump_info;
        SolutionStorage Z;
        // slots_per_item = 1 + epsilon
};


CLASS_HEADER
BUILD_FUNCTION_HEADER
void 
METHOD_HEADER::build(Iterator start, Iterator stop)
{
    using key_type = typename Iterator::value_type::first_type;
    using mapped_type = typename Iterator::value_type::second_type;
    // using value_type = std::pair<key_type, mapped_type>;
    using hash_pair_type = std::pair<typename Hasher::hash_t, mapped_type>;

    bump_info = bit::packed::vector<std::size_t>(2);
    const std::size_t run_id = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    auto ribbon_width = ::bit::size<Hasher::hash_t>();
    auto a = std::ceil(std::log2(ribbon_width));
    assert(a >= 0);
    auto b = 4 * a; // 2-bit variant
    auto c = std::floor(std::log2(((ribbon_width * ribbon_width) / b))); // w * w / (factor * log(w))
    assert(c >=0);
    bucket_size = std::size_t(1) << c; // round-down to next power of two
    {
    auto [low, up] = get_thresholds(ribbon_width, option_bundle.epsilon, bucket_size);
    lower = low;
    upper = up;
    }

    emem::external_memory_vector<hash_pair_type> hashed_keys(option_bundle.max_ram, option_bundle.tmp_dir, util::get_name("hashes", run_id));
    std::transform(start, stop, std::back_inserter(hashed_keys), [this](auto &v) {return std::make_pair(Hasher::hash(v.first, option_bundle.seed), v.second);});
    
    std::size_t layer = 0;
    do {
        hashed_keys = build_layer(run_id, hashed_keys, coefficients, results);
        ++layer;
    } while(hashed_keys.size() != 0 and layer < option_bundle.layers);
    assert(layer <= option_bundle.layers);
    if (layer < option_bundle.layers) option_bundle.layers = layer;

    if (hashed_keys.size() != 0) {
        // TODO save last items into fallback data structure
    }
}

CLASS_HEADER
BUILD_FUNCTION_HEADER
void 
METHOD_HEADER::build(Iterator start, std::size_t n)
{
    build(iterators::size_iterator(start, 0), iterators::size_iterator(start, n));
}



CLASS_HEADER
template <typename KeyType, typename ValueType> // Here KeyType are the hashed keys
emem::external_memory_vector<std::pair<KeyType, ValueType>> 
METHOD_HEADER::build_layer(
    const std::size_t run_id, 
    const emem::external_memory_vector<std::pair<KeyType, ValueType>>& pairs, 
    std::vector<KeyType>& coefficients, 
    std::vector<ValueType>& results)
{
    using hash_pair_type = std::pair<KeyType, ValueType>;
    auto ribbon_width = ::bit::size<KeyType>(); // recomputed here to avoid passing one additional argument
    std::size_t m = (1 + option_bundle.epsilon) * pairs.size(); // m = (1+e)*n (previously called num_slots)
    m = ((m + ribbon_width - 1) / ribbon_width) * ribbon_width; // round up to next multiple of ribbon_width for interleaved storage

    std::vector<typename Hasher::hash_t> coefficients(m);
    std::vector<mapped_type> results(m);
    emem::external_memory_vector<hash_pair_type> bumped_items(option_bundle.max_ram, option_bundle.tmp_dir, util::get_name("bumped", run_id));

    //----------------------------------------------------------------------------------------------------------------

    if (option_bundle.check and 
        std::adjacent_find(
            pairs.cbegin(), 
            pairs.cend(), 
            [](const auto &a, const auto &b) {
                return a.first == b.first;
            }
        ) != end
    ) throw std::runtime_error("[Build] Hashes are not unique");

    std::size_t prev_bucket = 0;
    std::vector<std::pair<std::size_t, hash_pair_type>> bump_cache;
    for(auto itr = pairs.cbegin(); itr != pairs.cend(); ++itr, ++i) {
        auto [bucket, cval] = hash_to_cval((*itr).first);

        if (bucket < prev_bucket) throw std::runtime_error("[build_layer] new bucket index < previous one");

        if (bucket != last_bucket) { // moving to next bucket
            if (thresh == bumping::NONE) {
                // sLOG << "Bucket" << last_bucket << "has no bumped items";
                bump_info[last_bucket] = thresh;
            }
            all_good = true;
            last_bucket = bucket;
            thresh = bumping::NONE; // maximum == "no bumpage"
            last_cval = cval;
            bump_cache.clear();
        } else if (!all_good) { // direct hard bump
            bumped_items.push_back(*itr);
            continue;
        } else if (cval != last_cval) { // clear bump cache
            bump_cache.clear();
            last_cval = cval;
        }

        const auto do_bump = [&coefficients, &results, &bumped_items](auto &vec) {
            for (auto [row, elem] : vec) {
                coefficients[row] = 0;
                results[row] = 0;
                bumped_items.push_back(elem);
            }
            vec.clear();
        };

        auto [success, row] = BandingAdd<kFCA1>(bs, start, hash, (*itr).second);
        if (!success) {
            assert(all_good);
            if (option_bundle.layers == 1) return false; // bumping disabled, abort!
            // if we got this far, this is the first failure in this bucket,
            // and we need to undo insertions with the same cval
            thresh = cval;
            bump_info[bucket] = thresh;
            all_good = false;
            do_bump(bump_cache);
            bumped_items.push_back(*itr);
        } else {
            bump_cache.emplace_back(row, *itr);
        }
    }

    // set final threshold
    if (thresh == bumping::NONE) {
        bump_info[last_bucket] = thresh;
    }

    total_empty_slots += m - pairs.size() + bumped_items.size();

    // TODO
    // The final result is a (interleaved) matrix Z + bumping information and the last bumped elements
    // Z must depend on concrete types, independent of the template paramters
    // Z is generated from two local vectors using back substitution

    Z.backsubstitution(coefficients, results);
}

CLASS_HEADER
std::pair<std::size_t, typename METHOD_HEADER::bumping> 
METHOD_HEADER::hash_to_bucket(typename Hasher::hash_t hkey, std::size_t m) const noexcept
{
    std::pair<bumping_t, std::size_t> toret;
    const auto hash = Hasher::hash(hkey, 0); // rehash for current run. IMPROVEMENT: maybe use a fast XOR rehasher
    assert(hash); // hash must contain at least one set bit for ribbon
    const std::size_t start = toolbox::fastrange(hash, m - ::bit::size<Hasher::hash_t>() + 1);
    assert(bucket_size != 0); // just check if bucket size is initialized
    const std::size_t sortpos = start ^ (bucket_size - 1);
    const std::size_t toret.first = sortpos / bucket_size;
    const std::size_t val = sortpos % bucket_size;
    
    if (val >= bucket_size) toret.second = bumping::NONE; // none bumped
    else if (val > upper) toret.second = bumping::SOME; // some bumped
    else if (val > lower) toret.second = bumping::MOST; // most bumped
    else toret.second = bumping::ALL; // all bumped
    return toret;
}

CLASS_HEADER
template <typename KeyType, typename ValueType>
ValueType 
METHOD_HEADER::at (const KeyType& x) const noexcept
{
    auto ribbon_width = ::bit::size<Hasher::hash_t>();
    

}

// General retrieval query a key from InterleavedSolutionStorage.
template <typename InterleavedSolutionStorage, typename Hasher>
std::pair<bool, typename InterleavedSolutionStorage::ResultRow>
InterleavedRetrievalQuery(const typename HashTraits<Hasher>::mhc_or_key_t &key,
                          const Hasher &hasher,
                          const InterleavedSolutionStorage &iss) {
    using Hash = typename Hasher::Hash;
    using Index = typename InterleavedSolutionStorage::Index;
    using CoeffRow = typename InterleavedSolutionStorage::CoeffRow;
    using ResultRow = typename InterleavedSolutionStorage::ResultRow;

    static_assert(sizeof(Index) == sizeof(typename Hasher::Index),
                  "must be same");
    static_assert(sizeof(CoeffRow) == sizeof(typename Hasher::CoeffRow),
                  "must be same");

    constexpr bool debug = false;
    constexpr auto kCoeffBits = static_cast<Index>(sizeof(CoeffRow) * 8U);
    constexpr Index num_columns = InterleavedSolutionStorage::kResultBits;

    // don't query an empty ribbon, please
    assert(iss.GetNumSlots() >= kCoeffBits);
    const Hash hash = hasher.GetHash(key);
    const Index start_slot = hasher.GetStart(hash, iss.GetNumStarts());
    const Index bucket = hasher.GetBucket(start_slot);

    const Index start_block_num = start_slot / kCoeffBits;
    Index segment = start_block_num * num_columns;
    iss.PrefetchQuery(segment);

    const Index val = hasher.GetIntraBucketFromStart(start_slot),
                cval = hasher.Compress(val);

    if (CheckBumped(val, cval, bucket, hasher, iss)) {
        sLOG << "Item was bumped, hash" << hash << "start" << start_slot
             << "bucket" << bucket << "val" << val << cval << "thresh"
             << (size_t)iss.GetMeta(bucket);
        return std::make_pair(true, 0);
    }

    sLOG << "Searching in bucket" << bucket << "start" << start_slot << "val"
         << val << cval << "below thresh =" << (size_t)iss.GetMeta(bucket);

    const Index start_bit = start_slot % kCoeffBits;
    const CoeffRow cr = hasher.GetCoeffs(hash);

    ResultRow sr = 0;
    const CoeffRow cr_left = cr << start_bit;
    for (Index i = 0; i < num_columns; ++i) {
        sr ^= rocksdb::BitParity(iss.GetSegment(segment + i) & cr_left) << i;
    }

    if (start_bit > 0) {
        segment += num_columns;
        const CoeffRow cr_right = cr >> (kCoeffBits - start_bit);
        for (Index i = 0; i < num_columns; ++i) {
            sr ^= rocksdb::BitParity(iss.GetSegment(segment + i) & cr_right) << i;
        }
    }

    return std::make_pair(false, sr);
}

} // namespace ribbon

#endif // BURR_HPP
