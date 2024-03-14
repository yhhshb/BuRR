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

CLASS_HEADER
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
        std::tuple<std::size_t, std::size_t> 
        get_thresholds(std::size_t ribbon_width, double eps, std::size_t bucket_size);

        template <typename KeyType, typename ValueType, typename Hasher>
        emem::external_memory_vector<std::pair<KeyType, ValueType>> 
        build_layer(
            const std::size_t run_id, 
            const emem::external_memory_vector<std::pair<KeyType, ValueType>>& pairs, 
            std::vector<KeyType>& coefficients, 
            std::vector<ValueType>& results
        );

        option_t option_bundle;
        std::size_t bucket_size;
        std::size_t lower;
        std::size_t upper;
        std::size_t m;
        bit::packed::vector<std::size_t> bump_info;
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
template <typename KeyType, typename ValueType, typename Hasher>
emem::external_memory_vector<std::pair<KeyType, ValueType>> 
METHOD_HEADER::build_layer(
    const std::size_t run_id, 
    const emem::external_memory_vector<std::pair<KeyType, ValueType>>& pairs, 
    std::vector<KeyType>& coefficients, 
    std::vector<ValueType>& results)
{
    using hash_pair_type = std::pair<KeyType, ValueType>;
    auto ribbon_width = ::bit::size<KeyType>(); // recomputed here to avoid passing one additional argument
    m = (1 + option_bundle.epsilon) * pairs.size(); // m = (1+e)*n
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
        const auto hash = Hasher::hash((*itr).first, 0); // rehash for current run. IMPROVEMENT: maybe use a fast XOR rehasher
        assert(hash); // hash must contain at least one set bit for ribbon
        const std::size_t start = hasher.GetStart(hash, num_starts);
        const std::size_t sortpos = Hasher::SortToStart(start);
        const std::size_t bucket = Hasher::GetBucket(sortpos);
        const std::size_t val = Hasher::GetIntraBucket(sortpos);
        const std::size_t cval = hasher.Compress(val);

        if (bucket < prev_bucket) throw std::runtime_error("[build_layer] new bucket index < previous one");

        if (bucket != last_bucket) { // moving to next bucket
            if (thresh == Hasher::NoBumpThresh()) {
                // sLOG << "Bucket" << last_bucket << "has no bumped items";
                bump_info[last_bucket] = thresh;
            }
            all_good = true;
            last_bucket = bucket;
            thresh = Hasher::NoBumpThresh(); // maximum == "no bumpage"
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
    if (thresh == Hasher::NoBumpThresh()) {
        bump_info[last_bucket] = thresh;
    }

    total_empty_slots += m - pairs.size() + bumped_items.size();

    // TODO
    // The final result is a (interleaved) matrix Z + bumping information and the last bumped elements
    // Z must depend on concrete types, independent of the template paramters
    // Z is generated from two local vectors using back substitution

    // TODO Z = backsubstitution(coefficients, results);
}

} // namespace ribbon

#endif // BURR_HPP