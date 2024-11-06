namespace retrieval::ribbon {

#define CLASS_HEADER template <class SolutionStorage, class Hasher>
#define METHOD_HEADER BuRR<threshold_t::TWOBITS, SolutionStorage, Hasher>
#define BUILD_FUNCTION_HEADER template <class Iterator>

CLASS_HEADER
class BuRR<threshold_t::TWOBITS, SolutionStorage, Hasher> // Standard IR data structure, 2-bit variant
{
    public:
        using size_type = std::size_t; // Unsigned integer type (usually std::size_t)
        using difference_type = std::ptrdiff_t; // Signed integer type (usually std::ptrdiff_t)

        BuRR(const option_t& build_options);

        BUILD_FUNCTION_HEADER
        void build(Iterator start, Iterator stop);

        BUILD_FUNCTION_HEADER
        void build(Iterator start, std::size_t n);

        template <typename KeyType, typename ValueType>
        ValueType at (const KeyType& x) const noexcept;

    private:
        enum class bumping : uint8_t {
            ALL = 0, 
            MOST, 
            LITTLE, 
            NONE
        };
        class layer
        {
            public:
                layer(const METHOD_HEADER* parent);
                
                template <typename KeyType, typename ValueType>
                emem::external_memory_vector<std::pair<KeyType, ValueType>> 
                build(
                    const std::size_t run_id, 
                    const emem::external_memory_vector<std::pair<KeyType, ValueType>>& pairs
                );

                template <typename ValueType>
                std::optional<ValueType>
                at(const typename Hasher::hash_type hval) const noexcept;

            private:
                std::tuple<std::size_t, std::size_t, typename METHOD_HEADER::bumping, typename Hasher::hash_type> 
                hash_to_bucket(typename Hasher::hash_type hkey, std::size_t m) const noexcept;

                template <typename ValueType>
                std::optional<std::size_t>
                insert(
                    const std::size_t start, 
                    const typename Hasher::hash_type ribbon_hash, 
                    const ValueType value,
                    std::vector<typename Hasher::hash_type>& coefficients,
                    std::vector<ValueType>& results
                );

                bool
                bumped(bumping info) const noexcept;

                const METHOD_HEADER* parent;
                bit::packed::vector<std::size_t> bump_info;
                std::size_t m; // size of the layer
                SolutionStorage Z; // the actual data
        };

        std::tuple<std::size_t, std::size_t> 
        get_thresholds(std::size_t ribbon_width, double eps, std::size_t bucket_size);

        option_t option_bundle;
        std::size_t total_empty_slots;
        std::size_t bucket_size;
        std::size_t lower;
        std::size_t upper;
        std::vector<layer> layers;
        // slots_per_item = 1 + epsilon
};

CLASS_HEADER
METHOD_HEADER::BuRR(const option_t& build_options)
    : 
    option_bundle(build_options), 
    total_empty_slots(0),
    bucket_size(0),
    lower(0),
    upper(0)
{}

CLASS_HEADER
BUILD_FUNCTION_HEADER
void 
METHOD_HEADER::build(Iterator start, Iterator stop)
{
    // using key_type = typename Iterator::value_type::first_type;
    using mapped_type = typename Iterator::value_type::second_type;
    // using value_type = std::pair<key_type, mapped_type>;
    using hash_pair_type = std::pair<typename Hasher::hash_type, mapped_type>;

    const std::size_t run_id = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    auto ribbon_width = ::bit::size<Hasher::hash_type>();
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
    if (option_bundle.check and 
        std::adjacent_find(
            hashed_keys.cbegin(), 
            hashed_keys.cend(), 
            [](const auto &a, const auto &b) {
                return a.first == b.first;
            }
        ) != hashed_keys.cend()
    ) throw std::runtime_error("[Build] Hashes are not unique");
    
    do {
        layers.emplace_back();
        hashed_keys = layers.back.build(run_id, hashed_keys);
    } while(hashed_keys.size() != 0 and layers.size() < option_bundle.nlayers);
    assert(layers.size() <= option_bundle.nlayers);
    // if (layers.size() < option_bundle.nlayers) option_bundle.nlayers = layers.size(); // just use layers.size()

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
template <typename KeyType, typename ValueType>
ValueType 
METHOD_HEADER::at(const KeyType& x) const noexcept
{
    const auto ribbon_width = ::bit::size<Hasher::hash_type>();
    auto hash = Hasher::hash(x, option_bundle.seed);
    std::optional<ValueType> ret = std::nullopt;
    std::size_t i = 0;
    while(!ret and i < layers.size()) {
        ret = layers.at(i).at(hash);
    }
    assert(ret.has_value());
    return ret.value_or(static_cast<ValueType>(std::numeric_limits<ValueType>::max()));
}

CLASS_HEADER
std::tuple<std::size_t, std::size_t>
METHOD_HEADER::get_thresholds(std::size_t ribbon_width, double eps, std::size_t bucket_size)
{
    auto compute_thresholds = [&eps, &bucket_size](std::tuple<double, double, double, double> coefficients) {
        auto [la, lb, ua, ub] = coefficients;
        return std::make_tuple(
            static_cast<std::size_t>((la * eps + lb) * bucket_size), 
            static_cast<std::size_t>((ua * eps + ub) * bucket_size)
        );
    };
    std::tuple<double, double, double, double> coefficients;
    if (ribbon_width <= 16) coefficients = std::make_tuple(0.5, 0.6, 0.5, 0.82);
    else if (ribbon_width <= 32) coefficients = std::make_tuple(0.5, 0.7, 0.5, 0.87);
    else coefficients = std::make_tuple(1.30, 0.78, 0.75, 0.91);
    return compute_thresholds(coefficients);
};

CLASS_HEADER
METHOD_HEADER::layer::layer(const METHOD_HEADER* p) 
    : 
    parent(p),
    bump_info(bit::packed::vector<std::size_t>(2)) // 2-bit variant
{}

CLASS_HEADER
template <typename KeyType, typename ValueType> // Here KeyType are the hashed keys
emem::external_memory_vector<std::pair<KeyType, ValueType>> 
METHOD_HEADER::layer::build(
    const std::size_t run_id, 
    const emem::external_memory_vector<std::pair<KeyType, ValueType>>& pairs
)
{
    using hash_pair_type = std::pair<KeyType, ValueType>;
    bump_info = bit::packed::vector<std::size_t>(2); // 2-bit thresholds
    const auto ribbon_width = ::bit::size<KeyType>(); // recomputed here to avoid passing one additional argument
    m = (1 + parent->option_bundle.epsilon) * pairs.size(); // m = (1+e)*n (previously called num_slots)
    m = ((m + ribbon_width - 1) / ribbon_width) * ribbon_width; // round up to next multiple of ribbon_width for interleaved storage

    std::vector<typename Hasher::hash_type> coefficients(m);
    std::vector<ValueType> results(m);
    emem::external_memory_vector<hash_pair_type> bumped_items(parent->option_bundle.max_ram, parent->option_bundle.tmp_dir, util::get_name("bumped", run_id));
    const auto do_bump = [&coefficients, &results, &bumped_items](auto &vec) {
        for (auto [row, elem] : vec) {
            coefficients[row] = 0;
            results[row] = 0;
            bumped_items.push_back(elem);
        }
        vec.clear();
    };

    bool all_good = true;
    std::size_t prev_bucket = 0;
    auto thresh = bumping::NONE;
    auto last_cval = -1;
    std::vector<std::pair<std::size_t, hash_pair_type>> bump_cache;
    for(auto itr = pairs.cbegin(); itr != pairs.cend(); ++itr, ++itr) {
        auto [offset, bucket, cval, hash] = hash_to_bucket((*itr).first);
        if (bucket < prev_bucket) throw std::runtime_error("[build_layer] new bucket index < previous one");
        if (bucket != prev_bucket) { // moving to next bucket
            if (thresh == bumping::NONE) bump_info[prev_bucket] = thresh;
            all_good = true;
            prev_bucket = bucket;
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

        auto row = insert(offset, hash, (*itr).second, coefficients, results);
        if (row) {
            bump_cache.emplace_back(row, *itr);
        } else {
            assert(all_good);
            if (parent->option_bundle.nlayers == 1) return false; // bumping disabled, abort!
            // if we got this far, this is the first failure in this bucket,
            // and we need to undo insertions with the same cval
            thresh = cval;
            bump_info[bucket] = thresh;
            all_good = false;
            do_bump(bump_cache);
            bumped_items.push_back(*itr);
        }
    }
    if (thresh == bumping::NONE) bump_info[prev_bucket] = thresh; // set final threshold
    parent->total_empty_slots += m - pairs.size() + bumped_items.size();

    // The final result is a (interleaved) matrix Z + bumping information and the last bumped elements
    // Z must depend on concrete types, independent of the template paramters
    // Z is generated from two local vectors using back substitution
    Z.backsubstitution(coefficients, results);
    return bumped_items;
}

CLASS_HEADER
template <typename ValueType>
std::optional<std::size_t>
METHOD_HEADER::layer::insert(
    const std::size_t start,
    const typename Hasher::hash_type ribbon_hash,
    const ValueType value,
    std::vector<typename Hasher::hash_type>& coefficients,
    std::vector<ValueType>& results
) 
{
    std::size_t pos = start;
    // if constexpr (!kFCA1) {
    //     int tz = bit::lsbll(ribbon_hash);
    //     pos += static_cast<size_t>(tz);
    //     ribbon_hash >>= tz;
    // }

    while (true) {
        assert(pos < m);
        assert((ribbon_hash & 1) == 1);

        auto other = coefficients.at(pos);
        if (other == 0) { // found an empty slot, insert
            coefficients[pos] = ribbon_hash;
            results[pos] = value;
            return pos;
        }

        assert((other & 1) == 1);
        ribbon_hash ^= other;
        value ^= results.at(pos);
        if (ribbon_hash == 0) return std::nullopt; // linearly dependent!

        // move to now-leading coefficient
        auto tz = bit::lsbll(ribbon_hash);
        pos += tz;
        ribbon_hash >>= tz;
    }
    throw std::runtime_error("[layer::insert] This should never happen");
}

CLASS_HEADER
std::tuple<
    std::size_t, // offset
    std::size_t, // position inside the bucket
    typename METHOD_HEADER::bumping, // bump threshold
    typename Hasher::hash_type // (re)hashed value (the key)
> 
METHOD_HEADER::layer::hash_to_bucket(typename Hasher::hash_type hkey, std::size_t m) const noexcept
{
    const auto bucket_size = parent->bucket_size;
    const auto hash = Hasher::hash(hkey, 0); // rehash for current run. IMPROVEMENT: maybe use a fast XOR rehasher
    assert(hash); // hash must contain at least one set bit for ribbon
    const std::size_t start = toolbox::fastrange(hash, m - ::bit::size<Hasher::hash_type>() + 1);
    assert(bucket_size != 0); // just check if bucket size is initialized
    const std::size_t sortpos = start ^ (bucket_size - 1); // ???

    const std::size_t val = sortpos % bucket_size;
    bumping thr_type;
    if (val >= bucket_size) thr_type = bumping::NONE; // none bumped
    else if (val > parent->upper) thr_type = bumping::SOME; // some bumped
    else if (val > parent->lower) thr_type = bumping::MOST; // most bumped
    else thr_type = bumping::ALL; // all bumped
    return std::make_tuple(start, sortpos / bucket_size, thr_type, hash);
}

CLASS_HEADER
template <typename ValueType>
std::optional<ValueType>
METHOD_HEADER::layer::at(const typename Hasher::hash_type hval) const noexcept
{
    // iss.PrefetchQuery(segment);
    auto [offset, bucket, cval, hash] = hash_to_bucket(hval, m);
    if (cval >= bump_info.at(bucket)) return std::nullopt; // this works since bumping threshold are assigned increasing integer codes (the enum)
    return Z.template at<typename Hasher::hash_type, ValueType>(offset, hash);
}

#undef CLASS_HEADER
#undef METHOD_HEADER
#undef BUILD_FUNCTION_HEADER

} // namespace retrieval::ribbon