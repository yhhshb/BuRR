namespace retrieval::ribbon {

#define CLASS_HEADER template <class SolutionStorage, class Hasher>
#define METHOD_HEADER BuRR<threshold_t::ONEBIT, SolutionStorage, Hasher>
#define BUILD_FUNCTION_HEADER template <class Iterator>

CLASS_HEADER
class BuRR<threshold_t::ONEBIT, SolutionStorage, Hasher> // Standard IR data structure, 2-bit variant
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
        class layer
        {
            public:
                layer(const METHOD_HEADER* parent);
                
            private:
        };

        std::size_t get_threshold(std::size_t ribbon_width, double eps, std::size_t bucket_size);

        option_t option_bundle;
        std::size_t total_empty_slots;
        std::size_t bucket_size;
        std::size_t threshold;
        std::vector<layer> layers;
        // slots_per_item = 1 + epsilon
};

CLASS_HEADER
METHOD_HEADER::BuRR(const option_t& build_options)
    : 
    option_bundle(build_options), 
    total_empty_slots(0),
    bucket_size(0),
    threshold(0)
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
    threshold = get_threshold(ribbon_width, option_bundle.epsilon, bucket_size);

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
std::size_t
METHOD_HEADER::get_threshold(std::size_t ribbon_width, double eps, std::size_t bucket_size)
{
    return (1 + 2 * option_bundle.epsilon) * bucket_size - 0.5 * sqrt(bucket_size / (1 + option_bundle.epsilon));
};

#undef CLASS_HEADER
#undef METHOD_HEADER
#undef BUILD_FUNCTION_HEADER

} // namespace retrieval::ribbon