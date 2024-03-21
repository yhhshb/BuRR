#pragma once

#include <optional>
#include "../bundled/biolib/include/packed_vector.hpp"

namespace ribbon::storage {

class interleaved
{
    public:
        interleaved();

        template <typename KeyType, typename ValueType>
        void backsubstitution(
            const std::size_t value_width, 
            const std::vector<KeyType>& coefficients, 
            const std::vector<ValueType>& results
        );

        std::size_t value_width() const noexcept;

        template <typename ValueType>
        std::optional<ValueType>
        at(const std::size_t offset) const noexcept;

    private:
        std::size_t valwidth;
        bit::packed::vector<std::size_t> data;
};

template <typename KeyType, typename ValueType>
void 
backsubstitution(
    const std::size_t value_width, // this can be < than ValueType's bit size
    const std::vector<KeyType>& coefficients, 
    const std::vector<ValueType>& results
)
{
    if (coefficients.size() != results.size()) throw std::runtime_error("[Backpropagation] different number of coefficients and results");
    valwidth = value_width;

    // using CoeffRow = typename BandingStorage::CoeffRow;
    // constexpr auto kCoeffBits = static_cast<Index>(sizeof(CoeffRow) * 8U);
    const auto ribbon_width = ::bit::size<KeyType>();
    assert(coefficient.size() % ribbon_width == 0); // m *MUST* be a multiple of ribbon_width
    // constexpr auto kResultBits = SolutionStorage::kResultBits;

    // TODO: consider fixed-column specializations with stack-allocated state

    // A column-major buffer of the solution matrix, 
    // containing enough recently-computed solution data for computing the next row
    // (also based on banding data).
    std::vector<KeyType> state;
    state.resize(value_width);

    std::size_t block = coefficients.size() / ribbon_width;
    std::size_t segment = num_blocks * value_width;
    assert(segment == block * value_width); // We should be utilizing all available segments
    data = bit::packed::vector(ribbon_width);
    while (block > 0) {
        --block;
        std::size_t i = (block + 1) * ribbon_width;
        while (i > block * ribbon_width) {
            --i;
            KeyType cr = bs.GetCoeffs(i);
            ValueType rr = bs.GetResult(i);
            for (std::size_t j = 0; j < num_columns; ++j) {
                // Compute next solution bit at row i, column j (see derivation below)
                KeyType tmp = state.at(j) << 1;
                auto bit = bit::parity(tmp & cr) ^ ((rr >> j) & static_cast<ValueType>(1));
                assert(bit >= 0);
                tmp |= static_cast<KeyType>(bit);

                // Now tmp is solution at column j from row i for next ribbon_width more rows. 
                // Thus, for valid solution, the dot product of the solution column with the 
                // coefficient row has to equal the result at that column,
                // bit::parity(tmp & cr) == ((rr >> j) & 1)

                // Update state.
                state[j] = tmp;
            }
        }
        segment -= value_width;
        for (std::size_t i = 0; i < value_width; ++i) {
            data[segment + i] = state.at(i); 
            // memcpy(data_.get() + (segment + i) * sizeof(CoeffRow), state.at(i), sizeof(KeyType));
        }
    }
    assert(block == 0);
    assert(segment == 0);
}

template <typename ValueType>
std::optional<ValueType> 
interleaved::at(const std::tuple<std::size_t, std::size_t, typename METHOD_HEADER::bumping, typename Hasher::hash_t>& pack) const noexcept
{
    const auto ribbon_width = ::bit::size<Hasher::hash_t>();
    const auto [offset, bucket, cval, hash] = pack;

    const std::size_t start_block_num = offset / ribbon_width;
    const std::size_t segment = start_block_num * Z.value_width();
    const std::size_t start_bit = offset % ribbon_width;

    const std::size_t cr = hasher.GetCoeffs(hash); // TODO check usage of GetCoeff from hasher.hpp
    const std::size_t cr_left = cr << start_bit;

    ValueType retrieved = 0;
    for (std::size_t i = 0; i < value_width(); ++i) {
        retrieved ^= bit::parity(iss.GetSegment(segment + i) & cr_left) << i;
    }

    if (start_bit > 0) {
        segment += value_width();
        const std::size_t cr_right = cr >> (kCoeffBits - start_bit);
        for (std::size_t i = 0; i < value_width(); ++i) {
            retrieved ^= bit::parity(iss.GetSegment(segment + i) & cr_right) << i;
        }
    }
}

} // namespace ribbon