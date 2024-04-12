#pragma once

#include <optional>
#include "../bundled/biolib/include/packed_vector.hpp"

namespace ribbon::storage {

class interleaved
{
    public:
        interleaved();

        template <typename KeyType, typename ValueType>
        void 
        backsubstitution(
            const std::size_t value_width, 
            const std::vector<KeyType>& coefficients, 
            const std::vector<ValueType>& results
        );

        std::size_t 
        value_width() const noexcept;

        template <typename KeyType, typename ValueType>
        std::optional<ValueType>
        at(const std::size_t offset, const KeyType hash) const noexcept;

        template <typename KeyType>
        KeyType 
        get_coefficients(const KeyType hkey);

    private:
        static const std::size_t salt = 0xc28f82822b650bedULL;
        std::size_t valwidth;
        bit::packed::vector<std::size_t> data;
};

template <typename KeyType>
KeyType 
interleaved::get_coefficients(const KeyType hkey)
{
    std::size_t target_pcnt = 0;
    if constexpr (::bit::size<KeyType>() == 128) target_pcnt = 16;
    else if (::bit::size<KeyType>() == 16) target_pcnt = 4;
    else target_pcnt = 8;
    assert(target_pcnt != 0);
    assert(::bit::size<KeyType>() / target_pcnt == 0);
    const std::size_t step = ::bit::size<KeyType>() / target_pcnt;
    const KeyType mask = step - 1;
    const std::size_t shift = std::ceil(std::log2(step));//tlx::integer_log2_ceil(step_);
    auto a = hkey * salt;
    KeyType cr = 0;
    for (std::size_t i = 0; i < target_pcnt; ++i) {
        std::size_t pos = (a & mask) + step * i;
        cr |= static_cast<KeyType>(1) << pos;
        a >>= shift;
    }

    // Now ensure the value is non-zero
    // if constexpr (kFirstCoeffAlwaysOne) {
        cr |= 1;
    //} ...
    return cr;
}

template <typename KeyType, typename ValueType>
void 
interleaved::backsubstitution(
    const std::size_t value_width, // this can be < than ValueType's bit size
    const std::vector<KeyType>& coefficients, // FIXME put coefficient and result in a pair for efficiency
    const std::vector<ValueType>& results
)
{
    if (coefficients.size() != results.size()) throw std::runtime_error("[Backpropagation] different number of coefficients and results");
    valwidth = value_width;
    const auto ribbon_width = ::bit::size<KeyType>();
    assert(coefficients.size() % ribbon_width == 0); // m *MUST* be a multiple of ribbon_width

    // A column-major buffer of the solution matrix, 
    // containing enough recently-computed solution data for computing the next row
    // (also based on banding data).
    std::vector<KeyType> state;
    state.resize(value_width);

    std::size_t block = coefficients.size() / ribbon_width; // here block starts as the total number of blocks
    std::size_t segment = block * value_width;
    assert(segment == block * value_width); // We should be utilizing all available segments
    data = bit::packed::vector(ribbon_width);
    while (block > 0) {
        --block;
        std::size_t i = (block + 1) * ribbon_width;
        while (i > block * ribbon_width) {
            --i;
            KeyType cr = coefficients.at(i);
            ValueType rr = results.at(i);
            for (std::size_t j = 0; j < value_width; ++j) {
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
        }
    }
    assert(block == 0);
    assert(segment == 0);
}

template <typename KeyType, typename ValueType>
std::optional<ValueType> 
interleaved::at(const std::size_t offset, const KeyType hash) const noexcept
{
    const auto ribbon_width = ::bit::size<KeyType>();
    const std::size_t start_block_num = offset / ribbon_width;
    const std::size_t start_bit = offset % ribbon_width;
    const std::size_t cr = get_coefficients(hash);
    const std::size_t cr_left = cr << start_bit;

    ValueType retrieved = 0;
    std::size_t segment = start_block_num * value_width();
    for (std::size_t i = 0; i < value_width(); ++i) {
        retrieved ^= bit::parity(static_cast<std::size_t>(data.at(segment + i)) & cr_left) << i;
    }

    if (start_bit > 0) {
        segment += value_width();
        const std::size_t cr_right = cr >> (ribbon_width - start_bit);
        for (std::size_t i = 0; i < value_width(); ++i) {
            retrieved ^= bit::parity(static_cast<std::size_t>(data.at(segment + i)) & cr_right) << i;
        }
    }
    return retrieved;
}

} // namespace ribbon