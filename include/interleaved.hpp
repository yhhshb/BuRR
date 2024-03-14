#pragma once

#include "../bundled/biolib/include/packed_vector.hpp"

namespace ribbon::storage {

class interleaved
{
    public:
        interleaved();

        template <typename KeyType, typename ValueType>
        void backpropagation(
            const std::size_t value_width, 
            const std::vector<KeyType>& coefficients, 
            const std::vector<ValueType>& results
        );

    private:
        bit::packed::vector<std::size_t> data;
};

template <typename KeyType, typename ValueType>
void backpropagation( 
    const std::size_t value_width, // this can be < than ValueType's bit size
    const std::vector<KeyType>& coefficients, 
    const std::vector<ValueType>& results
)
{
    if (coefficients.size() != results.size()) throw std::runtime_error("[Backpropagation] different number of coefficients and results");
    const auto ribbon_width = ::bit::size<KeyType>();
    // using CoeffRow = typename BandingStorage::CoeffRow;

    // constexpr auto kCoeffBits = static_cast<Index>(sizeof(CoeffRow) * 8U);
    // constexpr auto kResultBits = SolutionStorage::kResultBits;
    assert(coefficient.size() % ribbon_width == 0); // m *MUST* be a multiple of ribbon_width

    sol->Prepare(coefficients.size());
    const std::size_t num_blocks = sol->GetNumBlocks();
    const std::size_t num_segments = sol->GetNumSegments();

    // We should be utilizing all available segments
    assert(num_segments == num_blocks * kResultBits);

    // sLOG 
    //     << "Backsubstitution: have" << num_blocks << "blocks," 
    //     << num_segments << "segments for" 
    //     << coefficients.size() << "slots.\n" 
    //     << "kResultBits=" << kResultBits;

    // TODO: consider fixed-column specializations with stack-allocated state

    // A column-major buffer of the solution matrix, 
    // containing enough recently-computed solution data for computing the next row
    // (also based on banding data).
    std::vector<KeyType> state;
    state.resize(value_width);

    std::size_t block = num_blocks;
    std::size_t segment = num_segments;
    while (block > 0) {
        --block;
        // sLOG << "Backsubstituting block" << block << "segment" << segment;
        BackSubstBlock(state, value_width, bs, block * ribbon_width);
        segment -= value_width;
        for (std::size_t i = 0; i < value_width; ++i) {
            sol->SetSegment(segment + i, state.at(i));
        }
    }
    assert(block == 0);
    assert(segment == 0);
}

// A helper for InterleavedBackSubst.
template <typename BandingStorage>
inline void BackSubstBlock(typename BandingStorage::CoeffRow *state,
                           typename BandingStorage::Index num_columns,
                           const BandingStorage &bs,
                           typename BandingStorage::Index start_slot) {
    using CoeffRow = typename BandingStorage::CoeffRow;
    using Index = typename BandingStorage::Index;
    using ResultRow = typename BandingStorage::ResultRow;

    constexpr auto kCoeffBits = static_cast<Index>(sizeof(CoeffRow) * 8U);

    for (Index i = start_slot + kCoeffBits; i > start_slot;) {
        --i;
        CoeffRow cr = bs.GetCoeffs(i);
        ResultRow rr = bs.GetResult(i);
        for (Index j = 0; j < num_columns; ++j) {
            // Compute next solution bit at row i, column j (see derivation below)
            CoeffRow tmp = state[j] << 1;
            int bit = rocksdb::BitParity(tmp & cr) ^ ((rr >> j) & 1);
            tmp |= static_cast<CoeffRow>(bit);

            // Now tmp is solution at column j from row i for next kCoeffBits
            // more rows. Thus, for valid solution, the dot product of the
            // solution column with the coefficient row has to equal the result
            // at that column,
            //   BitParity(tmp & cr) == ((rr >> j) & 1)

            // Update state.
            state[j] = tmp;
        }
    }
}

template <typename BandingStorage, typename SolutionStorage>
void InterleavedBackSubst(const BandingStorage &bs, SolutionStorage *sol) {
    using CoeffRow = typename BandingStorage::CoeffRow;
    using Index = typename BandingStorage::Index;

    static_assert(sizeof(Index) == sizeof(typename SolutionStorage::Index), "must be same");
    static_assert(sizeof(CoeffRow) == sizeof(typename SolutionStorage::CoeffRow), "must be same");

    constexpr auto kCoeffBits = static_cast<Index>(sizeof(CoeffRow) * 8U);
    constexpr auto kResultBits = SolutionStorage::kResultBits;

    const Index num_slots = bs.GetNumSlots();
    // num_slots *MUST* be a multiple of kCoeffBits
    assert(num_slots >= kCoeffBits && num_slots % kCoeffBits == 0);
    sol->Prepare(num_slots);

    const Index num_blocks = sol->GetNumBlocks();
    const Index num_segments = sol->GetNumSegments();

    // We should be utilizing all available segments
    assert(num_segments == num_blocks * kResultBits);

    sLOG << "Backsubstitution: have" << num_blocks << "blocks," << num_segments
         << "segments for" << num_slots << "slots, kResultBits=" << kResultBits;

    // TODO: consider fixed-column specializations with stack-allocated state

    // A column-major buffer of the solution matrix, containing enough
    // recently-computed solution data to compute the next solution row
    // (based also on banding data).
    std::unique_ptr<CoeffRow[]> state{new CoeffRow[kResultBits]()};

    Index block = num_blocks;
    Index segment = num_segments;
    while (block > 0) {
        --block;
        sLOG << "Backsubstituting block" << block << "segment" << segment;
        BackSubstBlock(state.get(), kResultBits, bs, block * kCoeffBits);
        segment -= kResultBits;
        for (Index i = 0; i < kResultBits; ++i) {
            sol->SetSegment(segment + i, state[i]);
        }
    }
    // Verify everything processed
    assert(block == 0);
    assert(segment == 0);
}

} // namespace ribbon