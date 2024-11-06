#ifndef BURR_HPP
#define BURR_HPP

#include <cmath>
#include <chrono>
#include <optional>
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

namespace retrieval::ribbon {

enum class threshold_t {NORMAL, ONEBIT, TWOBITS};

class option_t 
{
    public:
        double epsilon;
        std::size_t nlayers;
        uint64_t seed;
        std::string tmp_dir;
        std::size_t max_ram;
        std::size_t nthreads;
        std::size_t verbose;
        bool check;
};

template <enum threshold_t, class SolutionStorage, class Hasher>
class BuRR {}; // dummy class for separate specializations

} // namespace retrieval::ribbon

#include "BuRR1bit.hpp"
#include "BuRR2bits.hpp"

#endif // BURR_HPP
