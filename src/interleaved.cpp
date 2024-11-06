#include "../include/interleaved.hpp"

namespace retrieval::ribbon::storage {

interleaved::interleaved()
    : valwidth(0), data(0)
{}

std::size_t 
interleaved::value_width() const noexcept
{
    return valwidth;
}

} // namespace retrieval::ribbon::storage