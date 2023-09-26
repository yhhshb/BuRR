#include "../include/burr_config.hpp"

namespace std {
ostream &operator<<(ostream &os, const __uint128_t &v) {
    return os << "0x" << ::std::hex << static_cast<uint64_t>(v >> 64)
              << ::std::setw(16) << static_cast<uint64_t>(v) << ::std::dec;
}
} // namespace std