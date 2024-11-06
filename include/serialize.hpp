/*
This method should be replaced by biolib's own endianness-aware functions
*/

template <typename T>
bool bswap_generic(T& val) {
    if constexpr (sizeof(T) == 1) {
        // Do nothing
    } else if constexpr (sizeof(T) == 2) {
        val = __builtin_bswap16(val);
    } else if constexpr (sizeof(T) == 4) {
        val = __builtin_bswap32(val);
    } else if constexpr (sizeof(T) == 8) {
        val = __builtin_bswap64(val);
    } else if constexpr (sizeof(T) == 16) {
        val = __builtin_bswap128(val);
    } else {
        return false;
    }
    return true;
}