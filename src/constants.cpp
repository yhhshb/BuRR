#include "../include/constants.hpp"

namespace ribbon::util {

std::string get_name(const std::string& prefix, uint64_t run_id)
{
    return prefix + "_" + std::to_string(run_id);
}

}