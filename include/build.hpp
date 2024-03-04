#ifndef BUILD_HPP
#define BUILD_HPP

#include <argparse/argparse.hpp>

namespace ribbon {

argparse::ArgumentParser get_parser_build();
int build_main(const argparse::ArgumentParser& parser);

} // namespace ribbon

#endif // BUILD_HPP