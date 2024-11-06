#ifndef BUILD_HPP
#define BUILD_HPP

#include <argparse/argparse.hpp>

argparse::ArgumentParser get_parser_build();
int build_main(const argparse::ArgumentParser& parser);

#endif // BUILD_HPP