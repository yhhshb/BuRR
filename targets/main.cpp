#include <iostream>
#include "../include/build.hpp"

int main(int argc, char* argv[])
{
    auto build_parser = get_parser_build();
    argparse::ArgumentParser program(argv[0]);
    program.add_subparser(build_parser);
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << "\n";
        std::cerr << program;
        return 1;
    }
    if (program.is_subcommand_used(build_parser)) return build_main(build_parser);
    else std::cerr << program << "\n";
    return 0;
}
