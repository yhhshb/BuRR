#include "../include/build.hpp"
#include "../include/BuRR.hpp"
#include "../include/interleaved.hpp"
#include "../bundled/biolib/include/hash.hpp"

using namespace retrieval::ribbon;

option_t check_args(const argparse::ArgumentParser& parser);

int build_main(const argparse::ArgumentParser& parser)
{
    auto opts = check_args(parser);
    BuRR<threshold_t::TWOBITS, storage::interleaved, hash::hash64> rds(opts);
    // rds.build(start, stop);
    return 0;
}

argparse::ArgumentParser get_parser_build()
{
    argparse::ArgumentParser parser("build");
    parser.add_description("Test BuRR data structure");
    parser.add_argument("-e", "--epsilon")
        .help("over-dimensioning parameter for the linear system. m = (1+e)n.")
        .scan<'f', double>()
        .required();
    parser.add_argument("-l", "--layers")
        .help("number of layers for bumping [4]")
        .scan<'d', std::size_t>()
        .default_value(std::size_t(4));
    parser.add_argument("-s", "--seed")
        .help("random seed")
        .scan<'d', uint64_t>()
        .default_value(uint64_t(42));
    parser.add_argument("-d", "--tmp-dir")
        .help("temporary directory")
        .default_value(std::string("."));
    parser.add_argument("-g", "--max-ram")
        .help("RAM limit (GB) [4]")
        .scan<'d', std::size_t>()
        .default_value(std::size_t(4));
    parser.add_argument("-C", "--check")
        .help("check correctness of BuRR data structures")
        .implicit_value(true)
        .default_value(false);
    parser.add_argument("-v", "--verbose")
        .help("verbosity level [0]")
        .scan<'d', std::size_t>()
        .default_value(std::size_t(0));
    return parser;
}

option_t check_args(const argparse::ArgumentParser& parser)
{
    option_t opts;
    opts.epsilon = parser.get<double>("-e");
    opts.nlayers = parser.get<std::size_t>("-l");
    opts.seed = parser.get<uint64_t>("-s");
    opts.tmp_dir = parser.get<std::string>("-d");
    opts.max_ram = parser.get<std::size_t>("-g");
    opts.nthreads = parser.get<std::size_t>("-t");
    opts.check = parser.get<bool>("-C");
    opts.verbose = parser.get<std::size_t>("-v");
    if (opts.epsilon < 0 or opts.epsilon > 1) throw std::invalid_argument("epsilon must be in [0, 1]");
    if (opts.nlayers > 4) std::cerr << "Warning: are " << opts.nlayers << " layers really necessary?\n";
    return opts;
}