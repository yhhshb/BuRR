//  Copyright (c) Lorenz Hübschle-Schneider
//  All Rights Reserved.  This source code is licensed under the Apache 2.0
//  License (found in the LICENSE file in the root directory).



#include <atomic>
#include <cstdlib>
#include <numeric>
#include <thread>
#include <vector>


#include "../include/ribbon.hpp"
#include "../include/rocksdb/stop_watch.h"

#include "../bundled/tlx/include/cmdline_parser.hpp"
#include "../bundled/tlx/include/logger.hpp"

#define DO_EXPAND(VAL) VAL##1
#define EXPAND(VAL)    DO_EXPAND(VAL)

#if !defined(RIBBON_BITS) || (EXPAND(RIBBON_BITS) == 1)
#undef RIBBON_BITS
#define RIBBON_BITS 8
#endif

namespace ribbon {

bool no_queries = false;

template <uint8_t depth, typename Config>
void run(size_t num_slots, double eps, size_t seed, unsigned num_threads) {
    IMPORT_RIBBON_CONFIG(Config);

    const double slots_per_item = eps + 1.0;
    const size_t num_items = num_slots / slots_per_item;
    LOG1 << "Running simple test with " << num_slots << " slots, eps=" << eps
         << " -> " << num_items << " items, seed=" << seed
         << " config: L=" << kCoeffBits << " B=" << kBucketSize
         << " r=" << kResultBits;

    rocksdb::StopWatchNano timer(true);

    auto input = std::make_unique<int[]>(num_items);
    std::iota(input.get(), input.get() + num_items, 0);
    LOG1 << "Input generation took " << timer.ElapsedNanos(true) / 1e6 << "ms";

    ribbon_filter<depth, Config> r(num_slots, slots_per_item, seed);

    LOG1 << "Allocation took " << timer.ElapsedNanos(true) / 1e6 << "ms\n";

    LOG1 << "Adding rows to filter....";
    r.AddRange(input.get(), input.get() + num_items);
    LOG1 << "Insertion took " << timer.ElapsedNanos(true) / 1e6 << "ms in total\n";

    input.reset();

    r.BackSubst();
    LOG1 << "Backsubstitution took " << timer.ElapsedNanos(true) / 1e6
         << "ms in total\n";

    const size_t bytes = r.Size();
    const double relsize = (bytes * 8 * 100.0) / (num_items * Config::kResultBits);
    LOG1 << "Ribbon size: " << bytes << " Bytes = " << (bytes * 1.0) / num_items
         << " Bytes per item = " << relsize << "%\n";

    std::atomic<bool> ok = true;
    auto pos_query = [&r, &ok, &num_items, &num_threads](unsigned id) {
        bool my_ok = true;
        size_t start = num_items / num_threads * id;
        // don't do the same queries on all threads
        for (size_t v = start; v < num_items; v++) {
            bool found = r.QueryFilter((int)v);
            assert(found);
            my_ok &= found;
        }
        for (size_t v = 0; v < start; v++) {
            bool found = r.QueryFilter((int)v);
            assert(found);
            my_ok &= found;
        }
        if (!my_ok)
            ok = false;
    };
    std::vector<std::thread> threads;
    for (unsigned i = 0; i < num_threads && !no_queries; i++) {
        threads.emplace_back(pos_query, i);
    }
    for (auto& t : threads)
        t.join();

    const auto check_nanos = timer.ElapsedNanos(true);
    LOG1 << "Parallel check with " << num_threads << " threads "
         << (ok ? "successful" : "FAILED") << " and took " << check_nanos / 1e6
         << "ms = " << check_nanos * 1.0 / num_items << "ns per key";
    // r.PrintStats();

    std::atomic<size_t> found = 0;
    auto neg_query = [&r, &found, &num_items, &num_threads](unsigned id) {
        size_t my_found = 0;
        // offset queries between threads
        size_t start = num_items + num_items / num_threads * id;
        for (size_t v = start; v < 2 * num_items; v++) {
            my_found += r.QueryFilter((int)v);
        }
        for (size_t v = num_items; v < start; v++) {
            my_found += r.QueryFilter((int)v);
        }
        found.fetch_add(my_found);
    };
    threads.clear();
    for (unsigned i = 0; i < num_threads && !no_queries; i++) {
        threads.emplace_back(neg_query, i);
    }
    for (auto& t : threads)
        t.join();

    const auto negq_nanos = timer.ElapsedNanos(true);
    const double fprate = found * 1.0 / (num_threads * num_items),
                 ratio = fprate * (1ul << Config::kResultBits);
    LOG1 << "Negative queries took " << negq_nanos / 1e6
         << "ms = " << negq_nanos * 1.0 / num_items << "ns per key, " << found
         << " FPs = " << fprate * 100 << "%, expecting "
         << 100.0 / (1ul << Config::kResultBits) << "% -> ratio = " << ratio;
    // r.PrintStats();
    auto [tl_bumped, tl_empty_slots, tl_frac_empty, tl_thresh_bytes] =
        r.GetStats();
    LOG1 << "RESULT n=" << num_items << " m=" << num_slots << " eps=" << eps
         << " d=" << (int)depth << dump_config<Config>() << " bytes=" << bytes
         << " tlempty=" << tl_empty_slots << " tlbumped=" << tl_bumped
         << " tlemptyfrac=" << tl_frac_empty
         << " tlthreshbytes=" << tl_thresh_bytes << " overhead=" << relsize - 100
         << " ok=" << ok << " tpos=" << check_nanos
         << " tpospq=" << (check_nanos * 1.0 / num_items) << " tneg=" << negq_nanos
         << " tnegpq=" << (negq_nanos * 1.0 / num_items) << " fps=" << found
         << " fpr=" << fprate << " ratio=" << ratio << " threads=" << num_threads;
}


template <ThreshMode mode, uint8_t depth, size_t L, size_t r, bool interleaved,
          bool cls, bool sparse, typename... Args>
void dispatch_shift(int shift, Args&... args) {
    switch (shift) {
        case 0:
            run<depth, RConfig<L, r, mode, sparse, interleaved, cls, 0>>(args...);
            break;
        case -1:
            run<depth, RConfig<L, r, mode, sparse, interleaved, cls, -1>>(args...);
            break;
        case 1:
            run<depth, RConfig<L, r, mode, sparse, interleaved, cls, 1>>(args...);
            break;
        default: LOG1 << "Unsupported bucket size shift: " << shift;
    }
}


template <ThreshMode mode, uint8_t depth, size_t L, size_t r, bool interleaved,
          bool cls, typename... Args>
void dispatch_sparse(bool sparse, Args&... args) {
    if (sparse) {
        if constexpr (interleaved) {
            LOG1 << "Sparse coefficients + interleaved sol doesn't make sense";
        } else {
            dispatch_shift<mode, depth, L, r, interleaved, cls, true>(args...);
        }
    } else {
        dispatch_shift<mode, depth, L, r, interleaved, cls, false>(args...);
    }
}

template <ThreshMode mode, uint8_t depth, size_t L, size_t r, typename... Args>
void dispatch_storage(bool cls, bool interleaved, Args&... args) {
    assert(!cls || !interleaved);
    if (cls) {
        // dispatch_sparse<mode, depth, L, r, false, true>(args...);
        LOG1 << "Cache-Line Storage is currently disabled";
    } else if (interleaved) {
        dispatch_sparse<mode, depth, L, r, true, false>(args...);
    } else {
        dispatch_sparse<mode, depth, L, r, false, false>(args...);
    }
}

template <ThreshMode mode, uint8_t depth, typename... Args>
void dispatch_width(size_t band_width, Args&... args) {
    static constexpr size_t r = RIBBON_BITS;
    switch (band_width) {
        // case 16: dispatch_storage<mode, depth, 16, r>(args...); break;
        case 32: dispatch_storage<mode, depth, 32, r>(args...); break;
        case 64: dispatch_storage<mode, depth, 64, r>(args...); break;
        case 128: dispatch_storage<mode, depth, 128, r>(args...); break;
        default: LOG1 << "Unsupported band width: " << band_width;
    }
}

template <ThreshMode mode, typename... Args>
void dispatch_depth(unsigned depth, Args&... args) {
    switch (depth) {
        case 0: dispatch_width<mode, 0>(args...); break;
        case 1: dispatch_width<mode, 1>(args...); break;
        case 2: dispatch_width<mode, 2>(args...); break;
        case 3: dispatch_width<mode, 3>(args...); break;
        case 4: dispatch_width<mode, 4>(args...); break;
        default: LOG1 << "Unsupported recursion depth: " << depth;
    }
}

template <typename... Args>
void dispatch(ThreshMode mode, Args&... args) {
    switch (mode) {
        case ThreshMode::onebit:
            dispatch_depth<ThreshMode::onebit>(args...);
            break;
        case ThreshMode::twobit:
            dispatch_depth<ThreshMode::twobit>(args...);
            break;
        case ThreshMode::normal:
            dispatch_depth<ThreshMode::normal>(args...);
            break;
        default:
            LOG1 << "Unsupported threshold compression mode: " << (int)mode;
    }
}

} // namespace ribbon

//----------------------------------------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include "../include/build.hpp"

using namespace ribbon;

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
