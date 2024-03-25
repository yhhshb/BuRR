#include "../include/BuRR.hpp"

namespace ribbon {

CLASS_HEADER
METHOD_HEADER::BuRR(const option_t& build_options)
    : 
    option_bundle(build_options), 
    total_empty_slots(0),
    bucket_size(0),
    lower(0),
    upper(0)
{}

CLASS_HEADER
METHOD_HEADER::layer::layer() : bump_info(bit::packed::vector<std::size_t>(2)) // 2-bit variant
{}

CLASS_HEADER
std::tuple<std::size_t, std::size_t>
METHOD_HEADER::get_thresholds(std::size_t ribbon_width, double eps, std::size_t bucket_size)
{
    auto compute_thresholds = [&eps, &bucket_size](std::tuple<double, double, double, double> coefficients) {
        auto [la, lb, ua, ub] = coefficients;
        return std::make_tuple(
            static_cast<std::size_t>((la * eps + lb) * bucket_size), 
            static_cast<std::size_t>((ua * eps + ub) * bucket_size)
        );
    };
    std::tuple<double, double, double, double> coefficients;
    if (ribbon_width <= 16) coefficients = std::make_tuple(0.5, 0.6, 0.5, 0.82);
    else if (ribbon_width <= 32) coefficients = std::make_tuple(0.5, 0.7, 0.5, 0.87);
    else coefficients = std::make_tuple(1.30, 0.78, 0.75, 0.91);
    return compute_thresholds(coefficients);
};

} // namespace ribbon