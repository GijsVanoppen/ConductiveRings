#ifndef CALC_JUNCTIONS_LAST_H
#define CALC_JUNCTIONS_LAST_H
#include <array>
#include <valarray>
std::array<std::array<double, 2>, 4> calc_junctions_last(std::valarray<std::array<double, 3>> circles, double r);
#endif