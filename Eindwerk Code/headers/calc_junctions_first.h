#ifndef CALC_JUNCTIONS_FIRST_H
#define CALC_JUNCTIONS_FIRST_H
#include <array>
#include <valarray>
std::array<std::array<double, 2>, 4> calc_junctions_first(std::valarray<std::array<double, 3>> circles, double r);
#endif