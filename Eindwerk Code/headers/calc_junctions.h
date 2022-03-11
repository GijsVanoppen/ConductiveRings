#ifndef CALC_JUNCTIONS_H
#define CALC_JUNCTIONS_H
#include "dist_circles.h"
#include <array>
#include <valarray>
std::array<std::array<double, 2>, 4> calc_junctions(int i, std::valarray<std::array<double, 3>> circles, double r);
#endif