#pragma once
#include <valarray>
#include <array>
#include <random>
class T {
    public:
    std::valarray<std::array<double, 3>> generate_circles(int N, double a, double r);
};