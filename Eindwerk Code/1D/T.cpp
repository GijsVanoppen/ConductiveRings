#include "T.h"
std::valarray<std::array<double, 3>> T::generate_circles(int N, double a, double r) {
    
    //rng
    using seed_dist_t = std::uniform_int_distribution<size_t>;
    std::random_device dev;
    seed_dist_t seed_distr(0, std::numeric_limits<size_t>::max());
    auto seed = seed_distr(dev);
    std::mt19937_64 engine(seed);
    auto distr = std::bind(std::uniform_real_distribution<double>(-1.0, 1.0),engine);
    
    double deviation_coefficient = -a/2 + 0.5*sqrt(8*pow(r,2) - pow(a,2));  //determines the max value the circles may deviate from regular lattice
    
    std::valarray<std::array<double, 3>> circles(N);
    for (int i {1}; i < N-1; i++) {
        circles[i][0] = a*i + distr()*deviation_coefficient;      //x values of circle center, with a small deviation from regular lattice (distr())
        circles[i][1] = distr()*deviation_coefficient;            //y values of circle center, with a small deviation from regular lattice (distr())
        circles[i][2] = r;                                        //circle radius
    } 
    circles[0][2] = r;              //first circle, no deviations from regular lattice
    circles[N-1][0] = a*(N-1);      //last circle, no deviations from regular lattice
    circles[N-1][2] = r;
    return circles;
}
