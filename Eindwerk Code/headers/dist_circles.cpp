#include "dist_circles.h"
#include <array>
#include <math.h>
double dist_circles(std::array<double,3> c_1, std::array<double,3> c_2){
    //returns the distance between two circle centers. the arrays contain {x_center, y_center, r}
    return sqrt(pow(c_1.front()-c_2.front(),2) + pow(c_1[1]-c_2[1], 2));
}
