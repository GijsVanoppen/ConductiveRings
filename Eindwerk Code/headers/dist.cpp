#include "generate_circles.h"
double dist(std::array<double,2> point_1, std::array<double,2> point_2){
    //returns the distance between two points
    return sqrt(pow(point_1.front()-point_2.front(),2) + pow(point_1.back()-point_2.back(), 2));
}