#include "calc_wire_resistance.h"
#include "dist.h"
#include <math.h>
using namespace std;
double calc_wire_resistance(array<double,3> circle, array<double, 2> point_1, array<double, 2> point_2, double R) {
    double r = circle.back();
    return (R*acos(1 - 0.5*pow(dist(point_1, point_2)/r,2)));
}