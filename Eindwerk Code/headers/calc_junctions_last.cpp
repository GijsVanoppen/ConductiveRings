#include "calc_junctions_last.h"
#include "dist_circles.h"
using namespace std;

array<array<double, 2>, 4> calc_junctions_last(valarray<array<double, 3>> circles, double r) {
    //same as function above, but now only for the last circle
    
    //determine circle center coörds
    int N = circles.size();
    double x_i = circles[N-1].front();    //get x value of circle center
    double y_i = circles[N-1][1];         //get y value of circle center
    double x_prev = circles[N-2].front();
    double y_prev = circles[N-2][1];
    
    //calculate junctions
    double d = dist_circles(circles[N-2], circles[N-1]);
    double alpha = asin((y_i - y_prev)/d);
    double M_x = x_prev + cos(alpha)*d*0.5;
    double M_y = y_prev + sin(alpha)*d*0.5;

    double junction_1_x = M_x - sin(alpha)*sqrt(pow(r,2) - 0.25*pow(d,2));
    double junction_1_y = M_y + cos(alpha)*sqrt(pow(r,2) - 0.25*pow(d,2));
    double junction_2_x = M_x + sin(alpha)*sqrt(pow(r,2) - 0.25*pow(d,2));
    double junction_2_y = M_y - cos(alpha)*sqrt(pow(r,2) - 0.25*pow(d,2));
    double junction_3_x = circles[N-1].front();
    double junction_3_y = r;   
    double junction_4_x = junction_3_x;
    double junction_4_y = -r;

    //put the coörds in arrays
    array<double,2> junction_1 {junction_1_x, junction_1_y};
    array<double,2> junction_2 {junction_2_x, junction_2_y};
    array<double,2> junction_3 {junction_3_x, junction_3_y};
    array<double,2> junction_4 {junction_4_x, junction_4_y};

    //put the arrays in an array, to be returned
    array<array<double, 2>, 4> junctions {junction_1, junction_2, junction_3, junction_4};
    return junctions;
}