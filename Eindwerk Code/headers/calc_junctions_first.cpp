#include "calc_junctions_first.h"
#include "dist_circles.h"
using namespace std;

array<array<double, 2>, 4> calc_junctions_first(valarray<array<double, 3>> circles, double r) {
    //same as function above, but now only for the first circle
    
    //determine circle center coörds
    double x_i = circles[0].front();    //get x value of circle center
    double y_i = circles[0][1];         //get y value of circle center
    double x_next = circles[1].front();
    double y_next = circles[1][1];
    
    //calculate junctions
    double d = dist_circles(circles[0], circles[1]);
    double alpha = asin((y_next - y_i)/d);
    double M_x = x_i + cos(alpha)*d*0.5;
    double M_y = y_i + sin(alpha)*d*0.5;

    double junction_1_x = 0;
    double junction_1_y = r;   
    double junction_2_x = 0;
    double junction_2_y = -r;
    double junction_3_x = M_x - sin(alpha)*sqrt(pow(r,2) - 0.25*pow(d,2));
    double junction_3_y = M_y + cos(alpha)*sqrt(pow(r,2) - 0.25*pow(d,2));
    double junction_4_x = M_x + sin(alpha)*sqrt(pow(r,2) - 0.25*pow(d,2));
    double junction_4_y = M_y - cos(alpha)*sqrt(pow(r,2) - 0.25*pow(d,2));

    //put the coörds in arrays
    array<double,2> junction_1 {junction_1_x, junction_1_y};
    array<double,2> junction_2 {junction_2_x, junction_2_y};
    array<double,2> junction_3 {junction_3_x, junction_3_y};
    array<double,2> junction_4 {junction_4_x, junction_4_y};

    //put the arrays in an array, to be returned
    array<array<double, 2>, 4> junctions {junction_1, junction_2, junction_3, junction_4};
    return junctions;
}
