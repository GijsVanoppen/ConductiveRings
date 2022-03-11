#include "calc_junctions.h"
std::array<std::array<double, 2>, 4> calc_junctions(int i, std::valarray<std::array<double, 3>> circles, double r) {
    //returns the coördinates of the junctions that the i-th circle with radius r makes with the other circles in the valarray 
    //junctions are ordered: left up (1), left down (2), right up (3), right down (4) as seen form the circles center.
    
    //determine circle center coörds
    double x_i = circles[i].front();    //get x value of circle center
    double y_i = circles[i][1];         //get y value of circle center
    double x_prev = circles[i-1].front();
    double y_prev = circles[i-1][1];
    double x_next = circles[i+1].front();
    double y_next = circles[i+1][1];

    //between first 2 circles
    double d = dist_circles(circles[i-1], circles[i]);
    double alpha = asin((y_i - y_prev)/d);
    double M_x = x_prev + cos(alpha)*d*0.5; 
    double M_y = y_prev + sin(alpha)*d*0.5;

    double junction_1_x = M_x - sin(alpha)*sqrt(pow(r,2) - 0.25*pow(d,2));
    double junction_1_y = M_y + cos(alpha)*sqrt(pow(r,2) - 0.25*pow(d,2));
    double junction_2_x = M_x + sin(alpha)*sqrt(pow(r,2) - 0.25*pow(d,2));
    double junction_2_y = M_y - cos(alpha)*sqrt(pow(r,2) - 0.25*pow(d,2));
    
    //between last 2 circles
    d = dist_circles(circles[i], circles[i+1]);
    alpha = asin((y_next - y_i)/d);
    M_x = x_i + cos(alpha)*d*0.5;
    M_y = y_i + sin(alpha)*d*0.5;

    double junction_3_x = M_x - sin(alpha)*sqrt(pow(r,2) - 0.25*pow(d,2));
    double junction_3_y = M_y + cos(alpha)*sqrt(pow(r,2) - 0.25*pow(d,2));
    double junction_4_x = M_x + sin(alpha)*sqrt(pow(r,2) - 0.25*pow(d,2));
    double junction_4_y = M_y - cos(alpha)*sqrt(pow(r,2) - 0.25*pow(d,2));

    //put the coörds in arrays
    std::array<double,2> junction_1 {junction_1_x, junction_1_y};
    std::array<double,2> junction_2 {junction_2_x, junction_2_y};
    std::array<double,2> junction_3 {junction_3_x, junction_3_y};
    std::array<double,2> junction_4 {junction_4_x, junction_4_y};

    //put the arrays in an array, to be returned
    std::array<std::array<double, 2>, 4> junctions {junction_1, junction_2, junction_3, junction_4};
    return junctions;
}
