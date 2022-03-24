#include <iostream>
#include <valarray>
using namespace std;
int main(){
    valarray<double> x(1);
    valarray<valarray<double>> y;
    y[1] = x;
    cout << x[0];

}





