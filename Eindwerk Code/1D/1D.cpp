#include <iostream>
#include <array>
#include <valarray>
#include <cmath>
#include <ctgmath>
#include <C:\Python\Python38-32\Scripts\Toolbox\eigen\eigen-3.4.0\Eigen\Sparse>
#include <fstream>
#include <string>
#include <random>


// #include "..\headers\Mat.cpp"
#include "..\headers\generate_circles.cpp"
#include "..\headers\dist.h"
#include "..\headers\dist_circles.h"
#include "..\headers\calc_junctions.h"
#include "..\headers\calc_junctions_first.h"
#include "..\headers\calc_junctions_last.h"
#include "..\headers\calc_wire_resistance.h"
#include "..\headers\write_circles_to_file.cpp"
#include "..\headers\write_results_to_file.h"
//#include "..\headers\Mat.h"

using namespace std;





/*
class Mat {
    private:
    int nr_rows_;
    int nr_cols_;
    std::valarray<double> row;
    std::valarray<std::valarray<double>> matrix_;
    std::valarray<double> I;

    public:
    Mat(const int nr_rows, const int nr_cols):
        nr_rows_ {nr_rows},
        nr_cols_ {nr_cols},
        row (nr_cols),
        matrix_ (nr_rows),
        I (nr_rows) {};
    
    double& operator()(const int row, const int col) { return matrix_[row][col];}
    int nr_rows() const { return nr_rows_; }
    int nr_cows() const { return nr_cols_; }

    void print(){
        for (int i {0}; i < nr_cols_; i++) {
            for (int j {0}; j < nr_rows_; j++){
                std::cout << matrix_[i][j] << " ";
            }
            std::cout << "\n";
        }
    }

    Eigen::VectorXd solve_matrix(std::valarray<double> b_) {
        Eigen::SparseMatrix<double> A(nr_rows_, nr_cols_);
        Eigen::VectorXd x (nr_cols_);
        Eigen::VectorXd b (nr_cols_);
        for (int i {0}; i < nr_rows_; i++){
            for (int j {0}; j < nr_cols_; j++){
                A.insert(i, j) = matrix_[i][j];
            }
        } 
        for (int i {0}; i < b_.size(); i++) {
            b(i) = b_[i];
        }
        A.makeCompressed();
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;   
        solver.analyzePattern(A);
        solver.factorize(A);
        x = solver.solve(b);
        return x;
    }


    void build_matrix(std::valarray<std::array<double, 3>> circles, double r, double R, double R_j){
        double R_ver;   //'vertical' resistance
        double R_hor;   //'horizontal' resistance
        //---FIRST CIRCLE---

        auto junctions = calc_junctions_first(circles, r);
        //node right up
        R_ver = calc_wire_resistance(circles[0], junctions[2], junctions[3], R);
        R_hor = calc_wire_resistance(circles[0], junctions[0], junctions[2], R);
        matrix_[0][0] = R_ver*R_j + R_hor*R_j + R_hor*R_ver;
        matrix_[0][1] = -R_hor*R_j;
        matrix_[0][2] = -R_hor*R_ver;

        //node right down
        R_hor = calc_wire_resistance(circles[0], junctions[1], junctions[3], R);
        matrix_[1][0] = -R_hor*R_j;
        matrix_[1][1] = R_ver*R_j + R_hor*R_j + R_hor*R_ver;
        matrix_[1][3] = -R_hor*R_ver;
        
        I[0] = R_ver*R_j;
        I[1] = R_ver*R_j;


        //---OTHER CIRCLES---

        for (int i = 1; i < circles.size()-1; i++) {
            junctions = calc_junctions(i, circles, r);
            //node left up
            R_ver = calc_wire_resistance(circles[i], junctions[0], junctions[1], R);
            R_hor = calc_wire_resistance(circles[i], junctions[0], junctions[2], R);
            matrix_[i*4-2][i*4-4] = -R_ver*R_hor;
            matrix_[i*4-2][i*4-2] = R_ver*R_hor + R_j*R_hor + R_j*R_ver;
            matrix_[i*4-2][i*4-1] = -R_j*R_hor;
            matrix_[i*4-2][i*4] -R_j*R_ver;

            //node left down
            R_hor = calc_wire_resistance(circles[i], junctions[1], junctions[3], R);
            matrix_[i*4-1][i*4-3] = -R_ver*R_hor;
            matrix_[i*4-1][i*4-2] = -R_j*R_hor;
            matrix_[i*4-1][i*4-1] = R_ver*R_hor + R_j*R_hor + R_j*R_ver;
            matrix_[i*4-1][i*4+1] = -R_j*R_ver;

            //node right up
            R_ver = calc_wire_resistance(circles[i], junctions[2], junctions[3], R);
            R_hor = calc_wire_resistance(circles[i], junctions[0], junctions[2], R);
            matrix_[i*4][i*4-2] = -R_ver*R_j;
            matrix_[i*4][i*4] = R_ver*R_j + R_hor*R_j + R_hor*R_ver;
            matrix_[i*4][i*4+1] = -R_hor*R_j;
            matrix_[i*4][i*4+2] = -R_hor*R_ver;

            //node right down
            R_hor = calc_wire_resistance(circles[i], junctions[1], junctions[3], R);
            matrix_[i*4+1][i*4-1] = -R_ver*R_j;
            matrix_[i*4+1][i*4] = -R_hor*R_j;
            matrix_[i*4+1][i*4+1] = R_ver*R_j + R_hor*R_j + R_hor*R_ver;
            matrix_[i*4+1][i*4+3] = -R_hor*R_ver;
        }
        //---FINAL CIRCLE---

        junctions = calc_junctions_last(circles, r);
        //node left up
        R_ver = calc_wire_resistance(circles[-1], junctions[0], junctions[1], R);
        R_hor = calc_wire_resistance(circles[-1], junctions[0], junctions[2], R);
        std::cout << R_ver << " " << R_hor << std::endl;
        for (const auto& i: junctions) {
            for (const auto& j: i){
                std::cout << j << " ";
            }
            std::cout << "\n";
        }
        matrix_[nr_rows_-2][nr_rows_-4] = -R_ver*R_hor;
        matrix_[nr_rows_-2][nr_rows_-2] = R_ver*R_hor + R_j*R_hor + R_j*R_ver;
        matrix_[nr_rows_-2][nr_rows_-1] = -R_j*R_hor;
        //node left down
        R_hor = calc_wire_resistance(circles[-1], junctions[1], junctions[3], R);
        matrix_[nr_rows_-1][nr_rows_-3] = -R_ver*R_hor;
        matrix_[nr_rows_-1][nr_rows_-2] = -R_j*R_hor;
        matrix_[nr_rows_-1][nr_rows_-1] = R_ver*R_hor + R_j*R_hor + R_j*R_ver;
    }
};
*/

double calc_wire_resistance(array<double,3> circle, array<double, 2> point_1, array<double, 2> point_2, double R) {
    double r = circle.back();
    return (R*acos(1 - 0.5*pow(dist(point_1, point_2)/r,2)));
}

class Mat {
    private:
    int nr_rows_;
    int nr_cols_;
    std::valarray<double> row;
    std::valarray<std::valarray<double>> matrix_;
    std::valarray<double> I;

    public:
    Mat(const int nr_rows, const int nr_cols):
        nr_rows_ {nr_rows},
        nr_cols_ {nr_cols},
        row (nr_cols),
        matrix_ (nr_rows),
        I (nr_rows) {};  

    double& operator()(const int row, const int col) { return matrix_[row][col];}
    int nr_rows() const { return nr_rows_; }
    int nr_cows() const { return nr_cols_; }
    
    void print(){
        for (int i {0}; i < nr_cols_; i++) {
            for (int j {0}; j < nr_rows_; j++){
                std::cout << matrix_[i][j] << " ";
            }
            std::cout << "\n";
        }
    }

    Eigen::VectorXd solve_matrix(std::valarray<double> b_){
        Eigen::SparseMatrix<double> A(nr_rows_, nr_cols_);
        Eigen::VectorXd x (nr_cols_);
        Eigen::VectorXd b (nr_cols_);
        for (int i {0}; i < nr_rows_; i++){
            for (int j {0}; j < nr_cols_; j++){
                A.insert(i, j) = matrix_[i][j];
            }
        } 
        for (int i {0}; i < b_.size(); i++) {
            b(i) = b_[i];
        }
        A.makeCompressed();
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;   
        solver.analyzePattern(A);
        solver.factorize(A);
        x = solver.solve(b);
        return x;
    }

    void build_matrix(std::valarray<std::array<double, 3>> circles, double r, double R, double R_j){
        double R_ver;   //'vertical' resistance
        double R_hor;   //'horizontal' resistance
        //---FIRST CIRCLE---

        auto junctions = calc_junctions_first(circles, r);
        //node right up
        R_ver = calc_wire_resistance(circles[0], junctions[2], junctions[3], R);
        R_hor = calc_wire_resistance(circles[0], junctions[0], junctions[2], R);
        matrix_[0][0] = R_ver*R_j + R_hor*R_j + R_hor*R_ver;
        matrix_[0][1] = -R_hor*R_j;
        matrix_[0][2] = -R_hor*R_ver;

        //node right down
        R_hor = calc_wire_resistance(circles[0], junctions[1], junctions[3], R);
        matrix_[1][0] = -R_hor*R_j;
        matrix_[1][1] = R_ver*R_j + R_hor*R_j + R_hor*R_ver;
        matrix_[1][3] = -R_hor*R_ver;
        
        I[0] = R_ver*R_j;
        I[1] = R_ver*R_j;


        //---OTHER CIRCLES---

        for (int i = 1; i < circles.size()-1; i++) {
            junctions = calc_junctions(i, circles, r);
            //node left up
            R_ver = calc_wire_resistance(circles[i], junctions[0], junctions[1], R);
            R_hor = calc_wire_resistance(circles[i], junctions[0], junctions[2], R);
            matrix_[i*4-2][i*4-4] = -R_ver*R_hor;
            matrix_[i*4-2][i*4-2] = R_ver*R_hor + R_j*R_hor + R_j*R_ver;
            matrix_[i*4-2][i*4-1] = -R_j*R_hor;
            matrix_[i*4-2][i*4] -R_j*R_ver;

            //node left down
            R_hor = calc_wire_resistance(circles[i], junctions[1], junctions[3], R);
            matrix_[i*4-1][i*4-3] = -R_ver*R_hor;
            matrix_[i*4-1][i*4-2] = -R_j*R_hor;
            matrix_[i*4-1][i*4-1] = R_ver*R_hor + R_j*R_hor + R_j*R_ver;
            matrix_[i*4-1][i*4+1] = -R_j*R_ver;

            //node right up
            R_ver = calc_wire_resistance(circles[i], junctions[2], junctions[3], R);
            R_hor = calc_wire_resistance(circles[i], junctions[0], junctions[2], R);
            matrix_[i*4][i*4-2] = -R_ver*R_j;
            matrix_[i*4][i*4] = R_ver*R_j + R_hor*R_j + R_hor*R_ver;
            matrix_[i*4][i*4+1] = -R_hor*R_j;
            matrix_[i*4][i*4+2] = -R_hor*R_ver;

            //node right down
            R_hor = calc_wire_resistance(circles[i], junctions[1], junctions[3], R);
            matrix_[i*4+1][i*4-1] = -R_ver*R_j;
            matrix_[i*4+1][i*4] = -R_hor*R_j;
            matrix_[i*4+1][i*4+1] = R_ver*R_j + R_hor*R_j + R_hor*R_ver;
            matrix_[i*4+1][i*4+3] = -R_hor*R_ver;
        }
        //---FINAL CIRCLE---

        junctions = calc_junctions_last(circles, r);
        //node left up
        R_ver = calc_wire_resistance(circles[-1], junctions[0], junctions[1], R);
        R_hor = calc_wire_resistance(circles[-1], junctions[0], junctions[2], R);
        std::cout << R_ver << " " << R_hor << std::endl;
        for (const auto& i: junctions) {
            for (const auto& j: i){
                std::cout << j << " ";
            }
            std::cout << "\n";
        }
        matrix_[nr_rows_-2][nr_rows_-4] = -R_ver*R_hor;
        matrix_[nr_rows_-2][nr_rows_-2] = R_ver*R_hor + R_j*R_hor + R_j*R_ver;
        matrix_[nr_rows_-2][nr_rows_-1] = -R_j*R_hor;
        //node left down
        R_hor = calc_wire_resistance(circles[-1], junctions[1], junctions[3], R);
        matrix_[nr_rows_-1][nr_rows_-3] = -R_ver*R_hor;
        matrix_[nr_rows_-1][nr_rows_-2] = -R_j*R_hor;
        matrix_[nr_rows_-1][nr_rows_-1] = R_ver*R_hor + R_j*R_hor + R_j*R_ver;
    }

};
    


string extract_parameter_from_string(string str) {
    //returns the part of the string after the ':' sign. Used for handle_input
    string parameter {""};
    int start_index = str.find(':') ;
    for(int i {start_index+1}; i < str.size(); i++) {
        parameter += str[i];
    }
    return parameter;
}

tuple<int, double, double, double, double> handle_input(string input_file_name){
    ifstream input_file(input_file_name);
    string line;
    int N;
    double a;
    double r;
    double R;
    double R_j;
    int counter {0};
    while (getline (input_file, line)) {
        cout << line << endl;
        if (counter == 0) {
            N = stoi(extract_parameter_from_string(line));
        } else if (counter == 1){
            a = stof(extract_parameter_from_string(line));
        } else if (counter == 2) {
            r = stof(extract_parameter_from_string(line));
        } else if (counter == 3) {
            R = stof(extract_parameter_from_string(line));
        } else if (counter == 4) {
            R_j = stof(extract_parameter_from_string(line));
        }
        counter++;
    }
    auto parameters = make_tuple(N, a, r, R, R_j);
    input_file.close();
    return parameters;
}


int main(){
    cout << "START OF MAIN" << endl;
    
    //---INITIALISATIONS---

    double R_ver;   //'vertical' resistance
    double R_hor;   //'horizontal' resistance

    auto parameters = handle_input("1D_input.txt");
    int N = get<0>(parameters);
    double a = get<1>(parameters);
    double r = get<2>(parameters);
    double R = get<3>(parameters);
    double R_j = get<4>(parameters);

    valarray<array<double, 3>> circles = generate_circles(N, a, r);
    write_circles_to_file(circles, "../1D/circles.txt");

    Mat G(4*(N-1), 4*(N-1));
    valarray<double> I(4*(N-1));

    cout << "Building matrix...\n";
    
    G.build_matrix(circles, r, R, R_j);

    //---SOLVING THE MATRIX---
    G.print();
    auto V = G.solve_matrix(I);
    cout << endl << "Solutions:\n" << V << endl;

    cout << "END OF MAIN" << endl;
    return 0;
}