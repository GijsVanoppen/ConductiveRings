#include "Mat.h"
#include "calc_junctions.h"
#include "calc_junctions_first.h"
#include "calc_junctions_last.h"
#include "calc_wire_resistance.h"



void Mat::print(){
    for (int i {0}; i < nr_cols_; i++) {
        for (int j {0}; j < nr_rows_; j++){
            std::cout << matrix_[i][j] << " ";
        }
        std::cout << "\n";
    }
}

Eigen::VectorXd Mat::solve_matrix(std::valarray<double> b_) {
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

void Mat::build_matrix(std::valarray<std::array<double, 3>> circles, double r, double R, double R_j){
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

Mat::Mat(const int nr_rows, const int nr_cols):
    nr_rows_ {nr_rows},
    nr_cols_ {nr_cols},
    row (nr_cols),
    matrix_ (nr_rows),
    I (nr_rows) {};
    