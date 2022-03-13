#include <iostream>
#include <array>
#include <valarray>
#include <cmath>
#include <ctgmath>
#include <C:\Python\Python38-32\Scripts\Toolbox\eigen\eigen-3.4.0\Eigen\Sparse>
#include <fstream>
#include <string>
#include <random>






using namespace std;





class Mat {
    private:
    int nr_rows_;
    int nr_cols_;
    valarray<valarray<double>> matrix_;
    valarray<double> matrix_row_;

    public:
    std::valarray<double> I;
    
    Mat(const int nr_rows, const int nr_cols){
        nr_rows_ = nr_rows;
        nr_cols_ = nr_cols;
        matrix_row_.resize(nr_cols_,0);
        matrix_.resize(nr_rows_);
        for (int i = 0; i < matrix_.size(); i++){
            matrix_[i] = matrix_row_;
        }    
        I.resize(nr_rows);
    }  
    
    double get(const int row, const int col) {
        return matrix_[row][col];
    }

    void set(const int row, const int col, const double value) {
        matrix_[row, col] = value;
    }

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
};

int main(){
    Mat G(
}