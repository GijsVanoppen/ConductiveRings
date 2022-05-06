#include <iostream>
#include <array>
#include <valarray>
#include <cmath>
#include <ctgmath>
#include <C:\Python\Python38-32\Scripts\Toolbox\eigen\eigen-3.4.0\Eigen\Sparse>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
#include <vector>

using namespace std;

class Ring{
    private:

    double x_;
    double y_;
    double r_;
    double R_;
    int ring_index_;

    public:

    Ring(const double x, const double y, const double r, const double R, const int ring_index) {
        x_ = x;
        y_ = y;
        r_ = r;
        R_ = R;
        ring_index_ = ring_index;
    }

    vector<tuple<double,double,double,int>> junctions;     //vector filled with angle between horizontal and center-junction line, index of overlapping circle, x_junction, y_junction


    int cluster_index = -1;
    
    double get_x() {return x_;}
    double get_y() {return y_;}
    double get_r() {return r_;}
    double get_R() {return R_;} 
    int get_ring_index() {return ring_index_;}
    void set_ring_index(int ring_index) {ring_index_ = ring_index;}
    
    double distance_between_rings(Ring ring2) {
        return sqrt(pow(get_x()-ring2.get_x(),2) + pow(get_y()-ring2.get_y(),2));
    }

    bool check_overlap(Ring ring2) {
        //check if the two Ring objects are not the same
        if (get_ring_index() == ring2.get_ring_index()) {
            return false;
        } else {
            if (distance_between_rings(ring2) < (get_r() + ring2.get_r())) {
                return true;
            } else {
                return false;
            }
        }
    }

    // int update_cluster_index(vector<vector<int>> clusters) {   //returns index of clusters in which vector the ring index is located. Returns -1 if this ring is not in a cluster
    //     for (int cluster_index_ = 0; cluster_index_ < clusters.size(); cluster_index_++) {
    //         for (const auto& value: clusters[cluster_index_]) {
    //             if (value == get_ring_index()) {
    //                 return cluster_index_;
    //             }
    //         }
    //     }
    //     return -1;
    // }

    array<double,4> calc_junctions(Ring ring) {
        const double x1 = get_x();            //x of center of this circle   
        const double y1 = get_y();            //y of center of this circle
        const double r1 = get_r();            //radius of this circle
        const double x2 = ring.get_x();       //x of center of other circle
        const double y2 = ring.get_y();       //y of center of other circle
        const double r2 = ring.get_r();       //radius of other circle
        const double d = distance_between_rings(ring);
        const double const1 = 0.5*(pow(r1,2) - pow(r2,2))/pow(d,2);
        const double const2 = sqrt(0.5*(pow(r1,2) + pow(r2,2))/pow(d,2) - pow(const1,2) - 0.25);

        const double junction_1_x = 0.5*(x1+x2) + const1*(x2-x1) + const2*(y2-y1);    
        const double junction_2_x = 0.5*(x1+x2) + const1*(x2-x1) - const2*(y2-y1);    
        const double junction_1_y = 0.5*(y1+y2) + const1*(y2-y1) + const2*(x1-x2);
        const double junction_2_y = 0.5*(y1+y2) + const1*(y2-y1) - const2*(x1-x2);
        array<double,4> junctions_ = {junction_1_x, junction_1_y, junction_2_x, junction_2_y};
        return junctions_;
    }

    double calc_junction_angle(double x_junction, double y_junction) {
        if (x_junction > get_x()) { 
            if (y_junction > get_y()) { //first quadrant
                cout << "kwadrant 1\n";
                return acos((x_junction - get_x())/get_r());
            } else {                    //fourth quadrant
                cout << "kwadrant 4\n";
                return (2*M_PI - acos((x_junction - get_x())/get_r()));
            }
        } else {
            if (y_junction < get_y()) { //third quadrant
                cout << "kwadrant 3\n";
                return (M_PI + acos((get_x() - x_junction)/get_r()));
            } else {                    //second quadrant
                cout << "kwadrant 2\n";
                return (M_PI - acos((get_x() - x_junction)/get_r()));
            }
        }
    }
    vector<Ring> return_neighbours(vector<Ring> rings) {     //returns series of Rings (or better, their indices) that overlap with this Ring 
        vector<Ring> neighbours;
        for (auto& ring: rings) {
            if (check_overlap(ring)) {
                neighbours.push_back(ring);
            }
        }
        return neighbours;
    }

    bool sort_junctions(tuple<double,double,double,int> junction_info1, tuple<double,double,double,int> junction_info2){ //used to sort
        return (get<2>(junction_info1) < get<2>(junction_info2));
    }

    void calc_all_junctions(vector<Ring> rings, double box_width) {
        vector<Ring> neighbours = return_neighbours(rings);
        for (auto& ring2: neighbours) {
            //get info about the pair of junctions            
            array<double, 4> new_junctions = calc_junctions(ring2);
            double x_junction_1 = new_junctions[0];
            double y_junction_1 = new_junctions[1];               
            double x_junction_2 = new_junctions[2];               
            double y_junction_2 = new_junctions[3]; 
            // add the info of the pair of junctions to this circle    
            tuple<double,double,double,int> junction_info;          //info: x,y,angle,ring_index
            if ((x_junction_1 > 0) && (x_junction_1 < box_width)) {
                junction_info = {x_junction_1, y_junction_1, calc_junction_angle(x_junction_1, y_junction_1), ring2.get_ring_index()};
                junctions.push_back(junction_info);
            }
            if ((x_junction_2 > 0) && (x_junction_2 < box_width)) {
                junction_info = {x_junction_2, y_junction_2, calc_junction_angle(x_junction_2, y_junction_2), ring2.get_ring_index()};
                junctions.push_back(junction_info);
            }
            // add the info of the pair of junctions to the neighbour circle           
            // junction_info = {x_junction_1, y_junction_1, ring2.calc_junction_angle(x_junction_1, y_junction_1), get_ring_index()};
            // ring2.junctions.push_back(junction_info);
            // junction_info = {x_junction_2, y_junction_2, ring2.calc_junction_angle(x_junction_2, y_junction_2), get_ring_index()};
            // ring2.junctions.push_back(junction_info);
        }
        //now sorting the junction based on increasing angle
        // sort(junctions.begin(), junctions.end(), sort_junctions);
        // cout << "Sorted junctions: \n";
        // for (auto x: junctions) {
        //     cout << get<2>(x) <<" ";
        // }


    }

    double calc_resistance(double angle1, double angle2) {
        return (abs(angle1-angle2)*get_r()*get_R());
    }

};






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
    
    double get_element(const int row, const int col) {
        return matrix_[row][col];
    }

    void set_element(const int row, const int col, const double value) {
        matrix_[row][col] = value;
    }

    int nr_rows() const { return nr_rows_; }
    int nr_cols() const { return nr_cols_; }
    
    void print(){
        double value;
        for (int i {0}; i < nr_rows(); i++) {
            for (int j {0}; j < nr_cols(); j++){
                value = matrix_[i][j];
                /*if (value < 0) {
                    cout.precision(4);
                } else {
                    cout.precision(5);
                }
                cout << value << "\t";
                */
                cout << value << " ";
            }   
            cout << "\n";
        }
    }




};

int main() {
    for (int i= 0; i< 0; i++) {
        cout << "y";
    }
}