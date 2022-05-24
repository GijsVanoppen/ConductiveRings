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





using namespace std;



double dist(std::array<double,2> point_1, std::array<double,2> point_2){
    //returns the distance between two points
    return sqrt(pow(point_1.front()-point_2.front(),2) + pow(point_1.back()-point_2.back(), 2));
}

double dist_circles(std::array<double,3> c_1, std::array<double,3> c_2){
    //returns the distance between two circle centers. the arrays contain {x_center, y_center, r}
    return sqrt(pow(c_1.front()-c_2.front(),2) + pow(c_1[1]-c_2[1], 2));
}

valarray<array<double, 3>> generate_circles(int N, double a, double r, bool wiggle) {
    std::valarray<std::array<double, 3>> circles(N);
    if (wiggle) {
        //rng
        using seed_dist_t = std::uniform_int_distribution<size_t>;
        std::random_device dev;
        seed_dist_t seed_distr(0, std::numeric_limits<size_t>::max());
        auto seed = seed_distr(dev);
        std::mt19937_64 engine(seed);
        auto distr = bind(std::uniform_real_distribution<double>(-1.0, 1.0),engine);
        
        double deviation_coefficient = -a/4 + 0.25*sqrt(8*pow(r,2) - pow(a,2));  //determines the max value the circles may deviate from regular lattice
        
        for (int i {1}; i < N-1; i++) {
            circles[i][0] = a*i + distr()*deviation_coefficient;      //x values of circle center, with a small deviation from regular lattice (distr())
            circles[i][1] = distr()*deviation_coefficient;            //y values of circle center, with a small deviation from regular lattice (distr())
            circles[i][2] = r;                                        //circle radius
        } 
        circles[0][2] = r;              //first circle, no deviations from regular lattice
        circles[N-1][0] = a*(N-1);      //last circle, no deviations from regular lattice
        circles[N-1][2] = r;
    } else {
        for (int i = 0; i < N; i++) {   //no deviations, regular pattern and y = 0
            circles[i][0] = a*i;
            circles[i][2] = r;
        }
    }
    return circles;
}

double calc_wire_resistance(array<double,3> circle, array<double, 2> point_1, array<double, 2> point_2, double R) {
    double r = circle.back();
    return (r*R*acos(1 - 0.5*pow(dist(point_1, point_2)/r,2)));
}

array<array<double, 2>, 4> calc_junctions(int i, valarray<array<double, 3>> circles, double r) {
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

class File_Handler {
    public:
    string extract_parameter_from_string(string str) {
        //returns the part of the string after the ':' sign. Used for handle_input
        string parameter {""};
        int start_index = str.find(':') ;
        for(int i {start_index+1}; i < str.size(); i++) {
            parameter += str[i];
        }
        return parameter;
    }

    tuple<int, double, double, double, double, bool, bool, bool, bool, int> handle_input(string input_file_name){
        ifstream input_file(input_file_name);
        string line;
        int N;
        double a;
        double r;
        double R;
        double R_j;
        bool print_G;
        bool print_V;
        bool wiggle;
        bool stretch;
        int a_iterations;
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
            } else if (counter == 5) {
                if (extract_parameter_from_string(line) == "1"){
                    print_G = true;
                } else {
                    print_G = false;
                }
            } else if (counter == 6) {
                if (extract_parameter_from_string(line) == "1"){
                    print_V = true;
                } else {
                    print_V = false;
                }
            } else if (counter == 7) {
                if (extract_parameter_from_string(line) == "1"){
                    wiggle = true;
                } else {
                    wiggle = false;
                }
            } else if (counter == 8) {
                if (extract_parameter_from_string(line) == "1"){
                    stretch = true;
                } else {
                    stretch = false;
                }
            } else if (counter == 9){
                a_iterations = stoi(extract_parameter_from_string(line));
            }
            counter++;
        }
        auto parameters = make_tuple(N, a, r, R, R_j, print_G, print_V, wiggle, stretch, a_iterations);
        input_file.close();
        return parameters;
    }

    void write_circles_to_file(std::valarray<std::array<double,3>> circles, std::string file_name) {
        std::ofstream circles_file(file_name);
        for (const auto& circle: circles) {
            circles_file << circle.front() << " " << circle[1] << " " << circle.back() << "\n";
        }
        circles_file.close();
    }

    void write_results_to_file(std::string file_name, std::valarray<std::valarray<double>> results_valarray) {
        std::ofstream results_file(file_name);
        for (const auto& result: results_valarray) {
            for (const auto& value: result) {
                results_file << value;
                if (value == result[result.size()-1]) {
                    results_file << "\n";
                } else {
                    results_file << " ";
                }
            }
        }
        results_file.close();
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
    
    double get(const int row, const int col) {
        return matrix_[row][col];
    }

    void set(const int row, const int col, const double value) {
        matrix_[row, col] = value;
    }

    int nr_rows() const { return nr_rows_; }
    int nr_cows() const { return nr_cols_; }
    
    void print(){
        double value;
        for (int i {0}; i < nr_cols_; i++) {
            for (int j {0}; j < nr_rows_; j++){
                value = matrix_[i][j];
                if (value < 0) {
                    cout.precision(4);
                } else {
                    cout.precision(5);
                }
                cout << value << "\t";
            }   
            cout << "\n";
        }
    }



    array<double, 2> build_matrix(std::valarray<std::array<double, 3>> circles, double r, double R, double R_j){
        double R_ver;   //'vertical' resistance
        double R_hor;   //'horizontal' resistance
        double R_1;     //for stretch
        double R_2;     //for stretch
        //---FIRST CIRCLE---

        auto junctions = calc_junctions_first(circles, r);


        //node right up
        R_ver = calc_wire_resistance(circles[0], junctions[2], junctions[3], R);
        R_hor = calc_wire_resistance(circles[0], junctions[0], junctions[2], R);
        matrix_[0][0] = R_ver*R_j + R_hor*R_j + R_hor*R_ver;
        matrix_[0][1] = -R_hor*R_j;
        matrix_[0][2] = -R_hor*R_ver;
        R_1 = R_hor;

        //node right down
        R_hor = calc_wire_resistance(circles[0], junctions[1], junctions[3], R);
        matrix_[1][0] = -R_hor*R_j;
        matrix_[1][1] = R_ver*R_j + R_hor*R_j + R_hor*R_ver;
        matrix_[1][3] = -R_hor*R_ver;
        R_2 = R_hor;
        
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
            matrix_[i*4-2][i*4] = -R_j*R_ver;

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
        R_ver = calc_wire_resistance(circles[circles.size()-1], junctions[0], junctions[1], R);
        R_hor = calc_wire_resistance(circles[circles.size()-1], junctions[0], junctions[2], R);

        matrix_[nr_rows_-2][nr_rows_-4] = -R_ver*R_hor;
        matrix_[nr_rows_-2][nr_rows_-2] = R_ver*R_hor + R_j*R_hor + R_j*R_ver;
        matrix_[nr_rows_-2][nr_rows_-1] = -R_j*R_hor;
        //node left down
        R_hor = calc_wire_resistance(circles[circles.size()-1], junctions[1], junctions[3], R);
        matrix_[nr_rows_-1][nr_rows_-3] = -R_ver*R_hor;
        matrix_[nr_rows_-1][nr_rows_-2] = -R_j*R_hor;
        matrix_[nr_rows_-1][nr_rows_-1] = R_ver*R_hor + R_j*R_hor + R_j*R_ver;

        array<double, 2> resistances_array {R_1, R_2}; //for stretch
        return resistances_array;
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
        
        //---WRITING RESULTS TO .txt FILE---
        valarray<valarray<double>> results(x.size());
        
        for (int index=0; index < x.size(); index++) {
            valarray<double> result (1);
            result[0] = x[index];
            results[index] = result;
        }

        File_Handler file_handler;
        file_handler.write_results_to_file("V.txt", results);
        return x;
    }

    void stretch(const int N, const double R, const double R_j, const int a_iterations, const double r){
        double R_tot;
        double R_1;
        double R_2;
        const double a_delta = r/((double)a_iterations+1);

        valarray<valarray<double>> results(a_iterations);   //will be filled with values of R_tot for different values of a  
        valarray<double> result(2);                         //will become the elements of results


        double a = r;

        for (int iteration = 0; iteration < a_iterations; iteration++) {
            a += a_delta;
            auto circles = generate_circles(N, a, r, false);
            auto resistances_array = build_matrix(circles, r, R, R_j);
            R_1 = resistances_array[0];    //two first resistances returned from build_matrix, used to calc R_tot
            R_2 = resistances_array[1];    //two first resistances returned from build_matrix, used to calc R_tot
            auto V = solve_matrix(I);
    
            //---CALCULATING R_tot AND COLLECTING RESULTS---
            R_tot = 1/((1-V[0])/R_1 + (1-V[1])/R_2);     //ohm's law onto the first 2 resistors, then added to get current. Lastly, used R = 1V/I
            result[0] = a;
            result[1] = R_tot;
            results[iteration] = result;
        }
        //---SENDING RESULTS TO .txt FILE---
        File_Handler file_handler;
        file_handler.write_results_to_file("results_stretch.txt", results);    
    
    }
};


void timing(chrono::time_point<chrono::system_clock> time1, chrono::time_point<chrono::system_clock> time2){
	//prints time between time1 and time 2
	const auto duration_us = chrono::duration_cast<chrono::microseconds>(time2 - time1);
	const int min = duration_us.count()/60000000;
	const int s = duration_us.count()/1000000 - min*60;
	const int ms = duration_us.count()/1000 - min*60000 - s*1000;
	const int um = duration_us.count() - min*60000000 - s*1000000 - ms*1000;
	cout << min << " min, " << s << " s, " << ms << " ms, " << um << " um \n";
}





int main(){
    const auto time_start = chrono::high_resolution_clock::now();
    cout << "START OF MAIN" << endl;
    //---INITIALISATIONS---
    
    File_Handler file_handler;

    double R_ver;   //'vertical' resistance
    double R_hor;   //'horizontal' resistance

    auto parameters = file_handler.handle_input("input.txt");
    int N = get<0>(parameters);
    double a = get<1>(parameters);
    double r = get<2>(parameters);
    double R = get<3>(parameters);
    double R_j = get<4>(parameters);
    bool print_G = get<5>(parameters);
    bool print_V = get<6>(parameters);
    bool wiggle = get<7>(parameters);
    bool stretch = get<8>(parameters);
    int a_iterations = get<9>(parameters);
    
    
    
    
    Mat G(4*(N-1), 4*(N-1));

    if (stretch) {
    
        G.stretch(N, R, R_j, a_iterations, r);
    
    } else {

        valarray<array<double, 3>> circles = generate_circles(N, a, r, wiggle);
        file_handler.write_circles_to_file(circles, "circles.txt");


        //---BUILDING THE MATRIX---
        cout << "Building matrix...\n";
        
        G.build_matrix(circles, r, R, R_j);
        const auto time_build = chrono::high_resolution_clock::now();
        cout << "Time to build matrix: \n";
        timing(time_start, time_build);
        if (print_G){
            G.print();
        }

        //---SOLVING THE MATRIX---
        auto V = G.solve_matrix(G.I);
        
        if (print_V) {
            cout << endl << "Solutions:\n" << V << endl;
        }
        cout << "Time to solve matrix: \n";
        const auto time_end = chrono::high_resolution_clock::now();
        timing(time_build, time_end);
        cout << "END OF MAIN" << endl;
        cout << "Total time: \n";
        timing(time_start, time_end);
    }
    return 0;
}

/*
for (const auto& junction: junctions){
    cout << junction[0] << " " << junction[1] << endl;
}
*/