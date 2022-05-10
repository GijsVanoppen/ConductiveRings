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


bool is_in_bounds(const double pt_x, const double box_width){
    return ((pt_x >= 0) && (pt_x <= box_width));
}

void timing(chrono::time_point<chrono::system_clock> time1, chrono::time_point<chrono::system_clock> time2){
	//prints time between time1 and time 2
	const auto duration_us = chrono::duration_cast<chrono::microseconds>(time2 - time1);
	const int min = duration_us.count()/60000000;
	const int s = duration_us.count()/1000000 - min*60;
	const int ms = duration_us.count()/1000 - min*60000 - s*1000;
	const int um = duration_us.count() - min*60000000 - s*1000000 - ms*1000;
	cout << min << " min, " << s << " s, " << ms << " ms, " << um << " um \n";
}

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
        //first, check if the two Ring objects are not the same
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
                return acos((x_junction - get_x())/get_r());
            } else {                    //fourth quadrant
                return (2*M_PI - acos((x_junction - get_x())/get_r()));
            }
        } else {
            if (y_junction < get_y()) { //third quadrant
                return (M_PI + acos((get_x() - x_junction)/get_r()));
            } else {                    //second quadrant
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
        }
        
    }

    double calc_resistance(double angle1, double angle2) {  //calculates resistance of wire from angle1 to angle2 (counter clockwise)
        double a = angle2-angle1;
        if (a < 0) {
            a += M_PI*2;
        }
        return (a*get_r()*get_R());
    }

    int calc_junction_index_before_crit_angle(double crit_angle) {
        if (junctions.size() == 1) {
            return 0;
        } else {
            for (int junction_index = 0; junction_index < junctions.size(); junction_index++) {
                if (get<2>(junctions[junction_index]) > crit_angle) {
                    return junction_index-1;
                } 
            }
            return junctions.size()-1;     
        }
    }

    array<double, 2> calc_angle_of_intersection_with_bound(double x_bound) {
        array<double,2> a;
        //first, we need to calculate the difference in y between the intersection and the ring center
        double delta_y = sqrt(pow(get_r(),2) - abs(get_x() - x_bound));
        a[0] = calc_junction_angle(x_bound, get_y() + delta_y);     //first junction (higher y)
        a[1] = calc_junction_angle(x_bound, get_y() - delta_y);     //second junction (lower y)
        return a;
    }

    int return_junction_index_from_other_ring(vector<Ring> rings, tuple<double, double, double, int> junction) {
        //first, get the neighbouring ring
        int wanted_ring_index = get<3>(junction);
        
        int index = 0;
        for (auto& ring: rings) {
            if (ring.get_ring_index() == wanted_ring_index) {
                break;
            }
            index++;
        }
        Ring ring2 = rings[index];


        for (int junction_index2 = 0; junction_index2 < ring2.junctions.size(); junction_index2++) {   //loop over its junctions
            auto junction2 = ring2.junctions[junction_index2];
            if ((get<0>(junction2) == get<0>(junction)) && (get<1>(junction2) == get<1>(junction))) {   //x and y are the same, we found the correct junction
                return junction_index2;
            }
        }
        cout << "Error: didnt find the correct junction index!!!!\n";
        return -1;
    }

};



bool sort_junctions(tuple<double,double,double,int> junction_info1, tuple<double,double,double,int> junction_info2){ //used to sort junctions by angle
    return (get<2>(junction_info1) < get<2>(junction_info2));
}

bool sort_rings(Ring ring1, Ring ring2){ //used to sort rings by x_value of center
    return (ring1.get_x() < ring2.get_x());
}

int node_nr_at_junction(vector<Ring> rings, int ring_index, int junction_index){
    int node_nr = 0;
    for (int ring_index2 = 0; rings[ring_index2].get_ring_index() != ring_index; ring_index2++){
        Ring ring = rings[ring_index2];
        node_nr += ring.junctions.size();
    }
    return (node_nr + junction_index);
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

    auto handle_input(string input_file_name){
        //open input file
        ifstream input_file(input_file_name);
        string line;
        int N;
        double box_width_min;
        double box_width_max;
        double box_width_iterations;
        double box_length_min;
        double box_length_max;
        double box_length_iterations;
        double V_diff;
        double r_min;
        double r_max;
        double R_min;
        double R_max;
        double R_j;
        bool print_G;
        bool print_V;
        bool import_rings;
        string rings_file_name;

        //read input parameters
        int counter {0};
        while (getline (input_file, line)) {
            cout << line << endl;
            if (counter == 0) {
                N = stoi(extract_parameter_from_string(line));
            } else if (counter == 1){
                box_width_min = stof(extract_parameter_from_string(line));
            } else if (counter == 2){
                box_width_max = stof(extract_parameter_from_string(line));
            } else if (counter == 3){
                box_width_iterations = stof(extract_parameter_from_string(line));
            } else if (counter == 4) {
                box_length_min = stof(extract_parameter_from_string(line));
            } else if (counter == 5) {
                box_length_max = stof(extract_parameter_from_string(line));
            } else if (counter == 6) {
                box_length_iterations = stof(extract_parameter_from_string(line));
            } else if (counter == 7) {
                V_diff = stof(extract_parameter_from_string(line));
            } else if (counter == 8) {
                r_min = stof(extract_parameter_from_string(line));
            } else if (counter == 9) {
                r_max = stof(extract_parameter_from_string(line));
            } else if (counter == 10) {
                R_min = stof(extract_parameter_from_string(line));
            } else if (counter == 11) {
                R_max = stof(extract_parameter_from_string(line));
            } else if (counter == 12) {
                R_j = stof(extract_parameter_from_string(line));
            } else if (counter == 13) {
                if (extract_parameter_from_string(line) == "1"){
                    print_G = true;
                } else {
                    print_G = false;
                }
            } else if (counter == 14) {
                if (extract_parameter_from_string(line) == "1"){
                    print_V = true;
                } else {
                    print_V = false;
                }
            } else if (counter == 15) {
                if (extract_parameter_from_string(line) == "1"){
                    import_rings = true;
                } else {
                    import_rings = false;
                }
            } else if (counter == 16) {
                rings_file_name = extract_parameter_from_string(line);
            }
            counter++;
        }
        //combine input parameters in tuple and close file
        auto parameters = make_tuple(N, box_width_min, box_width_max, box_width_iterations, box_length_min, box_length_max, box_length_iterations, V_diff, r_min, r_max, R_min, R_max, R_j, print_G, print_V, import_rings, rings_file_name);
        input_file.close();

        
        return parameters;
    }

    void write_rings_to_file(vector<Ring> rings, string file_name) {
        ofstream rings_file(file_name);
        for (auto& ring: rings) {
            rings_file << ring.get_ring_index() << " "  <<ring.get_x() << " " << ring.get_y() << " " << ring.get_r() << " " << ring.get_R() << "\n";
        }
        rings_file.close();
    }

    void write_junctions_to_file(vector<Ring> rings, string file_name) {
        ofstream junctions_file(file_name);
        for(auto ring: rings) {
            for (auto junction: ring.junctions) {
                junctions_file << get<0>(junction) << " " << get<1>(junction) << endl; 
            }
        }
        junctions_file.close();
    }

    void write_results_to_file(Eigen::VectorXd V, string file_name) {
        ofstream results_file(file_name);
        for (auto element: V) {
            results_file << element << endl;
        }
        results_file.close();
    }

    int nth_occurence(const string str, const char character , const int n) {
        int answer;
        int counter = 0;
        for (int i = 0; i < str.size(); i++) {
            if (str[i] == character) {
                counter++;
                if (counter == n) {
                    answer = i;
                    break;
                }
            }
        }
        return answer;
    } 

    vector<Ring> read_rings_from_file(string file_name) {
        string line;
        ifstream rings_file(file_name);
        double x;
        double y;
        double r;
        double R;
        int ring_index;
        vector<Ring> rings;
        while (getline (rings_file, line)) {
            int space_index_1 = nth_occurence(line, ' ', 1);            
            int space_index_2 = nth_occurence(line, ' ', 2);            
            int space_index_3 = nth_occurence(line, ' ', 3);                        
            int space_index_4 = nth_occurence(line, ' ', 4);

            int x_len = space_index_2-space_index_1;
            int y_len = space_index_3-space_index_2;
            int r_len = space_index_4-space_index_3;
            int R_len = line.size() - space_index_4;

            ring_index = stoi(line.substr(0, space_index_1));
            x = stod(line.substr(space_index_1+1, x_len-1));
            y = stod(line.substr(space_index_2+1, y_len-1));
            r = stod(line.substr(space_index_3+1, r_len-1));
            R = stod(line.substr(space_index_4+1, R_len-1));
            Ring ring(x,y,r,R,ring_index);
            rings.push_back(ring);

        
        }
        return rings;
    }

    void write_resistance_to_file(const string file_name, const double total_resistance, const int width_iteration_count, const int length_iteration_count) {
        
        if ((width_iteration_count==0) && (length_iteration_count==0)) {//first iteration, so delete previous content of the file
           ofstream resistance_file;
           resistance_file.open(file_name, ofstream::trunc);
           resistance_file.close();
        }
       
        ofstream resistance_file;
        resistance_file.open(file_name, ofstream::app);
        resistance_file << total_resistance << endl;
        resistance_file.close();
    }

};

vector<Ring> generate_rings(const int N, const double box_width, const double box_length, const double r_min, const double r_max, const double R_min, const double R_max){
    //random number generator
    using seed_dist_t = std::uniform_int_distribution<size_t>;
    seed_dist_t seed_distr(0, std::numeric_limits<size_t>::max());
    random_device dev;
    auto seed = seed_distr(dev);
    mt19937_64 engine(seed);
    auto distr = bind(std::uniform_real_distribution<double>(0.0, 1.0),engine);
    
    vector<Ring> rings;         //create vector, will be filled with Ring objects
    for (int i=0; i < N; i++) {
        //choose variables
        double r = r_min + distr()*(r_max-r_min);   //r between [r_min, r_max[
        double x = -r + distr()*(box_width+2*r);      //x between [-r, box_witdh+r[
        double y = r + distr()*(box_length-2*r);      //y between [r, box_length-r[
        double R = R_min + distr()*(R_max-R_min);     //R between [R_min, Rmax[
        Ring ring(x, y, r, R, i);          //create Ring object
        rings.push_back(ring);          //filling of vector
    }  
    return rings;
}


class Mat {
    private:
    int nr_rows_;
    int nr_cols_;
    valarray<valarray<double>> matrix_;
    valarray<double> matrix_row_;

    public:
    valarray<double> I;
    
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
                if (value < 0) {
                    cout.precision(4);
                } else {
                    cout.precision(5);
                }
                cout << value << "\t";
                
                // cout << value << " ";
            }   
            cout << "\n";
        }
    }

    void build_matrix(vector<Ring> rings, double R_j, double V_diff, double box_width, vector<tuple<int, double>>&  current_contributions){
        auto time_start = chrono::high_resolution_clock::now();
        cout << "Building matrix...\n";
        double R_p;
        double R_n;
        int node_counter = 0;
        for (int index = 0; index < rings.size(); index++) {    //loop over all rings
        
            auto& ring = rings[index];
            // cout << "ring " << ring.get_ring_index() << endl;
            for (int junction_index = 0; junction_index < ring.junctions.size(); junction_index++) {    //loop over all junctions of the ring
                
                if ((ring.get_x() <= ring.get_r()) &&  (junction_index == ring.calc_junction_index_before_crit_angle(M_PI))) {                                             //hits left boundary with first junction 
                    // cout << "Case1\n";

                    R_n = ring.calc_resistance(get<2>(ring.junctions[junction_index]), ring.calc_angle_of_intersection_with_bound(0)[0]);
                    current_contributions.push_back(make_tuple(node_counter, R_n));
                    if (junction_index == 0) {
                        if (ring.junctions.size() == 1) {   //only one junction at ring
                            R_p = ring.calc_resistance(ring.calc_angle_of_intersection_with_bound(0)[1], get<2>(ring.junctions[junction_index])); 
                            I[node_counter] += R_n*R_j*V_diff ;
                            current_contributions.push_back(make_tuple(node_counter, R_p));
                        } else {
                            R_p = ring.calc_resistance(get<2>(ring.junctions.back()), get<2>(ring.junctions[junction_index]));   
                            set_element(node_counter, node_counter+ring.junctions.size()-1, -R_j*R_n);                                                                        //previous element
                        }
                    } else {
                        R_p = ring.calc_resistance(get<2>(ring.junctions[junction_index-1]), get<2>(ring.junctions[junction_index]));   
                        set_element(node_counter, node_counter-1, -R_j*R_n);                                                                        //previous element
                    }
                    
                    set_element(node_counter, node_counter, R_n*R_p + R_p*R_j + R_j*R_n);                                                       //diagonal element
                    set_element(node_counter, node_nr_at_junction(rings, get<3>(ring.junctions[junction_index]), ring.return_junction_index_from_other_ring(rings, ring.junctions[junction_index])), -R_n*R_p);      //element corresponding to junction
                    I[node_counter] += R_p*R_j*V_diff;

                } else if ((ring.get_x() <= ring.get_r())  && (junction_index == ((ring.calc_junction_index_before_crit_angle(M_PI)+1) % ring.junctions.size()))) {        //hits left boundary with second junction 
                    // cout << "Case2\n";

                    //side note: here we dont have to make seperate cases for when there is only one junction present, as the previous if-statement would handle this case
                    R_p = ring.calc_resistance(ring.calc_angle_of_intersection_with_bound(0)[1], get<2>(ring.junctions[junction_index]));
                    R_n = ring.calc_resistance(get<2>(ring.junctions[junction_index]), get<2>(ring.junctions[(junction_index+1) % ring.junctions.size()]));     
                    current_contributions.push_back(make_tuple(node_counter, R_p));

                    if (junction_index == (ring.junctions.size()-1)){ //if this is the last junction, the next element should be put at the colomn of the first junction of this ring
                        R_n = ring.calc_resistance(get<2>(ring.junctions[junction_index]), get<2>(ring.junctions.front()));
                        set_element(node_counter, node_counter+1-ring.junctions.size(), -R_j*R_p);        //next element
                    } else {
                        R_n = ring.calc_resistance(get<2>(ring.junctions[junction_index]), get<2>(ring.junctions[junction_index+1]));
                        set_element(node_counter, node_counter+1, -R_j*R_p);        //next element
                    }
               
                    set_element(node_counter, node_counter, R_n*R_p + R_p*R_j + R_j*R_n);                                                       //diagonal element
                    set_element(node_counter, node_nr_at_junction(rings, get<3>(ring.junctions[junction_index]), ring.return_junction_index_from_other_ring(rings, ring.junctions[junction_index])), -R_n*R_p);      //element corresponding to junction
                    I[node_counter] = R_n*R_j*V_diff;
                    

                } else if ((abs(box_width - ring.get_x()) < ring.get_r()) && (junction_index == 0)) {                               //hits right boundary with first junction
                    // cout << "Case3\n";

                    R_p = ring.calc_resistance(ring.calc_angle_of_intersection_with_bound(box_width)[0], get<2>(ring.junctions[junction_index]));
                    if (ring.junctions.size() == 1) {
                        R_n = ring.calc_resistance(get<2>(ring.junctions[junction_index]), ring.calc_angle_of_intersection_with_bound(box_width)[1]); 
                    } else {
                        R_n = ring.calc_resistance(get<2>(ring.junctions[junction_index]), get<2>(ring.junctions[junction_index+1])); 
                        set_element(node_counter, node_counter+1, -R_j*R_p);                                                                        //next element  
                    }
                    set_element(node_counter, node_counter, R_n*R_p + R_p*R_j + R_j*R_n);                                                       //diagonal element
                    set_element(node_counter, node_nr_at_junction(rings, get<3>(ring.junctions[junction_index]), ring.return_junction_index_from_other_ring(rings, ring.junctions[junction_index])), -R_n*R_p);      //element corresponding to junction

                } else if ((abs(box_width - ring.get_x()) < ring.get_r()) && (junction_index == (ring.junctions.size()-1))) {       //hits right boundary with second junction
                    // cout << "Case4\n";

                    //side note: here we dont have to make seperate cases for when there is only one junction present, as the previous else if-statement would handle this case
                    R_p = ring.calc_resistance(get<2>(ring.junctions[junction_index-1]), get<2>(ring.junctions[junction_index]));   
                    R_n = ring.calc_resistance(get<2>(ring.junctions[junction_index]), ring.calc_angle_of_intersection_with_bound(box_width)[1]);

                    set_element(node_counter, node_counter, R_n*R_p + R_p*R_j + R_j*R_n);                                                       //diagonal element
                    set_element(node_counter, node_counter-1, -R_j*R_n);                                                                        //previous element
                    set_element(node_counter, node_nr_at_junction(rings, get<3>(ring.junctions[junction_index]), ring.return_junction_index_from_other_ring(rings, ring.junctions[junction_index])), -R_n*R_p);      //element corresponding to junction
                
                
                } else {
                    // cout << "Case5\n";

                    if (junction_index == 0) { //first junction
                        R_p = ring.calc_resistance(get<2>(ring.junctions.back()), get<2>(ring.junctions[junction_index]));  
                        R_n = ring.calc_resistance(get<2>(ring.junctions[junction_index]), get<2>(ring.junctions[junction_index+1]));

                        set_element(node_counter, node_counter-1+ring.junctions.size(), -R_j*R_n);                                                                        //previous element
                        set_element(node_counter, node_counter+1, -R_p*R_j);                                                                        //next element

                    } else if (junction_index == (ring.junctions.size()-1)){ //last junction
                        R_p = ring.calc_resistance(get<2>(ring.junctions[junction_index-1]), get<2>(ring.junctions[junction_index]));
                        R_n = ring.calc_resistance(get<2>(ring.junctions[junction_index]), get<2>(ring.junctions.front())); 
                        
                        set_element(node_counter, node_counter+1-ring.junctions.size(), -R_p*R_j);                                                                        //next element
                        set_element(node_counter, node_counter-1, -R_j*R_n);                                                                                              //previous element
                    } else {
                        R_p = ring.calc_resistance(get<2>(ring.junctions[junction_index-1]), get<2>(ring.junctions[junction_index]));   
                        R_n = ring.calc_resistance(get<2>(ring.junctions[junction_index]), get<2>(ring.junctions[junction_index+1]));    
                        set_element(node_counter, node_counter-1, -R_j*R_n);                                                                        //previous element
                        set_element(node_counter, node_counter+1, -R_p*R_j);                                                                        //next element
                    }
                    // cout << R_p << " " << R_n << " " << R_j << " " <<  R_n*R_p + R_p*R_j + R_j*R_n <<endl;

                    set_element(node_counter, node_counter, R_n*R_p + R_p*R_j + R_j*R_n);                                                       //diagonal element
                    set_element(node_counter, node_nr_at_junction(rings, get<3>(ring.junctions[junction_index]), ring.return_junction_index_from_other_ring(rings, ring.junctions[junction_index])), -R_n*R_p);      //element corresponding to junction
                }
                node_counter++;
            }
        }
        auto time_end = chrono::high_resolution_clock::now();
        cout << "Succesfully built matrix in ";
        timing(time_start, time_end);
    }
   
    Eigen::VectorXd solve_matrix(valarray<double> b_, bool print_V){
        auto time_start = chrono::high_resolution_clock::now();
        cout << "Solving matrix...\n";
        //first, insert data in special matrix and vectors from Eigen
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
        //solving the matrix using Eigen library
        A.makeCompressed();
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;   
        solver.analyzePattern(A);
        solver.factorize(A);
        x = solver.solve(b);
        //printing solutions if requested
        if (print_V) {
            cout << "Solutions:\n" << x << endl;
        }
        auto time_end = chrono::high_resolution_clock::now();
        cout << "Succesfully solved matrix in ";
        timing(time_start, time_end);
        return x;
    }
};

double calc_total_resistance(vector<tuple<int, double>> current_contributions, Eigen::VectorXd V, double V_diff) {
    double total_current = 0;
    for (auto element: current_contributions) {
        total_current += (V_diff - V[get<0>(element)])/get<1>(element);  //add each contribution to the total current (Ohm's law I = V/R)
    }
    return (V_diff/total_current); //using Ohm's law again, this gives the total resistance of the circuit
}


int main() {
    cout << "START OF MAIN \n"; 
    //---INITIALISATION---
    File_Handler file_handler;
    cout << "Input parameters:\n";
    auto input_pars = file_handler.handle_input("input.txt");

    const int N = get<0>(input_pars);
    const double box_width_min = get<1>(input_pars);
    const double box_width_max = get<2>(input_pars);
    const double box_width_iterations = get<3>(input_pars);
    const double box_length_min = get<4>(input_pars);
    const double box_length_max = get<5>(input_pars);
    const double box_length_iterations = get<6>(input_pars);
    const double V_diff = get<7>(input_pars);
    const double r_min = get<8>(input_pars);
    const double r_max = get<9>(input_pars);
    const double R_min = get<10>(input_pars);
    const double R_max = get<11>(input_pars);
    const double R_j = get<12>(input_pars);
    const bool print_G = get<13>(input_pars);
    const bool print_V = get<14>(input_pars);
    const bool import_rings = get<15>(input_pars);
    const string rings_file_name = get<16>(input_pars);

    cout << endl << endl;



    // ========================
    // ===START MAIN PROGRAM===
    // ========================
    double box_width;
    double box_length;
    
    for (int width_iteration_count = 0; width_iteration_count < box_width_iterations; width_iteration_count++) {    //iterate over different widths of box
        if (box_width_iterations==1) {
            box_width = box_width_max;
        } else {
            box_width = box_width_min + (box_width_max-box_width_min)/(box_width_iterations-1)*width_iteration_count;
        }
        
        for (int length_iteration_count = 0; length_iteration_count < box_length_iterations; length_iteration_count++) {    //iterate over different lengths of box
            if (box_length_iterations==1) {
                box_length = box_length_max;
            } else {
                box_length = box_length_min + (box_length_max-box_length_min)/(box_length_iterations-1)*length_iteration_count;
            }

            vector<tuple<int,double>> current_contributions; //to keep track of current going through circuit

            bool running = true;
            while (running) {   //this while loop enables the program to run again if no network is formed that spans the box
                
                //---GENERATE OR IMPORT RINGS---
                vector<Ring> rings;
                if (import_rings) {
                    rings = file_handler.read_rings_from_file(rings_file_name);
                } else {
                    rings = generate_rings(N, box_width, box_length, r_min, r_max, R_min, R_max);
                }
                file_handler.write_rings_to_file(rings, "rings_all.txt");






                //---FIND CLUSTERS---
                vector<int> rings_to_check;
                int amount_of_clusters = -1;

                for (auto& ring: rings) {
                    if (ring.cluster_index == -1) {
                        rings_to_check.push_back(ring.get_ring_index());
                        amount_of_clusters++;
                        ring.cluster_index = amount_of_clusters;
                    }

                    while (rings_to_check.size()!=0) {
                        int ring_index = rings_to_check.back();
                        Ring& ring = rings[ring_index];
                        rings_to_check.pop_back();

                        for (auto& ring2: rings) {
                            
                            if ((ring2.cluster_index == -1) && ring.check_overlap(ring2)) {
                                auto junctions = ring.calc_junctions(ring2);
                                if (is_in_bounds(junctions[0], box_width) || is_in_bounds(junctions[2], box_width)) {
                                    rings_to_check.push_back(ring2.get_ring_index());
                                    ring2.cluster_index = amount_of_clusters;
                                }
                            }
                        }
                    }
                }

                amount_of_clusters++;
           

                vector<vector<int>> clusters; //resize clusters
                clusters.resize(amount_of_clusters);
                for (auto ring: rings) {    //fill clusters
                    clusters[ring.cluster_index].push_back(ring.get_ring_index());
                }

            
                //---FILTER OUT CLUSTERS THAT DON'T HIT EDGES---
                vector<vector<int>> clusters_;
                int cluster_index_main;
                for (auto cluster: clusters) {
                    bool hit_left_edge = false;
                    bool hit_right_edge = false;
                    for (auto index: cluster) {
                        if (abs(rings[index].get_x()) < rings[index].get_r()) {
                            hit_left_edge = true;
                        }
                        if (abs(rings[index].get_x()-box_width) < rings[index].get_r()) {
                            hit_right_edge = true;
                        }
                    }
                    if (hit_left_edge && hit_right_edge) {
                        clusters_.push_back(cluster);
                        cluster_index_main = rings[cluster.front()].cluster_index; //get the cluster index of this cluster
                        break;
                    }
                    
                }
                
                vector<Ring> rings_;
                for (auto ring: rings) {
                    if (ring.cluster_index == cluster_index_main) {
                        rings_.push_back(ring);
                    }
                }
                rings = rings_;

                if (clusters_.size() == 0) {
                    cout << "No clusters span the box, trying again...\n";
                    continue;
                }
                cout << clusters_.size() << " cluster(s) span(s) the box.\n";
                
                clusters = clusters_;
                
                


                //---GIVE RINGS NEW INDICES AND FILTER RINGS---
                sort(rings.begin(), rings.end(), sort_rings);
                int index = 0;
                for (auto& ring: rings) {
                    if (ring.cluster_index == cluster_index_main) {
                        ring.set_ring_index(index);
                        index++;  
                    } 
                }



                //---MAKE SURE THE FIRST AND THE LAST RING ARE AT THE BOUNDARY---
                bool first_ring_found = false;
                bool last_ring_found = false; 
                for (auto ring: rings) {
                    if (not first_ring_found) {
                        if (abs(ring.get_x()) <= ring.get_r()) {
                            first_ring_found = true;
                            rings[ring.get_ring_index()] = rings.front();
                            rings.front() = ring;
                        }
                    }
                    if (not last_ring_found) {
                        if (abs(ring.get_x() - box_width) <= ring.get_r()) {
                            last_ring_found = true;
                            rings[ring.get_ring_index()] = rings.back();
                            rings.back() = ring;
                        }
                    }
                    if (first_ring_found && last_ring_found) {
                        break;
                    }
                }


                
                file_handler.write_rings_to_file(rings, "rings_main.txt");
            
                //---CALCULATE JUNCTIONS---
                int amount_of_nodes = 0;
                for (auto& ring: rings) {
                    ring.calc_all_junctions(rings, box_width);
                    //now sorting the junction based on increasing angle
                    sort(ring.junctions.begin(), ring.junctions.end(), sort_junctions);
                    amount_of_nodes += ring.junctions.size();
                }    
                file_handler.write_junctions_to_file(rings, "junctions.txt");

                //---BUILDING MATRIX---
                cout << "matrix size: " << amount_of_nodes << endl;
                Mat G(amount_of_nodes,amount_of_nodes);
                G.build_matrix(rings, R_j, V_diff, box_width, current_contributions);
                if (print_G) {
                    G.print();
                }


                
                //---SOLVING THE MATRIX---
                Eigen::VectorXd V = G.solve_matrix(G.I, print_V);
                file_handler.write_results_to_file(V, "results.txt");
                
                const double total_resistance = calc_total_resistance(current_contributions, V, V_diff);
                file_handler.write_resistance_to_file("resistances.txt", total_resistance, width_iteration_count, length_iteration_count);


                
                running = false;
            }
        }
    }


    cout << "END OF MAIN";
    return 0;
}








    // cout << "Printing junctions:\n"; 
    // for (auto ring: rings) {
    //     cout << "Junctions of ring " <<ring.get_ring_index() << endl;
    //     for (auto junction: ring.junctions) {
    //         cout << get<0>(junction) << " " << get<1>(junction) << " " << get<2>(junction) <<" " << get<3>(junction) << endl; 
    //     }
    //     cout << endl;
    // }
