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

    double calc_junction_angle(const double x_junction, const double y_junction) {
        double angle = atan((y_junction-get_y())/(x_junction-get_x()));
        if (x_junction < get_x()) {
            angle += M_PI;
        } else if (y_junction < get_y()){
            angle += 2*M_PI;
        }
        return angle;
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

    void calc_all_junctions(vector<Ring> rings) {
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
            junction_info = {x_junction_1, y_junction_1, calc_junction_angle(x_junction_1, y_junction_1), ring2.get_ring_index()};
            junctions.push_back(junction_info);
            junction_info = {x_junction_2, y_junction_2, calc_junction_angle(x_junction_2, y_junction_2), ring2.get_ring_index()};
            junctions.push_back(junction_info);
            // add the info of the pair of junctions to the neighbour circle           
            // junction_info = {x_junction_1, y_junction_1, ring2.calc_junction_angle(x_junction_1, y_junction_1), get_ring_index()};
            // ring2.junctions.push_back(junction_info);
            // junction_info = {x_junction_2, y_junction_2, ring2.calc_junction_angle(x_junction_2, y_junction_2), get_ring_index()};
            // ring2.junctions.push_back(junction_info);
        }
    }
};

vector<Ring> generate_rings(const int N, const double box_width, const double box_length, const double r_min, const double r_max, const double R){
    //rng
    using seed_dist_t = std::uniform_int_distribution<size_t>;
    seed_dist_t seed_distr(0, std::numeric_limits<size_t>::max());
    random_device dev;
    auto seed = seed_distr(dev);
    mt19937_64 engine(seed);
    auto distr = bind(std::uniform_real_distribution<double>(0.0, 1.0),engine);
    
    vector<Ring> rings;         //create vector, will be filled with Ring objects
    for (int i=0; i < N; i++) {
        //choose variables
        double r = r_min + distr()*(r_max-r_min);   //r between [r_min, r_max]
        double x = r + distr()*(box_width-2*r);      //x between [r, box_witdh-r]
        double y = r + distr()*(box_length-2*r);      //y between [r, box_length-r]
        Ring ring(x, y, r, R, i);          //create Ring object
        rings.push_back(ring);          //filling of vector
    }   
    return rings;

}



void write_rings_to_file(vector<Ring> rings, std::string file_name) {
    std::ofstream rings_file(file_name);
    for (auto& ring: rings) {
        rings_file << ring.get_x() << " " << ring.get_y() << " " << ring.get_r() << " " << ring.get_R() << "\n";
    }
    rings_file.close();
}



int main() {
    cout << "\n\n\nSTART OF MAIN \n"; 
    int N = 10;
    //---GENERATE RINGS---
    vector<Ring> rings = generate_rings(N, 10, 10, 1, 2, 1); 
    write_rings_to_file(rings, "rings.txt");


    //---FIND CLUSTERS---
    vector<vector<Ring>> clusters;
    for (auto& ring: rings) {
        vector<Ring> rings_to_check; //will be filled with rings that still need to be checked for overlapping rings
        bool done;

        vector<Ring> cluster;
        if (ring.cluster_index == -1){    //ring is not yet in a cluster, so we add a new, empty cluster to clusters
            cluster.push_back(ring);  //add this ring to the new cluster
            clusters.push_back(cluster);    // add the new cluster to clusters
            ring.cluster_index = clusters.size()-1;    // update the new cluster_index for the ring
            rings_to_check.push_back(ring);
            done = false;
        } else {
            done = true;
        }
        
        while (not done) {
            auto& ring_ = rings_to_check.back();   //we now look for rings that overlap with the last ring in the list
            ring_.cluster_index = clusters.size()-1;    //first I want to add the cluster_index to the ring, which is why I used a reference in the previous line
            rings_to_check.pop_back();  
            
            auto ring = ring_;  //we continue with a copy of ring_, because ring_ was defined as the last element of a vector. So if we chance that vector in the following lines, ring_ would also change, which is not what I want
            for (auto& ring2: rings){  //look if the ring has overlap with other rings
                if (ring.check_overlap(ring2)){      //found a connecting ring that is not in the cluster
                    if (ring2.cluster_index == -1){                  
                        clusters[ring.cluster_index].push_back(ring2);     //add ring2 to the cluster   
                        ring2.cluster_index = clusters.size()-1;           //assign the cluster_index to ring2
                        rings_to_check.push_back(ring2);                   //add the ring to rings_to_check
                    }
                }        
            }
            if (rings_to_check.size() == 0) {   //the whole cluster has been searched through
                done = true;
            }
        }


        cout << "CLUSTERS:\n";
        for (auto cluster: clusters) {
            for (auto ring: cluster) {
                // cout << ring.get_ring_index() << ": (" << ring.get_x() << " "<< ring.get_y() <<") ";
                cout << ring.get_ring_index() << " ";
            }
            cout << endl;
        }
    }


    //---CALCULATE JUNCTIONS---
    for (auto& ring: rings) {
        ring.calc_all_junctions(rings);
    }    







    cout << "end of main";

}
