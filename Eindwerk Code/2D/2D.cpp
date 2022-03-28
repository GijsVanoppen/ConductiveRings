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

    vector<array<double, 4>> junctions;
    double get_x() {return x_;}
    double get_y() {return y_;}
    double get_r() {return r_;}
    double get_R() {return R_;}
    int get_ring_index() {return ring_index_;}

    double distance_between_rings(Ring ring1, Ring ring2) {
        return sqrt(pow(ring1.get_x() - ring2.get_x(),2) + pow(ring1.get_y() - ring2.get_y(),2));
    }

    bool check_overlap(Ring ring1, Ring ring2) {
        //check if the two Ring objects are not the same
        if (ring1.get_x() == ring2.get_x() && ring1.get_y() == ring2.get_y() && ring1.get_r() == ring2.get_r() && ring1.get_R() == ring2.get_R()) {
            return false;
        }

        if (distance_between_rings(ring1, ring2) < (ring1.get_r() + ring2.get_r())) {
            return true;
        }
        else {
            return false;
        }
    }


    int cluster_index = -1;

    int calc_cluster_index(vector<vector<int>> clusters) {   //returns index of clusters in which vector the ring index is located. Returns -1 if this ring is not in a cluster
        for (int cluster_index_ = 0; cluster_index_ < clusters.size(); cluster_index_++) {
            for (const auto& value: clusters[cluster_index_]) {
                if (value == get_ring_index()) {
                    return cluster_index_;
                    cout << "this ring is already in cluster\n ";
                }
            }
        }
        return -1;
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
        rings.push_back(ring);          //filling vector
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
    int N = 1000;
    vector<Ring> rings = generate_rings(N, 10, 10, 1, 2, 1); 
    write_rings_to_file(rings, "rings.txt");


    //---FIND CLUSTERS---
    vector<vector<int>> clusters;
    
    for (int ring_index=0; ring_index<N; ring_index++) {
        auto ring = rings[ring_index];  //this is the ring that is currently focussed
        ring.cluster_index = ring.calc_cluster_index(clusters);

        vector<int> cluster;
        if (ring.cluster_index == -1){ //ring is not yet in a cluster, so we add a new, empty cluster to clusters
            cluster.push_back(ring_index);  //add the index of this ring to the new cluster
            clusters.push_back(cluster);    // add the new cluster to clusters
            ring.cluster_index = clusters.size()-1;     // update the new cluster_index for the ring
        } else {
            cluster = clusters[ring.cluster_index];     //ring is already in a cluster, so we want to expand that cluster with connecting rings
        }
        
        for (int ring_index2 = 0; ring_index2 < N; ring_index2++){  //look if the ring has overlap with other rings
            if (ring.check_overlap(ring, rings[ring_index2])){      //found a connecting ring
                if (find(cluster.begin(), cluster.end(), ring_index2) == cluster.end()) {   //the connected ring is not yet in the cluster
                    clusters[ring.cluster_index].push_back(ring_index2);     //add the ring_index to the cluster    
                }
            } 
    
        }
    }




}
