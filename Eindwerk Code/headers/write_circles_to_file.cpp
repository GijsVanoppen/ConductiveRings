#include "write_circles_to_file.h"
#include <fstream>
void write_circles_to_file(std::valarray<std::array<double,3>> circles, std::string file_name) {
    std::ofstream circles_file(file_name);
    for (const auto& circle: circles) {
        circles_file << circle.front() << " " << circle[1] << " " << circle.back() << "\n";
    }
    circles_file.close();
}