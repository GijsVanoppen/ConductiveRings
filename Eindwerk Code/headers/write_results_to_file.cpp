#include "write_results_to_file.h"
#include <fstream>
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