#include "extract_parameter_from_string.h"
std::string extract_parameter_from_string(std::string str) {
    //returns the part of the string after the ':' sign. Used for handle_input
    std::string parameter {""};
    int start_index = str.find(':') ;
    for(int i {start_index+1}; i < str.size(); i++) {
        parameter += str[i];
    }
    return parameter;
}