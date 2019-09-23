#ifndef CSV_PARSE_H
#define CSV_PARSE_H

// LIBRARIES
#include <string>

// FUNCTIONS
void load_geom(std::string csv_path, Geometry &geom, std::string test_tag);
void save_results(std::string out_dir, std::string time_tag, Geometry geom, std::string load_tag);

#endif