//
// IO C++ utilities
//
#include <R_ext/Error.h>
#include <fstream>
#include <vector>

// Writes the content of a vector of vector of int to a binary file
void write_binary_vect_vect_int(const std::vector<std::vector<int>>& vv, const std::string& out_file) {
    std::ofstream out(out_file, std::ios::binary);
    if (!out.is_open()) {
        Rf_error("Cannot open file for writing: %s", out_file.c_str());
    }

    // Write the outer vector size
    size_t outer_size = vv.size();
    out.write(reinterpret_cast<const char*>(&outer_size), sizeof(outer_size));

    for (const auto& inner_vector : vv) {
        // Write the inner vector size
        size_t inner_size = inner_vector.size();
        out.write(reinterpret_cast<const char*>(&inner_size), sizeof(inner_size));

        // Write the inner vector elements
        out.write(reinterpret_cast<const char*>(inner_vector.data()), inner_size * sizeof(int));
    }

    out.close();
}

// Reads the content of a binary file to a vector of vector of int
void read_binary_vect_vect_int(std::vector<std::vector<int>>& vv, const std::string& in_file) {
    std::ifstream in(in_file, std::ios::binary);
    if (!in.is_open()) {
        Rf_error("Cannot open file for reading: %s", in_file.c_str());
    }

    // Read the outer vector size
    size_t outer_size;
    in.read(reinterpret_cast<char*>(&outer_size), sizeof(outer_size));
    vv.resize(outer_size);

    for (auto& inner_vector : vv) {
        // Read the inner vector size
        size_t inner_size;
        in.read(reinterpret_cast<char*>(&inner_size), sizeof(inner_size));
        inner_vector.resize(inner_size);

        // Read the inner vector elements
        in.read(reinterpret_cast<char*>(inner_vector.data()), inner_size * sizeof(int));
    }

    in.close();
}
