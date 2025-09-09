#!/bin/bash

# Replace simple std::cout << std::endl patterns
sed -i '' 's/std::cout << std::endl;/Rprintf("\\n");/g' cpp_utils.hpp
sed -i '' 's/std::cout << std::endl;/Rprintf("\\n");/g' ext_ulm_priority_queue.hpp
sed -i '' 's/std::cout << std::endl;/Rprintf("\\n");/g' grid_vertex_path_model_priority_queue.hpp

# Replace std::cout << name << std::endl patterns
sed -i '' 's/std::cout << name << std::endl;/Rprintf("%s\\n", name.c_str());/g' cpp_utils.hpp
sed -i '' 's/std::cout << name << std::endl;/Rprintf("%s\\n", name.c_str());/g' ext_ulm_priority_queue.hpp
sed -i '' 's/std::cout << name << std::endl;/Rprintf("%s\\n", name.c_str());/g' grid_vertex_path_model_priority_queue.hpp

# Replace std::cout << name << ":" << std::endl patterns
sed -i '' 's/std::cout << name << ":" << std::endl;/Rprintf("%s:\\n", name.c_str());/g' cpp_utils.hpp

# Replace simple string literals
sed -i '' 's/std::cout << "  ";/Rprintf("  ");/g' cpp_utils.hpp
sed -i '' 's/std::cout << ", ";/Rprintf(", ");/g' cpp_utils.hpp
sed -i '' 's/std::cout << " ";/Rprintf(" ");/g' cpp_utils.hpp
sed -i '' 's/std::cout << "";/Rprintf("");/g' cpp_utils.hpp

