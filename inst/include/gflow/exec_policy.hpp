#pragma once
#include <execution>

#if defined(_WIN32)
// Windows: avoid oneTBB dependency
#  define GFLOW_EXEC_POLICY std::execution::seq
#else
// Non-Windows: keep your fast path
#  define GFLOW_EXEC_POLICY std::execution::par_unseq
#endif
