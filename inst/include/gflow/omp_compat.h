#pragma once

#ifdef _OPENMP
  #include <omp.h>
#else
  // Minimal stubs so code compiles without OpenMP
  inline int  omp_get_max_threads() { return 1; }
  inline int  omp_get_thread_num()  { return 0; }
  inline void omp_set_num_threads(int) {}
#endif
