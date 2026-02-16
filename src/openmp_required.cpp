#if !defined(_OPENMP) && !defined(GFLOW_ALLOW_NO_OPENMP)
#error "gflow dev profile requires OpenMP. Configure an OpenMP-enabled toolchain or set GFLOW_BUILD_PROFILE=cran-safe."
#endif

extern "C" void gflow_openmp_requirement_sentinel(void) {}
