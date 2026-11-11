# Header include order

For gflow-owned C++ sources, use this order:

1. The source file's own header, where applicable, to check that it is self-contained.
2. Required C and C++ standard-library headers.
3. Third-party C++ headers, such as Eigen.
4. R C API headers, with `R.h` before `Rinternals.h`.

Place the OpenMP compatibility header before any header that can include
`omp.h`; see [OpenMP header placement](omp_compat_placement.md).

A header that uses standard containers, Eigen and the R C API can use:

```cpp
#pragma once

// (optional) Windows macro guard to avoid min/max conflicts if you ever include <windows.h>
#ifdef _WIN32
  #ifndef NOMINMAX
    #define NOMINMAX
  #endif
#endif

// 1) C++ standard library
#include <cstddef>     // must precede any 'using std::size_t;'
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <fstream>

// 2) Third-party C++ (header-only)
#include <Eigen/Dense> // Eigen::VectorXd, Eigen::MatrixXd

// 3) R C API (order matters: R.h before Rinternals.h)
extern "C" {
  #include <R.h>
  #include <Rinternals.h>
  #include <R_ext/Print.h>  // Rprintf
}

// 4) Narrow using-declarations (avoid broad ones in headers)
using std::size_t;     // if you keep this pattern; otherwise prefer std::size_t in headers
```

Including R headers after template-heavy libraries limits macro collisions;
`R_NO_REMAP` additionally suppresses R's compatibility remapping macros.
The example uses explicit C linkage for the R API declarations.

Include `<cstddef>` before referring to `std::size_t`. Prefer qualified names
in public headers and avoid `using namespace std`. Include only the headers
needed by the declarations or implementation; for example, use `Eigen/Core`
when decomposition facilities from `Eigen/Dense` are unnecessary.
