# There isn’t a single “law,” but for C++ code in R packages a sane, low-drama include order is:

1. **This file’s own header** (if any) — catches missing prerequisites early.
2. **C/C++ standard library headers** — e.g. `<cstddef>`, `<vector>`, `<algorithm>`, `<map>`, `<string>`, `<fstream>`, `<cmath>`, etc.
3. **Third-party C++ headers** — e.g. Eigen.
4. **R C API headers** — `R.h` → `Rinternals.h` → extras like `R_ext/Print.h`; wrap in `extern "C"`.

Applied to your snippet:

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

### Why this order?

* **Self-contained headers:** Including your own header first (in a .cpp) ensures it compiles with only what it includes.
* **Reduce macro collisions:** Keeping R headers **last** minimizes macro pollution affecting template-heavy libs. With `-DR_NO_REMAP` you’re mostly safe either way, but this order avoids surprises.
* **C linkage:** R’s headers declare C functions; the `extern "C"` block ensures correct linkage in C++ (R headers handle this on many platforms, but wrapping is harmless and explicit).

### A few extra tips

* Put `using std::size_t;` **after** `<cstddef>`. Never before.
* In **public headers**, many teams avoid any `using` and write `std::size_t` explicitly. If changing all uses is too much right now, your narrow `using std::size_t;` is acceptable—just keep it to headers you control and avoid `using namespace std;`.
* Only include what you use; don’t pull in heavy headers (e.g., `<Eigen/Dense>`) from widely included headers if a forward-decl or a slimmer Eigen header would do (Eigen is header-only, so forward-decls are limited, but you can prefer `<Eigen/Core>` when you don’t need factorizations).

If you want, paste one of your headers and I’ll reorganize it to this pattern with minimal edits.


