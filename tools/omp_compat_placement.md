# omp_compat.h include placement

Keep `omp_compat.h` **early**, not last.

### Why

The whole point of `omp_compat.h` is to **undef any troublesome macros (like `match`) *before* `<omp.h>` is parsed**, and then include `<omp.h>` safely (or provide stubs when OpenMP is off). If you include it late, some other header (or your own) might already have included `<omp.h>` or used OpenMP pragmas — at which point the macro collision has already happened.

### Recommended pattern

1. **Never include `<omp.h>` directly anywhere.** Only include your wrapper.
2. Put this at (or very near) the **top of any .cpp/.hpp that uses OpenMP**:

```cpp
#include "omp_compat.h"   // must come before anything that might include <omp.h>
```

3. In `omp_compat.h`:

```cpp
#pragma once

// Kill problematic macros before <omp.h> is seen
#ifdef match
#  undef match
#endif
#ifdef check       // some SDKs define this one too; harmless to guard
#  undef check
#endif

#if defined(_OPENMP)
  #include <omp.h>
#else
  // tiny stubs so you can call these without #ifdefs elsewhere
  inline int  omp_get_max_threads() { return 1; }
  inline int  omp_get_thread_num()  { return 0; }
  inline void omp_set_num_threads(int) {}
#endif
```

### Sanity checks

* Make sure nothing includes `<omp.h>` directly:

  ```bash
  grep -RIn '^[[:space:]]*#include[[:space:]]*<omp\.h>' src inst/include
  ```

  If you find hits, replace them with `#include "omp_compat.h"`.

* If you want to be extra safe, you can have `omp_compat.h` define a marker and fail builds that include `<omp.h>` directly afterwards, but the grep is usually enough.

### TL;DR

Your instinct to include project headers first is fine — **as long as `omp_compat.h` is included before any code that might include `<omp.h>` or use OpenMP pragmas.** Don’t move it to the end; keep it at the top (or make other headers that rely on OpenMP include it themselves at the top).

