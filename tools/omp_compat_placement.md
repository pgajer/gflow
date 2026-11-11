# OpenMP compatibility header placement

Include `inst/include/gflow/omp_compat.h` before any header that may include
`omp.h`. The compatibility header removes conflicting `match` and `check`
macros before including the OpenMP declarations.

In gflow-owned sources, use:

```cpp
#include "omp_compat.h"
```

The header provides `gflow_get_max_threads()`, `gflow_get_thread_num()` and
other `gflow_*` wrappers. When OpenMP is disabled, its serial implementations
report one thread and make thread-setting calls no-ops. Use these wrappers
for code that must also compile without OpenMP.

To inspect direct includes:

```sh
rg -n '#include[[:space:]]*<omp\.h>' src inst/include/gflow
```

For gflow-owned callers, route direct OpenMP includes through the compatibility
header. Keep the compatibility header itself as the point that includes
`omp.h`; bundled third-party headers follow their own include conventions.
