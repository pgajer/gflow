#ifndef MEMORY_UTILS_HPP
#define MEMORY_UTILS_HPP

#include <cstddef>  // for size_t
#include <R.h>     // for Rprintf and R_FlushConsole

// Platform-specific includes for memory tracking
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>
#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>
#endif
#endif

// Function to get current RSS (Resident Set Size)
size_t get_current_rss();  // Platform-specific implementation

/**
 * @brief Memory usage tracking utility class
 *
 * Provides functionality to track and report memory usage changes between
 * initialization and subsequent checkpoints. Memory usage is measured using
 * Resident Set Size (RSS).
 *
 * Example usage:
 * @code
 *     {
 *         memory_tracker_t tracker("MyOperation");
 *         // ... perform memory-intensive operations ...
 *         tracker.report(); // Reports current usage and difference
 *     }
 * @endcode
 *
 * @note Memory measurements are in bytes but reported in megabytes (MB).
 *       The report includes both absolute usage and relative change.
 */
struct memory_tracker_t {
    size_t initial_usage; ///< Initial RSS measurement taken at construction time
    const char* context;  ///< Context string used to identify the tracking session in reports

    /**
     * @brief Constructs a new memory tracker
     *
     * @param ctx Context string to identify this tracking session in reports
     *
     * @note Takes an initial RSS measurement upon construction
     */
    memory_tracker_t(const char* ctx) : context(ctx) {
        initial_usage = get_current_rss();
    }

    /**
     * @brief Reports current memory usage and change since construction
     *
     * Prints a formatted message to R's console showing:
     *   - Current RSS memory usage in MB
     *   - Change in RSS memory usage since construction in MB (with sign)
     *
     * Output format: "<context> Memory: Current=X MB, Diff=±Y MB"
     *
     * @note Automatically flushes the R console after printing
     */
    void report() {
        size_t current = get_current_rss();
        size_t diff = current - initial_usage;
        Rprintf("%s Memory: Current=%zu MB, Diff=%+zd MB\n",
                context, current/1024/1024, diff/1024/1024);
        R_FlushConsole();
    }
};

#endif // MEMORY_UTILS_HPP
