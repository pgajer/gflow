/**
 * @brief Platform-independent memory tracking utility function
 *
 * Retrieves the current Resident Set Size (RSS) of the running process.
 * RSS represents the portion of a process's memory that is held in RAM.
 *
 * @note Platform-specific implementations:
 *   - Windows: Uses GetProcessMemoryInfo() to get WorkingSetSize
 *   - macOS: Uses task_info() with MACH_TASK_BASIC_INFO
 *   - Linux: Reads from /proc/self/statm and multiplies by page size
 *
 * @return size_t The current RSS in bytes. Returns 0 if measurement fails
 *         or if the platform is unsupported.
 */
#include "memory_utils.hpp"

#ifdef _WIN32
size_t get_current_rss() {
    PROCESS_MEMORY_COUNTERS_EX pmc;
    HANDLE h = GetCurrentProcess();
    if (GetProcessMemoryInfo(h, (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc))) {
        return (size_t)pmc.WorkingSetSize; // bytes
    }
    return 0;
}
#elif defined(__APPLE__) && defined(__MACH__)
size_t get_current_rss() {
    task_basic_info_data_t info;
    mach_msg_type_number_t count = TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&info, &count) == KERN_SUCCESS) {
        return (size_t)info.resident_size;
    }
    return 0;
}
#else
size_t get_current_rss() {
    // Linux/Unix: prefer /proc/self/statm if available; fall back to getrusage
    // 1) /proc path (Linux)
    FILE* f = std::fopen("/proc/self/statm", "r");
    if (f) {
        long pages = 0;
        if (std::fscanf(f, "%*s %ld", &pages) == 1) {
            std::fclose(f);
            long page_size = sysconf(_SC_PAGESIZE);
            if (page_size > 0 && pages > 0) return (size_t)pages * (size_t)page_size;
        } else {
            std::fclose(f);
        }
    }
    // 2) getrusage (portable-ish fallback)
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) == 0) {
    #ifdef __linux__
        return (size_t)ru.ru_maxrss * 1024ULL; // Linux: ru_maxrss in kilobytes
    #else
        return (size_t)ru.ru_maxrss;           // BSD/macOS: ru_maxrss in bytes (not used on mac path)
    #endif
    }
    return 0;
}
#endif
