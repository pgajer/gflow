#include "memory_utils.hpp"

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
size_t get_current_rss() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
        return pmc.WorkingSetSize;
    }
    return 0;

#elif defined(__APPLE__) && defined(__MACH__)
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                  (task_info_t)&info, &infoCount) == KERN_SUCCESS) {
        return info.resident_size;
    }
    return 0;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    FILE* fp = fopen("/proc/self/statm", "r");
    if (fp) {
        long rss;
        if (fscanf(fp, "%*s%ld", &rss) == 1) {
            fclose(fp);
            return rss * sysconf(_SC_PAGESIZE);
        }
        fclose(fp);
    }
    return 0;

#else
    return 0;
#endif
}
