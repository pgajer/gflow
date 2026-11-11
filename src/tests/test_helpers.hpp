// src/tests/test_helpers.hpp
#pragma once
#include <iostream>
#include <cmath>
#include <string>

#define TEST(name) \
    void test_##name(); \
    namespace { struct Register_##name { \
        Register_##name() { register_test(#name, test_##name); } \
    } reg_##name; } \
    void test_##name()

#define ASSERT_NEAR(a, b, tol) \
    if (std::abs((a) - (b)) > (tol)) { \
        std::cerr << "FAIL: " << #a << " != " << #b \
                  << " (" << (a) << " vs " << (b) << ")\n"; \
        throw std::runtime_error("Assertion failed"); \
    }

#define ASSERT_TRUE(cond) \
    if (!(cond)) { \
        std::cerr << "FAIL: " << #cond << "\n"; \
        throw std::runtime_error("Assertion failed"); \
    }

// Global test registry
inline std::vector<std::pair<std::string, void(*)()>>& get_tests() {
    static std::vector<std::pair<std::string, void(*)()>> tests;
    return tests;
}

inline void register_test(const char* name, void(*func)()) {
    get_tests().push_back({name, func});
}

inline int run_all_tests() {
    int passed = 0, failed = 0;
    for (auto& [name, func] : get_tests()) {
        try {
            func();
            std::cout << "✓ " << name << "\n";
            passed++;
        } catch (const std::exception& e) {
            std::cerr << "✗ " << name << ": " << e.what() << "\n";
            failed++;
        }
    }
    std::cout << "\n" << passed << " passed, " << failed << " failed\n";
    return failed;
}
