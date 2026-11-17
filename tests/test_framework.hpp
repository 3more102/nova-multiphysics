#pragma once

#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nova::test {

struct TestCase {
    std::string name;
    void (*function)();
};

inline std::vector<TestCase>& registry() {
    static std::vector<TestCase> tests;
    return tests;
}

class Registrar {
public:
    Registrar(std::string name, void (*function)()) {
        registry().push_back({std::move(name), function});
    }
};

inline void require_near(
    const double actual,
    const double expected,
    const double tolerance,
    const char* expression,
    const char* file,
    const int line
) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(actual - expected) > tolerance) {
        std::ostringstream message;
        message << file << ':' << line << ": " << expression
                << " (actual=" << actual << ", expected=" << expected
                << ", tolerance=" << tolerance << ')';
        throw std::runtime_error(message.str());
    }
}

inline int run_all_tests() {
    std::size_t failures = 0;

    for (const TestCase& test : registry()) {
        try {
            test.function();
            std::cout << "[PASS] " << test.name << '\n';
        } catch (const std::exception& error) {
            ++failures;
            std::cerr << "[FAIL] " << test.name << ": " << error.what() << '\n';
        } catch (...) {
            ++failures;
            std::cerr << "[FAIL] " << test.name << ": unknown exception\n";
        }
    }

    std::cout << "\nExecuted " << registry().size() << " tests; "
              << failures << " failed.\n";
    return failures == 0 ? 0 : 1;
}

} // namespace nova::test

#define NOVA_TEST(name) \
    static void name(); \
    static const ::nova::test::Registrar name##_registrar(#name, &name); \
    static void name()

#define NOVA_REQUIRE(condition) \
    do { \
        if (!(condition)) { \
            throw std::runtime_error( \
                std::string(__FILE__) + ":" + std::to_string(__LINE__) + \
                ": requirement failed: " #condition \
            ); \
        } \
    } while (false)

#define NOVA_REQUIRE_NEAR(actual, expected, tolerance) \
    ::nova::test::require_near( \
        static_cast<double>(actual), \
        static_cast<double>(expected), \
        static_cast<double>(tolerance), \
        #actual " ~= " #expected, \
        __FILE__, \
        __LINE__ \
    )
