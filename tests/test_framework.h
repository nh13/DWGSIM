/*
 * Minimal Unit Test Framework for DWGSIM
 *
 * Usage:
 *   #include "test_framework.h"
 *
 *   TEST(test_name) {
 *       ASSERT(condition);
 *       ASSERT_EQ(expected, actual);
 *       ASSERT_STR_EQ(expected, actual);
 *   }
 *
 *   int main(void) {
 *       RUN_TEST(test_name);
 *       TEST_SUMMARY();
 *       return test_failures;
 *   }
 */

#ifndef TEST_FRAMEWORK_H
#define TEST_FRAMEWORK_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* Global test counters */
static int test_count = 0;
static int test_failures = 0;
static int test_assertions = 0;
static const char *current_test = NULL;

/* Test function signature */
typedef void (*test_fn)(void);

/* Define a test function */
#define TEST(name) static void name(void)

/* Run a test */
#define RUN_TEST(name) do { \
    current_test = #name; \
    test_count++; \
    int prev_failures = test_failures; \
    name(); \
    if (test_failures == prev_failures) { \
        printf("  [PASS] %s\n", #name); \
    } \
} while (0)

/* Print test summary */
#define TEST_SUMMARY() do { \
    printf("\n----------------------------------------\n"); \
    printf("Tests: %d | Passed: %d | Failed: %d | Assertions: %d\n", \
           test_count, test_count - test_failures, test_failures, test_assertions); \
    printf("----------------------------------------\n"); \
} while (0)

/* Basic assertion */
#define ASSERT(cond) do { \
    test_assertions++; \
    if (!(cond)) { \
        printf("  [FAIL] %s\n", current_test); \
        printf("         Assertion failed: %s\n", #cond); \
        printf("         at %s:%d\n", __FILE__, __LINE__); \
        test_failures++; \
        return; \
    } \
} while (0)

/* Assert with message */
#define ASSERT_MSG(cond, msg) do { \
    test_assertions++; \
    if (!(cond)) { \
        printf("  [FAIL] %s\n", current_test); \
        printf("         %s\n", msg); \
        printf("         at %s:%d\n", __FILE__, __LINE__); \
        test_failures++; \
        return; \
    } \
} while (0)

/* Assert equality for integers */
#define ASSERT_EQ(expected, actual) do { \
    test_assertions++; \
    long long _exp = (long long)(expected); \
    long long _act = (long long)(actual); \
    if (_exp != _act) { \
        printf("  [FAIL] %s\n", current_test); \
        printf("         Expected: %lld, Actual: %lld\n", _exp, _act); \
        printf("         at %s:%d\n", __FILE__, __LINE__); \
        test_failures++; \
        return; \
    } \
} while (0)

/* Assert inequality */
#define ASSERT_NE(not_expected, actual) do { \
    test_assertions++; \
    long long _nexp = (long long)(not_expected); \
    long long _act = (long long)(actual); \
    if (_nexp == _act) { \
        printf("  [FAIL] %s\n", current_test); \
        printf("         Expected NOT: %lld, but got: %lld\n", _nexp, _act); \
        printf("         at %s:%d\n", __FILE__, __LINE__); \
        test_failures++; \
        return; \
    } \
} while (0)

/* Assert string equality */
#define ASSERT_STR_EQ(expected, actual) do { \
    test_assertions++; \
    const char *_exp = (expected); \
    const char *_act = (actual); \
    if (_exp == NULL || _act == NULL || strcmp(_exp, _act) != 0) { \
        printf("  [FAIL] %s\n", current_test); \
        printf("         Expected: \"%s\", Actual: \"%s\"\n", \
               _exp ? _exp : "(null)", _act ? _act : "(null)"); \
        printf("         at %s:%d\n", __FILE__, __LINE__); \
        test_failures++; \
        return; \
    } \
} while (0)

/* Assert floating point equality within epsilon */
#define ASSERT_FLOAT_EQ(expected, actual, epsilon) do { \
    test_assertions++; \
    double _exp = (double)(expected); \
    double _act = (double)(actual); \
    double _eps = (double)(epsilon); \
    if (fabs(_exp - _act) > _eps) { \
        printf("  [FAIL] %s\n", current_test); \
        printf("         Expected: %g, Actual: %g (epsilon: %g)\n", _exp, _act, _eps); \
        printf("         at %s:%d\n", __FILE__, __LINE__); \
        test_failures++; \
        return; \
    } \
} while (0)

/* Assert pointer is not NULL */
#define ASSERT_NOT_NULL(ptr) do { \
    test_assertions++; \
    if ((ptr) == NULL) { \
        printf("  [FAIL] %s\n", current_test); \
        printf("         Expected non-NULL pointer\n"); \
        printf("         at %s:%d\n", __FILE__, __LINE__); \
        test_failures++; \
        return; \
    } \
} while (0)

/* Assert pointer is NULL */
#define ASSERT_NULL(ptr) do { \
    test_assertions++; \
    if ((ptr) != NULL) { \
        printf("  [FAIL] %s\n", current_test); \
        printf("         Expected NULL pointer\n"); \
        printf("         at %s:%d\n", __FILE__, __LINE__); \
        test_failures++; \
        return; \
    } \
} while (0)

/* Assert a < b */
#define ASSERT_LT(a, b) do { \
    test_assertions++; \
    long long _a = (long long)(a); \
    long long _b = (long long)(b); \
    if (!(_a < _b)) { \
        printf("  [FAIL] %s\n", current_test); \
        printf("         Expected %lld < %lld\n", _a, _b); \
        printf("         at %s:%d\n", __FILE__, __LINE__); \
        test_failures++; \
        return; \
    } \
} while (0)

/* Assert a <= b */
#define ASSERT_LE(a, b) do { \
    test_assertions++; \
    long long _a = (long long)(a); \
    long long _b = (long long)(b); \
    if (!(_a <= _b)) { \
        printf("  [FAIL] %s\n", current_test); \
        printf("         Expected %lld <= %lld\n", _a, _b); \
        printf("         at %s:%d\n", __FILE__, __LINE__); \
        test_failures++; \
        return; \
    } \
} while (0)

/* Assert a > b */
#define ASSERT_GT(a, b) do { \
    test_assertions++; \
    long long _a = (long long)(a); \
    long long _b = (long long)(b); \
    if (!(_a > _b)) { \
        printf("  [FAIL] %s\n", current_test); \
        printf("         Expected %lld > %lld\n", _a, _b); \
        printf("         at %s:%d\n", __FILE__, __LINE__); \
        test_failures++; \
        return; \
    } \
} while (0)

/* Assert a >= b */
#define ASSERT_GE(a, b) do { \
    test_assertions++; \
    long long _a = (long long)(a); \
    long long _b = (long long)(b); \
    if (!(_a >= _b)) { \
        printf("  [FAIL] %s\n", current_test); \
        printf("         Expected %lld >= %lld\n", _a, _b); \
        printf("         at %s:%d\n", __FILE__, __LINE__); \
        test_failures++; \
        return; \
    } \
} while (0)

/* Test suite name header */
#define TEST_SUITE(name) printf("\n=== %s ===\n", name)

#endif /* TEST_FRAMEWORK_H */
