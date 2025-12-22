/*
 * DWGSIM Unit Test Runner
 *
 * This file runs all unit tests for DWGSIM.
 * Individual test files are included below.
 */

#include "test_framework.h"

/* Include test files here as they are added */
/* #include "test_position.c" */
/* #include "test_memory.c" */
/* #include "test_parsing.c" */

/*
 * Placeholder tests to verify the test framework works
 */
TEST(test_framework_assert) {
    ASSERT(1 == 1);
    ASSERT(0 == 0);
}

TEST(test_framework_assert_eq) {
    ASSERT_EQ(42, 42);
    ASSERT_EQ(-1, -1);
    ASSERT_EQ(0, 0);
}

TEST(test_framework_assert_ne) {
    ASSERT_NE(1, 2);
    ASSERT_NE(-1, 1);
}

TEST(test_framework_assert_str_eq) {
    ASSERT_STR_EQ("hello", "hello");
    ASSERT_STR_EQ("", "");
}

TEST(test_framework_assert_float_eq) {
    ASSERT_FLOAT_EQ(3.14, 3.14, 0.001);
    ASSERT_FLOAT_EQ(0.1 + 0.2, 0.3, 0.0001);
}

TEST(test_framework_assert_comparisons) {
    ASSERT_LT(1, 2);
    ASSERT_LE(1, 1);
    ASSERT_LE(1, 2);
    ASSERT_GT(2, 1);
    ASSERT_GE(2, 2);
    ASSERT_GE(3, 2);
}

TEST(test_framework_assert_null) {
    int *ptr = NULL;
    int val = 42;
    ASSERT_NULL(ptr);
    ASSERT_NOT_NULL(&val);
}

int main(void) {
    printf("DWGSIM Unit Tests\n");
    printf("========================================\n");

    TEST_SUITE("Test Framework Verification");
    RUN_TEST(test_framework_assert);
    RUN_TEST(test_framework_assert_eq);
    RUN_TEST(test_framework_assert_ne);
    RUN_TEST(test_framework_assert_str_eq);
    RUN_TEST(test_framework_assert_float_eq);
    RUN_TEST(test_framework_assert_comparisons);
    RUN_TEST(test_framework_assert_null);

    /* Add test suites here as they are created */
    /* TEST_SUITE("Position Calculation Tests"); */
    /* Include tests from test_position.c */

    TEST_SUMMARY();

    return test_failures > 0 ? 1 : 0;
}
