#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include "Logger.h"


// Use this macro at the beginning of a TEST to log the start of the test case

#define LOG_STARTING_OF_TEST     \
    const ::testing::TestInfo* const test_info = ::testing::UnitTest::GetInstance()->current_test_info(); \
    Log.info("Running GTEST: " + std::string(test_info->test_case_name()) + ", " + std::string(test_info->name()));

#endif /* DEFINITIONS_H */
