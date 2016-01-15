#include "gtest/gtest.h"

#include "FileUtilities.h"

TEST(FileUtilitiesTest, fileExists) {
    
    ASSERT_TRUE(FileUtilities::fileExists("../testData/input.yaml"));

    ASSERT_FALSE(FileUtilities::fileExists("../wrongDirectory/input.yaml"));

}