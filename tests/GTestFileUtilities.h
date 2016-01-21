#include "gtest/gtest.h"

#include "FileUtilities.h"

TEST(FileUtilitiesTest, fileExists) {
    
    ASSERT_TRUE(FileUtilities::fileExists("../testData/input.yaml"));

    ASSERT_FALSE(FileUtilities::fileExists("../wrongDirectory/input.yaml"));

}


TEST(FileUtilitiesTest, isRelative) {
    
    ASSERT_TRUE(FileUtilities::isRelative("../testData/input.yaml"));

    ASSERT_FALSE(FileUtilities::isRelative("/Users/rik/Git/PlatoSim3/testData/input.yaml"));

}

