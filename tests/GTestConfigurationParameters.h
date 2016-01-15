#include <string>
#include <vector>
#include <list>

#include "gtest/gtest.h"

#include "ConfigurationParameters.h"
#include "Exceptions.h"

using namespace std;

TEST(ConfigurationParametersTest, Constructor) {

    ASSERT_THROW(ConfigurationParameters ip = ConfigurationParameters("input.yaml"), IOException);

}

TEST(ConfigurationParametersTest, readGlobalValues) {
    ConfigurationParameters ip = ConfigurationParameters("../testData/input.yaml");

    string description = ip.getString("Description");
    EXPECT_STREQ(description.c_str(), "YAML Input File for 3rd Generation PLATO Simulator");

    string author = ip.getString("Author");
    EXPECT_STREQ(author.c_str(), "Rik Huygen");
}

TEST(ConfigurationParametersTest, readGeneralValues) {
    ConfigurationParameters ip = ConfigurationParameters("../testData/input.yaml");

    string projectLocation = ip.getString("General/ProjectLocation");
    EXPECT_STREQ(projectLocation.c_str(), "/Users/rik/Work/PLATO");
}

TEST(ConfigurationParametersTest, readObservingValues) {
    ConfigurationParameters ip = ConfigurationParameters("../testData/input.yaml");

    int exposureTime = ip.getInteger("Observing/ExposureTime");
    EXPECT_EQ(23, exposureTime);

    string filename = ip.getString("Observing/StarCatalogueFileName");
    EXPECT_STREQ(filename.c_str(), "inputFiles/starField_RA180Dec-70.txt");

    filename = ip.getProjectFileName("Observing/StarCatalogueFileName");
    EXPECT_STREQ(filename.c_str(), "/Users/rik/Work/PLATO/inputFiles/starField_RA180Dec-70.txt");

    filename = ip.getProjectFileName("Observing/AbsoluteFileName");
    EXPECT_STREQ(filename.c_str(), "/Users/rik/Work/PLATO/inputFiles/someInputFile.txt");

    double area = ip.getDouble("Observing/LightCollectingArea");
    EXPECT_DOUBLE_EQ(0.1131, area);
}

TEST(ConfigurationParametersTest, readSpecialValues) {
    ConfigurationParameters ip = ConfigurationParameters("../testData/input.yaml");

    int zeroValue = ip.getInteger("Special Values/zero");
    EXPECT_EQ(0, zeroValue);

    int oneValue = ip.getInteger("Special Values/one");
    EXPECT_EQ(1, oneValue);

    int minusOneValue = ip.getInteger("Special Values/minus-one");
    EXPECT_EQ(-1, minusOneValue);

    bool booleanTrue = ip.getBoolean("Special Values/boolean-true");
    EXPECT_TRUE(booleanTrue);

    bool booleanFalse = ip.getBoolean("Special Values/boolean-false");
    EXPECT_FALSE(booleanFalse);

}

