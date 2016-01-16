#include <string>
#include <vector>
#include <list>

#include "gtest/gtest.h"

#include "ConfigurationParameters.h"
#include "Exceptions.h"

using namespace std;



TEST(ConfigurationParametersTest, Constructor)
{

    ASSERT_THROW(ConfigurationParameters ip = ConfigurationParameters("input.yaml"), IOException);

}





TEST(ConfigurationParametersTest, readGlobalValues)
{
    ConfigurationParameters ip = ConfigurationParameters("../testData/input.yaml");

    string description = ip.getString("Description");
    EXPECT_STREQ(description.c_str(), "YAML Input File for 3rd Generation PLATO Simulator");

    string author = ip.getString("Author");
    EXPECT_STREQ(author.c_str(), "Rik Huygen");
}





TEST(ConfigurationParametersTest, readGeneralValues)
{
    ConfigurationParameters ip = ConfigurationParameters("../testData/input.yaml");

    string projectLocation = ip.getString("General/ProjectLocation");
    EXPECT_STREQ(projectLocation.c_str(), "/Users/rik/Work/PLATO");
}





TEST(ConfigurationParametersTest, readObservingValues)
{
    ConfigurationParameters ip = ConfigurationParameters("../testData/input.yaml");

    int exposureTime = ip.getInteger("Observing/ExposureTime");
    EXPECT_EQ(23, exposureTime);

    string filename = ip.getString("Observing/StarCatalogueFileName");
    EXPECT_STREQ(filename.c_str(), "inputFiles/starField_RA180Dec-70.txt");

    filename = ip.getAbsoluteFileName("Observing/StarCatalogueFileName");
    EXPECT_STREQ(filename.c_str(), "/Users/rik/Work/PLATO/inputFiles/starField_RA180Dec-70.txt");

    filename = ip.getAbsoluteFileName("Observing/AbsoluteFileName");
    EXPECT_STREQ(filename.c_str(), "/Users/rik/Work/PLATO/inputFiles/someInputFile.txt");

    double area = ip.getDouble("Observing/LightCollectingArea");
    EXPECT_DOUBLE_EQ(0.1131, area);
}





TEST(ConfigurationParametersTest, readSpecialValues)
{
    ConfigurationParameters ip = ConfigurationParameters("../testData/input.yaml");

    int zeroValue = ip.getInteger("Special Values/zero");
    EXPECT_EQ(0, zeroValue);

    int oneValue = ip.getInteger("Special Values/one");
    EXPECT_EQ(1, oneValue);

    int minusOneValue = ip.getInteger("Special Values/minus-one");
    EXPECT_EQ(-1, minusOneValue);

    // A 0 or a 1 can not be converted into a Boolean - CHECK THIS!
    ASSERT_ANY_THROW(ip.getBoolean("Special Values/zero"));
    ASSERT_ANY_THROW(ip.getBoolean("Special Values/one"));

    bool booleanTrue = ip.getBoolean("Special Values/boolean-true");
    EXPECT_TRUE(booleanTrue);

    bool booleanFalse = ip.getBoolean("Special Values/boolean-false");
    EXPECT_FALSE(booleanFalse);

}





TEST(ConfigurationParametersTest, testConversions)
{
    ConfigurationParameters ip = ConfigurationParameters("../testData/input.yaml");

    // Can convert an integer value into a double
    double exposureTime = ip.getDouble("Observing/ExposureTime");
    EXPECT_DOUBLE_EQ(23.0, exposureTime);

    // Can convert an integer value into a string
    string exposureTimeString = ip.getString("Observing/ExposureTime");
    EXPECT_STREQ("23", exposureTimeString.c_str());

    // Can not convert a double value into an Integer
    ASSERT_ANY_THROW(ip.getInteger("Observing/LightCollectingArea"));

    // Conversion from Integer 23 to Boolean should also throw an exception
    ASSERT_ANY_THROW(ip.getBoolean("Observing/ExposureTime"));

}


