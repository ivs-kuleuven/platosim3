#include <string>
#include <vector>
#include <list>

#include "gtest/gtest.h"

#include "gtest_definitions.h"

#include "ConfigurationParameters.h"
#include "Exceptions.h"

using namespace std;



TEST(ConfigurationParametersTest, Constructor)
{

    LOG_STARTING_OF_TEST

    ASSERT_THROW(ConfigurationParameters cp = ConfigurationParameters("input.yaml"), IllegalArgumentException);

}





TEST(ConfigurationParametersTest, readGlobalValues)
{
    LOG_STARTING_OF_TEST

    ConfigurationParameters cp = ConfigurationParameters("../testData/input_ConfigurationParametersTest.yaml");

    string description = cp.getString("Description");
    EXPECT_STREQ(description.c_str(), "YAML Input File for 3rd Generation PLATO Simulator");

    string author = cp.getString("Author");
    EXPECT_STREQ(author.c_str(), "Rik Huygen");
}





TEST(ConfigurationParametersTest, readGeneralValues)
{
    LOG_STARTING_OF_TEST

    ConfigurationParameters cp = ConfigurationParameters("../testData/input_ConfigurationParametersTest.yaml");

    string projectLocation = cp.getString("General/ProjectLocation");
    EXPECT_STREQ(projectLocation.c_str(), "/Users/rik/Work/PLATO");

    projectLocation = cp.getString("General/ProjectLocation");
    EXPECT_STREQ(projectLocation.c_str(), "/Users/rik/Work/PLATO");
}





TEST(ConfigurationParametersTest, readObservingValues)
{
    LOG_STARTING_OF_TEST

    ConfigurationParameters cp = ConfigurationParameters("../testData/input_ConfigurationParametersTest.yaml");

    int exposureTime = cp.getInteger("Observing/ExposureTime");
    EXPECT_EQ(23, exposureTime);

    string filename = cp.getString("Observing/StarCatalogueFilename");
    EXPECT_STREQ(filename.c_str(), "inputFiles/starField_RA180Dec-70.txt");

    filename = cp.getAbsoluteFilename("Observing/StarCatalogueFilename");
    EXPECT_STREQ(filename.c_str(), "/Users/rik/Work/PLATO/inputFiles/starField_RA180Dec-70.txt");

    filename = cp.getAbsoluteFilename("Observing/AbsoluteFilename");
    EXPECT_STREQ(filename.c_str(), "/Users/rik/Work/PLATO/inputFiles/someInputFile.txt");

    double area = cp.getDouble("Observing/LightCollectingArea");
    EXPECT_DOUBLE_EQ(0.1131, area);
}





TEST(ConfigurationParametersTest, readSpecialValues)
{
    LOG_STARTING_OF_TEST

    ConfigurationParameters cp = ConfigurationParameters("../testData/input_ConfigurationParametersTest.yaml");

    int zeroValue = cp.getInteger("Special Values/zero");
    EXPECT_EQ(0, zeroValue);

    int oneValue = cp.getInteger("Special Values/one");
    EXPECT_EQ(1, oneValue);

    int minusOneValue = cp.getInteger("Special Values/minus-one");
    EXPECT_EQ(-1, minusOneValue);

    bool booleanTrue = cp.getBoolean("Special Values/boolean-true");
    EXPECT_TRUE(booleanTrue);

    bool booleanFalse = cp.getBoolean("Special Values/boolean-false");
    EXPECT_FALSE(booleanFalse);

    bool yes = cp.getBoolean("Special Values/boolean-yes");
    EXPECT_TRUE(yes);

    bool no = cp.getBoolean("Special Values/boolean-no");
    EXPECT_FALSE(no);

    bool booleanOne = cp.getBoolean("Special Values/boolean-one");
    EXPECT_TRUE(booleanOne);

    bool booleanZero = cp.getBoolean("Special Values/boolean-zero");
    EXPECT_FALSE(booleanZero);

    ASSERT_THROW(cp.getBoolean("Special Values/boolean-integer"), ConfigurationException);

    ASSERT_THROW(cp.getInteger("Special Values/integer-float"), ConfigurationException);

}





TEST(ConfigurationParametersTest, testConversions)
{
    LOG_STARTING_OF_TEST

    ConfigurationParameters cp = ConfigurationParameters("../testData/input_ConfigurationParametersTest.yaml");

    // Can convert an integer value into a double
    double exposureTime = cp.getDouble("Observing/ExposureTime");
    EXPECT_DOUBLE_EQ(23.0, exposureTime);

    // Can convert an integer value into a string
    string exposureTimeString = cp.getString("Observing/ExposureTime");
    EXPECT_STREQ("23", exposureTimeString.c_str());

    // Can not convert a double value into an Integer
    ASSERT_ANY_THROW(cp.getInteger("Observing/LightCollectingArea"));

    // Conversion from Integer 23 to Boolean should also throw an exception
    ASSERT_ANY_THROW(cp.getBoolean("Observing/ExposureTime"));

}


TEST(ConfigurationParametersTest, testNonExistingKey)
{
    LOG_STARTING_OF_TEST

    ConfigurationParameters cp = ConfigurationParameters("../testData/input_ConfigurationParametersTest.yaml");

    ASSERT_THROW(string unknown = cp.getString("UnknownNode"), IllegalArgumentException);
    ASSERT_THROW(string unknown = cp.getString("Special Values/UnknownSubNode"), IllegalArgumentException);

    try
    {
        string unknown = cp.getString("Special Values/UnknownSubNode");
        FAIL() << "This should never fail";
    }
    catch(IllegalArgumentException ex)
    {
        string expected = "ConfigurationParameters: The sub-field \"UnknownSubNode\"";
        EXPECT_EQ(expected, string(ex.what()).substr(0, expected.size()));
    }

    ASSERT_THROW(cp.getBoolean("UnknownNode"), IllegalArgumentException);
    ASSERT_THROW(cp.getBoolean("Special Values/UnknownSubNode"), IllegalArgumentException);

    ASSERT_THROW(cp.getInteger("UnknownNode"), IllegalArgumentException);
    ASSERT_THROW(cp.getInteger("Special Values/UnknownSubNode"), IllegalArgumentException);

    ASSERT_THROW(cp.getDouble("UnknownNode"), IllegalArgumentException);
    ASSERT_THROW(cp.getDouble("Special Values/UnknownSubNode"), IllegalArgumentException);

    ASSERT_THROW(string unknown = cp.getAbsoluteFilename("UnknownNode"), IllegalArgumentException);
    ASSERT_THROW(string unknown = cp.getAbsoluteFilename("Special Values/UnknownSubNode"), IllegalArgumentException);

}





// This test checks that new fields can be set or added to a node.
// 
// * Set a new fields to a value
// * Set the same field to another value
// 
// A log message will be issued as a warning that the field is overwritten.

TEST(ConfigurationParametersTest, testSetNode)
{
    LOG_STARTING_OF_TEST

    ConfigurationParameters cp = ConfigurationParameters();

    string value;

    cp.setParameter("One New Key", "AFDEEBACFD");
    ASSERT_NO_THROW(value = cp.getString("One New Key"));
    ASSERT_STREQ("AFDEEBACFD", value.c_str());

    cp.setParameter("One New Key", "DFCABEEDFA");
    ASSERT_NO_THROW(value = cp.getString("One New Key"));
    ASSERT_STREQ("DFCABEEDFA", value.c_str());

}






// This test checks that new fields can be set or added to a map.
// The map doesn't exist initially.
// 
// * Set a new fields in a non-existing map to a value
// * Set the same field in the now existing map to another value
// * Set another field in the map to a value
// 
// A log message will be issued as a warning that the field is overwritten.

TEST(ConfigurationParametersTest, testSetSubNode)
{
    LOG_STARTING_OF_TEST

    ConfigurationParameters cp = ConfigurationParameters();

    string value;

    cp.setParameter("New Keys/Another New Key", "AFDEEBACFD");
    ASSERT_NO_THROW(value = cp.getString("New Keys/Another New Key"));
    ASSERT_STREQ("AFDEEBACFD", value.c_str());

    cp.setParameter("New Keys/Another New Key", "DFCABEEDFA");
    ASSERT_NO_THROW(value = cp.getString("New Keys/Another New Key"));
    ASSERT_STREQ("DFCABEEDFA", value.c_str());

    cp.setParameter("New Keys/Yet Another New Key", "FDBACEFDCB");
    ASSERT_NO_THROW(value = cp.getString("New Keys/Yet Another New Key"));
    ASSERT_STREQ("FDBACEFDCB", value.c_str());

}




TEST(ConfigurationParametersTest, Sequences)
{
    LOG_STARTING_OF_TEST

    ConfigurationParameters cp = ConfigurationParameters("../testData/input_ConfigurationParametersTest.yaml");

    vector<double> values = cp.getDoubleVector("Sequences/Polynomial/Coefficients");

    double expected[] = {-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5};

    for(std::vector<double>::size_type idx = 0; idx < values.size(); idx++) 
    {
        EXPECT_DOUBLE_EQ(expected[idx], values[idx]);
    }

}






