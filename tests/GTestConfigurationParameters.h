//
//  GTestInputParameters.h
//  New PLATO Simulator
//
//  Created by Rik Huygen on 26/11/15.
//  Copyright © 2015 KU Leuven. All rights reserved.
//

#include <string>
#include <vector>
#include <list>

#include "gtest/gtest.h"

#include "ConfigurationParameters.h"
#include "Exceptions.h"


TEST(ConfigurationParametersTest, Constructor) {

    ASSERT_THROW(ConfigurationParameters ip = ConfigurationParameters("input.json"), IOException);

    ASSERT_THROW(ConfigurationParameters ip = ConfigurationParameters("input.yaml"), IOException);

}

TEST(ConfigurationParametersTest, readGlobalValuesJSon) {
    ConfigurationParameters ip = ConfigurationParameters("../testData/input.json");

    std::string description = ip.getString("Description");
    EXPECT_STREQ(description.c_str(), "JSON Input File for 3rd Generation PLATO Simulator");

    std::string author = ip.getString("Author");
    EXPECT_STREQ(author.c_str(), "Rik Huygen");

    int exposureTime = ip.getInteger("Observing/ExposureTime");
    EXPECT_EQ(23, exposureTime);

}

TEST(ConfigurationParametersTest, readGlobalValuesYaml) {
    ConfigurationParameters ip = ConfigurationParameters("../testData/input.yaml");

    std::string description = ip.getString("Description");
    EXPECT_STREQ(description.c_str(), "YAML Input File for 3rd Generation PLATO Simulator");

    std::string author = ip.getString("Author");
    EXPECT_STREQ(author.c_str(), "Rik Huygen");

    int exposureTime = ip.getInteger("Observing/ExposureTime");
    EXPECT_EQ(23, exposureTime);

}

TEST(ConfigurationParametersTest, readGeneralValues) {
    ConfigurationParameters ip = ConfigurationParameters("../testData/input.json");

//    ConfigurationParameters general = ip.getGroup("General");

//    std::string dir2 = general.getString("ProjectLocation");
//    EXPECT_STREQ(dir1.c_str(), dir2.c_str());

//    std::list<std::string> authors = general.getValues<std::string>("Author");
//    EXPECT_STREQ(authors.front().c_str(), "Joris De Ridder");

    ASSERT_NO_THROW(ip.getString("Description"));
//    general.getString("Description");
//    ASSERT_THROW(general.getString("Description"), IOException);
}
