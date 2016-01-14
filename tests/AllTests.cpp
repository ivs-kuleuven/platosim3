//
//  AllTests.cpp
//  New PLATO Simulator
//
//  Created by Rik Huygen on 27/11/15.
//  Copyright © 2015 KU Leuven. All rights reserved.
//


// Include all the tests here

#include "GTestConfigurationParameters.h"


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

