#include "gtest/gtest.h"

#include "gtest_definitions.h"

#include "Polynomial2D.h"

using namespace std;


TEST(PolynomialTest, Evaluation)
{

    LOG_STARTING_OF_TEST

    double coeff_3[] = { -0.0805112828134,
        4.34991909425, -0.00540040725843, 0.00112767910467,
        0.0236240922537, -0.00237537065657, 7.24631958452e-05,
        -0.00366561903111, 0.000586015203847, 0.000128792686889 
    };

    Polynomial2D p = Polynomial2D(3, coeff_3);

    EXPECT_NEAR(12.98885, p.evaluateAt(3, 3), 0.00001);
    EXPECT_NEAR(13.60543, p.evaluateAt(3.1415, 3.1415), 0.00001);

}
