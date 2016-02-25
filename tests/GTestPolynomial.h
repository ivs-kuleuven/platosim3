#include "gtest/gtest.h"

#include "gtest_definitions.h"

#include "Polynomial2D.h"

using namespace std;



TEST(PolynomialTest, EvaluationDeg1)
{

    LOG_STARTING_OF_TEST

    double coeff_1[] = {0.5, 3.0, 2.5};

    Polynomial2D p = Polynomial2D(1, coeff_1);

    EXPECT_NEAR( 6.0, p.evaluateAt(1, 1),   0.00001);
    EXPECT_NEAR(55.5, p.evaluateAt(10, 10),   0.00001);

}







TEST(PolynomialTest, EvaluationDeg2)
{

    LOG_STARTING_OF_TEST

    // These coeff where taken from the 2nd degree 2D polynomial that 
    // was fit to the field distortion data of the PLATO Camera optics.
    // The fit was done using astropy and the results where obtained
    // from the same fitted 2D polynomial.

    double coeff_2[] = {
         0.609105776591, 
         4.07897153594, 
         0.0262060627497, 
        -0.0922743981284,
         0.00485527246587, 
         0.00920092300214
    };

    Polynomial2D p = Polynomial2D(2, coeff_2);

    EXPECT_NEAR( 4.6360651726, p.evaluateAt(1, 1),   0.00001);
    EXPECT_NEAR( 8.7435490851, p.evaluateAt(2, 2),   0.00001);
    EXPECT_NEAR(21.5491479211, p.evaluateAt(5, 5),   0.00001);
    EXPECT_NEAR(44.5023029764, p.evaluateAt(10, 10), 0.00001);
    EXPECT_NEAR(85.4146259196, p.evaluateAt(18, 18), 0.00001);



}






TEST(PolynomialTest, EvaluationDeg3)
{

    LOG_STARTING_OF_TEST

    // These coeff where taken from the 3rd degree 2D polynomial that 
    // was fit to the field distortion data of the PLATO Camera optics.
    // The fit was done using astropy and the results where obtained
    // from the same fitted 2D polynomial.

    double coeff_3[] = { -0.0805112828134,
        4.34991909425, -0.00540040725843, 0.00112767910467,
        0.0236240922537, -0.00237537065657, 7.24631958452e-05,
        -0.00366561903111, 0.000586015203847, 0.000128792686889 
    };

    Polynomial2D p = Polynomial2D(3, coeff_3);

    EXPECT_NEAR( 4.2835054569, p.evaluateAt(1, 1),   0.00001);
    EXPECT_NEAR( 8.6361291039, p.evaluateAt(2, 2),   0.00001);
    EXPECT_NEAR(21.7405384999, p.evaluateAt(5, 5),   0.00001);
    EXPECT_NEAR(44.4257310788, p.evaluateAt(10, 10), 0.00001);
    EXPECT_NEAR(86.1042429791, p.evaluateAt(18, 18), 0.00001);

}




TEST(PolynomialTest, EvaluationDeg4)
{

    LOG_STARTING_OF_TEST

    double coeff_4[15];

    ASSERT_THROW(Polynomial2D p = Polynomial2D(4, coeff_4), UnsupportedException);

}


