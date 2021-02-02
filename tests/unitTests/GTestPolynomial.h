#include "gtest/gtest.h"

#include "gtest_definitions.h"

#include "Polynomial1D.h"
#include "Polynomial2D.h"

using namespace std;



TEST(Polynomial1DTest, defaultConstructor)
{
    LOG_STARTING_OF_TEST

    Polynomial1D p = Polynomial1D();

    EXPECT_NEAR(1.0, p(42), 0.00001);    
    EXPECT_NEAR(1.0, p(-42), 0.00001);    
}

TEST(Polynomial1DTest, EvaluationDeg0)
{
    LOG_STARTING_OF_TEST

    vector<double> coeff {42.30};

    Polynomial1D p = Polynomial1D(0, coeff);

    EXPECT_NEAR(42.30, p(23.0), 0.00001);    
}



TEST(Polynomial1DTest, EvaluationDeg1)
{
    LOG_STARTING_OF_TEST

    vector<double> coeff = {42.30, 16.7};

    Polynomial1D p = Polynomial1D(1, coeff);

    EXPECT_NEAR(  42.30, p(  0.0), 0.00001);
    EXPECT_NEAR( 426.40, p( 23.0), 0.00001);
    EXPECT_NEAR(-341.80, p(-23.0), 0.00001);
}




TEST(Polynomial1DTest, EvaluationDeg2)
{
    LOG_STARTING_OF_TEST

    vector<double> coeff = {42.30, 16.7, 0.02554};

    Polynomial1D p = Polynomial1D(2, coeff);

    EXPECT_NEAR(  42.30000, p(  0.0), 0.00001);
    EXPECT_NEAR( 439.91066, p( 23.0), 0.00001);
    EXPECT_NEAR(-328.28934, p(-23.0), 0.00001);
}




TEST(Polynomial2DTest, defaultConstructor)
{
    LOG_STARTING_OF_TEST

    Polynomial2D p = Polynomial2D();

    EXPECT_NEAR(1.0, p( 1.0,  1.0), 0.00001);
    EXPECT_NEAR(1.0, p(42.0, 23.0), 0.00001);

}

TEST(Polynomial2DTest, EvaluationDeg1)
{

    LOG_STARTING_OF_TEST

    vector<double> coeff = {0.5, 3.0, 2.5};

    Polynomial2D p = Polynomial2D(1, coeff);

    EXPECT_NEAR( 6.0, p(1, 1),   0.00001);
    EXPECT_NEAR(55.5, p(10, 10),   0.00001);

}







TEST(Polynomial2DTest, EvaluationDeg2)
{

    LOG_STARTING_OF_TEST

    // These coeff where taken from the 2nd degree 2D polynomial that 
    // was fit to the field distortion data of the PLATO Camera optics.
    // The fit was done using astropy and the results where obtained
    // from the same fitted 2D polynomial.

    vector<double> coeff = {
         0.609105776591, 
         4.07897153594, 
         0.0262060627497, 
        -0.0922743981284,
         0.00485527246587, 
         0.00920092300214
    };

    Polynomial2D p = Polynomial2D(2, coeff);

    EXPECT_NEAR( 4.6360651726, p(1, 1),   0.00001);
    EXPECT_NEAR( 8.7435490851, p(2, 2),   0.00001);
    EXPECT_NEAR(21.5491479211, p(5, 5),   0.00001);
    EXPECT_NEAR(44.5023029764, p(10, 10), 0.00001);
    EXPECT_NEAR(85.4146259196, p(18, 18), 0.00001);



}






TEST(Polynomial2DTest, EvaluationDeg3)
{

    LOG_STARTING_OF_TEST

    // These coeff where taken from the 3rd degree 2D polynomial that 
    // was fit to the field distortion data of the PLATO Camera optics.
    // The fit was done using astropy and the results where obtained
    // from the same fitted 2D polynomial.

    vector<double> coeff = { 
        -0.0805112828134,
         4.34991909425, 
        -0.00540040725843, 
         0.00112767910467,
         0.0236240922537, 
        -0.00237537065657, 
         7.24631958452e-05,
        -0.00366561903111,
         0.000586015203847,
         0.000128792686889 
    };

    Polynomial2D p = Polynomial2D(3, coeff);

    EXPECT_NEAR( 4.2835054569, p(1, 1),   0.00001);
    EXPECT_NEAR( 8.6361291039, p(2, 2),   0.00001);
    EXPECT_NEAR(21.7405384999, p(5, 5),   0.00001);
    EXPECT_NEAR(44.4257310788, p(10, 10), 0.00001);
    EXPECT_NEAR(86.1042429791, p(18, 18), 0.00001);

}




TEST(Polynomial2DTest, EvaluationDeg4)
{

    LOG_STARTING_OF_TEST

    vector<double> coeff(15);

    ASSERT_THROW(Polynomial2D p = Polynomial2D(2, coeff), IllegalArgumentException);
    ASSERT_THROW(Polynomial2D p = Polynomial2D(4, coeff), UnsupportedException);

}


