#ifndef MAPPEDDISTORTION_H
#define MAPPEDDISTORTION_H

#include <vector>

class MappedDistortion
{
public:
    MappedDistortion(const std::vector<double> &x1,
                     const std::vector<double> &x2,
                     const std::vector<double> &z1,
                     const std::vector<double> &z2,
                     double focalLength);                                 // The constructor that we will use when we try to estimate the Wang distortion model.

    MappedDistortion(const std::vector<double> &x1,
                     const std::vector<double> &x2,
                     const std::vector<double> &z1,
                     const std::vector<double> &z2);                      // The constructor that we will use when we try to estimate the Polynomial distortion model.

    std::vector<double> getParameters();                                  // Returns the coefficients for the Wang distortion model.
    std::vector<double> getParametersX();                                 // Returns the coefficients for the polynomial distortion model for the x coordinates.
    std::vector<double> getParametersY();                                 // Returns the coefficients for the polynomial distortion model for the y coordinates.

    // These functions are only used to help debuggin
    bool resultIsSensible();
    std::pair<double, double> getRMS();

private:
    // Where we store the input data
    std::vector<double> x1;
    std::vector<double> x2;
    std::vector<double> y1;
    std::vector<double> y2;


    bool isPolynomial;                         // Is false if we estimate the Wang model.

    // These values are used fill matrix A and B for the Wang model.
    double focalLength;
    std::vector<double> r;
    std::vector<double> cos;
    std::vector<double> sin;
    std::vector<double> deltaX;
    std::vector<double> deltaY;

    // The matrices that are used to solve the equation A*x = B.
    // For polynomial model we use have two B matrices B and By.
    std::vector<std::vector<double>> A;
    std::vector<double> B;
    std::vector<double> By;

    // L(ower diagonal) matrix and U(pper diagonal) matrix so that LU = A.
    std::vector<std::vector<double>> L;
    std::vector<std::vector<double>> U;


    void addMatrixToA(std::vector<std::vector<double>> &matrixToAdd);
    void addVectorToB(std::vector<double> &vectorToAdd);


    void constructMatrix();                                         // Initializes A and B.
    void constructMatrixToFindPolynomial();                         // Initializes A, B and By.


    void luDecomposition();                                                // Decomposes the matrix A in L and U so that LU = A;
    void solveL(std::vector<double>& B, std::vector<double>& T);           // Solves L*T=B;
    void solveU(std::vector<double>& T, std::vector<double>& output);      // Solves U*output = T

    // Functions used for debugging purposes
    double applyDistortion(double x, double y, std::vector<double> &coefficients);
    void printMatrix(std::vector<std::vector<double>> &a);
    std::vector<std::vector<double>> multiplyMatrices(std::vector<std::vector<double>> &a,
                                                      std::vector<std::vector<double>> &b);
    std::vector<double> dotProduct(std::vector<std::vector<double>> A, std::vector<double> B);
    bool matricesAreEqual(std::vector<std::vector<double>> &a,
                          std::vector<std::vector<double>> &b);
    bool vectorsAreEqual(std::vector<double> &a,
                         std::vector<double> &b);
};
#endif
