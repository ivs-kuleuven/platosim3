#include "MappedDistortion.h"
#include <H5Ppublic.h>
#include <vector>
#include <iostream>
#include <cmath>











/**
 * \brief Estimate the closest analytical Wang model that approximates the distortion
 *        given by the input parameters.
 *
 * \param x1  Vector that contains the undistorted x focal plane coordinates.
 *
 * \param x2  Vector that contains the undistorted y focal plane coordinates.
 *
 * \param y1  Vector that contains the distorted x focal plane coordinates.
 *
 * \param y2 Vector that contains the distorted y focal plane coordinates.
 *
 * \param focLength double that contains the value of the focal length.
 *
 * \note  This procedure finds the 7 coeficients that describe the analytic Wang
 *        distortion model best. It does this by solving a set of linear equations
 *        (A*x = B) by decomposing the linear matrix (A) into an upper diagonal (U)
 *         and lower diagonal (L) matrix. (eq. LU decomposition). [ A = LU ]
 *
 */
MappedDistortion::MappedDistortion(const std::vector<double> &x1,
                                   const std::vector<double> &x2,
                                   const std::vector<double> &y1,
                                   const std::vector<double> &y2,
                                   double focLength)
:x1(x1), x2(x2), y1(y1), y2(y2)
{

    unsigned int length = x1.size();
    focalLength = focLength;
    isPolynomial = false;

    // Initialize the upper and lower diagonal matrix
    L.resize(7, std::vector<double>(7, 0.0));
    U.resize(7, std::vector<double>(7, 0.0));

    // Initialize dimensions of the matrix A and vector B.
    for (int j = 0; j < 7; j++)
    {
        std::vector<double> a(7, 0);
        A.push_back(a);
        B.push_back(0);
    }

    // We fill the matrix A and vector B.
    constructMatrix();

    // We decompose the matrix A so that A = L*U.
    luDecomposition();
}
















/**
 * \brief Estimate distortion by finding the closest 5th order 2D polynomial model
 *         that approximates the distortion.
 *
 * \param x1  Vector that contains the undistorted x focal plane coordinates.
 *
 * \param x2  Vector that contains the undistorted y focal plane coordinates.
 *
 * \param z1  Vector that contains the distorted x focal plane coordinates.
 *
 * \param z2 Vector that contains the distorted y focal plane coordinates.
 *
 * \note  This procedure finds the 36 coeficients that describe a polynomial
 *        distortion model best. It does this by solving a set of linear equations
 *        (A*x = B) by decomposing the linear matrix (A) into an upper diagonal (U)
 *         and lower diagonal (L) matrix. (eq. LU decomposition). [ A = LU ]
 *
 */
MappedDistortion::MappedDistortion(const std::vector<double> &x1,
                                   const std::vector<double> &x2,
                                   const std::vector<double> &y1,
                                   const std::vector<double> &y2)
:x1(x1), x2(x2), y1(y1), y2(y2)
{
    unsigned int length = x1.size();
    bool isPolynomial = true;

    // Initialize the upper and lower diagonal matrix
    L.resize(36, std::vector<double>(36, 0.0));
    U.resize(36, std::vector<double>(36, 0.0));

    // Initialize dimensions of the matrix A and vector B.
    for (int j = 0; j < 36; j++)
    {
        std::vector<double> a(36, 0);
        A.push_back(a);
        B.push_back(0);
        By.push_back(0);
    }

    // We fill the matrix A and vector B.
    constructMatrixToFindPolynomial();


    // We decompose the matrix A so that A = L*U.
    luDecomposition();
}














/**
 * \brief Fills the matrix A and vector B and By with their respective values,
 *        so that the system that we are trying to solve is A*x = B for distortion
 *        in the x-direction and A*x = By for distortion in the y-direction.
 *        (where we try and solve for vector x)
 */
void MappedDistortion::constructMatrixToFindPolynomial()
{

    for (int t = 0; t < x1.size(); t++)
    {
        // Add to the matrix A.
        for (int i_col = 0; i_col < 6; i_col++)
        {
            for (int j_col = 0; j_col < 6; j_col++)
            {
                for (int i_row = 0; i_row < 6; i_row++)
                {
                    for (int j_row = 0; j_row < 6; j_row++)
                    {
                        A[i_col+6*j_col][i_row+6*j_row] += std::pow( x1[t], i_col + i_row) * std::pow( x2[t], j_col + j_row);
                    }
                }
                B[i_col+6*j_col]  += y1[t]*std::pow(x1[t],i_col)*std::pow(x2[t],j_col);
                By[i_col+6*j_col] += y2[t]*std::pow(x1[t],i_col)*std::pow(x2[t],j_col);
            }
        }
    }
}








/**
 * \brief Fills the matrix A and vector B with their respective values,
 *        so that the system that we are trying to solve is A*x = B.
 *        (where we try and solve for vector x)
 */
void MappedDistortion::constructMatrix()
{

    for (unsigned int i = 0; i < x1.size(); i++)
    {
        // Initialize values that help with constructing the vector A and B.
        double r = std::sqrt(x1[i] * x1[i] + x2[i] * x2[i]) / focalLength;
        double cos, sin;
        double deltaX, deltaY;
        if (r == 0.0) {
          cos = 0;
          sin = 1;
        } else
        {
            cos = x1[i] / (r*focalLength);
            sin = x2[i] / (r*focalLength);
        }
        deltaX = y1[i] - x1[i];
        deltaY = y2[i] - x2[i];

        // Add to the matrix A.
        std::vector<std::vector<double>> a = {
            {std::pow(r, 6), std::pow(r, 8), std::pow(r, 10),
             std::pow(r, 5) * cos, std::pow(r, 5) * sin,
             std::pow(r, 5) * cos, std::pow(r, 5) * sin},
            {std::pow(r, 8), std::pow(r, 10), std::pow(r, 12),
             std::pow(r, 7) * cos, std::pow(r, 7) * sin,
             std::pow(r, 7) * cos, std::pow(r, 7) * sin},
            {std::pow(r, 10), std::pow(r, 12), std::pow(r, 14),
             std::pow(r, 9) * cos, std::pow(r, 9) * sin,
             std::pow(r, 9) * cos, std::pow(r, 9) * sin},
            {std::pow(r, 5) * cos, std::pow(r, 7) * cos,
             std::pow(r, 9) * cos, std::pow(r, 4), 0,
             std::pow(r, 4) * std::pow(cos, 2),
             std::pow(r, 4) * sin * cos},
            {std::pow(r, 5) * sin, std::pow(r, 7) * sin,
             std::pow(r, 9) * sin, 0, std::pow(r, 4),
             std::pow(r, 4) * cos * sin,
             std::pow(r, 4) * std::pow(sin, 2)},
            {std::pow(r, 5) * cos, std::pow(r, 7) * cos,
             std::pow(r, 9) * cos,
             std::pow(r, 4) * std::pow(cos, 2),
             std::pow(r, 4) * sin * cos,
             std::pow(r, 4) * std::pow(cos, 2),
             std::pow(r, 4) * sin * cos},
            {std::pow(r, 5) * sin, std::pow(r, 7) * sin,
             std::pow(r, 9) * sin,
             std::pow(r, 4) * sin * cos,
             std::pow(r, 4) * std::pow(sin, 2),
             std::pow(r, 4) * sin * cos,
             std::pow(r, 4) * std::pow(sin, 2)}
        };
        addMatrixToA(a);

        // Add to the vector B.
        std::vector<double> b =
            {
            (deltaX*cos + deltaY*sin)*std::pow(r, 3) / focalLength,
            (deltaX*cos + deltaY*sin)*std::pow(r, 5) / focalLength,
            (deltaX*cos + deltaY*sin)*std::pow(r, 7) / focalLength,
            deltaX*std::pow(r, 2) / focalLength,
            deltaY*std::pow(r, 2) / focalLength,
            (deltaX*cos + deltaY*sin)*cos*std::pow(r,2) / focalLength,
            (deltaX*cos + deltaY*sin)*sin*std::pow(r,2) / focalLength
            };

        addVectorToB(b);
    }
}









/**
 * \brief This method solved the linear set of equations given by A*x = B.
 *        The solution to this problem contains the coeficients with which the
 *        analytic distortion most closely matches the input values.
 *
 * \note  The way we solve this problem is by decomposing the matrix A in the
 *        product of lower diagonal matrix (L) and upper diagonal matrix (U), so
 *        that A * x = L*U*x = B. It is then straightforward to solve L*X' = B,
 *        using the fact that L is lower diagonal, and then U*x = X' using the fact
 *        that U is an upper diagonal matrix.
 */
std::vector<double> MappedDistortion::getParameters()
{
    int size_b = B.size();
    std::vector<double> output;
    std::vector<double> intermediate;

    output.resize(size_b, 0);
    intermediate.resize(size_b, 0);

    // Solve the linear equations L * intermediate = B, where L is a lower diagonal
    // matrix.
    solveL(B, intermediate);

    // Finally we solve the linear equatins U * output = intermediate, where U is an
    // upper diagonal matrix.
    solveU(intermediate, output);

    return output;
}






/**
 * \brief This method solved the linear set of equations given by A*x = B.
 *        The solution to this problem contains the coeficients with which the
 *        analytic distortion most closely matches the input values.
 *
 * \note  This is an alias for the method getParameters so that the user is able
 *        to use getParametersX and getParametersY for the polynomial model and
 *        getParameters for the Wang model.
 */
std::vector<double> MappedDistortion::getParametersX()
{
    return getParameters();
}












/**
 * \brief This method solved the linear set of equations given by A*x = By.
 *        The solution to this problem contains the coeficients with which the
 *        analytic distortion most closely matches the input values.
 */
std::vector<double> MappedDistortion::getParametersY()
{
    int size_b = By.size();
    std::vector<double> output;
    std::vector<double> intermediate;

    output.resize(size_b, 0);
    intermediate.resize(size_b, 0);

    // Solve the linear equations L * intermediate = B, where L is a lower diagonal
    // matrix.
    solveL(By, intermediate);

    // Finally we solve the linear equatins U * output = intermediate, where U is an
    // upper diagonal matrix.
    solveU(intermediate, output);

    return output;
}











/**
 * \brief Adds the input vector to the vector B.
 *
 * \param vector    The vector to add to B.
 *
 */
void MappedDistortion::addVectorToB(std::vector<double> &vector)
{
    for (int i = 0; i < 7; i++)
    {
        B[i] += vector[i];
    }
}










/**
 * \brief Adds the input matrix to the matrix A.
 *
 * \param matrix      The matrix to add to A.
 */
void MappedDistortion::addMatrixToA(std::vector<std::vector<double>> &matrix)
{
    for (int i = 0; i < 7; i++)
    {
        for (int j = 0; j < 7; j++)
        {
            (A[i])[j] += (matrix[i])[j];
        }
    }
}











/**
 * \brief This follows the LU decomposition alghoritmn as described in
 *        [William H. Press et all., Numerical Recipes in C : the Art of Scientific
 *        Computing. Cambridge [Cambridgeshire] ; New York :Cambridge University
 *        Press, 1992]
 */
void MappedDistortion::luDecomposition()
{

    int size_L = L.size();
    for (int i=0; i<size_L; i++)
    {
        L[i][i] = 1.0;
        for (int j=0; j<size_L; j++)
        {
            if (i<=j)
            {
                double sum = 0.0;
                for (int k = 0; k < i; k++)
                {
                  sum += L[i][k] * U[k][j];
                }
                U[i][j] = A[i][j] - sum;
            } else
            {
                double sum = 0.0;
                for (int k = 0; k < j; k++)
                {
                    sum += L[i][k] * U[k][j];
                }
                L[i][j] = (A[i][j] - sum)/ U[j][j];
            }
        }
    }
}















/**
 * \brief Solves the equation L*x = B.
 *
 * \note We exploit the fact that L is a lower diagonal matrix.
 *
 */
void MappedDistortion::solveL(std::vector<double> &B1, std::vector<double> &intermediate)
{
    int intermediate_size = intermediate.size();
    intermediate[0] = B1[0];
    for (int i = 1; i < intermediate_size; i++)
    {
        intermediate[i] = B1[i];
        for (int j = i-1; j >= 0; j--)
        {
            intermediate[i] -= L[i][j]*intermediate[j];
        }
    }

}










/**
 * \brief Solves the equation U*x = intermediate.
 *
 * \note We use the fact that U is an upper diagonal matrix.
 *
 */
void MappedDistortion::solveU(std::vector<double>& interm, std::vector<double>& output)
{

    int size_output = output.size();
    double sommen = 0;
    for (int i=size_output-1; i>=0; i--)
    {
        output[i] = interm[i];
        for (int j = size_output-1; j > i; j--)
        {
            output[i] -= U[i][j]*output[j];
        }

        if (U[i][i] == 0 && output[i] != 0) {
            throw std::invalid_argument(
                "We can not solve these equations!");
        } else if (output[i] == 0) {
            output[i] = 0;
        }
        else {
            output[i] = output[i] / U[i][i];
        }
    }
}










/**
 * \brief print out the input matrix to the console.
 *
 * \param a the matrix that should be printed out.
 *
 * \Note This function is only used for debugging purpose.
 *
 */
void MappedDistortion::printMatrix(std::vector<std::vector<double>> &a)
{
    int dim = a.size();
    for (int i=0; i<dim; i++)
    {
        for (int j = 0; j < dim; j++) {
            std::cout << "   " << a[i][j];
        }
        std::cout << " " << std::endl;
    }
}










/**
 * \brief multiply a matrix with a vector and returns the corresponding vector.
 *
 * \param a The matrix with whith we multply the input vector.
 *
 * \param b The vector that is used to multiply the matrix a with.
 *
 * \Note This function is only used for debugging purpose.
 *
 */
std::vector<std::vector<double>> MappedDistortion::multiplyMatrices(std::vector<std::vector<double>> &a,
                                                                    std::vector<std::vector<double>> &b)
{
    int dim = a.size();
    std::vector<std::vector<double>> c (dim);

    for (int i=0; i<dim; i++)
    {
        std::vector<double> row (dim, 0);
        for (int j = 0; j < dim; j++)
        {
            double sum = 0;
            for (int k = 0; k < dim; k++)
            {
                sum += a[i][k]*b[k][j];
            }
            row[j] = sum;
        }
        c[i] = row;
    }

    return c;
}










/**
 * \brief Calculate the dot product between two vectors and return the result.
 *
 * \param A, B The two vectors that are used to calculate the dot product.
 *
 * \Note This function is only used for debugging purpose.
 *
 */
std::vector<double> MappedDistortion::dotProduct(std::vector<std::vector<double>> A,
                             std::vector<double> B)
{

    std::vector<double> C(B.size(), 0);
    for (int i = 0; i < B.size(); i++)
    {
        double sum = 0;
        for (int j = 0; j < B.size(); j++) {
            sum += A[i][j]*B[j];
        }
        C[i] = sum;
    }
    return C;
}











/**
 * \brief Check that the two input matrices are equal.
 *
 * \param a, b The two matrix that we compare.
 *
 * \Note This function is only used for debugging purpose.
 *
 */bool MappedDistortion::matricesAreEqual(std::vector<std::vector<double>> &a,
                                        std::vector<std::vector<double>> &b)
{
    for (int i = 0; i < 7; i++) {
        for (int j = 0; j < 7; j++) {
            if (abs(a[i][j] - b[i][j]) > 0.1) {
                return false;
            }
        }
    }
    return true;
}










/**
 * \brief Check that the two input vectors are equal.
 *
 * \param a, b The two vectos that we compare.
 *
 * \Note This function is only used for debugging purpose.
 */
bool MappedDistortion::vectorsAreEqual( std::vector<double> &a,
                                        std::vector<double> &b)
{
    int a_size = a.size();
    if (a_size != b.size()){ return false;}
    for (int i=0; i<a_size; i++) {
        if (2*abs(a[i] - b[i]) / (a[i] + b[i]) > 0.01) {
            return false;
        }
    }
    return true;
}










/**
 * \brief Does a check on the intermediate results to make sure what we do seems
 *        sensible.
 *
 * \Note This function is only used for debugging purpose.
 */
bool MappedDistortion::resultIsSensible()
{
    std::cout << "Test1: LU = A" << std::endl;
    std::vector<std::vector<double>> x = multiplyMatrices(L, U);
    std::cout << matricesAreEqual(x, A) << std::endl;
    std::cout << " " << std::endl;
    if (!matricesAreEqual(x, A)) {
        return false;}

    std::cout << "Test2: L * intermediate= B" << std::endl;

    if (!isPolynomial)
    {
        std::vector<double> output;
        output = getParameters();

        std::vector<double> q1 = dotProduct(A, output);
        std::cout << vectorsAreEqual(q1, B) << std::endl;
        std::cout << " " << std::endl;

        if (!vectorsAreEqual(q1, B))
        {
            return false;
        }
    }
    else
    {
        std::vector<double> output1;
        std::vector<double> output2;

        output1 = getParametersX();
        output2 = getParametersY();

        std::vector<double> q1 = dotProduct(A, output1);
        std::vector<double> q2 = dotProduct(A, output2);

        std::cout << vectorsAreEqual(q1, B) << std::endl;
        std::cout << vectorsAreEqual(q1, By) << std::endl;

        if (!vectorsAreEqual(q1, B) || !vectorsAreEqual(q2, By))
        {
            return false;
        }
    }
    return true;
}











/**
 * \brief Return the root mean square (RMS) difference between the polynomial distortion and the true
 *        distortion as given by the values in y1 and y2.
 *
 * \output RMSx : RMS in the x-direction
 *         RMSy : RMS in the y-direction
 *
 * \Note This function is only used for debugging purpose.
 *
 * \Note Only useable for the polynomial model
 */
std::pair<double, double> MappedDistortion::getRMS()
{
    std::vector<double> output1;
    std::vector<double> output2;

    double RMSx = 0;
    double RMSy = 0;
    output1 = getParametersX();
    output2 = getParametersY();

    for (int i = 0; i < x1.size(); i++)
    {
        double xE, yE;
        xE = applyDistortion(x1[i], x2[i], output1);
        yE = applyDistortion(x1[i], x2[i], output2);
        RMSx += std::pow((xE - y1[i]), 2);
        RMSy += std::pow((yE - y2[i]), 2);
    }
    return std::make_pair(std::sqrt(RMSx / x1.size()), std::sqrt(RMSy / x2.size()));

}


















/**
 * \brief Apply the polynomial distortion on the points x, y with polynomial
 *        coefficients.
 *
 * \input: x:  input x-coordinate
 *         y:  input y-coordinate
 *         coefficients: polynomial coefficient that describe the distortion
 *
 * \output output: distorted coordinate
 *
 * \Note This function is only used for debugging purpose.
 */
double MappedDistortion::applyDistortion(double x, double y, std::vector<double> &coefficients)
{
    double output = 0;
    for (int i=0; i<6; i++)
    {
        for (int j=0; j<6; j++)
        {
            output += coefficients[i+6*j]*std::pow(x, i)*std::pow(y, j);
        }
    }
    return output;
}

