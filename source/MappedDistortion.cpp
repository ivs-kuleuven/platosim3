#include "MappedDistortion.h"
#include <vector>
#include <iostream>
#include <cmath>











/**
 * \brief Estimate the closest analytical model that approximates the distortion
 *        given by the input parameters.
 *
 * \param x1  Vector that contains the undistorted x focal plane coordinates.
 *
 * \param x2  Vector that contains the undistorted y focal plane coordinates.
 *
 * \param z1  Vector that contains the distorted x focal plane coordinates.
 *
 * \param z2 Vector that contains the distorted y focal plane coordinates.
 *
 * \param focLength double that contains the value of the focal length.
 *
 * \note  This procedure finds the 7 coeficients that describe the analytic
 *        distortion model best. It does this by solving a set of linear equations
 *        (A*x = B) by decomposing the linear matrix (A) into an upper diagonal (U)
 *         and lower diagonal (L) matrix. (eq. LU decomposition). [ A = LU ]
 *
 */
MappedDistortion::MappedDistortion(const std::vector<double> &x1,
                                   const std::vector<double> &x2,
                                   const std::vector<double> &z1,
                                   const std::vector<double> &z2,
                                   double focLength)
{
    unsigned int length = x1.size();
    focalLength = focLength;

    // Initialize the upper and lower diagonal matrix
    L.resize(7, std::vector<double>(7, 0.0));
    U.resize(7, std::vector<double>(7, 0.0));

    // Initialize the vectors that will be used to store the coeficients and
    // the intermidiate results (those that solve L x = B).
    intermediate.resize(7, 0);
    output.resize(7,0);

    // Calculate the values that are needed in construct the vector A and B.
    for (unsigned int i=0; i<length; i++)
    {
        double r1 = std::sqrt(x1[i] * x1[i] + x2[i] * x2[i]) / focalLength;
        r.push_back(r1);
        if (r1 == 0.0) {
          cos.push_back(0);
          sin.push_back(1);
        } else
        {
            cos.push_back(x1[i] / (r1*focalLength));
            sin.push_back(x2[i] / (r1*focalLength));
        }
        deltaX.push_back(z1[i] - x1[i]);
        deltaY.push_back(z2[i] - x2[i]);
    }

    // Initialize dimensions of the matrix A and vector B.
    for (int j = 0; j < 7; j++)
    {
        std::vector<double> a(7, 0);
        A.push_back(a);
        B.push_back(0);
    }

    // We fill the matrix A and vector B.
    constructMatrix();
}










/**
 * \brief Fills the matrix A and vector B with their respective values,
 *        so that the system that we are trying to solve is A*x = B.
 *        (where we try and solve for vector x)
 */
void MappedDistortion::constructMatrix()
{


    for (int i = 0; i < r.size(); i++)
    {
        // Add to the matrix A.
        std::vector<std::vector<double>> a = {
            {std::pow(r[i], 6), std::pow(r[i], 8), std::pow(r[i], 10),
             std::pow(r[i], 5) * cos[i], std::pow(r[i], 5) * sin[i],
             std::pow(r[i], 5) * cos[i], std::pow(r[i], 5) * sin[i]},
            {std::pow(r[i], 8), std::pow(r[i], 10), std::pow(r[i], 12),
             std::pow(r[i], 7) * cos[i], std::pow(r[i], 7) * sin[i],
             std::pow(r[i], 7) * cos[i], std::pow(r[i], 7) * sin[i]},
            {std::pow(r[i], 10), std::pow(r[i], 12), std::pow(r[i], 14),
             std::pow(r[i], 9) * cos[i], std::pow(r[i], 9) * sin[i],
             std::pow(r[i], 9) * cos[i], std::pow(r[i], 9) * sin[i]},
            {std::pow(r[i], 5) * cos[i], std::pow(r[i], 7) * cos[i],
             std::pow(r[i], 9) * cos[i], std::pow(r[i], 4), 0,
             std::pow(r[i], 4) * std::pow(cos[i], 2),
             std::pow(r[i], 4) * sin[i] * cos[i]},
            {std::pow(r[i], 5) * sin[i], std::pow(r[i], 7) * sin[i],
             std::pow(r[i], 9) * sin[i], 0, std::pow(r[i], 4),
             std::pow(r[i], 4) * cos[i] * sin[i],
             std::pow(r[i], 4) * std::pow(sin[i], 2)},
            {std::pow(r[i], 5) * cos[i], std::pow(r[i], 7) * cos[i],
             std::pow(r[i], 9) * cos[i],
             std::pow(r[i], 4) * std::pow(cos[i], 2),
             std::pow(r[i], 4) * sin[i] * cos[i],
             std::pow(r[i], 4) * std::pow(cos[i], 2),
             std::pow(r[i], 4) * sin[i] * cos[i]},
            {std::pow(r[i], 5) * sin[i], std::pow(r[i], 7) * sin[i],
             std::pow(r[i], 9) * sin[i],
             std::pow(r[i], 4) * sin[i] * cos[i],
             std::pow(r[i], 4) * std::pow(sin[i], 2),
             std::pow(r[i], 4) * sin[i] * cos[i],
             std::pow(r[i], 4) * std::pow(sin[i], 2)}
        };
        addMatrixToA(a);

        // Add to the vector B.
        std::vector<double> b = {
            (deltaX[i] * cos[i] + deltaY[i] * sin[i])*std::pow(r[i], 3) / focalLength,
            (deltaX[i] * cos[i] + deltaY[i] * sin[i])*std::pow(r[i], 5) / focalLength,
            (deltaX[i] * cos[i] + deltaY[i] * sin[i])*std::pow(r[i], 7) / focalLength,
            deltaX[i] * std::pow(r[i], 2) / focalLength,
            deltaY[i] * std::pow(r[i], 2) / focalLength,
            (deltaX[i] * cos[i] + deltaY[i] * sin[i])*cos[i]*std::pow(r[i], 2) / focalLength,
            (deltaX[i] * cos[i] + deltaY[i] * sin[i])*sin[i]*std::pow(r[i], 2) / focalLength};

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
    // We decompose the matrix A so that A = L*U.
    luDecomposition();

    // Solve the linear equations L * intermediate = B, where L is a lower diagonal
    // matrix.
    solveL();

    // Finally we solve the linear equatins U * output = intermediate, where U is an
    // upper diagonal matrix.
    solveU();

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
    for (int i=0; i<7; i++)
    {
        L[i][i] = 1.0;
        for (int j=0; j<7; j++)
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
void MappedDistortion::solveL()
{
    intermediate[0] = B[0];
    for (int i = 1; i < 7; i++)
    {
        intermediate[i] = B[i];
        for (int j = i-1; j >= 0; j--)
        {
            intermediate[i] -= L[i][j]*intermediate[j];
        }
    }

}










/**
 * \brief Solves the equation L*x = intermediate.
 *
 * \note We use the fact that U is an upper diagonal matrix.
 *
 */
void MappedDistortion::solveU()
{

    for (int i=6; i>=0; i--)
    {
        output[i] = intermediate[i];
        for (int j = 6; j > i; j--)
        {
            output[i] -= U[i][j]*output[j];

        }
        if (U[i][i] == 0 && intermediate[i] != 0) {
            throw std::invalid_argument(
                "We can not solve these equations!");
        } else if (intermediate[i] == 0) {
            output[0] = 0;
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
    for (int i=0; i<7; i++) {
        if (abs(a[i] - b[i]) > 0.1) {
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
        return false;
    }

    std::cout << "Test2: L * intermediate= B" << std::endl;
    std::vector<double> q = dotProduct(L, intermediate);
    std::cout << vectorsAreEqual(q, B) << std::endl;
    std::cout << " " << std::endl;
    if (!vectorsAreEqual(q, B)) {
        return false;
    }

    std::cout <<"Test3: U * output = output" << std::endl;
    std::vector<double> q3 = dotProduct(U, output);
    std::cout << vectorsAreEqual(q3, intermediate) << std::endl;;
    std::cout << " " << std::endl;
    if (!vectorsAreEqual(q3, intermediate)) {
        return false;
    }

    std::cout << "Test4: A * output2 = B" << std::endl;
    std::vector<double> q4 = dotProduct(A, output);
    std::cout << vectorsAreEqual(q4, B);
    std::cout << " " << std::endl;
    std::cout << " " << std::endl;
    if (!vectorsAreEqual(q4, B)) {
        return false;
    }

    std::cout << "Test5: Coordinates match those that we expect it to match" << std::endl;
    std::vector<double> theo_coef = {0.32419,  0.0232909,  0.407979,  0.00022463,  0.000217599,  0.000381958,  0.000963902};
    std::cout << vectorsAreEqual(theo_coef, output) << std::endl;
    if (!vectorsAreEqual(theo_coef, output)) {
        return false;
    }
    return true;
}

