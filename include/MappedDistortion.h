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
                     double focalLength);

    std::vector<double> getParameters();
    bool resultIsSensible();

private:
    double focalLength;

    std::vector<double> r;
    std::vector<double> cos;
    std::vector<double> sin;
    std::vector<double> deltaX;
    std::vector<double> deltaY;

    std::vector<std::vector<double>> L;
    std::vector<std::vector<double>> U;
    std::vector<double> intermediate;
    std::vector<double> output;

    std::vector<std::vector<double>> A;
    std::vector<double> B;


    void addMatrixToA(std::vector<std::vector<double>> &matrixToAdd);
    void addVectorToB(std::vector<double> &vectorToAdd);

    void constructMatrix();
    void luDecomposition();
    void solveL();
    void solveU();

    // Functions used for debugging purposes
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
