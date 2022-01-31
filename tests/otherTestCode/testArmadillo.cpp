
// compile with:
//  clang++ -o testArmadillo testArmadillo.cpp -I ../../dependencies/Installs/armadillo-6.500.4/include  -L../../dependencies/Installs/armadillo-6.500.4/lib -larmadillo -stdlib=libc++ -std=c++14 


#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

int main()
{
    Cube<double> A(3, 5, 2, fill::randu);

    Mat<double> x = A.slice(1);

    cout << "A: " << endl << A << endl;
    cout << "A.slice(1): " << endl << x << endl;

}
