
#ifndef CONVOLVER_H
#define CONVOLVER_H

#include <string>
#include "armadillo"
#include "fftw3.h"
#include "Logger.h"
#include "Exceptions.h"



// Convolver works with single-precision (float) matrices. The means that FFTW should be compiled
// for single-precision Fourier transforms (--enable-float), and be linked against "-l fftw3f"
// rather than "-l fftw".


typedef arma::Mat<float> matrix;


class Convolver
{
    public:

        Convolver();
        ~Convolver();
        void initialise(int Nrows, int Ncols, matrix &kernel);
        void convolve(matrix &in, matrix &out, double zeroThreshold = 0.0);

    protected:

        void free();
        void createExtendedMatrices(int Nrows, int Ncols, matrix &kernel);
        void createWrapAroundKernel(matrix &kernel);

    private:

        bool isInitialised;            // true if initialise() has been called, false otherwise

        int NrowsInOut;                // Number of rows on the input and output arrays
        int NcolsInOut;                // Number of columns on the input and output arrays

        int NrowsExtended;             // Number of rows on the zero-padded input and kernel arrays
        int NcolsExtended;             // Number of columns on the zero-padded input and kernel arrays

        matrix extendedKern;           // The zero-padded kernel in wrap-around format, column-major
        float *extendedOut;            // 2D row-major convolution result of convolving zero-padded arrays

        float *extendedIn;             // 2D zero-padded input array, row-major
        float *copyExtendedKern;       // 2D row-major copy of the column-major extendedKern

        fftwf_complex *fourierIn;      // Complex fourier transform of zero-padded input array
        fftwf_complex *fourierKern;    // Complex fourier transform of zero-padded and wrapped-around kernel
        fftwf_complex *fourierOut;     // Complex fourier transform of the convolution

        fftwf_plan forwardPlanIn;      // Fourier transformation of the zero-padded input array
        fftwf_plan forwardPlanKern;    // Fourier transformation of the zero-padded and wrapped-around kernel
        fftwf_plan backwardPlanOut;    // Inverse Fourier transformation of fourierOut
};


#endif
