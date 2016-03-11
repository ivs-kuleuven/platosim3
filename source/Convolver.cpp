
#include "Convolver.h"


/**
 * \brief  Default constructor
 * 
 */

Convolver::Convolver()
: isInitialised(false), extendedOut(NULL), extendedIn(NULL), copyExtendedKern(NULL), fourierIn(NULL), 
  fourierKern(NULL), fourierOut(NULL)
{
 
}







/**
 * \brief Destructor
 */

Convolver::~Convolver()
{
    if (isInitialised)
    {
        free();
    }
}








/**
 * \brief Free (deallocate) all arrays and FFTW plans that were previously allocated
 */

void Convolver::free()
{
    // The destruction should be done in the reverse order they were created
    // First the FFTW plans

    fftwf_destroy_plan(forwardPlanIn);
    fftwf_destroy_plan(forwardPlanKern);
    fftwf_destroy_plan(backwardPlanOut);

    // Then the extended real and complex arrays.

    fftwf_free(extendedIn);
    fftwf_free(copyExtendedKern);
    fftwf_free(extendedOut);

    fftwf_free(fourierIn);
    fftwf_free(fourierKern);
    fftwf_free(fourierOut);
}










/**
 * \brief    Initialise the Convolver
 * 
 * \details  Convolver convolves a (Nrows x Ncols) matrix with a kernel matrix, using FFTW.
 *           Care is taken that no wrap-around effects occur. The kernel matrix will be
 *           copied in wrap-around format and zero-padded. This needs to be done only once,
 *           hence the reason why the kernel is an argument of the constructor.
 * 
 * \param Nrows   Number of rows of the input matrix that will be convolved with the kernel
 * \param Ncols   Number of columns of the input matrix that will be convolved with the kernel
 * \param kernel  Kernel matrix to convolve with. Should be smaller in shape than input matrix.
 */

void Convolver::initialise(int Nrows, int Ncols, matrix &kernel)
{
    // First free any array allocations from a previous initialisation

    if (isInitialised)
    {
        free();
    }

    // Do sanity checks on the matrix dimensions

    if ((Nrows < 2) || (Ncols < 2))
    {
        string errorMessage = "Convolver: input matrix has a dimension < 2";
        Log.error(errorMessage);
        throw IllegalArgumentException(errorMessage);
    }

    if ((kernel.n_rows > Nrows) || (kernel.n_cols > Ncols))
    {
        string errorMessage = "Convolver: kernel dimensions (" + to_string(kernel.n_rows) + "," 
                            + to_string(kernel.n_cols) + ") should be smaller than input matrix dimensions "
                            + "(" + to_string(Nrows) + "," + to_string(Ncols) + ")";
        Log.error(errorMessage);
        throw IllegalArgumentException(errorMessage);
    }

    if ((kernel.n_rows < 2) || (kernel.n_cols < 2))
    {
        string errorMessage = "Convolver: kernel has a dimension < 2";
        Log.error(errorMessage);
        throw IllegalArgumentException(errorMessage);
    }

    // Keep the dimensions of the input and output matrices where
    //    Output = Input (x) kernel
    // where (x) denotes convolution.

    NrowsInOut = Nrows;
    NcolsInOut = Ncols;

    // To avoid boundary effects, the input and kernel matrices need to be zero-padded.
    // Make extended versions of the matrices that are large enough for zero-padding.

    createExtendedMatrices(Nrows, Ncols, kernel);

    // Create the FFTW plans. This can take some time, but needs to be done only once.
    // FFTW_ESTIMATE means that FFTW will try different algorithms, measure their performance,
    // and select the best one.
    // forwardPlanIn:   2D forward fourier transformation of the extended input matrix
    // forwardPlanKern: 2D forward fourier transformation of the extended kernel matrix
    // backwardPlanOut: 2D inverse fourier transformation of the extended output matrix

    forwardPlanIn   = fftwf_plan_dft_r2c_2d(NrowsExtended, NcolsExtended, extendedIn, fourierIn, FFTW_ESTIMATE);
    forwardPlanKern = fftwf_plan_dft_r2c_2d(NrowsExtended, NcolsExtended, copyExtendedKern, fourierKern, FFTW_ESTIMATE);
    backwardPlanOut = fftwf_plan_dft_c2r_2d(NrowsExtended, NcolsExtended, fourierOut, extendedOut, FFTW_ESTIMATE);

    // Compute the 2D fourier transform of the zero-padded wrap-around version of the kernel only once.

    fftwf_execute(forwardPlanKern);
 
    // Remember that all arrays and plans were initialised

    isInitialised = true;

}











/**
 * \brief Create larger versions of the input, kernel and output matrices that allow
 *        for zero-padding to avoid cyclic boundary effects.
 *
 * \param Nrows   Number of rows of the input matrix that will be convolved with the kernel
 * \param Ncols   Number of columns of the input matrix that will be convolved with the kernel
 * \param kernel  Kernel matrix to convolve with. Should be smaller in shape than input matrix.
 */


void Convolver::createExtendedMatrices(int Nrows, int Ncols, matrix &kernel)
{
    // The amount of zero-padding we need is half the size of the kernel, in either direction
    // (NrowsExtended, NcolsExtended) will be the shape of the extended input, the extended
    // kernel and the extended output matrices 

    NrowsExtended = Nrows + int(ceil(kernel.n_rows/2.0));
    NcolsExtended = Ncols + int(ceil(kernel.n_cols/2.0));

    // Create the extended zero-padded kernel
    // Copy the original kernel in wrap-around format in the extended kernel

    extendedKern = arma::zeros<matrix>(NrowsExtended, NcolsExtended);
    createWrapAroundKernel(kernel);

    // FFTW requires 2D matrices not in a [Nrows][Ncols] format, but in a [Nrows*Ncols]
    // format, in row-major. 
    // extendedIn:         2D zero-padded input array, row-major
    // copyExtendedKernel: 2D row-major copy of the column-major extendedKern.
    // extendedOut:        2D row-major convolution result of convolving zero-padded arrays

    extendedIn       = (float*) fftwf_malloc(sizeof(float) * NrowsExtended * NcolsExtended);
    copyExtendedKern = (float*) fftwf_malloc(sizeof(float) * NrowsExtended * NcolsExtended);
    extendedOut      = (float*) fftwf_malloc(sizeof(float) * NrowsExtended * NcolsExtended);

    // Copy the kernel from column-major to row-major format
    // Zero both the input and output arrays. 

    for (int i = 0; i < NrowsExtended; ++i)
    {
        for (int j = 0; j < NcolsExtended; ++j)
        {
            const int ij = i*NcolsExtended + j;
            copyExtendedKern[ij] = extendedKern(i,j);

            extendedIn[ij] = 0.0;
            extendedOut[ij] = 0.0;
        }
    }

    // fourierIn:   Complex fourier transform of zero-padded input array
    // fourierKern: Complex fourier transform of zero-padded and wrapped-around kernel
    // fourierOut:  Complex fourier transform of the convolution

    fourierIn   = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*NrowsExtended *(NcolsExtended/2+1));
    fourierKern = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*NrowsExtended *(NcolsExtended/2+1));
    fourierOut  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*NrowsExtended *(NcolsExtended/2+1));
}







/**
 * \brief Reshuffle the kernel in wrap-around format and zero-pad it. This is necessary
 *        to get rid of any cyclic boundary effects. 
 * 
 * \details The kernel matrix is divided into 4 quadrants, and reshuffled and zero-padded as:
 * 
 *                                    /  D 0 ... 0 C \
 *           /  A  |  B \             |  0    |    0 |
 *           |----------|     to      |  ... ---  ...|
 *           \  C  |  D /             |  0    |    0 |
 *                                    \  B 0 ... 0 A / 
 * 
 *           The quadrant boundaries are determined by the Nrows/2, and Ncols/2
 *           where the division is an integer division, and therefore rounded down for odd dimensions
 *           Care is taken that the method works for both even and odd dimensions.
 *  
 * \note  The extendedKern armadillo array was already set to zeros in initialise->createExtendedMatrices.
 *          
 * \param kernel:  Non-zero-padded kernel matrix
 */

void Convolver::createWrapAroundKernel(matrix &kernel)
{
    // Copy the lower right quarter of the original kernel into the upper left quarter of the extended kernel

    int rowStartOut = kernel.n_rows/2;            // Division is rounded downwards
    int colStartOut = kernel.n_cols/2;
    int rowEndOut   = kernel.n_rows-1;
    int colEndOut   = kernel.n_cols-1;
    int Nrows = rowEndOut - rowStartOut + 1;
    int Ncols = colEndOut - colStartOut + 1;

    int rowStartIn = 0;          
    int colStartIn = 0;
    int rowEndIn   = rowStartIn + Nrows-1;
    int colEndIn   = colStartIn + Ncols-1;

    extendedKern.submat(rowStartIn, colStartIn, rowEndIn, colEndIn) = kernel.submat(rowStartOut, colStartOut, rowEndOut, colEndOut);

    // Copy the upper right quarter of the original kernel into the lower left quarter of the extended kernel

    rowStartOut = 0;           
    colStartOut = kernel.n_cols/2;
    rowEndOut   = kernel.n_rows/2-1;
    colEndOut   = kernel.n_cols-1;
    Nrows = rowEndOut - rowStartOut + 1;
    Ncols = colEndOut - colStartOut + 1;    
    
    rowStartIn = extendedKern.n_rows - Nrows; 
    colStartIn = 0;
    rowEndIn   = rowStartIn + Nrows-1;
    colEndIn   = colStartIn + Ncols-1;
    
    extendedKern.submat(rowStartIn, colStartIn, rowEndIn, colEndIn) = kernel.submat(rowStartOut, colStartOut, rowEndOut, colEndOut);

    // Copy the lower left quarter of the original kernel into the upper right quarter of the extended kernel

    rowStartOut = kernel.n_rows/2;           
    colStartOut = 0;
    rowEndOut   = kernel.n_rows-1;
    colEndOut   = kernel.n_cols/2-1;
    Nrows = rowEndOut - rowStartOut + 1;
    Ncols = colEndOut - colStartOut + 1;    
 
    rowStartIn = 0; 
    colStartIn = extendedKern.n_cols - Ncols;
    rowEndIn   = rowStartIn + Nrows-1;
    colEndIn   = colStartIn + Ncols-1;

    extendedKern.submat(rowStartIn, colStartIn, rowEndIn, colEndIn) = kernel.submat(rowStartOut, colStartOut, rowEndOut, colEndOut);

    // Copy the upper left quarter of the original kernel into the lower right quarter of the extended kernel

    rowStartOut = 0;           
    colStartOut = 0;
    rowEndOut   = kernel.n_rows/2-1;
    colEndOut   = kernel.n_cols/2-1;
    Nrows = rowEndOut - rowStartOut + 1;
    Ncols = colEndOut - colStartOut + 1;    

    rowStartIn = extendedKern.n_rows - Nrows;
    colStartIn = extendedKern.n_cols - Ncols;
    rowEndIn   = rowStartIn + Nrows-1;
    colEndIn   = colStartIn + Ncols-1;

    extendedKern.submat(rowStartIn, colStartIn, rowEndIn, colEndIn) = kernel.submat(rowStartOut, colStartOut, rowEndOut, colEndOut); 
}







/**
 * \brief  Convolve the input matrix with the kernel to give the output matrix
 * 
 * \param in              Input matrix, should have the same shape as given in the constructor 
 * \param out             Output matrix, should have the same shape as the input matrix
 * \param zeroThreshold   To clip near-zero values of the output matrix to exactly zero.
 *                        if |value| < zeroThreshold -> value = 0.
 * 
 */

void Convolver::convolve(matrix &in, matrix &out, double zeroThreshold)
{
    // Check if both the input and output matrices have the same shape as was given in the constructor

    if ((out.n_rows != NrowsInOut) || (out.n_cols != NcolsInOut))
    {
        string errorMessage = "Convolver: output matrix shape (" + to_string(out.n_rows) + "," + to_string(out.n_cols)
                            + ") != expected shape (" + to_string(NrowsInOut) + "," + to_string(NcolsInOut) + ")";
        Log.error(errorMessage);
        throw IllegalArgumentException(errorMessage);
    }

   if ((in.n_rows != NrowsInOut) || (in.n_cols != NcolsInOut))
    {
        string errorMessage = "Convolver: input matrix shape (" + to_string(in.n_rows) + "," + to_string(in.n_cols)
                            + ") != expected shape (" + to_string(NrowsInOut) + "," + to_string(NcolsInOut) + ")";
        Log.error(errorMessage);
        throw IllegalArgumentException(errorMessage);
    }

    // First copy the column-major armadillo input matrix into a row-major zero-padded array

    for (int i = 0; i < in.n_rows; ++i)
    {
        for (int j = 0; j < in.n_cols; ++j)
        {
            const int ij = i*NcolsExtended + j;
            extendedIn[ij] = in(i,j);
        }
    }

    // Compute the 2D fourier transform of the input matrix 

    fftwf_execute(forwardPlanIn);

    // Multiply the complex fourier transforms of the input and kernel matrices
    // [0] contains the real part, [1] contains the imaginary part.

    for (int i = 0; i < NrowsExtended; ++i)
    {
        for (int j = 0; j < NcolsExtended/2+1; ++j) 
        {
            const int ij = i*(NcolsExtended/2+1) + j;
            fourierOut[ij][0] = (fourierIn[ij][0] * fourierKern[ij][0] - fourierIn[ij][1] * fourierKern[ij][1]);
            fourierOut[ij][1] = (fourierIn[ij][0] * fourierKern[ij][1] + fourierIn[ij][1] * fourierKern[ij][0]);
        }
    }

    // Inverse 2D-fourier-transform this product back to a real matrix. This is the convolved result.

    fftwf_execute(backwardPlanOut);

    // Copy the row-major extended convolution array into the user-given column-major armadillo array

    for (int i = 0; i < in.n_rows; ++i)
    {
        for (int j = 0; j < in.n_cols; ++j) 
        {
            const int ij = i*NcolsExtended + j;

            // FFTW does not scale the fourier transforms, so we have to do it manually

            const double value = extendedOut[ij] / (NrowsExtended * NcolsExtended);
            
            // If required, clip the values close to zero to exactly zero

            if (fabs(value) < zeroThreshold)
            {
                out(i,j) = 0.0;
            }
            else
            {
                out(i,j) = value;
            }
        }
    }
}








