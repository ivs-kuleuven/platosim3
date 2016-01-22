
#include "camera.h"


// Constructor 

Camera::Camera()
{
    
}





// Destructor

Camera::~Camera()
{

}






/**
 * @brief      Expose the subField to the Sky, i.e. add flux to the detectors, 
 *             add Background and convolve with the PSF.
 *
 * @param[in]  subField  { parameter_description }
 */
void Camera::exposeSubField(SubField subField)
{
    auto starCatalog = sky.getStarsWithinRadiusFrom(alpha, delta, radius);  
    double skyBackground = sky.getSkyBackground(alpha, delta)  

    double tickInterval = telescope.getTickInterval();

    while (currentTime < startingTime + exposureTime)
    {
        telescope.updatePointingCoordinates(raOpticalAxis, decOpticalAxis, tickInterval);
        currentTime += tickInterval;

        for (auto star : starCatalog)
        {
            computeFocalPlaneCoordinates(star, Xmm, Ymm)
            
            if (subField.containsPoint(Xmm, Ymm))
            {
                subField.addFlux(Xmm, Ymm, flux);
            }
        }
    }

    subField.add(skyBackground * exposureTime);
    subField.convolveWithPSF(psf);

    return;
}




