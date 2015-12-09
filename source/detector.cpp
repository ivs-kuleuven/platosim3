
#include "detector.h"


// Constructor

Detector::Detector()
{

}




// Destructor

Detector::~Detector()
{
    // Should also call the destructor of subField.
}




// Detector::takeExposure()
//
// PURPOSE:
// 
// INPUT:
//
// OUTPUT:
//

Detector::takeExposure()
{

    integrateLight();

    readOut();
}








// Detector::integrateLight()
//
// PURPOSE:
// 
// INPUT:
//
// OUTPUT:
//

Detector::integrateLight()
{
    // Get rid of the previous exposure, by zeroing the entire subfield.

    subField.reset();



    camera.exposeSubField(subField);

    // Apply the pixel response non-uniformity (flatfield) at subpixel level.
    // (Sub)pixel units before: [electrons]
    // (Sub)pixel units after: [electrons]

    subField.multiply(flatfield);

    // Rebin the subpixel map to the pixel map.

    pixelMap = subField.rebin();
}






// Detector::readOut()
//
// PURPOSE:
// 
// INPUT:
//
// OUTPUT:
//

Detector::readOut()
{
    // Apply quantum efficiency
    // Pixel units before: [photons]
    // Pixel units after: [electrons]

    applyQuantumEfficiency();

    // Apply poisson distributed photon noise. 
    // Pixel units before: [electrons]
    // Pixel units after: [electrons]

    addPhotonNoise();

    // Apply full well saturation. A pixel has a maximum capacity of electrons (the full well capacity). 
    // If photons free more electrons, the pixel saturates, and the electrons flow in the pixels above and below in 
    // the same column (potential barriers are smallest in that direction).
    // Pixel units before: [electrons]
    // Pixel units after: [electrons]

    applyFullWellSaturation();

    // Simulate the effects of the charge Transfer Inefficiency (CTI). When the 
    // CCD is read out, row after row, a part of the charge is always left behind
    // which then dribbles into the trailing pixels. This causes each star to have
    // a small "tail". Only visible when the CTI = 1 - CTE is poor.
    // Pixel units before: [electrons]
    // Pixel units after: [electrons]        

    applyChargeTransferSmearing();

    // Apply the effects of readout smearing due to an open shutter. Because there is no shutter,
    // the pixels are still receiving photons from the sky, while they are being transfered towards
    // the readout register.

    applyOpenShutterSmearing();

    // Each time the amplifier reads out a pixel, a tiny bit of noise is added.
    // Add the readout noise.
    // Pixel units before: [electrons]
    // Pixel units after: [electrons]
    
    addReadoutNoise();

    // Apply the gain, to increase the dynamic range of the detector.
    // Pixel units before: [electrons]
    // Pixel units after: [ADU]

    applyGain();

    // Take into account the bias level. I.e. add the constant "zero" level
    // introduced by the amplifier.
    // Pixel units before: [ADU]
    // Pixel units after: [ADU]

    addBias();

    // Take into acount digital saturation. If even after dividing by the gain
    // the number of ADUs in a pixel is still higher than the analog-digital 
    // converter (ADC) can represent with its fixed amount of bits, clip all 
    // values that are too high to the saturation level of the ADC.
    // Pixel units before: [ADU]
    // Pixel units after: [ADU]

    applyDigitalSaturation();

}
