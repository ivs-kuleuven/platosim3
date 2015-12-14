#ifndef SUBFIELD_H
#define SUBFIELD_H

#include <math.h>
#include "configurationparameters.h"
#include "starcatalog.h"

using namespace std;

class SubField
{
public:

	SubField(ConfigurationParameters configurationParameters,
			double detectorOriginOffsetX, double detectorOriginOffsetY,
			double detectorOrientation, double pixelSize);
	~SubField();

	double getDistanceFromOpticalAxisToFieldCenter();

	void addFlux(double xCoord, double yCoord, double flux);
	void addFlux(double);

	void convolveWithPsf(double **psf);
	void multiply(double **array);
	void rebin();
	void reset();

protected:

private:

	// Number of sub-pixels per pixel

	double numSubPixelsPerPixel;

	// Offset of the detector origin w.r.t. the centre of the focal plane in
	// both directions [mm]

	double detectorOriginOffsetX;
	double detectorOriginOffsetY;

	// Sub-pixel size [mm/pixel]

	double pixelSize;

	// Orientation of the detector w.r.t. the focal plane, measured
	// counterclockwise [radians]

	double detectorOrientation;

	// Position of the sub-field zeropoint w.r.t. the complete detector [pixels]

	double subFieldZeroPointX;
	double subFieldZeroPointY;

	// Size of the sub-field in both directions [pixels]

	double subFieldSizeX;
	double subFieldSizeY;

	// Number of pixels to extend the sub-field by on each side, to account for
	// the edge effect

	int numEdgePixels;

	// Sub-pixel map and its dimensions

	double **subPixelMap;

	int subPixelMapSizeX;
	int subPixelMapSizeY;

};

#endif
