from simfile import SimFile
import numpy as np
import math

def getOffset(serialPreScan, outlier_detection_k  = 10):
    """
    PURPOSE: Offset calculation as explained in PLATO-LESIA-PDC-DD-005 
             (PLATO: N_DPU Onboards Offset Calculation ATBD).
    INPUT:
        - serialPreScan: Serial pre-scan that is used to calculate the offset
        - outlier_detection-k: Number of largest and smallest values to flag as outliers
    
    OUTPUT:
        - offset_value_fc: Mean of the values in the serial pre-scan after discarding the 
                           outliers
        - offset_variance_fc: Variance of the values in the serial pre-scan after discarding 
                              the outliers
        - offset_pixels_error_number_fc: Number of pixels in the serial pre-scan that are 
                                         discarded as outliers
    """

    offset_rows = np.ravel(serialPreScan)                                   # 2D -> 1D
    offset_rows, n_useful = offset_outliers_detection(offset_outliers_detection, outlier_detection_k)

    # Approximation of the mean
    k = offset_rows[0]
    free_outliers_shifted_value = offset_rows - k
    
    pixel_value_sum = np.sum(free_outliers_shifted_value)
    square_pixel_value_sum = np.sum(np.power(free_outliers_shifted_value, 2))

    offset_value_fc = (pixel_value_sum + n_useful * k) / n_useful
    offset_variance_fc = (square_pixel_value_sum - (math.pow(pixel_value_sum, 2)) / n_useful) / (n_useful - 1)
    offset_pixels_error_number_fc = 2 * k


    # Return:
    #   - offset value
    #   - offset variance
    #   - number of discarded pixels (in case of outlier detection)

    return offset_value_fc, offset_variance_fc, offset_pixels_error_number_fc





def offset_outliers_detection(offset_rows, k = 10):

    """
    PURPOSE: Outlier detection in the serial pre-scan as explained in PLATO-MPSSR-PDC-DD-0002
             (PLATO On-Board Offset & Prescan Outlier Detection Algorithm Theoretical Baseline
             Document).
    
    INPUT:
        - offset_rows: 1D version of the serial pre-scan
        - k: Number of largest and smallest values to flag as outliers
    
    OUTPUT:
        - offset_rows: 1D version of the serial pre-scan after discarding the outliers
        - n_useful: Number of elements in the serial pre-scan after discarding the outliers
    """

    offset_rows = np.sort(offset_rows)  # Sort
    offset_rows = offset_rows[k : -k]   # Throw out the k largest and k smallest values

    n_useful = len(offset_rows) - 2 * k

    return offset_rows, n_useful





def getSmearing(simFile, exposure):

    return None
