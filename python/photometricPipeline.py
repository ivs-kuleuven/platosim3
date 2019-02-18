from simfile import SimFile
import numpy as np
import math

####################
# Offset calculation
####################

def getOffset(serialPreScan, outlier_detection_k  = 10):

    """
    PURPOSE: Offset calculation as explained in PLATO-LESIA-PDC-DD-005 
             (PLATO: N_DPU Onboards Offset Calculation ATBD).
    INPUT:
        - serialPreScan: Serial pre-scan that is used to calculate the offset
        - outlier_detection_k: Number of largest and smallest values to flag as outliers
    
    OUTPUT:
        - offset_value_fc: Mean of the values in the serial pre-scan after discarding the 
                           outliers
        - offset_variance_fc: Variance of the values in the serial pre-scan after discarding 
                              the outliers
        - offset_pixels_error_number_fc: Number of pixels in the serial pre-scan that are 
                                         discarded as outliers
    """

    offset_rows = np.ravel(serialPreScan)                                                               # 2D -> 1D
    offset_rows, n_useful = offset_outliers_detection(offset_outliers_detection, outlier_detection_k)   # Outlier detection

    
    k = offset_rows[0]                                  # Use the 1st non-flagged value as approximation of the mean
    free_outliers_shifted_value = offset_rows - k       # Shift the data
    
    pixel_value_sum = np.sum(free_outliers_shifted_value)                       # Sum of the shifted data
    square_pixel_value_sum = np.sum(np.power(free_outliers_shifted_value, 2))   # Sum of the squares of the shifted data

    offset_value_fc = (pixel_value_sum + n_useful * k) / n_useful                                               # Mean
    offset_variance_fc = (square_pixel_value_sum - (math.pow(pixel_value_sum, 2)) / n_useful) / (n_useful - 1)  # Variance
    offset_pixels_error_number_fc = 2 * k                                                                       # Number of outliers


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





######################
# Smearing calculation
######################

def getSmearing(parallelOverScan, offset_value, half_ccd_gain, a0_array, a_coefficients, b_coefficients, n0, outlier_detection_threshold, std_dev_previous, epsilon = 0.01):

    """
    PURPOSE: Smearing calculation as explained in PLATO-LESIA-PDC-TN- (Parallel overscan rows:
             correction of the CTI) and PLATO-LESIA-PDC-DD-006 (PLATO: N-DPU Onboard Smearing
             Calculation ATBD).

    INPUT:
        - parallelOverScan: Parallel over-scan that is used to calculate the smearing
        - offset_value: Electronic offset as calculated from the serial pre-scan
        - half_ccd_gain: CCD gain for the detector half to which the given parallel over-scan
                         corresponds
        - a0_array: Coefficient a0 for the CTI correction (one entry per column of the parallel
                    over-scan)
        - a_coefficients: Coefficients a1, a2, and a3 for the CTI correction
        - b_coefficients: Coefficients b0, b1, b2, b3 for the CTI correction
        - n0: Number of rows in the parallel over-scan that will be skipped for the CTI correction
              (as these may be affected by bright sources at the top of the detector)
        - outlier_detection_threshold: Threshold for outlier detection (number of median absolute 
                                       deviations used to flag outliers)
        - std_dev_previous: Standard deviation of the previous measurement of the parallel over-scan
                            (one entry per column of the parallel over-scan)
        - epsilon: Regularisation factor for the CTI correction

    OUTPUT:
        - smearing_pattern_fc: One smearing row for this CCD half
        - a0_array: Coefficients a1, a2, and a3 for the CTI correction, after updating
    """

    n1 = parallelOverScan.shape[0]

    smearing_rows = (parallelOverScan - offset_value) * half_ccd_gain

    Ic = np.zeros(parallelOverScan.shape)
    smearing_pattern_fc = np.zeros(parallelOverScan.shape[1])
    
    a1, a2, a3 = a_coefficients[0], a_coefficients[1], a_coefficients[2]    # a0 is not in the array
    b0, b1, b2, b3 = b_coefficients[0], b_coefficients[1], b_coefficients[2], b_coefficients[3]
    u0, u1, u2, u3 = math.exp(-b0), math.exp(-b1), math.exp(-b2), math.exp(-b3) # Eq. (12) in PLATO-LESIA-PDC-TN-

    for column in range(parallelOverScan.shape[1]):
        
        # CTI correction

        f0, f1, f2, f3 = 1, 1, 1, 1 # Will be updated iteratively -> initialisation for i = 0
        tau = (1 + a1 + a2 + a3)    # Will be updated iteratively -> initialisation for i = 0 -> Eq. (3) in PLATO-LESIA-PDC-TN-

        for i in range(n1):
            if i >= n0:
                # Correct the value
                Ic[i][column] = smearing_rows[i][column] - a0_array[column] * tau
                
            # Update f0, f1, f2, f3, and tau
            
            f0 *= u0
            f1 *= u1
            f2 *= u2
            f3 *= u3

            tau = f0 + a1 * f1 + a2 * f2 + a3 * f3

        # Outlier detection
        
        flag, n_useful = smearing_outlier_detection(Ic[n0:n1][column], outlier_detection_threshold, std_dev_previous[column])

        if(n_useful == 0):

            smearing_pattern_fc[column] = 0

        else:

            # Calculation of the mean smearing (first n0 measurements are excluded)

            smearing_pattern_fc[column] = np.mean(Ic[n0:][flag], axis = 0)


            # Update a0 for the current column

            chi, rho = 0, 0
            f0, f1, f2, f3 = 1, 1, 1, 1 # Will be updated iteratively -> initialisation for i = 0
            tau = (1 + a1 + a2 + a3)    # Will be updated iteratively -> initialisation for i = 0 -> Eq. (3) in PLATO-LESIA-PDC-TN-

            for i in range(n1):
            
                if (i >= n0) and (flag[i - n0] == 1):
                    chi += (Ic[i][column] - smearing_pattern_fc[column]) * tau
                    rho += pow(tau, 2)

                # Update f0, f1, f2, f3, and tau
            
                f0 *= u0
                f1 *= u1
                f2 *= u2
                f3 *= u3

                tau = f0 + a1 * f1 + a2 * f2 + a3 * f3

            a0_array[column] = (chi + epsilon * rho * a0_array[column]) / (rho * (1 + epsilon))

    # Return
    #   - one smearing row for the detector half
    #   - coefficients a0 for the CTI correction, after update
    
    return smearing_pattern_fc, a0_array





def smearing_outlier_detection(ctiCorrectedColumn, threshold, std_dev_previous):

    """
    PURPOSE: Outlier detection in the parallel over-scan as explained in PLATO-MPSSR-PDC-PT-0003
             (PLATO Onboards Overscan Outlier Detection Algorithm Theoretical Baseline Document).
    
    INPUT:
        - ctiCorrectedColumn: Column from the parallel over-scan, after CTI correction (and discarding
                              rows to avoid contamination by bright sources at the top of the detector)
        - threshold: Threshold for outlier detection (number of median absolute deviations used to flag
                     outliers)
        - std_dev_previous: Standard deviation of the previous measurement of this column of the parallel 
                            over-scan
    
    OUTPUT:
        - flag: Flag for the outliers in the given column of the parallel over-scan (0: outlier; 1: not an
                outlier)
        - n_useful: Number of elements in the given column of the parallel over-scan after discarding the outliers
    """

    median = np.median(ctiCorrectedColumn)

    flag = np.ones(len(ctiCorrectedColumn))
    flag[ctiCorrectedColumn - median >= threshold * std_dev_previous]

    n_useful = np.sum(flag)

    return flag, n_useful





###############################
# Background window calculation
###############################

def getBackground(offset_value, half_ccd_gain, smearing_pattern, outlier_detection_threshold):

    """
    PURPOSE: Background calculation as explained in PLATO-LESIA-PDC-DD-011 (PLATO Onboards Background
             Window Calculation ATBD).
    
    INPUT:
        -
        - offset_value: Electronic offset as calculated from the serial pre-scan
        - half_ccd_gain: CCD gain for the detector half to which the given parallel over-scan
                         corresponds
        - smearing_pattern: Smearing as calculated from the parallel over-scan
        - outlier_detection_threshold: Threshold for outlier detection
    
    OUTPUT:
        - background_value_fc: Mean of the values in the background window
        - background_variance_fc: Variance of the values in the background window
        - background_pixels_error_number_fc: Number of pixels in the background window that are 
                                             discarded as outliers
    """

    return None





def background_outlier_detection(backgroundWindow, threshold):

    """
    PURPOSE: Detection of outliers in the background window as explained in PLATO-MPSSR-PDC-DD-004
             (PLATO On-board Background Window Outlier Detection Algorithm Theoretical Baseline
             Document).

    INPUT:
        - backgroundWindow: Array of background window pixels.
        - threshold: Threshold for outlier detection (number of median absolute deviations used to flag
                     outliers)

    OUTPUT:
        - flag: Flag for the outliers in the background window (0: outlier; 1: not an outlier)
        - n_useful: Number of elements in the background window after discarding the outliers
    """

    return None





########################################################
# Flux & COB calculations using nominal & extended masks
########################################################

#########################
# Smearing time averaging
#########################

def smearingTimeAveraging(smearing_pattern_fc):

    """
    PURPOSE: Smearing time averaging for long cadence as explained in PLATO-LESIA-PDC-DD-007
             (PLATO Ob-board Smearing Pattern Time Averaging ATBD).

    INPUT:
        - smearing_pattern_fc: Smearing pattern computed every 25s.
    
    OUTPUT:
        - smearing _pattern_lc: Smearing pattern averaged over 600s.
    """

    numSamplesToAverage = 24

    numRowsOld = smearing_pattern_fc.shape[0]
    numRowsNew = numRowsOld // numSamplesToAverage
    numColumns = smearing_pattern_fc.shape[1]

    shape = (numRowsNew, numSamplesToAverage, numColumns, 1)
    smearing_pattern_lc = smearing_pattern_fc.reshape(shape).mean(-1).mean(1)

    # TODO Use time stamp of the 1st exposure of the time averaging cycle

    return smearing_pattern_lc






####################################
# Outlier detection over light curve
####################################

def outlier_detection_lightcurve(flux, threshold):

    """
    PURPOSE: Outlier detection over light curve as explained in PLATO-MPSSR-PDC-DD-0001
             (PLATO Onboard LC Outlier Detection ATBD).
    """

    return None





###########################
# Flux & COB time averaging
###########################

def fluxTimeAveraging():

    """
    PURPOSE: Flux time averaging for long cadence as explained in PLATO-LESIA-PDC-DD-010
             (PLATO Onboard Flux Long Cadence Time Averaging ATBD).
    """

    return None