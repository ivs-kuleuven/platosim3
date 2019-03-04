from simfile import SimFile
import numpy as np
import math

num_exposures_sc = 1
num_exposures_lc = 24


def processLeft(simFile, starID, window_dim):

    imagetteRadius = (window_dim - 1) / 2

    smearing_pattern_fc_array = np.array([])      # Smearing pattern for the different columns at fast cadence (25s)
    smearing_pattern_lc_array = np.array([])      # Smearing pattern for the difference columns at long cadence (600s)

    fx_fc_array = np.array([])      # Flux calculated with the nominal mask at fast cadence (25s)
    dfx_fc_array = np.array([])     # Difference in flux between the extended and the nominal mask at fast cadence (25s)
    ncob_fc_array = np.array([])    # COB calculated with the nominal mask at fast cadence (25s)
    ecob_fc_array = np.array([])    # COB calculated with the extende mask at fast cadence (25s)
    
    outlier_detection_b = 10
    nflag_fc = np.zeros(outlier_detection_b)
    eflag_fc = np.zeros(outlier_detection_b)

    fx_sc_array = np.array([])
    dfx_sc_array = np.array([])
    ncob_sc_array = np.array([])
    ecob_sc_array = np.array([])

    fx_lc_array = np.array([])
    fxvar_lc_array = np.array([])
    dfx_lc_array = np.array([])
    dfxvar_lc_array = np.array([])


    numExposures = simFile.getInputParameter("ObservingParameters", "NumExposures")

    for exposure in range(numExposures):

        # Calculate the offset for the current exposure (fast cadence)

        serialPreScan = simFile.getBiasMapLeft(exposure)
        outlier_detection_k = 10

        offset_value_fc = offset_calculation(serialPreScan, outlier_detection_k)[0]



        # Calculate the smearing pattern for the current exposure (fast cadence)

        parallelOverScan = simFile.getSmearingMap()

        half_ccd_gain = 1.0 / (simFile.getInputParameter("FEE", "Gain", "RefValueLeft") * simFile.getInputParameter("CCD", "Gain", "RefValueLeft")) # [e- / ADU]
        a0_array = np.empty(parallelOverScan.shape[1])
        a0_array.fill(101.3463)
        a_coefficients = np.array([0.9695, 4.6004, 5.9211])
        b_coefficients = np.array([2.1641, 0.2750, 0.0218, 0.0039])
        n0 = 10
        outlier_detection_threshold = 4.0
        epsilon = 0.01
        std_dev_previous = 9999

        smearing_pattern_fc, a0_array, std_dev_previous = smearing_calculation(parallelOverScan, offset_value_fc, half_ccd_gain, a0_array, a_coefficients, b_coefficients, n0, outlier_detection_threshold, std_dev_previous, epsilon)
        smearing_pattern_fc_array = np.append(smearing_pattern_fc_array, smearing_pattern_fc)



        # Calculate the background for the current exposure (fast cadence)
        # TODO
        # Real photometric pipeline: ~100 background windows (nominally 4x4) per CCD half,
        # for which we know the fluxes and the position on the CCD.  The latter is needed to 
        # subtract the correct smearing pattern from each background window.  For each of the
        # background windows the mean and variance are calculated, which are used in an
        # interpolation schema (based on radial basis functions) to determine the background
        # at the position of the target star.
        # It seems unfeasible to make (long-term) simulations for all these background windows, so
        # we will have to look for an alternative way to estimate the background at the target
        # location.



        # Calculate the flux & COB (nominal & extended mask) for the current exposure
        # -> add new datapoint to the light curve (fast cadence)
        # TODO
        # This will fill fx_fc (flux in the nominal window at fast cadence), 
        # dfx_fc (flux difference between the extended and the nominal mask), 
        # ncob_fc (COB as obtained in the nominal mask), and ecob_fc (COB as
        # obtained in the extended mask).

        star_window = simFile.getImagette(starID, exposure, imagetteRadius)
        fx_fc, dfx_fc, ncob_fc, ecob_fc = flux_cob_calculation(star_window, offset_value_fc, smearing_pattern_fc, half_ccd_gain, nmask, emask)
        fx_fc_array = np.append(fx_fc_array, fx_fc)
        dfx_fc_array = np.append(dfx_fc_array, dfx_fc)
        ncob_fc_array = np.append(ncob_fc_array, ncob_fc)
        ecob_fc_array = np.append(ecob_fc_array, ecob_fc)



        # Smearing time averaging (fast -> long cadence)

        if (exposure + 1) % num_exposures_lc == 0:
            smearing_pattern_lc = smearing_time_averaging(smearing_pattern_fc_array[-num_exposures_lc:], exposure)
            smearing_pattern_lc_array = np.append(smearing_pattern_lc_array, smearing_pattern_lc)



        # Light curve outlier detection
        # Note that you need b datapoints in the past and b datapoints in the future!

        outlier_detection_threshold = 5

        if (exposure + 1) > 2 * outlier_detection_b:
            
            last_nominal_fx_fc = fx_fc_array[-2 * outlier_detection_b - 1:]
            nflag_fc = flux_cob_outlier_detection(last_nominal_fx_fc, nflag_fc, outlier_detection_threshold)
            # TODO What at the end of the time series (when there are no future datapoints)?



        # Flux & COB time averaging (short cadence)

        if (exposure + 1) % num_exposures_sc == 0:

            # Nominal mask

            fx_sc, ncob_sc, n_useful = flux_cob_time_averaging_sc(fx_fc_array[-num_exposures_sc:], ncob_fc_array[-num_exposures_sc:], nflag_fc[-num_exposures_sc:])       # Select duration of short cadence (50s)
            if n_useful != 0:
                fx_sc_array = np.append(fx_sc_array, fx_sc)
                ncob_sc_array = np.append(ncob_sc_array, ncob_sc)

            # Extended mask

            dfx_sc, ecob_sc, n_useful = flux_cob_time_averaging_sc(dfx_fc_array[-num_exposures_sc:], ecob_fc_array[-num_exposures_sc:], eflag_fc[-num_exposures_sc:])     # Select duration of short cadence (50s)
            if n_useful != 0:
                dfx_sc_array = np.append(dfx_sc_array, dfx_sc)
                ecob_sc_array = np.append(ecob_sc_array, ecob_sc)

        # Flux & COB time averaging (long cadence)

        if (exposure + 1) % num_exposures_lc == 0:

            # Nominal mask

            fx_lc, fxvar_lc, n_useful = flux_cob_time_averaging_lc(fx_fc_array[-num_exposures_lc:], nflag_fc[-num_exposures_lc:])   # Select duration of long cadence (600s)
            if n_useful != 0:
                fx_lc_array = np.append(fx_lc_array, fx_lc)
                fxvar_lc_array = np.append(fxvar_lc_array, fxvar_lc)
            
            # Extended mask

            dfx_lc, dfxvar_lc, n_useful = flux_cob_time_averaging_lc(dfx_fc[-num_exposures_lc:], eflag_fc[-num_exposures_lc:])      # Select duration of long cadence (600s)
            if n_useful != 0:
                dfx_lc_array = np.append(dfx_lc_array, dfx_lc)
                dfxvar_lc_array = np.array(dfxvar_lc_array, dfxvar_lc)






####################
# Offset calculation
####################

def offset_calculation(offset_rows, outlier_detection_k  = 10):

    """
    PURPOSE: Offset calculation as explained in PLATO-LESIA-PDC-DD-005 
             (PLATO: N_DPU Onboards Offset Calculation ATBD).
             Algorithm name: ONB-OFFCAL-010.

    INPUT:
        - offset_rows: Serial pre-scan that is used to calculate the offset [ADU]
        - outlier_detection_k: Number of largest and smallest values to flag as outliers
    
    OUTPUT:
        - offset_value_fc: Mean of the values in the serial pre-scan after discarding the 
                           outliers [ADU]
        - offset_variance_fc: Variance of the values in the serial pre-scan after discarding 
                              the outliers [ADU]
    """

    # Outlier detection

    flag = offset_outliers_detection(offset_outliers_detection, outlier_detection_k)
    
    # Shifted data algorithm

    offset_value_fc, offset_variance_fc, n_useful = shifted_data_algorithm(offset_rows[~flag])

    return offset_value_fc, offset_variance_fc





def offset_outliers_detection(offset_rows, k = 10):

    """
    PURPOSE: Outlier detection in the serial pre-scan as explained in PLATO-MPSSR-PDC-DD-0002
             (PLATO On-Board Offset & Prescan Outlier Detection Algorithm Theoretical Baseline
             Document).
             Algorithm name: ONB-OFFOUTDET-010.
    
    INPUT:
        - offset_rows: Serial pre-scan that is used to calculate the offset [ADU]
        - k: Number of largest and smallest values to flag as outliers
    
    OUTPUT:
        - outliers_offset_array: Boolean truth values, where 0 means "no outlier" and 1 means
                                 "outlier".  The truth value 1 is equivalent to a flag.
    """

    offset_rows_1d = np.ravel(offset_rows)      # 2D -> 1D
    np.sort(offset_rows_1d)                     # Sort

    if len(offset_rows) < 2 * k:
        raise Exception("Number of entries, {0}, in the serial pre-scan (bias register map) should exceed 2k = {1}".format(len(offset_rows_1d), 2 * k))

    # Flag the k smallest values and the k largest values

    minOffset = offset_rows_1d[k]
    maxOffset = offset_rows_1d[-(k + 1)]

    outliers_offset_array = np.logical_or(offset_rows < minOffset, offset_rows > maxOffset)

    return outliers_offset_array





######################
# Smearing calculation
######################

def smearing_calculation(smearing_rows, offset_value_fc, half_ccd_gain, a0_array, a_coefficients, b_coefficients, n0, outlier_detection_threshold, std_dev_previous, epsilon = 0.01):

    """
    PURPOSE: Smearing calculation as explained in PLATO-LESIA-PDC-TN- (Parallel overscan rows:
             correction of the CTI) and PLATO-LESIA-PDC-DD-006 (PLATO: N-DPU Onboard Smearing
             Calculation ATBD).
             Algorithm name: ONB-SMRCAL-010.

    INPUT:
        - smearing_rows: Parallel over-scan that is used to calculate the smearing [ADU]
        - offset_value_fc: Electronic offset as calculated from the serial pre-scan [ADU]
        - half_ccd_gain: CCD gain for the detector half to which the given parallel over-scan
                         corresponds [e- / ADU]
        - a0_array: Coefficient a0 for the CTI correction (one entry per column of the parallel
                    over-scan)
        - a_coefficients: Coefficients a1, a2, and a3 for the CTI correction
        - b_coefficients: Coefficients b0, b1, b2, b3 for the CTI correction
        - n0: Number of rows in the parallel over-scan that will be skipped for the CTI correction
              (as these may be affected by bright sources at the top of the detector)
        - outlier_detection_threshold: Threshold for outlier detection (number of median absolute 
                                       deviations used to flag outliers)
        - std_dev_previous: Standard deviation of the previous measurement of the parallel over-scan
                            (one entry per column of the parallel over-scan) [e-]
        - epsilon: Regularisation factor for the CTI correction

    OUTPUT:
        - smearing_pattern_fc: One smearing row for this CCD half [e-]
        - a0_array: Coefficients a1, a2, and a3 for the CTI correction, after updating
        _ std_dev_previous: Standard deviation of the current measurement of this column of the parallel 
                            over-scan
    """

    n1 = smearing_rows.shape[0]

    smearing_rows = (smearing_rows - offset_value_fc) * half_ccd_gain   # [ADU] -> [e-]

    Ic = np.zeros(smearing_rows.shape)
    smearing_pattern_fc = np.zeros(smearing_rows.shape[1])
    
    a1, a2, a3 = a_coefficients[0], a_coefficients[1], a_coefficients[2]                            # a0 is not in the array (will be updated)
    b0, b1, b2, b3 = b_coefficients[0], b_coefficients[1], b_coefficients[2], b_coefficients[3]
    u0, u1, u2, u3 = math.exp(-b0), math.exp(-b1), math.exp(-b2), math.exp(-b3)                     # Eq. (12) in PLATO-LESIA-PDC-TN-

    for column in range(smearing_rows.shape[1]):
        
        # CTI correction

        f0, f1, f2, f3 = 1, 1, 1, 1     # Will be updated iteratively -> initialisation for i = 0
        tau = (1 + a1 + a2 + a3)        # Will be updated iteratively -> initialisation for i = 0 -> Eq. (3) in PLATO-LESIA-PDC-TN-

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
        
        flag = smearing_outlier_detection(Ic[n0:n1][column], outlier_detection_threshold, std_dev_previous[column])


        if(np.sum(flag) == np.size(flag)):

            smearing_pattern_fc[column] = 0

        else:

            # Calculation of the mean smearing (first n0 measurements are excluded)

            smearing_pattern_fc[column], std_dev_previous[column], n_useful = shifted_data_algorithm(Ic[n0:][column][~flag])

            # Update a0 for the current column

            chi, rho = 0, 0
            f0, f1, f2, f3 = 1, 1, 1, 1 # Will be updated iteratively -> initialisation for i = 0
            tau = (1 + a1 + a2 + a3)    # Will be updated iteratively -> initialisation for i = 0 -> Eq. (3) in PLATO-LESIA-PDC-TN-

            for i in range(n1):
            
                if (i >= n0) and (flag[i - n0] == 0):
                    chi += (smearing_rows[i][column] - smearing_pattern_fc[column]) * tau
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
    
    return smearing_pattern_fc, a0_array, std_dev_previous





def smearing_outlier_detection(ctiCorrectedColumn, threshold, std_dev_previous):

    """
    PURPOSE: Outlier detection in the parallel over-scan as explained in PLATO-MPSSR-PDC-PT-0003
             (PLATO Onboards Overscan Outlier Detection Algorithm Theoretical Baseline Document).
             Algorithm name: ONB-OVEROUTDET-010.
    
    INPUT:
        - ctiCorrectedColumn: Column from the parallel over-scan, after CTI correction (and discarding
                              rows to avoid contamination by bright sources at the top of the detector)
        - threshold: Threshold for outlier detection (number of median absolute deviations used to flag
                     outliers)
        - std_dev_previous: Standard deviation of the previous measurement of this column of the parallel 
                            over-scan
    
    OUTPUT:
        - raw_overscan_flags: Boolean truth values, where 0 means "no outlier" and 1 means "outlier".  The 
                              truth value 1 is equivalent to a flag.
    """

    median = np.median(ctiCorrectedColumn)

    raw_overscan_flags = np.ones(len(ctiCorrectedColumn))
    raw_overscan_flags[ctiCorrectedColumn - median >= threshold * std_dev_previous]

    return raw_overscan_flags                                           





###############################
# Background window calculation
###############################

def processBackgroundWindow(background_window, offset_value_fc, smearing_pattern_fc, half_ccd_gain, outlier_detection_threshold):

    """
    PURPOSE: Background calculation as explained in PLATO-LESIA-PDC-DD-011 (PLATO Onboard Background
             Window Calculation ATBD).
             Algorithm name: ONB-BKGCAL-010.

    INPUT:
        - background_window: values in the background window [ADU]
        - offset_value:_fc Electronic offset as calculated from the serial pre-scan
        - smearing_pattern_fc: Smearing as calculated from the parallel over-scan
        - half_ccd_gain: CCD gain for the detector half to which the given parallel over-scan
                         corresponds
        - outlier_detection_threshold: Threshold for outlier detection
    
    OUTPUT:
        - background_value_fc: Mean of the values in the background window [e-]
        - background_variance_fc: Variance of the values in the background window [e-]
    """

    window_smearing_pattern_fc = window_smearing_pattern(smearing_pattern_fc)

    # Subtract offset, multiply with gain, and subtract smearing pattern

    background_window = background_window - offset_value_fc

    background_window *= half_ccd_gain

    for column in range(background_window.shape[1]):

        background_window[:][column] = background_window[:][column] - window_smearing_pattern_fc[column]
    


    # Outlier detection

    flag = background_outlier_detection(background_window, outlier_detection_threshold)[0]

    # Shifted data algorithm

    background_value_fc, background_variance_fc, n_useful = shifted_data_algorithm(background_window[~flag])

    return background_value_fc, background_variance_fc





def background_outlier_detection(background_window, threshold):

    """
    PURPOSE: Detection of outliers in the background window as explained in PLATO-MPSSR-PDC-DD-004
             (PLATO On-board Background Window Outlier Detection Algorithm Theoretical Baseline
             Document).
             Algorithm name: ONB-BACKOUTDET-010.

    INPUT:
        - background_window: Array of background window pixels [-e]
        - threshold: Threshold for outlier detection (number of median absolute deviations used to flag
                     outliers)

    OUTPUT:
        - raw_back_flags: Boolean truth values, where 0 means "no outlier" and 1 means "outlier".  The
                          truth value 1 is equivalent to a flag
    """

    raw_back_flags = mad_median_clipping(background_window, threshold)

    return raw_back_flags





########################################################
# Flux & COB calculations using nominal & extended masks
########################################################

def flux_cob_calculation(star_window, offset_value_fc, smearing_pattern_fc, half_ccd_gain, nmask, emask):

    """
    PURPOSE: Flux and COB calculation as explained in PLATO-LESIA-PDC-DD-008 
             (PLATO Onboards Flux & COB Calculation ATBD).
             Algorithm name: ONB-FXCOBCAL-010.

    INPUT:
        - star_window: Array of star window pixels [e-]
        - offset_value_fc: Electronic offset as calculated from the serial pre-scan
        - smearing_pattern_fc: Smearing as calculated from the parallel over-scan
        - half_ccd_gain: CCD gain for the detector half to which the given parallel over-scan
                         corresponds
        - nmask: Nominal mask (computed on-ground and uploaded)
        - emask: Extended mask (computed on-ground and uploaded)

    OUTPUT:
        - fx_fc: Flux computed using the nominal mask
        - dfx_fc: Flux difference between the extended and the nominal mask
        - ncob_fc: COB computed using the nominal mask
        - ecob_fc: COB computed using the extended mask
    """

    window_smearing_pattern_fc = window_smearing_pattern(smearing_pattern_fc)
    #numRows = star_window.shape[0]
    #numColumns = star_window.shape[1]

    rowRange = np.arange(star_window.shape[0])
    columnRange = np.arange(star_window.shape[1])

    pixel_center_offset = 0.5

    # Variable name conventions from Sect. 3.6.1 in PLATO-LESIA-PDC-DD-008

    nmask_column_integral = np.sum(nmask, axis = 0)
    emask_column_integral = np.sum(emask, axis = 0)

    nmask_integral = np.sum(nmask_column_integral)
    emask_integral = np.sum(emask_column_integral)

    nmask_smearing_integral = np.dot(window_smearing_pattern_fc, nmask_column_integral)
    emask_smearing_integral = np.dot(window_smearing_pattern_fc, emask_column_integral)

    nmasked_flux_pixel = np.multiply(star_window, nmask)
    emasked_flux_pixel = np.multiply(star_window, emask)

    nmask_flux_column_integral = np.sum(nmasked_flux_pixel, axis = 0)
    emask_flux_column_integral = np.sum(emasked_flux_pixel, axis = 0)

    nmasked_flux_integral = np.sum(nmask_flux_column_integral)
    emasked_flux_integral = np.sum(emask_flux_column_integral)

    nmasked_Xcob_sum = np.dot(columnRange, nmask_flux_column_integral)
    emasked_Xcob_sum = np.dot(columnRange, emask_flux_column_integral)

    nmasked_Ycob_sum = np.dot(np.sum(nmasked_flux_pixel, axis = 1), rowRange)
    emasked_Ycob_sum = np.dot(np.sum(emasked_flux_pixel, axis = 1), rowRange)

    nmask_column_Xweighted = np.multiply(columnRange, nmask_column_integral)
    emask_column_Xweighted = np.multiply(columnRange, emask_column_integral)

    nmask_column_Yweighted = np.sum(np.multiply(nmask, columnRange), axis = 0)
    emask_column_Yweighted = np.sum(np.multiply(emask, columnRange), axis = 0)

    nmask_Xweighted_integral = np.sum(nmask_column_Xweighted)
    emask_Xweighted_integral = np.sum(emask_column_Xweighted)

    nmask_Yweighted_integral = np.sum(nmask_column_Yweighted)
    emask_Yweighted_integral = np.sum(emask_column_Yweighted)

    nmask_smearing_Xweighted_integral = np.dot(window_smearing_pattern, nmask_column_Xweighted)
    emask_smearing_Xweighted_integral = np.dot(window_smearing_pattern, emask_column_Xweighted)

    nmask_smearing_Yweighted_integral = np.dot(window_smearing_pattern, nmask_column_Yweighted)
    emask_smearing_Yweighted_integral = np.dot(window_smearing_pattern, emask_column_Yweighted)

    # Calculate the flux in the nominal mask & the difference in flux between the extended and the nominal mask

    fx_fc = (nmasked_flux_integral - offset_value_fc * nmask_integral) * half_ccd_gain - nmask_smearing_integral    # Flux in nominal mask
    dfx_fc = (emasked_flux_integral - offset_value_fc * emask_integral) * half_ccd_gain - emask_smearing_integral   # Flux difference between extended & nominal mask

    # Calculate the COB for the nominal and for the extended mask

    nmasked_Xcob_sum_corrected = (nmasked_Xcob_sum - offset_value_fc * nmask_Xweighted_integral) * half_ccd_gain - nmask_smearing_Xweighted_integral
    nmasked_Ycob_sum_corrected = (nmasked_Ycob_sum - offset_value_fc * nmask_Yweighted_integral) * half_ccd_gain - nmask_smearing_Yweighted_integral

    ncob_fc = [nmasked_Xcob_sum_corrected / fx_fc + pixel_center_offset, nmasked_Ycob_sum_corrected / fx_fc + pixel_center_offset]

    emasked_Xcob_sum_corrected = (emasked_Xcob_sum - offset_value_fc * emask_Xweighted_integral) * half_ccd_gain - emask_smearing_Xweighted_integral
    emasked_Ycob_sum_corrected = (emasked_Ycob_sum - offset_value_fc * emask_Yweighted_integral) * half_ccd_gain - emask_smearing_Yweighted_integral

    masked_Xcob_sum_corrected = nmasked_Xcob_sum_corrected + emasked_Xcob_sum_corrected
    masked_Ycob_sum_crrected = nmasked_Ycob_sum_corrected + emasked_Ycob_sum_corrected
    efx_fc = fx_fc + dfx_fc

    ecob_fc = [masked_Xcob_sum_corrected / efx_fc + pixel_center_offset, masked_Ycob_sum_crrected / efx_fc + pixel_center_offset]
    
    return fx_fc, dfx_fc, ncob_fc, ecob_fc





#########################
# Smearing time averaging
#########################

def smearing_time_averaging(smearing_pattern_fc_array, exposure):

    """
    PURPOSE: Smearing time averaging for long cadence as explained in PLATO-LESIA-PDC-DD-007
             (PLATO Ob-board Smearing Pattern Time Averaging ATBD).
             Algorithm name: ONB-SMRAVG-010

    INPUT:
        - smearing_pattern_fc_array: Smearing pattern computed every 25s (fast cadence) [e-], for a 
                               duration of 600s
    
    OUTPUT:
        - smearing _pattern_lc: Smearing pattern averaged over the given 24 samples (i.e. 600s) [e-]
    """

    shape = (1, num_exposures_lc, smearing_pattern_fc_array.shape[1], 1)
    smearing_pattern_lc = smearing_pattern_fc_array.reshape(shape).mean(-1).mean(1)

    return smearing_pattern_lc






####################################
# Outlier detection over light curve
####################################

def flux_cob_outlier_detection(fx_fc_array, fx_exposure_error_array, threshold):

    """
    PURPOSE: Outlier detection over light curve as explained in PLATO-MPSSR-PDC-DD-0001
             (PLATO Onboard LC Outlier Detection ATBD).
             Algorithm name: ONB-LCOUTDET-010.

    INPUT:
        - flux_fc_array: New segment in the flux time series (of length 2 * b + 1)
        - fx_exposure_error_array: Flag for the time series processed so far
        - threshold: Threshold for outlier detection

    OUTPUT:
        - fx_exposure_error_array: Flag for the time series processed so far, including
                                   the given time series segment
    """

    b = len(fx_fc_array + 1) / 2

    # Step 1 - 4

    isOutlier = mad_median_clipping(fx_fc_array, threshold, b)

    # Step 5

    if fx_exposure_error_array[-2] and fx_exposure_error_array[-1]:

        fx_exposure_error_array[-2] = False

        if not isOutlier:
            fx_exposure_error_array[-1] = False

    fx_exposure_error_array = np.append(fx_exposure_error_array, isOutlier)

    return fx_exposure_error_array





###########################
# Flux & COB time averaging
###########################

def flux_cob_time_averaging_sc(fx_fc_array, cob_fc_array, fx_exposure_error_array):

    """
    PURPOSE: Flux time averaging for short cadence (50s) as explained in PLATO-LESIA-PDC-DD-009
             (PLATO Onboard Flux and COB Short Cadence Time Averaging ATBD).
             Algorithm name: ONB-FXAGV-011.

    INPUT:
        - fx_fc_array: Flux time series at fast cadence (25s) for a duration of 50s.  This can
                       either be the flux obtained in the nominal mask (fx_fc) or the difference 
                       in flux between the extended and the nominal mask (dfx_sc).
        - cob_fc_array: COB time series at fast cadence (25s) for a duration of 50s.  This can
                        either be the COB obtained in the nominal mask (ncob_fc) or in the extended
                        mask (ecob_fc).
        - fx_exposure_error_array: Outlier detection flags for the time series at fast cadence (25s) 
                                   for a duration of 50s
    
    OUTPUT:
        - fx_sc: Mean of the given flux time series (over 50s), calculated with the non-flagged datapoints only
        - cob_sc: Mean of the given COB time series (over 50s), calculated with the non-flagged datapoints only
        - n_useful_sc: Number of non-flagged datapoints in the time series
    """

    n_useful_sc = len(fx_exposure_error_array) - np.sum(fx_exposure_error_array)

    fx_sc = np.sum(fx_fc_array[~fx_exposure_error_array]) / n_useful_sc

    cob_sc_x = np.sum(cob_fc_array[~fx_exposure_error_array][0]) / n_useful_sc
    cob_sc_y = np.sum(cob_fc_array[~fx_exposure_error_array][1]) / n_useful_sc
    cob_sc = [cob_sc_x, cob_sc_y]

    return fx_sc, cob_sc, n_useful_sc





def flux_cob_time_averaging_lc(fx_fc_array, fx_exposure_error_array):

    """
    PURPOSE: Flux time averaging for long cadence (600s) as explained in PLATO-LESIA-PDC-DD-010
             (PLATO Onboard Flux and COB Long Cadence Time Averaging ATBD).
             Algorithm name: ONB-FXAGV-012.

    INPUT:
        - fx_fc_array: Flux time series at fast cadence (25s) for a duration of 600s.  This can
                       either be the flux obtained in the nominal mask (fx_fc) or the difference 
                       in flux between the extended and the nominal mask (dfx_sc).
        - fx_exposure_error_array: Outlier detection flags for the time series at fast cadence (25s) 
                                   for a duration of 600s
    OUTPUT:
        - fx_lc: Mean of the given flux time series (over 600s), calculated with the non-flagged datapoints only
        - fxvar_lc: Variance of the given flux time series (over 600s), calculated with the non-flagged datapoints only
        - n_useful_lc: Number of non-flagged datapoints in the time series
    """

    fx_lc, fxvar_lc, n_useful_lc = shifted_data_algorithm(fx_fc_array[~fx_exposure_error_array])

    return fx_lc, fxvar_lc, n_useful_lc





#####################
# Auxiliary functions
#####################

def shifted_data_algorithm(data):

    """
    PURPOSE: Shifted data algorithm to compute the mean and variance for the given data.  This method
             avoids loss of significance with big numbers when computing the variance.

    INPUT:
        - data: Data array for which to calculated the mean and the variance

    OUTPUT:
        - mean: Mean value of the given data array, as computed with the shifted data algorithm
        - variance: Variance of the given data array, as computed with the shifted data algorithm
        - n_useful: Number of non-flagged datapoints
    """

    k = data[0]         # Use the 1st non-flagged value as approximation of the mean (arbitrary choice)
    data = data - k     # Shift the data

    data_shifted_sum = np.sum(data)                             # Sum of the shifted data
    data_shifted_squared_sum = np.sum(np.power(data, 2))        # Sum of the squares of the shifted data

    n_useful = len(data)

    mean = (data_shifted_sum + n_useful * k) / n_useful                                                 # Mean
    variance = (data_shifted_squared_sum - (pow(data_shifted_sum, 2) / n_useful)) / (n_useful - 1)      # Variance

    return (mean, variance, n_useful)





def mad_median_clipping(data, threshold, index = None):

    """
    PURPOSE: Outlier detection, based on Median Absolute Deviation (MAD) clipping around the median.

    INPUT:
        - data: Data for which to flag outliers
        - threshold: Threshold for MAD clipping around the median
        - index: If specified, it is calculated whether the datapoint at this position is an outlier.

    OUTPUT:
        - in case index was specified: boolean indicating whether or not the datapoint at the given
          position is an outlier; otherwise: boolean truth values, where 0 means "no outlier" and 1 means
                                 "outlier".  The truth value 1 is equivalent to a flag.
    """

    median = np.median(data)         # Median
    y = np.fabs(data - median)       # Subtract median from background window
    mad = np.median(y)               # MAD

    if index == None:

        flag = np.where(y > threshold * mad)
        return flag

    else:
        return y[index] / mad > threshold





def window_smearing_pattern(smearing_pattern_fc):

    # TODO We will need extra parameters (spanned columns) to extract the required smearing
    # pattern for the window.

    return smearing_pattern_fc
