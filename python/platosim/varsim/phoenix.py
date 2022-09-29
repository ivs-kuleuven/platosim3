#!/usr/bin/env python

"""
This is a script holding all relevant PHOENIX download and read features.
"""

import pathlib
import zipfile
import urllib.request
import numpy as np
from astropy.io import fits
from platosim.utilities import errorcode
from platonium.var.utilities import find_nearest


class Phoenix(object):
    """
    This class provides a convenient handling of PHOENIX spectra.
    For a given set of effective Temperature, log g, metalicity and alpha,
    it downloads the spectrum from PHOENIX ftp server.

    Webpage : https://phoenix.astro.physik.uni-goettingen.de/
    Download: http://phoenix.astro.physik.uni-goettingen.de/data/
    """

    def __init__(self):
        """
        PURPOSE: Initialize and prepare data structure
        """
        # Create data directories if they do not exist

        self.path = pathlib.Path(__file__).parent.resolve()
        pathlib.Path(self.path.joinpath('data')).mkdir(parents=False, exist_ok=True)


    def parameter_space(self, library):
        """
        PURPOSE: Return parameter space
        """

        if library == 'Atmos':
            self.valid_t = np.array([*list(range(2300, 7000, 100)),
                            *list((range(7000, 12200, 200)))])
            self.valid_g = np.array([*list(np.arange(0, 6, 0.5))])
            self.valid_z = np.array([*list(np.arange(-4, -2, 1)),
                            *list(np.arange(-2.0, 1.5, 0.5))])
            self.valid_a = np.array([-0.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])

        if library == 'HiRes':
            self.valid_t = np.array([*list(range(2300, 7000, 100)),
                            *list((range(7000, 12200, 200)))])
            self.valid_g = np.array([*list(np.arange(0, 6, 0.5))])
            self.valid_z = np.array([*list(np.arange(-4, -2, 1)),
                            *list(np.arange(-2.0, 1.5, 0.5))])
            self.valid_a = np.array([-0.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])

        if library == 'MedRes':
            self.valid_t = np.array([*list(range(2300, 7000, 100)),
                            *list((range(7000, 12200, 200)))])
            self.valid_g = np.array([*list(np.arange(-5.0, 1.0, 0.5))])
            self.valid_z = np.array([*list(np.arange(-6, 0.5, 0.5))])
            self.valid_a = np.array([0])

        if library == 'SpecInt':
            self.valid_t = np.array([*list(range(2300, 5000, 100)),
                                     *list(range(5100, 7000, 100)),
                                     *list(range(7000, 12200, 200))])
            self.valid_g = np.array([*list(np.arange(0, 6.5, 0.5))])
            self.valid_z = np.array([-4.0, -3.0, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0])
            self.valid_a = np.array([0])




            
    def nearest_parameters(self, Teff, logg, Z, alpha, library):
        """
        PURPOSE: Return compliant parameter space
        """

        # Find indices for closest valid parameter

        self.parameter_space(library)
        dex_t = find_nearest(self.valid_t, Teff)
        dex_g = find_nearest(self.valid_g, logg)
        dex_z = find_nearest(self.valid_z, Z)
        dex_a = find_nearest(self.valid_a, alpha)

        # Select valid parameter

        Teff  = self.valid_t[dex_t]
        logg  = self.valid_g[dex_g]
        Z     = self.valid_z[dex_z]
        alpha = self.valid_a[dex_a]

        return Teff, logg, Z, alpha




    
    def destructor(self):
        """
        PURPOSE: Destructor
        """
        string_1 = 'Invalid parameter space! Valid values are:'
        string_2 = '\nTeff: {0} \nlogg: {1} \nZ: {2} \nalpha: {3}\n'.format(
            np.sort(self.valid_t),
            np.sort(self.valid_g),
            np.sort(self.valid_z),
            np.sort(self.valid_a)
        )
        errorcode('error', string_1 + string_2)




        
    def getAtmosFITS(self, Teff, logg, Z, alpha):
        """
        PURPOSE: Download and load PHOENIX atmospheric SED models
        """

        # Create data folder if it do not exist

        dataDir = 'data/PHOENIX_AtmosFITS'
        pathlib.Path(self.path.joinpath(dataDir)).mkdir(parents=False, exist_ok=True)

        # Make sure code do not crash

        Teff, logg, Z, alpha = self.nearest_parameters(Teff, logg, Z, alpha, 'Atmos')

        # Accept only valid parameter space

        if Teff in self.valid_t and logg in self.valid_g \
           and Z in self.valid_z and alpha in self.valid_a:

            # Fetch zip url
            url_zip = (
                "ftp://phoenix.astro.physik.uni-goettingen.de/"
                "AtmosFITS/PHOENIX-ACES-AGSS-COND-2011_AtmosFITS_"
                "Z-{0:{1}2.1f}{2}{3}.zip".format(
                    Z,
                    "+" if Z > 0 else "-",
                    "" if alpha == 0 else ".Alpha=",
                    "" if alpha == 0 else "{:+2.2f}".format(alpha),
                )
            )

            # Select proper file name
            lte_file = (
                "lte{0:05}-{1:2.2f}-{2:{3}2.1f}{4}{5}."
                "PHOENIX-ACES-AGSS-COND-2011.ATMOS.fits".format(
                    Teff,
                    logg,
                    Z,
                    "+" if Z > 0 else "-",
                    "" if alpha == 0 else ".Alpha=",
                    "" if alpha == 0 else "{:+2.2f}".format(alpha),
                )
            )

            # Correct urls
            url_zip  = url_zip.replace("--", "-")
            lte_file = lte_file.replace("--", "-")

            # Download zip file
            phoenix_zip_file = self.path.joinpath(dataDir).joinpath(url_zip.split('/')[-1])
            if not phoenix_zip_file.is_file():
                print(f"Download Phoenix spectrum from {url_zip}")
                with urllib.request.urlopen(url_zip) as response, open(phoenix_zip_file, "wb") as out_file:
                    data = response.read()
                    out_file.write(data)

            # Extract only requested spectrum from zip
            phoenix_file = self.path.joinpath(dataDir).joinpath(lte_file)
            if not phoenix_file.is_file():
                print(f'Extracting file from zip {lte_file}')
                with zipfile.ZipFile(phoenix_zip_file, 'r') as zip:
                    # Extract all files
                    zip.extractall(self.path.joinpath(dataDir))
                    # Remove zip file
                    phoenix_zip_file.unlink(missing_ok=False)

            # Fetch spectrum
            spectrum = fits.getdata(phoenix_file)

        # Raise an error if parameters is outside limits
        else: self.destructor()

        # Finito!
        return spectrum




    

    def getHiResFITS(self, Teff, logg, Z, alpha):
        """
        PURPOSE: Download and load PHOENIX High Resolution SED models
        """

        # Prepare data structure
        dataDir = 'data/PHOENIX_HiResFITS'
        pathlib.Path(self.path.joinpath(dataDir)).mkdir(parents=False, exist_ok=True)
        phoenix_path = self.path.joinpath('data').joinpath('WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')

        # Accept only valid parameter space
        if Teff in self.valid_t and logg in self.valid_g \
           and Z in self.valid_z and alpha in self.valid_a:

            # FETCH WAVELENGTH DATA

            if not phoenix_path.is_file():
                print("Download PHOENIX wavelength file...")
                url = "ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
                with urllib.request.urlopen(url) as response, open(phoenix_path, "wb") as out_file:
                    data = response.read()
                    out_file.write(data)
            wvl = fits.getdata(phoenix_path)   # [AA]

            # FETCH SPECTRUM

            baseurl = (
                "ftp://phoenix.astro.physik.uni-goettingen.de/"
                "HiResFITS/PHOENIX-ACES-AGSS-COND-2011/"
                "Z-{0:{1}2.1f}{2}{3}/".format(
                    Z,
                    "+" if Z > 0 else "-",
                    "" if alpha == 0 else ".Alpha=",
                    "" if alpha == 0 else "{:+2.2f}".format(alpha),
                )
            )

            url = (
                baseurl + "lte{0:05}-{1:2.2f}-{2:{3}2.1f}{4}{5}."
                "PHOENIX-ACES-AGSS-COND-2011-HiRes.fits".format(
                    Teff,
                    logg,
                    Z,
                    "+" if Z > 0 else "-",
                    "" if alpha == 0 else ".Alpha=",
                    "" if alpha == 0 else "{:+2.2f}".format(alpha),
                )
            )

            # Fix url
            url = url.replace("--", "-")

            # Download spectrum
            spectrum_path = self.path.joinpath(dataDir).joinpath(url.split("/")[-1])
            if not spectrum_path.is_file():
                print(f"Download Phoenix spectrum from {url}")
                with urllib.request.urlopen(url) as response, open(spectrum_path, "wb") as out_file:
                    data = response.read()
                    out_file.write(data)

            # Fetch spectrum from file [ergs/s/cm^2/cm] -> convert with *0.1 to [uW/m^2/um]
            flux = fits.getdata(spectrum_path) / 1e8  # [ergs/s/cm2/AA]  

        # Raise an error if parameters is outside limits
        else: self.destructor()

        # Finito!
        return wvl, flux



    

    def getSpecIntFITS(self, Teff, logg, Z):
        """
        PURPOSE: Download and load PHOENIX spectral intensities
        """

        # Create data folder if it do not exist

        dataDir = 'data/PHOENIX_SpecIntFITS'
        pathlib.Path(self.path.joinpath(dataDir)).mkdir(parents=False, exist_ok=True)

        # Make sure that all paramters are valid

        Teff, logg, Z, alpha = self.nearest_parameters(Teff, logg, Z, 0, 'SpecInt')

        if Teff in self.valid_t and logg in self.valid_g and Z in self.valid_z:

            # DOWNLOAD SPECTRUM

            baseurl = (
                "ftp://phoenix.astro.physik.uni-goettingen.de/"
                "SpecIntFITS/PHOENIX-ACES-AGSS-COND-SPECINT-2011/"
                "Z{0}{1:2.1f}/".format(
                    "+" if Z > 0 else "-",
                    Z,
                )
            )

            url = (
                baseurl + "lte{0:05}{1}{2:2.2f}{3}{4:2.1f}."
                "PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits".format(
                    Teff,
                    "+" if logg == 0 else "-",
                    logg,
                    "+" if Z > 0 else "-",
                    Z
                )
            )

            # Fix url

            url = url.replace("--", "-")

            # Download spectrum

            phoenix_file = self.path.joinpath(dataDir).joinpath(url.split("/")[-1])

            if not phoenix_file.is_file():

                print(f"Download Phoenix spectrum from {url}")
                with urllib.request.urlopen(url) as response, open(phoenix_file, "wb") as out_file:
                    data = response.read()
                    out_file.write(data)

            # FETCH WAVELENGTH AND INTENSITIES

            data, hdr = fits.getdata(phoenix_file, header=True)
            mu  = fits.getdata(phoenix_file, 'MU')
            wvl = np.arange(hdr["CRVAL1"], hdr["CRVAL1"]+hdr["NAXIS1"]*hdr["CDELT1"], hdr["CDELT1"])

        # Raise an error if parameters is outside limits

        else: self.destructor()

        # Finito!

        return wvl, mu, data
