#!/usr/bin/env python3

"""
This script is an integrated part of PlatoSim's toolkit PLATOnium.
The usage consist of three main query methods:
  1) Create a customized PIC catalogue
  2) Create a single-star catalogue quering Simbad
  3) Create a PLATO pointing field catalogue for each camera FOV

General information:
  - One of the three query methods needs to be used to run this script.
  - All output files are directly useable by PALTOnium.
  - The arguments "--dmag" and "--dist" only apply to method 1 and 2.

1) Customized PIC catalogue (--pic):
  This option uses the PIC input catalog to create a smaller custimized
  stellar catalogue both for the PIC targets and contaminants. Moreover,
  this script can also be used to re-plot a old "--pic" catalog produced
  by this software. Usage example:

  $ picsim --pic 100 P1 LOPS2 --project <project_name> -p

2) Single-star cataloge (--simbad):
  This option can be used to feth a single target star by name, which
  is recognizable by the CDS Simbad query. By default is select stars
  within 45 arcsec (2 pixel) from the target star: Usage examples:

  $ picsim --simbad Mizar --project <project_name> -p

3) PLATO pointing field catalogue (--vizier):
  This option creates a coherent PLATO catalogue that covers the
  FoV of each camera group. It uses the Gaia DR3 catalogue and
  makes query of 16 grid points distribution around the requested
  PLATO pointing field. The output catalogues (one for each group)
  is recognized be by "platonium --fullframe" for the generation
  of full-frame CCD images. Usage examples:

  $ picsim --vizier LOPS2 --project <project_name> -p

"""

# Built-in
import os
import shutil
import random
import argparse
import datetime
import warnings
from pathlib import Path

# PlatoSim standard
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io.votable import parse
from astropy.coordinates import SkyCoord
from astropy import units as u
from tqdm import tqdm
from prettytable import PrettyTable

# PLATOnium extra
import ligo.skymap.plot

# PlatoSim functions
import platosim.plot      as pt
import platosim.utilities as ut
import platosim.starquery as sq
import platosim.referenceFrames as rf
from platosim.utilities  import errorcode
from platosim.simulation import Simulation


#==============================================================#
#                         BEGIN CLASS                          #
#==============================================================#


class PicSim(object):

    """Class to generate customised (PIC) catalogues.
    """
    
    def __init__(self, args):

        # I/O PARAMETERS

        # Plotting flag
        self.plot = args.plot

        # Verbosity (a.k.a log level) -> Identical to PlatoSim usage
        # verbose = 0: Cluster mode: Disabling print, warnings, and saving of logs
        # verbose = 1: Default mode: Print details to bash but do not save log files
        # verbose = 3: Debug mode  : Print details to bash and saves all log files
        if args.verbose == 0:
            self.verbose = 0            
            warnings.filterwarnings("ignore")
        elif args.verbose is None or args.verbose == 1:
            self.verbose = 1
            warnings.filterwarnings("ignore")
        else:
            self.verbose = 3

        # Output directory
        if (args.outdir is not None) or (args.project is not None):

            # NOTE "--outdir" overwrites "--project"
            if args.outdir is not None:
                self.outputDir = Path(args.outdir).resolve()
            elif args.project is not None:
                self.outputDir = Path(os.getenv('PLATO_WORKDIR')) / args.project / 'input'

            # Create directory if is doesn't exist
            self.outputDir.mkdir(parents=True, exist_ok=True)
        else:
            self.outputDir = args.outdir

            
        # GENERIC INPUT PARAMETERS
        
        # Contaminant magnitude limit
        if args.dmag is None:
            self.dmagConLimit = 5
        else:
            self.dmagConLimit = args.dmag

        # Warning for single star query
        if args.simbad and self.dmagConLimit > 21:
            errorcode('warning', 'Gaia is only complete to G < 21 mag!')

        # Contaminant distance limit (default 2 pixel)
        if args.dist is None:
            self.disConLimit = 45
        elif args.simbad or (args.pic is not None) and (args.dist in [30, 45, 60]):
            self.disConLimit = args.dist
        else:
            errorcode('warning', 'Not a valid contaminant-to-target distance! ' +
                      'Use {30, 45, 60} arcsec')
        

        

        
    #--------------------------------------------------------------#
    #                        PIC OF DESTINY                        #
    #--------------------------------------------------------------#            


    def printNotesPIC(self):

        """Function to show the notes about the PIC.
        """

        errorcode('message', '\nOverview notes for the PLATO Input Catalogue (PIC)')
        print("""
Notes on parsed argument "--pic":
  "star"   : Number of targets saved in catalog (Select "all" for all stars)
  "sample" : PIC samples [P1, P2, P4, P5, all]  (Select "all" for all samples)
  "field"  : PLATO pointing field [LOPS2, LOPN1, SPF, NPF]
        
Notes on PLATO fields:
  Number of stars in PIC200 LOPS2: (t:179,564, c:105,865,120)
  Number of stars in PIC200 LOPN1: (t:175,597, c:112,917,958)
  Number of stars in PIC110 SPF  : (t:163,772, c: )
  Number of stars in PIC110 NPF  : (t:156,971, c: )

Notes on Sample flag:
  P1: LOPS2:   8,835, LOPN1:   9,367, SPF:   6,817, NPF:   6,892
  P2: LOPS2:     678, LOPN1:     680, SPF:     717, NPF:     668
  P4: LOPS2:  12,026, LOPN1:  12,026, SPF:  16,866, NPF:  16,166
  P5: LOPS2: 157,543, LOPN1: 152,818, SPF: 132,571, NPF: 140,046
  Note that P2 is a bright sub-sample of P1, and P1 of P5.

Notes on PIC catalogue creation:
    PIC200:
      - Magnitude completeness of contaminants    : Pmag < 17 mag
      - Maximum radial contaminat-target distance : Rmax < 45 arcsec
      - Calibration methods: The Gaia DR3 colour information (BP-RP)
        and extinction maps are sued to derive the PLATO magnitudes.
    PIC110:
      - Magnitude completeness of contaminants    : Vmag < 17 mag
      - Maximum radial contaminat-target distance : Rmax < 60 arcsec
      - Calibration method using the de-reddened Gaia V magnitude obtained
        from Gaia photometry using a calibration technique being valid for: 
        0.5 < Bp-Rp < 5.17
        """)




    
    def initPIC(self):

        """Initialise the PIC input parameters.
        """

        if self.verbose > 0:
            errorcode('software', '\nPIC of Destiny')

        # Optinal parameters without checks
        self.saveAscii     = args.save
        self.inputFiles    = args.incat
        self.oldCatalogue  = args.unique
        self.numTimeseries = args.ntime

        
        # MANDATORY PARAMETERS
            
        # Mandatory parameters
        self.stars, self.sample, self.field = args.pic[0]
            
        # Check if sample is correct
        if self.sample == 'all':
            self.sample = None
        elif self.sample in ['P1', 'P2', 'P4', 'P5']:
            self.sample = self.sample
        else:
            errorcode('error', 'Not a valid PIC sample! Use {P1, P2, P4, P5, all}')

        # Check the PLATO field
        if self.field in ['SPF', 'NPF']:
            self.pic = 'PIC110'
            if self.field == 'SPF': self.numPIC = 137052
            if self.field == 'NPF': self.numPIC = 144507
        elif self.field in ['LOPS2', 'LOPN1']:
            self.pic  = 'PIC200'
            if self.field == 'LOPS2': self.numPIC = 179564
            if self.field == 'LOPN1': self.numPIC = 175325 # Original: 175,597            
        else:
            errorcode('error', 'Not valid pointing! Use {LOPS2, LOPN1, SPF, NPF}')        

        # Check number of stars
        self.numTargets = self.stars
        if self.numTargets[0].isdigit():
            self.numTargets = int(self.numTargets)
            if self.numTargets > self.numPIC:
                errorcode('error', 'More stars selected than available in the PIC! ' +
                          'See -h for PIC notes')
        elif self.stars == 'all':
            if self.numTimeseries is not None:
                self.numTargets = int(self.numTimeseries/6.) + 1
            else:
                self.numTargets = self.numPIC


        # OBSERVATIONAL PARAMETERS
        
        # Select magnitude range or select all stars
        if args.mag is None:
            self.magRange = [0, 21]
        else:
            self.magRange = ut.convertMagnitudeRange(magRange)
        
        # Visibility by N-cams
        if args.ncams in [None, 6, 12, 18, 24]:
            self.ncams = args.ncams
        else:
            errorcode('Not valid N-CAM visibility! Use {6, 12, 18, 24}')

        # Check cam-group entry
        if args.group in [None, 1, 2, 3, 4]:
            self.group = args.group
        else:
            errorcode('error', 'Not valid group number! Use {1, 2, 3, 4}')

        # Check evolutionary stage 
        if args.lum in [None, 'V', 'IV']:
            self.lumClass = args.lum
        else:
            errorcode('error', 'Not valid evolutionary stage! Use {V, IV}')

        # Check spectral type
        if args.spec in [None, 'F', 'G', 'K']:
            self.specType = args.spec
        else:
            errorcode('error', 'Not valid spectral type! Use {F, G, K}')
                                

        # EXTRA I/O PARAMETERS
        
        # Input directory        
        self.inputDir = Path(os.getenv("PLATO_PROJECT_HOME")) / 'inputfiles/data_picsim'    
        if not self.inputDir.is_dir():
            errorcode('message', 'Inuaguration: Welcome to PIC of Destiny!')

        # Name space for output files
        # Default is to use the latest PIC catalog saved
        # Else several old catalogs can be parsed for replotting
        # E.g. use: --intar starcat**targets.txt --incon starcat**contaminants.txt
        if self.outputDir:

            # Prefix name of output files
            self.outputPrefix = 'starcat'
            if self.sample: self.outputPrefix += '_'      + self.sample
            if self.field:  self.outputPrefix += '_'      + self.field
            if self.group:  self.outputPrefix += '_Group' + str(self.group)
            if self.ncams:  self.outputPrefix += '_Ncam'  + str(self.ncams)
            if self.inputFiles is not None: self.outputPrefix += '_NewCat'
            self.outputPrefixTar = self.outputPrefix + '_targets.ftr'
            self.outputPrefixCon = self.outputPrefix + '_contaminants.ftr'

            # Save either a PlatoSim (ascii) or a PLATOnium (feather) calogue
            if self.saveAscii:
                self.outputFileCat = self.outputDir / f'{self.outputPrefix}.txt'
            else:
                self.outputFileTar = self.outputDir / self.outputPrefixTar
                self.outputFileCon = self.outputDir / self.outputPrefixCon
                

        # PLOT PARAMETERS
            
        # Generic title to use for plotting
        self.title = f'{self.pic}, {self.field}, {self.sample} sample'
        
        # Add latex font if catalogue is saved
        if self.outputDir is None:
            from platosim.matplotlibrc import setup
            setup()
        else:
            from platosim.matplotlibrc import latex
            latex()
         

            

            
    def loadPIC(self):

        """Fetch PIC targets from feather files.
        """
        
        # Files names
        inputFileTar = self.inputDir / f'{self.pic}_{self.field}_targets.ftr'
        inputFileCon = self.inputDir / f'{self.pic}_{self.field}_contaminants.ftr'

        # Fetch PIC catalogue from FTP server
        if not inputFileTar.is_file() or not inputFileCon.is_file():
            errorcode('message', '\nInuaguration! Welcome to PIC of Destiny!')
            print(f'Downloading {self.pic} catalogue..')
            ut.downloadFromFTP(inputFileTar.name, self.inputDir, 'plato')
            ut.downloadFromFTP(inputFileCon.name, self.inputDir, 'plato')

        # Load catalogues
        if self.verbose > 0:
            print('Loading stellar catalogues..')

        self.df0 = pd.read_feather(inputFileTar)
        self.dc0 = pd.read_feather(inputFileCon)
        

    

        
    def loadOldPIC(self):

        """Load old PIC catalog to re-plot it.
        """

        # Check if one or more catalogs are parsed
        # NOTE using *.ftr will order in alphabetic order
        if len(self.inputFiles) > 2:
            inputFileTar = self.inputFiles[1::2]
            inputFileCon = self.inputFiles[0::2]
        else:
            inputFileTar = [self.inputFiles[1]]
            inputFileCon = [self.inputFiles[0]]

        # Append each catalog if several are parsed
        N = int(len(self.inputFiles)/2.)
        for i in range(N):
            if i == 0:
                self.df0 = pd.read_feather(inputFileTar[i])
                self.dc0 = pd.read_feather(inputFileCon[i])
            else:
                self.df0 = pd.concat(self.df, pd.read_feather(inputFileTar[i]))
                self.dc0 = pd.concat(self.dc, pd.read_feather(inputFileCon[i]))




                
    def getStellarClass(self, df):

        """Classifier of spectral type and luminosity class.
        """
        
        # Check stellar sample
        if self.sample is not None:
            df = df[df['sample'] == self.sample]

        # Seperate dwarf (MS) and sub-gaint (post MS) stars
        ds = df[df.R < ut.getMainSequenceLimit(df.Teff)]
        sg = df[df.R > ut.getMainSequenceLimit(df.Teff)]
        
        # Seperate dwarf spectral types
        dK = ds[ ds.Teff < 5300]
        dG = ds[(ds.Teff > 5300) & (ds.Teff < 5900)]
        dF = ds[ ds.Teff > 5900]
        
        # Seperate sub-giants spectral types
        sgK = sg[ sg.Teff < 5300]
        sgG = sg[(sg.Teff > 5300) & (sg.Teff < 5900)]
        sgF = sg[ sg.Teff > 5900]

        return df, ds, sg, dK, dG, dF, sgK, sgG, sgF



                

    def queryTargetsPIC(self):

        """Fetch PIC targets from feather files.
        """
        if self.verbose > 0:
            errorcode('module', '\nPIC targets')
        
        # Check parsing of old catalogue to select new unique targets
        if self.oldCatalogue:
            try:
                df_old = pd.read_feather(self.oldCatalogue)
            except:
                errorcode('error', 'Old picsim target catalogue do not exist!')
            else:
                cond = df['PIC'].isin(df_old['PIC']) # Boolen
                df   = df.drop(df[cond].index)
        else:
            df = self.df0

            
        # QUERY CUTS
            
        # Fetch stellar classifications
        df, ds, sg, dK, dG, dF, sgK, sgG, sgF = self.getStellarClass(df)
                
        # Check evolutionary stage 
        if self.lumClass == 'V':  df = ds
        if self.lumClass == 'IV': df = sg

        # Check spectral type
        if self.specType is not None:
            if self.specType == 'K':
                df = dK if self.lumClass == 'V' else (sgK if self.lumClass == 'IV' else pd.concat([dK, sgK]))
            elif self.specType == 'G':
                df = dG if self.lumClass == 'V' else (sgG if self.lumClass == 'IV' else pd.concat([dG, sgG]))
            elif self.specType == 'F':
                df = dF if self.lumClass == 'V' else (sgF if self.lumClass == 'IV' else pd.concat([dF, sgF]))

        # Select stars depending on camera visibility
        # NOTE this will be done properly later but here it is done to catch
        #      the failure when the stars requested are larger than catalogue
        if self.ncams is not None:
            df = df[df['ncams'] == int(self.ncams)]

        # Allow to select only one cam-group
        if self.group:
            sim = Simulation('picsim')
            # Fetch PLATO pointing field [deg]
            alpha, delta, kappa = ut.getPointingField(self.field, unit='deg')
            dex = sim.getStarsWithinCameraGroup(df.ra.to_numpy(), df.dec.to_numpy(),
                                                alpha, delta, kappa,
                                                self.group)

            # Select max number of stars if too many is requested
            max_stars = np.sum(dex[0])
            if self.numTargets > max_stars:
                self.numTargets = max_stars
            df = df[dex[0]]

        # Check P passband magnitude range
        if self.magRange is not None:
            df = df[(df['mag'] > self.magRange[0]) & (df['mag'] < self.magRange[1])]

        # Check if too many stars are selected
        if len(df) < self.numTargets:
            if (self.stars != 'all'):
                errorcode('warning', f'Only {len(df)} stars are available after cuts!')
            self.numTargets = len(df)
            

        # RANDOM CATALOGUE SELECTION

        # Only after the above cuts can we select randomly a sub-sample (if applicable)
        # Use random.choices(data, weights=some weight, k=numTargets) if a weigthing is needed
        # NOTE args.stars being a string is needed to select "all" in one of the samples too
        # NOTE we need to shuffle the sample if a specific number of timeseries are needed
        if (self.stars != 'all') or (self.numTimeseries is not None):
            df = df.sample(n = self.numTargets)

        # Redefine the number fo targets
        numTargets = len(df['PIC'])

        # Add ID column to df
        df  = df.reset_index()
        df  = df.drop(columns=['index'])
        ids = np.arange(1, len(df.ra)+1)
        IDs = [f'{i}'.zfill(9) for i in ids]
        df  = ut.pdAddColumn(df, IDs, 'ID')

        # Randomly choose sample until the number of timeseries are met
        if self.numTimeseries:
            numImagettes = 0
            for i in range(numTargets):
                # Add camera visibilities
                numImagettes += df['ncams'].to_numpy()[i]
                # The data it cut by the index that reached the limit
                if numImagettes >= self.numTimeseries:
                    print(i)
                    df = df.drop(df.index[:i])
                    break

                
        # MAKE AN OVERVIEW TABLE

        cameras = [6, 12, 18, 24]
        countStars = [df[df['ncams'] == 6].count()[0],  df[df['ncams'] == 12].count()[0],
                      df[df['ncams'] == 18].count()[0], df[df['ncams'] == 24].count()[0]]
        countImagettes = [countStars[i] * cameras[i] for i in range(4)]

        # Last collection of stars
        #if fileFormat == '.txt': numTargets = len(ncams)
        numTargets = np.sum(countStars)

        # Make table
        self.t = PrettyTable(['N-cams', 'Targets', 'Time series'])
        for i in range(4):
            self.t.add_row([str(cameras[i]), countStars[i], countImagettes[i]])
        if self.verbose > 0:
            print(f'User catalog contains {numTargets} stars of ' +
                  f'{np.sum(countImagettes)} time series drawn from:')
            print(self.t)

        # Store copy of modified df
        self.df = df

        



    def queryContaminantsPIC(self):
        
        """Fetch and select PIC contaminants.

        Function to select stellar contaminants belonging to specific 
        PIC target stars within a certain radial distance. Note that
        there are more than 8 million contaminants stars in the original
        catalog, hence, this step can take some time if the large
        catalog are to be used.
        """        
        if self.verbose > 0:
            errorcode('module', '\nPIC contaminants')
            print(f'Fetching contaminants within {self.disConLimit} arcsec ' +
                  f'and {self.dmagConLimit} mag from each target:')
            
        df = self.df
        dc = self.dc0

        # Create empty data-frame to append to 
        dc1 = pd.DataFrame()
        numConPerTar = np.zeros(self.numTargets)

        for starNo in tqdm(range(self.numTargets), bar_format=ut.tqdmBar()):

            # Find corresponding contaminants to our target
            dcc = dc[dc['PIC'] == df['PIC'].iloc[starNo]]

            # We will only add stars within a user-defined distance of our target
            dcc = dcc[dcc['dis'] < self.disConLimit]

            # Add only contaminants brigther than threshold
            dcc = dcc[dcc['mag'] < df['mag'].iloc[starNo] + self.dmagConLimit]

            # Add contaminats to contaminant list
            if dc1.empty and not dcc.empty:
                dc1 = dcc                
            elif not dc1.empty and not dcc.empty:
                dc1 = pd.concat([dc1, dcc])

            # Store number of contaminats per target for statistics below
            numConPerTar[starNo] = len(dcc)

        # Add the number of contaminats to the target catalog
        self.df['ncon'] = numConPerTar
        self.dc = dc1

    

        
            
    def plotTargetsPIC(self):

        """Plot function for PIC targets.
        """
        
        df = self.df
        
        # Plot targets in aitoff Galactic sky projection
        self.fig0, ax = pt.drawStarsInSkyAitoff(df.ra, df.dec, df.mag)
        if self.plot: plt.show()

        # Plot Zoom-in on FOV
        if df.shape[0] > 1000: mag = None
        else: mag = df.mag.to_numpy()
        self.fig1, ax = pt.plotPlatoFOV(self.field, df.ra, df.dec, magStars=mag,
                                        showGroups=False, title=self.title)
        if self.plot: plt.show()

        # Plot sample distribution in Teff vs. Radius
        _, ds, sg, dK, dG, dF, sgK, sgG, sgF = self.getStellarClass(self.df0)
        self.fig2, ax = pt.plotTeffvsRadius(ds, dK, dG, dF, sg, sgK, sgG, sgF, df,
                                            self.title)
        if self.plot: plt.show()



        
        
    def plotContaminantsPIC(self):

        """Plot function for PIC contaminants.
        """
        
        df = self.df
        dc = self.dc
        
        # Necessary to define again in case a PIC sample is defined
        magRange = [df.mag.min(), df.mag.max()]
        self.fig3 = plt.subplots(2, 2, figsize=(14,8))
        pt.plotStellarSampleDistributions(self.fig3, magRange,
                                          df.mag.to_numpy(),  dc.mag.to_numpy(),
                                          df.ncon.to_numpy(), dc.dis.to_numpy())
        if self.plot: plt.show()




        
    def prologuePIC(self):

        if self.verbose > 0:
            errorcode('module', '\nPrologue')
        
        if self.outputDir is not None:

            # Copy the YAML file to the project if it doesn't exist
            if self.verbose > 0:
                print(f"Copying YAML configuration file")
            ut.copyInputYAML(self.field, self.outputDir)

            # Save textfile with information about PIC stars
            log = ('PIC Catalogue include\n' +
                   '----------------------------------\n' +
                   f'Catalogue             : {self.pic}\n' +
                   f'PLATO field           : {self.field}\n' +
                   f'PLATO sample          : {self.sample}\n' +
                   f'PLATO camera-group    : {self.group}\n' +
                   f'Magnitude range       : {self.magRange[0]:.2f}-{self.magRange[1]:.2f}\n' +
                   f'Luminosity class      : {self.lumClass}\n' +
                   f'Spectral type         : {self.specType}\n' +
                   f'Contaminant distance  : {self.disConLimit} arcsec\n' +
                   f'Contaminant magnitude : {self.dmagConLimit} mag\n' +
                   self.t.get_string() + '\n')
            with open(self.outputDir / f'{self.outputPrefix}.log','w') as file:
                file.write(log)

            # Save figures
            resolution = 200
            if self.verbose > 0:
                print('Saving all plots to {0}'.format(self.outputDir))
            self.fig0.savefig(self.outputDir / f'plot_{self.outputPrefix}_allsky.png',
                              bbox_inches='tight', dpi=resolution)
            self.fig1.savefig(self.outputDir / f'plot_{self.outputPrefix}_pointing.png',
                              bbox_inches='tight', dpi=resolution)
            self.fig2.savefig(self.outputDir / f'plot_{self.outputPrefix}_TeffvsRadius.png',
                              bbox_inches='tight', dpi=resolution)
            if not self.fig3 is False:
                self.fig3[0].savefig(self.outputDir / f'plot_{self.outputPrefix}_distribution.png',
                                     bbox_inches='tight', dpi=resolution)            
                
            # Save ascii catalog (PlatoSim) or feathers (PLATOnium)
            if self.saveAscii:
                if self.verbose > 0:
                    print(f'Saving file {self.outputFileCat}')
                df0        = pd.concat([self.df.ra,  self.dc.ra])
                df0['dec'] = pd.concat([self.df.dec, self.dc.dec])
                df0['mag'] = pd.concat([self.df.mag, self.dc.mag])
                df0.to_csv(self.outputFileCat, sep=' ', header=False, float_format='%.6f')
            
            else:
                # We reset the index in order to save to feather
                df = self.df.reset_index(drop=True)
                dc = self.dc.reset_index(drop=True)
                if self.verbose > 0:
                    print(f'Saving file {self.outputFileTar}')
                df.to_feather(self.outputFileTar)
                if self.verbose > 0:
                    print(f'Saving file {self.outputFileCon}')
                dc.to_feather(self.outputFileCon)

                # Output file name of HPC data
                outputFileVarSim = f'{self.outputDir}/cluster_{self.outputPrefix}.data'
                if self.verbose > 0:
                    print(f'Saving file {outputFileVarSim}')
                headerVar = 'ID,PIC,M,R,Teff'
                df_hpc = pd.DataFrame()
                df_hpc['ID']   = df.ID
                df_hpc['PIC']  = df.PIC
                df_hpc['Teff'] = df.Teff
                df_hpc['R']    = df.R
                df_hpc['M']    = df.M
                df_hpc = df_hpc.reset_index(drop=True)
                df_hpc.to_csv(outputFileVarSim, sep=',', index=False)
                              #float_format=['%s', '%i', '%i', '%0.3f', '%0.3f'])
                


        
    def cameraObservability(self):
        
        # CAMERA OBSERVABILITY

        # Since the PIC are bigger than actual FOV (18.89 deg radius) we need to
        # Loop over each camera-group and avoid stars that are close to the FOV edge
        # NOTE the following also checks if the subfield can be placed on a CCD
        group = False
        if group:

            # CONFIGURE SPACECRAFT

            # # Set the telescope group ID
            # sim["Telescope/GroupID"] = camGroup        

            # # Solar panel orientation: 0, 90, 180, and 270 degrees for Q1, Q2, Q3, and Q4
            # solarPanelOrientation = sim["Platform/SolarPanelOrientation"] = math.fmod(self.quarter * 90., 360.)
            # self.solarPanelOrientation = np.deg2rad(float(solarPanelOrientation))

            # # Select PLATO pointing field from inputfile (e.g. SPF, NPF, etc.)
            # pointingField  = sim["ObservingParameters/StarCatalogFile"]
            # raPlatformDeg  = sim["ObservingParameters/RApointing"]  = alpha
            # decPlatformDeg = sim["ObservingParameters/DecPointing"] = delta

            # # Camera-group Alt (tilt) and Az
            # tiltTelescope    = sim["CameraGroups/TiltAngle"][self.group-1]
            # azimuthTelescope = sim["CameraGroups/AzimuthAngle"][self.group-1]

            # # Include spacecraft Pointing Repeatability Error (PRE) between consecutive quarters
            # # NOTE: Will only be included if the file "PRE.txt" is available in the input folder 
            # if os.path.exists(self.inputDir + '/PRE.txt'):
            #     PRE = np.loadtxt(self.inputDir + '/PRE.txt')
            #     dex = np.where(PRE[:,0] == self.quarter)[0]
            #     raPlatformDeg         += PRE[dex, 1][0]
            #     decPlatformDeg        += PRE[dex, 2][0]
            #     solarPanelOrientation += np.deg2rad(PRE[dex, 3][0])

            # # Include Absolute Pointing Error (APE) due to camera misalignments
            # # NOTE: Will only be included if the file "APE.txt" is available in the input folder
            # if os.path.exists(self.inputDir + '/APE.txt'):
            #     APE = np.loadtxt(self.inputDir + '/APE.txt')
            #     dex = (self.group - 1) * 6 + self.camera - 1
            #     tiltTelescope    += APE[dex, 0]
            #     azimuthTelescope += APE[dex, 1]

            print('Checking star visibility for each camera group..')
            dexGroup1, d1 = getStarsWithinCamGroup(1, raPF, decPF, df['ra'].to_numpy(), df['dec'].to_numpy())
            dexGroup2, d2 = getStarsWithinCamGroup(2, raPF, decPF, df['ra'].to_numpy(), df['dec'].to_numpy())
            dexGroup3, d3 = getStarsWithinCamGroup(3, raPF, decPF, df['ra'].to_numpy(), df['dec'].to_numpy())
            dexGroup4, d4 = getStarsWithinCamGroup(4, raPF, decPF, df['ra'].to_numpy(), df['dec'].to_numpy())

            # NOTE check how many stars that are observable in each group
            #-----------------
            # starGroup = np.sum([dexGroup1*1, dexGroup2*1, dexGroup3*1, dexGroup4*1], axis=1)
            # print(starGroup)
            #-----------------

            # NOTE the following is a small table for checking radial distance to OA
            # A ghost within 8 deg from the OA will produce a point-like ghost
            #----------------
            # Pmag = ut.passbandConversionV2P(mag, Teff)
            # t2 = PrettyTable(['PIC', 'V', 'P', 'd_G1 (deg)', 'd_G2 (deg)', 'd_G3 (deg)', 'd_G4 (deg)'])
            # for i in range(len(ID)):
            #     t2.add_row([int(ID[i]), '{:.2f}'.format(mag[i]), '{:.2f}'.format(Pmag[i]), '{:.2f}'.format(d1[i]),
            #                 '{:.2f}'.format(d2[i]), '{:.2f}'.format(d3[i]), '{:.2f}'.format(d4[i])])
            # print('\nUser have requested data table of all targets')
            # print(t2); exit()
            #----------------

            # Find the actual N-Cam visibility
            ncams = (dexGroup1*1 + dexGroup2*1 + dexGroup3*1 + dexGroup4*1) * 6

            # Remove stars that are not observable at all
            dex = dexGroup1 + dexGroup2 + dexGroup3 + dexGroup4
            PIC   = PIC[dex].astype(int)
            ra    = ra[dex]
            dec   = dec[dex]
            mag   = mag[dex]
            Teff  = Teff[dex]
            R     = R[dex]
            M     = M[dex]
            ncams = ncams[dex]

            # Again choose stars from N-Cam visibility
            if ncams:
                dex = ncams == int(ncams)
                PIC   = PIC[dex].astype(int)
                ra    = ra[dex]
                dec   = dec[dex]
                mag   = mag[dex]
                Teff  = Teff[dex]
                R     = R[dex]
                M     = M[dex]
                ncams = ncams[dex]

            # Allow to select only one cam-group
            if group:
                dex = starsWithinCamGroup(group, raPF, decPF, ra, dec)
                PIC   = PIC[dex].astype(int)
                ra    = ra[dex]
                dec   = dec[dex]
                mag   = mag[dex]
                Teff  = Teff[dex]
                R     = R[dex]
                M     = M[dex]
                ncams = ncams[dex]




    #--------------------------------------------------------------#
    #                   SINGLE GAIA STAR QUERY                     #
    #--------------------------------------------------------------#            

    
    def initSimbad(self):

        """Initialise the Simbad input parameters.
        """

        if self.verbose > 0:
            errorcode('software', '\nSimbad target query')
        
        # Arguments for GaiaDR3 star query
        self.simbad = args.simbad




        
    def querySimbad(self):

        """Function to query a Gaia DR3 star and its contaminants.
        """

        # Query Gaia DR3 star
        df = sq.simbadQuery(self.simbad, radius=self.disConLimit)
        if self.verbose > 0:
            print(f'\nCatalogue for {self.simbad}:')
            print(df)
            
        # Save output catalogue
        if self.outputDir is not None:
            da = df.loc[:, ['ra', 'dec', 'Pmag']]
            da.to_csv(self.outputDir / f"starcat_{self.simbad}.txt", sep=' ', header=False)
            df.to_feather(self.outputDir / f"starcat_{self.simbad}.ftr")


    
    #--------------------------------------------------------------#
    #                       PLATO CAMERA FOV                       #
    #--------------------------------------------------------------#            


    def initVizier(self):

        """Initialise the Simbad input parameters.
        """

        if self.verbose > 0:
            errorcode('software', '\nVizier PLATO FOV query')

        # Arguments for GaiaDR3 star query        
        self.field  = args.vizier
        self.bright = args.yale_stars
        self.astro  = args.gaia_astro

        # Magnitude limits
        if args.magmin is None:
            self.magmin = 0
        else:
            self.magmin = args.magmin
            
        if args.magmax is None:
            self.magmax = 15
        else:
            self.magmax = args.magmax
            
        # Check if output folder exist
        if self.outputDir:
            self.filename = self.outputDir / f'starcat_GaiaDR3_{self.field}'
        else:
            errorcode('error', 'Option --vizier needs a output destination!')

        # Check for inconsitent query
        if self.bright and self.astro:
            errorcode('error', 'Flags "yale_stars" and "gaia_astro" cannot be parsed at once!')
            
        # Constants [deg]
        self.rGroup = 19.0                       # Max radius to query stars within a group
        self.cGroup = 9.2                        # Diagonal opening angle of groups
        self.aGroup = np.sqrt(self.cGroup**2/2)  # Horizontal/vertical opening angle of group
        
        # Get pointing of platform [rad] 
        self.alpha, self.delta, self.kappa = ut.getPointingField(self.field, unit='rad')

        # Find the camera group pointings [rad]
        raGroups, decGroups = rf.getCameraGroupCoordinates(self.alpha, self.delta, self.kappa)
        self.raGroups  = np.rad2deg(np.append(raGroups,  self.alpha))
        self.decGroups = np.rad2deg(np.append(decGroups, self.delta))

        # Make a grid in azimuth and tilt angles and find the sky coordinates [deg]

        # Grid constants [deg]
        r = 7.3                    # Radius of grid point search
        a = 10.2                   # Distance between equidistant grid points
        c = np.sqrt(2*a**2)        # Diagonal distance of grid points
        b = np.sqrt(a**2+(2*a)**2) # Semi-diagonal distance of grid points
              
        gridCC = [(-2*a, +2*a), (-1*a, +2*a), (0, +2*a), (+1*a, +2*a), (+2*a, +2*a),
                  (-2*a, +1*a), (-1*a, +1*a), (0, +1*a), (+1*a, +1*a), (+2*a, +1*a),
                  (-2*a, +0*a), (-1*a, +0*a), (0, +0*a), (+1*a, +0*a), (+2*a, +0*a),
                  (-2*a, -1*a), (-1*a, -1*a), (0, -1*a), (+1*a, -1*a), (+2*a, -1*a),
                  (-2*a, -2*a), (-1*a, -2*a), (0, -2*a), (+1*a, -2*a), (+2*a, -2*a)]
        
        gridEQ = [( -45.0, 2*c), ( -22.5, b), (  0, 2*a), ( 67.5, b), ( 45.5, 2*c),
                  ( -67.5,   b), ( -45.0, c), (  0,   a), ( 45.0, c), ( 22.5,   b),
                  ( -90.0, 2*a), ( -90.0, a), (  0,   0), ( 90.0, a), ( 90.0, 2*a),
                  (-112.5,   b), (-135.0, c), (180,   a), (135.0, c), (112.5,   b), 
                  (-135.0, 2*c), (-157.5, b), (180, 2*a), (157.5, b), (135.0, 2*c)]

        # Get grid points in sky coordinates [deg]   
        x = np.zeros(len(gridEQ))
        y = np.zeros(len(gridEQ))
        for i in range(len(gridEQ)):
            x[i], y[i] = rf.platformToTelescopePointingCoordinates(self.alpha,
                                                                   self.delta,
                                                                   self.kappa,
                                                                   np.deg2rad(gridEQ[i][0]),
                                                                   np.deg2rad(gridEQ[i][1]))
        self.raGrid  = np.rad2deg(x)
        self.decGrid = np.rad2deg(y)

        # Shorten names
        self.r, self.a, self.c = r, a, c
        self.gridCC, self.gridEQ = gridCC, gridEQ

        # Plot grid used for query
        if self.plot:
            self.plotVizier()
        
        



    def plotVizier(self):

        # PLOT GRID IN CARTESIAN COORDINATES
        
        # Shorten names
        r, a, c = self.r, self.a, self.c
        rg, ag  = self.rGroup, self.aGroup
        gridCC = self.gridCC
        gridEQ = self.gridEQ
        
        # Plot the grid in cartesian coordinates
        fig, ax = plt.subplots(figsize=(7,7))
        for i in range(len(gridCC)):
            ax.add_artist(plt.Circle((gridCC[i][0], gridCC[i][1]), r, color='k', alpha=.2))
        ax.add_artist(plt.Circle(( ag, ag), rg, color='b', alpha=.2))
        ax.add_artist(plt.Circle(( ag,-ag), rg, color='b', alpha=.2))
        ax.add_artist(plt.Circle((-ag,-ag), rg, color='b', alpha=.2))
        ax.add_artist(plt.Circle((-ag, ag), rg, color='b', alpha=.2))

        # Settings
        ax.set_title('Cartesian query grid')
        ax.set_xlabel('x [deg]')
        ax.set_ylabel('y [deg]')
        ax.set_xlim(-30, 30)
        ax.set_ylim(-30, 30)
        ax.set_aspect('equal')
        plt.show()

        # PLOT GRID IN EQUATORIAL COORDINATES        
        
        fig = plt.figure(figsize=(7,7))

        # Plot the grid points on the sky
        platform = SkyCoord(np.rad2deg(self.alpha), np.rad2deg(self.delta),
                            frame='icrs', unit=u.deg)
        ax = plt.axes(projection='astro degrees zoom', center=platform,
                      radius='35 deg', rotate='180 deg')

        # Plot pointing of platform
        ax.plot(platform.ra.deg, platform.dec.deg, '*', c='k', mfc='magenta', ms=25,
                    transform=ax.get_transform('world'))

        # Plot the pointing of each camera group
        colors = ['b', 'limegreen', 'yellow', 'r']
        for i, c in zip(range(4), colors):
            ax.plot(self.raGroups[i], self.decGroups[i], 'o', ms=13, c=c, mec='k',
                    transform=ax.get_transform('world'))

        # Plot grid points
        grid = SkyCoord(self.raGrid, self.decGrid, frame='icrs', unit=u.deg)
        ax.plot(grid.ra.deg, grid.dec.deg, 'o', ms=10, c='k',mec='k',
                transform=ax.get_transform('world'))      

        # Settings
        ax.scalebar((0.05, 0.05), 10 * u.deg).label()
        ax.compass(0.95, 0.05, 0.1)
        ax.grid(color='gray')
        ax.set_xlabel('RA [deg]')
        ax.set_ylabel('Dec [deg]')
        plt.show()




        
    def queryVizier(self):

        """Query function for Gaia DR3.

        This function query a circular area from the Gaia DR3 given a
        equatorial pointing of the PLATO spacecraft. It uses a grid of 
        9 circular regions to search for stars (in order not to exceed
        the RAM memory) and concatenates these grids into a final Gaia
        custom catalogue spanning each of the camera groups.

        Authors: Juan Cabrera & Nicholas Jannsen
        https://www.cosmos.esa.int/web/gaia-users/archive/programmatic-access
        """
                
        if self.verbose > 0:
            print('\nStart Gaia DR3 query using magnitude limits: ' +
                  f'{self.magmin} - {self.magmax}')
            
        # Query stars within the FOV of each grid
        for i in tqdm(range(len(self.raGrid)), bar_format=ut.tqdmBar()):
            
            df0 = sq.gaiaRegionQuery(self.raGrid[i], self.decGrid[i], radius=self.r,
                                     maglim_min=self.magmin, maglim_max=self.magmax,
                                     flag_astro=self.astro, ofile=f'{self.filename}.vot')

            # Contatenate catalogue
            if i == 0: df = df0
            else:      df = pd.concat([df, df0])
        
        # Remove duplicate stars (from overlapping grid)
        df = df.drop_duplicates(subset=['gaiaDR3'])

        # Replace missing Gaia colors assuming M0 dwarfs
        if df.BP_RP.isna().sum() > 0:
            df.BP_RP[df.BP_RP.isna()] = 2.0
            
        # Convert Gmag to Pmag
        df['Pmag']  = ut.passbandConversionG2P(df.Gmag, df.BP_RP)
        df['PBmag'] = ut.passbandConversionG2P(df.Gmag, df.BP_RP, camera='fast_blue')
        df['PRmag'] = ut.passbandConversionG2P(df.Gmag, df.BP_RP, camera='fast_red')
            
        # If requested, add bright stars not available in the Gaia catalogue
        if self.bright:
            sirius  = {'gaiaDR3':'Sirius', 'ra':101.2871667, 'dec':-16.7161167,
                       'Pmag':  ut.passbandConversionV2P(-1.46, 9940),
                       'PBmag': ut.passbandConversionV2P(-1.46, 9940),
                       'PRmag': ut.passbandConversionV2P(-1.46, 9940)}
            canopus = {'gaiaDR3':'Canopus', 'ra': 95.9879167, 'dec':-52.6956611,
                       'Pmag':  ut.passbandConversionV2P(-0.72, 7400),
                       'PBmag': ut.passbandConversionV2P(-0.59, 7400),
                       'PRmag': ut.passbandConversionV2P(-0.96, 7400)}
            epscma  = {'gaiaDR3':'epsCMa', 'ra':104.6564583, 'dec':-28.9720861,
                       'Pmag':  ut.passbandConversionV2P(1.50, 22900),
                       'PBmag': ut.passbandConversionV2P(1.29, 22900),
                       'PRmag': ut.passbandConversionV2P(1.59, 22900)}
            df_sirius  = pd.DataFrame([sirius])
            df_canopus = pd.DataFrame([canopus])
            df_epscma  = pd.DataFrame([epscma])
            df = pd.concat([df, df_sirius, df_canopus, df_epscma])

        # Keep only stars within the camera group FOV
        if self.verbose > 0:
            print(f'\nCreating catalogue for each camera group')

        for i in range(5):

            # Calculate angular distance [deg]
            dOA = ut.radialDistance(self.raGroups[i], self.decGroups[i],
                                    df.ra.to_numpy(), df.dec.to_numpy())
            df0 = df[dOA < 19]
            
            # Select output filename
            ofile = f'{self.filename}_group{i+1}.ftr'

            # Save new catalogue
            df0.reset_index(drop=True, inplace=True)
            df0.to_feather(ofile)

            

        
        
#==============================================================#
#               PARSING COMMAND-LINE ARGUMENTS                 #
#==============================================================#

parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

out_group = parser.add_argument_group('I/O PARAMETERS')
out_group.add_argument('-p', '--plot',    action='store_true',      help='Flag to plot each step of the software')
out_group.add_argument('-v', '--verbose', type=int, metavar='INT',  help='Verbosity level [0, 1, 3] (Default: 1)')
out_group.add_argument('-o', '--outdir',  type=str, metavar='STR',  help='Output directory to save catalogue and plots')
out_group.add_argument('--project',       type=str, metavar='NAME', help='Name of project folder within $PLATO_WORKDIR')

pic_group = parser.add_argument_group('PIC QUERY (PLATO FIELDS)')
pic_group.add_argument('--pic',    type=str, nargs=3, metavar=('STARS', 'SAMPLE', 'FIELD'), action="append", help='Mandatory arguments (Check "--notes")')
pic_group.add_argument('--ncams',  type=int, metavar='INT',   help='N-CAM visibility [6, 12, 18, 24]')
pic_group.add_argument('--group',  type=int, metavar='INT',   help='Camera-group visibility [1, 2, 3, 4]')
pic_group.add_argument('--mag',    type=str, metavar='RANGE', help='P magnitude range of PIC targets (float or dash-range)')
pic_group.add_argument('--spec',   type=str, metavar='STR',   help='Spectral type [F, G, K]')
pic_group.add_argument('--lum',    type=str, metavar='STR',   help='Luminosity class [V, IV]')
pic_group.add_argument('--ntime',  type=int, metavar='INT',   help='Number of individual timeseries to be generated')
pic_group.add_argument('--incat',  type=str, nargs='*',       help='PIC target and contaminant input files [*.ftr, *.txt]')
pic_group.add_argument('--unique', type=str, metavar='STR',   help='Parse old target catalogue to select new unique stars [*ftr])')
pic_group.add_argument('--save',   action='store_true',       help='Flag to save stars into ascii PlatoSim-like catalogue')
pic_group.add_argument('--notes',  action='store_true',       help='Flag to show the notes of the available PIC catalogues')

bad_group = parser.add_argument_group('SIMBAD QUERY (SKY REGION)')
bad_group.add_argument('--simbad', type=str, metavar='NAME',  help='Simbad target name')

viz_group = parser.add_argument_group('VIZIER QUERY (PLATO FOV)')
viz_group.add_argument('--vizier', type=str,   metavar='FIELD', help='PLATO pointing field')
viz_group.add_argument('--magmin', type=float, metavar='MAG',   help='Min magnitude to query (Default: 0 mag)')
viz_group.add_argument('--magmax', type=float, metavar='MAG',   help='Max magnitude to query (Default: 15 mag)')
viz_group.add_argument('--yale_stars', action='store_true',     help='Flag to add the Yale bright stars catalogue')
viz_group.add_argument('--gaia_astro', action='store_true',     help='Flag to add parameters from the Astrophysical Gaia table')

que_group = parser.add_argument_group('GENERIC QUERY OPTIONS')
que_group.add_argument('--dmag', type=int, metavar='MAG',    help='Delta magnitude target-to-contaminant limit (Default: 5 mag)')
que_group.add_argument('--dist', type=int, metavar='ARCSEC', help='Radial distance target-to-contaminant limit (Default: 45 arcsec)')

args = parser.parse_args()


#--------------------------------------------------------------#
#                            WORKFLOW                          #
#--------------------------------------------------------------#

# Start time tracking
tic = datetime.datetime.now()

# Initialize instance of class
p = PicSim(args)

# Print notes only
if args.notes:
    p.printNotesPIC()
    exit()

# Query single star using Simbad
if args.simbad:
    p.initSimbad()
    p.querySimbad()

# Query a PLATO pointing field
elif args.vizier:
    p.initVizier()
    p.queryVizier()
    
# Query stars from the PIC     
elif args.pic:
    
    # Old picsim catalogue    
    if args.incat:
        p.loadOldPIC()
        p.plotTargetsPIC()
        p.plotContaminantsPIC()
        p.prologuePIC()
        
    # New picsim catalogue
    else:
        p.initPIC()
        p.loadPIC()
        p.queryTargetsPIC()
        p.plotTargetsPIC()
        p.queryContaminantsPIC()
        p.plotContaminantsPIC()
        p.prologuePIC()
    
# Finish with output
if (args.verbose is None) or (args.verbose > 0):
    toc = datetime.datetime.now()
    print(f'\nTotal execution time: {toc-tic} [hh:mm:ss]')
    print('')
