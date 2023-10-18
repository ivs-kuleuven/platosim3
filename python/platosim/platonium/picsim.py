#!/usr/bin/env python3

"""
This script is an integrated part of PlatoSim's toolkit PLATOnium.
 takes the PIC input catalog and fetch the 
correct entries for target stars and their contaminants. Notice that
the catalog will only be saved if an ouput directory and prefix (-o)
are given as input. This script can also be used to re-plot a old 
output ascii catalog produced by this software. Note that the PIC
 contaminant catalog is complete up to 20 mag.

Notes on PLATO fields:
  Number of stars in PIC200 LOPS2: (t:179,564, c:105,865,120)
  Number of stars in PIC200 LOPN1: (t:175,597, c:112,917,958)
  Number of stars in PIC110 SPF  : (t:163,772, c: )
  Number of stars in PIC110 NPF  : (t:156,971, c: )

Note on Sample flag:
  Number of targets in P1:  15,094 (SPF:   6,817, NPF:   6,892)
  Number of targets in P2:   1,385 (SPF:     717, NPF:     668)
  Number of targets in P4:  33,032 (SPF:  16,866, NPF:  16,166)
  Number of targets in P5: 272,617 (SPF: 132,571, NPF: 140,046)
  Note that P2 is a bright sub-sample of P1.
  Note that P5 is on-board processed lightcures only.

Notes on Calibration methods:
  PIC110 uses the de-reddened Gaia V magnitude obtained from Gaia
  photometry using a calibration technique being valid for:
  0.5 < Bp-Rp < 5.17
  PIC200

Usage examples:
  >> picsim 100 P1 SPF -p
  >> picsim 1000 P5 NPF --ncams 24 --project <project_name>
  >> picsim all all SPF --ncams 6 --group 1 --mag 9.5-12.2 -o </path/to/outdir>
"""

# Python standard
import os
import shutil
import random
import warnings
import argparse
import datetime

# PlatoSim standard
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from prettytable import PrettyTable
from tqdm import tqdm

# PlatoSim functions
import platosim.plot      as pt
import platosim.utilities as ut
import platosim.starquery as sq
from platosim.utilities  import errorcode
from platosim.simulation import Simulation


#==============================================================#
#                         BEGIN CLASS                          #
#==============================================================#


class PicSim(object):

    """Class to generate customised (PIC) catalogues.
    """
    
    def __init__(self, args):

        # Arguments for PIC query
        self.plot          = args.plot
        self.saveAscii     = args.save
        self.inputFiles    = args.incat
        self.oldCatalogue  = args.unique
        self.numTimeseries = args.ntime

        # Arguments for GaiaDR3 star query
        self.simbad = args.simbad

        
        # MANDATOEY PARAMETERS
        
        # Check if sample is correct
        self.sample = args.sample
        if self.sample == 'all':
            self.samplePIC = None
        elif self.sample in ['P1', 'P2', 'P4', 'P5']:
            self.samplePIC = self.sample
        elif self.simbad:
            pass
        else:
            errorcode('error', 'Not a valid PIC sample! Use {P1, P2, P4, P5}')

        # Check the PLATO field
        self.field  = args.field
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
        if args.stars[0].isdigit():
            self.numTargets = int(args.stars)
            if self.numTargets > self.numPIC:
                errorcode('error', 'More stars selected than available in the PIC! ' +
                          'See -h for PIC notes')
        elif args.stars == 'all':
            if self.numTimeseries is not None:
                self.numTargets = int(self.numTimeseries/6.) + 1
            else:
                self.numTargets = self.numPIC

            
        # OBS PARAMETERS            
            
        # Select magnitude range or select all stars
        if args.mag is None:
            self.magRange = [0, 21]
        else:
            self.magRange = ut.convertMagnitudeRange(magRange)

        # Visibility by N-cams
        if args.ncams in [None, 6, 12, 18, 24]:
            self.ncamVisible = args.ncams
        else:
            errorcode('Not valid N-CAM visibility! Use {6, 12, 18, 24}')

        # Check cam-group entry
        if args.group in [None, 1, 2, 3, 4]:
            self.camGroup = args.group
        else:
            errorcode('error', 'Not valid group number! Use {1, 2, 3, 4}')

        # Contaminant magnitude limit
        if args.dmag is None:
            self.dmagConLimit = 5
        elif args.dmag > 21:
            errorcode('warning', 'PIC (Gaia) is only complete to G < 21 mag!')

        # Contaminant distance limit (default 2 pixel)
        if args.dist is None:
            self.disConLimit = 45
        elif args.dist in [30, 45, 60]:
            self.disConLimit = args.dist
        else:
            errorcode('error', 'Not a valid contaminant-to-target distance! ' +
                      'Use {30, 45, 60} arcsec')

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

            
        # I/O PARAMETERS

        # Activate plot if `-t` is parsed
        if args.targ:
            self.preview = True
            self.plot    = True
        else:
            self.preview = False
        
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
        
        # Input directory        
        self.inputDir = Path(os.getenv("PLATO_PROJECT_HOME")) / 'inputfiles/data_picsim'    
        if not self.inputDir.is_dir():
            errorcode('message', 'Inuaguration: Welcome to PIC of Destiny!')
        
        # Prefix name of output files
        self.outputPrefix = 'starcat'
        if args.sample: self.outputPrefix += '_' + args.sample
        if args.field:  self.outputPrefix += '_' + args.field
        if args.group:  self.outputPrefix += '_Group' + str(args.group)
        if args.ncams:  self.outputPrefix += '_Ncam'  + str(args.ncams)
        if self.inputFiles is not None: self.outputPrefix += '_NewCat'
        self.outputPrefixTar = self.outputPrefix + '_targets.ftr'
        self.outputPrefixCon = self.outputPrefix + '_contaminants.ftr'

        # Name space for output files
        # Default is to use the latest PIC catalog saved
        # Else several old catalogs can be parsed for replotting
        # E.g. use: --intar starcat**targets.txt --incon starcat**contaminants.txt
        if (args.outdir is not None) or (args.project is not None):

            # Resolve absolute output directory
            # NOTE "--outdir" overwrites "--project"
            if args.outdir is not None:
                self.outputDir = Path(args.outdir).resolve()
            elif args.project is not None:
                self.outputDir = Path(os.getenv('PLATO_WORKDIR')) / args.project / 'input'

            # Create directory if is doesn't exist
            self.outputDir.mkdir(parents=True, exist_ok=True)
                
            # Save either a PlatoSim (ascii) or a PLATOnium (feather) calogue
            if self.saveAscii:
                self.outputFileCat = self.outputDir / self.outputPrefix/ '.txt'
            else:
                self.outputFileTar = self.outputDir / self.outputPrefixTar
                self.outputFileCon = self.outputDir / self.outputPrefixCon
        else:
            self.outputDir = args.outdir

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

            
        # OBS PARAMETERS

        # Fetch PLATO pointing field [deg]
        self.alpha, self.delta, self.kappa = ut.getPointingField(self.field, unit='deg')

            


            
    def loadPIC(self):

        """Fetch PIC targets from feather files.
        """
        
        # Files names
        inputFileTar = self.inputDir / f'{self.pic}_{self.field}_targets.ftr'
        inputFileCon = self.inputDir / f'{self.pic}_{self.field}_contaminants.ftr'

        # Fetch PIC catalogue from FTP server
        if not inputFileTar.is_file() or not inputFileCon.is_file():
            print(f'Downloading {self.pic} catalogue..')
            ut.downloadFromFTP(inputFileTar.name, self.inputDir, 'plato')
            ut.downloadFromFTP(inputFileCon.name, self.inputDir, 'plato')

        # Load catalogues
        if self.verbose > 0:
            print('Loading stellar catalogues..')
        df = pd.read_feather(inputFileTar)
        dc = pd.read_feather(inputFileCon)

        # Remove any star with NaN value
        self.df0 = df.dropna()
        self.dc0 = dc.dropna()

        # Removing dublicate stars (if any)
        self.df0 = self.df0.drop_duplicates(subset=['PIC'])

        

    

        
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
                self.df0 = self.df.append(pd.read_feather(inputFileTar[i]))
                self.dc0 = self.dc.append(pd.read_feather(inputFileCon[i]))





    
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
        if self.ncamVisible is not None:
            df = df[df['ncams'] == int(self.ncamVisible)]

        # Allow to select only one cam-group
        if self.camGroup:
            sim = Simulation('picsim')
            dex = sim.getStarsWithinCameraGroup(df.ra.to_numpy(), df.dec.to_numpy(),
                                                self.alpha, self.delta, self.kappa,
                                                self.camGroup)

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
            errorcode('warning', f'Only {len(df)} stars are available after cuts!')
            self.numTargets = len(df)
            

        # RANDOM CATALOGUE SELECTION

        # Only after the above cuts can we select randomly a sub-sample (if applicable)
        # Use random.choices(data, weights=some weight, k=numTargets) if a weigthing is needed
        # NOTE args.stars being a string is needed to select "all" in one of the samples too
        # NOTE we need to shuffle the sample if a specific number of timeseries are needed
        if args.stars != 'all' or self.numTimeseries is not None:
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
            print(f'Fetching contaminants within {self.disConLimit} arcsec from each target')
            
        df = self.df
        dc = self.dc0

        # Create empty data-frame to append to 
        dc1 = pd.DataFrame(columns=['PIC', 'ra', 'dec', 'mag', 'dis'])
        numConPerTar = np.zeros(self.numTargets)

        for starNo in tqdm(range(self.numTargets), bar_format=ut.tqdmBar()):

            # Find corresponding contaminants to our target
            dcc = dc[dc['PIC'] == df['PIC'].iloc[starNo]]

            # We will only add stars within a user-defined distance of our target
            dcc = dcc[dcc['dis'] < self.disConLimit]

            # Add only contaminants brigther than threshold
            dcc = dcc[dcc['mag'] < df['mag'].iloc[starNo] + self.dmagConLimit]

            # Append contaminats to contaminant list
            dc1 = dc1.append(dcc)

            # Store number of contaminats per target for statistics below
            numConPerTar[starNo] = len(dcc)

        # Add the number of contaminats to the target catalog
        self.df['ncon'] = numConPerTar
        self.dc = dc1

    

        
            
    def plotsTargetsPIC(self):

        """Plot function for PIC targets.
        """
        
        df = self.df
        
        # Plot targets in aitoff Galactic sky projection
        self.fig0, ax = pt.drawStarsInSkyAitoff(df.ra, df.dec, df.mag)
        if self.plot: plt.show()

        # Plot Zoom-in on FOV
        if df.shape[0] > 200: mag = None
        else: mag = df.mag.to_numpy()
        self.fig1, ax = pt.plotPlatoFOV(self.field, df.ra, df.dec, magStars=mag,
                                        showGroups=False, title=self.title)
        if self.plot: plt.show()

        # Plot sample distribution in Teff vs. Radius
        _, ds, sg, dK, dG, dF, sgK, sgG, sgF = self.getStellarClass(self.df0)
        self.fig2, ax = pt.plotTeffvsRadius(ds, dK, dG, dF, sg, sgK, sgG, sgF, df,
                                            self.title)
        if self.plot: plt.show()

        # Stop here when flag -t
        if self.preview: exit()


        

        
    def plotsContaminantsPIC(self):

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




        
    def prologuePIC(self, tictoc):

        if self.verbose > 0:
            errorcode('module', '\nPrologue')
        
        if self.outputDir is not None:

            # Copy the YAML file to the project if it doesn't exist
            if self.verbose > 0:
                print(f"Copying YAML configuration file")
            ut.copyInputYAML(self.field, self.outputDir)

            # Save textfile with information about PIC stars
            log = ('PIC Catalogue include\n' +
                   f'Catalogue            : {self.pic}\n' +
                   f'PLATO field          : {self.field}\n' +
                   f'PLATO sample         : {self.sample}\n' +
                   f'PLATO camera-group   : {self.camGroup}\n' +
                   f'Magnitude range      : {self.magRange[0]:.2f}-{self.magRange[1]:.2f} \n' +
                   f'Luminosity class     : {self.lumClass}\n' +
                   f'Spectral type        : {self.specType}\n' +
                   f'Contaminant distance : {self.disConLimit} arcsec\n' +
                   f'Contaminant magnitude: {self.dmagConLimit} mag\n' +
                   '\n' + self.t.get_string())
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
                df0        = pd.concat([df.ra,  dfc.ra])
                df0['dec'] = pd.concat([df.dec, dfc.dec])
                df0['mag'] = pd.concat([df.mag, dfc.mag])
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
                outputFileVarSim = f'{self.outputDir}/cluster_{self.outputPrefix}.txt'
                if self.verbose > 0:
                    print(f'Saving file {outputFileVarSim}')
                headerVar = 'ID,PIC,M,R,Teff'
                df_hpc         = df.ID
                df_hpc['PIC']  = df.PIC
                df_hpc['Teff'] = df.Teff
                df_hpc['R']    = df.R
                df_hpc['M']    = df.M
                df_hpc.to_csv(outputFileVarSim, sep=',',
                              float_format=['%s', '%i', '%0.3f', '%0.3f', '%i'])
                



    #--------------------------------------------------------------#
    #                   SINGLE GAIA STAR QUERY                     #
    #--------------------------------------------------------------#            

    def queryStar(self):

        """Function to query a Gaia DR3 star and its contaminants.
        """

        # Query Gaia DR3 star
        df = sq.simbadQuery(self.simbad, radius=self.disConLimit)
        if self.verbose > 0:
            print(f'\nCatalogue for {self.simbad}:')
            print(df)
            
        # Save output catalogue
        if self.outputDir is not None:
            da = df.loc[:, ['ra', 'dec', 'mag']]
            da.to_csv(self.outputDir / f"starcat_{self.simbad}.txt", sep=' ', header=False)
            df.to_feather(self.outputDir / f"starcat_{self.simbad}.ftr")





    #--------------------------------------------------------------#
    #                      PIC-VARSIM CATALOGS                     #
    #--------------------------------------------------------------#


    def loadNumpyTargetsPIC110(inputFileTar):

        """Function to load target catalogue. 

        PICidDR1   : PIC-ID-DR1 from Gaia DR2
        ra         : ICRS RA
        decl       : ICRS Dec
        gaiaV      : De-reddened V mag from Gaia colour photometry
        sampleFlag : Bitmaskdefining PIC samples
        teff       : Stellar effective temperature [K]
        radius     : Stellar radius [R_sun]
        mass       : Stellar mass [M_sun]
        nCameraObs : EOL number of cameras seeing the star
        """
        df = pd.DataFrame()
        pic_tar = np.load(inputFileTar)
        df['PIC']    = pic_tar[:,0].astype(float).astype(int)
        df['ra']     = pic_tar[:,1].astype(np.float64)
        df['dec']    = pic_tar[:,2].astype(np.float64)
        df['mag']    = pic_tar[:,3].astype(np.float64)
        df['sample'] = pic_tar[:,4].astype(float).astype(int)
        df['Teff']   = pic_tar[:,5].astype(np.float64)
        df['R']      = pic_tar[:,6].astype(np.float64)
        df['M']      = pic_tar[:,7].astype(np.float64)
        df['ncams']  = pic_tar[:,8].astype(float).astype(int)
        df['field']  = pic_tar[:,9].astype(str)
        df = df.iloc[1:]

        return df.reset_index()





    def loadNumpyContaminantsPIC110(inputFileCon):

        """Function to load target catalogue. 
        """

        df = pd.DataFrame()
        pic_con = np.load(inputFileCon)
        df['PIC'] = pic_con[:,0].astype(float).astype(int)
        df['ra']  = pic_con[:,1].astype(float)
        df['dec'] = pic_con[:,2].astype(float)
        df['mag'] = pic_con[:,4].astype(float)
        df['dis'] = pic_con[:,3].astype(float)

        return df





    def convertPIC110(path):

        """Function to create PIC110 target file.

        This function loads the original PIC110 ascii file
        containing the PIC targets, selects the right columns
        and saves them to a binary feather file.
        """

        # PIC TARGETS

        # Load ascii catalogue
        data = np.genfromtxt(f'{path}/filename.csv', delimiter=',',
                            usecols=[0, 3, 5, 53, 55, 56, 58, 60, 71])
        df = pd.DataFrame()
        df['PIC']    = data[:,0].astype(float).astype(int)
        df['ra']     = data[:,1].astype(np.float64)
        df['dec']    = data[:,2].astype(np.float64)
        df['mag']    = data[:,3].astype(np.float64)
        df['Teff']   = data[:,5].astype(np.float64)
        df['R']      = data[:,6].astype(np.float64)
        df['M']      = data[:,7].astype(np.float64)
        df['ncams']  = data[:,8].astype(float).astype(int)
        df['sample'] = data[:,4].astype(float).astype(int)

        # String field needs to be loaded seperately: PLATO field: N=North, S=South
        df['field'] = np.loadtxt(inputFileTar.with_suffix('.csv'), delimiter=',',
                                 usecols=[68], dtype=str)

        # Drop nan rows
        df = df.dropna()

        # Store columns
        col_sample = df['sample']
        df = df.drop(columns=['sample'])
        df['sample'] = col_sample

        # Select PIC sample
        # NOTE 2 is stated in documentation but 3 is correct..
        df['sample'] = df['sample'].replace([1, 3, 4, 8], ['P1', 'P2', 'P4', 'P5'])

        # Change camera numbers
        df['ncams'] = df['ncams'].replace([5, 11, 16, 17, 22], [6, 12, 18, 18, 24])

        # Convert V Jonhson-Cousin to P passband
        df['mag'] = ut.passbandConversionV2P(df.mag, df.Teff)

        # Select catalogues
        ds = df[df.field == 'S']
        dn = df[df.field == 'N']

        # Drop field before saving
        ds = ds.drop(columns=['field'])
        dn = dn.drop(columns=['field'])

        # Reset indices
        ds = ds.reset_index(drop=True)
        dn = ds.reset_index(drop=True)

        # Save to feather files
        ds.to_feather('PIC110_SPF_targets.ftr')
        dn.to_feather('PIC110_NPF_targets.ftr')

        
        # PIC CONTAMINATS

        # Load data
        # NOTE The PIC is the target to which the contaminants refers to
        data = np.loadtxt(inputFileCon.with_suffix('.csv'), delimiter=',',
                          skiprows=1, usecols=[2, 6, 8, 5, 19])
        dc =pd.DataFrame()
        dc['PIC'] = data[:,0].astype(float).astype(int)
        dc['ra']  = data[:,1].astype(np.float64)
        dc['dec'] = data[:,2].astype(np.float64)
        dc['mag'] = data[:,4].astype(np.float64)
        dc['dis'] = data[:,3].astype(float).astype(int)

        # Convert Vmag to Pmag using host star Teff
        # NOTE assumption needed for PlatoSim!
        for i in tqdm(range(len(dc)), bar_format=ut.tqdmBar()):
            df_i = df[df.PIC == dc.PIC.iloc[i]]
            dc.mag.iloc[i] = ut.passbandConversionV2P(dc.mag.iloc[0], df_i.Teff)
        
        # Sort after pointing
        ds = df[(ds.dec < 0)]
        dn = df[(df.dec > 0)]

        # Save to feather files
        ds.to_feather('PIC110_SPF_contaminants.ftr')
        dn.to_feather('PIC110_NPF_contaminants.ftr')





    def createPIC200(path, field):

        import pandas as pd
        from astropy.io.votable import parse
        from platosim.utilities import votable2pandas
        
        # TARGETS
            
        # Load targets
        tfile = f'p{field}PICtarget2.0.0.1-t.vot'
        votable = parse(f'{path}/{tfile}')
        df0 = votable2pandas(votable)

        # Write relevant columns to df
        df = pd.DataFrame()
        df['PIC']   = df0.PICid
        df['ra']    = df0.RAdeg
        df['dec']   = df0.DEdeg
        df['mag']   = df0.PlatoMagNCAM
        df['BPmag'] = df0.PlatoMagFCAMb
        df['RPmag'] = df0.PlatoMagFCAMr
        df['Teff']  = df0.Teff
        df['R']     = df0.Radius
        df['M']     = df0.Mass
        df['ncams'] = df0.BOLnCameraObs

        # Remove NaNs
        df = df.dropna()

        # Remove stellar duplicates (if any)
        df = df.drop_duplicates(subset=['PIC'])
        
        # Save to feather files
        df = df.reset_index(drop=True)
        df.to_feather(f'PIC200_{field}_targets.ftr')

        # CONTAMINANTS

        # File is huge hence read only one column at the time
        cfile = f'p{field}PICcontaminant2001t.csv'

        # Create data frame
        dc = pd.DataFrame()
        dc['PIC']   = pd.read_csv(f'{path}/{cfile}', usecols=['PICcontaminantId'])
        dc['ra']    = pd.read_csv(f'{path}/{cfile}', usecols=['RAdeg'])
        dc['dec']   = pd.read_csv(f'{path}/{cfile}', usecols=['DEdeg'])
        dc['mag']   = pd.read_csv(f'{path}/{cfile}', usecols=['Gmag'])
        dc['BPmag'] = pd.read_csv(f'{path}/{cfile}', usecols=['BPmag'])
        dc['RPmag'] = pd.read_csv(f'{path}/{cfile}', usecols=['RPmag'])

        # Remove NaNs
        dc = dc.dropna()

        # Use colour proxy for bandpass conversion
        #Teff = 
        #df['mag'] = ut.passbandConversionV2P(df.mag, Teff)
        
        # Save to feather files
        dc = dc.reset_index(drop=True)
        dc.to_feather(f'PIC200_{field}_contaminants.ftr')




        
    def cameraObservability(self):
        
        # CAMERA OBSERVABILITY

        # Since the PIC are bigger than actual FOV (18.89 deg radius) we need to
        # Loop over each camera-group and avoid stars that are close to the FOV edge
        # NOTE the following also checks if the subfield can be placed on a CCD
        camGroup = False
        if camGroup:

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
            if ncamVisible:
                dex = ncams == int(ncamVisible)
                PIC   = PIC[dex].astype(int)
                ra    = ra[dex]
                dec   = dec[dex]
                mag   = mag[dex]
                Teff  = Teff[dex]
                R     = R[dex]
                M     = M[dex]
                ncams = ncams[dex]

            # Allow to select only one cam-group
            if camGroup:
                dex = starsWithinCamGroup(camGroup, raPF, decPF, ra, dec)
                PIC   = PIC[dex].astype(int)
                ra    = ra[dex]
                dec   = dec[dex]
                mag   = mag[dex]
                Teff  = Teff[dex]
                R     = R[dex]
                M     = M[dex]
                ncams = ncams[dex]




                
    def massLuminosityRelation(R, Teff):

        """Calculate mass from proxy luminosity.

        Correct for missing masses using mass-luminosity relation.
        Method valid for (0.43 < M/Msun < 2) and reference from:
        https://en.wikipedia.org/wiki/Mass%E2%80%93luminosity_relation
        """
        return R**(1/2) * Teff[i]/5778.


                
        
#==============================================================#
#               PARSING COMMAND-LINE ARGUMENTS                 #
#==============================================================#

parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=errorcode('software', '\nPIC of Destiny'))

pic_group = parser.add_argument_group('PIC INFO (required)')
pic_group.add_argument('stars',  type=str, help='Number of targets saved in catalog (Select "all" for all stars)')
pic_group.add_argument('sample', type=str, help='PIC samples [P1, P2, P4, P5, all]  (Select "all" for all samples)')
pic_group.add_argument('field',  type=str, help='PLATO fields [LOPS2, LOPN1, SPF, NPF]')

out_group = parser.add_argument_group('I/O PARAMETERS')
out_group.add_argument('-p', '--plot', action='store_true', help='Flag to plot target and contaminant catalog creation')
out_group.add_argument('-t', '--targ', action='store_true', help='Flag to only show catalogue illustration (mainly testing)')
parser.add_argument('-v', '--verbose', metavar='NUM',  type=int, help='Verbosity level [0, 1, 2] (Default: 1)')
out_group.add_argument('-s', '--save', action='store_true', help='Flag to save target and contaminant catalog into a PlatoSim like starcatalog format')
out_group.add_argument('-o', '--outdir', type=str, metavar='STR', help='Output directory to save catalogue and plots')
out_group.add_argument('--project',     type=str, metavar='NAME', help='Name project folder within $PLATO_WORKDIR')
out_group.add_argument('-i', '--incat',  type=str, nargs='*',     help='PIC input target- and contaminant files (Format: csv, txt, ftr)')
out_group.add_argument('-u', '--unique', type=str, metavar='STR', help='Parse old picsim target catalogue to select new unique stars (Format: ftr)')

que_group = parser.add_argument_group('QUERY OPTIONS (optional)')
que_group.add_argument('--simbad', type=str, metavar='NAME',  help='Simbad target for query using Gaia DR2')
que_group.add_argument('--dmag',  type=int, metavar='MAG',    help='Delta magnitude inclusion threhold of target-to-contaminatn (Default: 5 mag)')
que_group.add_argument('--dist',  type=int, metavar='ARCSEC', help='Radial distance to fecth contaminants [30, 45, 60] (Default: 45 arcsec)')

obs_group = parser.add_argument_group('OBSERVABLES (optional)')
obs_group.add_argument('--ncams', metavar='NUM',   type=int, help='N-CAM visibility [6, 12, 18, 24]')
obs_group.add_argument('--group', metavar='NUM',   type=int, help='Camera-group visibility [1, 2, 3, 4]')
obs_group.add_argument('--mag',   metavar='RANGE', type=str, help='P magnitude range of PIC targets (float or dash-range)')
obs_group.add_argument('--spec',  metavar='TYPE',  type=str, help='Spectral type [F, G, K]')
obs_group.add_argument('--lum',   metavar='CLASS', type=str, help='Luminosity class [V, IV]')
obs_group.add_argument('--ntime', metavar='NUM',   type=int, help='Number of individual timeseries to be generated')

args = parser.parse_args()

#--------------------------------------------------------------#
#                            WORKFLOW                          #
#--------------------------------------------------------------#

# Start time tracking
tic = datetime.datetime.now()

# Initialize instance of class
p = PicSim(args)

# Query single star using Simbad
if args.simbad:
    p.queryStar()

# Old picsim catalogue
elif args.incat:
    p.loadOldPIC()
    p.plotsTargetsPIC()
    p.plotsContaminantsPIC()
    p.prologuePIC()
    
# New picsim catalogue
else:
    p.loadPIC()
    p.queryTargetsPIC()
    p.plotsTargetsPIC()
    p.queryContaminantsPIC()
    p.plotsContaminantsPIC()
    p.prologuePIC()
    
# Finish with output
toc = datetime.datetime.now()
print(f'\nTotal execution time: {toc-tic} [hh:mm:ss]')
print('')

