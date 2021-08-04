from starPositionOnCCD                             import StarPositionOnCCD
from stellarVariability                            import StellarVariability
from StellarAberration.absoluteAberration          import AbsoluteAberration
from fieldDistortion                               import FieldDistortion
from ThermoElasticDrift.tedFromFile                import TedFromFile
from ThermoElasticDrift.tedYawPitchRoll            import TedYawPitchRoll
from ThermoElasticDrift.tedOnCCDs                  import TedOnCCDs
from Jitter.jitterYawPitchRoll                     import JitterYawPitchRoll
from Jitter.jitterFromFile                         import JitterFromFile
from Jitter.jitterOnCCDs                           import JitterOnCCDs
from Jitter.jitterOnCameras                        import JitterOnCameras
from skyBackground                                 import SkyBackGround
from Convolution.mappedGaussianPSF                 import MappedGaussianPSF
from Convolution.analyticNonGaussian               import AnalyticNonGaussianPSF
from Convolution.analyticGaussian                  import AnalyticGaussianPSF
from prnu                                          import PRNU
from rebinning                                     import Rebinning
from ThroughputEfficiency.transmissionEfficiency   import TransmissionEfficiency
from ThroughputEfficiency.vignetting               import Vignetting
from ThroughputEfficiency.polarization             import Polarization
from ThroughputEfficiency.quantumEfficiency        import QuantumEfficiency
from ThroughputEfficiency.particulateContamination import ParticulateContamination
from ThroughputEfficiency.molecularContamination   import MolecularContamination
from DarkSignal.shotNoise                          import ShotNoise
from DarkSignal.darkSignalNonUniformity            import DarkSignalNonUniformity
from DarkSignal.tempVariationOfCCD                 import TempVariationOfCCD 
from brighterFatterEffect                          import BrighterFatterEffect
from cosmics                                       import Cosmics
from openShutterSmearing                           import OpenShutterSmearing
from ChargeTransferInefficiency.simpleCTI          import SimpleCTI
from photonNoise                                   import PhotonNoise
from readOutNoise                                  import ReadoutNoise
from fullWellSaturation                            import FullWellSaturation
from Quantisation.gain                             import Gain
from Quantisation.electronicOffset                 import ElectronicOffset
from Quantisation.flooring                         import Flooring
from Quantisation.digitalSaturation                import DigitalSaturation

from contextlib import contextmanager
import sys, os


@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout



testMessages = []

try:
    from termcolor import colored
    successMessage = colored('Success', 'green')
    failureMessage = colored('Failure', 'red', attrs=['reverse', 'blink'])
except ImportError:
    successMessage = 'Success'
    failureMessage = 'Failure'

success = {True: successMessage, False: failureMessage}




with suppress_stdout():
    test1 = StarPositionOnCCD()
    out1  = test1.run()

name = "Star position on CCD"
testMessages.append("{:<9}  {:^42}:{}".format("Test1:", name, success[out1]))
print(testMessages[-1])


with suppress_stdout():
    test2   = StellarVariability()
    out2    = test2.run()

name = "Stellar Variability"
testMessages.append("{:<9}  {:^42}:{}".format("Test2:", name, success[out2]))
print(testMessages[-1])


with suppress_stdout():
    test3   = AbsoluteAberration()
    out3    = test3.run()

name = "Absolute Aberration"
testMessages.append("{:<9}  {:^42}:{}".format("Test3:", name, success[out3]))
print(testMessages[-1])

#with suppress_stdout():
#    test4   = FieldDistortion()
#    out4    = test4.run()

#name = "Field Distortion"
#testMessages.append("{:<9}  {:^42}:{}".format("Test4:", name, success[out4]))
#print(testMessages[-1])


with suppress_stdout():
    test5_1 = TedYawPitchRoll()
    out5_1  = test5_1.run()

name = "Thermo-Elastic drift from noise"
testMessages.append("{:<9}  {:^42}:{}".format("Test5.1:", name, success[out5_1]))
print(testMessages[-1])


with suppress_stdout():
    test5_2 = TedFromFile()
    out5_2  = test5_2.run()

name = "Thermo-Elastic drift from file"
testMessages.append("{:<9}  {:^42}:{}".format("Test5.2:", name, success[out5_2]))
print(testMessages[-1])


with suppress_stdout():
    test5_3 = TedOnCCDs()
    out5_3 = test5_3.run()

name = "TED on different CCDs"
testMessages.append("{:<9}  {:^42}:{}".format("Test5.3:", name, success[out5_3]))
print(testMessages[-1])


with suppress_stdout():
    test6_1 = JitterYawPitchRoll()
    out6_1  = test6_1.run()

name = "Jitter from red noise"
testMessages.append("{:<9}  {:^42}:{}".format("Test6.1:", name, success[out6_1]))
print(testMessages[-1])


with suppress_stdout():
    test6_2 = JitterFromFile()
    out6_2  = test6_2.run()

name = "Jitter from file"
testMessages.append("{:<9}  {:^42}:{}".format("Test6.2:", name, success[out6_2]))
print(testMessages[-1])


with suppress_stdout():
    test6_3 = JitterOnCCDs()
    out6_3  = test6_3.run()

name = "Jitter on different CCDs"
testMessages.append("{:<9}  {:^42}:{}".format("Test6.3:", name, success[out6_3]))
print(testMessages[-1])


with suppress_stdout():
    test6_4 = JitterOnCameras()
    out6_4 = test6_4.run()

name = "Jitter on different cameras"
testMessages.append("{:<9}  {:^42}:{}".format("Test6.4:", name, success[out6_4]))
print(testMessages[-1])


with suppress_stdout():
    test7   = SkyBackGround()
    out7    = test7.run()

name = "Sky Background"
testMessages.append("{:<9}  {:^42}:{}".format("Test7:", name, success[out7]))
print(testMessages[-1])


with suppress_stdout():
    test8_1 = MappedGaussianPSF()
    out8_1  = test8_1.run()

name = "Mapped Gaussian PSF and Charge Diffusion"
testMessages.append("{:<9}  {:^42}:{}".format("Test8.1:", name, success[out8_1]))
print(testMessages[-1])


with suppress_stdout():
    test8_2 = AnalyticNonGaussianPSF()
    out8_2  = test8_2.run()

name = "Analytic non Gaussian PSF"
testMessages.append("{:<9}  {:^42}:{}".format("Test8.2:", name, success[out8_2]))
print(testMessages[-1])


with suppress_stdout():
    test8_3 = AnalyticGaussianPSF()
    out8_3  = test8_3.run()

name = "Analytic Gaussian PSF"
testMessages.append("{:<9}  {:^42}:{}".format("Test8.3:", name, success[out8_3]))
print(testMessages[-1])


with suppress_stdout():
    test9 = PRNU()
    out9  = test9.run()

name = "Pixel-Responsivity Non Uniformity"
testMessages.append("{:<9}  {:^42}:{}".format("Test9:", name, success[out9]))
print(testMessages[-1])


with suppress_stdout():
    test10 = Rebinning()
    out10  = test10.run()

name = "Rebinning Subpixel"
testMessages.append("{:<9}  {:^42}:{}".format("Test10:", name, success[out10]))
print(testMessages[-1])


with suppress_stdout():
    test11_1 = TransmissionEfficiency()
    out11_1  = test11_1.run()

name = "Transmission Efficiency"
testMessages.append("{:<9}  {:^42}:{}".format("Test11.1:", name, success[out11_1]))
print(testMessages[-1])


with suppress_stdout():
    test11_2 = Vignetting()
    out11_2  = test11_2.run()

name = "Vignetting"
testMessages.append("{:<9}  {:^42}:{}".format("Test11.2:", name, success[out11_2]))
print(testMessages[-1])


with suppress_stdout():
    test11_3 = Polarization()
    out11_3  = test11_3.run()

name = "Polarization"
testMessages.append("{:<9}  {:^42}:{}".format("Test11.3:", name, success[out11_3]))
print(testMessages[-1])


with suppress_stdout():
    test11_4 = QuantumEfficiency()
    out11_4  = test11_4.run()

name = "Quantum Efficiency"
testMessages.append("{:<9}  {:^42}:{}".format("Test11.4:", name, success[out11_4]))
print(testMessages[-1])


with suppress_stdout():
    test11_5 = ParticulateContamination()
    out11_5  = test11_5.run()

name = "Particulate Contamination"
testMessages.append("{:<9}  {:^42}:{}".format("Test11.5:", name, success[out11_5]))
print(testMessages[-1])


with suppress_stdout():
    test11_6 = MolecularContamination()
    out11_6  = test11_6.run()

name = "Molecular Contamination"
testMessages.append("{:<9}  {:^42}:{}".format("Test11.6:", name, success[out11_6]))
print(testMessages[-1])


with suppress_stdout():
    test12_1 = ShotNoise()
    out12_1  = test12_1.run()

name = "Shot Noise"
testMessages.append("{:<9}  {:^42}:{}".format("Test12.1:", name, success[out12_1]))
print(testMessages[-1])


with suppress_stdout():
    test12_2 = DarkSignalNonUniformity()
    out12_2  = test12_2.run()

name = "Dark Signal Non Uniformity"
testMessages.append("{:<9}  {:^42}:{}".format("Test12.2:", name, success[out12_2]))
print(testMessages[-1])


with suppress_stdout():
    test12_3 = TempVariationOfCCD()
    out12_3  = test12_3.run()

name = "Temperature Variation of CCD"
testMessages.append("{:<9}  {:^42}:{}".format("Test12.3:", name, success[out12_3]))
print(testMessages[-1])


with suppress_stdout():
    test13 = BrighterFatterEffect()
    out13  = test13.run()

name = "Brighter-Fatter Effect"
testMessages.append("{:<9}  {:^42}:{}".format("Test13:", name, success[out13]))
print(testMessages[-1])

with suppress_stdout():
    test14 = Cosmics()
    out14  = test14.run()

name = "Cosmics"
testMessages.append("{:<9}  {:^42}:{}".format("Test14:", name, success[out14]))
print(testMessages[-1])

with suppress_stdout():
    test15 = OpenShutterSmearing()
    out15  = test15.run()

name = "Open-Shutter Smearing"
testMessages.append("{:<9}  {:^42}:{}".format("Test15:", name, success[out15]))
print(testMessages[-1])

with suppress_stdout():
    test16_1 = SimpleCTI()
    out16_1  = test16_1.run()

name = "Simple CTI model"
testMessages.append("{:<9}  {:^42}:{}".format("Test16.1:", name, success[out16_1]))
print(testMessages[-1])

with suppress_stdout():
    test17 = PhotonNoise()
    out17  = test17.run()

name = "Photon Noise"
testMessages.append("{:<9}  {:^42}:{}".format("Test17:", name, success[out17]))
print(testMessages[-1])

with suppress_stdout():
    test18 = ReadoutNoise()
    out18  = test18.run()

name = "Readout Noise"
testMessages.append("{:<9}  {:^42}:{}".format("Test18:", name, success[out18]))
print(testMessages[-1])

with suppress_stdout():
    test19 = FullWellSaturation()
    out19  = test19.run()

name = "Full-Well Saturation"
testMessages.append("{:<9}  {:^42}:{}".format("Test19:", name, success[out19]))
print(testMessages[-1])

with suppress_stdout():
    test20_1 = Gain()
    out20_1  = test20_1.run()

name = "Gain"
testMessages.append("{:<9}  {:^42}:{}".format("Test20.1:", name, success[out20_1]))
print(testMessages[-1])

with suppress_stdout():
    test20_2 = ElectronicOffset()
    out20_2  = test20_2.run()

name = "Electronic Offset"
testMessages.append("{:<9}  {:^42}:{}".format("Test20.2:", name, success[out20_2]))
print(testMessages[-1])

with suppress_stdout():
    test20_3 = Flooring()
    out20_3  = test20_3.run()

name = "Flooring"
testMessages.append("{:<9}  {:^42}:{}".format("Test20.3:", name, success[out20_3]))
print(testMessages[-1])

with suppress_stdout():
    test20_4 = DigitalSaturation()
    out20_4  = test20_4.run()

name = "Digital Saturation"
testMessages.append("{:<9}  {:^42}:{}".format("Test20.4:", name, success[out20_4]))
print(testMessages[-1])










with open(os.environ["PLATO_PROJECT_HOME"] + '/tests/validationTests/ioFiles/result.txt', 'w') as outputFile:
    for line in testMessages:
        outputFile.write(line+'\n')
