from starPositionOnCCD                    import StarPositionOnCCD
from stellarVariability                   import StellarVariability
from StellarAberration.absoluteAberration import AbsoluteAberration
from fieldDistortion                      import FieldDistortion
from ThermoElasticDrift.tedFromFile       import TedFromFile
from ThermoElasticDrift.tedYawPitchRoll   import TedYawPitchRoll
from Jitter.jitterYawPitchRoll            import JitterYawPitchRoll
from Jitter.jitterFromFile                import JitterFromFile
from skyBackground                        import SkyBackGround
from Convolution.mappedGaussianPSF        import MappedGaussianPSF

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



L = []
succes = {True: "Success", False: "Failure"}




with suppress_stdout():
    test1 = StarPositionOnCCD()
    out1  = test1.run()

name = "Star position on CCD"
L.append("{:<9}  {:^36}:{}\n".format("Test1:", name, succes[out1]))
print("{:<9}  {:^36}:{}".format("Test1:", name, succes[out1]))


with suppress_stdout():
    test2   = StellarVariability()
    out2    = test2.run()

name = "Stellar Variability"
L.append("{:<9}  {:^36}:{}\n".format("Test2:", name, succes[out2]))
print("{:<9}  {:^36}:{}".format("Test2:", name, succes[out2]))


with suppress_stdout():
    test3   = AbsoluteAberration()
    out3    = test3.run()

name = "Absolute Aberration"
L.append("{:<9}  {:^36}:{}\n".format("Test3:", name, succes[out3]))
print("{:<9}  {:^36}:{}".format("Test3:", name, succes[out3]))


with suppress_stdout():
    test4   = FieldDistortion()
    out4    = test4.run()

name = "Field Distortion"
L.append("{:<9}  {:^36}:{}\n".format("Test4:", name, succes[out4]))
print("{:<9}  {:^36}:{}".format("Test4:", name, succes[out4]))


with suppress_stdout():
    test5_1 = TedYawPitchRoll()
    out5_1  = test5_1.run()

name = "Thermo-Elastic drift from noise"
L.append("{:<9}  {:^36}:{}\n".format("Test5.1:", name, succes[out5_1]))
print("{:<9}  {:^36}:{}".format("Test5.1:", name, succes[out5_1]))


with suppress_stdout():
    test5_2 = TedFromFile()
    out5_2  = test5_2.run()

name = "Thermo-Elastic drift from file"
L.append("{:<9}  {:^36}:{}\n".format("Test5.2:", name, succes[out5_2]))
print("{:<9}  {:^36}:{}".format("Test5.2:", name, succes[out5_2]))


with suppress_stdout():
    test6_1 = JitterYawPitchRoll()
    out6_1  = test6_1.run()

name = "Jitter from red noise"
L.append("{:<9}  {:^36}:{}\n".format("Test6.1:", name, succes[out6_1]))
print("{:<9}  {:^36}:{}".format("Test6.1:", name, succes[out6_1]))


with suppress_stdout():
    test6_2 = JitterFromFile()
    out6_2  = test6_2.run()

name = "Jitter from file"
L.append("{:<9}  {:^36}:{}\n".format("Test6.2:", name, succes[out6_2]))
print("{:<9}  {:^36}:{}".format("Test6.2:", name, succes[out6_2]))


with suppress_stdout():
    test7   = SkyBackGround()
    out7    = test7.run()

name = "Sky Background"
L.append("{:<9}  {:^36}:{}\n".format("Test7:", name, succes[out7]))
print("{:<9}  {:^36}:{}".format("Test7:", name, succes[out7]))


with suppress_stdout():
    test8_1 = MappedGaussianPSF()
    out8_1  = test8_1.run()

name = "Mapped Gaussian PSF"
L.append("{:<9}  {:^36}:{}\n".format("Test8.1:", name, succes[out8_1]))
print("{:<9}  {:^36}:{}".format("Test8.1:", name, succes[out8_1]))








f = open('../ioFiles/result.txt', 'w')
f.writelines(L)
f.close
