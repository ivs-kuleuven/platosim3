from starPositionOnCCD  import StarPositionOnCCD
from stellarVariability import StellarVariability
from StellarAberration.absoluteAberation import AbsoluteAberration
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
succes = {True: "Succes", False: "Failure"}
with suppress_stdout():
    test1 = StarPositionOnCCD()
    out1  = test1.run()

L.append("Test1: Star Position on CCD:  {0}\n".format(succes[out1]))

with suppress_stdout():
    test2 = StellarVariability()
    out2  = test2.run()

L.append("Test2: Stellar Variability: {}\n".format(succes[out2]))

with suppress_stdout():
    test3 = AbsoluteAberration()
    out3  = test3.run()

L.append("Test3: Absolute Aberration: {}\n".format(succes[out3]))


f = open('result.txt', 'w')
f.writelines(L)
f.close
