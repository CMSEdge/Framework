import ROOT, math, optparse, copy
import os


def LoadPoissonError():

    os.system('root -l -b -q PoissonError.C++')
    #ROOT.gROOT.LoadMacro('include/PoissonError.C+')
#LoadPoissonError()

if __name__ == "__main__":
    print '!!!!!!!!'

