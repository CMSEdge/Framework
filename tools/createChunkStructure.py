#####################################################################
######                                                              #
###### 8=========D         OO                        8==========D   #  
###### OO                  ||                                ,88   #
###### ||                  ||                              ,88'    #  
###### ||O---\     ,adPPYb,|| ,adPPYb,d8  ,adPPYba,      ,88'      #
###### ||O---/    a8'    `Y||a8'    `Y88 a8P_____88    ,88'        #
###### ||         8b       ||8b       88 8PP'''''''  ,88'          #
###### \/         '8a,   ,d||'8a,   ,d88 '8b,   ,aa 88'            #
###### 8=========D `'8bbdP'\/  `'YbbdP'Y8  `'Ybbd8'' 8==========D   #
######                       aa,    ,88                             #
######                         'Y8bbdP'                             #
######                                                              #
#####################################################################

import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats
import math, sys, optparse, copy, subprocess




class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'



def processChunk(sampleName, eosdir, number):

    fileName = "heppyOutput_" + number + ".tgz"
    chunkName = sampleName + "_Chunk" + number
    subprocess.call("/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select cp /eos/cms" + eosdir + fileName + " .", shell=True)
    subprocess.call("tar xvfz " + fileName, shell=True)
    subprocess.call("mv Output " + chunkName, shell=True)
    
    rootFileName = "/eos/cms" + eosdir + "tree_" + number + ".root"
    dest = chunkName + "/treeProducerSusyMultilepton/tree.root"
    subprocess.call("/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select cp " + rootFileName + " " + dest, shell=True)


##Main body of the analysis
if __name__ == '__main__':

    print bcolors.HEADER 
    print '#######################################################################'
    print '             Starting the creation of chunk structures                 ' 
    print '#######################################################################' + bcolors.ENDC

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-d', '--directory', action='store', type=str, dest='eosdir', default='.', help='EOS directory containing the root and tgz files')
    parser.add_option('-n', '--name', action='store', type=str, dest='name', default='.', help='Name of the sample')
    (opts, args) = parser.parse_args()

    eosdir = "/eos/cms" + opts.eosdir
    sampleName = opts.name

    outputA = subprocess.Popen("/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select ls " + eosdir, stdout=subprocess.PIPE, shell=True)
    outputB = outputA.communicate()[0]
    output = outputB.split("\n")
  
    for line in output:
        if(line.find("tgz") != -1):
             pos1 = line.find('_')
             pos2 = line.find('.')
             number = line[pos1+1:pos2] 
             print bcolors.OKBLUE + "Processing chunk number " + number + bcolors.ENDC
             processChunk(sampleName, opts.eosdir, number)





 
