import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TH1F, TGraph
import math, sys, optparse, copy, re, array, itertools, os


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Tables     as Tables




if __name__ == "__main__":
   

    gROOT.SetBatch(True)
    gROOT.ProcessLine('.L include/tdrstyle.C')
    r.setTDRStyle() 
    r.gStyle.SetPaintTextFormat("4.2f")
    
    
    global lumi
    lumi= 35
    mcDatasets = ['TT_pow_ext34']
    dyDatasets = ['DYJetsToLL_M10to50','DYJetsToLL_M50_HT100to200_ext','DYJetsToLL_M50_HT200to400_ext',
                  'DYJetsToLL_M50_HT400to600_ext','DYJetsToLL_M50_HT600toInf']

    treeTT = Sample.Tree(helper.selectSamples('samples.dat', mcDatasets, 'MC'), 'MC', 0, isScan = 0)
    treeDY = Sample.Tree(helper.selectSamples('samples.dat', dyDatasets, 'MC'), 'MC', 0, isScan = 0)

    cuts = CutManager.CutManager()
    mllCut = {'low': 'lepsMll_Edge < 81',
              'onz': 'lepsMll_Edge > 81 && lepsMll_Edge < 101',
              'high1': 'lepsMll_Edge > 101 && lepsMll_Edge < 200',
              'high2': 'lepsMll_Edge > 200 && lepsMll_Edge < 300',
              'high3': 'lepsMll_Edge > 300'}
    nllCut = {'ttbar': 'nll_Edge < 21',
              'non-ttbar': 'nll_Edge > 21'}
    mt2Cut = {'highmt2': 'mt2_Edge > 80.',
              'lowmt2' : 'mt2_Edge < 80.'}
    Cuts = [cuts.SignalRegionBaseLineNoTrigger, cuts.SF]
    print 'yields for %4.1f /fb'%lumi
    for process in ['TT']:
        tree = treeTT if 'TT' in process else treeDY
        print 'Process', process
        for mt2 in ['highmt2','lowmt2']:
            print mt2
            print '\t & \\multicolumn{3}{c}{ttbar-like} & \t  \\multicolumn{3}{c}{Non-ttbar-like} \\\\'
            for mll in ['low','onz','high1','high2','high3']:
                yields = []
                for nll in ['ttbar','non-ttbar']:
                    theCuts = Cuts + [nllCut[nll]] + [mllCut[mll]] + [mt2Cut[mt2]]
                    yields.append(tree.getYields(lumi,'0.5',0,1,cuts.AddList(theCuts)))
                print '%s\t & %4.1f & +/- & %4.1f\t & %4.1f & +/- &%4.1f \\\\'%(mll,
                                                                                yields[0][0], yields[0][1],
                                                                                yields[1][0], yields[1][1])
