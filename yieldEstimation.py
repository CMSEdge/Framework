import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TH1F, TGraph
import math, sys, optparse, copy, re, array, itertools, os


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Tables     as Tables
import include.SMSReweight as SMSReweight



if __name__ == "__main__":
   

    gROOT.SetBatch(True)
    gROOT.ProcessLine('.L include/tdrstyle.C')
    r.setTDRStyle() 
    r.gStyle.SetPaintTextFormat("4.2f")
    
    
    global lumi
    lumi= 35.
    mcDatasets = ['TT_pow_ext34']
    dyDatasets = ['DYJetsToLL_M10to50','DYJetsToLL_M50_HT100to200_ext','DYJetsToLL_M50_HT200to400_ext',
                  'DYJetsToLL_M50_HT400to600_ext','DYJetsToLL_M50_HT600toInf']

    siDatasets = ['SMS_T6bbllslepton_mSbottom400to575_mLSP150to550',
                  'SMS_T6bbllslepton_mSbottom600to775_mLSP150to725',
                  'SMS_T6bbllslepton_mSbottom800to950_mLSP150to900']

    treeTT = Sample.Tree(helper.selectSamples('samples.dat', mcDatasets, 'MC'), 'MC', 0, isScan = 0)
    treeDY = Sample.Tree(helper.selectSamples('samples.dat', dyDatasets, 'MC'), 'MC', 0, isScan = 0)
    treeSI = Sample.Tree(helper.selectSamples('samples.dat', siDatasets, 'SI'), 'SI', 0, isScan = 1)

    cuts = CutManager.CutManager()
    mt2Cut = {'highmt2': 'mt2_Edge > 80.',
              'lowmt2' : 'mt2_Edge < 80.'}
    Cuts = [cuts.SignalRegionBaseLineNoTrigger, cuts.SF]
    SMSRW =  SMSReweight.SMSReweight('T6bbllslepton')

    print 'yields for %4.1f /fb'%lumi
    for process in ['T6bbllslepton_900_200','T6bbllslepton_900_400','T6bbllslepton_900_600','TT']:
        print 'Process', process
        for mt2 in ['highmt2','lowmt2']:
            print mt2
            theCuts = Cuts + [mt2Cut[mt2]]
            if 'T6bbllslepton' in process:
                theCuts = theCuts + ['GenSusyMScan1_Edge == %s && GenSusyMScan2_Edge == %s'%(process.split('_')[1],
                                         process.split('_')[2])]
                tree = treeSI
        
            else:
                tree = treeTT if 'TT' in process else treeDY
        
            theCuts = cuts.AddList(theCuts)
            reWeight = 1
            if 'T6bbllslepton' in process:
                theCuts = theCuts.replace(cuts.twoLeptons, 'nPairLep_Edge > 0')
                reWeight = SMSRW.getReWeightFactor(float(process.split('_')[1]),float(process.split('_')[2]))


            h = tree.getTH2F(lumi,'yields%s'%process,'lepsMll_Edge:nll_Edge',[0,21.,100000.],0,0,[20.,81.,101.,200.,300.,10000.],0,0,theCuts,"","","")
            h.Scale(reWeight)
            print '\t & \\multicolumn{3}{c}{ttbar-like} & \t  \\multicolumn{3}{c}{Non-ttbar-like} \\\\'
            for i,mll in enumerate(['low','onz','high1','high2','high3']):
                print '%s\t & %4.1f & +/- & %4.1f\t & %4.1f & +/- &%4.1f \\\\'%(mll,
                                                                                h.GetBinContent(1,i+1), h.GetBinError(1,i+1),
                                                                                h.GetBinContent(2,i+1), h.GetBinError(2,i+1))
