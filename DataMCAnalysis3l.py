#####################################################################
######                                                              #
###### 88888888888         88                        888888888888   #  
###### 88                  88                                 ,88   #
###### 88                  88                               ,88"    #  
###### 88aaaaa     ,adPPYb,88  ,adPPYb,d8  ,adPPYba,      ,88"      #
###### 88"""""    a8"    `Y88 a8"    `Y88 a8P_____88    ,88"        #
###### 88         8b       88 8b       88 8PP"""""""  ,88"          #
###### 88         "8a,   ,d88 "8a,   ,d88 "8b,   ,aa 88"            #
###### 88888888888 `"8bbdP"Y8  `"YbbdP"Y8  `"Ybbd8"' 888888888888   #
######                       aa,    ,88                             #
######                         "Y8bbdP"                             #
######                                                              #
#####################################################################

import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership, TLatex
import math, sys, optparse, array, copy
import gc, inspect

import include.LeptonSF
import include.nll
import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample



class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def dumpObjects():
    gc.collect()
    oo = gc.get_objects()
    for o in oo:
        if getattr(o, "__class__", None):
            name = o.__class__.__name__
            if name not in exclude:
                filename = inspect.getabsfile(o.__class__)            
                print "Object of class:", name, "...",
                print "defined in file:", filename                


def makePlot(lumi, lumi_str, treeMC, sigs, var, name, nbin, xmin, xmax, theCut, cuts,labelx, logx, names):
    normalized = False
    theCutMC = cuts.AddList([theCut, cuts.goodLepton])
    theCutSIG = cuts.AddList([theCut, cuts.goodLeptonSig])
    MC   = treeMC.getTH1F(lumi, "hMC_%s"%(name), var, nbin, xmin, xmax, theCutMC, '', labelx)
    MCS  = treeMC.getStack(lumi, "hMCS_%s"%(name), var, nbin, xmin, xmax, theCutMC, "", labelx)
    SIG1 = sigs[0].getTH1F(lumi, "hSIG1_%s"%(name), var, nbin, xmin, xmax, theCutSIG, '', labelx)
    SIG2 = sigs[1].getTH1F(lumi, "hSIG2_%s"%(name), var, nbin, xmin, xmax, theCutSIG, '', labelx)
    SIG3 = sigs[2].getTH1F(lumi, "hSIG3_%s"%(name), var, nbin, xmin, xmax, theCutSIG, '', labelx)
    SIG4 = sigs[3].getTH1F(lumi, "hSIG4_%s"%(name), var, nbin, xmin, xmax, theCutSIG, '', labelx)

    for _i,_h in enumerate(MCS.GetHists()):
        print "sample ", _h.GetName(), " has integral ", _h.Integral()
        
    if normalized:
        hists = []
        for _i,_h in enumerate(MCS.GetHists()):
            yie_e = r.Double()
            yie = _h.IntegralAndError(1, _h.GetNbinsX()+1, yie_e)
            _h.SetName(_h.GetName().split('_')[-2]+' %.1f +- %.1f'%(yie, yie_e))
            _h.Scale(1./_h.Integral())
            _h.SetLineColor(_h.GetFillColor())
            _h.SetFillColor(0)
            _h.SetLineWidth(2)
            hists.append(_h)
        yie_e = r.Double()
        yie = MC.IntegralAndError(1, MC.GetNbinsX()+1, yie_e)
        MC.SetName('Full MC')
        MC  .Scale(1./MC.Integral())
        MC  .SetLineColor(r.kRed)
        MC  .SetLineWidth(2)         
        hists.append(MC)
    SIG1.SetLineStyle(2);SIG2.SetLineStyle(2);SIG3.SetLineStyle(2);SIG4.SetLineStyle(2);
    SIG1.SetLineWidth(2);SIG2.SetLineWidth(2);SIG3.SetLineWidth(2);SIG4.SetLineWidth(2);
    print "SIG1 ", SIG1.Integral()
    maxVal = MC.GetMaximum() #GetBinContent(MC.GetMaximumBin())
    if not logx:
        MCS.SetMaximum(1.5*maxVal)
    else:
        MCS.SetMaximum(10000*maxVal)

    print name
    SetOwnership(MC, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(MCS, 0 )   # 0 = release (not keep), 1 = keep
    plot = Canvas.Canvas('DataMC/%s/plot_%s'%(lumi_str,name), 'png,pdf,root', 0.65, 0.4, 0.85, 0.9)
    if not normalized:
        plot.addStack(MCS, "HIST", 1, 1)
       # plot.addHisto(SIG1, "HIST,SAME", "m_{l}: %s, m_{{\Chi_{1}^{0}}} : %s" %(names[0][0], names[0][1]), "L", r.kGreen-9, 1, 0)
       # plot.addHisto(SIG2, "HIST,SAME", "m_{l}: %s, m_{{\Chi_{1}^{0}}} : %s" %(names[1][0], names[1][1]), "L", r.kCyan-7, 1, 0)
       # plot.addHisto(SIG3, "HIST,SAME", "m_{l}: %s, m_{{\Chi_{1}^{0}}} : %s" %(names[2][0], names[2][1]), "L", r.kPink-4, 1, 0)
       # plot.addHisto(SIG4, "HIST,SAME", "m_{l}: %s, m_{{\Chi_{1}^{0}}} : %s" %(names[3][0], names[3][1]), "L", r.kYellow-9, 1, 0)
    else:
        for h in hists:
            plot.addHisto(h, 'hist,same' if not h.GetName() == 'data' else 'p,same', h.GetName(), 'PL', h.GetLineColor(), 1, 0)
    plot.saveRatio(1, 1, logx, lumi, MC, MC)


if __name__ == "__main__":

    print bcolors.HEADER
    print '#######################################################################'
    print '                  Starting r_SFOF analysis...                          '
    print '#######################################################################' + bcolors.ENDC
    jet = "0"
    print "Doing ", jet, " jet plots"  
    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples_'+jet+'jetSkim.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientsFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    (opts, args) = parser.parse_args()

    print bcolors.HEADER + '[Data - MC comparisons] ' + bcolors.OKBLUE + 'Loading DATA and MC trees...' + bcolors.ENDC

    mcDatasets = ['ZZTo4L', 'WZTo3LNu', 'WWW', 'WWZ',  'WZZ', 'ZZZ', 'DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO','TTTo2L2Nu_1',  'WWTo2L2Nu', 'ZZTo2L2Nu', 'TTHnobb_pow', 'VHToNonbb', 'TWZ', 'WZTo2L2Q',  'TBar_tch_powheg', 'T_tch_powheg',  'TTJets_SingleLeptonFromT', 'TTTT',  'TTZToQQ', 'TTWToQQ', 'TTLLJets_m1to10', 'TTWToLNu_ext2', 'WJetsToLNu_LO', 'GluGluToContinToZZTo2e2tau', 'GluGluToContinToZZTo2mu2tau', 'GluGluToContinToZZTo2mu2nu', 'GluGluToContinToZZTo4e', 'GluGluToContinToZZTo4mu', 'GluGluToContinToZZTo4tau']

    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0, isScan = 0)
    print bcolors.HEADER + '[Data - MC comparisons] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    lumi = 35.9 ; maxrun = 999999; lumi_str = '35.9invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    cuts = CutManager.CutManager()
    labelmll = 'm_{ll} [GeV]'
    labelmet = 'E_{T}^{miss} [GeV]'
    labelnjet = "N. Jets"
    njetCut = "nJetSel_Edge == "+jet
    signals = ['SMS_400_150__Chunk1978901673', 'SMS_300_10__Chunk1681815997', 'SMS_250_70__Chunk771063984', 'SMS_150_40__Chunk9117347303']
    signames = [[400, 150], [300, 10], [250, 70], [150, 40]]
    xsec = [2*(1.310+ 0.500), 2*(4.430+ 1.680),2*(9.210+3.470),2*(63.340+23.320) ]
    treeSIG1 = Sample.Tree(helper.selectSamples(opts.sampleFile, [signals[0]], 'SI'), 'SI', 0, isScan = xsec[0])
    treeSIG2 = Sample.Tree(helper.selectSamples(opts.sampleFile, [signals[1]], 'SI'), 'SI', 0, isScan = xsec[1])
    treeSIG3 = Sample.Tree(helper.selectSamples(opts.sampleFile, [signals[2]], 'SI'), 'SI', 0, isScan = xsec[2])
    treeSIG4 = Sample.Tree(helper.selectSamples(opts.sampleFile, [signals[3]], 'SI'), 'SI', 0, isScan = xsec[3])
    sigs = [treeSIG1, treeSIG2, treeSIG3, treeSIG4]
    makePlot(lumi, lumi_str, treeMC, sigs, "Lep1_pt_Edge", "lep1_pt_"+jet+"_jet_3l", 20,0, 200, cuts.AddList([cuts.slep3lCR, njetCut]), cuts, "Lep 1 p_{T} [GeV] ("+ jet+" jet)",  0, signames)
    makePlot(lumi, lumi_str, treeMC, sigs, "Lep2_pt_Edge", "lep2_pt_"+jet+"_jet_3l", 20,0, 200, cuts.AddList([cuts.slep3lCR, njetCut]), cuts, "Lep 2 p_{T} [GeV] ("+ jet+" jet)",  0, signames)
    makePlot(lumi, lumi_str, treeMC, sigs, "minMT_Edge", "minMT_"+jet+"_jet_3l", 40,0, 400, cuts.AddList([cuts.slep3lCR, njetCut]), cuts, "min m_{T} [GeV] ("+ jet+" jet)",  0, signames)
    makePlot(lumi, lumi_str, treeMC, sigs, "met_Edge","met_"+jet+"_jet_3l", 30, 100, 400, cuts.AddList([cuts.slep3lCR, njetCut]), cuts, "E_{T}^{miss} [GeV] ("+ jet+" jet)",  0, signames)
    makePlot(lumi, lumi_str, treeMC, sigs, "lepsMll_Edge", "mll_"+jet+"_jet_3l", 35, 1, 351, cuts.AddList([cuts.slep3lCR, njetCut]), cuts, "m_{ll} [GeV] ("+ jet+" jet)",  0, signames)
    makePlot(lumi, lumi_str, treeMC, sigs, "abs(lepsDPhi_Edge)", "lepsDPhi_"+jet+"_jet_3l", 32, 0, 3.2, cuts.AddList([cuts.slep3lCR, njetCut]), cuts, "|#Delta#Phi(l,l)| ("+ jet+" jet)",  0, signames)
    makePlot(lumi, lumi_str, treeMC, sigs, "abs(metl1DPhi_Edge)", "dPhi_lep1_met_"+jet+"_jet_3l", 32, 0, 3.2, cuts.AddList([cuts.slep3lCR, njetCut]),cuts,"|#Delta#Phi(E_{T}^{miss},l_{1})| ("+ jet+" jet)",  0, signames)
    makePlot(lumi, lumi_str, treeMC, sigs, "abs(metl2DPhi_Edge)", "dPhi_lep2_met_"+jet+"_jet_3l", 32, 0, 3.2, cuts.AddList([cuts.slep3lCR, njetCut]),cuts,"|#Delta#Phi(E_{T}^{miss},l_{2})| ("+ jet+" jet)",  0, signames)
    makePlot(lumi, lumi_str, treeMC, sigs, "lepsMETRec_Edge", "recoil_"+jet+"_jet_3l", 30, 0, 300, cuts.AddList([cuts.slep3lCR, njetCut]), cuts, "recoil [GeV] ("+ jet+" jet)",  0, signames)
    makePlot(lumi, lumi_str, treeMC, sigs, "mt2_Edge", "mt2_"+jet+"_jet_3l", 30, 100, 400, cuts.AddList([cuts.slep3lCR, njetCut]), cuts, "M_{T2} [GeV] ("+ jet+" jet)",  0, signames)
    makePlot(lumi, lumi_str, treeMC, sigs, "nLepLoose_Edge", "nLepLoose_"+jet+"_jet_3l", 3, 1.5, 4.5, cuts.AddList([njetCut,cuts.slep3lCR]), cuts, "number of loose leptons ("+ jet+" jet)",  0, signames)
    if jet == "1":
        makePlot(lumi, lumi_str, treeMC, sigs, "JetSel_Edge_pt", "jetPT_"+jet+"_jet_3l", 40, 0, 400, cuts.AddList([cuts.slep3lCR, njetCut]), cuts, "Jet p_{T} [GeV] ("+ jet+" jet)",  0, signames)
        makePlot(lumi, lumi_str, treeMC, sigs, "abs(j1MetDPhi_Edge)", "j1MetDPhi_"+jet+"_jet_3l", 32, 0, 3.2, cuts.AddList([cuts.slep3lCR, njetCut]), cuts, "|#Delta#Phi(jet, E_{T}^{miss})| ("+ jet+" jet)",  0, signames)



