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
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership
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

def compareMMEEPlots(plotmm, plotee):
    for _i,h in enumerate(plotmm.histos):
        if 'hDATA' in h.GetName():
            ind = _i
    h_mm = plotmm.histos[_i]
    h_ee = plotee.histos[_i]
    name = plotmm.name.split('/')[-1].replace('plot_','').replace('mm_','')
    lumistr = plotmm.name.split('/')[-2]

    plot = Canvas.Canvas('DataMC/%s/plot_mmeeComparison_%s'%(lumistr,name), 'png,pdf,root', 0.7, 0.7, 0.9, 0.9)
    plot.addHisto(h_ee, "PE1,SAME", "elel", "P", r.kYellow+1, 1, 1)
    plot.addHisto(h_mm, "PE1,SAME", "mumu", "P", r.kBlack   , 1, 0)
    plot.saveRatio(1, 1, 0, 4, h_ee, h_mm)
    

def makeSummaryTable3l4l(plot_3l, plot_4l, plot_ttZ):
    wz3l, wz3l_e, zz4l, zz4l_e, ttZ, ttZ_e = 0., 0., 0., 0.,0., 0.
    bkg3l, bkg3l_e, bkg4l, bkg4l_e, bkgttZ, bkgttZ_e = 0., 0., 0., 0.,0., 0.
    obs3l, obs4l, obsttZ = 0, 0, 0
    rat3l, rat3l_e, rat4l, rat4l_e, ratttZ, ratttZ_e = 0., 0., 0., 0.,0., 0.
    
    for i, h_3l in enumerate(plot_3l.histos):
        if not i: continue
        tmp_double = r.Double()
        if 'WZ' in h_3l.GetName():
            wz3l = h_3l.IntegralAndError(1, h_3l.GetNbinsX()+1, tmp_double)
            wz3l_e = tmp_double
        elif not 'DATA' in h_3l.GetName(): 
            tmp_yield = h_3l.IntegralAndError(1, h_3l.GetNbinsX()+1, tmp_double)
            bkg3l  += tmp_yield
            bkg3l_e = math.sqrt(bkg3l_e**2 + tmp_double**2)
        elif 'DATA' in h_3l.GetName():
            obs3l = h_3l.Integral()
        
    for i, h_4l in enumerate(plot_4l.histos):
        if not i: continue
        tmp_double = r.Double()
        print "h_4l.GetName() ", h_4l.GetName()
        if 'ZZ' in h_4l.GetName():
            zz4l = h_4l.IntegralAndError(1, h_4l.GetNbinsX()+1, tmp_double)
            zz4l_e = tmp_double
        elif not 'DATA' in h_4l.GetName(): 
            tmp_yield = h_4l.IntegralAndError(1, h_4l.GetNbinsX()+1, tmp_double)
            bkg4l  += tmp_yield
            bkg4l_e = math.sqrt(bkg4l_e**2 + tmp_double**2)
        elif 'DATA' in h_4l.GetName():
            obs4l = h_4l.Integral()                                                      
    
    for i, h_ttZ in enumerate(plot_ttZ.histos):
        if not i: continue
        tmp_double = r.Double()
        print "h_ttZ.GetName() ", h_ttZ.GetName()
        if 'TTZ' in h_ttZ.GetName():
            ttZ = h_ttZ.IntegralAndError(1, h_ttZ.GetNbinsX()+1, tmp_double)
            ttZ_e = tmp_double
        elif not 'DATA' in h_ttZ.GetName(): 
            print "now here!!! "
            tmp_yield = h_ttZ.IntegralAndError(1, h_ttZ.GetNbinsX()+1, tmp_double)
            bkgttZ  += tmp_yield
            bkgttZ_e = math.sqrt(bkgttZ_e**2 + tmp_double**2)
        elif 'DATA' in h_ttZ.GetName():
            print "here obs!!!"
            obsttZ = h_ttZ.Integral()     


    obsSub3l, obsSub3l_e, obsSub4l, obsSub4l_e = obs3l-bkg3l, math.sqrt(bkg3l_e**2 + obs3l), obs4l-bkg4l, math.sqrt(bkg4l_e**2 + obs4l)
    obsSubttZ, obsSubttZ_e = obsttZ-bkgttZ, math.sqrt(bkgttZ_e**2 + obsttZ)
    rat3l, rat3l_e = helper.ratioError(obsSub3l, obsSub3l_e, wz3l, wz3l_e)
    rat4l, rat4l_e = helper.ratioError(obsSub4l, obsSub4l_e,zz4l, zz4l_e)
    ratttZ, ratttZ_e = helper.ratioError(obsSubttZ, obsSubttZ_e,ttZ, ttZ_e)

    table3l4lComparisonString = '''\\documentclass[12pt,a4paper]{{article}}
\\begin{{document}}
\\begin{{table}}[hbtp] 
\\begin{{center}} 
\\bgroup 
\\def\\arraystretch{{1.2}} 
\\caption{{3-lepton, 4-lepton and ttZ control regions. Signal MC is WZ-$>$3lnu in the 3-lepton region, ZZ-$>$4l in the 4-lepton region and ttZ-$>$2l2nu and ttZ-$>$qq in the ttZ region}} 
\\label{{tab:tab3l4lcontrol}} 
\\begin{{tabular}}{{l| c c c}} 
              & 3-lepton region & 4-lepton region &  ttZ region\\\\ \\hline \\hline
signal MC        & {wz3l:.2f}     $\\pm$  {wz3l_e:.2f}       & {zz4l:.2f}     $\\pm$  {zz4l_e:.2f}    & {ttZ:.2f}     $\\pm$  {ttZ_e:.2f}     \\\\ \\hline
bkg. MC          & {bkg3l:.2f}  $\\pm$  {bkg3l_e:.2f}        & {bkg4l:.2f}  $\\pm$  {bkg4l_e:.2f}    &   {bkgttZ:.2f}  $\\pm$  {bkgttZ_e:.2f}\\\\ \\hline \\hline
\\textbf{{data}}             & \\textbf{{{obs3l}}}                               & \\textbf{{{obs4l}}}    & \\textbf{{{obsttZ}}}   \\\\ \\hline
data-bkg.        & {obsSub3l:.2f}   $\\pm$  {obsSub3l_e:.2f} & {obsSub4l:.2f}   $\\pm$  {obsSub4l_e:.2f}   & {obsSubttZ:.2f}   $\\pm$  {obsSubttZ_e:.2f}    \\\\ \\hline \\hline
(data-bkg.)/sig. & {rat3l:.2f}   $\\pm$  {rat3l_e:.2f}       & {rat4l:.2f}   $\\pm$  {rat4l_e:.2f}     & {ratttZ:.2f}   $\\pm$  {ratttZ_e:.2f}    \\\\ \\hline

\\end{{tabular}} 
\\egroup 
\\end{{center}} 
\\end{{table}} 
\\end{{document}}'''.format(wz3l=wz3l, wz3l_e=wz3l_e, zz4l=zz4l, zz4l_e=zz4l_e,ttZ=ttZ, ttZ_e=ttZ_e,
                            bkg3l=bkg3l, bkg3l_e=bkg3l_e, bkg4l=bkg4l,bkg4l_e=bkg4l_e, bkgttZ=bkgttZ,bkgttZ_e=bkgttZ_e,
                            obs3l=int(obs3l),obs4l=int(obs4l),obsttZ=int(obsttZ),
                            obsSub3l=obsSub3l, obsSub3l_e=obsSub3l_e, obsSub4l=obsSub4l, obsSub4l_e=obsSub4l_e,obsSubttZ=obsSubttZ, obsSubttZ_e=obsSubttZ_e,
                            rat3l=rat3l, rat3l_e=rat3l_e, rat4l=rat4l, rat4l_e=rat4l_e, ratttZ=ratttZ, ratttZ_e=ratttZ_e)
    compTableFile = open('plots/DataMC/12.9invfb/tables/controlRegion3l4l_comparison.txt','w')
    compTableFile.write(table3l4lComparisonString)
    compTableFile.close()

def makePlot(lumi, lumi_str, treeDA, treeMC, var, name, nbin, xmin, xmax, theCut, cuts,labelx, logx, returnplot = False, addHistos = False , scaleToData = False, normalized = False, cumulative = False, onlyMC = False):


    theCutTau = cuts.AddList([theCut, cuts.bothTau ])
    theCutNoTau = cuts.AddList([theCut, cuts.bothNoTau ])
    MCTAU   = treeMC.getTH1F(lumi, "hMC_Tau%s"%(name), var, nbin, xmin, xmax, theCutTau, '', labelx)
    MCNOTAU = treeMC.getTH1F(lumi, "hMC_NoTau%s"%(name), var, nbin, xmin, xmax, theCutNoTau, '', labelx)
    #MClep1Tau2 = treeMC.getTH1F(lumi, "hMC_lep1Tau1%s"%(name), var, nbin, xmin, xmax, theCutlep1Tau1, '', labelx)
    mcMAX = MCTAU.GetMaximum() 
    mcMAX1 = MCNOTAU.GetMaximum() 
    #maxMClep1Tau2 = MClep1Tau2.GetMaximum()
    maxVal = max(mcMAX1,  mcMAX)
    MCNOTAU.GetYaxis().SetRangeUser(0.1, maxVal*1.5)
    plot = Canvas.Canvas('DataMC/%s/plot_%s_TauMatch'%(lumi_str,name), 'png,pdf,root', 0.65, 0.7, 0.85, 0.9)
    plot.addHisto(MCNOTAU, 'hist', 'both no tau', 'L', r.kGreen+3, 1, 0)
    #plot.addHisto(MClep1Tau2, 'hist same', '1 lep from tau ', 'PL', r.kOrange+5, 1, 0)
    plot.addHisto(MCTAU, 'hist same', 'both tau', 'L', r.kOrange+5, 1, 0)
    plot.save(1, 1, logx, lumi)

    del plot


if __name__ == "__main__":

    print bcolors.HEADER
    print '#######################################################################'
    print '                  Starting r_SFOF analysis...                          '
    print '#######################################################################' + bcolors.ENDC
    jet = "1"
    print "Doing ", jet, " jet plots"  
    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples_'+jet+'jetSkim.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientsFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    (opts, args) = parser.parse_args()

    print bcolors.HEADER + '[Data - MC comparisons] ' + bcolors.OKBLUE + 'Loading DATA and MC trees...' + bcolors.ENDC

    #mcDatasets = ['ZZTo4L','GGHZZ4L','VVTo2L2Nu_ext',   'WZTo3LNu', 'WWW', 'WWZ',  'WZZ', 'ZZZ', 'DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO', 'TTTo2L2Nu','tZq_ll','TTHnobb_pow', 'VHToNonbb', 'TWZ', 'WZTo2L2Q',  'TBar_tch_powheg', 'T_tch_powheg',  'TTJets_SingleLeptonFromT', 'TTTT',  'TTZToQQ', 'TTWToQQ', 'TTLLJets_m1to10', 'TTZToLLNuNu_ext2' ,'TTWToLNu_ext2', 'WJetsToLNu_LO']
    #mcDatasets = ['ZZTo2L2Nu']
    mcDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO']
    #mcDatasets = ['ZZTo4L', 'WZTo3LNu', 'WWW', 'WWZ',  'WZZ', 'ZZZ', 'DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO', 'TTTo2L2Nu_1','WWTo2L2Nu', 'ZZTo2L2Nu', 'TTHnobb_pow', 'VHToNonbb', 'TWZ', 'WZTo2L2Q',  'TBar_tch_powheg', 'T_tch_powheg',  'TTJets_SingleLeptonFromT', 'TTTT',  'TTZToQQ', 'TTWToQQ', 'TTLLJets_m1to10', 'TTWToLNu_ext2', 'WJetsToLNu_LO', 'GluGluToContinToZZTo2e2tau', 'GluGluToContinToZZTo2mu2tau', 'GluGluToContinToZZTo2mu2nu', 'GluGluToContinToZZTo4e', 'GluGluToContinToZZTo4mu', 'GluGluToContinToZZTo4tau']
    #mcDatasets = ['ZZTo4L','GGHZZ4L',  'WZTo3LNu', 'WWW', 'WWZ',  'WZZ', 'ZZZ', 'DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO', 'TTTo2L2Nu','tZq_ll','WWTo2L2Nu', 'ZZTo2L2Nu', 'TTHnobb_pow', 'VHToNonbb', 'TWZ', 'WZTo2L2Q',  'TBar_tch_powheg', 'T_tch_powheg',  'TTJets_SingleLeptonFromT', 'TTTT',  'TTZToQQ', 'TTWToQQ', 'TTLLJets_m1to10', 'TTZToLLNuNu_ext2' ,'TTWToLNu_ext2', 'WJetsToLNu_LO']

    daDatasetsB = ['DoubleEG_Run2016B_03Feb2017_ver2_v2_runs_273150_275376',                            
                   'DoubleMuon_Run2016B_03Feb2017_ver2_v2_runs_273150_275376',
                   'MuonEG_Run2016B_03Feb2017_ver2_v2_runs_273150_275376']    
                                                                              
    daDatasetsC = ['DoubleEG_Run2016C_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016C_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016C_03Feb2017_v1_runs_271036_284044']    
    
    daDatasetsD = ['DoubleEG_Run2016D_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016D_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016D_03Feb2017_v1_runs_271036_284044']    
                                                                              
    daDatasetsE = ['DoubleEG_Run2016E_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016E_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016E_03Feb2017_v1_runs_271036_284044']    
                                                                              
    daDatasetsF = ['DoubleEG_Run2016F_03Feb2017_v1_runs_271036_284044',
                  'DoubleMuon_Run2016F_03Feb2017_v1_runs_271036_284044',
                  'MuonEG_Run2016F_03Feb2017_v1_runs_271036_284044']  
                                                                              
    daDatasetsG = ['DoubleEG_Run2016G_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016G_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016G_03Feb2017_v1_runs_271036_284044']    
                                                                              
    daDatasetsH = ['DoubleEG_Run2016H_03Feb2017_ver2_v1_runs_281085_284035',
                   'DoubleEG_Run2016H_03Feb2017_ver3_v1_runs_284036_284044',
                   'DoubleMuon_Run2016H_03Feb2017_ver2_v1_runs_281085_284035',
                   'DoubleMuon_Run2016H_03Feb2017_ver3_v1_runs_284036_284044',
                   'MuonEG_Run2016H_03Feb2017_ver2_v1_runs_281085_284035', 
                   'MuonEG_Run2016H_03Feb2017_ver3_v1_runs_284036_284044']    

 
    daDatasets = daDatasetsB + daDatasetsC + daDatasetsD +daDatasetsE + daDatasetsF + daDatasetsG + daDatasetsH       


    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0, isScan = 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1, isScan = 0)
    #treeSI = Sample.Tree(helper.selectSamples(opts.sampleFile, siDatasets, 'SI'), 'SI', 0, isScan = 1)
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
    #makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_OFSF_nll", [20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300], 1, 1, cuts.AddList([cuts.goodLepton, cuts.SignalRegionBaseLine]), cuts, labelmll, 0)

##########################################Stuff from Leonora##########################    
    ## EWino signal region and table
    ## =================================================================
    #ewino_SR = makePlot(12.9, '12.9invfb', treeDA, treeMC, "met_Edge", "met_ewino_SR", 12, 0, 300, cuts.AddList([cuts.ewinoSR, 'run_Edge <= 999999']), cuts, labelmet, 0, True)
    #makeSummaryEWino(ewino_SR)
    
    #plot_nll_sf     = makePlot(lumi, lumi_str, treeDA, treeMC, "nll_Edge", "nll_edgeBaseline"    , 14, 13., 27, cuts.AddList([cuts.goodLepton, cuts.SignalRegionBaseLine, cuts.SF, cuts.Zveto ]), cuts, 'NLL', 0, True, False, True, True, True)
    #plot_nll_sf_onz = makePlot(lumi, lumi_str, treeDA, treeMC, "nll_Edge", "nll_edgeBaseline_onZ", 14, 13., 27, cuts.AddList([cuts.goodLepton, cuts.SignalRegionBaseLine, cuts.SF, cuts.Zmass ]), cuts, 'NLL', 0, True, False, True, True, True)
    #plot_nll_sf_lm  = makePlot(lumi, lumi_str, treeDA, treeMC, "nll_Edge", "nll_edgeBaseline_loM", 14, 13., 27, cuts.AddList([cuts.goodLepton, cuts.SignalRegionBaseLine, cuts.SF, cuts.loMass]), cuts, 'NLL', 0, True, False, True, True, True)
    #plot_nll_sf_hm  = makePlot(lumi, lumi_str, treeDA, treeMC, "nll_Edge", "nll_edgeBaseline_hiM", 14, 13., 27, cuts.AddList([cuts.goodLepton, cuts.SignalRegionBaseLine, cuts.SF, cuts.hiMass]), cuts, 'NLL', 0, True, False, True, True, True)

    #makePlot(lumi, lumi_str, treeSI, treeMC, "bestMjj_Edge", "mjj_ewinoSR", 20,  0,  300, cuts.AddList([cuts.ewinoSR]), cuts, 'm_{jj}', 0, True)
    #makePlot(lumi, lumi_str, treeSI, treeMC, "bestMjj_Edge", "mjj", 20,  0,  300, cuts.AddList([cuts.goodLepton]), cuts, 'm_{jj}', 0, True)
    ## 3l and 4l plots
    ## =================================================================
    #plot_ttZ = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttZregion_AF_allSamples", 10,  0, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.regionttZ]), cuts, labelmll, 0, True)
    #plot_ttZ_met = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttZregion_AF_allSamples", 20,  0, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.regionttZ, cuts.Zmass]), cuts, labelmet, 0, True)
    #plot_3l = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge"    , "met_3lregion", 14, 60, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region3l]), cuts, labelmet, 0, True)
    #plot_3l_Zmass = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge"    , "met_3lregion_Zmass", 14, 60, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region3l, cuts.Zmass]), cuts, labelmet, 0, True)
    #plot_4l = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_4lregion", 20,  6, 206, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region4l]), cuts, labelmll, 0, True)
    #plot_4l_Zmass = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_4lregion_Zmass", 20,  6, 206, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region4l, cuts.Zmass]), cuts, labelmll, 0, True)
    #plot_ttZ = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttZregion", 20,  6, 206, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.regionttZ]), cuts, labelmll, 0, True)
    #plot_ttZ_Zmass = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttZregion_Zmass", 20,  6, 206, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.regionttZ, cuts.Zmass]), cuts, labelmll, 0, True)
    ## ## ## plot_4l = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_4lregion_met0to30_AF", 10,  0, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region4l, 'met_Edge < 30']), cuts, labelmll, 0, True)
    #makeSummaryTable3l4l(plot_3l_Zmass, plot_4l_Zmass, plot_ttZ_Zmass)
    #print "made the table"
    #ttbar_region = cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj2, 'run_Edge <= %d'%maxrun, 'met_Edge >150'])
    #makePlot(lumi, lumi_str, treeDA, treeMC, "nBJetMedium35_Edge", "ttbar_region_met150_2jets_of_dataScaled", 4, 0, 4, ttbar_region, cuts, 'n_{b-jets}', 0, False, True)
##########################################End Stuff from Leonora##########################    
    #makePlot(lumi, lumi_str, treeMC, treeMC, "mT", "mT_"+jet, 30, 0, 300, cuts.AddList([cuts.goodLepton, cuts.slep, cuts.ThirdLeptonVeto]), cuts, "m_{T} [GeV] ("+ jet+" jet)",  0, True, True)
    #bins =[100, 120, 140, 160, 180, 200, 240, 280]
    bins =[20, 40, 60 , 80, 100, 120, 140, 160, 180, 200, 240, 280]
    makePlot(lumi, lumi_str, treeMC, treeMC, "mt2_Edge", "mt2_"+jet+"_onlyDY", bins, 1, 1, cuts.AddList([cuts.goodLepton, njetCut, cuts.slep, cuts.ThirdLeptonVeto]), cuts, "M_{T2} [GeV] ("+ jet+" jet)",  0, True, False)
   # makePlot(lumi, lumi_str, treeMC, treeMC, "lepsMll_Edge", "mll_SR", 30, 0, 300, cuts.AddList([cuts.goodLepton, cuts.SignalRegion,  cuts.SF]), cuts, 'mll',  0, True)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll", 30, 0, 300, cuts.AddList([cuts.goodLepton, cuts.DYControlRegionNoMll,  cuts.ee]), cuts, 'nll',  0, True)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "mbb_Edge", "mbb_SF", 20, 0, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.ewinoNeuNeuNomt2bbmbb]), cuts, 'mbb',  0, True)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "mt2bb_Edge", "mt2bb_SF", 20, 0, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.ewinoNeuNeuNomt2bbmbb]), cuts, 'mt2bb',  0, True)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_OF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.DYControlRegionNoMll]), cuts, labelmll,  0, True)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_SF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.DYControlRegionNoMll]), cuts, labelmll,  0, True)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_ee", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.DYControlRegionNoMll]), cuts, labelmll,  0, True)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_mm", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.DYControlRegionNoMll]), cuts, labelmll,  0, True)
   # 
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_SF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll,  0, True)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_OF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll,  0, True)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_ee", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll,  0, True)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_mm", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll,  0, True)
 
   # makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_SF", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.SF]), cuts, labelnjet, 0, True)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_OF", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.OF]), cuts, labelnjet, 0, True)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_ee", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.ee]), cuts, labelnjet, 0, True)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_mm", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.mm]), cuts, labelnjet, 0, True)                                              


   # #This one needs protection
 
   # makePlot(17.0, '17invfb', treeDA, treeMC, "met_Edge", "met_ttbarcontrol_SF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.RSFOFDirectControlRegionNoMET, cuts.protection]), cuts, labelmet, 0, True)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge",  "met_ttbarcontrol_OF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet,  0, True)
   # makePlot(17.0, '17invfb', treeDA, treeMC, "met_Edge", "met_ttbarcontrol_ee", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.RSFOFDirectControlRegionNoMET, cuts.protection]), cuts, labelmet,  0, True)
   # makePlot(17.0, '17invfb', treeDA, treeMC, "met_Edge", "met_ttbarcontrol_mm", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.RSFOFDirectControlRegionNoMET, cuts.protection]), cuts, labelmet,  0, True)













    # =================================================================

    # r0b1b_of_cuts = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge >= 0', cuts.nj2, 'j1MetDPhi_Edge > 1.', 'lepsMll_Edge > 61 && lepsMll_Edge < 121', 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.'])
    # r0b1b_of_den_cuts_e1b = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge == 1', cuts.nj2, 'j1MetDPhi_Edge > 1.', 'lepsMll_Edge > 61 && lepsMll_Edge < 121', 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.'])
    # #makePlot(lumi, lumi_str, treeDA, treeMC, "nBJetMedium25_Edge", "r0b1b_of_nb_btagscaledPOWHEG", 4, -0.5, 3.5, r0b1b_of_cuts, cuts, 'n_{b-jets}', 0, True, True)
    # #makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "r0b1b_of_den_e1b", 14, 60, 130, r0b1b_of_den_cuts_e1b, cuts, labelmll, 0, True)

    # fmll_of_1b_cuts = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge >= 1', cuts.nj2, 'j1MetDPhi_Edge > 1.', 'lepsMll_Edge > 61 && lepsMll_Edge < 121', 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.'])
    # fmll_of_0b_cuts = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge == 0', cuts.nj2, 'j1MetDPhi_Edge > 1.', 'lepsMll_Edge > 61 && lepsMll_Edge < 121', 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.'])
    # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "fmll_of_1b", 14, 60, 200, fmll_of_1b_cuts, cuts, labelmll, 0, True)
    # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "fmll_of_0b", 14, 60, 200, fmll_of_0b_cuts, cuts, labelmll, 0, True)

    #ewino_of_eq1b = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge == 1', cuts.nj2, 'j1MetDPhi_Edge > 1.', 'lepsMll_Edge > 61 && lepsMll_Edge < 121', 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.'])
    #a =makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "ewino_of_1b_eq1b", 14, 60, 130, ewino_of_eq1b, cuts, labelmll, 0, True)
    # ewino_sf = cuts.AddList([cuts.goodLepton, cuts.SF, 'nBJetMedium25_Edge >= 1', cuts.nj2, 'j1MetDPhi_Edge > 1.', 'lepsMll_Edge > 61 && lepsMll_Edge < 121', cuts.Zveto, 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.'])
    # c = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "ewino_sf_1b_zveto", 14, 60, 130, ewino_sf, cuts, labelmll, 0, True)
    # ewino_of = cuts.AddList([cuts.goodLepton, cuts.OF, 'nBJetMedium25_Edge >= 1', cuts.nj2, 'j1MetDPhi_Edge > 1.', 'lepsMll_Edge > 61 && lepsMll_Edge < 121', 'run_Edge <= %d'%maxrun, 'mt2_Edge > 80.'])
    # b = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "ewino_of_1b", 14, 60, 130, ewino_of, cuts, labelmll, 0, True)

   # nll = '(nll_Edge > 21)'
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_SF_nll", [20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300], 1, 1, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.SignalRegionBaseLine, nll]), cuts, labelmll, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_OF_nll", [20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300], 1, 1, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.SignalRegionBaseLine, nll]), cuts, labelmll, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_OFSF_nll", [20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300], 1, 1, cuts.AddList([cuts.goodLepton, cuts.SignalRegionBaseLine]), cuts, labelmll, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_ee", [20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300], 1, 1, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.SignalRegionBaseLine]), cuts, labelmll, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_mm", [20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300], 1, 1, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.SignalRegionBaseLine]), cuts, labelmll, 0)
   # tightCharge = '(Lep1_tightCharge_Edge > 0 && Lep2_tightCharge_Edge > 0)'
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_SF_nll_tightCharge_ttpow", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.SignalRegionBaseLine, nll, tightCharge]), cuts, labelmll, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_OF_nll_tightCharge_ttpow", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.SignalRegionBaseLine, nll, tightCharge]), cuts, labelmll, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_ee_nll_tightCharge_ttpow", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.SignalRegionBaseLine, nll, tightCharge]), cuts, labelmll, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_mm_nll_tightCharge_ttpow", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.SignalRegionBaseLine, nll, tightCharge]), cuts, labelmll, 0)
   # makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_signal_SF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj2]), cuts, labelmet, 1)



