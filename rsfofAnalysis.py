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

import include.LeptonSF
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




def saveInFile(theFile, measuredValueMC, measuredValueUncMC, measuredValueData, measuredValueUncData):

    foutput = open(theFile + "_aux", 'w')
    for line in open(theFile).readlines():
        if line.find("rsfof") != -1 and line.find("direct") != -1:
            if line.find("DATA") != -1:
                foutput.write('rsfof       direct          DATA        %.4f      %0.4f       %.4f\n'%(measuredValueData, measuredValueUncData, measuredValueUncMC))
            else:
                foutput.write('rsfof       direct          MC          %.4f      %0.4f       %.4f\n'%(measuredValueMC, measuredValueUncMC, 0))
        else:
            foutput.write(line)

    foutput.close()
    subprocess.call(['mv ' + theFile + "_aux " + theFile], shell=True)                                                                                                        
  

def makeTable(DATAee, DATAmm, DATASF, DATAOF, MCee, MCmm, MCSF, MCOF,DATArsfof, DATArsfof_e, MCrsfof, MCrsfof_e):
    line0 = '  \hline'
    line1 = '  ee &' 
    line2 = '  \mu\mu   &' 
    line3 = '  OF   &' 
    line4 = '  SF   &' 
    line5 = '  R_{\mathrm{SF/OF}}   &' 
    line0 += ' & Data  & MC   \\\\'
    line1 += ' %.f & %.f   \\\ ' %(DATAee.Integral(), MCee.Integral())
    line2 += ' %.f & %.f   \\\ ' %(DATAmm.Integral(), MCmm.Integral())
    line3 += ' %.f & %.f   \\\ ' %(DATAOF.Integral(), MCOF.Integral())
    line4 += ' %.f & %.f   \\\ ' %(DATASF.Integral(), MCSF.Integral())
    line5 += ' %.3f \pm %.3f & %.3f \pm %.3f  \\\ ' %(DATArsfof,DATArsfof_e, MCrsfof, MCrsfof_e)
    line0 += '\\hline'; line4 += '\\hline'; line5 += '\\hline';
                                                                                                                                                                                     
    helper.ensureDirectory('plots/rsfof/%s/'%lumi_str)
    helper.ensureDirectory('plots/rsfof/%s/tables/'%lumi_str)
    compTableFile = open('plots/rsfof/%s/tables/resultTable_%s%s.txt'%(lumi_str, str(lumi).replace('.','p'), "rt"),'w')
    compTableFile.write(line0+'\n')
    compTableFile.write(line1+'\n')
    compTableFile.write(line2+'\n')                                                                                             
    compTableFile.write(line3+'\n')                                                                                             
    compTableFile.write(line4+'\n')                                                                                             
    compTableFile.write(line5+'\n')                                                                                             
    print line0
    print line1
    print line2                                                                                                                                                                      
    print line3                                                                                                                                             
    print line4                                                                                                                                                      
    print line5                                                                                                                                                     


                                                                                                                                                                              
def make_rsfof(histo_sf, histo_of, dataMC):

    ratio = histo_sf.Clone('rsfof_' + histo_sf.GetName())
    ratio.Divide(histo_of)
    ratio.GetYaxis().SetTitle('r_{SFOF}')

    doFit = 0
    if doFit:
        fit = TF1('myfit','pol0', ratio.GetXaxis().GetXmin(), ratio.GetXaxis().GetXmax())
        fit.SetLineColor(r.kBlack if dataMC  == 'DATA' else r.kRed+1)
        fit.SetLineStyle(2)
        
        ratio.Fit('myfit','0')

    ratio.GetYaxis().SetRangeUser(0,3)

    return ratio


def runAnalysis(lumi, treeDA, treeMC, cuts, specialcut, tag, save, ingredientsFile):

    labelx = "m_{ll} [GeV]"
    labelmet = "E_{T}^{miss} [GeV]"
    labelnjet = "N. Jets"
    labelmt2 = "mt2 [GeV]"

    #####Main mll control plot
    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Starting mll plot' + bcolors.ENDC
    bins = [20,60,70,110,150,200,300,400] 
    MCControlSF =      treeMC.getTH1F(lumi, "MCControlSF", "lepsMll_Edge", bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectControlRegion, cuts.SF]), '', labelx)
    MCControlee =      treeMC.getTH1F(lumi, "MCControlee", "lepsMll_Edge", bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectControlRegion, cuts.ee]), '', labelx)
    MCControlmm =      treeMC.getTH1F(lumi, "MCControlmm", "lepsMll_Edge", bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectControlRegion, cuts.mm]), '', labelx)
    MCControlOF =      treeMC.getTH1F(lumi, "MCControlOF", "lepsMll_Edge", bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectControlRegion, cuts.OF]), '', labelx)
    MCSignalSF =       treeMC.getTH1F(lumi, "MCSignalSF", "lepsMll_Edge", bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectSignalRegion, cuts.SF]), '', labelx)
    MCSignalOF =       treeMC.getTH1F(lumi, "MCSignalOF", "lepsMll_Edge", bins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectSignalRegion, cuts.OF]), '', labelx)
    MCControlSFvalue = treeMC.getTH1F(lumi, "MCControlSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectControlRegion, cuts.SF]), '', labelx)
    MCControlOFvalue = treeMC.getTH1F(lumi, "MCControlOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectControlRegion, cuts.OF]), '', labelx)
    MCSignalSFvalue =  treeMC.getTH1F(lumi, "MCSignalSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectSignalRegion, cuts.SF]), '', labelx)
    MCSignalOFvalue =  treeMC.getTH1F(lumi, "MCSignalOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.RSFOFDirectSignalRegion, cuts.OF]), '', labelx)
    DataControlee =    treeDA.getTH1F(lumi, "DataControlee", "lepsMll_Edge", bins, 1, 1, cuts.AddList([specialcut, cuts.trigger, cuts.goodLepton,cuts.RSFOFDirectControlRegion, cuts.ee]), '', labelx)
    DataControlmm=treeDA.getTH1F(lumi,"DataControlmm","lepsMll_Edge",bins,1,1,cuts.AddList([specialcut,cuts.trigger,cuts.goodLepton,cuts.RSFOFDirectControlRegion, cuts.mm]), '', labelx)
    DataControlSF=treeDA.getTH1F(lumi,"DataControlSF","lepsMll_Edge",bins,1,1,cuts.AddList([specialcut,cuts.trigger,cuts.goodLepton,cuts.RSFOFDirectControlRegion, cuts.SF]), '', labelx)
    DataControlOF=treeDA.getTH1F(lumi,"DataControlOF","lepsMll_Edge",bins,1,1,cuts.AddList([specialcut,cuts.trigger,cuts.goodLepton,cuts.RSFOFDirectControlRegion, cuts.OF]), '', labelx)
    DataControlSFvalue= treeDA.getTH1F(lumi,"DataControlSFvalue","lepsMll_Edge",[20,1000], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.trigger, cuts.RSFOFDirectControlRegion, cuts.SF]), '', labelx)
    DataControlOFvalue= treeDA.getTH1F(lumi,"DataControlOFvalue","lepsMll_Edge",[20,1000], 1, 1, cuts.AddList([specialcut, cuts.goodLepton, cuts.trigger, cuts.RSFOFDirectControlRegion, cuts.OF]), '', labelx)
    MCControl =          make_rsfof(MCControlSF, MCControlOF, "MC")
    MCSignal =           make_rsfof(MCSignalSF, MCSignalOF, "MC")
    MCControlvalue =     make_rsfof(MCControlSFvalue, MCControlOFvalue, "MC")
    MCSignalvalue =      make_rsfof(MCSignalSFvalue, MCSignalOFvalue, "MC")
    DataControl =        make_rsfof(DataControlSF, DataControlOF, "DATA")
    DataControlvalue =   make_rsfof(DataControlSFvalue, DataControlOFvalue, "DATA")

    measuredValueMC =         MCControlvalue.GetBinContent(1)
    measuredValueUncMC =      MCControlvalue.GetBinError(1)
    measuredValueData =       DataControlvalue.GetBinContent(1)
    if(measuredValueUncMC < measuredValueData*0.04):
        measuredValueUncMC = measuredValueData*0.04
    measuredValueUncData =    DataControlvalue.GetBinError(1)
    measuredValueUncTotData = math.sqrt(measuredValueUncData**2 + measuredValueUncMC**2)
    
    print "MC:   " , measuredValueMC
    print "Data: " ,measuredValueData
    makeTable(DataControlee, DataControlmm, DataControlSF, DataControlOF, MCControlee, MCControlmm, MCControlSF, MCControlOF,measuredValueData, measuredValueUncData, measuredValueMC, measuredValueUncMC)
    plot_rsfof = Canvas.Canvas('rsfof/%s_%s/plot_rsfof_mll'%(lumi_str,tag), 'png,pdf',  0.45, 0.6, 0.7, 0.8)
    plot_rsfof.addHisto(MCControl, 'PE', 'Control region - MC', 'PL', r.kRed+1 , 1, 0)
    #plot_rsfof.addHisto(MCControlvalue, 'PE,SAME', 'ttjets measurement - MC', 'PL', r.kBlue , 1, 0)
    plot_rsfof.addHisto(MCSignal, 'PE,SAME', 'Signal region - MC', 'PL', r.kGreen+1 , 1, 0)
    plot_rsfof.addHisto(DataControl, 'PE,SAME', 'Control region - DATA', 'PL', r.kBlack , 1, 0)
    #plot_rsfof.addHisto(DataControlvalue, 'PE,SAME', 'ttjets measurement - DATA', 'PL', r.kGreen , 1, 0)
    plot_rsfof.addBand (MCControl.GetXaxis().GetXmin(), measuredValueData-measuredValueUncTotData, MCControl.GetXaxis().GetXmax(), measuredValueData+measuredValueUncTotData, r.kGreen, 0.2)
    plot_rsfof.addLine (MCControl.GetXaxis().GetXmin(), measuredValueData, MCControl.GetXaxis().GetXmax(), measuredValueData, r.kBlack)
    plot_rsfof.save(1, 1, 0, lumi, 0.5, 1.8)
    
    if save==True:  
        "saved!!!"
        saveInFile(ingredientsFile, measuredValueMC, measuredValueUncMC, measuredValueData, measuredValueUncData)    

    ######Met and JET plots 
    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Starting jet and met plot' + bcolors.ENDC
    MCSignalSFMET =treeMC.getTH1F(lumi,"MCSignalSFMET", "met_Edge", [50, 100, 150, 200, 300], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFDirectSignalRegionNoMET, cuts.SF]), '', labelmet)
    DASignalSFMET =treeDA.getTH1F(lumi,"DASignalSFMET", "met_Edge", [50, 100, 150, 200, 300], 1, 1, cuts.AddList([cuts.trigger, cuts.goodLepton, cuts.RSFOFDirectSignalRegionNoMET, cuts.SF]), '', labelmet)
    MCSignalOFMET =treeMC.getTH1F(lumi,"MCSignalOFMET", "met_Edge", [50, 100, 150, 200, 300], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFDirectSignalRegionNoMET, cuts.OF]), '', labelmet)
    DASignalOFMET =treeDA.getTH1F(lumi,"DASignalOFMET", "met_Edge", [50, 100, 150, 200, 300], 1, 1, cuts.AddList([cuts.trigger, cuts.goodLepton, cuts.RSFOFDirectSignalRegionNoMET, cuts.OF]), '', labelmet)
    MCSignalSFJet =treeMC.getTH1F(lumi,"MCSignalSFJets", "nJetSel_Edge", 8, 0, 8, cuts.AddList([cuts.goodLepton, cuts.RSFOFDirectSignalRegionNoJet, cuts.SF]), '', labelnjet)
    DASignalSFJet =treeDA.getTH1F(lumi,"DASignalSFJets", "nJetSel_Edge", 8, 0, 8, cuts.AddList([cuts.trigger, cuts.goodLepton, cuts.RSFOFDirectSignalRegionNoJet, cuts.SF]), '', labelnjet)
    MCSignalOFJet =treeMC.getTH1F(lumi,"MCSignalOFJets", "nJetSel_Edge", 8, 0, 8, cuts.AddList([cuts.goodLepton, cuts.RSFOFDirectSignalRegionNoJet, cuts.OF]), '', labelnjet)
    DASignalOFJet =treeDA.getTH1F(lumi,"DASignalOFJets", "nJetSel_Edge", 8, 0, 8, cuts.AddList([cuts.trigger, cuts.goodLepton, cuts.RSFOFDirectSignalRegionNoJet, cuts.OF]), '', labelnjet)
    MCSignalSFmt2 =treeMC.getTH1F(lumi,"MCSignalSFmt2","mt2_Edge", [0, 20, 40, 60, 80, 100,  130], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFDirectSignalRegionNoJet, cuts.SF]), '', labelmt2)
    DASignalSFmt2 =treeDA.getTH1F(lumi,"DASignalSFmt2","mt2_Edge", [0, 20, 40, 60, 80, 100,  130], 1, 1, cuts.AddList([cuts.trigger, cuts.goodLepton, cuts.RSFOFDirectSignalRegionNoJet, cuts.SF]), '', labelmt2)
    MCSignalOFmt2 =treeMC.getTH1F(lumi,"MCSignalOFmt2","mt2_Edge", [0, 20, 40, 60, 80, 100,  130], 1, 1, cuts.AddList([cuts.goodLepton, cuts.RSFOFDirectSignalRegionNoJet, cuts.OF]), '', labelmt2)
    DASignalOFmt2 =treeDA.getTH1F(lumi,"DASignalOFmt2","mt2_Edge", [0, 20, 40, 60, 80, 100,  130], 1, 1, cuts.AddList([cuts.trigger, cuts.goodLepton, cuts.RSFOFDirectSignalRegionNoJet, cuts.OF]), '', labelmt2)

    MCSignalMET   =make_rsfof(MCSignalSFMET, MCSignalOFMET, "MC")
    DASignalMET   =make_rsfof(DASignalSFMET, DASignalOFMET, "DATA")
    MCSignalJet   =make_rsfof(MCSignalSFJet, MCSignalOFJet, "MC")
    DASignalJet   =make_rsfof(DASignalSFJet, DASignalOFJet, "DATA")
    MCSignalmt2   =make_rsfof(MCSignalSFmt2, MCSignalOFmt2, "MC")
    DASignalmt2   =make_rsfof(DASignalSFmt2, DASignalOFmt2, "DATA")

    plot_rsfofmet = Canvas.Canvas('rsfof/%s_%s/plot_rsfof_met'%(lumi_str,tag), 'png,pdf', 0.45, 0.6, 0.7, 0.8)
    plot_rsfofmet.addHisto(MCSignalMET, 'PE', 'ttjets region - MC', 'PL', r.kRed+1 , 1, 0)
    plot_rsfofmet.addHisto(DASignalMET, 'PE, SAME', 'ttjets region - Data', 'PL', r.kBlack , 1, 0)
    plot_rsfofmet.addBand(MCSignalMET.GetXaxis().GetXmin(), measuredValueMC-measuredValueUncMC, MCSignalMET.GetXaxis().GetXmax(), measuredValueMC+measuredValueUncMC, r.kGreen, 0.2)
    plot_rsfofmet.addLine(MCSignalMET.GetXaxis().GetXmin(), measuredValueMC, MCSignalMET.GetXaxis().GetXmax(), measuredValueMC, r.kBlack)
    plot_rsfofmet.save(1, 1, 0, lumi, 0.5, 1.8)

    plot_rsfofjet = Canvas.Canvas('rsfof/%s_%s/plot_rsfof_jet'%(lumi_str,tag), 'png,pdf', 0.45, 0.6, 0.7, 0.8)
    plot_rsfofjet.addHisto(MCSignalJet, 'PE', 'ttjets region - MC', 'PL', r.kRed+1 , 1, 0)
    #plot_rsfofjet.addHisto(DASignalJet, 'PE, SAME', 'ttjets region - Data', 'PL', r.kBlack , 1, 0)
    plot_rsfofjet.addBand(MCSignalJet.GetXaxis().GetXmin(), measuredValueMC-measuredValueUncMC, MCSignalJet.GetXaxis().GetXmax(), measuredValueMC+measuredValueUncMC, r.kGreen, 0.2)
    plot_rsfofjet.addLine(MCSignalJet.GetXaxis().GetXmin(), measuredValueMC, MCSignalJet.GetXaxis().GetXmax(), measuredValueMC, r.kBlack)
    plot_rsfofjet.save(1, 1, 0, lumi, 0.5, 1.8)                                                                                                                                           
    
    plot_rsfofmt2 = Canvas.Canvas('rsfof/%s_%s/plot_rsfof_mt2'%(lumi_str,tag), 'png,pdf', 0.45, 0.6, 0.7, 0.8)
    plot_rsfofmt2.addHisto(MCSignalmt2, 'PE', 'ttjets region - MC', 'PL', r.kRed+1 , 1, 0)
    #plot_rsfofmt2.addHisto(DASignalmt2, 'PE, SAME', 'ttjets region - Data', 'PL', r.kBlack , 1, 0)
    plot_rsfofmt2.addBand(MCSignalmt2.GetXaxis().GetXmin(), measuredValueMC-measuredValueUncMC, MCSignalmt2.GetXaxis().GetXmax(), measuredValueMC+measuredValueUncMC, r.kGreen, 0.2)
    plot_rsfofmt2.addLine(MCSignalmt2.GetXaxis().GetXmin(), measuredValueMC, MCSignalmt2.GetXaxis().GetXmax(), measuredValueMC, r.kBlack)
    plot_rsfofmt2.save(1, 1, 0, lumi, 0.5, 1.8)                                                                                                                                           


##Main body of the analysis
if __name__ == '__main__':

    print bcolors.HEADER 
    print '#######################################################################'
    print '                  Starting r_SFOF analysis...                          ' 
    print '#######################################################################' + bcolors.ENDC

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientsFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    (opts, args) = parser.parse_args()

    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Loading DATA and MC trees...' + bcolors.ENDC

    DYDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO']
    ttDatasets = ['TTJets_DiLepton']
    mcDatasets = ['TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT',  'T_tch_powheg', 'TBar_tch_powheg', 'WWTo2L2Nu',  'WZTo3LNu','WZTo2L2Q', 'ZZTo4L', 'ZZTo2L2Nu', 'ZZTo2L2Q', 'WWW', 'WWZ', 'WZZ', 'ZZZ', 'TWZ', 'tZq_ll', 'TTZToLLNuNu' , 'TTZToQQ', 'TTLLJets_m1to10', 'TTWToLNu','TTWToQQ',  'TTTT', 'TTHnobb_pow', 'VHToNonbb',  'GGHZZ4L',  'WJetsToLNu_LO']
    mcDatasets += ttDatasets
    mcDatasets += DYDatasets

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
    
    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)
  
    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    maxrun = 999999
    lumi = 35.9 ; maxrun = 999999; lumi_str = '35.9'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 

    cuts = CutManager.CutManager()
    
    runAnalysis(lumi, treeDA, treeMC, cuts, '', 'nocut', True, opts.ingredientsFile)
 
