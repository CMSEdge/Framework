############################################################
############################################################
##          _          _            _              _      ##
##         /\ \       /\ \         /\ \           /\ \    ##
##        /  \ \     /  \ \____   /  \ \         /  \ \   ##
##       / /\ \ \   / /\ \_____\ / /\ \_\       / /\ \ \  ##
##      / / /\ \_\ / / /\/___  // / /\/_/      / / /\ \_\ ##
##     / /_/_ \/_// / /   / / // / / ______   / /_/_ \/_/ ##
##    / /____/\  / / /   / / // / / /\_____\ / /____/\    ##
##   / /\____\/ / / /   / / // / /  \/____ // /\____\/    ##
##  / / /______ \ \ \__/ / // / /_____/ / // / /______    ##
## / / /_______\ \ \___\/ // / /______\/ // / /_______\   ##
## \/__________/  \/_____/ \/___________/ \/__________/   ##
############################################################
############################################################
                                                       

import ROOT as r
from ROOT import gROOT, TCanvas, TFile, TF1, TGraphAsymmErrors
import math,sys,optparse, copy, subprocess
import math



import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Rounder    as rounder


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\032.4m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'




def makeTable(DATAnumeratoree, DATAnumeratormm, DATAnumeratorOF, DATAdenominatoree, DATAdenominatormm, DATAdenominatorOF, effMlleevalue, effMllmmvalue , effMllOFvalue, eevale, mmvale, emvale):
    line0 = '  \hline'
    line1 = '  ee &' 
    line2 = '  \mu\mu   &' 
    line3 = '  e\mu   &' 
    line0 += ' & nominator  & denominator & $\epsilon_{\mathrm{trigger}}$ $\pm \sigma_{stat}}$ \\\\'
    line1 += ' %.f & %.f &    %.3f $\\pm$ %.3f       %s' %(DATAnumeratoree.Integral(), DATAdenominatoree.Integral(), effMlleevalue  , eevale  ,   '\\\\')
    line2 += ' %.f & %.f &    %.3f $\\pm$ %.3f       %s' %(DATAnumeratormm.Integral(), DATAdenominatormm.Integral(), effMllmmvalue  , mmvale  ,   '\\\\')
    line3 += ' %.f & %.f &    %.3f $\\pm$ %.3f       %s' %(DATAnumeratorOF.Integral(), DATAdenominatorOF.Integral(), effMllOFvalue  , emvale  ,   '\\\\')
    line0 += '\\hline'; line3 += '\\hline';
                                                                                                                                                                                     
    helper.ensureDirectory('plots/rt/%s/'%lumi_str)
    helper.ensureDirectory('plots/rt/%s/tables/'%lumi_str)
    compTableFile = open('plots/rt/%s/tables/resultTable_%s%s.txt'%(lumi_str, str(lumi).replace('.','p'), "rt"),'w')
    compTableFile.write(line0+'\n')
    compTableFile.write(line1+'\n')
    compTableFile.write(line2+'\n')                                                                                             
    compTableFile.write(line3+'\n')                                                                                             
    print line0
    print line1
    print line2                                                                                                                                                                      
    print line3                                                                                                                                                                      


############################################################
def getRT(eff_ee, unc_ee, eff_mm, unc_mm, eff_em, unc_em):
    if(eff_em == 0):
        print "Division by zero"
        return [0, 0, 0]
    RT = math.sqrt(eff_ee*eff_mm)/eff_em
    uncRTee = (math.sqrt(eff_mm)/eff_em)*0.5*(unc_ee/math.sqrt(eff_ee)) 
    uncRTmm = (math.sqrt(eff_ee)/eff_em)*0.5*(unc_mm/math.sqrt(eff_mm)) 
    uncRTem = math.sqrt(eff_ee*eff_mm)*((unc_em)/(eff_em*eff_em)) 
    uncRT = math.sqrt(uncRTee * uncRTee + uncRTmm * uncRTmm + uncRTem * uncRTem)
    uncsysRTee = (math.sqrt(eff_mm)/eff_em)*0.5*(0.05*eff_ee/math.sqrt(eff_ee)) 
    uncsysRTmm = (math.sqrt(eff_ee)/eff_em)*0.5*(0.05*eff_mm/math.sqrt(eff_mm)) 
    uncsysRTem = math.sqrt(eff_ee*eff_mm)*((0.05*eff_em)/(eff_em*eff_em)) 
    uncsysRT = math.sqrt(uncsysRTee * uncsysRTee + uncsysRTmm * uncsysRTmm + uncsysRTem * uncsysRTem)

    return [RT, uncRT, uncsysRT]


def RT(ref, eff_ee, eff_mm, eff_em):

    RTratio = copy.deepcopy(ref)
    RTratio.SetName("RT_" + ref.GetName())
    RTratioee = copy.deepcopy(ref)
    RTratioee.SetName("RT_ee" + ref.GetName())
    RTratiomm = copy.deepcopy(ref)
    RTratiomm.SetName("RT_mm" + ref.GetName())
    RTratioem = copy.deepcopy(ref)
    RTratioem.SetName("RT_em" + ref.GetName())

    RTratio.GetYaxis().SetTitle("RT")
    RTratio.GetXaxis().SetTitle(ref.GetXaxis().GetTitle())
    RTratioee.GetYaxis().SetTitle("RT")
    RTratioee.GetXaxis().SetTitle(ref.GetXaxis().GetTitle())
    RTratiomm.GetYaxis().SetTitle("RT")
    RTratiomm.GetXaxis().SetTitle(ref.GetXaxis().GetTitle())
    RTratioem.GetYaxis().SetTitle("RT")
    RTratioem.GetXaxis().SetTitle(ref.GetXaxis().GetTitle())
    RTratioee.Reset()
    RTratiomm.Reset()
    RTratioem.Reset()

    xee = eff_ee.GetX()
    yee = eff_ee.GetY()
    eyhee = eff_ee.GetEYhigh()
    eylee = eff_ee.GetEYlow()
    for i in range(0, eff_ee.GetN()):
        j = RTratioee.FindBin(xee[i])
        RTratioee.SetBinContent(j, yee[i])
        RTratioee.SetBinError(j, max(eyhee[i], eylee[i]))

    xmm = eff_mm.GetX()
    ymm = eff_mm.GetY()
    eyhmm = eff_mm.GetEYhigh()
    eylmm = eff_mm.GetEYlow()
    for i in range(0, eff_mm.GetN()):
        j = RTratiomm.FindBin(xmm[i])
        RTratiomm.SetBinContent(j, ymm[i])
        RTratiomm.SetBinError(j, max(eyhmm[i], eylmm[i]))

    xem = eff_em.GetX()
    yem = eff_em.GetY()
    eyhem = eff_em.GetEYhigh()
    eylem = eff_em.GetEYlow()
    for i in range(0, eff_em.GetN()):
        j = RTratioem.FindBin(xem[i])
        RTratioem.SetBinContent(j, yem[i])
        RTratioem.SetBinError(j, max(eyhem[i], eylem[i]))

    for i in range(1, RTratioee.GetNbinsX()+2):
        [rt, uncrt, uncsys] = getRT(RTratioee.GetBinContent(i), RTratioee.GetBinError(i), RTratiomm.GetBinContent(i), RTratiomm.GetBinError(i), RTratioem.GetBinContent(i), RTratioem.GetBinError(i)) 
        RTratio.SetBinContent(i, rt)
        RTratio.SetBinError(i, uncrt)
    
    #SetOwnership(RTratio, 0 )   # 0 = release (not keep), 1 = keep
    #SetOwnership(RTratioee, 0 )   # 0 = release (not keep), 1 = keep
    #SetOwnership(RTratiomm, 0 )   # 0 = release (not keep), 1 = keep
    #SetOwnership(RTratioem, 0 )   # 0 = release (not keep), 1 = keep
    #del RTratioee
    #del RTratiomm
    #del RTratioem
    return RTratio


def getTriggerEffs(num, den):

    trig = den.Clone("trig_" + den.GetName())
    trig.Reset()



    errs = TGraphAsymmErrors(num, den, 'v')
    #x = errs.GetX()
    #y = errs.GetY()
    #eyh = errs.GetEYhigh()
    #eyl = errs.GetEYlow()
    #for i in range(1, num.GetNbinsX()+1):
    #    if(x[i] != 0 or y[i] != 0):
    #        trig.SetBinContent(i, y[i])
    #        trig.SetBinError(i, max(eyh[i], eyl[i]))

    return errs


def saveInFile(theFile, measuredValueMC, measuredValueUncMC, measuredValueSystMC, measuredValueData, measuredValueUncData, measuredValueSystData):


    foutput = open(theFile + "_aux", 'w')
    for line in open(theFile).readlines():
        if line.find("rt") != -1:
            if line.find("DATA") != -1:
                foutput.write('rt          region          DATA        %.4f      %0.4f       %.4f\n'%(measuredValueData, measuredValueUncData, measuredValueSystData))
            else:
                foutput.write('rt          region          MC          %.4f      %0.4f       %.4f\n'%(measuredValueMC, measuredValueUncMC, measuredValueSystMC))
        else:
            foutput.write(line)

    foutput.close()
    subprocess.call(['mv ' + theFile + "_aux " + theFile], shell=True)



##Main body of the analysis
if __name__ == '__main__':

    print bcolors.HEADER
    print '#######################################################################'
    print '                  Starting RT analysis...                          '
    print '#######################################################################' + bcolors.ENDC

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientsFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    (opts, args) = parser.parse_args()

    theFile = opts.ingredientsFile

    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Loading DATA and MC trees...' + bcolors.ENDC

    #mcDatasets = ['TT_pow_ext34', 'DYJetsToLL_M10to50', 'DYJetsToLL_M50']
    #daDatasets = ['JetHT_Run2016B-PromptReco-v2_runs_271036_276097', 'HTMHT_Run2016B-PromptReco-v2_runs_271036_276097',
    #        'JetHT_Run2016C-PromptReco-v2_runs_271036_276811', 'HTMHT_Run2016C-PromptReco-v2_runs_271036_276811',
    #        'JetHT_Run2016D-PromptReco-v2_runs_271036_276811', 'HTMHT_Run2016D-PromptReco-v2_runs_271036_276811']
    daDatasets = ['JetHT_Run2016B-PromptReco-v2_runs_271036_276097',
            'JetHT_Run2016C-PromptReco-v2_runs_271036_276811',
            'JetHT_Run2016D-PromptReco-v2_runs_271036_276811']

    #treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)

    print bcolors.HEADER + '[RSFOFAnalysis] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    
    lumi = 12.9 ; maxrun = 999999; lumi_str = '12.9invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    cuts = CutManager.CutManager()
    labelx = 'm_{ll} [GeV]'
    labelmet = 'E_{T}^{miss} [GeV]'
    labelmt2 = 'mT2 [GeV]'
    labelpt1 = 'p_{T} leading [GeV]'
    labelpt2 = 'p_{T} subleading [GeV]'
    labeleta1 = '#eta leading'
    labeleta2 = '#eta subleading'
    mllbins = [20, 70, 81, 101, 111, 300]
    metbins = [0, 50, 100, 150]
    mt2bins = [0, 20, 40, 60, 80, 100, 160]
    ptbins = [ 25, 30, 40, 50, 70, 100, 150, 200]
    etabins = [-3.0, -2.4, -1.8, -1.2, -0.6, 0, 0.6, 1.2, 1.8, 2.4, 3.0]
    specialcut = ''

    DATAdenominatorMllee = treeDA.getTH1F(lumi, "DATAdenominatoree", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelx)
    DATAnumeratorMllee =   treeDA.getTH1F(lumi, "DATAnumeratoree", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labelx)
    DATAdenominatorMllmm = treeDA.getTH1F(lumi, "DATAdenominatormm", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelx)
    DATAnumeratorMllmm =   treeDA.getTH1F(lumi, "DATAnumeratormm", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labelx)
    DATAdenominatorMllSF = treeDA.getTH1F(lumi, "DATAdenominatorSF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelx)
    DATAnumeratorMllSF =   treeDA.getTH1F(lumi, "DATAnumeratorSF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labelx)
    DATAdenominatorMllOF = treeDA.getTH1F(lumi, "DATAdenominatorOF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelx)
    DATAnumeratorMllOF =   treeDA.getTH1F(lumi, "DATAnumeratorOF", "lepsMll_Edge", mllbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labelx)
    
    DATAdenominatorpt1ee = treeDA.getTH1F(lumi, "DATAdenominatorpt1ee", "Lep1_pt_Edge", ptbins , 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelpt1)
    DATAnumeratorpt1ee =   treeDA.getTH1F(lumi, "DATAnumeratorpt1ee", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labelpt1)
    DATAdenominatorpt1mm = treeDA.getTH1F(lumi, "DATAdenominatorpt1mm", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelpt1)
    DATAnumeratorpt1mm =   treeDA.getTH1F(lumi, "DATAnumeratorpt1mm", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labelpt1)
    DATAdenominatorpt1SF = treeDA.getTH1F(lumi, "DATAdenominatorpt1SF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelpt1)
    DATAnumeratorpt1SF =   treeDA.getTH1F(lumi, "DATAnumeratorpt1SF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labelpt1)
    DATAdenominatorpt1OF = treeDA.getTH1F(lumi, "DATAdenominatorpt1OF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelpt1)
    DATAnumeratorpt1OF =   treeDA.getTH1F(lumi, "DATAnumeratorpt1OF", "Lep1_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labelpt1)
 
    DATAdenominatorpt2ee = treeDA.getTH1F(lumi, "DATAdenominatorpt2ee", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelpt2)
    DATAnumeratorpt2ee =   treeDA.getTH1F(lumi, "DATAnumeratorpt2ee", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labelpt2)
    DATAdenominatorpt2mm = treeDA.getTH1F(lumi, "DATAdenominatorpt2mm", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelpt2)
    DATAnumeratorpt2mm =   treeDA.getTH1F(lumi, "DATAnumeratorpt2mm", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labelpt2)
    DATAdenominatorpt2SF = treeDA.getTH1F(lumi, "DATAdenominatorpt2SF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelpt2)
    DATAnumeratorpt2SF =   treeDA.getTH1F(lumi, "DATAnumeratorpt2SF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labelpt2)
    DATAdenominatorpt2OF = treeDA.getTH1F(lumi, "DATAdenominatorpt2OF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelpt2)
    DATAnumeratorpt2OF =   treeDA.getTH1F(lumi, "DATAnumeratorpt2OF", "Lep2_pt_Edge", ptbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labelpt2)
 
    DATAdenominatoreta1ee = treeDA.getTH1F(lumi, "DATAdenominatoreta1ee", "Lep1_eta_Edge", etabins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labeleta1)
    DATAnumeratoreta1ee =   treeDA.getTH1F(lumi, "DATAnumeratoreta1ee", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labeleta1)
    DATAdenominatoreta1mm = treeDA.getTH1F(lumi, "DATAdenominatoreta1mm", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labeleta1)
    DATAnumeratoreta1mm =   treeDA.getTH1F(lumi, "DATAnumeratoreta1mm", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labeleta1)
    DATAdenominatoreta1SF = treeDA.getTH1F(lumi, "DATAdenominatoreta1SF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labeleta1)
    DATAnumeratoreta1SF =   treeDA.getTH1F(lumi, "DATAnumeratoreta1SF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labeleta1)
    DATAdenominatoreta1OF = treeDA.getTH1F(lumi, "DATAdenominatoreta1OF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labeleta1)
    DATAnumeratoreta1OF =   treeDA.getTH1F(lumi, "DATAnumeratoreta1OF", "Lep1_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labeleta1)
 
    DATAdenominatoreta2ee = treeDA.getTH1F(lumi, "DATAdenominatoreta2ee", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labeleta2)
    DATAnumeratoreta2ee =   treeDA.getTH1F(lumi, "DATAnumeratoreta2ee", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labeleta2)
    DATAdenominatoreta2mm = treeDA.getTH1F(lumi, "DATAdenominatoreta2mm", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labeleta2)
    DATAnumeratoreta2mm =   treeDA.getTH1F(lumi, "DATAnumeratoreta2mm", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labeleta2)
    DATAdenominatoreta2SF = treeDA.getTH1F(lumi, "DATAdenominatoreta2SF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labeleta2)
    DATAnumeratoreta2SF =   treeDA.getTH1F(lumi, "DATAnumeratoreta2SF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labeleta2)
    DATAdenominatoreta2OF = treeDA.getTH1F(lumi, "DATAdenominatoreta2OF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labeleta2)
    DATAnumeratoreta2OF =   treeDA.getTH1F(lumi, "DATAnumeratoreta2OF", "Lep2_eta_Edge", etabins, 1, 1,cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labeleta2)
 
    DATAdenominatorMlleevalue =   treeDA.getTH1F(lumi, "DATAdenominatoreevalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelx)
    DATAnumeratorMlleevalue =     treeDA.getTH1F(lumi, "DATAnumeratoreevalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labelx)
    DATAdenominatorMllmmvalue =   treeDA.getTH1F(lumi, "DATAdenominatormmvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelx)
    DATAnumeratorMllmmvalue =     treeDA.getTH1F(lumi, "DATAnumeratormmvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labelx)
    DATAdenominatorMllSFvalue =   treeDA.getTH1F(lumi, "DATAdenominatorSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelx)
    DATAnumeratorMllSFvalue =     treeDA.getTH1F(lumi, "DATAnumeratorSFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labelx)
    DATAdenominatorMllOFvalue =   treeDA.getTH1F(lumi, "DATAdenominatorOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelx)
    DATAnumeratorMllOFvalue =     treeDA.getTH1F(lumi, "DATAnumeratorOFvalue", "lepsMll_Edge", [20, 1000], 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labelx)
    
    DATAdenominatorMETee =   treeDA.getTH1F(lumi, "DATAdenominatoreevalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelmet)
    DATAnumeratorMETee =     treeDA.getTH1F(lumi, "DATAnumeratoreevalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labelmet)
    DATAdenominatorMETmm =   treeDA.getTH1F(lumi, "DATAdenominatormmvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelmet)
    DATAnumeratorMETmm =     treeDA.getTH1F(lumi, "DATAnumeratormmvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labelmet)
    DATAdenominatorMETSF =   treeDA.getTH1F(lumi, "DATAdenominatorSFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelmet)
    DATAnumeratorMETSF =     treeDA.getTH1F(lumi, "DATAnumeratorSFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labelmet)
    DATAdenominatorMETOF =   treeDA.getTH1F(lumi, "DATAdenominatorOFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelmet)
    DATAnumeratorMETOF =     treeDA.getTH1F(lumi, "DATAnumeratorOFvalue", "met_Edge", metbins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labelmet)

    DATAdenominatormt2ee =   treeDA.getTH1F(lumi, "DATAdenominatoreevalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.ee]), '', labelmt2)
    DATAnumeratormt2ee =     treeDA.getTH1F(lumi, "DATAnumeratoreevalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.ee]), '', labelmt2)
    DATAdenominatormt2mm =   treeDA.getTH1F(lumi, "DATAdenominatormmvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.mm]), '', labelmt2)
    DATAnumeratormt2mm =     treeDA.getTH1F(lumi, "DATAnumeratormmvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.mm]), '', labelmt2)
    DATAdenominatormt2SF =   treeDA.getTH1F(lumi, "DATAdenominatorSFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.SF]), '', labelmt2)
    DATAnumeratormt2SF =     treeDA.getTH1F(lumi, "DATAnumeratorSFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.SF]), '', labelmt2)
    DATAdenominatormt2OF =   treeDA.getTH1F(lumi, "DATAdenominatorOFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.denominator, cuts.OF]), '', labelmt2)
    DATAnumeratormt2OF =     treeDA.getTH1F(lumi, "DATAnumeratorOFvalue", "mt2_Edge", mt2bins, 1, 1, cuts.AddList([specialcut, cuts.numerator, cuts.trigger, cuts.OF]), '', labelmt2)

    
    ######################## Calculation of total values #############################
    effMlleevalue = getTriggerEffs(DATAnumeratorMlleevalue, DATAdenominatorMlleevalue)
    effMllmmvalue = getTriggerEffs(DATAnumeratorMllmmvalue, DATAdenominatorMllmmvalue)
    effMllSFvalue = getTriggerEffs(DATAnumeratorMllSFvalue, DATAdenominatorMllSFvalue)
    effMllOFvalue = getTriggerEffs(DATAnumeratorMllOFvalue, DATAdenominatorMllOFvalue)

    eeval_ = effMlleevalue.GetY()
    eevaleh_ = effMlleevalue.GetEYhigh()
    eevalel_ = effMlleevalue.GetEYlow()
    eeval = eeval_[0]
    eevale = max(eevaleh_[0], eevalel_[0])

    mmval_ = effMllmmvalue.GetY()
    mmvaleh_ = effMllmmvalue.GetEYhigh()
    mmvalel_ = effMllmmvalue.GetEYlow()
    mmval = mmval_[0]
    mmvale = max(mmvaleh_[0], mmvalel_[0])

    emval_ = effMllOFvalue.GetY()
    emvaleh_ = effMllOFvalue.GetEYhigh()
    emvalel_ = effMllOFvalue.GetEYlow()
    emval = emval_[0]
    emvale = max(emvaleh_[0], emvalel_[0])

    [rt, uncrt, systrt] = getRT(eeval, eevale, mmval, mmvale, emval, emvale)



    makeTable(DATAnumeratorMlleevalue, DATAnumeratorMllmmvalue, DATAnumeratorMllOFvalue, DATAdenominatorMlleevalue, DATAdenominatorMllmmvalue, DATAdenominatorMllOFvalue, eeval, mmval, emval, eevale, mmvale, emvale)
    ######################## Calculation of total values #############################

    effMETee = getTriggerEffs(DATAnumeratorMETee, DATAdenominatorMETee)
    effMETmm = getTriggerEffs(DATAnumeratorMETmm, DATAdenominatorMETmm)    
    effMETSF = getTriggerEffs(DATAnumeratorMETSF, DATAdenominatorMETSF)    
    effMETOF = getTriggerEffs(DATAnumeratorMETOF, DATAdenominatorMETOF)
   
    effmt2ee = getTriggerEffs(DATAnumeratormt2ee, DATAdenominatormt2ee)
    effmt2mm = getTriggerEffs(DATAnumeratormt2mm, DATAdenominatormt2mm)
    effmt2SF = getTriggerEffs(DATAnumeratormt2SF, DATAdenominatormt2SF)
    effmt2OF = getTriggerEffs(DATAnumeratormt2OF, DATAdenominatormt2OF)

    effMllee = getTriggerEffs(DATAnumeratorMllee, DATAdenominatorMllee)
    effMllmm = getTriggerEffs(DATAnumeratorMllmm, DATAdenominatorMllmm)
    effMllSF = getTriggerEffs(DATAnumeratorMllSF, DATAdenominatorMllSF)
    effMllOF = getTriggerEffs(DATAnumeratorMllOF, DATAdenominatorMllOF)
    
    effpt1ee = getTriggerEffs(DATAnumeratorpt1ee, DATAdenominatorpt1ee)
    effpt1mm = getTriggerEffs(DATAnumeratorpt1mm, DATAdenominatorpt1mm)
    effpt1SF = getTriggerEffs(DATAnumeratorpt1SF, DATAdenominatorpt1SF)
    effpt1OF = getTriggerEffs(DATAnumeratorpt1OF, DATAdenominatorpt1OF)
    
    effpt2ee = getTriggerEffs(DATAnumeratorpt2ee, DATAdenominatorpt2ee)
    effpt2mm = getTriggerEffs(DATAnumeratorpt2mm, DATAdenominatorpt2mm)
    effpt2SF = getTriggerEffs(DATAnumeratorpt2SF, DATAdenominatorpt2SF)
    effpt2OF = getTriggerEffs(DATAnumeratorpt2OF, DATAdenominatorpt2OF)
    
    effeta1ee = getTriggerEffs(DATAnumeratoreta1ee, DATAdenominatoreta1ee)
    effeta1mm = getTriggerEffs(DATAnumeratoreta1mm, DATAdenominatoreta1mm)
    effeta1SF = getTriggerEffs(DATAnumeratoreta1SF, DATAdenominatoreta1SF)
    effeta1OF = getTriggerEffs(DATAnumeratoreta1OF, DATAdenominatoreta1OF)
    
    effeta2ee = getTriggerEffs(DATAnumeratoreta2ee, DATAdenominatoreta2ee)
    effeta2mm = getTriggerEffs(DATAnumeratoreta2mm, DATAdenominatoreta2mm)
    effeta2SF = getTriggerEffs(DATAnumeratoreta2SF, DATAdenominatoreta2SF)
    effeta2OF = getTriggerEffs(DATAnumeratoreta2OF, DATAdenominatoreta2OF)


    RTMET = RT(DATAnumeratorMETSF, effMETee, effMETmm, effMETOF)
    RTmt2 = RT(DATAnumeratormt2SF, effmt2ee, effmt2mm, effmt2OF)
    RTMll = RT(DATAnumeratorMllSF, effMllee, effMllmm, effMllOF)
    RTpt1 = RT(DATAnumeratorpt1SF, effpt1ee, effpt1mm, effpt1OF)
    RTpt2 = RT(DATAnumeratorpt2SF, effpt2ee, effpt2mm, effpt2OF)
    RTeta1 = RT(DATAnumeratoreta1SF, effeta1ee, effeta1mm, effeta1OF)
    RTeta2 = RT(DATAnumeratoreta2SF, effeta2ee, effeta2mm, effeta2OF)

    effRTMll = Canvas.Canvas('rt/%s/plot_rt_mll'%(lumi_str), 'png,pdf',  0.6, 0.3, 0.8, 0.5)
    h_auxrtMll = r.TH1F("h_auxRTMll", "", 1, 0, 250)
    h_auxrtMll.GetYaxis().SetRangeUser(0.8, 1.2)
    h_auxrtMll.GetXaxis().SetRangeUser(0, 250)
    h_auxrtMll.GetXaxis().SetTitle(labelx)
    effRTMll.addHisto(h_auxrtMll, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
    effRTMll.addHisto(RTMll, 'PE,SAME', 'R_{T}', 'PL', r.kBlack , 1, 0)
    effRTMll.addBand(h_auxrtMll.GetXaxis().GetXmin(), rt-systrt, h_auxrtMll.GetXaxis().GetXmax(), rt+systrt, r.kOrange+6, 0.2)
    effRTMll.addLine(h_auxrtMll.GetXaxis().GetXmin(), rt, h_auxrtMll.GetXaxis().GetXmax(), rt ,r.kBlue)
    effRTMll.addLatex (0.6, 0.25, 'Mean R_{T}: %.2f '%(rt))
    effRTMll.save(1, 1, 0, lumi, 0.8, 1.2)                                                            
    
    effRTMET = Canvas.Canvas('rt/%s/plot_rt_met'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
    h_auxrtMET = r.TH1F("h_auxRTMET", "", 1, 0, 200)
    h_auxrtMET.GetYaxis().SetRangeUser(0.8, 1.2)
    h_auxrtMET.GetXaxis().SetRangeUser(0, 150)
    h_auxrtMET.GetXaxis().SetTitle(labelmet)
    effRTMET.addHisto(h_auxrtMET, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
    effRTMET.addHisto(RTMET, 'PE,SAME', 'R_{T}', 'PL', r.kBlack , 1, 0)
    effRTMET.addBand(h_auxrtMET.GetXaxis().GetXmin(), rt-systrt, h_auxrtMET.GetXaxis().GetXmax(), rt+systrt, r.kOrange+6, 0.2)
    effRTMET.addLine(h_auxrtMET.GetXaxis().GetXmin(), rt, h_auxrtMET.GetXaxis().GetXmax(), rt ,r.kBlue)
    effRTMET.addLatex (0.6, 0.25, 'Mean R_{T}: %.2f '%(rt))
    effRTMET.save(1, 1, 0, lumi, 0.8, 1.2)                                                                                                  

    
    effRTmt2 = Canvas.Canvas('rt/%s/plot_rt_mt2'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
    h_auxrtmt2 = r.TH1F("h_auxRTMET", "", 1, 0, 160)
    h_auxrtmt2.GetYaxis().SetRangeUser(0.8, 1.2)
    h_auxrtmt2.GetXaxis().SetRangeUser(0, 160)
    h_auxrtmt2.GetXaxis().SetTitle(labelmt2)
    effRTmt2.addHisto(h_auxrtmt2, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
    effRTmt2.addHisto(RTmt2, 'PE,SAME', 'R_{T}', 'PL', r.kBlack , 1, 0)
    effRTmt2.addBand(h_auxrtmt2.GetXaxis().GetXmin(), rt-systrt, h_auxrtmt2.GetXaxis().GetXmax(), rt+systrt, r.kOrange+6, 0.2)
    effRTmt2.addLine(h_auxrtmt2.GetXaxis().GetXmin(), rt, h_auxrtmt2.GetXaxis().GetXmax(), rt ,r.kBlue)
    effRTmt2.addLatex (0.6, 0.25, 'Mean R_{T}: %.2f '%(rt))
    effRTmt2.save(1, 1, 0, lumi, 0.8, 1.2)                                                                                                 

    effRTpt1 = Canvas.Canvas('rt/%s/plot_rt_pt1'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
    h_auxrtpt1 = r.TH1F("h_auxRTpt1", "", 1, 20, 150)
    h_auxrtpt1.GetYaxis().SetRangeUser(0.8, 1.2)
    h_auxrtpt1.GetXaxis().SetRangeUser(20, 150)
    h_auxrtpt1.GetXaxis().SetTitle(labelpt1)
    effRTpt1.addHisto(h_auxrtpt1, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
    effRTpt1.addHisto(RTpt1, 'PE,SAME', 'RT', 'PL', r.kBlack , 1, 0)
    effRTpt1.addBand(h_auxrtpt1.GetXaxis().GetXmin(), rt-systrt, h_auxrtpt1.GetXaxis().GetXmax(), rt+systrt, r.kOrange+6, 0.2)
    effRTpt1.addLine(h_auxrtpt1.GetXaxis().GetXmin(), rt, h_auxrtpt1.GetXaxis().GetXmax(), rt ,r.kBlue)
    effRTpt1.addLatex (0.6, 0.25, 'Mean R_{T}: %.2f '%(rt))
    effRTpt1.save(1, 1, 0, lumi, 0.8, 1.2)

    effRTpt2 = Canvas.Canvas('rt/%s/plot_rt_pt2'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
    h_auxrtpt2 = r.TH1F("h_auxRTpt2", "", 1, 20, 150)
    h_auxrtpt2.GetYaxis().SetRangeUser(0.8, 1.2)
    h_auxrtpt2.GetXaxis().SetRangeUser(20, 150)
    h_auxrtpt2.GetXaxis().SetTitle(labelpt2)
    effRTpt2.addHisto(h_auxrtpt2, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
    effRTpt2.addHisto(RTpt2, 'PE,SAME', 'R_{T}', 'PL', r.kBlack , 1, 0)
    effRTpt2.addBand(h_auxrtpt2.GetXaxis().GetXmin(), rt-systrt, h_auxrtpt2.GetXaxis().GetXmax(), rt+systrt, r.kOrange+6, 0.2)
    effRTpt2.addLine(h_auxrtpt2.GetXaxis().GetXmin(), rt, h_auxrtpt2.GetXaxis().GetXmax(), rt ,r.kBlue)
    effRTpt2.addLatex (0.6, 0.25, 'Mean R_{T}: %.2f '%(rt))
    effRTpt2.save(1, 1, 0, lumi, 0.8, 1.2)

    effRTeta1 = Canvas.Canvas('rt/%s/plot_rt_eta1'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
    h_auxrteta1 = r.TH1F("h_auxRTeta1", "", 1, -2.4, 2.4)
    h_auxrteta1.GetYaxis().SetRangeUser(0.8, 1.2)
    h_auxrteta1.GetXaxis().SetRangeUser(-2.4, 2.4)
    h_auxrteta1.GetXaxis().SetTitle(labeleta1)
    effRTeta1.addHisto(h_auxrteta1, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
    effRTeta1.addHisto(RTeta1, 'PE,SAME', 'RT', 'PL', r.kBlack , 1, 0)
    effRTeta1.addBand(h_auxrteta1.GetXaxis().GetXmin(), rt-systrt, h_auxrteta1.GetXaxis().GetXmax(), rt+systrt, r.kOrange+6, 0.2)
    effRTeta1.addLine(h_auxrteta1.GetXaxis().GetXmin(), rt, h_auxrteta1.GetXaxis().GetXmax(), rt ,r.kBlue)
    effRTeta1.addLatex (0.6, 0.25, 'Mean R_{T}: %.2f '%(rt))
    effRTeta1.save(1, 1, 0, lumi, 0.8, 1.2)

    effRTeta2 = Canvas.Canvas('rt/%s/plot_rt_eta2'%(lumi_str), 'png,pdf', 0.6, 0.3, 0.8, 0.5)
    h_auxrteta2 = r.TH1F("h_auxRTeta2", "", 1, -2.4, 2.4)
    h_auxrteta2.GetYaxis().SetRangeUser(0.8, 1.2)
    h_auxrteta2.GetXaxis().SetRangeUser(-2.4, 2.4)
    h_auxrteta2.GetXaxis().SetTitle(labeleta2)
    effRTeta2.addHisto(h_auxrteta2, 'h', '', 'R_{T}', r.kRed+1, 1, 0)
    effRTeta2.addHisto(RTeta2, 'PE,SAME', 'RT', 'PL', r.kBlack , 1, 0)
    effRTeta2.addBand(h_auxrteta2.GetXaxis().GetXmin(), rt-systrt, h_auxrteta2.GetXaxis().GetXmax(), rt+systrt, r.kOrange+6, 0.2)
    effRTeta2.addLine(h_auxrteta2.GetXaxis().GetXmin(), rt, h_auxrteta2.GetXaxis().GetXmax(), rt ,r.kBlue)
    effRTeta2.addLatex (0.6, 0.25, 'Mean R_{T}: %.2f '%(rt))
    effRTeta2.save(1, 1, 0, lumi, 0.8, 1.2)
   
    effMll = Canvas.Canvas('rt/%s/plot_eff_mll'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
    h_auxMll = r.TH1F("h_auxMll", "", 1, 0, 250)
    h_auxMll.GetYaxis().SetRangeUser(0, 2)
    h_auxMll.GetXaxis().SetRangeUser(0, 250)
    h_auxMll.GetXaxis().SetTitle(labelx)
    effMll.addHisto(h_auxMll, 'h', '', '', r.kRed+1, 1, 0)
    effMll.addHisto(effMllee, 'PE,SAME', 'Double Electron', 'PL', r.kRed+1 , 1, 0)
    effMll.addHisto(effMllmm, 'PE,SAME', 'Double Muon', 'PL', r.kBlue+1 , 1, 0)
    effMll.addHisto(effMllSF, 'PE,SAME', 'Same Flavor', 'PL', r.kGreen+1 , 1, 0)
    effMll.addHisto(effMllOF, 'PE,SAME', 'Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
    effMll.save(1, 1, 0, lumi, 0.2, 1.8)                                                                    

    effmt2 = Canvas.Canvas('rt/%s/plot_eff_mt2'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
    h_auxmt2 = r.TH1F("h_auxmt2", "", 1, 0, 160)
    h_auxmt2.GetYaxis().SetRangeUser(0, 2)
    h_auxmt2.GetXaxis().SetRangeUser(0, 160)
    h_auxmt2.GetXaxis().SetTitle(labelx)
    effmt2.addHisto(h_auxmt2, 'h', '', '', r.kRed+1, 1, 0)
    effmt2.addHisto(effmt2ee, 'PE,SAME', 'Double Electron', 'PL', r.kRed+1 , 1, 0)
    effmt2.addHisto(effmt2mm, 'PE,SAME', 'Double Muon', 'PL', r.kBlue+1 , 1, 0)
    effmt2.addHisto(effmt2SF, 'PE,SAME', 'Same Flavor', 'PL', r.kGreen+1 , 1, 0)
    effmt2.addHisto(effmt2OF, 'PE,SAME', 'Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
    effmt2.save(1, 1, 0, lumi, 0.2, 1.8)                                                                    







    effpt1 = Canvas.Canvas('rt/%s/plot_eff_pt1'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
    h_auxpt1 = r.TH1F("h_auxpt1", "", 1, 0, 150)
    h_auxpt1.GetYaxis().SetRangeUser(0, 2)
    h_auxpt1.GetXaxis().SetRangeUser(0, 150)
    h_auxpt1.GetXaxis().SetTitle(labelpt1)
    effpt1.addHisto(h_auxpt1, 'h', '', '', r.kRed+1, 1, 0)
    effpt1.addHisto(effpt1ee, 'PE,SAME', 'Double Electron', 'PL', r.kRed+1 , 1, 0)
    effpt1.addHisto(effpt1mm, 'PE,SAME', 'Double Muon', 'PL', r.kBlue+1 , 1, 0)
    effpt1.addHisto(effpt1SF, 'PE,SAME', 'Same Flavor', 'PL', r.kGreen+1 , 1, 0)
    effpt1.addHisto(effpt1OF, 'PE,SAME', 'Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
    effpt1.save(1, 1, 0, lumi, 0.2, 1.8)

    effpt2 = Canvas.Canvas('rt/%s/plot_eff_pt2'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
    h_auxpt2 = r.TH1F("h_auxpt2", "", 1, 0, 150)
    h_auxpt2.GetYaxis().SetRangeUser(0, 2)
    h_auxpt2.GetXaxis().SetRangeUser(0, 150)
    h_auxpt2.GetXaxis().SetTitle(labelpt2)
    effpt2.addHisto(h_auxpt2, 'h', '', '', r.kRed+1, 1, 0)
    effpt2.addHisto(effpt2ee, 'PE,SAME', 'Double Electron', 'PL', r.kRed+1 , 1, 0)
    effpt2.addHisto(effpt2mm, 'PE,SAME', 'Double Muon', 'PL', r.kBlue+1 , 1, 0)
    effpt2.addHisto(effpt2SF, 'PE,SAME', 'Same Flavor', 'PL', r.kGreen+1 , 1, 0)
    effpt2.addHisto(effpt2OF, 'PE,SAME', 'Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
    effpt2.save(1, 1, 0, lumi, 0.2, 1.8)
 
    effeta1 = Canvas.Canvas('rt/%s/plot_eff_eta1'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
    h_auxeta1 = r.TH1F("h_auxeta1", "", 1, -2.4, 2.4)
    h_auxeta1.GetYaxis().SetRangeUser(0, 2)
    h_auxeta1.GetXaxis().SetRangeUser(-2.4, 2.4)
    h_auxeta1.GetXaxis().SetTitle(labeleta1)
    effeta1.addHisto(h_auxeta1, 'h', '', '', r.kRed+1, 1, 0)
    effeta1.addHisto(effeta1ee, 'PE,SAME', 'Double Electron', 'PL', r.kRed+1 , 1, 0)
    effeta1.addHisto(effeta1mm, 'PE,SAME', 'Double Muon', 'PL', r.kBlue+1 , 1, 0)
    effeta1.addHisto(effeta1SF, 'PE,SAME', 'Same Flavor', 'PL', r.kGreen+1 , 1, 0)
    effeta1.addHisto(effeta1OF, 'PE,SAME', 'Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
    effeta1.save(1, 1, 0, lumi, 0.2, 1.8)

    effeta2 = Canvas.Canvas('rt/%s/plot_eff_eta2'%(lumi_str), 'png,pdf', 0.4, 0.2, 0.65, 0.4)
    h_auxeta2 = r.TH1F("h_auxeta2", "", 1, -2.4, 2.4)
    h_auxeta2.GetYaxis().SetRangeUser(0, 2)
    h_auxeta2.GetXaxis().SetRangeUser(-2.4, 2.4)
    h_auxeta2.GetXaxis().SetTitle(labeleta2)
    effeta2.addHisto(h_auxeta2, 'h', '', '', r.kRed+1, 1, 0)
    effeta2.addHisto(effeta2ee, 'PE,SAME', 'Double Electron', 'PL', r.kRed+1 , 1, 0)
    effeta2.addHisto(effeta2mm, 'PE,SAME', 'Double Muon', 'PL', r.kBlue+1 , 1, 0)
    effeta2.addHisto(effeta2SF, 'PE,SAME', 'Same Flavor', 'PL', r.kGreen+1 , 1, 0)
    effeta2.addHisto(effeta2OF, 'PE,SAME', 'Opposite Flavor', 'PL', r.kBlack+1 , 1, 0)
    effeta2.save(1, 1, 0, lumi, 0.2, 1.8)
 
 
    print 'Measured RT value', rt, ' +/- ', uncrt, ' +/- ', systrt
    saveInFile(theFile, 0, 0, 0, rt, uncrt, systrt)

