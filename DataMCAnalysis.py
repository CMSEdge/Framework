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

    plot = Canvas.Canvas('DataMC/%s/plot_mmeeComparison_%s'%(lumistr,name), 'png,pdf,root', 0.65, 0.7, 0.85, 0.9)
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

def makePlot(lumi, lumi_str, treeDA, treeMC, var, name, nbin, xmin, xmax, theCut, cuts,labelx, logx, returnplot = False, scaleToData = False, normalized = False, cumulative = False, onlyMC = False, treeSI=None):
    print theCut
    theCut = theCut.replace(cuts.trigger17,'1')
    if treeSI:
        print cuts.AddList( [theCut.replace(cuts.filters17,'1'), cuts.mass_600_25])
        SI = treeSI.getTH1F(lumi,"hSI_%s"%(name), var, nbin, xmin, xmax, cuts.AddList( [theCut.replace(cuts.filters17,'1'), cuts.mass_600_25]), "", labelx,'1','ZZpt')
        print 'Signal yield', SI.Integral()

    MCS  = treeMC.getStack(lumi, "hMCS_%s"%(name), var, nbin, xmin, xmax, theCut, "", labelx, "1", 'ZZpt')
    if not onlyMC and treeDA:
        DATA = treeDA.getTH1F(lumi, "hDATA_%s"%(name), var, nbin, xmin, xmax, theCut, '', labelx, "1", 'ZZpt')
    else: DATA = None

    f = 0
    for _i,_h in enumerate(MCS.GetHists()):
        if not f: MC = copy.deepcopy(_h)
        else:MC.Add(_h, 1.)               
        f = 1
        print "sample ", _h.GetName(), " has integral ", _h.Integral(), 'and overflow', _h.GetBinContent(_h.GetNbinsX())
        


    maxValmc = MC.GetMaximum() #GetBinContent(MC.GetMaximumBin())
    maxValdata = DATA.GetMaximum() if DATA else 0#GetBinContent(DATA.GetMaximumBin())
    maxVal = max(maxValmc, maxValdata)
    if not logx:
        MCS.SetMaximum(1.3*maxVal)
        if scaleToData:
            newStack.SetMaximum(1.3*maxVal)
        if DATA:
            DATA.GetYaxis().SetRangeUser(0.1, 1.3*maxVal)
    else:
        MCS.SetMaximum(2.0*maxVal)
        if scaleToData:
            newStack.SetMaximum(2.0*maxVal)
        if DATA:
            DATA.GetYaxis().SetRangeUser(0.1, 2.0*maxVal)

   
    print name
    SetOwnership(MC, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(MCS, 0 )   # 0 = release (not keep), 1 = keep
    if DATA: SetOwnership(DATA, 0 )   # 0 = release (not keep), 1 = keep
    plot = Canvas.Canvas('DataMC/%s/plot_%s'%(lumi_str,name), 'png,pdf,root', 0.7, 0.55, 0.9, 0.9)
    if not normalized:
        plot.addStack(MCS if not scaleToData else newStack, "HIST", 1, 1)
        if DATA:
            plot.addHisto(DATA, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
        if treeSI:
            SI.SetLineWidth(2)
            plot.addHisto(SI, "L,SAME", "Signal (600,25)", "L", r.kRed, 1,1)
    else:
        for h in hists:
            plot.addHisto(h, 'hist,same' if not h.GetName() == 'data' else 'p,same', h.GetName(), 'PL', h.GetLineColor(), 1, 0)
    plot.saveRatio(1, 1, logx, lumi, DATA if DATA else MC, MC)

    if cumulative:
        cum_plot = Canvas.Canvas('DataMC/%s/plot_cumulative_%s'%(lumi_str,name), 'png,pdf,root', 0.65, 0.15, 0.85, 0.45)
        if not normalized:
            cum_plot.addStack(MCS if not scaleToData else newStack, "HIST", 1, 1)
            if DATA: cum_plot.addHisto(DATA, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
        else:
            for h in hists:
                if onlyMC and not h.GetName() == 'data':
                    cum_plot.addHisto(h.GetCumulative(), 'hist,same' if not h.GetName() == 'w dPhi cut' else 'p,same', h.GetName(), 'PL' if h.GetName() == 'w dPhi cut' else 'L', h.GetLineColor(), 1, 0)
        cum_plot.save(1, 1, logx, lumi)
    
    if returnplot:
        return plot
    else:
        del plot
        return

def rebinEWino(histo):
    #a= [0,50,100,150,225,300,325]
    a= [50,100,150,225,300,325]
    _arr = array.array('d', a)
    h_ret = r.TH1F(histo.GetName(),histo.GetTitle(), len(a)-1, _arr)
    for ii,ib in enumerate(a):
        if ib == a[-1]: continue
        tmp_err = r.Double()
        if not ib == 300:
            tmp_cont = histo.IntegralAndError(histo.FindBin(ib), histo.FindBin(a[ii+1]-1.), tmp_err )
        else:
            tmp_cont = histo.IntegralAndError(histo.FindBin(ib), histo.GetNbinsX()+1, tmp_err )
        h_ret.SetBinContent(ii+1, tmp_cont)
        h_ret.SetBinError  (ii+1, tmp_err )
    return copy.deepcopy(h_ret)

def makeSummaryEWino(plot_EWino):
    f_ttbar = r.TFile('plots/ttbarMET/lumi12.9invfbTTPOW/closure_highStats_withDataCorrected.root','READ')
    #f_vince = r.TFile('plots/ttbarMET/lumi7.65invfbTTPOW/data_SR_EWK_hists_4fb.root','READ')
    #f_vince = r.TFile('h_ewk_prediction_7p65invfb.root','READ')
    #f_vince = r.TFile('h_ewk_prediction_12p9invfb.root','READ')
    f_vince = r.TFile('h_ewk_prediction_12p9invfb_EWKsubtracted.root','READ')
    c_ttbar = f_ttbar.Get('c1_n2')
    h_ttPredMC = c_ttbar.FindObject("estimated ttbar")
    h_ttPredDA = c_ttbar.FindObject("data est. corrected")
    #h_vince = copy.deepcopy(f_vince.Get('h_templ_met'))
    h_vince = copy.deepcopy(f_vince.Get('methist_data_SR_EWK'))
    h_vince.Rebin(25)
    
    h_vincePred = h_ttPredMC.Clone('vincePredictionFinal')
    h_vincePred.Reset()
    for ib in range(1,h_vincePred.GetNbinsX()+2):
        h_vincePred.SetBinContent(ib, h_vince.GetBinContent(ib))
        h_vincePred.SetBinError  (ib, h_vince.GetBinError  (ib))
        if ib == h_vincePred.FindBin(301.):
            tmp_err = r.Double()
            tmp_cont = h_vince.IntegralAndError(ib, h_vince.GetNbinsX()+1, tmp_err)
            h_vincePred.SetBinContent(ib, tmp_cont)
            h_vincePred.SetBinError  (ib, tmp_err )
            break

    for i,h_ew in enumerate(plot_EWino.histos):
        if 'WZ'   in h_ew.GetName():
            h_wz = copy.deepcopy(h_ew)
        if 'ZZ'   in h_ew.GetName():
            h_zz = copy.deepcopy(h_ew)
        if 'TTZ'   in h_ew.GetName():
            h_ttZ = copy.deepcopy(h_ew)
        if 'rare' in h_ew.GetName():
            h_rare = copy.deepcopy(h_ew)
        if 'DATA' in h_ew.GetName():
            h_data = copy.deepcopy(h_ew)
    
    h_vincePred = rebinEWino(h_vincePred); h_vincePred.SetFillColor(r.kAzure+5)
    h_wz        = rebinEWino(h_wz       ); h_wz       .SetFillColor(r.kRed+2); h_wz.GetXaxis().SetTitle('MET (GeV)')
    h_zz        = rebinEWino(h_zz       ); h_zz       .SetFillColor(r.kGreen+1)
    h_rare      = rebinEWino(h_rare     ); h_rare     .SetFillColor(r.kGray)
    h_ttPredDA  = rebinEWino(h_ttPredDA ); h_ttPredDA .SetFillColor(r.kYellow+1)
    h_data      = rebinEWino(h_data     );

    h_vincePred .SetName('DY')
    h_wz        .SetName('WZ')
    h_zz        .SetName('ZZ')
    h_rare      .SetName('Other rares')
    h_ttPredDA  .SetName('ttbar')
    h_data      .SetName('data')

    e_vincePred = h_vincePred .Clone('err_vincePred')
    e_wz        = h_wz        .Clone('err_wz')
    e_zz        = h_zz        .Clone('err_zz')
    e_rare      = h_rare      .Clone('err_rare')
    e_ttPredDA  = h_ttPredDA  .Clone('err_ttPredDA')

    for i in range(1,e_vincePred.GetNbinsX()+1):
        e_vincePred .SetBinContent(i,1.  ); e_vincePred .SetBinError(i,0.15)
        e_wz        .SetBinContent(i,1.  ); e_wz        .SetBinError(i,0.30)
        e_zz        .SetBinContent(i,1.35); e_zz        .SetBinError(i,0.50)
        e_rare      .SetBinContent(i,1.  ); e_rare      .SetBinError(i,0.50)
        e_ttPredDA  .SetBinContent(i,1.  ); e_ttPredDA  .SetBinError(i,0.30)
        if e_vincePred.GetBinCenter(i) > 150: e_vincePred .SetBinError(i,math.sqrt(0.05**2+0.10**2))
        if e_vincePred.GetBinCenter(i) > 295: e_vincePred .SetBinError(i,math.sqrt(0.05**2+0.10**2))
        if e_vincePred.GetBinCenter(i) > 300: e_vincePred .SetBinError(i,math.sqrt(0.05**2+0.15**2))

    h_vincePred .Multiply(e_vincePred );
    h_wz        .Multiply(e_wz        )
    h_zz        .Multiply(e_zz        )
    h_rare      .Multiply(e_rare      )
    h_ttPredDA  .Multiply(e_ttPredDA  )

    bin75 = h_data.FindBin(75.)
    bkg = h_wz.GetBinContent(bin75)+h_zz.GetBinContent(bin75)+h_rare.GetBinContent(bin75)+h_ttPredDA.GetBinContent(bin75)
    data = h_data.GetBinContent(bin75)
    dy = h_vincePred.GetBinContent(bin75)

    h_vincePred.Scale((data-bkg)/dy)

    ## UGLY UGLY UGLY
    h_vincePred.SetBinError(e_vincePred.GetNbinsX(),1.56); h_vincePred.SetBinError(e_vincePred.GetNbinsX()-1,2.32)

    stack = r.THStack()
    stack.Add(h_wz       )
    stack.Add(h_zz       )
    stack.Add(h_rare     )
    stack.Add(h_ttPredDA )
    stack.Add(h_vincePred)
    can_aux = r.TCanvas("tmpcanvas")
    can_aux.cd()
    stack.Draw()
    stack.GetXaxis().SetTitle('E_{T}^{miss} [GeV]')
    stack.GetYaxis().SetTitle('Events')
    del can_aux


    totalBKG = h_vincePred.Clone('totalBkg')
    totalBKG.Add(h_wz)
    totalBKG.Add(h_zz)
    totalBKG.Add(h_rare)
    totalBKG.Add(h_ttPredDA)
    totalBKG.SetFillColor(1)
    totalBKG.SetMarkerSize(0.)
    totalBKG.SetFillStyle(3004)

    res = {}
    for i in range(1,totalBKG.GetNbinsX()+1):
        res['ra%d'%i]                  = h_rare     .GetBinContent(i)
        res['ra%d_e'%i]                = h_rare     .GetBinError  (i)
        res['wz%d'%i]                  = h_wz       .GetBinContent(i)
        res['wz%d_e'%i]                = h_wz       .GetBinError  (i)
        res['zz%d'%i]                  = h_zz       .GetBinContent(i)
        res['zz%d_e'%i]                = h_zz       .GetBinError  (i)
        res['dy%d'%i]                  = h_vincePred.GetBinContent(i)
        res['dy%d_e'%i]                = h_vincePred.GetBinError  (i)
        res['tt%d'%i]                  = h_ttPredDA .GetBinContent(i)
        res['tt%d_e'%i]                = h_ttPredDA .GetBinError  (i)
        res['to%d'%i]                  = totalBKG .GetBinContent(i)
        res['to%d_e'%i]                = totalBKG .GetBinError  (i)
        res['ob%d'%i]                  = h_data     .GetBinContent(i)
        res['ob%d_e'%i]                = h_data     .GetBinError  (i)
        res['di%d'%i], res['di%d_e'%i] = helper.ratioError(h_data.GetBinContent(i), h_data.GetBinError(i), totalBKG.GetBinContent(i), totalBKG.GetBinError(i) )

    stack.SetMaximum(2.0*max(totalBKG.GetMaximum(),h_data.GetMaximum()) )
    helper.ensureDirectory('ttbarMET/%s/'%('lumi12.9invfb'))
    plot = Canvas.Canvas('ttbarMET/%s/plot_resultFancyEWino'%('lumi12.9invfb'), 'png,pdf,root', 0.7, 0.65, 0.90, 0.9)
    plot.addStack(stack, "HIST", 1, 1)
    plot.addHisto(totalBKG, "E2,SAME", "", "", r.kBlack, 1, -1)
    plot.addHisto(h_data, "E1,SAME", "Data", "PL", r.kBlack, 1, 0)
    plot.changeLabelsToNames()
    plot.saveRatio(1, 1, 0, 12.9, h_data, totalBKG)

    
    tableResultsEWino = '''\\documentclass{{article}}
\\usepackage[a4paper,margin=1in,landscape]{{geometry}}
\\begin{{document}}
\\begin{{table}}[hbtp] 
\\begin{{center}} 
\\bgroup 
\\def\\arraystretch{{1.2}} 
\\caption{{Result for the EWKino search for 12.9 fb$^{{-1}}$.}} 
\\label{{tab:ewinoResult}} 
\\begin{{tabular}}{{l| c c c c c}} 
 MET region      & 50 -- 100 GeV                   & 100 -- 150 GeV                  & 150 -- 225 GeV                  & 225 -- 300 GeV                 & $\\geq$ 300 GeV \\\\ \\hline
rare             & {ra1:.2f} $\\pm$ {ra1_e:.2f} & {ra2:.2f} $\\pm$ {ra2_e:.2f} & {ra3:.2f} $\\pm$ {ra3_e:.2f} & {ra4:.2f} $\\pm$ {ra4_e:.2f}& {ra5:.2f} $\\pm$ {ra5_e:.2f}           \\\\
WZ3l             & {wz1:.2f} $\\pm$ {wz1_e:.2f} & {wz2:.2f} $\\pm$ {wz2_e:.2f} & {wz3:.2f} $\\pm$ {wz3_e:.2f} & {wz4:.2f} $\\pm$ {wz4_e:.2f}& {wz5:.2f} $\\pm$ {wz5_e:.2f}           \\\\
ZZ               & {zz1:.2f} $\\pm$ {zz1_e:.2f} & {zz2:.2f} $\\pm$ {zz2_e:.2f} & {zz3:.2f} $\\pm$ {zz3_e:.2f} & {zz4:.2f} $\\pm$ {zz4_e:.2f}& {zz5:.2f} $\\pm$ {zz5_e:.2f}           \\\\
DY prediction    & {dy1:.2f} $\\pm$ {dy1_e:.2f} & {dy2:.2f} $\\pm$ {dy2_e:.2f} & {dy3:.2f} $\\pm$ {dy3_e:.2f} & {dy4:.2f} $\\pm$ {dy4_e:.2f}& {dy5:.2f} $\\pm$ {dy5_e:.2f}           \\\\
tt prediction    & {tt1:.2f} $\\pm$ {tt1_e:.2f} & {tt2:.2f} $\\pm$ {tt2_e:.2f} & {tt3:.2f} $\\pm$ {tt3_e:.2f} & {tt4:.2f} $\\pm$ {tt4_e:.2f}& {tt5:.2f} $\\pm$ {tt5_e:.2f}           \\\\ \\hline
total bkg        & {to1:.1f} $\\pm$ {to1_e:.1f} & {to2:.1f} $\\pm$ {to2_e:.1f} & {to3:.1f} $\\pm$ {to3_e:.1f} & {to4:.1f} $\\pm$ {to4_e:.1f}& {to5:.1f} $\\pm$ {to5_e:.1f}           \\\\ \\hline
observed         & {ob1:.0f}                    & {ob2:.0f}                    & {ob3:.0f}                    & {ob4:.0f}                   & {ob5:.0f}                          \\\\
obs./pred.       & {di1:.2f} $\\pm$ {di1_e:.2f} & {di2:.2f} $\\pm$ {di2_e:.2f} & {di3:.2f} $\\pm$ {di3_e:.2f} & {di4:.2f} $\\pm$ {di4_e:.2f}& {di5:.2f} $\\pm$ {di5_e:.2f}      \\\\
\\end{{tabular}} 
\\egroup 
\\end{{center}} 
\\end{{table}} 
\\end{{document}}'''.format( ## try calling this with **res. this should work! amazing!
ra1=res['ra1'],  ra1_e=res['ra1_e'], ra2=res['ra2'], ra2_e=res['ra2_e'], ra3=res['ra3'], ra3_e=res['ra3_e'], ra4=res['ra4'], ra4_e=res['ra4_e'], ra5=res['ra5'], ra5_e=res['ra5_e'],
wz1=res['wz1'],  wz1_e=res['wz1_e'], wz2=res['wz2'], wz2_e=res['wz2_e'], wz3=res['wz3'], wz3_e=res['wz3_e'], wz4=res['wz4'], wz4_e=res['wz4_e'], wz5=res['wz5'], wz5_e=res['wz5_e'],
zz1=res['zz1'],  zz1_e=res['zz1_e'], zz2=res['zz2'], zz2_e=res['zz2_e'], zz3=res['zz3'], zz3_e=res['zz3_e'], zz4=res['zz4'], zz4_e=res['zz4_e'], zz5=res['zz5'], zz5_e=res['zz5_e'],
dy1=res['dy1'],  dy1_e=res['dy1_e'], dy2=res['dy2'], dy2_e=res['dy2_e'], dy3=res['dy3'], dy3_e=res['dy3_e'], dy4=res['dy4'], dy4_e=res['dy4_e'], dy5=res['dy5'], dy5_e=res['dy5_e'],
tt1=res['tt1'],  tt1_e=res['tt1_e'], tt2=res['tt2'], tt2_e=res['tt2_e'], tt3=res['tt3'], tt3_e=res['tt3_e'], tt4=res['tt4'], tt4_e=res['tt4_e'], tt5=res['tt5'], tt5_e=res['tt5_e'],
to1=res['to1'],  to1_e=res['to1_e'], to2=res['to2'], to2_e=res['to2_e'], to3=res['to3'], to3_e=res['to3_e'], to4=res['to4'], to4_e=res['to4_e'], to5=res['to5'], to5_e=res['to5_e'],
ob1=res['ob1'],                      ob2=res['ob2'],                     ob3=res['ob3'],                     ob4=res['ob4'],                     ob5=res['ob5'],                   
di1=res['di1'],  di1_e=res['di1_e'], di2=res['di2'], di2_e=res['di2_e'], di3=res['di3'], di3_e=res['di3_e'], di4=res['di4'], di4_e=res['di4_e'], di5=res['di5'], di5_e=res['di5_e'])
    helper.ensureDirectory('plots/ttbarMET/lumi12.9invfb/tables/')
    compTableFile = open('plots/ttbarMET/lumi12.9invfb/tables/resultTable12p9invfb.tex','w')
    compTableFile.write(tableResultsEWino)
    compTableFile.close()

    return copy.deepcopy(h_vincePred)
    
    

if __name__ == "__main__":

    print bcolors.HEADER
    print '#######################################################################'
    print '                  Starting r_SFOF analysis...                          '
    print '#######################################################################' + bcolors.ENDC

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientsFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    parser.add_option('-d', '--do', action='store', type=str, dest='do', default='1', help='do')
    (opts, args) = parser.parse_args()

    print bcolors.HEADER + '[Data - MC comparisons] ' + bcolors.OKBLUE + 'Loading DATA and MC trees...' + bcolors.ENDC


    
    mcDatasets = 



    daDatasetsB = ['DoubleEG_Run2017B_17Nov2017_v1_runs_297046_299329',                                             
                   'DoubleMuon_Run2017B_17Nov2017_v1_runs_297046_299329',
                   'MuonEG_Run2017B_17Nov2017_v1_runs_297046_299329',    
                   'MET_Run2017B_17Nov2017_v1_runs_297046_299329',    
                   'JetHT_Run2017B_17Nov2017_v1_runs_297046_299329',    
                   'SingleElectron_Run2017B_17Nov2017_v1_runs_297046_299329',   
                   'SingleMuon_Run2017B_17Nov2017_v1_runs_297046_299329']    
                                                                              
    daDatasetsC = ['DoubleEG_Run2017C_17Nov2017_v1_runs_299368_302029',
                   'DoubleMuon_Run2017C_17Nov2017_v1_runs_299368_302029',
                   'MuonEG_Run2017C_17Nov2017_v1_runs_299368_302029',    
                   'MET_Run2017C_17Nov2017_v1_runs_299368_302029',    
                   'JetHT_Run2017C_17Nov2017_v1_runs_299368_302029',     
                   'SingleElectron_Run2017C_17Nov2017_v1_runs_299368_302029',   
                   'SingleMuon_Run2017C_17Nov2017_v1_runs_299368_302029']   
    
    daDatasetsD = ['DoubleEG_Run2017D_17Nov2017_v1_runs_302030_303434',
                   'DoubleMuon_Run2017D_17Nov2017_v1_runs_302030_303434',
                   'MuonEG_Run2017D_17Nov2017_v1_runs_302030_303434',     
                   'MET_Run2017D_17Nov2017_v1_runs_302030_303434',     
                   'JetHT_Run2017D_17Nov2017_v1_runs_302030_303434',    
                   'SingleElectron_Run2017D_17Nov2017_v1_runs_302030_303434',     
                   'SingleMuon_Run2017D_17Nov2017_v1_runs_302030_303434']     
                                                                              
    daDatasetsE = ['DoubleEG_Run2017E_17Nov2017_v1_runs_303824_304797',
                   'DoubleMuon_Run2017E_17Nov2017_v1_runs_303824_304797',
                   'MuonEG_Run2017E_17Nov2017_v1_runs_303824_304797',    
                   'MET_Run2017E_17Nov2017_v1_runs_303824_304797',    
                   'JetHT_Run2017E_17Nov2017_v1_runs_303824_304797',    
                   'SingleElectron_Run2017E_17Nov2017_v1_runs_303824_304797',    
                   'SingleMuon_Run2017E_17Nov2017_v1_runs_303824_304797']    
                                                                              
    daDatasetsF = ['DoubleEG_Run2017F_17Nov2017_v1_runs_305040_306462',
                  'DoubleMuon_Run2017F_17Nov2017_v1_runs_305040_306462',
                  'MuonEG_Run2017F_17Nov2017_v1_runs_305040_306462',           
                  'MET_Run2017F_17Nov2017_v1_runs_305040_306462',           
                  'JetHT_Run2017F_17Nov2017_v1_runs_305040_306462',           
                  'SingleElectron_Run2017F_17Nov2017_v1_runs_305040_306462',           
                  'SingleMuon_Run2017F_17Nov2017_v1_runs_305040_306462']           
  
    daDatasets16B = ['DoubleEG_Run2016B_03Feb2017_ver2_v2_runs_273150_275376',
                   'DoubleMuon_Run2016B_03Feb2017_ver2_v2_runs_273150_275376', 
                   'MuonEG_Run2016B_03Feb2017_ver2_v2_runs_273150_275376']     
                                                                                 
    daDatasets16C = ['DoubleEG_Run2016C_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016C_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016C_03Feb2017_v1_runs_271036_284044']    
    
    daDatasets16D = ['DoubleEG_Run2016D_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016D_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016D_03Feb2017_v1_runs_271036_284044']    
                                                                                 
    daDatasets16E = ['DoubleEG_Run2016E_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016E_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016E_03Feb2017_v1_runs_271036_284044']    
                                                                                 
    daDatasets16F = ['DoubleEG_Run2016F_03Feb2017_v1_runs_271036_284044',
                  'DoubleMuon_Run2016F_03Feb2017_v1_runs_271036_284044',
                  'MuonEG_Run2016F_03Feb2017_v1_runs_271036_284044']  
                                                                                 
    daDatasets16G = ['DoubleEG_Run2016G_03Feb2017_v1_runs_271036_284044',
                   'DoubleMuon_Run2016G_03Feb2017_v1_runs_271036_284044',
                   'MuonEG_Run2016G_03Feb2017_v1_runs_271036_284044']    
                                                                                 
    daDatasets16H = ['DoubleEG_Run2016H_03Feb2017_ver2_v1_runs_281085_284035',
                   'DoubleEG_Run2016H_03Feb2017_ver3_v1_runs_284036_284044',
                   'DoubleMuon_Run2016H_03Feb2017_ver2_v1_runs_281085_284035',
                   'DoubleMuon_Run2016H_03Feb2017_ver3_v1_runs_284036_284044',
                   'MuonEG_Run2016H_03Feb2017_ver2_v1_runs_281085_284035', 
                   'MuonEG_Run2016H_03Feb2017_ver3_v1_runs_284036_284044']       



    daDatasets16 = daDatasets16B + daDatasets16C + daDatasets16D +daDatasets16E + daDatasets16F + daDatasets16G + daDatasets16H       
    daDatasets = daDatasetsB + daDatasetsC + daDatasetsD +daDatasetsE + daDatasetsF  

    
    datasetTChiWZ = ['SMS_TChiWZ']

    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0, isScan = 0)
    treeSI = Sample.Tree(helper.selectSamples(opts.sampleFile, datasetTChiWZ,'SI'),'SI', 0, (600,25,'TChiWZ')  ) 
    # treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1, isScan = 0)
    # treeDAB = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasetsB, 'DA'), 'DATA', 1, isScan = 0)
    # treeDAC = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasetsC, 'DA'), 'DATA', 1, isScan = 0)
    # treeDAD = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasetsD, 'DA'), 'DATA', 1, isScan = 0)
    # treeDAE = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasetsE, 'DA'), 'DATA', 1, isScan = 0)
    # treeDAF = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasetsF, 'DA'), 'DATA', 1, isScan = 0)
    #treeDA16 = Sample.Tree(helper.selectSamples("samplesEdge.dat", daDatasets16, 'DA'), 'DATA', 1)
    
    lumiB = 4.80
    lumiC = 9.76
    lumiD = 4.3 
    lumiE = 9.4
    lumiF = 13.5
    lumi_strB = 'runB'
    lumi_strC = 'runC'
    lumi_strD = 'runD'
    lumi_strE = 'runE'
    lumi_strF = 'runF'                                                                                              
    
    print bcolors.HEADER + '[Data - MC comparisons] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    lumi = 140. ; maxrun = 999999; lumi_str = '140.0invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    cuts = CutManager.CutManager()
    mll = 'm_{ll} [GeV]'
    met = 'p_{T}^{miss} [GeV]'
    njet = "N. Jets"
    dy = "(DY control region)"
    tt = "(ttbar control region)"
    incl = "(inclusive selection)"
    ee = ", ee"
    mm = ", #mu#mu"
    SF = ", SF"
    OF = ", OF"

    if opts.do == 'dy':
        makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DY_OF", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.OF, cuts.DYControlRegionNoMll]), cuts, mll+dy+OF,  0, True)
        makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DY_SF", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.SF, cuts.DYControlRegionNoMll]), cuts, mll+dy+SF,  0, True)
        makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DY_ee", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.ee, cuts.DYControlRegionNoMll]), cuts, mll+dy+ee,  0, True)
        makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DY_mm", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.mm, cuts.DYControlRegionNoMll]), cuts, mll+dy+mm,  0, True)
        makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DY_OF", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.OF, cuts.DYControlRegionNoMET]), cuts, met+dy+OF,  0, True)
        makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DY_SF", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.SF, cuts.DYControlRegionNoMET]), cuts, met+dy+SF,  0, True)
        makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DY_ee", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.ee, cuts.DYControlRegionNoMET]), cuts, met+dy+ee,  0, True)
        makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DY_mm", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.mm, cuts.DYControlRegionNoMET]), cuts, met+dy+mm,  0, True)
        makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "nJet_DY_OF", 7, 0, 7,  cuts.AddList([cuts.goodLepton17, cuts.OF, cuts.DYControlRegionNoJet]), cuts, njet+dy+OF,  0, True)
        makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "nJet_DY_SF", 7, 0, 7,  cuts.AddList([cuts.goodLepton17, cuts.SF, cuts.DYControlRegionNoJet]), cuts, njet+dy+SF,  0, True)
        makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "nJet_DY_ee", 7, 0, 7,  cuts.AddList([cuts.goodLepton17, cuts.ee, cuts.DYControlRegionNoJet]), cuts, njet+dy+ee,  0, True)
        makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "nJet_DY_mm", 7, 0, 7,  cuts.AddList([cuts.goodLepton17, cuts.mm, cuts.DYControlRegionNoJet]), cuts, njet+dy+mm,  0, True)
    if opts.do == 'ttbar':
        makePlot(lumi,lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbar_SF", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.SF, cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+OF,0,True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbar_OF", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.OF, cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+SF,0,True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbar_ee", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.ee, cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+ee,0,True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbar_mm", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.mm, cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+mm,0,True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "met_Edge", "met_ttbar_SF", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.SF, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+OF,  0, True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "met_Edge", "met_ttbar_OF", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.OF, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+SF,  0, True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "met_Edge", "met_ttbar_ee", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.ee, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+ee,  0, True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "met_Edge", "met_ttbar_mm", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.mm, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+mm,  0, True)    
        makePlot(lumi,lumi_str, treeDA, treeMC, "nJetSel_Edge", "nJet_ttbar_SF", 7, 0, 7,cuts.AddList([cuts.goodLepton17, cuts.SF, cuts.RSFOFDirectCRNoJet]), cuts, njet+tt+OF,  0, True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "nJetSel_Edge", "nJet_ttbar_OF", 7, 0, 7,cuts.AddList([cuts.goodLepton17, cuts.OF, cuts.RSFOFDirectCRNoJet]), cuts, njet+tt+SF,  0, True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "nJetSel_Edge", "nJet_ttbar_ee", 7, 0, 7,cuts.AddList([cuts.goodLepton17, cuts.ee, cuts.RSFOFDirectCRNoJet]), cuts, njet+tt+ee,  0, True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "nJetSel_Edge", "nJet_ttbar_mm", 7, 0, 7,cuts.AddList([cuts.goodLepton17, cuts.mm, cuts.RSFOFDirectCRNoJet]), cuts, njet+tt+mm,  0, True)
    if opts.do == 'incl':
        makePlot(lumi,lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_incl_SF", 7, 0, 7, cuts.AddList([cuts.goodLepton17, cuts.SF]), cuts, njet+incl+SF, 0, True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_incl_OF", 7, 0, 7, cuts.AddList([cuts.goodLepton17, cuts.OF]), cuts, njet+incl+OF, 0, True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_incl_ee", 7, 0, 7, cuts.AddList([cuts.goodLepton17, cuts.ee]), cuts, njet+incl+ee, 0, True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_incl_mm", 7, 0, 7, cuts.AddList([cuts.goodLepton17, cuts.mm]), cuts, njet+incl+mm, 0, True)                               
        makePlot(lumi,lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_incl_SF", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.SF]), cuts, mll+incl+SF,  0, True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_incl_OF", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.OF]), cuts, mll+incl+OF,  0, True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_incl_ee", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.ee]), cuts, mll+incl+ee,  0, True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_incl_mm", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.mm]), cuts, mll+incl+mm,  0, True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "met_Edge", "met_incl_SF", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.SF]), cuts, met+incl+SF,  0, True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "met_Edge", "met_incl_OF", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.OF]), cuts, met+incl+OF,  0, True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "met_Edge", "met_incl_ee", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.ee]), cuts, met+incl+ee,  0, True)
        makePlot(lumi,lumi_str, treeDA, treeMC, "met_Edge", "met_incl_mm", 20, 20, 300, cuts.AddList([cuts.goodLepton17, cuts.mm]), cuts, met+incl+mm,  0, True)        
    
    if opts.do == 'dy_splittedRuns':
        makePlot(lumiB, lumi_strB, treeDAB, treeMC, "lepsMll_Edge", "mll_DY_OF_runB", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+OF+", runB",1,0)
        makePlot(lumiB, lumi_strB, treeDAB, treeMC, "lepsMll_Edge", "mll_DY_SF_runB", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+SF+", runB",1,0)
        makePlot(lumiB, lumi_strB, treeDAB, treeMC, "lepsMll_Edge", "mll_DY_ee_runB", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+ee+", runB",1,0)
        makePlot(lumiB, lumi_strB, treeDAB, treeMC, "lepsMll_Edge", "mll_DY_mm_runB", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+mm+", runB",1,0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC, "lepsMll_Edge", "mll_DY_OF_runC", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+OF+", runC",1,0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC, "lepsMll_Edge", "mll_DY_SF_runC", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+SF+", runC",1,0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC, "lepsMll_Edge", "mll_DY_ee_runC", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+ee+", runC",1,0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC, "lepsMll_Edge", "mll_DY_mm_runC", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+mm+", runC",1,0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC, "lepsMll_Edge", "mll_DY_OF_runD", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+OF+", runD",1,0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC, "lepsMll_Edge", "mll_DY_SF_runD", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+SF+", runD",1,0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC, "lepsMll_Edge", "mll_DY_ee_runD", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+ee+", runD",1,0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC, "lepsMll_Edge", "mll_DY_mm_runD", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+mm+", runD",1,0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC, "lepsMll_Edge", "mll_DY_OF_runE", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+OF+", runE",1,0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC, "lepsMll_Edge", "mll_DY_SF_runE", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+SF+", runE",1,0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC, "lepsMll_Edge", "mll_DY_ee_runE", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+ee+", runE",1,0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC, "lepsMll_Edge", "mll_DY_mm_runE", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+mm+", runE",1,0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC, "lepsMll_Edge", "mll_DY_OF_runF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+OF+", runF",1,0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC, "lepsMll_Edge", "mll_DY_SF_runF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+SF+", runF",1,0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC, "lepsMll_Edge", "mll_DY_ee_runF", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+ee+", runF",1,0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC, "lepsMll_Edge", "mll_DY_mm_runF", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17,cuts.DYControlRegionNoMll]), cuts, mll+dy+mm+", runF",1,0)
        
        makePlot(lumiB, lumi_strB, treeDAB, treeMC, "met_Edge", "met_DY_OF_runB", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+OF+", runB",  1, 0)
        makePlot(lumiB, lumi_strB, treeDAB, treeMC, "met_Edge", "met_DY_SF_runB", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+SF+", runB",  1, 0)
        makePlot(lumiB, lumi_strB, treeDAB, treeMC, "met_Edge", "met_DY_ee_runB", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+ee+", runB",  1, 0)
        makePlot(lumiB, lumi_strB, treeDAB, treeMC, "met_Edge", "met_DY_mm_runB", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+mm+", runB",  1, 0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC, "met_Edge", "met_DY_OF_runC", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+OF+", runC",  1, 0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC, "met_Edge", "met_DY_SF_runC", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+SF+", runC",  1, 0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC, "met_Edge", "met_DY_ee_runC", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+ee+", runC",  1, 0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC, "met_Edge", "met_DY_mm_runC", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+mm+", runC",  1, 0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC, "met_Edge", "met_DY_OF_runD", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+OF+", runD",  1, 0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC, "met_Edge", "met_DY_SF_runD", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+SF+", runD",  1, 0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC, "met_Edge", "met_DY_ee_runD", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+ee+", runD",  1, 0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC, "met_Edge", "met_DY_mm_runD", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+mm+", runD",  1, 0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC, "met_Edge", "met_DY_OF_runE", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+OF+", runE",  1, 0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC, "met_Edge", "met_DY_SF_runE", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+SF+", runE",  1, 0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC, "met_Edge", "met_DY_ee_runE", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+ee+", runE",  1, 0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC, "met_Edge", "met_DY_mm_runE", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+mm+", runE",  1, 0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC, "met_Edge", "met_DY_OF_runF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+OF+", runF",  1, 0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC, "met_Edge", "met_DY_SF_runF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+SF+", runF",  1, 0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC, "met_Edge", "met_DY_ee_runF", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+ee+", runF",  1, 0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC, "met_Edge", "met_DY_mm_runF", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17, cuts.DYControlRegionNoMET]), cuts, met+dy+mm+", runF",  1, 0)
    if opts.do == 'ttbar_splittedRuns':
        makePlot(lumiB, lumi_strB, treeDAB, treeMC,"lepsMll_Edge","mll_ttbar_OF_runB", 20, 20, 300, cuts.AddList([cuts.OF,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+OF+", runB",1,0)
        makePlot(lumiB, lumi_strB, treeDAB, treeMC,"lepsMll_Edge","mll_ttbar_SF_runB", 20, 20, 300, cuts.AddList([cuts.SF,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+SF+", runB",1,0)
        makePlot(lumiB, lumi_strB, treeDAB, treeMC,"lepsMll_Edge","mll_ttbar_ee_runB", 20, 20, 300, cuts.AddList([cuts.ee,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+ee+", runB",1,0)
        makePlot(lumiB, lumi_strB, treeDAB, treeMC,"lepsMll_Edge","mll_ttbar_mm_runB", 20, 20, 300, cuts.AddList([cuts.mm,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+mm+", runB",1,0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC,"lepsMll_Edge","mll_ttbar_OF_runC", 20, 20, 300, cuts.AddList([cuts.OF,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+OF+", runC",1,0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC,"lepsMll_Edge","mll_ttbar_SF_runC", 20, 20, 300, cuts.AddList([cuts.SF,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+SF+", runC",1,0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC,"lepsMll_Edge","mll_ttbar_ee_runC", 20, 20, 300, cuts.AddList([cuts.ee,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+ee+", runC",1,0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC,"lepsMll_Edge","mll_ttbar_mm_runC", 20, 20, 300, cuts.AddList([cuts.mm,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+mm+", runC",1,0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC,"lepsMll_Edge","mll_ttbar_OF_runD", 20, 20, 300, cuts.AddList([cuts.OF,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+OF+", runD",1,0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC,"lepsMll_Edge","mll_ttbar_SF_runD", 20, 20, 300, cuts.AddList([cuts.SF,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+SF+", runD",1,0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC,"lepsMll_Edge","mll_ttbar_ee_runD", 20, 20, 300, cuts.AddList([cuts.ee,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+ee+", runD",1,0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC,"lepsMll_Edge","mll_ttbar_mm_runD", 20, 20, 300, cuts.AddList([cuts.mm,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+mm+", runD",1,0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC,"lepsMll_Edge","mll_ttbar_OF_runE", 20, 20, 300, cuts.AddList([cuts.OF,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+OF+", runE",1,0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC,"lepsMll_Edge","mll_ttbar_SF_runE", 20, 20, 300, cuts.AddList([cuts.SF,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+SF+", runE",1,0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC,"lepsMll_Edge","mll_ttbar_ee_runE", 20, 20, 300, cuts.AddList([cuts.ee,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+ee+", runE",1,0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC,"lepsMll_Edge","mll_ttbar_mm_runE", 20, 20, 300, cuts.AddList([cuts.mm,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+mm+", runE",1,0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC,"lepsMll_Edge","mll_ttbar_OF_runF", 20, 20, 300, cuts.AddList([cuts.OF,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+OF+", runF",1,0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC,"lepsMll_Edge","mll_ttbar_SF_runF", 20, 20, 300, cuts.AddList([cuts.SF,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+SF+", runF",1,0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC,"lepsMll_Edge","mll_ttbar_ee_runF", 20, 20, 300, cuts.AddList([cuts.ee,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+ee+", runF",1,0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC,"lepsMll_Edge","mll_ttbar_mm_runF", 20, 20, 300, cuts.AddList([cuts.mm,cuts.goodLepton17,cuts.RSFOFDirectCRNoMll]), cuts, mll+tt+mm+", runF",1,0)
        
        makePlot(lumiB, lumi_strB, treeDAB, treeMC,"met_Edge", "met_ttbar_OF_runB", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+OF+", runB",1, 0)
        makePlot(lumiB, lumi_strB, treeDAB, treeMC,"met_Edge", "met_ttbar_SF_runB", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+SF+", runB",1, 0)
        makePlot(lumiB, lumi_strB, treeDAB, treeMC,"met_Edge", "met_ttbar_ee_runB", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+ee+", runB",1, 0)
        makePlot(lumiB, lumi_strB, treeDAB, treeMC,"met_Edge", "met_ttbar_mm_runB", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+mm+", runB",1, 0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC,"met_Edge", "met_ttbar_OF_runC", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+OF+", runC",1, 0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC,"met_Edge", "met_ttbar_SF_runC", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+SF+", runC",1, 0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC,"met_Edge", "met_ttbar_ee_runC", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+ee+", runC",1, 0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC,"met_Edge", "met_ttbar_mm_runC", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+mm+", runC",1, 0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC,"met_Edge", "met_ttbar_OF_runD", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+OF+", runD",1, 0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC,"met_Edge", "met_ttbar_SF_runD", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+SF+", runD",1, 0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC,"met_Edge", "met_ttbar_ee_runD", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+ee+", runD",1, 0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC,"met_Edge", "met_ttbar_mm_runD", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+mm+", runD",1, 0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC,"met_Edge", "met_ttbar_OF_runE", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+OF+", runE",1, 0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC,"met_Edge", "met_ttbar_SF_runE", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+SF+", runE",1, 0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC,"met_Edge", "met_ttbar_ee_runE", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+ee+", runE",1, 0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC,"met_Edge", "met_ttbar_mm_runE", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+mm+", runE",1, 0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC,"met_Edge", "met_ttbar_OF_runF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+OF+", runF",1, 0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC,"met_Edge", "met_ttbar_SF_runF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+SF+", runF",1, 0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC,"met_Edge", "met_ttbar_ee_runF", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+ee+", runF",1, 0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC,"met_Edge", "met_ttbar_mm_runF", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17, cuts.RSFOFDirectCRNoMET]), cuts, met+tt+mm+", runF",1, 0)

    if opts.do == 'incl_splittedRuns':
        makePlot(lumiB, lumi_strB, treeDAB, treeMC,"lepsMll_Edge","mll_incl_OF_runB", 20, 20, 300, cuts.AddList([cuts.OF,cuts.goodLepton17]), cuts, mll+incl+OF+", runB",1,0)
        makePlot(lumiB, lumi_strB, treeDAB, treeMC,"lepsMll_Edge","mll_incl_SF_runB", 20, 20, 300, cuts.AddList([cuts.SF,cuts.goodLepton17]), cuts, mll+incl+SF+", runB",1,0)
        makePlot(lumiB, lumi_strB, treeDAB, treeMC,"lepsMll_Edge","mll_incl_ee_runB", 20, 20, 300, cuts.AddList([cuts.ee,cuts.goodLepton17]), cuts, mll+incl+ee+", runB",1,0)
        makePlot(lumiB, lumi_strB, treeDAB, treeMC,"lepsMll_Edge","mll_incl_mm_runB", 20, 20, 300, cuts.AddList([cuts.mm,cuts.goodLepton17]), cuts, mll+incl+mm+", runB",1,0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC,"lepsMll_Edge","mll_incl_OF_runC", 20, 20, 300, cuts.AddList([cuts.OF,cuts.goodLepton17]), cuts, mll+incl+OF+", runC",1,0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC,"lepsMll_Edge","mll_incl_SF_runC", 20, 20, 300, cuts.AddList([cuts.SF,cuts.goodLepton17]), cuts, mll+incl+SF+", runC",1,0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC,"lepsMll_Edge","mll_incl_ee_runC", 20, 20, 300, cuts.AddList([cuts.ee,cuts.goodLepton17]), cuts, mll+incl+ee+", runC",1,0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC,"lepsMll_Edge","mll_incl_mm_runC", 20, 20, 300, cuts.AddList([cuts.mm,cuts.goodLepton17]), cuts, mll+incl+mm+", runC",1,0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC,"lepsMll_Edge","mll_incl_OF_runD", 20, 20, 300, cuts.AddList([cuts.OF,cuts.goodLepton17]), cuts, mll+incl+OF+", runD",1,0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC,"lepsMll_Edge","mll_incl_SF_runD", 20, 20, 300, cuts.AddList([cuts.SF,cuts.goodLepton17]), cuts, mll+incl+SF+", runD",1,0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC,"lepsMll_Edge","mll_incl_ee_runD", 20, 20, 300, cuts.AddList([cuts.ee,cuts.goodLepton17]), cuts, mll+incl+ee+", runD",1,0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC,"lepsMll_Edge","mll_incl_mm_runD", 20, 20, 300, cuts.AddList([cuts.mm,cuts.goodLepton17]), cuts, mll+incl+mm+", runD",1,0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC,"lepsMll_Edge","mll_incl_OF_runE", 20, 20, 300, cuts.AddList([cuts.OF,cuts.goodLepton17]), cuts, mll+incl+OF+", runE",1,0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC,"lepsMll_Edge","mll_incl_SF_runE", 20, 20, 300, cuts.AddList([cuts.SF,cuts.goodLepton17]), cuts, mll+incl+SF+", runE",1,0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC,"lepsMll_Edge","mll_incl_ee_runE", 20, 20, 300, cuts.AddList([cuts.ee,cuts.goodLepton17]), cuts, mll+incl+ee+", runE",1,0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC,"lepsMll_Edge","mll_incl_mm_runE", 20, 20, 300, cuts.AddList([cuts.mm,cuts.goodLepton17]), cuts, mll+incl+mm+", runE",1,0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC,"lepsMll_Edge","mll_incl_OF_runF", 20, 20, 300, cuts.AddList([cuts.OF,cuts.goodLepton17]), cuts, mll+incl+OF+", runF",1,0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC,"lepsMll_Edge","mll_incl_SF_runF", 20, 20, 300, cuts.AddList([cuts.SF,cuts.goodLepton17]), cuts, mll+incl+SF+", runF",1,0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC,"lepsMll_Edge","mll_incl_ee_runF", 20, 20, 300, cuts.AddList([cuts.ee,cuts.goodLepton17]), cuts, mll+incl+ee+", runF",1,0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC,"lepsMll_Edge","mll_incl_mm_runF", 20, 20, 300, cuts.AddList([cuts.mm,cuts.goodLepton17]), cuts, mll+incl+mm+", runF",1,0)
        
        makePlot(lumiB, lumi_strB, treeDAB, treeMC,"met_Edge", "met_incl_OF_runB", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17]), cuts, met+incl+OF+", runB",1, 0)
        makePlot(lumiB, lumi_strB, treeDAB, treeMC,"met_Edge", "met_incl_SF_runB", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17]), cuts, met+incl+SF+", runB",1, 0)
        makePlot(lumiB, lumi_strB, treeDAB, treeMC,"met_Edge", "met_incl_ee_runB", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17]), cuts, met+incl+ee+", runB",1, 0)
        makePlot(lumiB, lumi_strB, treeDAB, treeMC,"met_Edge", "met_incl_mm_runB", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17]), cuts, met+incl+mm+", runB",1, 0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC,"met_Edge", "met_incl_OF_runC", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17]), cuts, met+incl+OF+", runC",1, 0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC,"met_Edge", "met_incl_SF_runC", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17]), cuts, met+incl+SF+", runC",1, 0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC,"met_Edge", "met_incl_ee_runC", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17]), cuts, met+incl+ee+", runC",1, 0)
        makePlot(lumiC, lumi_strC, treeDAC, treeMC,"met_Edge", "met_incl_mm_runC", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17]), cuts, met+incl+mm+", runC",1, 0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC,"met_Edge", "met_incl_OF_runD", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17]), cuts, met+incl+OF+", runD",1, 0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC,"met_Edge", "met_incl_SF_runD", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17]), cuts, met+incl+SF+", runD",1, 0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC,"met_Edge", "met_incl_ee_runD", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17]), cuts, met+incl+ee+", runD",1, 0)
        makePlot(lumiD, lumi_strD, treeDAD, treeMC,"met_Edge", "met_incl_mm_runD", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17]), cuts, met+incl+mm+", runD",1, 0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC,"met_Edge", "met_incl_OF_runE", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17]), cuts, met+incl+OF+", runE",1, 0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC,"met_Edge", "met_incl_SF_runE", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17]), cuts, met+incl+SF+", runE",1, 0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC,"met_Edge", "met_incl_ee_runE", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17]), cuts, met+incl+ee+", runE",1, 0)
        makePlot(lumiE, lumi_strE, treeDAE, treeMC,"met_Edge", "met_incl_mm_runE", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17]), cuts, met+incl+mm+", runE",1, 0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC,"met_Edge", "met_incl_OF_runF", 20, 20, 300, cuts.AddList([cuts.OF, cuts.goodLepton17]), cuts, met+incl+OF+", runF",1, 0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC,"met_Edge", "met_incl_SF_runF", 20, 20, 300, cuts.AddList([cuts.SF, cuts.goodLepton17]), cuts, met+incl+SF+", runF",1, 0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC,"met_Edge", "met_incl_ee_runF", 20, 20, 300, cuts.AddList([cuts.ee, cuts.goodLepton17]), cuts, met+incl+ee+", runF",1, 0)
        makePlot(lumiF, lumi_strF, treeDAF, treeMC,"met_Edge", "met_incl_mm_runF", 20, 20, 300, cuts.AddList([cuts.mm, cuts.goodLepton17]), cuts, met+incl+mm+", runF",1, 0)




 
    if opts.do == "boosted_reg_sergio_tests":
        # ossf pair (lo mitico: pt cut, dr(ll), mll > x
        # at least 1 ak4 jet
        # at dphi(ak4[0], met) > 0.4
        # met > 80 
        # bveto
        # mt2 > 100 
        # third lepton veto, 86 < mll  < 106
        #mycuts = [cuts.BaselineNoTriggerNoNJet.replace(cuts.dPhiJET1MET,'deltaPhi(met_phi_Edge,FatJetSel_Edge_phi) > 0.8'), cuts.ThirdLeptonVetoOLD, cuts.ZmassLoose] 
        #makePlot(lumi, lumi_str, None, treeMC,"nFatJetSel_Edge", "nFatJetSel_Edge", 5,-0.5, 4.5, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)

        # esto es la prueba sin ak4 veto
        mycuts = [cuts.BaselineNoTriggerNoNJet.replace(cuts.dPhiJET1MET,'deltaPhi(met_phi_Edge,FatJetSel_Edge_phi) > 0.4'), cuts.ThirdLeptonVetoOLD, cuts.ZmassLoose, cuts.SF, 'nFatJetSel_Edge > 0',  'FatJetSel_Edge_softDropMass  > 60', 'FatJetSel_Edge_softDropMass  < 100', 'FatJetSel_Edge_tau2[0]/FatJetSel_Edge_tau1[0] < 0.6','nBJetMedium25_Edge == 0'] 
        makePlot(lumi, lumi_str, None, treeMC, "met_Edge", "met_boosted", [0,100,200,300,400,500],0,0, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)        
        print kk 

        
        # estos son los finales 
        mycuts = [cuts.BaselineNoTriggerNoNJet.replace(cuts.dPhiJET1MET,'deltaPhi(met_phi_Edge,FatJetSel_Edge_phi) > 0.4'), cuts.ThirdLeptonVetoOLD, cuts.ZmassLoose, cuts.SF, 'nFatJetSel_Edge > 0', 'nJetSel_Edge == 1',  'FatJetSel_Edge_softDropMass  > 60', 'FatJetSel_Edge_softDropMass  < 100', 'FatJetSel_Edge_tau2[0]/FatJetSel_Edge_tau1[0] < 0.6'] 
        makePlot(lumi, lumi_str, None, treeMC, "met_Edge", "met", [0,50,100,200,300,400],0,0, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)        


        

        

        # mycuts = mycuts + ['met_Edge > 300']
        # makePlot(lumi, lumi_str, None, treeMC, "FatJetSel_Edge_softDropMass", "SDmass", 50,0, 300, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)        
        # makePlot(lumi, lumi_str, None, treeMC,"FatJetSel_Edge_tau2[0]/FatJetSel_Edge_tau1[0]", "fat_jet_tau21", 30, 0, 2,cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)
        # makePlot(lumi, lumi_str, None, treeMC, "FatJetSel_Edge_mass", "mass", 50,0, 300, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)        
        
        #
        
        print kk
        #makePlot(lumi, lumi_str, None, treeMC,"met_Edge", "met", 50,0, 1000, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)        
        #makePlot(lumi, lumi_str, None, treeMC,"FatJetSel_Edge_mass", "mass", 50,0, 300, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)        
        #makePlot(lumi, lumi_str, None, treeMC,"FatJetSel_Edge_tau2[0]/FatJetSel_Edge_tau1[0]", "fat_jet_tau21", 100, 0, 2, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)

        #makePlot(lumi, lumi_str, None, treeMC,"deltaPhi(FatJetSel_Edge_phi[0],met_phi_Edge)", "deltaPhiFatMet", 100,-4.15, 3.15, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)        

        #makePlot(lumi, lumi_str, None, treeMC,"FatJetSel_Edge_tau2[0]/FatJetSel_Edge_tau1[0]", "fat_jet_tau21", 100, 0, 2, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)
        #makePlot(lumi, lumi_str, None, treeMC,"FatJetSel_Edge_pt[0]", "fat_jet_pt", 100, 0, 2000, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)
        #makePlot(lumi, lumi_str, None, treeMC,"lepsZPt_Edge", "zpt", 100, 0, 2000, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)


        # makePlot(lumi, lumi_str, None, treeMC,"met_Edge", "fat_jet_met", 100, 0, 2000, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)
        # makePlot(lumi, lumi_str, None, treeMC,"FatJetSel_Edge_pt[0]", "fat_jet_pt", 20, 150, 1000, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)
        # makePlot(lumi, lumi_str, None, treeMC,"FatJetSel_Edge_softDropMass", "fat_jet_sdmass", 50, 0, 150, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)
        # makePlot(lumi, lumi_str, None, treeMC,"FatJetSel_Edge_mass", "fat_jet_mass", 50, 0, 150, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)
        
        # makePlot(lumi, lumi_str, None, treeMC,"met_Edge", "fat_jet_met", 100, 0, 2000, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)


        # plots with some met cut
        #mycuts = [cuts.BaselineNoTriggerNoNJet, cuts.METg80, cuts.bveto,cuts.mT2_100,  cuts.ThirdLeptonVetoOLD, cuts.ZmassLoose, 'nFatJetSel_Edge > 0','nJetSel_Edge == 1','met_Edge > 100']
        #makePlot(lumi, lumi_str, None, treeMC,"FatJetSel_Edge_pt[0]", "fat_jet_pt_metcut", 20, 150, 1000, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)
        #makePlot(lumi, lumi_str, None, treeMC,"nJetSel_Edge", "n_jet_sel_metcut", 6, -.5, 5.5, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)
        #makePlot(lumi, lumi_str, None, treeMC,"FatJetSel_Edge_softDropMass", "fat_jet_sdmass_metcut", 100, 0, 500, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)
        #makePlot(lumi, lumi_str, None, treeMC,"FatJetSel_Edge_mass", "fat_jet_mass_metcut", 100, 0, 500, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)
        #makePlot(lumi, lumi_str, None, treeMC,"FatJetSel_Edge_tau2[0]/FatJetSel_Edge_tau1[0]", "fat_jet_tau21_metcut", 100, 0, 2, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)
        #makePlot(lumi, lumi_str, None, treeMC,"met_Edge", "fat_jet_met_metcut", 100, 0, 500, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)



         
    if opts.do == "resolved_reg_sergio_tests":
        mycuts = [cuts.SF, cuts.ewinoWZNoTrigger.replace(cuts.METg100,'1'), cuts.ewinoWZExtMll.replace(cuts.ThirdLeptonVeto,cuts.ThirdLeptonVetoOLD), 'nFatJetSel_Edge == 0']
        makePlot(lumi, lumi_str, None, treeMC, "met_Edge", "met_resolved", [0,50,100,150,250,350],0,0, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)        
        print kk 


    if opts.do == "regions_sergio_summary":
        mycuts = [cuts.BaselineNoTriggerNoNJet.replace(cuts.dPhiJET1MET,'deltaPhi(met_phi_Edge,FatJetSel_Edge_phi) > 0.4'), cuts.ThirdLeptonVetoOLD, cuts.ZmassLoose, cuts.SF, 'nFatJetSel_Edge > 0',  'FatJetSel_Edge_softDropMass  > 60', 'FatJetSel_Edge_softDropMass  < 100', 'FatJetSel_Edge_tau2[0]/FatJetSel_Edge_tau1[0] < 0.6','nBJetMedium25_Edge == 0'] 
        makePlot(lumi, lumi_str, None, treeMC, "met_Edge", "met_boosted_strategyA", [0,100,200,300,400,500],0,0, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)
        makePlot(lumi, lumi_str, None, treeMC, "met_Edge", "met_boosted_strategyB", [0,100,200,300,400,500],0,0, cuts.AddList(mycuts + ['nJetSel_Edge == 0']), cuts, '',1, 0,treeSI=treeSI)        


        mycuts = [cuts.SF, cuts.ewinoWZNoTrigger.replace(cuts.METg100,'1'), cuts.ewinoWZExtMll.replace(cuts.ThirdLeptonVeto,cuts.ThirdLeptonVetoOLD) ]
        makePlot(lumi, lumi_str, None, treeMC, "met_Edge", "met_resolved_strategyB", [0,50,100,150,250,350],0,0, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)
        makePlot(lumi, lumi_str, None, treeMC, "met_Edge", "met_resolved_strategyA", [0,50,100,150,250,350],0,0, cuts.AddList(mycuts + ['nFatJetSel_Edge == 0']), cuts, '',1, 0,treeSI=treeSI)        


    if opts.do == "regions_sergio_summary_1":
        mycuts = [cuts.BaselineNoTriggerNoNJet.replace(cuts.dPhiJET1MET,'deltaPhi(met_phi_Edge,FatJetSel_Edge_phi) > 0.4'), cuts.ThirdLeptonVetoOLD, cuts.ZmassLoose, cuts.SF, 'nFatJetSel_Edge > 0',  'FatJetSel_Edge_softDropMass  > 60', 'FatJetSel_Edge_softDropMass  < 100', 'FatJetSel_Edge_tau2[0]/FatJetSel_Edge_tau1[0] < 0.6','nBJetMedium25_Edge == 0'] 
        makePlot(lumi, lumi_str, None, treeMC, "met_Edge", "met_boosted_strategyB", [0,100,200,300,400,500],0,0, cuts.AddList(mycuts + ['nJetSel_Edge < 2']), cuts, '',1, 0,treeSI=treeSI)        


    if opts.do == "regions_sergio_summary_2":
        mycuts = [cuts.SF, cuts.ewinoWZNoTrigger.replace(cuts.METg100,'1'), cuts.ewinoWZExtMll.replace(cuts.ThirdLeptonVeto,cuts.ThirdLeptonVetoOLD) ]
        makePlot(lumi, lumi_str, None, treeMC, "met_Edge", "met_resolved_strategyB", [0,50,100,150,250,350],0,0, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=treeSI)

    if opts.do == "regions_sergio_summary_3":
        mycuts = [cuts.SF, cuts.ewinoWZNoTrigger.replace(cuts.METg100,'1'), cuts.ewinoWZExtMll.replace(cuts.ThirdLeptonVeto,cuts.ThirdLeptonVetoOLD) ]
        makePlot(lumi, lumi_str, None, treeMC, "met_Edge", "met_resolved_strategyA", [0,50,100,150,250,350],0,0, cuts.AddList(mycuts + ['nFatJetSel_Edge == 0']), cuts, '',1, 0,treeSI=treeSI)        


    if opts.do == "regions_sergio_hz":
        mycuts = [cuts.SF, cuts.ewinoZHExtMll.replace(cuts.ThirdLeptonVeto,cuts.ThirdLeptonVetoOLD) ]
        makePlot(lumi, lumi_str, None, treeMC, "met_Edge", "met_resolved_strategyA", [0,50,100,150,250,350],0,0, cuts.AddList(mycuts), cuts, '',1, 0,treeSI=None)        
