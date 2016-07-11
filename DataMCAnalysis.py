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

def makePlot(lumi, lumi_str, treeDA, treeMC, var, name, nbin, xmin, xmax, theCut, cuts,labelx, logx, returnplot = False, scaleToData = False, normalized = False, cumulative = False, onlyMC = False):

    print 'returnplot', returnplot
    print 'scaletodata', scaleToData
    print 'normalized', normalized
    print 'cumulative', cumulative
    print 'onlyMC', onlyMC

    theCutDATA = cuts.AddList([theCut, cuts.trigger ])
    theCutZ = cuts.AddList([theCut, cuts.Zmass ])
    MC   = treeMC.getTH1F(lumi, "hMC_%s"%(name), var, nbin, xmin, xmax, theCut, '', labelx)
    MCS  = treeMC.getStack(lumi, "hMCS_%s"%(name), var, nbin, xmin, xmax, theCut, "", labelx)
    MCSZ  = treeMC.getStack(lumi, "hMCZS_%s"%(name), var, nbin, xmin, xmax, theCutZ, "", labelx)
    DATA = treeDA.getTH1F(lumi, "hDATA_%s"%(name), var, nbin, xmin, xmax, theCutDATA, '', labelx)

    if scaleToData:
        newStack = r.THStack()
        mcInt = MC  .Integral()
        daInt = DATA.Integral()
        for _i,_h in enumerate(MCS.GetHists()):
            _hnew = _h.Clone(_h.GetName()+'_new')
            _tmpInt = _h.Integral()
            _hnew.Scale(daInt/mcInt)
            newStack.Add(_hnew)
        MC.Scale(daInt/mcInt)

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
        DATA.SetName('Data')
        DATA.Scale(1./DATA.Integral())
        DATA.SetMarkerSize(0.8)
        DATA.SetMarkerColor(r.kBlack)
        DATA.SetLineColor(r.kBlack)
        hists.append(MC)
        hists.append(DATA)                                                                   

    maxValmc = MC.GetMaximum() #GetBinContent(MC.GetMaximumBin())
    maxValdata = DATA.GetMaximum() #GetBinContent(DATA.GetMaximumBin())
    maxVal = max(maxValmc, maxValdata)
    if not logx:
        MC.GetYaxis().SetRangeUser(0.1, 1.3*maxVal)
        MCS.SetMaximum(1.3*maxVal)
        if scaleToData:
            newStack.SetMaximum(1.3*maxVal)
        DATA.GetYaxis().SetRangeUser(0.1, 1.3*maxVal)
    else:
        MC.GetYaxis().SetRangeUser(0.1, 2.0*maxVal)
        MCS.SetMaximum(2.0*maxVal)
        if scaleToData:
            newStack.SetMaximum(2.0*maxVal)
        DATA.GetYaxis().SetRangeUser(0.1, 2.0*maxVal)

    if normalized:
        MC.GetYaxis().SetRangeUser(0.01, 1.2)
   
    print name
    SetOwnership(MC, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(MCS, 0 )   # 0 = release (not keep), 1 = keep
    SetOwnership(DATA, 0 )   # 0 = release (not keep), 1 = keep
    plot = Canvas.Canvas('DataMC/%s/plot_%s'%(lumi_str,name), 'png,pdf,root', 0.65, 0.6, 0.85, 0.9)
    if not normalized:
        plot.addStack(MCS if not scaleToData else newStack, "HIST", 1, 1)
        plot.addHisto(DATA, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
        #for h in enumerate(MCSZ.GetHists()):
        #    plot.addLatex(0.6, 0.6, "%s : %.1f"%(h.GetName, h.Integral()) )
    else:
        for h in hists:
            plot.addHisto(h, 'hist,same' if not h.GetName() == 'data' else 'p,same', h.GetName(), 'PL', h.GetLineColor(), 1, 0)
    plot.saveRatio(1, 1, logx, lumi, DATA, MC)

    if cumulative:
        cum_plot = Canvas.Canvas('DataMC/%s/plot_cumulative_%s'%(lumi_str,name), 'png,pdf,root', 0.60, 0.15, 0.80, 0.45)
        if not normalized:
            cum_plot.addStack(MCS if not scaleToData else newStack, "HIST", 1, 1)
            cum_plot.addHisto(DATA, "E1,SAME", "Data", "P", r.kBlack, 1, 0)
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
    (opts, args) = parser.parse_args()

    print bcolors.HEADER + '[Data - MC comparisons] ' + bcolors.OKBLUE + 'Loading DATA and MC trees...' + bcolors.ENDC

    dyDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO']
    #dyDatasets = ['DYJetsToLL_M10to50','DYJetsToLL_M50_HT100to200', 'DYJetsToLL_M50_HT200to400', 'DYJetsToLL_M50_HT400to600', 'DYJetsToLL_M50_HT600toInf']
    ttDatasets = ['TTJets_DiLepton', 'TTJets_DiLepton_ext']
    mcDatasets = ['ZZTo4L', 'WZTo3LNu', 'WWW', 'WWZ','WZZ', 'ZZZ',  'TTZToLLNuNu' ,'TTWToLNu', 'T_tWch', 'TBar_tWch' , 'TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT', 'TToLeptons_sch', 'TToLeptons_tch_powheg', 'TBarToLeptons_tch_powheg', 'TTHnobb_pow', 'VHToNonbb', 'WJetsToLNu_LO']
    mcDatasets += ttDatasets
    mcDatasets += dyDatasets
    #siDatasets = ['SMS_TChiWZ_ZToLL']
    daDatasets = ['DoubleMuon_Run2016B-PromptReco-v2_runs_273150_275376', 'DoubleEG_Run2016B-PromptReco-v2_runs_273150_275376', 'MuonEG_Run2016B-PromptReco-v2_runs_273150_275376',
                  'DoubleMuon_Run2016C-PromptReco-v2_runs_275420_276283', 'DoubleEG_Run2016C-PromptReco-v2_runs_275420_276283', 'MuonEG_Run2016C-PromptReco-v2_runs_275420_276283',
                  'DoubleMuon_Run2016D-PromptReco-v2_runs_276315_276811', 'DoubleEG_Run2016D-PromptReco-v2_runs_276315_276811', 'MuonEG_Run2016D-PromptReco-v2_runs_276315_276811']

    treeMC = Sample.Tree(helper.selectSamples(opts.sampleFile, mcDatasets, 'MC'), 'MC'  , 0, isScan = 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1, isScan = 0)
    #treeSI = Sample.Tree(helper.selectSamples(opts.sampleFile, siDatasets, 'SI'), 'SI', 0, isScan = 1)
    print bcolors.HEADER + '[Data - MC comparisons] ' + bcolors.OKBLUE + 'Trees successfully loaded...' + bcolors.ENDC

    lumi = 12.9 ; maxrun = 276811; lumi_str = '12.9invfb'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    cuts = CutManager.CutManager()
    labelmll = 'm_{ll} [GeV]'
    labelmet = 'E_{T}^{miss} [GeV]'
    labelnjet = "N. Jets"
    #makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_signal_OFSF_nll", [20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300], 1, 1, cuts.AddList([cuts.goodLepton, cuts.SignalRegionBaseLine]), cuts, labelmll, 0)
    
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
    plot_3l = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge"    , "met_3lregion", 14, 60, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region3l]), cuts, labelmet, 0, True)
    plot_3l_Zmass = makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge"    , "met_3lregion_Zmass", 14, 60, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region3l, cuts.Zmass]), cuts, labelmet, 0, True)
    plot_4l = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_4lregion", 20,  6, 206, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region4l]), cuts, labelmll, 0, True)
    plot_4l_Zmass = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_4lregion_Zmass", 20,  6, 206, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region4l, cuts.Zmass]), cuts, labelmll, 0, True)
    plot_ttZ = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttZregion", 20,  6, 206, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.regionttZ]), cuts, labelmll, 0, True)
    plot_ttZ_Zmass = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttZregion_Zmass", 20,  6, 206, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.regionttZ, cuts.Zmass]), cuts, labelmll, 0, True)
    ## ## ## plot_4l = makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_4lregion_met0to30_AF", 10,  0, 200, cuts.AddList([cuts.goodLepton, cuts.AF, cuts.nj0, cuts.region4l, 'met_Edge < 30']), cuts, labelmll, 0, True)
    makeSummaryTable3l4l(plot_3l_Zmass, plot_4l_Zmass, plot_ttZ_Zmass)
    #print "made the table"
    #ttbar_region = cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj2, 'run_Edge <= %d'%maxrun, 'met_Edge >150'])
    #makePlot(lumi, lumi_str, treeDA, treeMC, "nBJetMedium35_Edge", "ttbar_region_met150_2jets_of_dataScaled", 4, 0, 4, ttbar_region, cuts, 'n_{b-jets}', 0, False, True)

 #   makePlot(lumi, lumi_str, treeDA, treeMC, "nVert_Edge", "nvtx_inclusive", 20,  0,  40, cuts.AddList([cuts.goodLepton]), cuts, 'n_{vertices}', 0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "nVert_Edge", "nvtx_inclusive_OF", 20,  0,  40, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj1]), cuts, 'n_{vertices}',  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "nVert_Edge", "nvtx_inclusive_SF", 20,  0,  40, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj1]), cuts, 'n_{vertices}',  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "nVert_Edge", "nvtx_inclusive_ee", 20,  0,  40, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.nj1]), cuts, 'n_{vertices}',  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "nVert_Edge", "nvtx_inclusive_mm", 20,  0,  40, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.nj1]), cuts, 'n_{vertices}',  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_inclusive_SF", 40, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj1]), cuts, labelmll,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_inclusive_OF", 40, 20, 300, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj1]), cuts, labelmll,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_inclusive_ee", 40, 20, 300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.nj1]), cuts, labelmll,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_inclusive_mm", 40, 20, 300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.nj1]), cuts, labelmll,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_inclusive_SF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.nj1]), cuts, labelmet,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_inclusive_OF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.nj1]), cuts, labelmet,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_inclusive_ee", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.nj1]), cuts, labelmet,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_inclusive_mm", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.nj1]), cuts, labelmet,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_SF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.DYControlRegionNoMll]), cuts, labelmll,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_OF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.DYControlRegionNoMll]), cuts, labelmll,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_ee", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.DYControlRegionNoMll]), cuts, labelmll,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_DYJets_mm", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.DYControlRegionNoMll]), cuts, labelmll,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DYJets_SF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.DYControlRegionNoMET]), cuts, labelmet,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DYJets_OF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.DYControlRegionNoMET]), cuts, labelmet,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DYJets_ee", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.DYControlRegionNoMET]), cuts, labelmet,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_DYJets_mm", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.DYControlRegionNoMET]), cuts, labelmet,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_SF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_OF", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_ee", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "lepsMll_Edge", "mll_ttbarcontrol_mm", 20, 20, 300, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.RSFOFDirectControlRegionNoMll]), cuts, labelmll,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbarcontrol_SF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.SF, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbarcontrol_OF", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.OF, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbarcontrol_ee", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.ee, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "met_Edge", "met_ttbarcontrol_mm", 25, 0, 200, cuts.AddList([cuts.goodLepton, cuts.mm, cuts.RSFOFDirectControlRegionNoMET]), cuts, labelmet,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_SF", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.SF]), cuts, labelnjet,  0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_OF", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.OF]), cuts, labelnjet, 0 , True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_ee", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.ee]), cuts, labelnjet, 0, True)
 #   makePlot(lumi, lumi_str, treeDA, treeMC, "nJetSel_Edge", "njet_inclusive_mm", 9, 0, 9, cuts.AddList([cuts.goodLepton, cuts.mm]), cuts, labelnjet, 0, True)                                              





















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
