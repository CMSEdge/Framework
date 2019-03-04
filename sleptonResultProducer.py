#####################################################################
######                                                              #
###### 8=========D         OO                        8==========D   #  
###### OO                  ||                                ,88    #
###### ||                  ||                              ,88'     #  
###### ||O---\     ,adPPYb,|| ,adPPYb,d8  ,adPPYba,      ,88'       #
###### ||O---/    a8'    `Y||a8'    `Y88 a8P_____88    ,88'         #
###### ||         8b       ||8b       88 8PP'''''''  ,88'           #
###### \/         '8a,   ,d||'8a,   ,d88 '8b,   ,aa 88'             #
###### 8=========D `'8bbdP'\/  `'YbbdP'Y8  `'Ybbd8'' 8==========D   #
######                       aa,    ,88                             #
######                         'Y8bbdP'                             #
######                                                              #
#####################################################################

import ROOT as r
from   ROOT import gROOT, TCanvas,TGraph,  TFile, TF1, TH1D, TPaveStats, TStyle, TH1, SetOwnership, TMath, TGraphAsymmErrors
import math, sys, optparse, copy, re, array, subprocess

import os
from array import array
import include.nll     
import include.LeptonSF 
import include.PoissonError as PE 
import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Tables     as Tables
import include.Scans      as Scans

def getPoissonError(n_obs):

    Lambda_1sigma_low = -1.;
    Lambda_1sigma_up = -1.;
    Lambda_max = n_obs + 6*math.sqrt( n_obs );
    Lambda_min = 0.;
    Nsteps_lambda = int(1e4)
    Integral_fraction_1sigma = r.Math.gaussian_cdf(-1.,1.,0.);
    h_pdf_full = r.TH1D('h_pdf_full','Pdf for lambda',Nsteps_lambda,Lambda_min,Lambda_max)
    Nobs_max = n_obs+100 
    for i_bin in range(1,Nsteps_lambda):
        lambdaP = h_pdf_full.GetBinCenter(i_bin)   
        Poisson_sum_low = 0.
        Poisson_sum_up = 0.
        for i_obs1 in range(int(n_obs),int(Nobs_max)):
            Poisson_sum_low = TMath.Poisson(i_obs1,lambdaP) + Poisson_sum_low
        if(Poisson_sum_low > Integral_fraction_1sigma and  Lambda_1sigma_low < 0.) :
       	    Lambda_1sigma_low = lambdaP;                                 
        for i_obs2 in range(0,int(n_obs+1)):
            Poisson_sum_up = TMath.Poisson(i_obs2,lambdaP) + Poisson_sum_up
        if(Poisson_sum_up < Integral_fraction_1sigma and  Lambda_1sigma_up < 0.) :
            Lambda_1sigma_up = lambdaP;                                 
    
    del h_pdf_full
    return [n_obs - Lambda_1sigma_low, Lambda_1sigma_up-n_obs]

def saveInFile(theFile, measuredValueMC, measuredValueUncMC,  measuredValueData, measuredValueUncData):

    foutput = open(theFile + "_aux", 'w')
    for line in open(theFile).readlines():
        if line.find("rsfof") != -1 and  line.find("final") != -1 :
            if line.find("DATA") != -1:
                foutput.write('rsfof        final           DATA        %.4f      %0.4f      0.0000\n'%(measuredValueData, measuredValueUncData))
            else:
                foutput.write('rsfof        final           MC          %.4f      %0.4f      0.0000\n'%(measuredValueMC, measuredValueUncMC))           
        else:
            foutput.write(line)

    foutput.close()
    subprocess.call(['mv ' + theFile + "_aux " + theFile], shell=True)                                                                                      

                                                                                                                      
def getFinalError(stat, syst):
    err = math.sqrt(stat**2 + syst**2)
    return err                           

def getFactor(x,  xstat, xsyst, y, ystat, ysyst):
    z = x*y
    err = z* math.sqrt((xstat/x)**2 + (xsyst/x)**2+ (ystat/y)**2 + (ysyst/y)**2)
    return z, err                           



def makeRatioTable(myregion, dataMC, eta, nbs): ## for this to make sense the region should be properly binned!!
    header= 'THIS IS THE TABLE FOR %s in %s for %s b-tags'%(dataMC, eta, str(nbs))
    line0 = ' %30s &           & '  %('')
    line1 = ' %30s & %s pred.  & ' %('\multirow{3}{*}{%s}' %(myregion.name), dataMC)
    line2 = ' %30s & %s obs.   & ' %('', dataMC)
    line3 = ' %30s &  ratio    & ' %('')
    tmp_histo_obs  = myregion.mll     .getHisto(dataMC, eta)
    tmp_histo_pred = myregion.mll_pred.getHisto(dataMC, eta)
    my_range = range(1,tmp_histo_obs.GetNbinsX()+1)
    for i in my_range:
        tmp_ratio   = tmp_histo_obs.GetBinContent(i) / tmp_histo_pred.GetBinContent(i)
        tmp_ratio_e = math.sqrt( (tmp_histo_obs.GetBinError(i)/tmp_histo_obs.GetBinContent(i))**2  +  (tmp_histo_pred.GetBinError(i)/tmp_histo_obs.GetBinContent(i))**2)  * tmp_ratio
        line0 += '%.0f $<$ \\mll $<$ %.0f ' %(tmp_histo_pred.GetXaxis().GetBinLowEdge(i), tmp_histo_pred.GetXaxis().GetBinUpEdge(i))
        line1 += '  %.2f $\\pm$ %.2f      ' %(tmp_histo_pred.GetBinContent(i), tmp_histo_pred.GetBinError(i))
        line2 += '  %.2f $\\pm$ %.2f      ' %(tmp_histo_obs.GetBinContent(i), tmp_histo_obs.GetBinError(i))
        line3 += '  %.2f $\\pm$ %.2f      ' %(tmp_ratio, tmp_ratio_e)
        if i != max(my_range):
            line0+=' & '
            line1+=' & '
            line2+=' & '
            line3+=' & '
        else:
            line0+=' \\\\ '
            line1+=' \\\\ '
            line2+=' \\\\ '
            line3+=' \\\\ '
    return header, line0, line1, line2, line3                                                                                                                                                      

def getSystError(tot, stat):
    if tot > 0 and stat > 0 and tot > stat:
        syst = math.sqrt(tot**2-stat**2)
        return syst
    else:
        return stat                               

def getTotErr(syst, stat):
    tot = math.sqrt(syst**2+stat**2)
    return tot                              


def getFinalTotErr(err1, err2, err3):
    tot = math.sqrt(err1**2+err2**2+err3**2)
    return tot                            

def check(check):
    if check < 0:
        return 0.
    else:
        return check


def scaleByRSFOF(histo, rsfof, rsfof_err):
    h_rsfof = copy.deepcopy(histo)
    h_rsfof.SetName('h_rsfof')
    for i in range(1, h_rsfof.GetNbinsX()+1):
        h_rsfof.SetBinContent(i,rsfof)
        h_rsfof.SetBinError  (i,rsfof_err)
    histo.Multiply(h_rsfof)
    return histo                                 

def makeResultsTableSig(da, fs,da_OF, fs_stat, fs_syst,  zz, zz_stat, zz_syst, wz, wz_stat, wz_syst, rare, rare_stat, dy,dy_stat, dy_syst,  mc_full, sig1, sig2, sig3,sig4, signames):
    
    line0 = '\\begin{tabular}{c c c c c}  '                                                               
    line1 = '\\\\ \hline\hline   '                                                               
    line2 = ' p_{T}^{miss} [GeV]& 100-150& 150-225 & 225-300& 300+  \\\\ \hline  '                                                               
    fsStatUp = [];fsStatDn = [];fsSyst = [];fsTotUp  = [];fsTotDn  = []; rest = []
    dyStatUp = [];dyStatDn = [];dySyst = [];dyTotUp  = [];dyTotDn  = [];
    for i in range(1,5):
        print i
        rest.append(math.sqrt((zz_syst.GetBinError(i)**2) + (zz_stat.GetBinError(i)**2)+(wz_syst.GetBinError(i)**2) + (wz_stat.GetBinError(i)**2) + (rare_stat.GetBinError(i)**2) + (0.5*check(rare.GetBinContent(i))**2)))
        fsStatUp.append(getPoissonError(da_OF.GetBinContent(i))[1]*fs.GetBinContent(i)/da_OF.GetBinContent(i))
        fsStatDn.append(getPoissonError(da_OF.GetBinContent(i))[0]*fs.GetBinContent(i)/da_OF.GetBinContent(i))
        fsSyst.append(fs_syst.GetBinError(i))
    dyStatUp.append(getPoissonError(dy.GetBinContent(1))[1]*0.14)
    dyStatDn.append(getPoissonError(dy.GetBinContent(1))[0]*0.14)
    dySyst.append(0.5*0.14*check(dy.GetBinContent(1)))
    for i in range(2,5):
        dyStatUp.append(1.8*0.14)
        dyStatDn.append(0)
        dySyst.append(0.14*1.8*0.5)
    for i in range(0, 4):
        fsTotUp.append(getTotErr(fsStatUp[i], fsSyst[i]))
        fsTotDn.append(getTotErr(fsStatDn[i], fsSyst[i]))
        dyTotUp.append(getTotErr(dyStatUp[i], dySyst[i]))
        dyTotDn.append(getTotErr(dyStatDn[i], dySyst[i]))
    metLabels = ['100-150 ', '150-225 ','225-300, 300+']
    print "fs bin content ", da_OF.GetBinContent(4)
    print "mc_full ", mc_full.GetBinError(4)
    print "rest ", rest[3]
    print "fsStat pure", getPoissonError(da_OF.GetBinContent(4))[1]
    print "rsfof", fs.GetBinContent(4)/da_OF.GetBinContent(4)
    print "fsStatUp ", fsStatUp[3]
    print "fsSyst ", fsSyst[3]
    print "fsTotUp ", fsTotUp[3]
    lines = []
    print line0                                                                                                                                                      
    print line1                                                                                                                                                      
    print line2                                                                                                                                                      
    print 'FS bkg.&%.2f$^{+%.2f}_{-%.2f}$ &%.2f$^{+%.2f}_{-%.2f}$ &%.2f$^{+%.2f}_{-%.2f}$ &%.2f$^{+%.2f}_{-%.2f}$\\\ ' %(fs.GetBinContent(1), fsTotUp[0],fsTotDn[0], fs.GetBinContent(2), fsTotUp[1],fsTotDn[1], fs.GetBinContent(3), fsTotUp[2],fsTotDn[2],   fs.GetBinContent(4), fsTotUp[3],fsTotDn[3]) 
    print 'ZZ &%.2f$\\pm$ %.2f&%.2f$\\pm$%.2f&%.2f$\\pm$%.2f&%.2f$\\pm$ %.2f\\\ ' %(zz.GetBinContent(1),getTotErr(zz_stat.GetBinError(1), zz_syst.GetBinError(1)),zz.GetBinContent(2),getTotErr(zz_stat.GetBinError(2), zz_syst.GetBinError(2)),zz.GetBinContent(3),getTotErr(zz_stat.GetBinError(3), zz_syst.GetBinError(3)),zz.GetBinContent(4), getTotErr(zz_stat.GetBinError(4), zz_syst.GetBinError(4)))
    print 'WZ &%.2f$\\pm$ %.2f&%.2f$\\pm$%.2f&%.2f$\\pm$%.2f&%.2f$\\pm$ %.2f\\\ ' %(wz.GetBinContent(1),getTotErr(wz_stat.GetBinError(1), wz_syst.GetBinError(1)),wz.GetBinContent(2),getTotErr(wz_stat.GetBinError(2), wz_syst.GetBinError(2)),wz.GetBinContent(3),getTotErr(wz_stat.GetBinError(3), wz_syst.GetBinError(3)),wz.GetBinContent(4), getTotErr(wz_stat.GetBinError(4), wz_syst.GetBinError(4)))
    #print 'WZ &%.2f$\\pm$ %.2f&%.2f$\\pm$%.2f&%.2f$\\pm$%.2f&%.2f$\\pm$ %.2f\\\ ' %(wz.GetBinContent(1),wz.GetBinError(1),wz.GetBinContent(2),wz.GetBinError(2),wz.GetBinContent(3),wz.GetBinError(3),wz.GetBinContent(4),wz.GetBinError(4))
    print 'DY+jets&%.2f$^{+%.2f}_{-%.2f}$ &%.2f$^{+%.2f}_{-%.2f}$ &%.2f$^{+%.2f}_{-%.2f}$ &%.2f$^{+%.2f}_{-%.2f}$\\\ ' %(dy.GetBinContent(1), dyTotUp[0],dyTotDn[0], dy.GetBinContent(2), dyTotUp[1],dyTotDn[1], dy.GetBinContent(3), dyTotUp[2],dyTotDn[2],   dy.GetBinContent(4), dyTotUp[3],dyTotDn[3]) 
    print 'Rare processes&%.2f$\\pm$%.2f&%.2f$\\pm$%.2f&%.2f$\\pm$%.2f&%.2f$\\pm$ %.2f \\\\ \hline ' %(check(rare.GetBinContent(1)),rare.GetBinError(1),check(rare.GetBinContent(2)),rare.GetBinError(2),check(rare.GetBinContent(3)),rare.GetBinError(3), check(rare.GetBinContent(4)),rare.GetBinError(4))
    print 'Total Prediction &%.2f$^{+%.2f}_{-%.2f}$ &%.2f$^{+%.2f}_{-%.2f}$ &%.2f$^{+%.2f}_{-%.2f}$ &%.2f$^{+%.2f}_{-%.2f}$\\\ '%(mc_full.GetBinContent(1), getFinalTotErr(fsTotUp[0],dyTotUp[0],rest[0]), getFinalTotErr(fsTotDn[0],dyTotDn[0],rest[0]), mc_full.GetBinContent(2), getFinalTotErr(fsTotUp[1],dyTotUp[1],rest[1]) , getFinalTotErr(fsTotDn[1],dyTotDn[1],rest[1]), mc_full.GetBinContent(3), getFinalTotErr(fsTotUp[2],dyTotUp[2],rest[2]), getFinalTotErr(fsTotDn[2],dyTotDn[2],rest[2]), mc_full.GetBinContent(4), getFinalTotErr(fsTotUp[3],dyTotUp[3],rest[3]), getFinalTotErr(fsTotDn[3],dyTotDn[3],rest[3])) 
    print 'Observed data &%d$ &%d &%d &%d\\\ ' %(da.GetBinContent(1), da.GetBinContent(2), da.GetBinContent(3), da.GetBinContent(4)) 
    print 'm(slep):%s,m(LSP):%s&%.2f$\\pm$%.2f& %.2f$\\pm$%.2f& %.2f$\\pm$%.2f& %.2f$\\pm$%.2f \\\  ' %( signames[0][0], signames[0][1],   sig1.GetBinContent(1),sig1.GetBinError(1),sig1.GetBinContent(2),sig1.GetBinError(2),sig1.GetBinContent(3),sig1.GetBinError(3), sig1.GetBinContent(4),sig1.GetBinError(4))
    print 'm(slep):%s,m(LSP):%s&%.2f$\\pm$%.2f& %.2f$\\pm$%.2f& %.2f$\\pm$%.2f& %.2f$\\pm$%.2f \\\  ' %( signames[1][0], signames[1][1],   sig2.GetBinContent(1),sig2.GetBinError(1),sig2.GetBinContent(2),sig2.GetBinError(2),sig2.GetBinContent(3),sig2.GetBinError(3), sig2.GetBinContent(4),sig2.GetBinError(4))
    print 'm(slep):%s,m(LSP):%s&%.2f$\\pm$%.2f& %.2f$\\pm$%.2f& %.2f$\\pm$%.2f& %.2f$\\pm$%.2f \\\  ' %( signames[2][0], signames[2][1],   sig3.GetBinContent(1),sig3.GetBinError(1),sig3.GetBinContent(2),sig3.GetBinError(2),sig3.GetBinContent(3),sig3.GetBinError(3), sig3.GetBinContent(4),sig3.GetBinError(4))
    print 'm(slep):%s,m(LSP):%s&%.2f$\\pm$%.2f& %.2f$\\pm$%.2f& %.2f$\\pm$%.2f& %.2f$\\pm$%.2f \\\  ' %( signames[3][0], signames[3][1],   sig4.GetBinContent(1),sig4.GetBinError(1),sig4.GetBinContent(2),sig4.GetBinError(2),sig4.GetBinContent(3),sig4.GetBinError(3), sig4.GetBinContent(4),sig4.GetBinError(4))
    print '\\hline \hline'
    print '\end{tabular}'                                                                                                                                                         

def makeResultsTable(da, fs, of, fs_stat, fs_syst,  zz, zz_stat, zz_syst, wz, wz_stat , wz_syst, rare, rare_stat,  mc):
    
    line0 = '\\begin{tabular}{c c c c c c c}  '                                                               
    line2 = ' $\mathrm{p_{T}^{miss}} [GeV]& FS & ZZ & WZ & Rares & Total & Data \\\\ \hline  '                                                               
    line1 = '\\\\ \hline\hline   '                                                               
    metLabels = ['100-150 ', '150-225 ','225-300', '+300']
    lines = []
    print line0                                                                                                                                                      
    print line1                                                                                                                                                      
    for bin,metLabel in enumerate(metLabels):
        totUp = math.sqrt((abs(getPoissonError(of.GetBinContent(bin+1))[1])*fs.GetBinContent(bin+1)/of.GetBinContent(bin+1))**2 + (fs_syst.GetBinError(bin+1))**2 +(zz_stat.GetBinError(bin+1))**2 + (zz_syst.GetBinError(bin+1))**2 +(wz_stat.GetBinError(bin+1))**2 +(wz_syst.GetBinError(bin+1))**2 + (rare_stat.GetBinError(bin+1))**2 + 0.5*check(rare.GetBinContent(bin+1))**2)
        totDn = math.sqrt((abs(getPoissonError(of.GetBinContent(bin+1))[0])*fs.GetBinContent(bin+1)/of.GetBinContent(bin+1))**2 + (fs_syst.GetBinError(bin+1))**2 +(zz_stat.GetBinError(bin+1))**2 + (zz_syst.GetBinError(bin+1))**2 +(wz_stat.GetBinError(bin+1))**2 +(wz_syst.GetBinError(bin+1))**2 + (rare_stat.GetBinError(bin+1))**2 + 0.5*check(rare.GetBinContent(bin+1))**2)
        print '%s &   %.2f^{+ %.2f}_{ - %.2f}\\pm$ %.2f & %.2f $\\pm$ %.2f \\pm$ %.2f & %.2f $\\pm$ %.2f \\pm$ %.2f & %.2f $\\pm$ %.2f \\pm$ %.2f &  %.2f^{+ %.2f}_{ - %.2f}   & %.2f  \\\  ' %(metLabel,
                                                                                                                      fs.GetBinContent(bin+1),
                                                                                                                      abs(getPoissonError(of.GetBinContent(bin+1))[1])*fs.GetBinContent(bin+1)/of.GetBinContent(bin+1),
                                                                                                                      abs(getPoissonError(of.GetBinContent(bin+1))[0])*fs.GetBinContent(bin+1)/of.GetBinContent(bin+1),
                                                                                                                      fs_syst.GetBinError(bin+1),
                                                                                                                      zz.GetBinContent(bin+1),
                                                                                                                      zz_stat.GetBinError(bin+1),
                                                                                                                      zz_syst.GetBinError(bin+1),
                                                                                                                      wz.GetBinContent(bin+1),
                                                                                                                      wz_stat.GetBinError(bin+1),
                                                                                                                      wz_syst.GetBinError(bin+1),
                                                                                                                      check(rare.GetBinContent(bin+1)),
                                                                                                                      rare_stat.GetBinError(bin+1),
                                                                                                                      0.5*check(rare.GetBinContent(bin+1)),
                                                                                                                      mc.GetBinContent(bin+1),
                                                                                                                      totUp,
                                                                                                                      totDn,
                                                                                                                      da.GetBinContent(bin+1)
                                                                                                                      )

    print '\\hline \hline'
    print '\end{tabular}'                                                                                                                                                         



def makeClosureTests(var, specialcut = '', scutstring = '', doCumulative = False, nbins=0, xmin=0, xmax=0, save=True):

    if var == 'met':
        treevar = 'met_Edge'
        nbins, xmin, xmax = 29, 16, 306
        xlabel = 'm_{ll} (GeV)'
    else: 
        treevar = var
        xlabel = ''
    
    rsfof_da = helper.readFromFileRsfofD("ingredients.dat", "DATA") 
    rsfof_mc = helper.readFromFileRsfofD("ingredients.dat", "MC")   
    rt_da = helper.readFromFileRT("ingredients.dat", "DATA")
    rt_mc = helper.readFromFileRT("ingredients.dat", "MC")  
    rmue_a_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffA", "DATA")
    rmue_a_mc = helper.readFromFileRmueCoeff("ingredients.dat","coeffA", "MC")
    rmue_b_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffB", "DATA")
    rmue_b_mc = helper.readFromFileRmueCoeff("ingredients.dat","coeffB", "MC")

                                                                                                                                
    mc_OF = treeFS.getTH1F(lint, var+"mc_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.Zveto, cuts.OF]), '', xlabel,"1", kf)
    mc_SF = treeFS.getTH1F(lint, var+"mc_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.Zveto, cuts.SF]), '', xlabel,"1", kf)
    dy_SF = treeDY.getTH1F(lint, var+"dy_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion, cuts.Zveto, cuts.SF]), '', xlabel,"1", kf)
    mc_OF_factor = treeFS.getTH1F(lint, var+"mc_OF_factor"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,cuts.Zveto, cuts.OF]), '', xlabel, '(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]),kf)
    mc_OF_factorUp = treeFS.getTH1F(lint, var+"mc_OF_factorUp"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,cuts.Zveto, cuts.OF]), '', xlabel,'(0.5*( ({a} + {b}/Lep2_pt_Edge)*1.1 + 1/(1.1*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]), kf)
    mc_OF_factorDn = treeFS.getTH1F(lint, var+"mc_OF_factorDn"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,cuts.Zveto, cuts.OF]), '', xlabel,'(0.5*( ({a} + {b}/Lep2_pt_Edge)*0.9 + 1/(0.9*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]), kf)

    mc_OF_err = copy.deepcopy(mc_OF)
    mc_OF_direct = copy.deepcopy(mc_OF) 
    mc_OF_direct = scaleByRSFOF(mc_OF_direct, rsfof_mc[0], rsfof_mc[1])
    mc_OF_factor = getRMueError(mc_OF_factor, mc_OF_factorUp, mc_OF_factorDn)
    mc_OF_factor = scaleByRSFOF(mc_OF_factor, rt_mc[0], getFinalError(rt_mc[1],rt_mc[2]))
    result = weightedAverage( mc_OF_factor, mc_OF_direct, mc_OF)
    return result


def makeClosureTestPlots(analysis, var, specialcut = '', scutstring = '', doCumulative = False, nbins=0, xmin=0, xmax=0, save=True):
    treevar = 'met_Edge'        ; nbins = [100, 125, 150, 175, 200, 225, 265, 305]; xmin =1 ; xmax = 1; xlabel = 'p_{T}^{miss} [GeV] (M_{T2} > '+scutstring+')'
    #treevar = 'met_Edge'        ; nbins = [100, 150, 225, 300]; xmin =1 ; xmax = 1; xlabel = 'p_{T}^{miss} [GeV]'
    scan = Scans.Scan(analysis) 
    rsfof_da = helper.readFromFileRsfofD("ingredients.dat", "DATA") 
    rsfof_mc = helper.readFromFileRsfofD("ingredients.dat", "MC")   
    rt_da = helper.readFromFileRT("ingredients.dat", "DATA")
    rt_mc = helper.readFromFileRT("ingredients.dat", "MC")  
    rmue_a_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffA", "DATA")
    rmue_a_mc = helper.readFromFileRmueCoeff("ingredients.dat","coeffA", "MC")
    rmue_b_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffB", "DATA")
    rmue_b_mc = helper.readFromFileRmueCoeff("ingredients.dat","coeffB", "MC")

    da_OF = treeDA.getTH1F(lint, var+"da_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.OF]), '', xlabel, "1", kf)
    da_OF_factor = treeDA.getTH1F(lint, var+"da_OF_factor"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.OF]), '', xlabel,'(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_da[0],b=rmue_b_da[0]), kf)
    da_OF_factorUp = treeDA.getTH1F(lint, var+"da_OF_factorUp"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline,  cuts.OF]), '', xlabel,'(0.5*( ({a} + {b}/Lep2_pt_Edge)*1.1 + 1/(1.1*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]), kf)
    da_OF_factorDn = treeDA.getTH1F(lint, var+"da_OF_factorDn"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline,  cuts.OF]), '', xlabel,'(0.5*( ({a} + {b}/Lep2_pt_Edge)*0.9 + 1/(0.9*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]), kf)
    print '(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_da[0],b=rmue_b_da[0])
    for i in range(1, da_OF.GetNbinsX()+1):
        if not  da_OF.GetBinContent(i): continue
    da_OF_direct = copy.deepcopy(da_OF)
    da_OF_direct = scaleByRSFOF(da_OF_direct, rsfof_da[0], rsfof_da[1])
    da_OF_factor = getRMueError(da_OF_factor, da_OF_factorUp, da_OF_factorDn)
    da_OF_factor = scaleByRSFOF(da_OF_factor, rt_da[0], getFinalError(rt_da[1],rt_da[2]))
    da_result = weightedAverage( da_OF_factor, da_OF_direct, da_OF)                                   
    ### 
    da_prediction = copy.deepcopy( da_OF )
    da_OF.SetBinErrorOption( TH1.kPoisson)
    for i in range(1, da_prediction.GetNbinsX()+1):
        da_prediction.SetBinError(i,da_result[3].GetBinError(i))
    da_prediction.Multiply(da_result[2])                       
                                                                                                                                
    for bin, label in scan.SRLabels.items():
        bin = da_prediction.FindBin(bin)
        # to get asymmetric errors (i dont know why it doesnt work another way)
        dummyHisto = r.TH1F()
        dummyHisto.SetBinErrorOption(TH1.kPoisson)
        for i in range(1, int(da_OF.GetBinContent(bin))+1):
            dummyHisto.Fill(0.5)
        dummyHisto.GetBinContent(1), '+/-', dummyHisto.GetBinErrorUp(1), dummyHisto.GetBinErrorLow(1)
        syst = '  {value:4.1f}^{{+ {errUp:4.1f}}}_{{- {errDn:4.1f}}}'.format(value = da_prediction.GetBinContent(bin),
                                                                             errUp = da_prediction.GetBinErrorUp(bin),
                                                                             errDn = da_prediction.GetBinErrorLow(bin))
        stat = '^{{+ {errUp:4.1f}}}_{{- {errDn:4.1f}}}'.format(errUp = dummyHisto.GetBinErrorUp(1) * da_result[2].GetBinContent(bin),
                                                               errDn = dummyHisto.GetBinErrorLow(1) * da_result[2].GetBinContent(bin))
        del dummyHisto                                                                                                                 

    ## ## mll distributions
    da_SF = treeDA.getTH1F(lint, var+"da_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.SF]), '', xlabel, "1", kf)
    mc_OF = treeFS.getTH1F(lint, var+"mc_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.OF]), '', xlabel, "1", kf)
    mc_SF = treeFS.getTH1F(lint, var+"mc_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.SF]), '', xlabel, "1", kf)
    mc_ee = treeFS.getTH1F(lint, var+"mc_ee"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.ee]), '', xlabel, "1", kf)
    mc_mm = treeFS.getTH1F(lint, var+"mc_mm"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.mm]), '', xlabel, "1", kf)
    dy_SF = treeDY.getTH1F(lint, var+"dy_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.SF]), '', xlabel, "1", kf)
    mc_OF_factor = treeFS.getTH1F(lint, var+"mc_OF_factor"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.OF]), '', xlabel,'(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]), kf)
    mc_OF_factorUp = treeFS.getTH1F(lint, var+"mc_OF_factorUp"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.OF]), '', xlabel,'(0.5*( ({a} + {b}/Lep2_pt_Edge)*1.1 + 1/(1.1*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]), kf)
    mc_OF_factorDn = treeFS.getTH1F(lint, var+"mc_OF_factorDn"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slepBaseline, cuts.OF]), '', xlabel,'(0.5*( ({a} + {b}/Lep2_pt_Edge)*0.9 + 1/(0.9*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_mc[0],b=rmue_b_mc[0]), kf)

    mc_OF_err = copy.deepcopy(mc_OF)
    mc_OF_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    mc_OF_err.SetFillStyle(3004); mc_OF_err.SetMarkerSize(0.)
    dy_SF.SetFillColorAlpha(r.kGreen+2,0.5)

    mc_OF_direct = copy.deepcopy(mc_OF) 
    mc_OF_direct = scaleByRSFOF(mc_OF_direct, rsfof_mc[0], getFinalError(rsfof_mc[1], rsfof_mc[2]))
    mc_OF_factor = getRMueError(mc_OF_factor, mc_OF_factorUp, mc_OF_factorDn)
    mc_OF_factor = scaleByRSFOF(mc_OF_factor, rt_mc[0], getFinalError(rt_mc[1],rt_mc[2]))
    result = weightedAverage( mc_OF_factor, mc_OF_direct, mc_OF)
    mc_OF_rsfofScaled = result[3] 
    mc_OF_rsfofScaled_err = copy.deepcopy( mc_OF_rsfofScaled ) 
    mc_OF_rsfofScaled_err.SetFillColorAlpha(r.kBlue+1, 0.8)
    mc_OF_rsfofScaled_err.SetFillStyle(3004); mc_OF_rsfofScaled_err.SetMarkerSize(0.)          

    maxCont = max(da_prediction.GetMaximum(), mc_OF_rsfofScaled.GetMaximum())
    mc_SF.GetYaxis().SetRangeUser(0.1, 1.50*maxCont)
    mc_OF_rsfofScaled.GetYaxis().SetRangeUser(0.1, 1.50*maxCont)

    maxrat = 0.5
    for ib in range(1,da_SF.GetNbinsX()+1):
        tmp_rat = da_SF.GetBinContent(ib)/( mc_OF_rsfofScaled.GetBinContent(ib) if mc_OF_rsfofScaled.GetBinContent(ib) > 0 else 1. )
        if tmp_rat > maxrat:
            maxrat = tmp_rat                                                                                          
    mc_SF.GetYaxis().SetRangeUser(0.01, 1.5*mc_SF.GetMaximum())
    print helper.bcolors.HEADER + '[MC only closure test not scaled by RSFOF] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC

    print helper.bcolors.HEADER + '[MC only closure test scaled by RSFOF] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
    plot_closure = Canvas.Canvas('closure/%s/plot_closure_met_mcPredmcObs%s'%(lint_str,  '' if not scutstring else '_mt2'+scutstring), 'png,pdf', 0.6, 0.6, 0.75, 0.8)
    plot_closure.addHisto(mc_SF                , 'PE'       , 'MC - SF', 'PL', r.kGreen+1  , 1,  0)
    plot_closure.addHisto(mc_ee                , 'PE, same'       , 'MC - ee', 'PL', r.kRed+1  , 1,  0)
    plot_closure.addHisto(mc_mm                , 'PE, same'       , 'MC - mm', 'PL', r.kCyan+1  , 1,  0)
    plot_closure.addHisto(mc_OF_rsfofScaled_err, 'e2,same'  , ''       , 'PL', r.kBlue+1 , 1, -1)
    plot_closure.addHisto(mc_OF_rsfofScaled    , 'hist,SAME', 'MC - OF', 'L' , r.kBlue+1 , 1,  1)
    plot_closure.addLatex (0.61, 0.82, 'R_{SFOF} scaled')
    plot_closure.saveRatio(1, 0, 0, lint, mc_SF, mc_OF_rsfofScaled,  0. , 2, "SF/OF")                                                                                              
                                                                                                                                                                             
def getRMueError(norm, up, dn):
    final = copy.deepcopy(norm)
    for i in range(1,norm.GetNbinsX()+1):
        stat = norm.GetBinError(i)
        syst = max( abs(norm.GetBinContent(i) - up.GetBinContent(i)), abs(norm.GetBinContent(i) - dn.GetBinContent(i)))
        final.SetBinError(i, math.sqrt(stat**2 + syst**2))

    return final                                                                                                                 
                                  

def weightedAverage( factor, direct, unwght):
    
    rsfof_factor = copy.deepcopy(factor)
    rsfof_direct = copy.deepcopy(factor)
    rsfof_final  = copy.deepcopy(factor)
    prediction_final    = copy.deepcopy(unwght)
    
    for i in range(1, unwght.GetNbinsX()+1):
        if not unwght.GetBinContent(i): 
            rsfof_factor.SetBinContent( i, 0.)             
            rsfof_direct.SetBinContent( i, 0.) 
            rsfof_final .SetBinContent( i, 0.) 
            prediction_final   .SetBinContent( i, 0.)    
            continue 
        rsfof_factor.SetBinContent(i,factor.GetBinContent(i) / unwght.GetBinContent(i))
        rsfof_factor.SetBinError  (i, math.sqrt(abs(factor.GetBinError(i)**2 - (rsfof_factor.GetBinContent(i)*unwght.GetBinError(i))**2)) / unwght.GetBinContent(i) )
        rsfof_direct.SetBinContent(i, direct.GetBinContent(i) / unwght.GetBinContent(i))
        rsfof_direct.SetBinError  (i, math.sqrt(abs(direct.GetBinError(i)**2 - (rsfof_direct.GetBinContent(i)*unwght.GetBinError(i))**2)) / unwght.GetBinContent(i) )
        
        rsfof_final .SetBinContent(i, (rsfof_factor.GetBinContent(i) / rsfof_factor.GetBinError(i)**2 + rsfof_direct.GetBinContent(i) / rsfof_factor.GetBinError(i)**2) / ( 1. / rsfof_factor.GetBinError(i)**2 + 1. / rsfof_factor.GetBinError(i)**2 ))
        rsfof_final .SetBinError  (i, 1 / math.sqrt(1/rsfof_factor.GetBinError(i)**2 + 1/rsfof_direct.GetBinError(i)**2))
        
        prediction_final.SetBinContent(i, unwght.GetBinContent(i)*rsfof_final.GetBinContent(i))
        prediction_final.SetBinError  (i, math.sqrt( (unwght.GetBinError(i)*rsfof_final.GetBinContent(i))**2 + (unwght.GetBinContent(i)*rsfof_final.GetBinError(i))**2))

    return [rsfof_factor, rsfof_direct, rsfof_final, prediction_final]                                                                                                                                                                                             

def weightedAverageEEMM( factor, direct, unwght):
    
    rsfof_direct = copy.deepcopy(factor)
    prediction_final    = copy.deepcopy(unwght)
    
    for i in range(1, unwght.GetNbinsX()+1):
        if not unwght.GetBinContent(i): 
            rsfof_direct.SetBinContent( i, 0.) 
            prediction_final.SetBinContent( i, 0.)    
            continue 
        rsfof_direct.SetBinContent(i, direct.GetBinContent(i) / unwght.GetBinContent(i))
        rsfof_direct.SetBinError  (i, math.sqrt(abs(direct.GetBinError(i)**2 - (rsfof_direct.GetBinContent(i)*unwght.GetBinError(i))**2)) / unwght.GetBinContent(i) )
        
        prediction_final.SetBinContent(i, unwght.GetBinContent(i)*rsfof_direct.GetBinContent(i))
        prediction_final.SetBinError  (i, math.sqrt( (unwght.GetBinError(i)*rsfof_direct.GetBinContent(i))**2 + (unwght.GetBinContent(i)*rsfof_direct.GetBinError(i))**2))

    return [rsfof_direct, rsfof_direct, rsfof_direct, prediction_final]                                                                                                                                                                                            



def makePred(scan, specialcut = '', scutstring = '', doCumulative = False, nbins=0, xmin=0, xmax=0):
    kf = "noKFactor"
    var = scan.srID
    if var == 'mll':
        treevar = 'lepsMll_Edge'
        xlabel = 'm_{ll} (GeV)'
    elif var == 'nll':
        treevar = 'nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge)'
        xlabel = 'NLL'
    elif var == 'nllMC':
        treevar = 'nll_mc_Edge'
        nbins, xmin, xmax = 26, 10, 36
        xlabel = 'NLL'
    elif var == 'nllMCSF':
        treevar = 'nll_mc_sf_Edge'
        nbins, xmin, xmax = 26, 10, 36
        xlabel = 'NLL'
    else: 
        treevar = var
        xlabel = ''
    rsfof_da = helper.readFromFileRsfofD("ingredients.dat", "DATA") 
    rsfofee_da = helper.readFromFileRsfofDee("ingredients.dat", "DATA") 
    rsfofmm_da = helper.readFromFileRsfofDmm("ingredients.dat", "DATA") 
    rt_da = helper.readFromFileRT("ingredients.dat", "DATA")
    rmue_a_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffA","DATA")
    rmue_b_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffB","DATA")

    da_OF = treeDA.getTH1F(lint, var+"da_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,  cuts.OF, cuts.Zveto]), '', xlabel, "1", kf)
    da_OF_factor = treeDA.getTH1F(lint, var+"da_OF_factor"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,  cuts.OF, cuts.Zveto]), '', xlabel,'(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_da[0],b=rmue_b_da[0]), kf)
    da_OF_factorUp = treeDA.getTH1F(lint, var+"da_OF_factorUp"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,  cuts.OF, cuts.Zveto]), '', xlabel,'(0.5*( ({a} + {b}/Lep2_pt_Edge)*1.1 + 1/(1.1*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]), kf)
    da_OF_factorDn = treeDA.getTH1F(lint, var+"da_OF_factorDn"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.SignalRegion,  cuts.OF, cuts.Zveto]), '', xlabel,'(0.5*( ({a} + {b}/Lep2_pt_Edge)*0.9 + 1/(0.9*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]), kf)
    print '(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_da[0],b=rmue_b_da[0])
    print 'trivial check'
    for i in range(1, da_OF.GetNbinsX()+1):
        if not  da_OF.GetBinContent(i): continue

    da_OF_direct = copy.deepcopy(da_OF)
    print scutstring
    if scutstring == '':
        da_OF_direct = scaleByRSFOF(da_OF_direct, rsfof_da[0], getFinalError(rsfof_da[1], rsfof_da[2]))
    if scutstring == 'ee':
        da_OF_direct = scaleByRSFOF(da_OF_direct, rsfofee_da[0], getFinalError(rsfofee_da[1], rsfofee_da[2]))
    if scutstring == 'mm':
        da_OF_direct = scaleByRSFOF(da_OF_direct, rsfofmm_da[0], getFinalError(rsfofmm_da[1], rsfofmm_da[2]))
    da_OF_factor = getRMueError(da_OF_factor, da_OF_factorUp, da_OF_factorDn)
    da_OF_factor = scaleByRSFOF(da_OF_factor, rt_da[0], getFinalError(rt_da[1],rt_da[2]))
    result = weightedAverage( da_OF_factor, da_OF_direct, da_OF)                                   
    ### 

    prediction = copy.deepcopy( da_OF )
    da_OF.SetBinErrorOption( TH1.kPoisson)
    for i in range(1, prediction.GetNbinsX()+1):
        prediction.SetBinError(i,result[3].GetBinError(i))
    prediction.Multiply(result[2])                       

    for bin, label in scan.SRLabels.items():
        bin = prediction.FindBin(bin)
        # to get asymmetric errors (i dont know why it doesnt work another way)
        dummyHisto = r.TH1F()
        dummyHisto.SetBinErrorOption(TH1.kPoisson)
        for i in range(1, int(da_OF.GetBinContent(bin))+1):
            dummyHisto.Fill(0.5)
        dummyHisto.GetBinContent(1), '+/-', dummyHisto.GetBinErrorUp(1), dummyHisto.GetBinErrorLow(1)
        print da_OF.GetBinContent(bin), '+', da_OF.GetBinErrorUp(bin), '-', da_OF.GetBinErrorLow(bin)
        syst = '  {value:4.1f}^{{+ {errUp:4.1f}}}_{{- {errDn:4.1f}}}'.format(value = prediction.GetBinContent(bin),
                                                                             errUp = prediction.GetBinErrorUp(bin),
                                                                             errDn = prediction.GetBinErrorLow(bin))
        stat = '^{{+ {errUp:4.1f}}}_{{- {errDn:4.1f}}}'.format(errUp = dummyHisto.GetBinErrorUp(1) * result[2].GetBinContent(bin),
                                                               errDn = dummyHisto.GetBinErrorLow(1) * result[2].GetBinContent(bin))
        del dummyHisto
        print label, syst, stat                                                                                                           
    
    return result                                                                                                                                                             

def makeRSOFTable(analysis):
    scan = Scans.Scan(analysis)
    rsfof_factor_mc, rsfof_direct_mc, rsfof_final_mc, nof_final_mc = makeClosureTests(scan.srID, specialcut = '', scutstring = '', doCumulative = False, nbins=scan.srIDMax+1, xmin=-0.5, xmax=scan.srIDMax+0.5,save=False)
    rsfof_factor_da, rsfof_direct_da, rsfof_final_da, prediction_final_da = makePred(scan, specialcut = '', scutstring = '', doCumulative = False, nbins=scan.srIDMax+1, xmin=-0.5, xmax=scan.srIDMax+0.5)
    print '                &  Data       & & & MC  & \\\ \\hline'
    print '     m_{ll}  &    $\mathrm{r_{SF/OF}^{fact}}$  & $\mathrm{r_{SF/OF}^{dict}}$  &  $\mathrm{r_{SF/OF}}$ &  $\mathrm{r_{SF/OF}^{fact}}$  & $\mathrm{r_{SF/OF}^{dict}}$  &  $\mathrm{r_{SF/OF}}$ \\\ \\hline'
    for bin, label in scan.SRLabels.items():
        print ' %s &      %4.3f $\pm$ %4.3f   & %4.3f $\pm$ %4.3f & %4.3f $\pm$ %4.3f & %4.3f $\pm$ %4.3f & %4.3f $\pm$ %4.3f & %4.3f $\pm$ %4.3f \\\\'%(label, 
                rsfof_factor_da.GetBinContent(rsfof_factor_da.FindBin(bin)), rsfof_factor_da.GetBinError(rsfof_factor_da.FindBin(bin)), 
                rsfof_direct_da.GetBinContent(rsfof_direct_da.FindBin(bin)), rsfof_direct_da.GetBinError(rsfof_direct_da.FindBin(bin)), 
                rsfof_final_da.GetBinContent(  rsfof_final_da.FindBin(bin)), rsfof_final_da.GetBinError(rsfof_final_da.FindBin(bin)), 
                rsfof_factor_mc.GetBinContent(rsfof_factor_mc.FindBin(bin)), rsfof_factor_mc.GetBinError(rsfof_factor_mc.FindBin(bin)), 
                rsfof_direct_mc.GetBinContent(rsfof_direct_mc.FindBin(bin)), rsfof_direct_mc.GetBinError(rsfof_direct_mc.FindBin(bin)), 
                rsfof_final_mc.GetBinContent(  rsfof_final_mc.FindBin(bin)), rsfof_final_mc.GetBinError(rsfof_final_mc.FindBin(bin)))


def makeResultData(analysis, var, signames, maxrun = 999999, lint = 35.9, specialcut = '', scutstring = '', _options = ''):
    returnplot, addRares, splitFlavor, makeTable, printIntegral = True, True, False, False, False
    if   var == 'met'      : treevar = 'met_Edge'        ; nbins = [100, 150, 225, 300]; xmin =1 ; xmax = 1; xlabel = 'p_{T}^{miss} [GeV]'
    if not specialcut:
        specialcut = specialcut 
    else:
        specialcut = specialcut 
    scan = Scans.Scan(analysis)
    mc_stack = r.THStack()
    newLumiString = '35.9invfb'
    ##get the ingredients
    rsfof_da = helper.readFromFileRsfofD("ingredients.dat", "DATA") 
    rsfofee_da = helper.readFromFileRsfofDee("ingredients.dat", "DATA") 
    rsfofmm_da = helper.readFromFileRsfofDmm("ingredients.dat", "DATA") 
    rt_da = helper.readFromFileRT("ingredients.dat", "DATA")
    rmue_a_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffA","DATA")
    rmue_b_da = helper.readFromFileRmueCoeff("ingredients.dat","coeffB","DATA")
    #get the fs prediction
    sig1 = treeSIG1.getTH1F(lint, var+"sig1"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLeptonSignal, cuts.slep0jet, cuts.SF]), '', xlabel, "1",kf)
    sig2 = treeSIG2.getTH1F(lint, var+"sig2"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLeptonSignal, cuts.slep0jet, cuts.SF]), '', xlabel, "1",kf)
    sig3 = treeSIG3.getTH1F(lint, var+"sig3"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLeptonSignal, cuts.slep0jet, cuts.SF]), '', xlabel, "1",kf)
    sig4 = treeSIG4.getTH1F(lint, var+"sig4"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLeptonSignal, cuts.slep0jet, cuts.SF]), '', xlabel, "1",kf)
    da_SF = treeDA.getTH1F(lint, var+"da_SF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slep0jet, cuts.SF]), '', xlabel, "1",kf)
    mc_FS = treeFS.getTH1F(lint, var+"mc_FS"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slep0jet, cuts.SF]), '', xlabel, "1",kf)
    mc_OF = treeFS.getTH1F(lint, var+"mc_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([specialcut, cuts.goodLepton, cuts.slep0jet, cuts.OF]), '', xlabel, "1",kf)
    
    sig1.SetLineStyle(2);sig2.SetLineStyle(2);sig3.SetLineStyle(2);sig4.SetLineStyle(2);
    sig1.SetLineWidth(2);sig2.SetLineWidth(2);sig3.SetLineWidth(2);sig4.SetLineWidth(2);
    da_OF_orig = treeDA.getTH1F(lint, var+"da_OF_orig"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([cuts.goodLepton, cuts.slep0jet, cuts.OF]), '', xlabel, "1",kf)
    da_OF = treeDA.getTH1F(lint, var+"da_OF"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([cuts.goodLepton, cuts.slep0jet, cuts.OF]), '', xlabel, "1",kf)
    da_OF_factor = treeDA.getTH1F(lint, var+"da_OF_factor"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([cuts.goodLepton, cuts.slep0jet, cuts.OF]), '', xlabel,'(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_da[0],b=rmue_b_da[0]), kf)
    da_OF_factorUp = treeDA.getTH1F(lint, var+"da_OF_factorUp"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([cuts.goodLepton,cuts.slep0jet,  cuts.OF]), '', xlabel, '(0.5*( ({a} + {b}/Lep2_pt_Edge)*1.1 + 1/(1.1*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]), kf)
    da_OF_factorDn = treeDA.getTH1F(lint, var+"da_OF_factorDn"+scutstring, treevar, nbins, xmin, xmax, cuts.AddList([cuts.goodLepton,cuts.slep0jet,  cuts.OF]), '', xlabel, '(0.5*( ({a} + {b}/Lep2_pt_Edge)*0.9 + 1/(0.9*({a} + {b}/Lep2_pt_Edge))))'.format(a=rmue_a_da[0],b=rmue_b_da[0]), kf)
    print '(0.5*({a} + {b}/Lep2_pt_Edge + 1/({a} + {b}/Lep2_pt_Edge)))'.format(a=rmue_a_da[0],b=rmue_b_da[0])
    for i in range(1, da_OF.GetNbinsX()+1):
        if not  da_OF.GetBinContent(i): continue
    
    da_OF_direct = copy.deepcopy(da_OF)
    print "doing this flavor!", scutstring
    if scutstring == 'SF':
        da_OF_direct = scaleByRSFOF(da_OF_direct, rsfof_da[0], getFinalError(rsfof_da[1], rsfof_da[2]))
    if scutstring == 'ee':
        da_OF_direct = scaleByRSFOF(da_OF_direct, rsfofee_da[0], getFinalError(rsfofee_da[1], rsfofee_da[2]))
    if scutstring == 'mm':
        da_OF_direct = scaleByRSFOF(da_OF_direct, rsfofmm_da[0], getFinalError(rsfofmm_da[1], rsfofmm_da[2]))
    
    da_OF_factor = getRMueError(da_OF_factor, da_OF_factorUp, da_OF_factorDn)
    da_OF_factor = scaleByRSFOF(da_OF_factor, rt_da[0], getFinalError(rt_da[1],rt_da[2]))
    
    result = weightedAverage( da_OF_factor, da_OF_direct, da_OF) 
    if scutstring == 'ee':
        result = weightedAverageEEMM( da_OF_direct, da_OF_direct, da_OF)
    if scutstring == 'mm':
        result = weightedAverageEEMM( da_OF_direct, da_OF_direct, da_OF) 
    if scutstring == 'SF':
        result = weightedAverage( da_OF_factor, da_OF_direct, da_OF) 

    ### 
    prediction = copy.deepcopy( da_OF )
    prediction_stat = copy.deepcopy( da_OF_orig )
    prediction_syst = copy.deepcopy( da_OF_orig )
    da_OF.SetTitle("data OF")  
    for i in range(1, prediction.GetNbinsX()+1):
        #prediction.SetBinError(i,result[3].GetBinError(i))
        prediction.SetBinError(i,math.sqrt(result[3].GetBinError(i)**2-da_OF.GetBinError(i)**2+(getPoissonError(prediction.GetBinContent(i))[1])**2) )
        prediction_syst.SetBinError(i,math.sqrt(abs(result[3].GetBinError(i)**2-da_OF_orig.GetBinError(i)**2)))
    prediction.Multiply(result[2])                      
    prediction_syst.Multiply(result[2])                      
    prediction_stat.Multiply(result[2])                      

        # to get asymmetric errors (i dont know why it doesnt work another way)
#    dummyHisto = r.TH1F()
#    dummyHisto.Sumw2(False);
#    dummyHisto.FillRandom("gaus",100);
#    dummyHisto.SetBinErrorOption(TH1.kPoisson)
#    for i in range(1, prediction.GetNbinsX()+1):
#    #for i in range(1, int(da_OF.GetBinContent(bin))+1):
#        print "doing ", prediction.GetBinContent(i)
#        err_low = dummyHisto.GetBinErrorLow(i)
#        err_up  = dummyHisto.GetBinErrorUp(i)
#        print "err_low", err_low
#        print "err_up ", err_up
    #dummyHisto.GetBinContent(1), '+/-', dummyHisto.GetBinErrorUp(1), dummyHisto.GetBinErrorLow(1)
    #syst = '  {value:4.1f}^{{+ {errUp:4.1f}}}_{{- {errDn:4.1f}}}'.format(value = da_prediction.GetBinContent(bin),
    #                                                                     errUp = da_prediction.GetBinErrorUp(bin),
    #                                                                     errDn = da_prediction.GetBinErrorLow(bin))
    #stat = '^{{+ {errUp:4.1f}}}_{{- {errDn:4.1f}}}'.format(errUp = dummyHisto.GetBinErrorUp(1) * da_result[2].GetBinContent(bin),
    #                                                       errDn = dummyHisto.GetBinErrorLow(1) * da_result[2].GetBinContent(bin))



    doPDFvariations = False
    zz_orig = treeZZ.getTH1F(lint, var+"zz_SF", "met_Edge",nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF, cuts.lepsFromZ]), '', xlabel, "1",'ZZpt')
    wz_orig = treeWZ.getTH1F(lint, var+"wz_SF", "met_Edge",nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF, cuts.lepsFromZ]), '', xlabel, "1",kf)
    dy_orig = treeDY.getTH1F(lint, var+"dy_SF", "met_Edge",nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jetOnZ,cuts.SF]), '', xlabel, "1",'ZZpt')
    zz      = treeZZ.getTH1F(lint,"zz_SF"+scutstring, treevar,nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF, cuts.lepsFromZ]), '', xlabel, "1",'ZZpt')
    zz_syst  = treeZZ.getTH1F(lint,"zz_systSF"+scutstring, treevar,nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF, cuts.lepsFromZ]), '', xlabel, "1",'ZZpt')
    wz      = treeWZ.getTH1F(lint,"wz_SF"+scutstring, treevar,nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF, cuts.lepsFromZ]), '', xlabel,"1", kf)
    wz_syst  = treeWZ.getTH1F(lint,"wz_systSF"+scutstring, treevar,nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF, cuts.lepsFromZ]), '', xlabel,"1", kf)
    dy      = treeDY.getTH1F(lint,"dy_SF"+scutstring, treevar,nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jetOnZ,cuts.SF]), '', xlabel, "1",'ZZpt')
    dy_syst  = treeDY.getTH1F(lint,"dy_systSF"+scutstring, treevar,nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jetOnZ,cuts.SF]), '', xlabel, "1",'ZZpt')
    dy_stat  = treeDY.getTH1F(lint,"dy_statSF"+scutstring, treevar,nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jetOnZ,cuts.SF]), '', xlabel, "1",'ZZpt')
    if doPDFvariations:
        arrZZ1 = [];arrZZ2 = [];arrZZ3 = []; arrZZ4 = [];arrWZ1 = [];arrWZ2 = [];arrWZ3 = []; arrWZ4 = [];valZZ = [];valWZ = []
        fPDF = TFile("pdfVariations.root", "UPDATE");
        for i in range(9, 109, 1): #these are the indices of the 100 pdf weights 
            #print "doing ", i
            zz_lhe= treeZZ.getTH1F(lint,"z"+str(i),"met_Edge",nbins,1,1,cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF,cuts.lepsFromZ]),'',xlabel,"LHEweight_wgt_Edge["+str(i)+"]","ZZpt")
            wz_lhe= treeWZ.getTH1F(lint,"w"+str(i),"met_Edge",nbins,1,1,cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF,cuts.lepsFromZ]),'',xlabel,"LHEweight_wgt_Edge["+str(i)+"]",kf)
            valZZ.append(1)
            valWZ.append(1)
            arrZZ1.append(zz_orig.GetBinContent(1) - zz_lhe.GetBinContent(1)) 
            arrWZ1.append(wz_orig.GetBinContent(1) - wz_lhe.GetBinContent(1)) 
            arrZZ2.append(zz_orig.GetBinContent(2) - zz_lhe.GetBinContent(2)) 
            arrWZ2.append(wz_orig.GetBinContent(2) - wz_lhe.GetBinContent(2)) 
            arrZZ3.append(zz_orig.GetBinContent(3) - zz_lhe.GetBinContent(3)) 
            arrWZ3.append(wz_orig.GetBinContent(3) - wz_lhe.GetBinContent(3))
            arrZZ4.append(zz_orig.GetBinContent(4) - zz_lhe.GetBinContent(4)) 
            arrWZ4.append(wz_orig.GetBinContent(4) - wz_lhe.GetBinContent(4))
        tZZ_met1 = TGraph(len(valZZ), array("f", arrZZ1), array("f", valZZ))
        tZZ_met2 = TGraph(len(valZZ), array("f", arrZZ2), array("f", valZZ))
        tZZ_met3 = TGraph(len(valZZ), array("f", arrZZ3), array("f", valZZ))
        tZZ_met4 = TGraph(len(valZZ), array("f", arrZZ4), array("f", valZZ))
        tWZ_met1 = TGraph(len(valWZ), array("f", arrWZ1), array("f", valWZ))
        tWZ_met2 = TGraph(len(valWZ), array("f", arrWZ2), array("f", valWZ))
        tWZ_met3 = TGraph(len(valWZ), array("f", arrWZ3), array("f", valWZ))
        tWZ_met4 = TGraph(len(valWZ), array("f", arrWZ4), array("f", valWZ))
        tZZ_met1.Write("tZZ_met1");tZZ_met2.Write("tZZ_met2");tZZ_met3.Write("tZZ_met3");tZZ_met4.Write("tZZ_met4");
        tWZ_met1.Write("tWZ_met1");tWZ_met2.Write("tWZ_met2");tWZ_met3.Write("tWZ_met3");tWZ_met4.Write("tWZ_met4");
        fPDF.Close()
        #print "tZZ_met1.GetRMS()",tZZ_met1.GetRMS() 
        #print "tWZ_met1.GetRMS()",tWZ_met1.GetRMS()
        #print "tZZ_met2.GetRMS()",tZZ_met2.GetRMS() 
        #print "tWZ_met2.GetRMS()",tWZ_met2.GetRMS()
        #print "tZZ_met3.GetRMS()",tZZ_met3.GetRMS() 
        #print "tWZ_met3.GetRMS()",tWZ_met3.GetRMS()
        #print "tZZ_met4.GetRMS()",tZZ_met4.GetRMS() 
        #print "tWZ_met4.GetRMS()",tWZ_met4.GetRMS() 
    
    ZZ_met1 = 0.220514311817
    WZ_met1 = 0.111147118377
    ZZ_met2 = 0.15735404677
    WZ_met2 = 0.0519250007685
    ZZ_met3 = 0.0498430447302
    WZ_met3 = 0.0219493486014
    ZZ_met4 = 0.0364173724046
    WZ_met4 = 0.00542837400671
    ZZ_met = [ZZ_met1, ZZ_met2, ZZ_met3, ZZ_met4]
    WZ_met = [WZ_met1, WZ_met2, WZ_met3, WZ_met4]
    pdfZZ = r.TH1F('pdfsZZ','pdfZZ',4,0,4)
    pdfWZ = r.TH1F('pdfsWZ','pdfWZ',4,0,4)
    pdfZZ.SetBinContent(1, ZZ_met1 );pdfWZ.SetBinContent(1, WZ_met1);
    pdfZZ.SetBinContent(2, ZZ_met2 );pdfWZ.SetBinContent(2, WZ_met2);
    pdfZZ.SetBinContent(3, ZZ_met3 );pdfWZ.SetBinContent(3, WZ_met3);
    pdfZZ.SetBinContent(4, ZZ_met4 );pdfWZ.SetBinContent(4, WZ_met4);
    #zz_mass     = treeZZ.getTH1F(lint,"zz_massFF"+scutstring, treevar,nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF, cuts.lepsFromZ]), '', xlabel, "1","ZZmass")
    zz_pt       = treeZZ.getTH1F(lint,"zz_ptFF"+scutstring, treevar,nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF, cuts.lepsFromZ]), '', xlabel, "1","ZZpt")
    zz_noKF     = treeZZ.getTH1F(lint,"zz_noKF"+scutstring, treevar,nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF, cuts.lepsFromZ]), '', xlabel, "1","noKFactor")
    zz_noScale  = treeZZ.getTH1F(lint,"zz_noScale"+scutstring, treevar,nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF, cuts.lepsFromZ]), '', xlabel, "1","noKFactor")
    zz_qcdInclNom = treeZZ.getTH1F(lint,"zz_SF_iqNom",treevar,nbins,1,1,cuts.AddList([specialcut,cuts.goodLepton,cuts.slepInclJet,cuts.SF,cuts.lepsFromZ]),'',xlabel,"1",'ZZpt')
    zz_qcdInclUp  = treeZZ.getTH1F(lint,"zz_SF_iqUp" ,treevar,nbins,1,1,cuts.AddList([specialcut,cuts.goodLepton,cuts.slepInclJet,cuts.SF,cuts.lepsFromZ]),'',xlabel,"LHEweight_wgt_Edge[4]",'ZZpt')
    zz_qcdInclDn  = treeZZ.getTH1F(lint,"zz_SF_iqDn" ,treevar,nbins,1,1,cuts.AddList([specialcut,cuts.goodLepton,cuts.slepInclJet,cuts.SF,cuts.lepsFromZ]),'',xlabel,"LHEweight_wgt_Edge[8]",'ZZpt')
    zz_qcdExclNom = treeZZ.getTH1F(lint,"zz_SF_eqNom",treevar,nbins,1,1,cuts.AddList([specialcut,cuts.goodLepton,cuts.slepExclJet,cuts.SF,cuts.lepsFromZ]),'',xlabel,"1",'ZZpt')
    zz_qcdExclUp  = treeZZ.getTH1F(lint,"zz_SF_eqUp" ,treevar,nbins,1,1,cuts.AddList([specialcut,cuts.goodLepton,cuts.slepExclJet,cuts.SF,cuts.lepsFromZ]),'',xlabel,"LHEweight_wgt_Edge[4]",'ZZpt')
    zz_qcdExclDn  = treeZZ.getTH1F(lint,"zz_SF_eqDn" ,treevar,nbins,1,1,cuts.AddList([specialcut,cuts.goodLepton,cuts.slepExclJet,cuts.SF,cuts.lepsFromZ]),'',xlabel,"LHEweight_wgt_Edge[8]",'ZZpt')
    wz_qcdInclNom = treeWZ.getTH1F(lint,"wz_SF_iqNom",treevar,nbins,1,1,cuts.AddList([specialcut,cuts.goodLepton,cuts.slepInclJet,cuts.SF,cuts.lepsFromZ]),'',xlabel,"1",kf)
    wz_qcdInclUp  = treeWZ.getTH1F(lint,"wz_SF_iqUp" ,treevar,nbins,1,1,cuts.AddList([specialcut,cuts.goodLepton,cuts.slepInclJet,cuts.SF,cuts.lepsFromZ]),'',xlabel,"LHEweight_wgt_Edge[4]", kf)
    wz_qcdInclDn  = treeWZ.getTH1F(lint,"wz_SF_iqDn" ,treevar,nbins,1,1,cuts.AddList([specialcut,cuts.goodLepton,cuts.slepInclJet,cuts.SF,cuts.lepsFromZ]),'',xlabel,"LHEweight_wgt_Edge[8]", kf)
    wz_qcdExclNom = treeWZ.getTH1F(lint,"wz_SF_eqNom",treevar,nbins,1,1,cuts.AddList([specialcut,cuts.goodLepton,cuts.slepExclJet,cuts.SF,cuts.lepsFromZ]),'',xlabel,"1",kf)
    wz_qcdExclUp  = treeWZ.getTH1F(lint,"wz_SF_eqUp" ,treevar,nbins,1,1,cuts.AddList([specialcut,cuts.goodLepton,cuts.slepExclJet,cuts.SF,cuts.lepsFromZ]),'',xlabel,"LHEweight_wgt_Edge[4]", kf)
    wz_qcdExclDn  = treeWZ.getTH1F(lint,"wz_SF_eqDn" ,treevar,nbins,1,1,cuts.AddList([specialcut,cuts.goodLepton,cuts.slepExclJet,cuts.SF,cuts.lepsFromZ]),'',xlabel,"LHEweight_wgt_Edge[8]", kf)
    zz_jecUp    = treeZZ.getTH1F(lint,"zz_SF_jecUp", "met_jecUp_Edge",nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF, cuts.lepsFromZ]), '', xlabel, "1",'ZZpt')
    zz_jecDn    = treeZZ.getTH1F(lint,"zz_SF_jecDn", "met_jecDn_Edge",nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF, cuts.lepsFromZ]), '', xlabel, "1",'ZZpt')
    wz_jecUp    = treeWZ.getTH1F(lint,"wz_SF_jecUp", "met_jecUp_Edge",nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF, cuts.lepsFromZ]), '', xlabel, "1",kf)
    wz_jecDn    = treeWZ.getTH1F(lint,"wz_SF_jecDn", "met_jecDn_Edge",nbins, 1,1, cuts.AddList([specialcut,cuts.goodLepton,cuts.slep0jet,cuts.SF, cuts.lepsFromZ]), '', xlabel, "1",kf)
    #print "#############zzQCD INCL #############################################"
    #print "zz_qcdUp 1", abs(zz_qcdInclUp.GetBinContent(1) - zz_qcdInclNom.GetBinContent(1))
    #print "zz_qcdDn 1", abs(zz_qcdInclDn.GetBinContent(1) - zz_qcdInclNom.GetBinContent(1))
    #print "zz_qcdUp 2", abs(zz_qcdInclUp.GetBinContent(2) - zz_qcdInclNom.GetBinContent(2))
    #print "zz_qcdDn 2", abs(zz_qcdInclDn.GetBinContent(2) - zz_qcdInclNom.GetBinContent(2))
    #print "zz_qcdUp 3", abs(zz_qcdInclUp.GetBinContent(3) - zz_qcdInclNom.GetBinContent(3))
    #print "zz_qcdDn 3", abs(zz_qcdInclDn.GetBinContent(3) - zz_qcdInclNom.GetBinContent(3))
    #print "zz_qcdUp 4", abs(zz_qcdInclUp.GetBinContent(4) - zz_qcdInclNom.GetBinContent(4))
    #print "zz_qcdDn 4", abs(zz_qcdInclDn.GetBinContent(4) - zz_qcdInclNom.GetBinContent(4))
    #print "#############zzQCD EXCL #############################################"
    #print "zz_qcdUp 1", abs(zz_qcdExclUp.GetBinContent(1) - zz_qcdExclNom.GetBinContent(1))
    #print "zz_qcdDn 1", abs(zz_qcdExclDn.GetBinContent(1) - zz_qcdExclNom.GetBinContent(1))
    #print "zz_qcdUp 2", abs(zz_qcdExclUp.GetBinContent(2) - zz_qcdExclNom.GetBinContent(2))
    #print "zz_qcdDn 2", abs(zz_qcdExclDn.GetBinContent(2) - zz_qcdExclNom.GetBinContent(2))
    #print "zz_qcdUp 3", abs(zz_qcdExclUp.GetBinContent(3) - zz_qcdExclNom.GetBinContent(3))
    #print "zz_qcdDn 3", abs(zz_qcdExclDn.GetBinContent(3) - zz_qcdExclNom.GetBinContent(3))
    #print "zz_qcdUp 4", abs(zz_qcdExclUp.GetBinContent(4) - zz_qcdExclNom.GetBinContent(4))
    #print "zz_qcdDn 4", abs(zz_qcdExclDn.GetBinContent(4) - zz_qcdExclNom.GetBinContent(4))
    #print "#############zzQCD TOT #############################################"
    #print "zz_qcdUp 1",math.sqrt( abs(zz_qcdInclUp.GetBinContent(1) - zz_qcdInclNom.GetBinContent(1))**2+ abs(zz_qcdExclUp.GetBinContent(1) - zz_qcdExclNom.GetBinContent(1))**2)
    #print "zz_qcdDn 1",math.sqrt( abs(zz_qcdInclDn.GetBinContent(1) - zz_qcdInclNom.GetBinContent(1))**2+ abs(zz_qcdExclDn.GetBinContent(1) - zz_qcdExclNom.GetBinContent(1))**2)
    #print "zz_qcdUp 2",math.sqrt( abs(zz_qcdInclUp.GetBinContent(2) - zz_qcdInclNom.GetBinContent(2))**2+ abs(zz_qcdExclUp.GetBinContent(2) - zz_qcdExclNom.GetBinContent(2))**2)
    #print "zz_qcdDn 2",math.sqrt( abs(zz_qcdInclDn.GetBinContent(2) - zz_qcdInclNom.GetBinContent(2))**2+ abs(zz_qcdExclDn.GetBinContent(2) - zz_qcdExclNom.GetBinContent(2))**2)
    #print "zz_qcdUp 3",math.sqrt( abs(zz_qcdInclUp.GetBinContent(3) - zz_qcdInclNom.GetBinContent(3))**2+ abs(zz_qcdExclUp.GetBinContent(3) - zz_qcdExclNom.GetBinContent(3))**2)
    #print "zz_qcdDn 3",math.sqrt( abs(zz_qcdInclDn.GetBinContent(3) - zz_qcdInclNom.GetBinContent(3))**2+ abs(zz_qcdExclDn.GetBinContent(3) - zz_qcdExclNom.GetBinContent(3))**2)
    #print "zz_qcdUp 4",math.sqrt( abs(zz_qcdInclUp.GetBinContent(4) - zz_qcdInclNom.GetBinContent(4))**2+ abs(zz_qcdExclUp.GetBinContent(4) - zz_qcdExclNom.GetBinContent(4))**2)
    #print "zz_qcdDn 4",math.sqrt( abs(zz_qcdInclDn.GetBinContent(4) - zz_qcdInclNom.GetBinContent(4))**2+ abs(zz_qcdExclDn.GetBinContent(4) - zz_qcdExclNom.GetBinContent(4))**2)
    
   # print "#############wzQCD#############################################"
   # print "wz_qcdUp 1", abs(wz_qcdUp.GetBinContent(1) - wz_qcdNom.GetBinContent(1))
   # print "wz_qcdDn 1", abs(wz_qcdDn.GetBinContent(1) - wz_qcdNom.GetBinContent(1))
   # print "wz_qcdUp 2", abs(wz_qcdUp.GetBinContent(2) - wz_qcdNom.GetBinContent(2))
   # print "wz_qcdDn 2", abs(wz_qcdDn.GetBinContent(2) - wz_qcdNom.GetBinContent(2))
   # print "wz_qcdUp 3", abs(wz_qcdUp.GetBinContent(3) - wz_qcdNom.GetBinContent(3))
   # print "wz_qcdDn 3", abs(wz_qcdDn.GetBinContent(3) - wz_qcdNom.GetBinContent(3))
   # print "wz_qcdUp 4", abs(wz_qcdUp.GetBinContent(4) - wz_qcdNom.GetBinContent(4))
   # print "wz_qcdDn 4", abs(wz_qcdDn.GetBinContent(4) - wz_qcdNom.GetBinContent(4))
#    print "#############zzQCD#############################################"
#    print "zz_jecUp 1", abs(zz_jecUp.GetBinContent(1) - zz.GetBinContent(1))
#    print "zz_jecDn 1", abs(zz_jecDn.GetBinContent(1) - zz.GetBinContent(1))
#    print "zz_jecUp 2", abs(zz_jecUp.GetBinContent(2) - zz.GetBinContent(2))
#    print "zz_jecDn 2", abs(zz_jecDn.GetBinContent(2) - zz.GetBinContent(2))
#    print "zz_jecUp 3", abs(zz_jecUp.GetBinContent(3) - zz.GetBinContent(3))
#    print "zz_jecDn 3", abs(zz_jecDn.GetBinContent(3) - zz.GetBinContent(3))
#    print "zz_jecUp 4", abs(zz_jecUp.GetBinContent(4) - zz.GetBinContent(4))
#    print "zz_jecDn 4", abs(zz_jecDn.GetBinContent(4) - zz.GetBinContent(4))
#    print "#############wzQCD##########################################"
#    print "wz_jecUp 1", abs(wz_jecUp.GetBinContent(1) - wz.GetBinContent(1))
#    print "wz_jecDn 1", abs(wz_jecDn.GetBinContent(1) - wz.GetBinContent(1))
#    print "wz_jecUp 2", abs(wz_jecUp.GetBinContent(2) - wz.GetBinContent(2))
#    print "wz_jecDn 2", abs(wz_jecDn.GetBinContent(2) - wz.GetBinContent(2))
#    print "wz_jecUp 3", abs(wz_jecUp.GetBinContent(3) - wz.GetBinContent(3))
#    print "wz_jecDn 3", abs(wz_jecDn.GetBinContent(3) - wz.GetBinContent(3))
#    print "wz_jecUp 4", abs(wz_jecUp.GetBinContent(4) - wz.GetBinContent(4))
#    print "wz_jecDn 4", abs(wz_jecDn.GetBinContent(4) - wz.GetBinContent(4))
#    print "#############PDF##########################################"
#    print "zz_Pdf 1", pdfZZ.GetBinContent(1)
#    print "zz_Pdf 2", pdfZZ.GetBinContent(2)
#    print "zz_Pdf 3", pdfZZ.GetBinContent(3)
#    print "zz_Pdf 4", pdfZZ.GetBinContent(4)
#    print "wz_Pdf 1", pdfWZ.GetBinContent(1)
#    print "wz_Pdf 2", pdfWZ.GetBinContent(2)
#    print "wz_Pdf 3", pdfWZ.GetBinContent(3) 
#    print "wz_Pdf 4", pdfWZ.GetBinContent(4) 
#    print "#############STAT##########################################"
#    print "zz_stat 1", zz.GetBinError(1)/zz.GetBinContent(1)
#    print "zz_stat 2", zz.GetBinError(2)/zz.GetBinContent(2)
#    print "zz_stat 3", zz.GetBinError(3)/zz.GetBinContent(3)
#    print "zz_stat 4", zz.GetBinError(4)/zz.GetBinContent(4)
#    print "wz_stat 1", wz.GetBinError(1)/wz.GetBinContent(1)
#    print "wz_stat 2", wz.GetBinError(2)/wz.GetBinContent(2)
#    print "wz_stat 3", wz.GetBinError(3)/wz.GetBinContent(3)
#    print "wz_stat 4", wz.GetBinError(4)/wz.GetBinContent(4)
#    print "#############ZZ noKF - ZZ pt #############################################"
    zz_noKF.Scale(1./zz_noKF.Integral())
    zz_pt.Scale(1./zz_pt.Integral())
    zz_noKF.Scale(zz_noScale.Integral()) 
    zz_pt.Scale(zz_noScale.Integral()) 
    #print "zz kfactoris now 1", abs(zz_noKF.GetBinContent(1) - zz_pt.GetBinContent(1))
    #print "zz kfactoris now 2", abs(zz_noKF.GetBinContent(2) - zz_pt.GetBinContent(2))
    #print "zz kfactoris now 3", abs(zz_noKF.GetBinContent(3) - zz_pt.GetBinContent(3))
    #print "zz kfactoris now 4", abs(zz_noKF.GetBinContent(4) - zz_pt.GetBinContent(4))




    kfactorZZ = r.TH1F('kfactorZZ','kfactorZZ',4,0,4)
    kfactorZZ.SetBinContent(1, abs(zz_noKF.GetBinContent(1) - zz_pt.GetBinContent(1)) );
    kfactorZZ.SetBinContent(2, abs(zz_noKF.GetBinContent(2) - zz_pt.GetBinContent(2)) );
    kfactorZZ.SetBinContent(3, abs(zz_noKF.GetBinContent(3) - zz_pt.GetBinContent(3)) );
    kfactorZZ.SetBinContent(4, abs(zz_noKF.GetBinContent(4) - zz_pt.GetBinContent(4)) );

    lumiUnc = 0.025;lepUnc = 0.05;triggerUnc = 0.03;
    for i in [0, 1, 2, 3]:
        zz.SetBinError(i, math.sqrt((zz.GetBinContent(i)*math.sqrt(0.07**2 + lepUnc**2+triggerUnc**2))**2+zz.GetBinError(i)**2 + max(math.sqrt(abs(zz_qcdExclUp.GetBinContent(i)-zz_qcdExclNom.GetBinContent(i))**2+ abs(zz_qcdInclUp.GetBinContent(i)-zz_qcdInclNom.GetBinContent(i))**2),math.sqrt(abs(zz_qcdExclDn.GetBinContent(i)-zz_qcdExclNom.GetBinContent(i))**2+ abs(zz_qcdInclUp.GetBinContent(i)-zz_qcdInclNom.GetBinContent(i))**2))+max(abs(zz_jecUp.GetBinContent(i)-zz.GetBinContent(i))**2,abs(zz_jecDn.GetBinContent(i)-zz.GetBinContent(i))**2)+abs(zz_noKF.GetBinContent(i)-zz_pt.GetBinContent(i))**2 + (ZZ_met[i])**2  ))
        zz_syst.SetBinError(i, math.sqrt((zz.GetBinContent(i)*math.sqrt(0.07**2 + lepUnc**2+triggerUnc**2))**2+ max(math.sqrt(abs(zz_qcdExclUp.GetBinContent(i)-zz_qcdExclNom.GetBinContent(i))**2+ abs(zz_qcdInclUp.GetBinContent(i)-zz_qcdInclNom.GetBinContent(i))**2),math.sqrt(abs(zz_qcdExclDn.GetBinContent(i)-zz_qcdExclNom.GetBinContent(i))**2+ abs(zz_qcdInclUp.GetBinContent(i)-zz_qcdInclNom.GetBinContent(i))**2))+max(abs(zz_jecUp.GetBinContent(i)-zz.GetBinContent(i))**2,abs(zz_jecDn.GetBinContent(i)-zz.GetBinContent(i))**2)+abs(zz_noKF.GetBinContent(i)-zz_pt.GetBinContent(i))**2 + (ZZ_met[i])**2  ))
        wz.SetBinError(i, math.sqrt((wz.GetBinContent(i)*math.sqrt(0.06**2 +0.05**2+  lepUnc**2+triggerUnc**2))**2+wz.GetBinError(i)**2 + max(math.sqrt(abs(wz_qcdExclUp.GetBinContent(i)-wz_qcdExclNom.GetBinContent(i))**2+ abs(wz_qcdInclUp.GetBinContent(i)-wz_qcdInclNom.GetBinContent(i))**2),math.sqrt(abs(wz_qcdExclDn.GetBinContent(i)-wz_qcdExclNom.GetBinContent(i))**2+ abs(wz_qcdInclUp.GetBinContent(i)-wz_qcdInclNom.GetBinContent(i))**2))+max(abs(wz_jecUp.GetBinContent(i)-wz.GetBinContent(i))**2,abs(wz_jecDn.GetBinContent(i)-wz.GetBinContent(i))**2)+ (WZ_met[i])**2  ))
        wz_syst.SetBinError(i, math.sqrt((wz.GetBinContent(i)*math.sqrt(0.06**2 +0.05**2+ lepUnc**2+triggerUnc**2))**2+ max(math.sqrt(abs(wz_qcdExclUp.GetBinContent(i)-wz_qcdExclNom.GetBinContent(i))**2+ abs(wz_qcdInclUp.GetBinContent(i)-wz_qcdInclNom.GetBinContent(i))**2),math.sqrt(abs(wz_qcdExclDn.GetBinContent(i)-wz_qcdExclNom.GetBinContent(i))**2+ abs(wz_qcdInclUp.GetBinContent(i)-wz_qcdInclNom.GetBinContent(i))**2))+max(abs(wz_jecUp.GetBinContent(i)-wz.GetBinContent(i))**2,abs(wz_jecDn.GetBinContent(i)-wz.GetBinContent(i))**2)+ (WZ_met[i])**2  ))
        dy.SetBinError(i, math.sqrt(dy.GetBinError(i)**2+(0.5*dy.GetBinContent(i))**2  ))
        dy_syst.SetBinError(i, math.sqrt((0.5*dy.GetBinContent(i))**2  ))

    others = treeOTHERS.getTH1F(lint, var+"others"+scutstring, treevar, nbins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton,  cuts.slep0jet, cuts.SF]), '', xlabel, "1",kf)
    others_orig = treeOTHERS.getTH1F(lint, var+"others_orig"+scutstring, treevar, nbins, 1, 1, cuts.AddList([specialcut, cuts.goodLepton,  cuts.slep0jet, cuts.SF]), '', xlabel, "1",kf)
    for bin, label in scan.SRLabels.items():
        others.SetBinError(bin, math.sqrt((others.GetBinContent(bin)*0.5)**2 + others.GetBinError(bin)**2))
    
    # aesthetics
    print "dy before rinout ", dy.Integral()
    dy.Scale(0.138)
    print "dy after rinout ", dy.Integral()
    wz.Scale(1.06);wz_jecUp.Scale(1.06);wz_jecDn.Scale(1.06)
    zz.Scale(0.94);zz_jecUp.Scale(0.94);zz_jecDn.Scale(0.94)
    rare = copy.deepcopy(others)
    #rare.Add(dy, 1.);
    rare.SetFillColorAlpha(r.kViolet+3, 1);rare.SetTitle("Rare processes");rare.SetLineColor(r.kBlack)
    dy.SetFillColorAlpha(r.kYellow-9, 1);dy.SetTitle("DY+jets");dy.SetLineColor(r.kBlack)
    zz.SetFillColorAlpha(r.kCyan+2, 1);zz.SetTitle("ZZ#rightarrow 2l");zz.SetLineColor(r.kBlack)
    wz.SetFillColorAlpha(r.kGreen-8, 1);wz.SetTitle("WZ#rightarrow 3l");wz.SetLineColor(r.kBlack)
    prediction.SetFillColorAlpha(r.kRed-9, 0.7);prediction.SetTitle("Flavor symmetric");
    da_SF.SetTitle("data SF")       
    #mc_FS.SetTitle("FS MC");mc_FS.SetFillColorAlpha(r.kRed-9, 0.7)       

    mc_stack.Add(dy      );                      
    mc_stack.Add(rare      );                      
    mc_stack.Add(wz        );                     
    mc_stack.Add(zz        );                     
    #mc_stack.Add(mc_FS);                      
    mc_stack.Add(prediction);                      
    mc_stack.Draw();                               
    mc_stack.GetXaxis().SetTitle('p_{T}^{miss} [GeV]'); 
    #mc_full = copy.deepcopy(mc_FS);         
    mc_full = copy.deepcopy(prediction);         
    mc_full.Add(dy, 1.);                       
    mc_full.Add(rare, 1.);                       
    mc_full.Add(wz, 1.);                      
    mc_full.Add(zz, 1.);                     
    print "mc_full error", mc_full.GetBinError(4)
    mc_full_e = copy.deepcopy(mc_full);          
    mc_full_e.SetFillColorAlpha(r.kBlue+3, 1);mc_full_e.SetFillStyle(3017); mc_full_e.SetMarkerSize(0.);
    
    maxrat = 0.5
    for ib in range(1,da_SF.GetNbinsX()+1):
        tmp_rat = da_SF.GetBinContent(ib)/( mc_full.GetBinContent(ib) if mc_full.GetBinContent(ib) > 0 else 1. )
        if tmp_rat > maxrat:
            maxrat = tmp_rat                                                                            
        
    maxCont = max(da_SF.GetMaximum(), mc_full.GetMaximum())
    da_SF.GetYaxis().SetRangeUser(0.5, 1.50*maxCont)
    mc_stack.SetMaximum(1.5*maxCont)
    mc_stack.SetMinimum(0.5)
    SetOwnership(mc_stack, 0);SetOwnership(da_SF, 0);SetOwnership(mc_full, 0);            

    print helper.bcolors.HEADER + '[result scaled by RSFOF for DATA] ' + helper.bcolors.OKBLUE + 'Producing plot...' + helper.bcolors.ENDC
    plot_result = Canvas.Canvas('results/%s/plot_result_%s_daPreddaObs%s'%(newLumiString, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.6, 0.59, 0.90, 0.87)
    plot_result.addStack(mc_stack, "HIST" , 1, 1)
    plot_result.addHisto(mc_full_e, 'E2,SAME'  , '' , 'PL', r.kBlack , 1, -1)
    if scutstring == 'ee':string = 'ee'
    if scutstring == 'mm':string = '#mu#mu'
    if scutstring == 'SF':string = 'ee+#mu#mu'
    plot_result.addHisto(da_SF    , 'E1,SAME'  , 'Observed data '+string, 'PL', r.kBlack  , 1,  0)
    #plot_result.saveRatio(1, 1, 0, lint, mc_full, mc_full, 0. , int(maxrat+1.0)) 
    #THIS IS TO UNBLIND!
    plot_result.saveRatio(1, 1, 0, lint, da_SF, mc_full, 0. , int(maxrat+1.0))                                                                                                       
    makeResultsTable(da_SF, prediction, da_OF_orig, prediction_stat, prediction_syst, zz, zz_orig, zz_syst,  wz, wz_orig, wz_syst, rare,  others_orig, mc_full)
     
    plot_resultSigs = Canvas.Canvas('results/%s/plot_result_%s_daPreddaObs%s_withSignal'%(newLumiString, var, '' if not scutstring else '_'+scutstring), 'png,pdf', 0.6, 0.59, 0.90, 0.87)
    plot_resultSigs.addStack(mc_stack, "HIST" , 1, 1)
    plot_resultSigs.addHisto(mc_full_e, 'E2,SAME'  , '' , 'PL', r.kBlack , 1, -1)
    plot_resultSigs.addHisto(sig1, "HIST,SAME", "m_{l}: %s, m_{{\Chi_{1}^{0}}} : %s" %(signames[0][0], signames[0][1]), "L", r.kGreen-9, 1, 0)
    plot_resultSigs.addHisto(sig2, "HIST,SAME", "m_{l}: %s, m_{{\Chi_{1}^{0}}} : %s" %(signames[1][0], signames[1][1]), "L", r.kCyan-7, 1, 0)
    plot_resultSigs.addHisto(sig3, "HIST,SAME", "m_{l}: %s, m_{{\Chi_{1}^{0}}} : %s" %(signames[2][0], signames[2][1]), "L", r.kPink-4, 1, 0)
    plot_resultSigs.addHisto(sig4, "HIST,SAME", "m_{l}: %s, m_{{\Chi_{1}^{0}}} : %s" %(signames[3][0], signames[3][1]), "L", r.kOrange-4, 1, 0)
    plot_resultSigs.addHisto(da_SF, 'E1,SAME'  , 'Observed data'+string, 'PL', r.kBlack  , 1,  0)
    plot_resultSigs.saveRatio(1, 1, 0, lint, da_SF, mc_full, 0. , int(maxrat+1.0)) 
    makeResultsTableSig(da_SF, prediction, da_OF_orig, prediction_stat, prediction_syst,  zz, zz_orig, zz_syst, wz, wz_orig, wz_syst, rare, others_orig, dy, dy_orig, dy_syst, mc_full, sig1, sig2, sig3,sig4, signames)

    tfile = r.TFile('datacards/forDatacards_slepton_2017%s_withDY.root'%scutstring,'recreate')
    tfile.cd()
    tfile.WriteTObject(dy_stat      , 'dy_stat')
    tfile.WriteTObject(rare      , 'rare')
    tfile.WriteTObject(wz     , 'wz_SF')
    tfile.WriteTObject(zz     , 'zz_SF')
    tfile.WriteTObject(da_SF     , 'da_SF')
    tfile.WriteTObject(da_OF     , 'da_OF')
    tfile.WriteTObject(result[2] , 'tf_CR_SR')
    tfile.WriteTObject(zz_jecUp     , 'zz_jecUp')
    tfile.WriteTObject(zz_jecDn     , 'zz_jecDn')
    tfile.WriteTObject(wz_jecUp     , 'wz_jecUp')
    tfile.WriteTObject(wz_jecDn     , 'wz_jecDn')
    tfile.WriteTObject(zz_qcdInclNom    , 'zz_qcdInclNom')
    tfile.WriteTObject(zz_qcdInclUp     , 'zz_qcdInclUp')
    tfile.WriteTObject(zz_qcdInclDn     , 'zz_qcdInclDn')  
    tfile.WriteTObject(zz_qcdExclNom    , 'zz_qcdExclNom')
    tfile.WriteTObject(zz_qcdExclUp     , 'zz_qcdExclUp')
    tfile.WriteTObject(zz_qcdExclDn     , 'zz_qcdExclDn')  
    tfile.WriteTObject(wz_qcdInclNom    , 'wz_qcdInclNom')
    tfile.WriteTObject(wz_qcdInclUp     , 'wz_qcdInclUp')
    tfile.WriteTObject(wz_qcdInclDn     , 'wz_qcdInclDn')  
    tfile.WriteTObject(wz_qcdExclNom    , 'wz_qcdExclNom')
    tfile.WriteTObject(wz_qcdExclUp     , 'wz_qcdExclUp')
    tfile.WriteTObject(wz_qcdExclDn     , 'wz_qcdExclDn')  
    tfile.WriteTObject(pdfZZ        , 'pdfZZ')
    tfile.WriteTObject(pdfWZ        , 'pdfWZ')
    tfile.WriteTObject(kfactorZZ    , 'kfactorZZ')
    tfile.Close()

    if returnplot:
        return plot_result                                                                                                                                                                         
        del mc_stack;da_SF;mc_stack_perGeV;da_perGeV;mc_full; mc_full_perGeV;                                                                                                                                                                                                                                                                                   


if __name__ == '__main__':

    print 'Starting to produce some good ol\' results...'
    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples FilenameWithIngredients', version='%prog 1.0')
    parser.add_option('-t', action='store_true', dest='onlyTTbar', default=False, help='just do OF closure test')
    parser.add_option('-c', action='store_true', dest='onlyClosure', default=False, help='just do the closure test. don\'t bother with data')
    parser.add_option('-b', '--nbs', action='store', type=int, dest='nbs', default=0, help='do this for different numbers of b\'s')
    parser.add_option('-i', '--ingredients', action='store', type=str, dest='ingredientsFile', default='ingredients.dat', help='the ingredients file. default \'ingredients.dat\'')
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samplesSlepton.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-l', '--loadShapes', action='store_true', dest='loadShapes', default=False, help='reload dy shapes. default is off since this takes a while')
    parser.add_option('-M', '--maxRun', action='store', type=int, dest='maxRun', default=999999, help='max run to use for analysis (run is included)')
    parser.add_option('-m', '--minRun', action='store', type=int, dest='minRun', default=-1    , help='min run to use for analysis (run not included)')
    ## make the options globa.. also the lumi
    global opts, lumi, lumi_str, dy_shapes, nbstring
    global ingMC, ingDA, onZ, treeDA, treeMC, treeDY, treeTT, treeRA
    (opts, args) = parser.parse_args()

    ingredientsFile = opts.ingredientsFile
    print 'running with these options \n'
    for key, value in opts.__dict__.items():
        print '%-20s : %-20s' %(key, value)
    print ' \n\n'
    kf = "noKFactor"
    print 'Going to load DATA and MC trees...'
    dyDatasets = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO']
    zzDatasets = ['ZZTo2L2Nu', 'GluGluToContinToZZTo2mu2nu', 'GluGluToContinToZZTo2e2nu']
    wzDatasets = ['WZTo3LNu']
    othersDatasets = ['WWZ', 'WZZ', 'ZZZ', 'TWZ', 'tZq_ll', 'TTZToLLNuNu_ext2', 'TTZToQQ', 'TTLLJets_m1to10', 'TTHnobb_pow', 'VHToNonbb', 'GGHZZ4L', 'QQHZZ4L']
    fsDatasets = ['TTTT',  'TTTo2L2Nu', 'TBar_tch_powheg', 'T_tch_powheg', 'WWTo2L2Nu', 'WWW', 'WWG', 'WWDouble', 'WpWpJJ', 'TTWToLNu_ext2',  'TTWToQQ', 'TTJets_SingleLeptonFromTbar', 'TTJets_SingleLeptonFromT']                       

    mcDatasets = fsDatasets+dyDatasets + othersDatasets + zzDatasets + wzDatasets
 
                                                                                     
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
    treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasets, 'DY'), 'DY'  , 0, isScan = 0)
    treeFS = Sample.Tree(helper.selectSamples(opts.sampleFile, fsDatasets, 'FS'), 'FS'  , 0, isScan = 0)
    treeOTHERS = Sample.Tree(helper.selectSamples(opts.sampleFile, othersDatasets, 'OTHERS'), 'OTHERS'  , 0, isScan = 0)
    treeWZ = Sample.Tree(helper.selectSamples(opts.sampleFile, wzDatasets, 'WZ'), 'WZ'  , 0, isScan = 0)
    treeZZ = Sample.Tree(helper.selectSamples(opts.sampleFile, zzDatasets, 'ZZ'), 'ZZ'  , 0, isScan = 0)
    treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1, isScan = 0)
    
    #signals = ['SMS_450_40', 'SMS_375_160', 'SMS_250_180', 'SMS_100_1']
    #signames = [[450, 40], [375, 160], [250, 180], [100, 1]]
    signals = ['SMS_450_20', 'SMS_400_20', 'SMS_350_20', 'SMS_100_1']
    signames = [[450, 20], [400, 20], [350, 20], [100, 1]]    
    xsecf = open('datacards/xsec_SUM_13tev_fit.txt', 'r')
    #xsecf = open('datacards/xsec_lLlL_13tev_fit.txt', 'r')
    xsecs = eval(xsecf.read())
    xsecf.close()                                             
    treeSIG1 = Sample.Tree(helper.selectSamples("samplesSlepton.dat", [signals[0]], 'SI'), 'SI', 0, isScan =  xsecs[signames[0][0]][0])
    treeSIG2 = Sample.Tree(helper.selectSamples("samplesSlepton.dat", [signals[1]], 'SI'), 'SI', 0, isScan =  xsecs[signames[1][0]][0])
    treeSIG3 = Sample.Tree(helper.selectSamples("samplesSlepton.dat", [signals[2]], 'SI'), 'SI', 0, isScan =  xsecs[signames[2][0]][0])
    treeSIG4 = Sample.Tree(helper.selectSamples("samplesSlepton.dat", [signals[3]], 'SI'), 'SI', 0, isScan =  xsecs[signames[3][0]][0])
    sigs = [treeSIG1, treeSIG2, treeSIG3, treeSIG4]                                                               

    print 'Trees successfully loaded...'

    isBlinded = True

    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()
    lint = 35.9  ; maxrun = 999999; lint_str = '35.9invfb'
    print 'Running with an integrated luminosity of %.2f fb-1' %(lint)

    ## ============================================================
    ## ========== set RSFOF globally ==============================
    ## ============================================================
    #rsfof, rsfof_e, rsfof_mc, rsfof_mc_e, rmue_a_da, rmue_a_mc, rmue_b_da, rmue_b_mc =  makeFactorsTable()
    ## result plots in different variables:
    for v in ['met']:
        makeResultData('slepton2017', v,signames, maxrun,lint,specialcut='' , scutstring = 'SF',    _options='returnplot,splitFlavor')
        #makeResultData('slepton2017ee', v,signames, maxrun,lint,specialcut='(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -121)' , scutstring = 'ee',    _options='returnplot,splitFlavor')
        #makeResultData('slepton2017mm', v,signames, maxrun,lint,specialcut='(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -169)' , scutstring = 'mm',    _options='returnplot,splitFlavor')
        #makeClosureTestPlots('slepton2017','met0jet','nJet25_Edge  == 0 && mt2_Edge < 90 && mt2_Edge > 60' , '60-90', True)
        #makeClosureTestPlots('slepton2017','met0jet','nJet25_Edge  == 0 && mt2_Edge > 90' , 'ht0', True)
        #makeClosureTestPlots('slepton2017','met0jet','nJet25_Edge  == 0 && mt2_Edge > 90' , '90', True)
        #makeClosureTestPlots('slepton2017','met0jet','nJet25_Edge  == 0 && mt2_Edge > 80' , '80', True)
        #makeClosureTestPlots('slepton2017','met0jet','nJet25_Edge  == 0 && mt2_Edge > 50' , '50', True)
        #makeClosureTestPlots('slepton2017','met0jet','nJet25_Edge  == 0 && mt2_Edge > 60' , '60', True)
        #makeClosureTestPlots('slepton2017','met0jet','nJet25_Edge  == 0 && mt2_Edge > 70' , '70', True)
        #makeClosureTestPlots('slepton2017',v,'', 'inclusive', True)
        #resultPlotLoNll = makeResultData('Edge_Moriond2017', v,maxrun,lint,specialcut='nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) < 21.' , scutstring = 'nllBelow21',    _options='returnplot,splitFlavor')
        #resultPlotHiNll = makeResultData('Edge_Moriond2017', v,maxrun,lint,specialcut='nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) >= 21.' , scutstring = 'nllAbove21',    _options='returnplot,splitFlavor')

   # makeRSOFTable('Edge_Moriond2017')
    
