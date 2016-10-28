import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats
import math, sys, optparse, copy, re, array, os, pickle


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Tables     as Tables
import include.Scans      as Scans

from multiprocessing import Pool

global _r
_r = r.TRandom3(42)

def smearNumber(num):
    rand = _r.Gaus(0,0.1)
    smeared = num*(1.+rand)
    print 'smeared', num, 'to', smeared
    return smeared

def getSRIDs(region):
    # region[0] is mll binning
    # region[1] is nll binning
    # region[2] is mt2 binning

    if region[0]   == 'belowZ'   : mllids = [1]
    elif region[0] == 'onZ'      : mllids = [2]
    elif region[0] == 'aboveZ'   : mllids = [3]
    elif region[0] == 'highMass' : mllids = [4]
    elif region[0] == 'vHighMass': mllids = [5]
    elif region[0] == 'highinc'  : mllids = [3,4,5] # for only three mll bins

    if   region[1]   == 'lowNll'   : nllids = [0]
    elif region[1]   == 'highNll'  : nllids = [1]
    
    if   region[2]   == 'lowMT2'   : mt2ids = [0]
    elif region[2]   == 'highMT2'  : mt2ids = [1]
    elif region[2]   == 'incMT2'   : mt2ids = [0,1]
    elif region[2]   == 'lowHT'    : mt2ids = [0]
    elif region[2]   == 'medHT'    : mt2ids = [1]
    elif region[2]   == 'highHT'   : mt2ids = [2]
    allSRs = []
    for mllid in mllids:
        for nllid in nllids:
            for mt2id in mt2ids:
                allSRs.append(100*mt2id + 10 *nllid + mllid)


    # print 'for region', '_'.join(region), 'returning srIDs', allSRs

    return allSRs
    
def bins(var):
    if var =='met': return [ 30,100,400]
    if var =='zpt': return [ 50,  0,250]
    if var =='maxjj': return [ 50,  0,500]
    if var =='minjj': return [ 50,  0,500]
    if var =='bestjj': return [ 50, 0,500]

def makePlots(var):
    dyDatasets = ['DYJetsToLL_M50_HT100to200',
                  'DYJetsToLL_M50_HT200to400',
                  'DYJetsToLL_M50_HT400to600',
                  'DYJetsToLL_M50_HT600toInf']
    samples = []
    samples.append( Sample.Tree(helper.selectSamples(opts.sampleFile, ['TTLep_pow'] , 'TT'), 'TT', 0) )
    samples.append( Sample.Tree(helper.selectSamples(opts.sampleFile, ['WZTo3L1Nu'     ] , 'WZ'), 'WZ', 0) )
    samples.append( Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasets    , 'DY'), 'DY', 0) )
    if var =='met': v = 'met_Edge'      
    if var =='zpt': v = 'lepsZPt_Edge'
    if var =='maxjj': v = 'maxMjj_Edge'
    if var =='minjj': v = 'minMjj_Edge'
    if var =='bestjj': v = 'bestMjj_Edge'
    b = bins(var)

    histos = []
    for sample in samples:
        histos.append(sample.getTH1F(lumi, '%s'%(sample.name), v, b[0], b[1], b[2], scan.cuts_norm, '', var) )
    histos.append( scan.tree.getTH1F(lumi, 'CN_350_20' , v, b[0], b[1], b[2], '('+scan.cuts_norm+'&& GenSusyMNeutralino2_Edge == 350 && GenSusyMNeutralino_Edge ==  20)', '', var) )
    histos.append( scan.tree.getTH1F(lumi, 'CN_350_100', v, b[0], b[1], b[2], '('+scan.cuts_norm+'&& GenSusyMNeutralino2_Edge == 350 && GenSusyMNeutralino_Edge == 100)', '', var) )

    plot = Canvas.Canvas('TChiNeuWZ/plot_%s'%var, 'png,pdf', 0.6, 0.65, 0.85, 0.82)
    lego = 0
    for i,h in enumerate(histos):
        h.SetMarkerStyle(20)
        h.SetMarkerSize(0.8*h.GetMarkerSize())
        h.SetMarkerColor(r.kBlack+i)
        h.Scale(1./h.Integral())
        plot.addHisto(h, 'p %s'%('' if not i else 'same')    , h.GetName(), 'PL' , r.kBlack+i, 1,  lego)
        lego+=1
    plot.save(1, 0, 0 , lumi)

def makeMCDatacards():
    print 'producing mc datacards'
    ttDatasets = ['TT_pow_ext34']
    dyDatasets = ['DYJetsToLL_M50_HT100to200_ext', 'DYJetsToLL_M50_HT200to400_ext',
                  'DYJetsToLL_M50_HT400to600_ext', 'DYJetsToLL_M50_HT600toInf_ext']
    print 'getting trees'
    treeTT = Sample.Tree(helper.selectSamples(opts.sampleFile, ttDatasets, 'TT'), 'TT', 0, isScan = 0)
    treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasets, 'DY'), 'DY', 0, isScan = 0)
    print 'getting yields'
    print cuts.AddList([scan.cuts_norm,cuts.goodLepton])
    tt = treeTT.getTH1F(lumi, 'tt_srs', scan.srID, 300, 0, 300, cuts.AddList([scan.cuts_norm,cuts.goodLepton]), '', '')
    c = r.TCanvas()
    tt.Draw()
    c.SaveAs("srs.pdf")
    dy = treeDY.getTH1F(lumi, 'dy_srs', scan.srID, 300, 0, 300, cuts.AddList([scan.cuts_norm,cuts.goodLepton]), '', '')
    print 'doing datacards'
    for region in scan.regions:
        totalTT = sum( tt.GetBinContent(tt.FindBin(i)) for i in getSRIDs(region) )
        totalDY = sum( dy.GetBinContent(dy.FindBin(i)) for i in getSRIDs(region) )
        mll = region[0]; mass = region[1]; nb = region[2]
        bin_name = '%s_%s_%s'%(mll, mass, nb)
        obs = int(totalTT+totalDY)
        fs  = totalTT
        fs_unc = 0.05
        onz = totalDY
        onz_e = 0.3*onz
        of_yield = int(totalTT)
        rsfof = 1.
        
        dc = '''# this is the datacard for bin {bin_name}
imax 1  number of channels
jmax 2  number of backgrounds
kmax *  number of nuisance parameters (sources of systematical uncertainties)
------------
# only one channel here with observation {obs}
bin            {bin_name}
observation    {obs}
------------
bin        {bin_name}     {bin_name}     {bin_name}
process    sig            FS             DY    
process    0              1              2
rate       XXRATEXX         {fs_bkg:.2f}       {dy_bkg:.2f}
------------
deltaS       lnN              1.20      -           -       
lumi         lnN              1.12      -           -       
FS_unc       lnN              -         {fs_unc:.2f}    -       
{bin_name}_fs_stat      gmN {of_yield}   -         {rsfof:.3f}        -       
DY_unc       lnN              -         -           {dy_unc:.2f}'''.format(bin_name=bin_name,
                                                                           obs=int(obs),#
                                                                           fs_bkg=fs,#
                                                                           fs_unc=1+fs_unc,#
                                                                           dy_bkg=max(0.01,onz),#
                                                                           dy_unc=1+onz_e/abs(max(0.01,onz)),
                                                                           of_yield=int(of_yield), #
                                                                           rsfof=rsfof)
        tmp_file = open('datacards/datacards_%s/mc_%s.txt'%(scan.name,bin_name),'w')
        tmp_file.write(dc)
        tmp_file.close()

# def makeMCDatacards():
#     daDatasets = ['DoubleMuon_Run2015C_25ns-05Oct_v1_runs_246908_260628' , 'DoubleEG_Run2015C_25ns-05Oct_v1_runs_246908_260628',
#                   'DoubleMuon_Run2015D-05Oct_v1_runs_246908_260628'      , 'DoubleEG_Run2015D-05Oct_v1_runs_246908_260628'     ,
#                   'DoubleMuon_Run2015D_v4_runs_246908_260628'            , 'DoubleEG_Run2015D_v4_runs_246908_260628'           ]

#     dyDatasets = ['DYJetsT0LL_M50']

#     dyDatasets = ['DYJetsToLL_M50_HT100to200',
#                   'DYJetsToLL_M50_HT200to400',
#                   'DYJetsToLL_M50_HT400to600',
#                   'DYJetsToLL_M50_HT600toInf']

#     global dic
#     samples = []; dic = {}
#     samples.append( Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets    , 'DA'), 'DA', 1) )
#     samples.append( Sample.Tree(helper.selectSamples(opts.sampleFile, ['TTLep_pow'] , 'TT'), 'TT', 0) )
#     samples.append( Sample.Tree(helper.selectSamples(opts.sampleFile, ['WZTo3L1Nu'     ] , 'WZ'), 'WZ', 0) )
#     samples.append( Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasets    , 'DY'), 'DY', 0) )

#     global nllDistributions, cumDistributions
#     nllDistributions = []
#     cumDistributions = []
#     for sample in samples:
#         all_cuts = scan.cuts_norm
#         #all_cuts = cuts.AddList([scan.cuts_norm, cuts.twoLeptons, cuts.trigger])
#         dic[sample.name] = sample.getTH1F(lumi, '%s_yields'%sample.name, 'srID_Edge', 200, 100, 300, all_cuts, '', 'SR ID')
#         nllDistributions.append( sample.getTH1F(lumi, '%s_nll'%sample.name, '-1.*TMath::Log(lh_ana_met_data_Edge*lh_ana_mlb_data_Edge*lh_ana_a3d_data_Edge*lh_ana_zpt_data_Edge)', 80, 12, 30, all_cuts, '', 'NLL') )

#     nllDistributions.append( scan.tree.getTH1F(lumi, 'signal_nll_350/100', '-1.*TMath::Log(lh_ana_met_data_Edge*lh_ana_mlb_data_Edge*lh_ana_a3d_data_Edge*lh_ana_zpt_data_Edge)', 80, 12, 30, '('+scan.cuts_norm+'&& GenSusyMNeutralino2_Edge == 350 && GenSusyMNeutralino_Edge == 100)', '', 'NLL') )
#     nllDistributions.append( scan.tree.getTH1F(lumi, 'signal_nll_350/20', '-1.*TMath::Log(lh_ana_met_data_Edge*lh_ana_mlb_data_Edge*lh_ana_a3d_data_Edge*lh_ana_zpt_data_Edge)', 80, 12, 30, '('+scan.cuts_norm+'&& GenSusyMNeutralino2_Edge == 350 && GenSusyMNeutralino_Edge == 20)', '', 'NLL') )
#     for i in nllDistributions:
#         #i.GetYaxis().SetRangeUser(0. )
#         #i.Scale(1./i.Integral())
#         cumDistributions.append( i.GetCumulative() )
    
#     plot = Canvas.Canvas('%s/plot_nll'%scan.name, 'png,pdf', 0.6, 0.65, 0.85, 0.82)
#     lego = 0
#     for i,h in enumerate(nllDistributions):
#         if 'DA' in h.GetName(): continue
#         if 'WZ' in h.GetName(): continue
#         if 'DY' in h.GetName(): continue
#         h.SetMarkerStyle(20)
#         h.SetMarkerSize(0.8*h.GetMarkerSize())
#         h.SetMarkerColor(r.kBlack+i)
#         plot.addHisto(h, 'p %s'%('' if not i else 'same')    , h.GetName(), 'PL' , r.kBlack+i, 1,  lego)
#         lego+=1
#     plot.save(1, 0, 0 , lumi)
#     ##
#     plot = Canvas.Canvas('%s/plot_nll_cum'%scan.name, 'png,pdf', 0.6, 0.65, 0.85, 0.82)
#     lego = 0
#     for i,h in enumerate(cumDistributions):
#         if 'DA' in h.GetName(): continue
#         if 'WZ' in h.GetName(): continue
#         if 'DY' in h.GetName(): continue
#         h.SetMarkerStyle(20)
#         h.SetMarkerSize(0.8*h.GetMarkerSize())
#         h.SetMarkerColor(r.kBlack+i)
#         plot.addHisto(h, 'p %s'%('' if not i else 'same')    , h.GetName(), 'PL' , r.kBlack+i, 1,  lego)
#         lego+=1
#     plot.save(1, 0, 0 , lumi)
#     print adfsdafsd

#     for region in scan.regions:
#         dayield = (lumi/2.3)*sum( dic['DA'].GetBinContent(dic['DA'].FindBin(i)) for i in getSRIDs(region) )
#         bkyields = []; bkstats = []; systs = []
#         for name, bkg in dic.items():
#             if not name  in ['DA','DY']:
#                 bkyields.append(sum( dic[name].GetBinContent(dic[name].FindBin(i))    for i in getSRIDs(region) ) )
#                 bkstats .append(math.sqrt(sum( dic[name].GetBinError(dic[name].FindBin(i))**2 for i in getSRIDs(region) ) ) )
#                 systs.append(1.15 if name == 'TT' else 1.2 if name == 'WZ' else 1.9)
#         tmp_dy = sum(dic['DY'].GetBinContent(dic['DY'].FindBin(i)) for i in getSRIDs(region))
#         print 'sum bkg', sum(bkyields)
#         print 'dy', tmp_dy
#         ## dyyield = math.sqrt( abs((dayield-sum(bkyields))) * tmp_dy )
#         dyyield = tmp_dy
#         bkyields.append(dyyield); bkstats.append(math.sqrt(dyyield));
#         systs.append(1.5)
#         tmp_dcf = open('datacards/empty_datacard.txt','r')
#         tmp_dcc = tmp_dcf.read()
#         nbkgs = len(list(d for d in dic.keys() if d != 'DA'))
#         tmp_dcc = tmp_dcc.replace('BINNAME'  , '_'.join(region) )
#         tmp_dcc = tmp_dcc.replace('BINS'     , ' '.join(['_'.join(region) for i in range(nbkgs)]) )
#         tmp_dcc = tmp_dcc.replace('BKGNAMES' , ' '.join(list(d for d in dic.keys() if d != 'DA')) )
#         tmp_dcc = tmp_dcc.replace('BKGNUMS'  , ' '.join(str(i) for i in range(1,nbkgs+1) )        )
#         tmp_dcc = tmp_dcc.replace('BKGYIELDS', ' '.join('%.2f'%i for i in bkyields )              )
#         tmp_dcc = tmp_dcc.replace('SYSDELTAS', ' '.join('-'    for i in range(nbkgs) )            )
#         tmp_dcc = tmp_dcc.replace('SYSLUMI'  , ' '.join(' -  ' for i in range(nbkgs) )            )
#         for i in range(nbkgs):
#             tmp_dcc+= 'systBKG{n} lnN - {before} {foo} {after} \n'.format(n=i+1, 
#                                                                           foo=systs[i], 
#                                                                           before=' '.join('-' for i in range(i)), 
#                                                                           after =' '.join('-' for i in range(nbkgs - i - 1)))
#         tmp_dcc = tmp_dcc.replace('SYSTS'    , ' '.join(str(i) for i in systs )                   )
#         #tmp_dcc = tmp_dcc.replace('SYSTATS'  , ' '.join('%.2f'%(1.+bkstats[i]/bkyields[i]) for i in range(nbkgs) ) )
#         tmp_dcc = tmp_dcc.replace('DATAOBS'  , str(int(dayield)))
#         tmp_dcc = tmp_dcc.replace('NBKG'     , str(nbkgs)       )
#         tmp_dc = open('datacards/datacards_{name}/mc_{reg}.txt'.format(name=scan.name, reg = '_'.join(region)), 'w')
#         tmp_dc.write(tmp_dcc)
#         tmp_dc.close()

    
#     return dic
    

def runCmd(cmd):
    os.system(cmd[0])
    os.chdir (cmd[1])
    os.system(cmd[2])
    return True


def produceLimits(whichb, njobs):
    basedir = os.path.abspath('.')+'/datacards/datacards_{name}/{name}/'.format(name=scan.name)
    subdirs = [i for i in os.listdir(basedir) if os.path.isdir(basedir+i)]
    pool = Pool(njobs)
    tasks = []
    for ind,d in enumerate(subdirs):
        mass = ''.join(i for i in d.split('_') if i.isdigit() )
        #if ind > 10: continue
        fd = basedir+d
        bstring = ''; b= ''
#        for i in whichb:
        bstring += ' {fd}/{d}*.txt'.format(fd=fd,b=i,d=d)
#            b+=i
        dc_name    = 'datacard_{massp}_full.txt'.format(massp=d)
        final_dc   = '{fd}/{dc_name}'.format(fd=fd,dc_name=dc_name)
        runcmd     = 'combineCards.py -S {bstr} > {final_dc}'.format(bstr=bstring,final_dc=final_dc)
        combinecmd = 'combine -m {mass} -M Asymptotic {dc_name}'.format(mass=mass,dc_name=dc_name)
        tasks.append([runcmd, fd, combinecmd])
    pool.map(runCmd, tasks)
    haddcmd = 'hadd -f {bd}/{name}_allLimits.root {bd}/*/higgs*.root'.format(name=scan.name,bd=basedir)
    #print haddcmd
    os.system(haddcmd)

def fillAndSaveDatacards(nbs):
    for region in scan.regions:
        mll = region[0]; mass = region[1]; nb = region[2]
        if not scan.makeMCDatacards:
            tmp_file = open('datacards/datacards_%s/%s_%s_%s.txt'%(scan.name,mll, mass, nb) ,'r')
        else:
            tmp_file = open('datacards/datacards_%s/mc_%s_%s_%s.txt'%(scan.name,mll, mass, nb) ,'r')
        tmp_dc  = tmp_file.read()
        tmp_histo = getattr(scan, 'eff_%s_%s_%s'%(mll, mass, nb))
        for i in range(1, tmp_histo.GetXaxis().GetNbins()+1):
            for j in range(1, tmp_histo.GetYaxis().GetNbins()+1):
                xmass = tmp_histo.GetXaxis().GetBinCenter(i)
                ymass = tmp_histo.GetYaxis().GetBinCenter(j)
                if ymass > xmass: continue
                if tmp_histo.GetBinContent(tmp_histo.GetBin(i,j)) == 0.:
                    continue
                mass_string = 'mSbottom_%.0f_mchi2_%.0f'%(xmass, ymass)
                helper.ensureDirectory('datacards/datacards_{name}/{name}'.format(name=scan.name))
                helper.ensureDirectory('datacards/datacards_{name}/{name}/{mass}'.format(name=scan.name,mass=mass_string))
                tmp_rate = scan.br*lumi*scan.xsecs[xmass][0]*tmp_histo.GetBinContent(tmp_histo.GetBin(i,j)) ## LUMI WEIGHTING HAPPENS HERE!!!!
                tmp_out = tmp_dc.replace('XXRATEXX', '%.3f'%(tmp_rate))
                print 'this is the fn', tmp_file.name.split('/')[-1]
                tmp_new = open('datacards/datacards_{name}/{name}/{mass}/{mass}_{fn}'.format(name=scan.name,mass=mass_string,fn=tmp_file.name.split('/')[-1]), 'w')
                tmp_new.write(tmp_out)
                tmp_new.close()
        tmp_file.close()

def adaptBinning(target, current):
    final = copy.deepcopy(target)
    final.Reset()
    final = final.Project3D('yx')
    current_yx = current.Project3D('yx')
    ret_histo = copy.deepcopy(target)
    ret_histo.Reset()
    for i in range(1,final.GetXaxis().GetNbins()+1):
        for j in range(1, final.GetYaxis().GetNbins()+1):
            xval = final.GetXaxis().GetBinCenter(i)
            yval = final.GetYaxis().GetBinCenter(j)
            cont = current_yx.GetBinContent(current_yx.FindBin(xval, yval))
            err  = current_yx.GetBinError  (current_yx.FindBin(xval, yval))
            final.SetBinContent(final.FindBin(xval, yval), cont)
            final.SetBinError  (final.FindBin(xval, yval), err )
            for k in range(1,ret_histo.GetZaxis().GetNbins()+1):
                ret_histo.SetBinContent(ret_histo.FindBin(xval, yval, ret_histo.GetZaxis().GetBinCenter(k)), cont)
                ret_histo.SetBinError  (ret_histo.FindBin(xval, yval, ret_histo.GetZaxis().GetBinCenter(k)), err )
    return final, ret_histo

def getSRYield(mll, nll, mt2):
    if mll   == 'belowZ'   : mllids = [1]
    elif mll == 'onZ'      : mllids = [2]
    elif mll == 'aboveZ'   : mllids = [3]
    elif mll == 'highMass' : mllids = [4]
    elif mll == 'vHighMass': mllids = [5]
    elif mll == 'highinc'  : mllids = [3,4,5] # for only three mll bins

    if   nll   == 'lowNll'   : nllids = [0]
    elif nll   == 'highNll'  : nllids = [1]
    
    if   mt2   == 'lowMT2'   : mt2ids = [0]
    elif mt2   == 'highMT2'  : mt2ids = [1]
    elif mt2   == 'incMT2'   : mt2ids = [0,1]
    elif mt2   == 'lowHT'    : mt2ids = [0]
    elif mt2   == 'medHT'    : mt2ids = [1]
    elif mt2   == 'highHT'   : mt2ids = [2]

    allSRs = []
    for mllid in mllids:
        for nllid in nllids:
            for mt2id in mt2ids:
                allSRs.append(100*mt2id + 10 *nllid + mllid)

    print 'at signal region',mll, nll, mt2
    print 'combining signal regions:',allSRs

    scan_sr_yield = scan.ngen_2d.Clone('yield_%s_%s_%s'%(mll, nll, mt2))
    scan_sr_yield.SetTitle('yield_%s_%s_%s'%(mll, nll, mt2))
    scan_sr_yield.SetName ('yield_%s_%s_%s'%(mll, nll, mt2))
    scan_sr_yield.Reset()
    for x in range(1, scan.norm.GetNbinsX()+1):
        for y in range(1, scan.norm.GetNbinsY()+1):
            if scan.norm.GetYaxis().GetBinCenter(y) > scan.norm.GetXaxis().GetBinCenter(x): continue
            tmp_yield  = 0.
            tmp_yield2 = 0.
            for z in allSRs:
                zbin = scan.norm.GetZaxis().FindBin(z)
                tmp_yield  +=  scan.norm.GetBinContent(scan.norm.GetBin(x, y, zbin))
                tmp_yield2 += (scan.norm.GetBinError  (scan.norm.GetBin(x, y, zbin)))**2
            scan_sr_yield.SetBinContent(scan_sr_yield.GetBin(x,y), tmp_yield)
            scan_sr_yield.SetBinError  (scan_sr_yield.GetBin(x,y), math.sqrt(tmp_yield2))
    scan_sr_eff = scan_sr_yield.Clone('eff_%s_%s_%s'%(mll, nll, mt2))
    scan_sr_eff.SetTitle('eff_%s_%s_%s'%(mll, nll, mt2))
    scan_sr_eff.GetZaxis().SetRangeUser(0., scan.zmaxEff)
    scan_sr_eff.Divide(scan.ngen_2d)
    setattr(scan, 'yield_%s_%s_%s'%(mll, nll, mt2), scan_sr_yield)
    setattr(scan, 'eff_%s_%s_%s'  %(mll, nll, mt2), scan_sr_eff  )

 
if __name__ == "__main__":

    r.gStyle.SetPaintTextFormat(".2f")
    r.gStyle.SetOptStat(0)

    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    parser.add_option('-r', '--reloadScan' , action='store_true', dest='reloadScan', help='reload scan. default %default')
    parser.add_option('-l', '--reloadLimits' , action='store_true', dest='reloadLimits', help='reload limits. default %default')
    parser.add_option('-n', '--scanName'     , action='store', type=str, dest='scanName', default='T6bbslepton', help='scan name. e.g. T6bbslepton (this is the default)')
    (opts, args) = parser.parse_args()

    print 'running with these options'
    for key, value in opts.__dict__.items():
        print '%-20s : %-20s' %(key, value)
    print ' \n\n'
    print 'Going to load the tree(s)...'

    global lumi, scan

    ##ttDatasets = ['TTLep_pow']
    ##treeTT = Sample.Tree(helper.selectSamples(opts.sampleFile, ttDatasets, 'TT'), 'TT'  , 0, isScan = False)
    cuts = CutManager.CutManager()
    lumi = 10.0

    ## have to think a way of reweighting the events with trigger and lepton SFs.
    ## weighting now done with isScan=True flag in samples.py

    ## this should be then sr-ID:m_slepton:m_sbottom for the final scan
    #cuts_norm = cuts.AddList([cuts.METJetsSignalRegion, cuts.GoodLeptonSFNoTrigger()]) ## trigger not available in fastsim
    #cuts_norm = cuts_norm.replace(cuts.twoLeptons, 't.nPairLep_Edge > 0') ## remove the filters, ugly, but it's a bit intricate in the samples

    ## ==============================================
    ## == try loading stuff from the pickled file ===
    ## ==============================================
    
    pickfile = 'datacards/datacards_{name}/{name}/{name}.pkl'.format(name=opts.scanName)
    if os.path.isfile(pickfile) and not opts.reloadScan:
        print 'getting the scan object from the pickled file'
        scan = pickle.load(open(pickfile,'r'))
    else:
        ## everything that takes long should be done here!
        scan = Scans.Scan(opts.scanName)
        scan.tree = Sample.Tree(helper.selectSamples(opts.sampleFile, scan.datasets, 'SIG'), 'SIG'  , 0, isScan = True)
        scan.norm = scan.tree.getTH3F(1., 'nPass_norm', scan.srID+':'+scan.yvar+':'+scan.xvar, scan.xbins.n+1, scan.xbins._min-scan.xbins.w/2., scan.xbins._max+scan.xbins.w/2.,  ## lumi set later for scans!!
                                      scan.ybins.n+1, scan.ybins._min-scan.ybins.w/2., scan.ybins._max+scan.ybins.w/2., 
                                      300, 0, 300, scan.cuts_norm, '', scan.xtitle, scan.ytitle, scan.ztitle)
        
        #for i in ['met', 'zpt', 'maxjj', 'minjj', 'bestjj']: makePlots(i)
        #print asfsdfs
        if scan.makeMCDatacards:
            print 'preparing datacards from MC'
            a = makeMCDatacards()
        #print asfsdfs

        print '=================================================='
        print '=================================================='
        print '===== this is the cut ============================\n', scan.cuts_norm
        print '=================================================='
        print '=================================================='

        ## we also need the number of generated events here!!
        if scan.has3DGen:
            scan.ngen = scan.tree.blocks[0].samples[0].smsCount ## take the first slice's ngen histo
            for ind,i in enumerate(scan.tree.blocks[0].samples):
                if ind: ## do not add the first one twice
                    scan.ngen.Add(i.smsCount, 1.) ## add all others
        else:
            scan.make3DGen()

        ## =====================================================
        ## == if anything takes a long time, load it with pickle
        ## =====================================================

        print 'this is the type of scan.ngen before', type(scan.ngen)
        newbinning   = adaptBinning(scan.norm, scan.ngen)
        scan.ngen_2d = newbinning[0]
        scan.ngen_3d = newbinning[1] ## this one has ngen in every single bin. for every SR. and it's 3D, so that's cool
        print 'this is the type of scan.ngen after', type(scan.ngen)

        tmp_ngen = copy.deepcopy(scan.ngen)

        #print adfasdf

        scan.eff_norm = scan.norm.Clone('efficiency_norm')
        scan.eff_norm.Divide(scan.ngen_3d)

        scan.xy = scan.eff_norm.Project3D('xy') # this means x versus y. so x is on the y-axis
        scan.xz = scan.eff_norm.Project3D('xz')
        scan.yx = scan.eff_norm.Project3D('yx') # that's the inclusive efficiency map
        scan.yz = scan.eff_norm.Project3D('yz')
        scan.zx = scan.eff_norm.Project3D('zx')
        scan.zy = scan.eff_norm.Project3D('zy')


        ## here we get the 2D efficiencies for all signal regions
        for region in scan.regions:
            mll = region[0]; nll = region[1]; mt2 = region[2]
            getSRYield(mll, nll, mt2)

        ## set the scan tree to 0 before pickling. it's huge (that's what she said)
        scan.tree = 0
        scan.ngen = tmp_ngen

    # now we take the default datacards and save them for each point into a subdirectory of the scan

    if opts.reloadLimits:
        print 'reloading limits and datacards'
        dobs = ['0b']#,'1b']
        fillAndSaveDatacards(dobs)
        produceLimits(dobs,5)

    os.system('mkdir -p mkdir -p makeExclusionPlot/config/%s/'%scan.paper)
    scan.makeExclusion()
    scan.makePrettyPlots()

    ## save the scan object in a pickle file to save time the second time around.
    pickle.dump(scan, open('datacards/datacards_{name}/{name}/{name}.pkl'.format(name=scan.name),'w') )
    print 'marc is stupid'
