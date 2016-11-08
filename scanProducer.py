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
    tt = treeTT.getTH1F(lumi, 'FSbkg', scan.srID, scan.srIDMax+1, -0.5, scan.srIDMax+0.5, cuts.AddList([scan.cuts_norm,cuts.goodLepton]), '', ''); tt.SetName('FSbkg')
    dy = treeDY.getTH1F(lumi, 'otherbkg', scan.srID, scan.srIDMax+1, -0.5, scan.srIDMax+0.5, cuts.AddList([scan.cuts_norm,cuts.goodLepton]), '', ''); dy.SetName('otherbkg')
    data = tt.Clone('data_obs'); data.Add(dy)
    helper.ensureDirectory('datacards/datacards_%s/%s'%(scan.name,scan.name))
    for SR, label in scan.shortLabels.items():
        datacard = '''imax 1 number of bins
jmax 2 number of processes minus 1
kmax *  number of nuisance parameters
----------------------------------------------------------------------------------------------------------------------------------
bin          {label}
observation  {obs}   
----------------------------------------------------------------------------------------------------------------------------------
bin          {label}      {label}        {label}
process      XXXSIGNAL     FSbkg          otherbkg
process      0             1              2   
rate         XXXSIGRATE    {fs}           {other}
----------------------------------------------------------------------------------------------------------------------------------
fs_stat_{label} gmN  {fs_int}  -            1.0             - 
fs_unc lnN             -            1.05            - 
jec   lnN           XXjecXX            -               -
El    lnN           XXElXX            -               -
Mu    lnN           XXMuXX            -               -
bHe   lnN           XXbHeXX            -               -
bLi   lnN           XXbLiXX            -               -
signalMCstats_{label} lnN    XXmcStatXX       -                -
lumi   lnN             1.2            -               -
'''.format(obs = tt.GetBinContent(tt.FindBin(SR)) + dy.GetBinContent(dy.FindBin(SR)), 
           fs  = tt.GetBinContent(tt.FindBin(SR)), other = dy.GetBinContent(dy.FindBin(SR)), 
           label = label, fs_int = int(tt.GetBinContent(tt.FindBin(SR))))
        outputFile = open('datacards/datacards_{scan}/{scan}/template_{sr}.txt'.format(scan=scan.name,sr=label),'w')
        outputFile.write(datacard)
        outputFile.close()

    # no shapes anymore :(
    # shapes = r.TFile('datacards/datacards_%s/shapes.root'%scan.name,'recreate'); 
    # tt.Write(); dy.Write(); data.Write();
    # shapes.Close()
    
    

def runCmd(cmd):
    os.chdir (cmd[2])
    os.system(cmd[0])
    os.system(cmd[1])
    return True


def produceLimits( njobs ):
    basedir = os.path.abspath('.')+'/datacards/datacards_{name}/{name}/'.format(name=scan.name)
    subdirs = [i for i in os.listdir(basedir) if os.path.isdir(basedir+i)]
    pool = Pool(njobs)
    tasks = []
    for ind,d in enumerate(subdirs):
        mass = ''.join(i for i in d.split('_') if i.isdigit() )
        fd = basedir+d
        srCards = ''
        for SR, label in scan.shortLabels.items():
            srCards += ' {fd}/datacard_{mass}_{label}.txt'.format(fd=fd,mass=d,label=label)
        dc_name    = '{fd}/datacard_{suffix}.txt'.format(fd=fd,suffix=d)
        runcmd = 'combineCards.py -S {bstr} > {final_dc}'.format(bstr=srCards,final_dc=dc_name)
        combinecmd = 'combine -m {mass} -M Asymptotic {dc_name}'.format(mass=mass,dc_name=dc_name)
        print runcmd
        tasks.append([runcmd,combinecmd, fd])
    pool.map(runCmd, tasks)
    print 'hadding everything'
    haddcmd = 'hadd -f {bd}/{name}_allLimits.root {bd}/*/higgs*.root'.format(name=scan.name,bd=basedir)
    os.system(haddcmd)


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


def getEffMapsSys(sys):
    print ('getting scan for sys', sys) if len(sys) > 0 else 'getting nominal scan'
    theCuts = scan.cuts_norm
    for rpl in replaceCutsForSys[sys]:
        theCuts = theCuts.replace(rpl[0],rpl[1])
                
    effMap = scan.tree.getTH3F(1., 'nPass_norm'+sys, scan.srID+':'+scan.yvar+':'+scan.xvar,
                               scan.xbins.n+1, scan.xbins._min-scan.xbins.w/2.,
                               scan.xbins._max+scan.xbins.w/2., scan.ybins.n+1,
                               scan.ybins._min-scan.ybins.w/2., 
                               scan.ybins._max+scan.ybins.w/2., 
                               scan.srIDMax+1, -0.5, scan.srIDMax+0.5, theCuts, '',
                               scan.xtitle, scan.ytitle, scan.ztitle, 
                               extraWeightsForSys[sys])
    effMap.Divide(scan.ngen_3d)
    return effMap

def getSREffMaps():
    if hasattr(scan,'maps'): 
        if '' in scan.maps: effmap =  scan.maps[''] 
    else: effmap = getEffMapsSys('')

    tmp = effmap.Project3D('yx'); tmp.Reset()

    scan.effMapSR = {}
    effmaps = r.TFile('effmaps.root','recreate')
    for SR, label in scan.SRLabels.items():
        scan.effMapSR[SR] = tmp.Clone()
        for i in range(1, scan.effMapSR[SR].GetXaxis().GetNbins()+1):
            for j in range(1, scan.effMapSR[SR].GetYaxis().GetNbins()+1):
                bin = effmap.FindBin(effmap.GetXaxis().GetBinCenter(i),
                                     effmap.GetYaxis().GetBinCenter(j),
                                     SR)
                scan.effMapSR[SR].SetBinContent(i, j, effmap.GetBinContent(bin))
        scan.effMapSR[SR].SetName(label)
        scan.effMapSR[SR].SetTitle(label)
        scan.effMapSR[SR].Write()
    effmaps.Close()

def PutHistosIntoRootFiles():
    
    for i in range(1, scan.norm.GetXaxis().GetNbins()+1):
        for j in range(1, scan.norm.GetYaxis().GetNbins()+1):
            xval = scan.norm.GetXaxis().GetBinCenter(i)
            yval = scan.norm.GetYaxis().GetBinCenter(j)
            if (yval > xval): continue
            massString = 'mSbottom_%.0f_mchi2_%.0f'%(xval, yval)
            
            sysHistos = {}
            for sys, histo in scan.maps.items():
                out = histo.ProjectionZ('%4.0f_%4.0f'%(xval,yval), i, i, j, j)
                out.Scale(scan.br*lumi*scan.xsecs[xval][0])
                out.SetName(massString + sys)
                sysHistos[sys] = out
            if sysHistos[''].Integral() == 0: continue # if no sensitivity


            helper.ensureDirectory('datacards/datacards_{scan}/{scan}/{mass}/'.format(scan=scan.name,mass=massString))
            for SR, label in scan.shortLabels.items():
                template = open('datacards/datacards_{scan}/{scan}/template_{sr}.txt'.format(scan=scan.name,sr=label),'r').read()
                template = template.replace('XXXSIGNAL',massString)
                template = template.replace('XXXSIGRATE', '%f'%sysHistos[''].GetBinContent(sysHistos[''].FindBin(SR)))
                for sys in scan.SysString.split():
                    up = sysHistos[sys+'Up']  .GetBinContent(sysHistos[sys+'Up']  .FindBin(SR))
                    dn = sysHistos[sys+'Down'].GetBinContent(sysHistos[sys+'Down'].FindBin(SR))
                    nom= sysHistos['']        .GetBinContent(sysHistos['']        .FindBin(SR))
                    template = template.replace('XX'+sys+'XX', '%4.2f/%4.2f'%(dn/nom,up/nom))
                mcStat = sysHistos[''].GetBinError(sysHistos[''].FindBin(SR)) / sysHistos[''].GetBinContent(sysHistos[''].FindBin(SR))
                template = template.replace('XXmcStatXX', '%f'%(1 + mcStat))
                card = open('datacards/datacards_{scan}/{scan}/{mass}/datacard_{mass}_{label}.txt'.format(scan=scan.name,
                                                                                                          mass=massString,
                                                                                                          label=label),'w')
                card.write(template)
                card.close()
                
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

    global replaceCutsForSys, extraWeightsForSys
    replaceCutsForSys = {'':      [],
                         'jecUp'  : [['nBJetMedium35_Edge','nBJetMedium35_jecUp_Edge'],
                                     ['mt2_Edge','mt2_jecUp_Edge'],
                                     ['nll_Edge', 'nll_jecUp_Edge'],
                                     ['met_Edge', 'met_jecUp_Edge'],
                                     ['nJetSel_Edge','nJetSel_jecUp_Edge']],
                         'jecDown': [['nBJetMedium35_Edge','nBJetMedium35_jecDn_Edge'],
                                     ['mt2_Edge','mt2_jecDn_Edge'],
                                     ['nll_Edge', 'nll_jecDn_Edge'],
                                     ['met_Edge', 'met_jecDn_Edge'],
                                     ['nJetSel_Edge','nJetSel_jecDn_Edge']],
                         'bHeUp'  : [],
                         'bHeDown': [],
                         'bLiUp'  : [],
                         'bLiDown': [],
                         'ElUp'   : [],
                         'ElDown' : [],
                         'MuUp'   : [],
                         'MuDown' : []}
    
    extraWeightsForSys = {''       : '1',
                          'jecUp'  : '1',
                          'jecDown': '1',
                          'bHeUp'  : 'weight_btagsf_heavy_UP_Edge/weight_btagsf_Edge',
                          'bHeDown': 'weight_btagsf_heavy_DN_Edge/weight_btagsf_Edge',
                          'bLiUp'  : 'weight_btagsf_light_UP_Edge/weight_btagsf_Edge',
                          'bLiDown': 'weight_btagsf_light_DN_Edge/weight_btagsf_Edge',
                          'ElUp'   : 'weight_LepSF_ElUp_Edge / weight_LepSF_Edge',
                          'ElDown' : 'weight_LepSF_ElDn_Edge / weight_LepSF_Edge',
                          'MuUp'   : 'weight_LepSF_MuUp_Edge / weight_LepSF_Edge',
                          'MuDown' : 'weight_LepSF_MuDn_Edge / weight_LepSF_Edge'}


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
        if scan.makeMCDatacards:
            print 'preparing datacards from MC'
            a = makeMCDatacards()


        ## Load the number of generated events to produce efficiency maps per systematic
        scan.dummy = scan.tree.getTH3F(1., 'dummy', '1:1:1',
                                       scan.xbins.n+1, scan.xbins._min-scan.xbins.w/2.,
                                       scan.xbins._max+scan.xbins.w/2., scan.ybins.n+1,
                                       scan.ybins._min-scan.ybins.w/2., 
                                       scan.ybins._max+scan.ybins.w/2., 
                                       scan.srIDMax+1, -0.5, scan.srIDMax+0.5, '1', '',
                                       scan.xtitle, scan.ytitle, scan.ztitle, 
                                       '1')
        if scan.has3DGen:
            scan.ngen = scan.tree.blocks[0].samples[0].smsCount ## take the first slice's ngen histo
            for ind,i in enumerate(scan.tree.blocks[0].samples):
                if ind: ## do not add the first one twice
                    scan.ngen.Add(i.smsCount, 1.) ## add all others
        else:
            scan.make3DGen()

        print 'this is the type of scan.ngen before', type(scan.ngen)
        newbinning   = adaptBinning(scan.dummy, scan.ngen)
        scan.ngen_2d = newbinning[0]
        scan.ngen_3d = newbinning[1] ## this one has ngen in every single bin. for every SR. and it's 3D, so that's cool
        print 'this is the type of scan.ngen after', type(scan.ngen)
        getSREffMaps()

        scan.SysString = 'El Mu jec bHe bLi'
        scan.maps = {}
        scan.maps['']   = getEffMapsSys('')
        for sys in scan.SysString.split():
            scan.maps[sys+'Up']   = getEffMapsSys(sys+'Up')
            scan.maps[sys+'Down'] = getEffMapsSys(sys+'Down')

        scan.norm = scan.maps['']

        scan.xy = scan.norm.Project3D('xy') # this means x versus y. so x is on the y-axis
        scan.xz = scan.norm.Project3D('xz')
        scan.yx = scan.norm.Project3D('yx') # that's the inclusive efficiency map
        scan.yz = scan.norm.Project3D('yz')
        scan.zx = scan.norm.Project3D('zx')
        scan.zy = scan.norm.Project3D('zy')
 
        PutHistosIntoRootFiles()

        

    if opts.reloadLimits:
        print 'reloading limits and datacards'
 #        fillAndSaveDatacards(dobs)
        produceLimits(6)

    os.system('mkdir -p mkdir -p makeExclusionPlot/config/%s/'%scan.paper)
    scan.makeExclusion()
    scan.makePrettyPlots()

    ## save the scan object in a pickle file to save time the second time around.
    pickle.dump(scan, open('datacards/datacards_{name}/{name}/{name}.pkl'.format(name=scan.name),'w') )
    print 'marc is stupid'
