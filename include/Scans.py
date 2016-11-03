import sys, copy
import Sample
import ROOT as r
import CutManager as CutManager

cuts = CutManager.CutManager()

class binning:
    def __init__(self, mi, ma, w):
        self._min = mi
        self._max = ma
        self.w  = w 
        self.n  = abs(ma-mi)/w

class Scan(object):
    def __init__(self, name):
        self.cuts_norm = '1'
        self.makeMCDatacards = False
        self.name  = name
        self.zminUL  = 1e-3; self.zmaxUL = 1e3
        self.xbins = binning(100,100,1000); self.ybins = binning(100,100,1000)
        self.xtitle = 'mother'; self.ytitle = 'daughter'; self.ztitle = 'SR-ID'
        self.br    = 1.
        self.xsecFile = 0
        self.zmaxEff = 0.5
        self.tree = 0
        self.paper = 'SUS66666'
        self.has3DGen = True
        self.xvar = 'GenSusyMScan1_Edge'
        self.yvar = 'GenSusyMScan2_Edge'
        self.loadData()
        self.loadXsecs()

    def __getstate__(self): 
        return self.__dict__
    def __setstate__(self, d): 
        self.__dict__.update(d)

    def loadData(self):
        if self.name == 'Edge_ICHEP':
            self.makeMCDatacards = True
            self.paper = 'SUS15011'
            self.datasets = ['SMS_T6bbllslepton_mSbottom400to575_mLSP150to550',
                             'SMS_T6bbllslepton_mSbottom600to775_mLSP150to725',
                             'SMS_T6bbllslepton_mSbottom800to950_mLSP150to900']
            self.xbins = binning(400,950,25)
            self.ybins = binning(200,900,25)
            self.xvar = 'GenSusyMScan1_Edge'
            self.yvar = 'GenSusyMScan2_Edge'
            self.cuts_norm = cuts.AddList([cuts.SignalRegionBaseLineNoTrigger, cuts.SF])
            self.cuts_norm = self.cuts_norm.replace(cuts.twoLeptons, 'nPairLep_Edge > 0') ## remove the filters, ugly.
            self.zminUL = 1e-3; self.zmaxUL = 1e3
            self.zmaxEff = 0.30
            self.xsecFile = ('datacards/sbottomXsec.txt')
            self.regions = []
            self.xtitle = 'm_{sbottom}'; self.ytitle = 'm_{neu2}'
            self.srID   = '0*(lepsMll_Edge < 81) + 1*(lepsMll_Edge > 101) + 0*(nll_Edge < 21.) + 2*(nll_Edge > 21.)'
            self.srIDMax = 3 # its good to avoid empty bins to be in the safe side
            self.shortLabels = { 0 : 'lowmll_lownll',
                                 1 : 'highmll_lownll',
                                 2 : 'lowmll_highnll',
                                 3 : 'highmll_lownll'}

                }
            self.SRLabels = { 0 : 'Low  m_{ll} / t#bar{t}-like ',
                              1 : 'High m_{ll} / t#bar{t}-like ',
                              2 : 'Low  m_{ll} / Non t#bar{t}-like ',
                              3 : 'High m_{ll} / Non t#bar{t}-like '}

# if self.name == 'Edge_ICHEP':

        #     self.makeMCDatacards = True
        #     self.paper = 'SUS15011'
        #     self.datasets = ['SMS_T6bbllslepton_mSbottom400to575_mLSP150to550',
        #                      'SMS_T6bbllslepton_mSbottom600to775_mLSP150to725',
        #                      'SMS_T6bbllslepton_mSbottom800to950_mLSP150to900']
        #     self.xbins = binning(400,950,25)
        #     self.ybins = binning(200,900,25)
        #     self.xvar = 'GenSusyMScan1_Edge'
        #     self.yvar = 'GenSusyMScan2_Edge'
        #     self.cuts_norm = cuts.AddList([cuts.SignalRegionBaseLineNoTrigger, cuts.SF])
        #     self.cuts_norm = self.cuts_norm.replace(cuts.twoLeptons, 'nPairLep_Edge > 0') ## remove the filters, ugly.
        #     self.zminUL = 1e-3; self.zmaxUL = 1e3
        #     self.zmaxEff = 0.30
        #     self.xsecFile = ('datacards/sbottomXsec.txt')
        #     self.regions = []
        #     # 4 signal regions
        #     for mll in ['belowZ','highinc']:
        #         for nll in ['lowNll', 'highNll']:
        #             for mt2 in ['incMT2']:
        #                 self.regions.append([mll,nll,mt2])
        #     self.xtitle = 'm_{sbottom}'; self.ytitle = 'm_{neu2}'
        #     self.srID   = '1*(lepsMll_Edge < 81) + 2*(lepsMll_Edge > 81)*(lepsMll_Edge < 101) + 3*(lepsMll_Edge > 101)*(lepsMll_Edge < 200.) + 4*(lepsMll_Edge > 200)*(lepsMll_Edge < 300) + 5*(lepsMll_Edge > 300) + 0*(nll_Edge < 21.) + 10*(nll_Edge > 21.) + 0*(mt2_Edge < 80.) + 100*(mt2_Edge > 80.)'

        # if self.name == 'Edge_MLLbins_MT2inc':
        #     self.makeMCDatacards = True
        #     self.paper = 'SUS15011'
        #     self.datasets = ['SMS_T6bbllslepton_mSbottom400to575_mLSP150to550',
        #                      'SMS_T6bbllslepton_mSbottom600to775_mLSP150to725',
        #                      'SMS_T6bbllslepton_mSbottom800to950_mLSP150to900']
        #     self.xbins = binning(400,950,25)
        #     self.ybins = binning(200,900,25)
        #     self.xvar = 'GenSusyMScan1_Edge'
        #     self.yvar = 'GenSusyMScan2_Edge'
        #     self.cuts_norm = cuts.AddList([cuts.SignalRegionBaseLineNoTrigger, cuts.SF])
        #     self.cuts_norm = self.cuts_norm.replace(cuts.twoLeptons, 'nPairLep_Edge > 0') ## remove the filters, ugly.
        #     self.zminUL = 1e-3; self.zmaxUL = 1e3
        #     self.zmaxEff = 0.30
        #     self.xsecFile = ('datacards/sbottomXsec.txt')
        #     self.regions = []
        #     # 10 signal regions
        #     for mll in ['belowZ','onZ','aboveZ','highMass', 'vHighMass' ]:
        #         for nll in ['lowNll', 'highNll']:
        #             for mt2 in ['incMT2']:
        #                 self.regions.append([mll,nll,mt2])
        #     self.xtitle = 'm_{sbottom}'; self.ytitle = 'm_{neu2}'
        #     self.srID   = '1*(lepsMll_Edge < 81) + 2*(lepsMll_Edge > 81)*(lepsMll_Edge < 101) + 3*(lepsMll_Edge > 101)*(lepsMll_Edge < 200.) + 4*(lepsMll_Edge > 200)*(lepsMll_Edge < 300) + 5*(lepsMll_Edge > 300) + 0*(nll_Edge < 21.) + 10*(nll_Edge > 21.) + 0*(mt2_Edge < 80.) + 100*(mt2_Edge > 80.)'

        # if self.name == 'Edge_MLLbins_MT2bins':
        #     self.paper = 'SUS15011'
        #     self.makeMCDatacards = True
        #     self.datasets = ['SMS_T6bbllslepton_mSbottom400to575_mLSP150to550',
        #                      'SMS_T6bbllslepton_mSbottom600to775_mLSP150to725',
        #                      'SMS_T6bbllslepton_mSbottom800to950_mLSP150to900']
        #     self.xbins = binning(400,950,25)
        #     self.ybins = binning(200,900,25)
        #     self.xvar = 'GenSusyMScan1_Edge'
        #     self.yvar = 'GenSusyMScan2_Edge'
        #     self.cuts_norm = cuts.AddList([cuts.SignalRegionBaseLineNoTrigger, cuts.SF])
        #     self.cuts_norm = self.cuts_norm.replace(cuts.twoLeptons, 'nPairLep_Edge > 0') ## remove the filters, ugly.
        #     self.zminUL = 1e-3; self.zmaxUL = 1e3
        #     self.zmaxEff = 0.30
        #     self.xsecFile = ('datacards/sbottomXsec.txt')
        #     self.regions = []
        #     # 20 signal regions
        #     for mll in ['belowZ','onZ','aboveZ','highMass', 'vHighMass' ]:
        #         for nll in ['lowNll', 'highNll']:
        #             for mt2 in ['lowMT2','highMT2']:
        #                 self.regions.append([mll,nll,mt2])
        #     self.xtitle = 'm_{sbottom}'; self.ytitle = 'm_{neu2}'
        #     self.srID   = '1*(lepsMll_Edge < 81) + 2*(lepsMll_Edge > 81)*(lepsMll_Edge < 101) + 3*(lepsMll_Edge > 101)*(lepsMll_Edge < 200.) + 4*(lepsMll_Edge > 200)*(lepsMll_Edge < 300) + 5*(lepsMll_Edge > 300) + 0*(nll_Edge < 21.) + 10*(nll_Edge > 21.) + 0*(mt2_Edge < 80.) + 100*(mt2_Edge > 80.)'

        # if self.name == 'Edge_MLLbins_MT2cut':
        #     self.paper = 'SUS15011'
        #     self.makeMCDatacards = True
        #     self.datasets = ['SMS_T6bbllslepton_mSbottom400to575_mLSP150to550',
        #                      'SMS_T6bbllslepton_mSbottom600to775_mLSP150to725',
        #                      'SMS_T6bbllslepton_mSbottom800to950_mLSP150to900']
        #     self.xbins = binning(400,950,25)
        #     self.ybins = binning(200,900,25)
        #     self.xvar = 'GenSusyMScan1_Edge'
        #     self.yvar = 'GenSusyMScan2_Edge'
        #     self.cuts_norm = cuts.AddList([cuts.SignalRegionBaseLineNoTrigger, cuts.SF])
        #     self.cuts_norm = self.cuts_norm.replace(cuts.twoLeptons, 'nPairLep_Edge > 0') ## remove the filters, ugly.
        #     self.zminUL = 1e-3; self.zmaxUL = 1e3
        #     self.zmaxEff = 0.30
        #     self.xsecFile = ('datacards/sbottomXsec.txt')
        #     self.regions = []
        #     # 20 signal regions
        #     for mll in ['belowZ','onZ','aboveZ','highMass', 'vHighMass' ]:
        #         for nll in ['lowNll', 'highNll']:
        #             for mt2 in ['highMT2']:
        #                 self.regions.append([mll,nll,mt2])
        #     self.xtitle = 'm_{sbottom}'; self.ytitle = 'm_{neu2}'
        #     self.srID   = '1*(lepsMll_Edge < 81) + 2*(lepsMll_Edge > 81)*(lepsMll_Edge < 101) + 3*(lepsMll_Edge > 101)*(lepsMll_Edge < 200.) + 4*(lepsMll_Edge > 200)*(lepsMll_Edge < 300) + 5*(lepsMll_Edge > 300) + 0*(nll_Edge < 21.) + 10*(nll_Edge > 21.) + 0*(mt2_Edge < 80.) + 100*(mt2_Edge > 80.)'

        if self.name == 'Edge_MLLbins_HTbins':
            self.paper = 'SUS15011'
            self.makeMCDatacards = True
            self.datasets = ['SMS_T6bbllslepton_mSbottom400to575_mLSP150to550',
                             'SMS_T6bbllslepton_mSbottom600to775_mLSP150to725',
                             'SMS_T6bbllslepton_mSbottom800to950_mLSP150to900']
            self.xbins = binning(400,950,25)
            self.ybins = binning(200,900,25)
            self.xvar = 'GenSusyMScan1_Edge'
            self.yvar = 'GenSusyMScan2_Edge'
            self.cuts_norm = cuts.AddList([cuts.SignalRegionBaseLineNoTrigger, cuts.SF])
            self.cuts_norm = self.cuts_norm.replace(cuts.twoLeptons, 'nPairLep_Edge > 0') ## remove the filters, ugly.
            self.zminUL = 1e-3; self.zmaxUL = 1e3
            self.zmaxEff = 0.30
            self.xsecFile = ('datacards/sbottomXsec.txt')
            self.regions = []
            # 20 signal regions
            for mll in ['belowZ','onZ','aboveZ','highMass', 'vHighMass' ]:
                for nll in ['lowNll', 'highNll']:
                    for ht in ['lowHT','medHT','highHT']:
                        self.regions.append([mll,nll,ht])
            self.xtitle = 'm_{sbottom}'; self.ytitle = 'm_{neu2}'
            self.srID   = '1*(lepsMll_Edge < 81) + 2*(lepsMll_Edge > 81)*(lepsMll_Edge < 101) + 3*(lepsMll_Edge > 101)*(lepsMll_Edge < 200.) + 4*(lepsMll_Edge > 200)*(lepsMll_Edge < 300) + 5*(lepsMll_Edge > 300) + 0*(nll_Edge < 21.) + 10*(nll_Edge > 21.) + 0*(htJet35j_Edge < 400.) + 100*(htJet35j_Edge > 400.)*(htJet35j_Edge < 800.) + 200*(htJet35j_Edge > 800.)'


        if self.name == 'T6bbslepton':
            self.paper = 'SUS15011'
            self.datasets = ['SMS_T6bbllslepton_mSbottom400To575_mLSP150To550',
                             'SMS_T6bbllslepton_mSbottom600To775_mLSP150To725',
                             'SMS_T6bbllslepton_mSbottom800To950_mLSP150To900']
            self.xbins = binning(400,900,25)
            self.ybins = binning(200,900,25)
            self.xvar = 'GenSusyMScan1_Edge'
            self.yvar = 'GenSusyMScan2_Edge'
            self.cuts_norm = cuts.AddList([cuts.METJetsSignalRegion, cuts.GoodLeptonSFNoTrigger()]) ## trigger not available in fastsim
            self.cuts_norm = self.cuts_norm.replace(cuts.twoLeptons, 'nPairLep_Edge > 0') ## remove the filters, ugly.
            self.zminUL = 1e-3; self.zmaxUL = 1e3
            self.zmaxEff = 0.30
            self.xsecFile = ('/afs/cern.ch/user/m/mdunser/edgeSW/Framework/datacards/sbottomXsec.txt')
            self.regions = []
            for eta in ['central','forward']:#,'inclusive']:
                for mass in ['lowMass', 'belowZ', 'onZ', 'aboveZ', 'highMass']:#,'allMass']:
                    for nb in ['incb', '0b', '1b', '2b']:
                        self.regions.append([eta,mass,nb])
            #self.regions.append([eta,mass,nb] for eta in ['central','forward'] for mass in ['allMass', 'lowMass', 'belowZ', 'onZ', 'aboveZ', 'highMass'] for nb in ['incb', '0b', '1b', '2b'] )
            self.xtitle = 'm_{sbottom}'; self.ytitle = 'm_{neu2}'

        if self.name == 'T6bbsleptonMET150':
            self.paper = 'SUS15011MET150'
            self.datasets = ['SMS_T6bbllslepton_mSbottom400To575_mLSP150To550',
                             'SMS_T6bbllslepton_mSbottom600To775_mLSP150To725',
                             'SMS_T6bbllslepton_mSbottom800To950_mLSP150To900']
            self.xbins = binning(400,900,25)
            self.ybins = binning(200,900,25)
            self.xvar = 'GenSusyMScan1_Edge'
            self.yvar = 'GenSusyMScan2_Edge'
            self.cuts_norm = cuts.AddList([cuts.METJetsSignalRegion, cuts.GoodLeptonSFNoTrigger(), 'met_Edge > 150']) ## trigger not available in fastsim
            self.cuts_norm = self.cuts_norm.replace(cuts.twoLeptons, 'nPairLep_Edge > 0') ## remove the filters, ugly.
            self.zminUL = 1e-3; self.zmaxUL = 1e3
            self.zmaxEff = 0.30
            self.xsecFile = ('/afs/cern.ch/user/m/mdunser/edgeSW/Framework/datacards/sbottomXsec.txt')
            self.regions = []
            for eta in ['central','forward']:#,'inclusive']:
                for mass in ['lowMass', 'belowZ', 'onZ', 'aboveZ', 'highMass']:#,'allMass']:
                    for nb in ['incb', '0b', '1b', '2b']:
                        self.regions.append([eta,mass,nb])
            #self.regions.append([eta,mass,nb] for eta in ['central','forward'] for mass in ['allMass', 'lowMass', 'belowZ', 'onZ', 'aboveZ', 'highMass'] for nb in ['incb', '0b', '1b', '2b'] )
            self.xtitle = 'm_{sbottom}'; self.ytitle = 'm_{neu2}'

        if self.name == 'TChiChiSleps':
            self.paper = 'SUS16777'
            self.datasets = ['TChiChi_slep_mCha600_mLSP50',
                             'TChiChi_slep_mCha350_mLSP200']
            self.xbins = binning(350,650,50)
            self.ybins = binning(  0,250,50)
            self.xvar = 'GenSusyMChargino1_Edge'
            self.yvar = 'GenSusyMNeutralino_Edge'
            self.cuts_norm = cuts.AddList(['met_Edge > 150', 'nJetSel_Edge < 2', cuts.GoodLeptonAFNoTrigger()]) ## trigger not available in fastsim
            self.cuts_norm = self.cuts_norm.replace(cuts.twoLeptons, 'nPairLep_Edge > 0') ## remove the filters, ugly.
            self.zminUL = 1e-3; self.zmaxUL = 1e3
            self.zmaxEff = 0.30
            self.xsecFile = ('/afs/cern.ch/user/m/mdunser/edgeSW/Framework/datacards/chaChaXsec_wino.txt')
            self.regions = []
            for eta in ['central','forward']:#,'inclusive']:
                for mass in ['lowMass', 'belowZ', 'onZ', 'aboveZ', 'highMass']:#,'allMass']:
                    for nb in ['0b']:
                        self.regions.append([eta,mass,nb])
            #self.regions.append([eta,mass,nb] for eta in ['central','forward'] for mass in ['allMass', 'lowMass', 'belowZ', 'onZ', 'aboveZ', 'highMass'] for nb in ['incb', '0b', '1b', '2b'] )
            self.xtitle = 'm_{sbottom}'; self.ytitle = 'm_{neu2}'
            self.has3DGen = False

        if self.name == 'TChiNeuWZ':
            self.makeMCDatacards = True
            self.paper = 'SUS16888'
            self.datasets = ['TChiNeu_WZ_mCha350_mLSP100',
                             'TChiNeu_WZ_mCha200_mLSP100',
                             'TChiNeu_WZ_mCha350_mLSP20' ]
            self.xbins = binning(150,400,50)
            self.ybins = binning(  0,140,20)
            self.xvar = 'GenSusyMNeutralino2_Edge'
            self.yvar = 'GenSusyMNeutralino_Edge'
            self.cuts_norm = cuts.AddList([cuts.METJetsSignalRegion, cuts.GoodLeptonSFNoTrigger()])#, 'met_Edge > 200 && maxMjj_Edge > 50 && maxMjj_Edge < 150 && -1.*TMath::Log(lh_ana_met_data_Edge*lh_ana_mlb_data_Edge*lh_ana_a3d_data_Edge*lh_ana_zpt_data_Edge) > 19.0']) ## trigger not available in fastsim
            self.cuts_norm = self.cuts_norm.replace(cuts.twoLeptons, 'nPairLep_Edge > 0') ## remove the filters, ugly.
            self.zminUL = 1e-3; self.zmaxUL = 1e3
            self.zmaxEff = 0.30
            self.xsecFile = ('/afs/cern.ch/user/m/mdunser/edgeSW/Framework/datacards/chaNeuXsec_wino.txt')
            self.regions = []
            #for eta in ['central','forward']:#,'inclusive']:
            for eta in ['inclusive']:
                for mass in ['onZ']:#,'allMass']:
                    for nb in ['0b']:
                        self.regions.append([eta,mass,nb])
            #self.regions.append([eta,mass,nb] for eta in ['central','forward'] for mass in ['allMass', 'lowMass', 'belowZ', 'onZ', 'aboveZ', 'highMass'] for nb in ['incb', '0b', '1b', '2b'] )
            self.xtitle = 'm_{Chargino=Neutralino2}'; self.ytitle = 'm_{LSP}'
            self.has3DGen = False
            self.br = 0.1 ## Z goes only to leptons

    def loadXsecs(self):
        if not self.xsecFile:
            estring = 'ERROR: no xsec file specified for scan {name}!!!\nExiting...'.format(name=self.name)
            sys.exit(estring)
        xsecf = open(self.xsecFile, 'r')
        self.xsecs = eval(xsecf.read())
        xsecf.close()
        keys = self.xsecs.keys()
        self.xsec_histo = r.TH1F('x-sections in fb for %s'%self.name,'xsec_histo', len(keys), min(keys), max(keys) )
        self.xsec_histo.Sumw2()
        for key, value in self.xsecs.items():
            self.xsecs[key][0] = self.xsecs[key][0]*1000.
            self.xsecs[key][1] = self.xsecs[key][0]*0.01*self.xsecs[key][1]
            self.xsec_histo.SetBinContent(self.xsec_histo.FindBin(key), value[0])
            self.xsec_histo.SetBinError  (self.xsec_histo.FindBin(key), value[1]) ## it's a percent value

    def makeExclusion(self):
        print 'making the exclusions!!'
        limitfile = r.TFile('datacards/datacards_{name}/{name}/{name}_allLimits.root'.format(name=self.name),'READ')
        limittree = limitfile.Get('limit')
        ## observed limit
        self.ex_obs = copy.deepcopy(self.yx)
        self.ex_obs.SetTitle("nominal exclusion plot observed"); self.ex_obs.SetName('ex_obs')
        self.ex_obs.Reset()
        ## observed + 1 sigma xsec limit
        self.ex_obs_p1s = copy.deepcopy(self.yx)
        self.ex_obs_p1s.SetTitle("nominal exclusion plot observed + 1 sig xsec"); self.ex_obs_p1s.SetName('ex_obs_p1s')
        self.ex_obs_p1s.Reset()
        ## observed - 1 sigma xsec limit
        self.ex_obs_m1s = copy.deepcopy(self.yx)
        self.ex_obs_m1s.SetTitle("nominal exclusion plot observed - 1 sig xsec"); self.ex_obs_m1s.SetName('ex_obs_m1s')
        self.ex_obs_m1s.Reset()
        ## expected limit
        self.ex_exp = copy.deepcopy(self.yx)
        self.ex_exp.SetTitle("nominal exclusion plot expected"); self.ex_exp.SetName('ex_exp')
        self.ex_exp.Reset()
        ## expected limit +1 sigma(exp)
        self.ex_exp_p1s = copy.deepcopy(self.yx)
        self.ex_exp_p1s.SetTitle("nominal exclusion plot expected + 1 sigma exp"); self.ex_exp_p1s.SetName('ex_exp_p1s')
        self.ex_exp_p1s.Reset()
        ## expected limit +1 sigma(exp)
        self.ex_exp_m1s = copy.deepcopy(self.yx)
        self.ex_exp_m1s.SetTitle("nominal exclusion plot expected - 1 sigma exp"); self.ex_exp_m1s.SetName('ex_exp_m1s')
        self.ex_exp_m1s.Reset()
    
        for point in limittree:
            mass      = str(int(point.mh))
            massx     = int(mass[:3]); massy = int(mass[3:])
            if point.quantileExpected == -1:
                self.ex_obs    .Fill(massx, massy, point.limit)
                self.ex_obs_p1s.Fill(massx, massy, point.limit*(self.xsecs[massx][0]+self.xsecs[massx][1])/self.xsecs[massx][0])
                self.ex_obs_m1s.Fill(massx, massy, point.limit*(self.xsecs[massx][0]-self.xsecs[massx][1])/self.xsecs[massx][0])
            elif 0.49 < point.quantileExpected < 0.51:
                self.ex_exp    .Fill(massx, massy, point.limit)
            elif 0.15 < point.quantileExpected < 0.17:
                self.ex_exp_p1s.Fill(massx, massy, point.limit)
            elif 0.83 < point.quantileExpected < 0.85:
                self.ex_exp_m1s.Fill(massx, massy, point.limit)
        zmax = self.ex_obs.GetMaximum()
        self.ex_obs    .GetZaxis().SetRangeUser(0.,10.)
        self.ex_obs_p1s.GetZaxis().SetRangeUser(0.,10.)
        self.ex_obs_m1s.GetZaxis().SetRangeUser(0.,10.)
        self.ex_exp    .GetZaxis().SetRangeUser(0.,10.)
        self.ex_exp_p1s.GetZaxis().SetRangeUser(0.,10.)
        self.ex_exp_m1s.GetZaxis().SetRangeUser(0.,10.)
        out  = r.TFile('limits_%s.root'%self.name,'recreate')
        self.ex_obs.Write()
        out.Close()
    
    def getSmoothedGraph(self,h_orig):
        tmp_name       = (h_orig.GetName()+'_smoothed').replace('Graph2D_from_','')
        tmp_name_graph = tmp_name+'_graph'
        global h_smoothed, smoothed_2dg, smoothed_gl, smoothed_g
        h_smoothed = h_orig.Clone(tmp_name)
        smoothedbins = []
        for i in range(h_smoothed.GetNbinsX()+1):
            for j in range(h_smoothed.GetNbinsY()+1):
                if j > 0 and not h_smoothed.GetBinContent(i,j):
                    h_smoothed.SetBinContent(i,j,h_smoothed.GetBinContent(i-1,j))
                    smoothedbins.append((i,j))
                    break
        h_smoothed.Smooth(1, 'k3a')
        for i,j in smoothedbins:
            h_smoothed.SetBinContent(i,j,0.)
        bruno = 0#'gay'
        if bruno:
            h_smoothed.SetMinimum(0)
            h_smoothed.SetMaximum(2)
            h_smoothed.SetContour(4)
            c= r.TCanvas()
            h_smoothed.Draw("contz list")
            r.gPad.Update()
            obj = r.gROOT.GetListOfSpecials().FindObject("contours")
            ls = obj.At(1)
            graph = []
            for l in range(ls.GetSize()):
                gr = ls.At(l).Clone()
                n_points = gr.GetN()
                graph.append((n_points,gr))
            graph.sort(reverse=True)
            name = h_smoothed.GetName()
            for i,g in enumerate(graph):
                g[1].SetName(name+(str(i)if i>0 else "")+'_graph')
            if len(graph)==1: graph.append(graph[0])
        else:
            smcopy = copy.deepcopy(h_smoothed)
            smoothed_2dg = r.TGraph2D(smcopy)
            xbinsize = 5.; ybinsize = 5.
            smoothed_2dg.SetNpx( int((smoothed_2dg.GetXmax() - smoothed_2dg.GetXmin())/xbinsize) )
            smoothed_2dg.SetNpy( int((smoothed_2dg.GetYmax() - smoothed_2dg.GetYmin())/ybinsize) )
            brunosucksballs = smoothed_2dg.GetHistogram() ## have to call this, otherwise root will freak out
            c = smoothed_2dg.GetContourList(1.0)
            smoothed_g   = c[max( (i.GetN(),j) for j,i in enumerate( c )  )[1]]
            smoothed_g.SetName(tmp_name_graph)
        setattr(self, tmp_name      , copy.deepcopy(h_smoothed) )
        setattr(self, tmp_name_graph, copy.deepcopy(smoothed_g) )#copy.deepcopy(graph[0][1]) )
    
    def makeULPlot(self):
        print 'making the UL histo with holes filled etc.'
        h_rhisto = copy.deepcopy(self.ex_obs)
        for i in range(h_rhisto.GetNbinsX()+1):
            for j in range(h_rhisto.GetNbinsY()+1):
                xs = self.xsecs[int(h_rhisto.GetXaxis().GetBinCenter(i))][0]
                h_rhisto.SetBinContent(i,j,h_rhisto.GetBinContent(i,j)*xs/1000.)
        gr2d = r.TGraph2D(h_rhisto)
        xbinsize = 12.5; ybinsize = 12.5
        gr2d.SetNpx( int((gr2d.GetXmax() - gr2d.GetXmin())/xbinsize) )
        gr2d.SetNpy( int((gr2d.GetYmax() - gr2d.GetYmin())/ybinsize) )
        tmp_2d_histo = gr2d.GetHistogram()
        tmp_2d_histo.SetName('ul_histo')
        setattr(self, 'ul_histo', copy.deepcopy(tmp_2d_histo))
        print '... done making the UL histo'
                
    def saveULGraphsInFile(self):
        print 'saving the UL histo and graphs in a file'
        f = r.TFile('makeExclusionPlot/config/%s/%s_results.root'%(self.paper,self.name),'RECREATE')
        f.cd()
        self.ul_histo.Write()
        self.ex_exp_smoothed_graph    .Write()
        self.ex_exp_p1s_smoothed_graph.Write()
        self.ex_exp_m1s_smoothed_graph.Write()
        self.ex_obs_smoothed_graph    .Write()
        self.ex_obs_p1s_smoothed_graph.Write()
        self.ex_obs_m1s_smoothed_graph.Write()
        f.Close()
        print '... done saving the results file'
    
    def makePrettyPlots(self):
        print 'starting to make pretty plots'
        for t in ['_obs', '_obs_p1s', '_obs_m1s', '_exp', '_exp_p1s', '_exp_m1s']:
            tmp_2d_graph = r.TGraph2D(getattr(self, 'ex%s'%t))
            xbinsize = 12.5; ybinsize = 12.5
            tmp_2d_graph.SetNpx( int((tmp_2d_graph.GetXmax() - tmp_2d_graph.GetXmin())/xbinsize) )
            tmp_2d_graph.SetNpy( int((tmp_2d_graph.GetYmax() - tmp_2d_graph.GetYmin())/ybinsize) )
            tmp_2d_histo = tmp_2d_graph.GetHistogram()
            tmp_graph_list = tmp_2d_graph.GetContourList(1.0)
            tmp_graph = tmp_graph_list[max( (i.GetN(),j) for j,i in enumerate( tmp_graph_list )  )[1]]
            setattr(self, 'ex%s_graph'%t, copy.deepcopy(tmp_graph    ) )
            setattr(self, 'ex%s_2dg'  %t, copy.deepcopy(tmp_2d_graph ) )
            setattr(self, 'ex%s_2dh'  %t, copy.deepcopy(tmp_2d_histo ) )
            self.getSmoothedGraph(getattr(self,'ex%s_2dh'%t))
        self.makeULPlot()
        self.saveULGraphsInFile()
        print '... done making pretty plots'

    def make3DGen(self):
        self.ngen = r.TH3F('3dngen', '3dngen', self.xbins.n, self.xbins._min, self.xbins._max, 
                                               self.ybins.n, self.ybins._min, self.ybins._max,
                                               1, 0 , 1)
        #self.ngen = self.norm.Clone('smsCountHisto') ##scan.tree.blocks[0].samples[0].smsCount ## take the first slice's ngen histo
        #self.ngen.Reset()
        for ind,i in enumerate(self.tree.blocks[0].samples):
            massx = int( ''.join(j for j in i.name.split('_')[-2] if j.isdigit() ) )
            massy = int( ''.join(j for j in i.name.split('_')[-1] if j.isdigit() ) )
            print 'at ', i, 'getting masses', massx, massy
            self.ngen.Fill(massx,massy, 0.5, i.count) ## add all others

