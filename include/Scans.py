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
        self.doTwoSigmas = False
    def __getstate__(self): 
        return self.__dict__
    def __setstate__(self, d): 
        self.__dict__.update(d)

    def loadData(self):
        if self.name == 'Edge_ICHEP2016':
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
                                 3 : 'highmll_highnll'}
            self.SRLabels = { 0 : 'Low  m_{ll} / t#bar{t}-like ',
                              1 : 'High m_{ll} / t#bar{t}-like ',
                              2 : 'Low  m_{ll} / Non t#bar{t}-like ',
                              3 : 'High m_{ll} / Non t#bar{t}-like '}

        if self.name == 'Edge_Moriond2017':
            self.makeMCDatacards = True; print 'cambiar esto'
            self.doTwoSigmas = True
            self.paper = 'SUS16034'
            self.datasets = ['SMS_T6bbllslepton_mSbottom400to575_mLSP150to550',
                             'SMS_T6bbllslepton_mSbottom600to775_mLSP150to725',
                             'SMS_T6bbllslepton_mSbottom800to950_mLSP150to900']
            self.xbins = binning(400,950,25)
            self.ybins = binning(200,900,25)
            self.xvar = 'GenSusyMScan1_Edge'
            self.yvar = 'GenSusyMScan2_Edge'
            self.cuts_norm = cuts.AddList([cuts.BaselineNoTrigger, cuts.SF, cuts.EdgeBaseline])
            self.cuts_norm = self.cuts_norm.replace(cuts.twoLeptons, 'nPairLep_Edge > 0') ## remove the filters, ugly.
            self.zminUL = 1e-3; self.zmaxUL = 1e3
            self.zmaxEff = 0.30
            self.xsecFile = ('datacards/sbottomXsec.txt')
            self.regions = []
            self.xtitle = 'm_{sbottom}'; self.ytitle = 'm_{neu2}'
            self.srID   = '0*(lepsMll_Edge < 60) + 1*(lepsMll_Edge > 60)*(lepsMll_Edge < 86) + 2*(lepsMll_Edge > 96)*(lepsMll_Edge < 150.) + 3*(lepsMll_Edge > 150)*(lepsMll_Edge < 200) + 4*(lepsMll_Edge > 200)*(lepsMll_Edge < 300) + 5*(lepsMll_Edge > 300)*(lepsMll_Edge < 400) + 6*(lepsMll_Edge > 400) + 0*(nll_Edge < 21.) + 7*(nll_Edge > 21.)'
            self.srIDMax = 13 # its good to avoid empty bins to be in the safe side
            self.shortLabels = { 0: 'lowmll_lownll',
                                 1: 'belowZ_lownll',
                                 2: 'aboveZ_lownll',
                                 3: 'highmll1_lownll',
                                 4: 'highmll2_lownll',
                                 5: 'highmll3_lownll',
                                 6: 'highmll4_lownll',
                                 7: 'lowmll_highnll',
                                 8: 'belowZ_highnll',
                                 9: 'aboveZ_highnll',
                                10: 'highmll1_highnll',
                                11: 'highmll2_highnll',
                                12: 'highmll3_highnll',
                                13: 'highmll4_highnll'}
            self.SRLabels = { 0: 'Low Z / t#bar{t}-like ',
                              1: 'Below Z / t#bar{t}-like',
                              2: 'Above Z / t#bar{t}-like',
                              3: 'High m_{ll} / t#bar{t}-like',
                              4: 'Very high m_{ll} / t#bar{t}-like',
                              5: 'Even higher m_{ll} / t#bar{t}-like',
                              6: 'Huge m_{ll} / t#bar{t}-like',
                              7: 'Low Z / Non t#bar{t}-like ',
                              8: 'Below Z / Non t#bar{t}-like',
                              9: 'Above Z / Non t#bar{t}-like',
                             10: 'High m_{ll} / Non t#bar{t}-like',
                             11: 'Very high m_{ll} / Non t#bar{t}-like',
                             12: 'Even higher m_{ll} / Non t#bar{t}-like',
                             13: 'Huge m_{ll} / Non t#bar{t}-like'}

        if self.name == 'CharNeu_Moriond2017':
            self.makeMCDatacards = False
            self.paper = 'SUS16034'
            self.datasets = ['work in progress :) ']
            self.xbins = binning(400,950,25)
            self.ybins = binning(200,900,25)
            self.xvar = 'GenSusyMScan1_Edge'
            self.yvar = 'GenSusyMScan2_Edge'
            self.cuts_norm = cuts.AddList([cuts.BaselineNoTrigger, cuts.SF, cuts.ewinoCharNeu])
            self.cuts_norm = self.cuts_norm.replace(cuts.twoLeptons, 'nPairLep_Edge > 0')
            self.zminUL = 1e-3; self.zmaxUL = 1e3
            self.zmaxEff = 0.30
            self.xsecFile = ('datacards/charneuXsec.txt')
            self.regions = []
            self.xtitle = 'm_{chi^{0}_{2}} = m_{chi^{+}_{1}}'; self.ytitle = 'm_{chi^{0}_{1}}'
            self.srID   = '0'
            self.srIDMax = 0 
            self.shortLabels = {0: 'only SR' }
            self.SRLabels    = {0: 'only SR' }

        if self.name == 'NeuNeu_Moriond2017':
            self.makeMCDatacards = False
            self.paper = 'SUS16034'
            self.datasets = ['work in progress :) ']
            self.xbins =  binning(400,950,25) # ### to do the proper binning
            self.ybins =  binning(200,900,25) # ### to do the proper binning
            self.xvar = 'GenSusyMScan1_Edge'
            self.yvar = 'GenSusyMScan2_Edge'
            self.cuts_norm = cuts.AddList([cuts.BaselineNoTrigger, cuts.SF, cuts.ewinoNeuNeu])
            self.cuts_norm = self.cuts_norm.replace(cuts.twoLeptons, 'nPairLep_Edge > 0')
            self.zminUL = 1e-3; self.zmaxUL = 1e3
            self.zmaxEff = 0.30
            self.xsecFile = ('datacards/neuneuXsec.txt')
            self.regions = []
            self.xtitle = 'm_{chi^{0}_{1}} = m_{chi^{0}_{1}}'; self.ytitle = 'm_{chi^{0}_{1}}'
            self.srID   = '0*(met_Edge > 100)*(met_Edge < 150) + 1*(met_Edge > 150)*(met_Edge < 250) + 2*(met_Edge > 250)'
            self.srIDMax = 2
            self.shortLabels = {0: 'lowmet',
                                1: 'medmet',
                                2: 'highmet'}
            self.SRLabels    = {0: '100 GeV < ME_{T} < 150 GeV',
                                1: '150 GeV < ME_{T} < 250 GeV',
                                2: '250 GeV > ME_{T}'}

        if self.name == 'StrongOnZ_Moriond2017':
            self.makeMCDatacards = False
            self.paper = 'SUS16034'
            self.datasets = ['T5ZZ work in progress :) ']
            self.xbins = binning(400,950,25) # ### to do the proper binning 
            self.ybins = binning(200,900,25) # ### to do the proper binning 
            self.xvar = 'GenSusyMScan1_Edge'
            self.yvar = 'GenSusyMScan2_Edge'
            self.cuts_norm = cuts.AddList([cuts.BaselineNoTrigger, cuts.SF,cuts.strongOnZBase])
            self.cuts_norm = self.cuts_norm.replace(cuts.twoLeptons, 'nPairLep_Edge > 0')
            self.zminUL = 1e-3; self.zmaxUL = 1e3
            self.zmaxEff = 0.30
            self.xsecFile = ('datacards/gluinoXsec.txt')
            self.regions = []
            self.xtitle = 'm_{#tilde{g}}'; self.ytitle = 'm_{chi^{0}_{1}}'
            self.srID   = '{bveto} * ( {j23} *( {ht500}*0 -10000*(1-{ht500}))  + {j45}*( {ht500}*1 -10000*(1-{ht500})) + {j6}*2) + {withb} * ({j23}*({ht200}*3 + -100000*(1-{ht200})) + {j45}*({ht200}*4 -10000*(1-{ht200})) + {j6}*5) -100000*(!{bveto} && !{withb}) + 0*(met_Edge > 100)*(met_Edge < 150) + 6*(met_Edge > 150)*(met_Edge < 225) + 12*(met_Edge > 225)*(met_Edge < 300) + 18*(met_Edge > 300)'.format(bveto = cuts.strongOnZBVeto, withb= cuts.strongOnZWithB, j23 = '(nJetSel_Edge == 2 || nJetSel_Edge == 3)', j45 = '(nJetSel_Edge == 4 || nJetSel_Edge ==5 )', j6 = '(nJetSel_Edge > 5)', ht200 = '(htJet35j_Edge > 200)', ht500 = '(htJet35j_Edge > 500)')
            # negative values are not inside the selection. they should fall on the underflow of the histogram
            self.srIDMax = 23
            self.shortLabels = { 0: 'bveto_23jets_lowmet',
                                 1: 'bveto_45jets_lowmet',
                                 2: 'bveto_6jets_lowmet',
                                 3: 'withb_23jets_lowmet',
                                 4: 'withb_45jets_lowmet',
                                 5: 'withb_6jets_lowmet',
                                 6: 'bveto_23jets_medmet',
                                 7: 'bveto_45jets_medmet',
                                 8: 'bveto_6jets_medmet',
                                 9: 'withb_23jets_medmet',
                                10: 'withb_45jets_medmet',
                                11: 'withb_6jets_medmet',
                                12: 'bveto_23jets_highmet',
                                13: 'bveto_45jets_highmet',
                                14: 'bveto_6jets_highmet',
                                15: 'withb_23jets_highmet',
                                16: 'withb_45jets_highmet',
                                17: 'withb_6jets_highmet',
                                18: 'bveto_23jets_vhighmet',
                                19: 'bveto_45jets_vhighmet',
                                20: 'bveto_6jets_vhighmet',
                                21: 'withb_23jets_vhighmet',
                                22: 'withb_45jets_vhighmet',
                                23: 'withb_6jets_vhighmet'
                                 }
            self.SRLabels    = { 0: 'bveto_23jets_lowmet',
                                 1: 'bveto_45jets_lowmet',
                                 2: 'bveto_6jets_lowmet',
                                 3: 'withb_23jets_lowmet',
                                 4: 'withb_45jets_lowmet',
                                 5: 'withb_6jets_lowmet',
                                 6: 'bveto_23jets_medmet',
                                 7: 'bveto_45jets_medmet',
                                 8: 'bveto_6jets_medmet',
                                 9: 'withb_23jets_medmet',
                                10: 'withb_45jets_medmet',
                                11: 'withb_6jets_medmet',
                                12: 'bveto_23jets_highmet',
                                13: 'bveto_45jets_highmet',
                                14: 'bveto_6jets_highmet',
                                15: 'withb_23jets_highmet',
                                16: 'withb_45jets_highmet',
                                17: 'withb_6jets_highmet',
                                18: 'bveto_23jets_vhighmet',
                                19: 'bveto_45jets_vhighmet',
                                20: 'bveto_6jets_vhighmet',
                                21: 'withb_23jets_vhighmet',
                                22: 'withb_45jets_vhighmet',
                                23: 'withb_6jets_vhighmet'
                                 }
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
        ## expected limit -1 sigma(exp)
        self.ex_exp_m1s = copy.deepcopy(self.yx)
        self.ex_exp_m1s.SetTitle("nominal exclusion plot expected - 1 sigma exp"); self.ex_exp_m1s.SetName('ex_exp_m1s')
        self.ex_exp_m1s.Reset()

        ## expected limit +2 sigma(exp)
        self.ex_exp_p2s = copy.deepcopy(self.yx)
        self.ex_exp_p2s.SetTitle("nominal exclusion plot expected + 2 sigma exp"); self.ex_exp_p2s.SetName('ex_exp_p2s')
        self.ex_exp_p2s.Reset()
        ## expected limit -2 sigma(exp)
        self.ex_exp_m2s = copy.deepcopy(self.yx)
        self.ex_exp_m2s.SetTitle("nominal exclusion plot expected - 2 sigma exp"); self.ex_exp_m2s.SetName('ex_exp_m2s')
        self.ex_exp_m2s.Reset()
    
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
            elif 0.97 < point.quantileExpected < 0.98:
                self.ex_exp_m2s.Fill(massx, massy, point.limit)
            elif 0.02 < point.quantileExpected < 0.03:
                self.ex_exp_p2s.Fill(massx, massy, point.limit)

        zmax = self.ex_obs.GetMaximum()
        self.ex_obs    .GetZaxis().SetRangeUser(0.,10.)
        self.ex_obs_p1s.GetZaxis().SetRangeUser(0.,10.)
        self.ex_obs_m1s.GetZaxis().SetRangeUser(0.,10.)
        self.ex_exp    .GetZaxis().SetRangeUser(0.,10.)
        self.ex_exp_p1s.GetZaxis().SetRangeUser(0.,10.)
        self.ex_exp_m1s.GetZaxis().SetRangeUser(0.,10.)
        self.ex_exp_p2s.GetZaxis().SetRangeUser(0.,10.)
        self.ex_exp_m2s.GetZaxis().SetRangeUser(0.,10.)
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
        if self.doTwoSigmas: 
            self.ex_exp_p2s_smoothed_graph.Write()
            self.ex_exp_m2s_smoothed_graph.Write()
        f.Close()
        print '... done saving the results file'
    
    def makePrettyPlots(self):
        print 'starting to make pretty plots'
        for t in ['_obs', '_obs_p1s', '_obs_m1s', '_exp', '_exp_p1s', '_exp_m1s', '_exp_p2s', '_exp_m2s']:
            if self.doTwoSigmas and '2' in t: continue
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

