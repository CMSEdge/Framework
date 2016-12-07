import ROOT, math, optparse, copy

import include.CutManager as CutManager
import include.Sample     as Sample

import itertools


def makeFramesNice(frames):
    for frame in frames:
        tit = 0
        legcoords = [0.5, 0.7, 0.88, 0.85]
        if   'lepsDPhi' in frame.GetName():
            name = 'ldp_frame_nice'
            leg = ROOT.TLegend(legcoords[0], legcoords[1], legcoords[2], legcoords[3])
            leg.SetNColumns(2)
            leg.AddEntry(frame.findObject('h_em_data_cuts_of_sr_met150_of'), 'data  (OF) ', 'LP')
            leg.AddEntry(frame.findObject('ldp_analyticalPDF_DA_Norm[lepsDPhi_Edge]'), 'fit data (OF)', 'L')
            leg.AddEntry(frame.findObject('h_tt_mc_cuts_of_sr_met150_of'  ), 't#bar{t}  (OF) ', 'LP')
            leg.AddEntry(frame.findObject('ldp_analyticalPDF_MC_Norm[lepsDPhi_Edge]'), 'fit t#bar{t} (OF)', 'L')
            leg.AddEntry(frame.findObject('h_tt_mc_cuts_sf_sr_met150_sf'  ), 't#bar{t}  (SF) ', 'LP')
            leg.AddEntry(frame.findObject('ldp_analyticalPDF_MC_SF_Norm[lepsDPhi_Edge]'), 'fit t#bar{t} (SF)', 'L')
            tit  = 'pdf of leptons #Delta#phi'
            atit = '#Delta#phi_{ll}'
        elif 'lepsZPt' in frame.GetName():
            name = 'zpt_frame_nice'
            leg = ROOT.TLegend(legcoords[0], legcoords[1], legcoords[2], legcoords[3])
            leg.SetNColumns(2)
            leg.AddEntry(frame.findObject('h_em_data_cuts_of_sr_met150_of'), 'data  (OF) ', 'LP')
            leg.AddEntry(frame.findObject('zpt_analyticalPDF_DA_Norm[lepsZPt_Edge]'), 'fit data (OF)', 'LP')
            leg.AddEntry(frame.findObject('h_tt_mc_cuts_of_sr_met150_of'  ), 't#bar{t}  (OF) ', 'LP')
            leg.AddEntry(frame.findObject('zpt_analyticalPDF_MC_Norm[lepsZPt_Edge]'), 'fit t#bar{t} (OF)', 'LP')
            leg.AddEntry(frame.findObject('h_tt_mc_cuts_sf_sr_met150_sf'  ), 't#bar{t}  (SF) ', 'LP')
            leg.AddEntry(frame.findObject('zpt_analyticalPDF_MC_SF_Norm[lepsZPt_Edge]'), 'fit t#bar{t} (SF)', 'LP')
            tit  = 'pdf of dilepton-p_{T}'
            atit = 'p_{T}^{ll} (GeV)'
        elif 'mlb' in frame.GetName():
            name = 'mlb_frame_nice'
            leg = ROOT.TLegend(legcoords[0], legcoords[1], legcoords[2], legcoords[3])
            leg.SetNColumns(2)
            leg.AddEntry(frame.findObject('h_em_data_cuts_of_sr_met150_of'), 'data  (OF) ', 'LP')
            leg.AddEntry(frame.findObject('mlb_analyticalPDF_DA_Norm[sum_mlb_Edge]'), 'fit data (OF)', 'LP')
            leg.AddEntry(frame.findObject('h_tt_mc_cuts_of_sr_met150_of'  ), 't#bar{t}  (OF) ', 'LP')
            leg.AddEntry(frame.findObject('mlb_analyticalPDF_MC_Norm[sum_mlb_Edge]'), 'fit t#bar{t} (OF)', 'LP')
            leg.AddEntry(frame.findObject('h_tt_mc_cuts_sf_sr_met150_sf'  ), 't#bar{t}  (SF) ', 'LP')
            leg.AddEntry(frame.findObject('mlb_analyticalPDF_MC_SF_Norm[sum_mlb_Edge]'), 'fit t#bar{t} (SF)', 'LP')
            tit  = 'pdf of sum-m_{lb}'
            atit = '#Sigma m_{lb} (GeV)'
        elif 'met' in frame.GetName():
            name = 'met_frame_nice'
            leg = ROOT.TLegend(legcoords[0], legcoords[1], legcoords[2], legcoords[3])
            leg.SetNColumns(2)
            leg.AddEntry(frame.findObject('h_em_data_cuts_of_sr_met150_of'), 'data  (OF) ', 'LP')
            leg.AddEntry(frame.findObject('met_analyticalPDF_DA_Norm[met_Edge]'), 'fit data (OF)', 'LP')
            leg.AddEntry(frame.findObject('h_tt_mc_cuts_of_sr_met150_of'  ), 't#bar{t}  (OF) ', 'LP')
            leg.AddEntry(frame.findObject('met_analyticalPDF_MC_Norm[met_Edge]'), 'fit t#bar{t} (OF)', 'LP')
            leg.AddEntry(frame.findObject('h_tt_mc_cuts_sf_sr_met150_sf'  ), 't#bar{t}  (SF) ', 'LP')
            leg.AddEntry(frame.findObject('met_analyticalPDF_MC_SF_Norm[met_Edge]'), 'fit t#bar{t} (SF)', 'LP')
            tit  = 'pdf of E_{T}^{miss}'
            atit = 'E_{T}^{miss} (GeV)'


        frame.SetAxisRange(0., 1.3*frame.findObject('h_em_data_cuts_of_sr_met150_of').GetYaxis().GetXmax(), 'Y')
        frame.GetXaxis().SetTitle(atit)
        frame.GetXaxis().SetLabelSize(0.04)
        frame.GetYaxis().SetLabelSize(0.04)
        frame.SetTitle('')
        canv_tmp = ROOT.TCanvas('foo', 'bar', 600, 600)
        canv_tmp.SetLeftMargin(0.13)
        canv_tmp.SetBottomMargin(0.13)
        frame.GetYaxis().SetTitleOffset(1.5)
        frame.GetXaxis().SetTitleOffset(1.3)
        frame.Draw('same')
        leg.SetBorderSize(0)
        leg.Draw()
        lat = ROOT.TLatex()
        lat.SetNDC()
        lat.SetTextSize(0.040)
        lat.SetTextFont(61)
        lat.DrawLatex(0.16, 0.84, "CMS")
        lat.SetTextFont(52)
        lat.SetTextSize(0.035)
        lat.DrawLatex(0.16, 0.80, "Preliminary")

        lat.SetTextFont(42)
        lat.DrawLatex(0.63, 0.92, '12.9 fb^{-1} (13 TeV)')
        canv_tmp.SaveAs(name+'.png')
        canv_tmp.SaveAs(name+'.pdf')
        canv_tmp.SaveAs(name+'.root')


def cleanCut(cut):
    cut = cut.replace(cuts.twoLeptons, 't.nPairLep_Edge')
    cut = cut.replace('t.', '')
    return cut

class pdfClass:

    def __init__(self, _name, _ttree, _skimVal=1, isMC=False):
        self.name = _name
        self.tree = _ttree
        self.setTreeBranches()
        self.pdfs     = []
        self.ana_pdfs = []
        self.datasets = []
        self.frames   = []
        self.mcWgt = isMC if isMC else 1
        self.makeVarSet()
        self.skim = _skimVal
        ## setting some variables properties for the fit and the frame
        self.setVarProperties()


    def setTreeBranches(self):
        self.tree.SetBranchStatus('*'               , 0); ## old self.tree.SetBranchStatus('evt'             , 1);
        self.tree.SetBranchStatus('Lep1_pt_Edge'    , 1); self.tree.SetBranchStatus('Lep2_pt_Edge'    , 1);
        self.tree.SetBranchStatus('Lep1_eta_Edge'   , 1); self.tree.SetBranchStatus('Lep2_eta_Edge'   , 1);
        self.tree.SetBranchStatus('Lep1_pdgId_Edge' , 1); self.tree.SetBranchStatus('Lep2_pdgId_Edge' , 1);
        self.tree.SetBranchStatus('lepsDR_Edge'     , 1); self.tree.SetBranchStatus('lepsMll_Edge'    , 1);
        self.tree.SetBranchStatus('nJetSel_Edge'    , 1); self.tree.SetBranchStatus('nPairLep_Edge'   , 1);
        self.tree.SetBranchStatus('met_Edge'        , 1); self.tree.SetBranchStatus('min_mlb1_Edge'   , 1); 
        self.tree.SetBranchStatus('sum_mlb_Edge'    , 1); self.tree.SetBranchStatus('st_Edge'         , 1);
        self.tree.SetBranchStatus('lepsZPt_Edge'    , 1); self.tree.SetBranchStatus('d3D_Edge'        , 1);
        self.tree.SetBranchStatus('lepsDPhi_Edge'   , 1); self.tree.SetBranchStatus('j1MetDPhi_Edge'  , 1);
        self.tree.SetBranchStatus('j2MetDPhi_Edge'  , 1)
        
    def makeVarSet(self):
        self.l1p = ROOT.RooRealVar('Lep1_pt_Edge'    , 'l1p',   0., 1000., 'GeV'); self.l2p = ROOT.RooRealVar('Lep2_pt_Edge'    , 'l2p',   0., 1000., 'GeV');
        self.l1e = ROOT.RooRealVar('Lep1_eta_Edge'   , 'l1e', -2.5,   2.5, ''   ); self.l2e = ROOT.RooRealVar('Lep2_eta_Edge'   , 'l2e', -2.5,   2.5, ''   );
        self.l1i = ROOT.RooRealVar('Lep1_pdgId_Edge' , 'l1i',  -15,    15, ''   ); self.l2i = ROOT.RooRealVar('Lep2_pdgId_Edge' , 'l2i',  -15,    15, ''   );
        self.ldr = ROOT.RooRealVar('lepsDR_Edge'     , 'ldr',  0.1,  5.8 , ''   ); self.njs = ROOT.RooRealVar('nJetSel_Edge'    , 'njs',   0 ,   25 , ''   );
        self.nle = ROOT.RooRealVar('nPairLep_Edge'   , 'nle',   -2,    3 , ''   ); self.met = ROOT.RooRealVar('met_Edge'        , 'met',   0., 1000., 'GeV');
        self.mlp = ROOT.RooRealVar('metl1DPhi_Edge'  , 'mlp',   0.,3.142 , ''   ); self.mll = ROOT.RooRealVar('lepsMll_Edge'    , 'mll',  20., 1000., 'GeV');
        self.mlb = ROOT.RooRealVar('sum_mlb_Edge'    , 'mlb',   0., 3000., 'GeV'); self.st  = ROOT.RooRealVar('st_Edge'         , 'st' , 100., 5000., 'GeV');
        self.zpt = ROOT.RooRealVar('lepsZPt_Edge'    , 'zpt',   0., 1000., 'GeV'); ## old self.evt = ROOT.RooRealVar('evt'             , 'evt',   0.,1e7   , ''   );
        self.a3d = ROOT.RooRealVar('d3D_Edge'        , 'a3d', 0.2 , 3.14 , ''   ); self.dp1 = ROOT.RooRealVar('j1MetDPhi_Edge'  , 'dp1', -3.2, 3.2, '');
        self.ldp = ROOT.RooRealVar('lepsDPhi_Edge'   , 'ldp', 0.  , 3.14 , ''   ); self.dp2 = ROOT.RooRealVar('j2MetDPhi_Edge'  , 'dp2', -3.2, 3.2, '');
        self.nlt = ROOT.RooRealVar('nLepTight_Edge'  , 'nlt',  0 , 4 , ''   );
        self.wgt = ROOT.RooRealVar('lumwgt'          , 'wgt', self.mcWgt, self.mcWgt, ''); self.wgt.setConstant(1);

        self.varSet = ROOT.RooArgSet(self.l1p, self.l2p, self.l1e, self.l2e, self.l1i, self.l2i, self.ldr, self.njs, self.nle)
        self.varSet.add(self.met); self.varSet.add(self.mlp); self.varSet.add(self.mll); self.varSet.add(self.mlp);
        self.varSet.add(self.mlb); self.varSet.add(self.st); self.varSet.add(self.zpt); ## old self.varSet.add(self.evt);
        self.varSet.add(self.wgt); self.varSet.add(self.a3d); self.varSet.add(self.nlt); self.varSet.add(self.ldp)
        self.varSet.add(self.dp1); self.varSet.add(self.dp2)

    def makeDatasets(self, cutlist):
        for cut, name in cutlist.items():
            print 'loading a dataset into the class for cut %s'%(name)
            ## old cut = cut + ' && evt%{aa} == 0'.format(aa=self.skim)
            tmp_ds = ROOT.RooDataSet('ds_'+name, 'ds_'+name, self.tree, self.varSet, cut, 'lumwgt')
            setattr(self, name, tmp_ds)
            self.datasets.append(tmp_ds)

    def setVarProperties(self, var=0, _min=0, _max=0, rho=0):
        if not var:
            setattr(self, 'mlb_min',   0.); setattr(self, 'mlb_max', 3000.); setattr(self, 'mlb_rho', 1.4);
            setattr(self, 'met_min', 150.); setattr(self, 'met_max', 1000.); setattr(self, 'met_rho', 1.5);
            setattr(self, 'ldr_min',  0.3); setattr(self, 'ldr_max',   5.8); setattr(self, 'ldr_rho', 2.0);
            setattr(self, 'zpt_min',   0.); setattr(self, 'zpt_max', 1000.); setattr(self, 'zpt_rho', 1.4);
            setattr(self, 'a3d_min',  0.2); setattr(self, 'a3d_max',  3.14); setattr(self, 'a3d_rho', 1.4);
            setattr(self, 'st_min' , 140.); setattr(self, 'st_max' , 5000.); setattr(self, 'st_rho' , 1.8);
            setattr(self, 'ldp_min',   0.); setattr(self, 'ldp_max',  3.14); setattr(self, 'ldp_rho', 1.8);
            self.propSet = True
        else:
            setattr(self, var+'_min', _min); setattr(self, var+'_max', _max); setattr(self, var+'_rho', rho);

    def addNDPDF(self, variable, dataset, option, alpha):
        ds = getattr(self, dataset)
        pdf_name = 'pdf_%s_%s'%(variable, ds.GetName())
        pdf_histo_name = 'pdf_histo_%s_%s'%(variable, ds.GetName())
        print 'getting NDKeysPDF for variable %s in dataset %s of sample %s' %(getattr(self,variable).getTitle(), ds.GetName(), self.name)
        self.setRange(variable)
        ret_pdf = ROOT.RooNDKeysPdf(pdf_name, pdf_name, ROOT.RooArgList(getattr(self,variable)), ds, option, alpha)
        setattr(self, pdf_name, ret_pdf)
        pdf_histo = ret_pdf.createHistogram(pdf_histo_name, getattr(self,variable), ROOT.RooFit.Binning(1000))
        pdf_histo.SetName(self.name+'_'+pdf_histo_name)
        pdf_histo.Scale(1./pdf_histo.Integral())
        setattr(self, pdf_histo_name, pdf_histo)
        self.pdfs.append(ret_pdf)
        return ret_pdf

    def setRange(self, variable):
        getattr(self, variable).setMin(getattr(self, variable+'_min'))
        getattr(self, variable).setMax(getattr(self, variable+'_max'))
        #if   variable == 'ldr':
        #    getattr(self, variable).setMin(0.); getattr(self, variable).setMax(5.8)
        #elif variable == 'mlb' or variable == 'met':
        #    getattr(self, variable).setMin(0.); getattr(self, variable).setMax(750)
        #elif variable == 'zpt':
        #    getattr(self, variable).setMin(0.); getattr(self, variable).setMax(600)
        #elif variable == 'st':
        #    getattr(self, variable).setMin(100); getattr(self, variable).setMax(1000)
        

    def makeFrame(self, variable, _min=0, _max=0, _nbins=100):
        if _min != _max:
            tmp_frame = getattr(self, variable).frame(_min, _max, _nbins)
        else:
            tmp_frame = getattr(self, variable).frame(getattr(self, variable+'_min'), getattr(self, variable+'_max'), _nbins)
        for i,ds in enumerate(self.datasets):
            ds.plotOn(tmp_frame, ROOT.RooFit.MarkerColor(1+i), ROOT.RooFit.MarkerSize(0.8))
            if hasattr(self, 'pdf_%s_%s'%(variable, ds.GetName()) ):
                getattr(self, 'pdf_%s_%s'%(variable, ds.GetName()) ).plotOn(tmp_frame, ROOT.RooFit.LineColor(1+i)) 
        self.frames.append(tmp_frame)
        setattr(self, 'frame_%s'%(variable), tmp_frame)

    def saveInFile(self):
        for obj in dir(self):
            if not obj.startswith('pdf_histo'): continue
            getattr(self, obj).Write()

def getFrame(var, xmin=0, xmax=0):
    ## old if   var == 'met': vartree = 'met_pt'
    if   var == 'met': vartree = 'met_Edge'
    elif var == 'mlb': vartree = "sum_mlb_Edge";
    elif var == 'ldr': vartree = 'lepsDR_Edge';
    elif var == 'zpt': vartree = 'lepsZPt_Edge';
    elif var == 'a3d': vartree = 'd3D_Edge';
    elif var == 'ldp': vartree = 'lepsDPhi_Edge';
    if xmax <= xmin:
        return w.var(vartree).frame()
    else:
        return w.var(vartree).frame(xmin,xmax,100)

def buildModels(w):
    for var,t in itertools.product(['a3d', 'met', 'mlb', 'ldr', 'zpt', 'ldp'],['_DA', '_MC','_MC_SF']):
        if var == 'met':
            ## old vartree = 'met_pt'
            vartree = 'met_Edge'
            w.factory('met_c0{t}[-0.028,-1,0.]'.format(t=t))
            w.factory('met_c1{t}[-0.01,-0.05,-0.005]'.format(t=t))
            w.var(vartree).setMin(150.); w.var(vartree).setMax(1000.);
            w.factory('SUM::{var}_analyticalPDF{t}(met_n1{t}[0,1e4]*Exponential::exp1{t}({vartree},met_c0{t}),met_n2{t}[0,1e4]*Exponential::exp2{t}({vartree},met_c1{t}))'.format(var=var, vartree=vartree, t=t))
        elif var == 'mlb':
            vartree = "sum_mlb_Edge";
            w.factory("mlb_peak{t}[170,120,220]".format(t=t))
            w.factory("mlb_sigma{t}[40.,1.,100.]".format(t=t))
            w.factory("mlb_alpha{t}[-1.,-2.5,-0.5]".format(t=t))
            w.factory("mlb_n{t}[1.5,0.,2.5]".format(t=t))
            w.var(vartree).setMin(0.); w.var(vartree).setMax(3000.);
            w.factory("CBShape::{var}_analyticalPDF{t}({vartree},mlb_peak{t},mlb_sigma{t},mlb_alpha{t},mlb_n{t})".format(var=var, vartree=vartree,t=t))
        elif var == 'ldr':
            vartree = 'lepsDR_Edge';
            w.factory("ldr_peak{t}[3.0,2.0,3.5]".format(t=t))
            w.factory("ldr_sigma{t}[3.,0.,10.]".format(t=t))
            w.factory("ldr_alpha{t}[1.,0.,4.5]".format(t=t))
            w.factory("ldr_n{t}[2.5,1.,15.0]".format(t=t))
            #w.factory('CBShape::{var}_analyticalPDF({vartree},ldr_peak,ldr_sigma,ldr_alpha,ldr_n)'.format(var=var, vartree=vartree))
            w.factory('SUM::{var}_analyticalPDF{t}(ldr_n1{t}[0,1]*CBShape::{var}_cb{t}({vartree},ldr_peak{t},ldr_sigma{t},ldr_alpha{t},ldr_n{t}),ldr_n2{t}[0,1]*Gaussian::{var}_gaus{t}({vartree},ldr_gm{t}[0,-3,3],ldr_gs{t}[0.5,0.,0.8]))'.format(var=var, vartree=vartree, t=t))
        elif var == 'zpt':
            vartree = 'lepsZPt_Edge';
            w.factory("zpt_peak{t}[60.,30.,100.]".format(t=t))
            w.factory("zpt_sigma{t}[25., 5.,80.]".format(t=t))
            w.factory("zpt_alpha{t}[-1.,-2.5,0.]".format(t=t))
            w.factory("zpt_n{t}[5.0,5.0,100.0]".format(t=t))
            #w.factory('CBShape::{var}_analyticalPDF({vartree},zpt_peak,zpt_sigma,zpt_alpha,zpt_n)'.format(var=var, vartree=vartree))
            w.factory('SUM::{var}_analyticalPDF{t}(zpt{t}[0,1]*CBShape::{var}_cb{t}({vartree},zpt_peak{t},zpt_sigma{t},zpt_alpha{t},zpt_n{t}))'.format(var=var, vartree=vartree, t=t))
        elif var == 'a3d':
            vartree = 'd3D_Edge'
            ## w.factory("a3d_offset[1.5,1.4,1.6]")
            ## w.factory("a3d_a0[0.,-10.,10.]")
            ## w.factory("a3d_a1[0.,-10.,10.]")
            ## w.factory("a3d_a2[-1.,-5.,0.]")
            w.factory('EXPR::{var}_analyticalPDF{t}("a3d_a2{t}*({vartree}-a3d_offset{t})*({vartree}-a3d_offset{t}) + a3d_a1{t}*({vartree}-a3d_offset{t}) + a3d_a0{t}",{{{vartree},a3d_offset{t}[1.5,1.2,1.8],a3d_a0{t}[3.,2.5,5.],a3d_a1{t}[0.,-10,10],a3d_a2{t}[-1.,-3.,-1.]}})'.format(var=var, vartree=vartree,t=t))
            ## double gaussian
            ## w.factory("a3d_mean[1.5,1.2,1.7]")
            ## w.factory("a3d_sigma[1.,0.,3.]")
            ## #w.factory('Gaussian::{var}_gauss({vartree},a3d_mean,a3d_sigma)'.format(var=var, vartree=vartree))
            ## w.factory("a3d_mean2[1.5,1.5,2.0]")
            ## w.factory("a3d_sigma2[1.,0.,3.]")
            ## w.factory('Gaussian::{var}_gauss({vartree},a3d_mean,a3d_sigma)'.format(var=var, vartree=vartree))
            ## w.factory('Gaussian::{var}_gauss2({vartree},a3d_mean2,a3d_sigma2)'.format(var=var, vartree=vartree))
            ## w.factory('SUM::{var}_analyticalPDF(n11[0,1000]*a3d_gauss2,n22[0,1000]*a3d_gauss)'.format(var=var, vartree=vartree))

            ##double sided crystal ball approach (not very good)
            ## w.factory("a3d_mean[1.5,1.0,5.0]")
            ## w.factory("a3d_sigma[1.2,1.0,2.0]")
            ## w.factory("a3d_alpha1[1.,0.,2.]")
            ## w.factory("a3d_alpha2[1.,0.,2.]")
            ## w.factory("a3d_n1[1.,0.,10.]")
            ## w.factory("a3d_n2[1.,0.,10.]")
            ## w.factory("DoubleCB::{var}_analyticalPDF({vartree},a3d_mean,a3d_sigma,a3d_alpha1,a3d_n1,a3d_alpha2,a3d_n2)".format(var=var, vartree=vartree))
        elif var == 'ldp':
            vartree = 'lepsDPhi_Edge'
            w.factory('EXPR::{var}_analyticalPDF{t}("ldp_a2{t}*({vartree})*({vartree}) + ldp_a0{t}",{{{vartree},ldp_a0{t}[0.5,0,2],ldp_a2{t}[0.8,0.,2.]}})'.format(var=var, vartree=vartree,t=t))
            ## w.factory("ldp_lambda{t}[-3.,3.,-1.]".format(t=t))
            ## w.factory('SUM::{var}_analyticalPDF{t}(ldp{t}[0,1]*RooExponential::{var}_cb{t}({vartree},ldp_lambda{t}))'.format(var=var, vartree=vartree, t=t))

if __name__ == '__main__':

    print 'Starting pdf producer...'
    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-s', '--skim' , action='store', type=int , dest='skimVal', default=1, help='take only every nth event into the samples.')
    parser.add_option('-t', '--dott' , action='store_true', dest='dott' , help='do ttbar MC as well (takes time). default is 0')
    parser.add_option('-b', '--batch', action='store_true', dest='batch', help='set batch mode to not open windows. default false.')
    (opts, args) = parser.parse_args()

    if opts.batch:
        ROOT.gROOT.SetBatch()

    print 'running with these options \n'
    for key, value in opts.__dict__.items():
        print '%-20s : %-20s' %(key, value)

    cuts = CutManager.CutManager()
#    cuts_sr_met150_of = cuts.AddList([cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.OF]); cuts_sr_met150_of = cleanCut(cuts_sr_met150_of)
    cuts_sr_met150_of = cuts.AddList([cuts.BaselineNoTrigger, cuts.OF, cuts.MET150]); cuts_sr_met150_of = cleanCut(cuts_sr_met150_of)
#    cuts_sr_met150_sf = cuts.AddList([cuts.goodLepton, cuts.SignalRegionBaseLineNoTrigger, cuts.SF]); cuts_sr_met150_sf = cleanCut(cuts_sr_met150_sf)
    cuts_sr_met150_sf = cuts.AddList([cuts.BaselineNoTrigger, cuts.SF, cuts.MET150]); cuts_sr_met150_sf = cleanCut(cuts_sr_met150_sf)
    
    ## 2016 data and MC
#    tt_tfile = ROOT.TFile('/afs/cern.ch/work/m/mdunser/public/edgeTrees/trees_80X_ICHEP/mc_jun17_miniaodv2_noLHE/friends/evVarFriend_TTJets_DiLepton_total.root')
#    tt_tree = tt_tfile.Get('sf/t')
    #em_tfile = ROOT.TFile('/afs/cern.ch/work/m/mdunser/public/edgeTrees/trees_80X_ICHEP/data_jun20_prompt/friends/evVarFriend_MuonEG_Run2016B-PromptReco-v2.root')
#    em_tfile = ROOT.TFile('/afs/cern.ch/work/m/mdunser/public/edgeTrees/trees_80X_ICHEP/data_jun23_prompt_4invfb/friends/evVarFriend_MuonEG_Run2016B-PromptReco-v2_runs_271036_275125.root')
#    em_tree = em_tfile.Get('sf/t')
    em_tree = ROOT.TChain('sf/t')
    em_tree.Add('/afs/cern.ch/work/s/sesanche/public/ntuples/Data/evVarFriend_MuonEG_Run2016B_23Sep2016_v3_runs_273150_275376_part1.root')
    em_tree.Add('/afs/cern.ch/work/s/sesanche/public/ntuples/Data/evVarFriend_MuonEG_Run2016B_23Sep2016_v3_runs_273150_275376_part2.root')
    em_tree.Add('/afs/cern.ch/work/s/sesanche/public/ntuples/Data/evVarFriend_MuonEG_Run2016B_23Sep2016_v3_runs_273150_275376_part3.root')
    em_tree.Add('/afs/cern.ch/work/s/sesanche/public/ntuples/Data/evVarFriend_MuonEG_Run2016B_23Sep2016_v3_runs_273150_275376_part4.root')
    em_tree.Add('/afs/cern.ch/work/s/sesanche/public/ntuples/Data/evVarFriend_MuonEG_Run2016C_23Sep2016_v1_runs_271036_284044.root')
    em_tree.Add('/afs/cern.ch/work/s/sesanche/public/ntuples/Data/evVarFriend_MuonEG_Run2016D_23Sep2016_v1_runs_271036_284044.root')
    em_tree.Add('/afs/cern.ch/work/s/sesanche/public/ntuples/Data/evVarFriend_MuonEG_Run2016E_23Sep2016_v1_runs_271036_284044.root')
    em_tree.Add('/afs/cern.ch/work/s/sesanche/public/ntuples/Data/evVarFriend_MuonEG_Run2016H-PromptReco-v2_runs_281613_284035.root')
    em_tree.Add('/afs/cern.ch/work/s/sesanche/public/ntuples/Data/evVarFriend_MuonEG_Run2016H-PromptReco-v3_runs_284036_284044.root')


   
    dss = []
    em_data = pdfClass('em_data', em_tree)
    dss.append(em_data)
    if opts.dott: 
        tt_mc   = pdfClass('tt_mc'   , tt_tree, opts.skimVal,    (1./1.25)*12.9*831760.0/(2.58500e+08)) #13.5*87314.8*2.6/(1.98059e+08))#1.10523e+08)) ## not sure why this normalization does not work
        dss.append(tt_mc)
    
    cutsAndDs = {#cuts_sr       : 'cuts_of_sr'       ,
                 #cuts_sr_met100: 'cuts_of_sr_met100'}
                 cuts_sr_met150_of: 'cuts_of_sr_met150_of',
                 cuts_sr_met150_sf: 'cuts_sf_sr_met150_sf'}
    w = ROOT.RooWorkspace('w')
    for ds in reversed(dss):
        ds.makeDatasets(cutsAndDs)
        getattr(w,'import')(getattr(ds, 'cuts_of_sr_met150_of'),ROOT.RooFit.Rename(ds.name+'_cuts_of_sr_met150_of'))
        if not 'em_data' in ds.name:
            getattr(w,'import')(getattr(ds, 'cuts_sf_sr_met150_sf'),ROOT.RooFit.Rename(ds.name+'_cuts_sf_sr_met150_sf'))
    
    buildModels(w)

    doNDKeys = False
    
    frs = []; pdf_histos = []; frames = []
    writeobjs = []
    for var in ['ldp', 'met', 'mlb', 'zpt']:
        print '====================================================================='
        print '====================================================================='
        print '====== AT VARIABLE %s ================================================'%(var)
        print '====================================================================='
        print '====================================================================='
        opt = 'a' if var != 'met' else 'am'
        fmin = (150 if var == 'met' else 0 if var == 'zpt' else 0 if var == 'mlb' else 0)
        fmax = (500 if var == 'met' else 600 if var == 'zpt' else 1000 if var == 'mlb' else 0)
        tmp_frame = getFrame(var, fmin, fmax); 
        dslist = w.allData()
        for ds in dslist:
        #for i, ds in enumerate(reversed(dss)):
            model = var+'_analyticalPDF'
            #if ds.name == 'tt_mc'  : model += '_MC'
            #if ds.name == 'em_data': model += '_DA'
            if 'tt_mc'   in ds.GetName() and '_of' in ds.GetName(): model += '_MC'    ; color = ROOT.kRed   ; linestyle = ROOT.kSolid
            if 'tt_mc'   in ds.GetName() and '_sf' in ds.GetName(): model += '_MC_SF' ; color = ROOT.kBlue  ; linestyle = ROOT.kSolid
            if 'em_data' in ds.GetName() and '_of' in ds.GetName(): model += '_DA'    ; color = ROOT.kBlack ; linestyle = ROOT.kSolid
            #if 'em_data' in ds.GetName() and '_sf' in ds.GetName(): continue
            print '====================================================================='
            print '====================================================================='
            print '====== AT DATASET %s ================================================'%(ds.GetName())
            print '====================================================================='
            print '====================================================================='
            #color = ROOT.kBlue if i==0 else ROOT.kRed if i == 1 else ROOT.kBlack
            #linestyle = ROOT.kSolid if i == 0 else ROOT.kDashed
            ##if doNDKeys: ds.addNDPDF(var, 'cuts_of_sr_met150', opt, getattr(ds, var+'_rho') )

            #tmp_fr=w.pdf(model).fitTo(w.data(ds.name),ROOT.RooFit.Verbose(1),ROOT.RooFit.PrintLevel(1),ROOT.RooFit.NumCPU(1,0),ROOT.RooFit.Save(1),ROOT.RooFit.Strategy(2))
            tmp_fr=w.pdf(model).fitTo(w.data(ds.GetName()),ROOT.RooFit.Verbose(1),ROOT.RooFit.PrintLevel(1),ROOT.RooFit.NumCPU(1,0),ROOT.RooFit.Save(1),ROOT.RooFit.Strategy(2))
            #tmp_fr=w.pdf(model).fitTo(w.data(ds.name),ROOT.RooFit.Save(1))
            frs.append(copy.deepcopy(tmp_fr))

            #tmp_histo = w.pdf(model).createHistogram(ds.name+'_fit_histo_'+var+'_ds_cuts_of_sr_met150', getattr(ds,var), ROOT.RooFit.Binning(1000) )
            tmp_histo = w.pdf(model).createHistogram(ds.GetName()+'_fit_histo_'+var, getattr(dss[0],var), ROOT.RooFit.Binning(1000) )
            #tmp_histo.SetName(tmp_histo.GetName()[:tmp_histo.GetName().rfind('_ds_cuts_of_sr_met150')+len('_ds_cuts_of_sr_met150')] )
            tmp_histo.SetName(tmp_histo.GetName()[:tmp_histo.GetName().rfind('__')] )
            tmp_histo.SetLineColor(color); tmp_histo.SetLineStyle(linestyle)
            #ds.ana_pdfs.append(tmp_histo); pdf_histos.append(tmp_histo)

            #w.data(ds.name).plotOn(tmp_frame,ROOT.RooFit.Normalization(w.data('em_data').numEntries(),ROOT.RooAbsReal.NumEvent),ROOT.RooFit.MarkerColor(color), ROOT.RooFit.MarkerSize(0.8))
            w.data(ds.GetName()).plotOn(tmp_frame,ROOT.RooFit.MarkerColor(color), ROOT.RooFit.MarkerSize(0.8))
            #w.data(ds.GetName()).plotOn(tmp_frame,ROOT.RooFit.Normalization(w.data('em_data').numEntries(),ROOT.RooAbsReal.NumEvent),ROOT.RooFit.MarkerColor(color), ROOT.RooFit.MarkerSize(0.8))
            w.pdf(model).plotOn(tmp_frame,ROOT.RooFit.LineColor(color), ROOT.RooFit.LineStyle(linestyle) )
            #if doNDKeys: getattr(ds, 'pdf_%s_ds_cuts_of_sr_met150'%var).plotOn(tmp_frame,ROOT.RooFit.LineColor(color), ROOT.RooFit.LineStyle(ROOT.kDashed) )

            ## # damn root takes the objects away if i write them now
            writeobjs.append(copy.deepcopy(tmp_histo))
            writeobjs.append(copy.deepcopy(w.pdf(model)))
            #if doNDKeys: writeobjs.append(copy.deepcopy(getattr(ds, 'pdf_histo_%s_ds_%s'%(var, 'cuts_of_sr_met150') ) ))
            #a = raw_input()

        frames.append(tmp_frame);
        c = ROOT.TCanvas()
        tmp_frame.Draw()
        c.SaveAs('frame_%s.pdf' %(var))
        c.SaveAs('frame_%s.root'%(var))

    #makeFramesNice(copy.deepcopy(frames))
    

    w.writeToFile("pdfs/pdfs_version6_80X_2016Data_savingTheWorkspace_withSFPDFs_12p9invfb_forFrameICHEP.root")

    outfile = ROOT.TFile('pdfs/pdfs_version5_80X_2016Data_withSFPDFs_12p9invfb_forFrameICHEP.root', 'RECREATE')
    outfile.cd()
    for i in writeobjs:
        i.Write()
    for i in frs:
        i.Write()
    outfile.Close()


