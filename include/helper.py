import math, sys, os, copy, re
import ROOT as r
from   ROOT import TGraphErrors, gROOT, TCanvas, TFile





def ensurePath( _path):
    d = os.path.dirname(_path)
    if not os.path.exists(d):
        os.makedirs(d)             





def ratioError(num, num_e, den, den_e, opt=''):
    tmp_n = r.TH1F('n', 'n', 1,0,1)
    tmp_d = r.TH1F('d', 'd', 1,0,1)
    tmp_r = r.TH1F('r', 'r', 1,0,1)
    tmp_n.SetBinContent(1,num); tmp_n.SetBinError(1,num_e)
    tmp_d.SetBinContent(1,den); tmp_d.SetBinError(1,den_e)
    #tmp_n.Divide(tmp_d)
    tmp_r.Divide(tmp_n, tmp_d, 1., 1., opt)
    return tmp_r.GetBinContent(1), tmp_r.GetBinError(1)


def readFromFileEWKFactors(theFile, dataMC, factor, reg): 

    for line in open(theFile).readlines():
        if line.find(factor) != -1 and line.find(reg) != -1:
            if line.find(dataMC) != -1:
                splitline = line.split(" ")
                x = 0
                if dataMC == 'MC': x = 2
                r = splitline[19+x]
                stat = splitline[25+x]
                arr = []
                arr.append(r)
                arr.append(stat)
                print "using this ",factor , " ", reg, " : ",  r
                if arr[0] == '':
                    print "warning probably not right rsfof!!!", arr
                return map(float, arr)                                       


def readFromFileRsfofD(theFile, dataMC): 

    for line in open(theFile).readlines():
        if line.find("rsfof") != -1 and line.find("direct") != -1:
            if line.find(dataMC) != -1:
                splitline = line.split(" ")
                x = 0
                if dataMC == 'MC': x = 2
                rsfof = splitline[25+x]
                stat = splitline[31+x]
                syst = splitline[38+x]
                arr = []
                arr.append(rsfof)
                arr.append(stat)
                arr.append(syst)
                print "using these direct rsfofs ", dataMC, ": ", arr
                if arr[0] == '':
                    print "warning probably not right rsfof!!!", arr
                return map(float, arr)                                                          

def readFromFileRsfofDee(theFile, dataMC): 

    for line in open(theFile).readlines():
        if line.find("rsfofee") != -1 and line.find("direct") != -1:
            if line.find(dataMC) != -1:
                splitline = line.split(" ")
                x = 0
                if dataMC == 'MC': x = 2
                rsfof = splitline[25+x]
                stat = splitline[31+x]
                syst = splitline[38+x]
                arr = []
                arr.append(rsfof)
                arr.append(stat)
                arr.append(syst)
                print "using these direct rsfofs for ee", dataMC, ": ", arr
                if arr[0] == '':
                    print "warning probably not right rsfof!!!", arr
                return map(float, arr)                                          

def readFromFileRsfofDmm(theFile, dataMC): 

    for line in open(theFile).readlines():
        if line.find("rsfofmm") != -1 and line.find("direct") != -1:
            if line.find(dataMC) != -1:
                splitline = line.split(" ")
                x = 0
                if dataMC == 'MC': x = 2
                rsfof = splitline[25+x]
                stat = splitline[31+x]
                syst = splitline[38+x]
                arr = []
                arr.append(rsfof)
                arr.append(stat)
                arr.append(syst)
                print "using these direct rsfofs for mm", dataMC, ": ", arr
                if arr[0] == '':
                    print "warning probably not right rsfof!!!", arr
                return map(float, arr)                                          









def readFromFileRinout(theFile, dataMC, reg): 

    for line in open(theFile).readlines():
        if line.find("rinout") != -1 and line.find(reg) != -1:
            if line.find(dataMC) != -1:
                splitline = line.split(" ")
                x = 0
                if dataMC == 'MC': x = 2
                rinout = splitline[25+x]
                stat = splitline[31+x]
                syst = splitline[38+x]
                arr = []
                arr.append(rinout)
                arr.append(stat)
                arr.append(syst)
                print "using these rinouts ", dataMC, ": ", arr
                if arr[0] == '':
                    print "warning probably not right rinout!!!", arr
                return map(float, arr)                                           


def readFromFileRmueCoeff(theFile, coeff, dataMC): 

    for line in open(theFile).readlines():
        if line.find("rmue") != -1 and line.find(coeff) != -1:
            if line.find(dataMC) != -1:
                splitline = line.split(" ")
                rmue = splitline[27]
                stat = splitline[33]
                arr = []
                arr.append(rmue)
                arr.append(stat)
                arr.append('0.00')
                print "using these rmues coeff", dataMC, ": ", arr
                if arr[0] == '':
                    print "warning probably not right rmue!!!", arr
                return map(float, arr)                                        


def readFromFileRmue(theFile, dataMC): 

    for line in open(theFile).readlines():
        if line.find("rmue") != -1 and line.find("factor") != -1:
            if line.find(dataMC) != -1:
                splitline = line.split(" ")
                x = 0
                if dataMC == 'MC': x = 2
                rmue = splitline[27+x]
                stat = splitline[33+x]
                arr = []
                arr.append(rmue)
                arr.append(stat)
                arr.append('0.00')
                print "using these rmues factor", dataMC, ": ", arr
                if arr[0] == '':
                    print "warning probably not right rmue!!!", arr
                return map(float, arr)                                               

def readFromFileRT(theFile, dataMC): 

    for line in open(theFile).readlines():
        if line.find("rt") != -1 and line.find("region") != -1:
            if line.find(dataMC) != -1:
                splitline = line.split(" ")
                x = 0
                if dataMC == 'MC': x = 2
                rt   = splitline[28+x]
                stat = splitline[34+x]
                syst = splitline[41+x]
                arr = []
                arr.append(rt)
                arr.append(stat)
                arr.append(syst)
                print "using these rts ", dataMC, ": ", arr
                if arr[0] == '':
                    print "warning probably not right rt!!!", arr
                return map(float, arr)                                       


def makeRatioPlot(hlist, name, opts = {}):

    norm = False; do = 'pe'; legco = [0.6,0.5,0.88,0.9];
    if opts.has_key('norm' ): norm = opts['norm']
    if opts.has_key('do'   ): do   = opts['do']
    if opts.has_key('legco'): legco= opts['legco']

    c1 = r.TCanvas(); c1.cd()
    pad1 = r.TPad("pad1", "pad1", 0, 0.25, 1, 1.0)
    pad1.SetBottomMargin(0.1)
    pad1.Draw()
    pad2 = r.TPad("pad2", "pad2", 0, 0.05, 1, 0.25)
    pad2.SetTopMargin(0.1);
    pad2.SetBottomMargin(0.3);
    pad2.Draw();

    leg = r.TLegend(legco[0], legco[1], legco[2], legco[3])
    leg.SetTextSize(0.04)
    ymax = 0.
    for ih,h in enumerate(hlist):
        tmp_max = h.GetMaximum()
        if tmp_max > ymax: ymax = tmp_max

    pad1.cd()
    for ih,h in enumerate(hlist):
        if norm: h.Scale(1./h.Integral())
        else: h.GetYaxis().SetRangeUser(0.01, ymax*1.2)
        h.Draw(do + ('' if not ih else ' same'))
        #leg.AddEntry(h, ' '.join(h.GetName().replace('_norm','').replace('_',' ').split()[:-1]), 'pl')
        leg.AddEntry(h, ' '.join(h.GetName().replace('_norm','').replace('_',' ').replace('observed','obs.').replace('predicted','pred.').replace('estimated','pred.').replace('corrected','(corr.)')), 'pl')

    leg.Draw('same')

    if opts.has_key('fitratio'): 
        fitratio = opts['fitratio']
        print 'found ratio fit key. fitting index %d with function %s'%(fitratio[0], fitratio[1])
    else: fitratio = [-1 , ''] 
    r_fiti = fitratio[0]; r_fitf = fitratio[1]
    ratios = []
    pad2.cd()
    for ih,h in enumerate(hlist):
        if not ih: 
            refHisto = copy.deepcopy(h)
        else:
            divh = copy.deepcopy(h)
            divh.Divide(refHisto)
            divh.GetYaxis().SetRangeUser(0.0,2.0)
            divh.GetYaxis().SetTitle('ratio')
            divh.GetYaxis().CenterTitle()
            divh.GetYaxis().SetTitleSize(0.14)
            divh.GetYaxis().SetLabelSize(0.12)
            divh.GetYaxis().SetNdivisions(4)
            divh.GetXaxis().SetLabelSize(0.12)
            divh.GetXaxis().SetTitle('')
            if ih == r_fiti:
                divh.Fit(r_fitf,'','',50,250)
                retval = copy.deepcopy(divh.GetFunction(fitratio[1]))
            #divh.Draw('pe same')#+ '' if ih == 1 else ' same')
            ratios.append(divh)

    for ratio in ratios:
        ratio.Draw('e0 same')

    l1 = r.TLine(divh.GetXaxis().GetXmin(), 1., divh.GetXaxis().GetXmax(), 1.)
    l1.SetLineStyle(2)
    l1.Draw('same')

    c1.SaveAs(name+'.png')
    c1.SaveAs(name+'.pdf')
    c1.SaveAs(name+'.root')
    pad1.cd()
    pad1.SetLogy(1)
    leg.SetX1NDC(0.35); leg.SetY1NDC(0.25); leg.SetX2NDC(0.58); leg.SetY2NDC(0.55)
    c1.SaveAs(name+'_log.png')
    c1.SaveAs(name+'_log.pdf')
    c1.SaveAs(name+'_log.root')

    if r_fiti >= 0:
        return retval

class vParams:
    def __init__(self,name):
        self.name = name
        self.setParams()
    def setParams(self):
        if   self.name == 'met': self.xmin =   0   ; self.xmax = 300 ; self.nbins = 12; self.ymax = 0.35; self.title = 'ME_{T} (GeV)'         ; self.vtree = 'met_Edge'            ;
        elif self.name == 'mt2': self.xmin =   0   ; self.xmax = 200 ; self.nbins = 25; self.ymax = 0.35; self.title = 'MT_{2} (GeV)'         ; self.vtree = 'mt2_Edge'            ;
        elif self.name == 'hjj': self.xmin =  50   ; self.xmax = 300 ; self.nbins = 25; self.ymax = 0.35; self.title = 'm_{jj, hard}'         ; self.vtree = 'hardMjj_Edge'        ;
        elif self.name == 'ldr': self.xmin = 0.1   ; self.xmax = 0.5 ; self.nbins =  4; self.ymax = 0.08; self.title = '#Delta R_{ll}'        ; self.vtree = 'lepsDR_Edge'         ;
        elif self.name == 'mll': self.xmin = 11.   ; self.xmax = 121.; self.nbins = 22; self.ymax = 0.08; self.title = 'm_{ll}'               ; self.vtree = 'lepsMll_Edge'        ;
        elif self.name == 'ldphi': self.xmin = 0.  ; self.xmax = 3.15; self.nbins = 25; self.ymax = 0.08; self.title = '#Delta #phi_{ll}'     ; self.vtree = 'lepsDPhi_Edge'       ;
        elif self.name == 'zpt': self.xmin = 0.0   ; self.xmax = 250 ; self.nbins = 25; self.ymax = 0.24; self.title = 'p_{T}^{ll} (GeV)'     ; self.vtree = 'lepsZPt_Edge'        ;
        elif self.name == 'ht': self.xmin = 0.0    ; self.xmax = 250 ; self.nbins = 25; self.ymax = 0.24; self.title = 'H_{T} (GeV)'          ; self.vtree = 'htJet35j_Edge'       ;
        elif self.name == 'mlb': self.xmin = 0.0   ; self.xmax = 500 ; self.nbins = 25; self.ymax = 0.30; self.title = '#Sigma m_{lb} (Gev)'  ; self.vtree = 'sum_mlb_Edge'        ;
        elif self.name == 'a3d': self.xmin = 0.0   ; self.xmax = 3.15; self.nbins = 25; self.ymax = 0.14; self.title = '#Delta 3D'            ; self.vtree = 'd3D_Edge'            ;
        elif self.name == 'jjp': self.xmin = 0.0   ; self.xmax = 3.15; self.nbins = 25; self.ymax = 0.14; self.title = '#Delta #phi_{jj,hard}'; self.vtree = 'abs(hardJJDphi_Edge)';
        elif self.name == 'j1mphi': self.xmin = 0.0; self.xmax = 3.15; self.nbins = 25; self.ymax = 0.14; self.title = '#Delta #phi_{j1,met}' ; self.vtree = 'abs(j1MetDPhi_Edge)' ;
        elif self.name == 'j2mphi': self.xmin = 0.0; self.xmax = 3.15; self.nbins = 25; self.ymax = 0.14; self.title = '#Delta #phi_{j2,met}' ; self.vtree = 'abs(j2MetDPhi_Edge)' ;
        elif self.name == 'jjdr': self.xmin = 0.0  ; self.xmax = 5.00; self.nbins = 25; self.ymax = 0.14; self.title = '#Delta R_{j1,j2}'     ; self.vtree = 'abs(hardJJDR_Edge)'  ;
        elif self.name == 'l1pt': self.xmin = 20.  ; self.xmax = 270.; self.nbins = 25; self.ymax = 0.14; self.title = 'p_{T}^{l1}'           ; self.vtree = 'Lep1_pt_Edge'        ;
        elif self.name == 'l2pt': self.xmin = 20.  ; self.xmax = 270.; self.nbins = 25; self.ymax = 0.14; self.title = 'p_{T}^{l2}'           ; self.vtree = 'Lep2_pt_Edge'        ;
        elif self.name == 'l1eta': self.xmin = -3. ; self.xmax = 3.  ; self.nbins = 25; self.ymax = 0.14; self.title = '#eta^{l1}'            ; self.vtree = 'Lep1_eta_Edge'       ;
        elif self.name == 'l2eta': self.xmin = -3. ; self.xmax = 3.  ; self.nbins = 25; self.ymax = 0.14; self.title = '#eta^{l2}'            ; self.vtree = 'Lep2_eta_Edge'       ;
        elif self.name == 'nj': self.xmin =  0.    ; self.xmax = 3.  ; self.nbins =  3; self.ymax = 0.14; self.title = 'n_{jets}'             ; self.vtree = 'nJetSel_Edge'        ;
        elif self.name == 'nb': self.xmin =  0.    ; self.xmax = 3.  ; self.nbins =  3; self.ymax = 0.14; self.title = 'n_{b-jets}'           ; self.vtree = 'nBJetMedium35_Edge'  ;
        elif self.name == 'metBIN': self.xmin = 100; self.xmax = 600 ; self.nbins = 10; self.ymax = 0.35; self.title = 'ME_{T} (GeV)'         ; self.vtree = 'met_Edge'            ;
        elif self.name == 'nll': 
            self.xmin = 12  ; self.xmax =  30  ; self.nbins =  6 ; self.ymax = 0.20 ; self.title = 'NLL'                 ;
            self.vtree = '-1.*TMath::Log(lh_ana_met_data_Edge*lh_ana_mlb_data_Edge*lh_ana_a3d_data_Edge*lh_ana_zpt_data_Edge)'; 
        else:
            sys.exit('TRYING TO CALL VPARAMS WITH UNKNOWN VARIABLE')

def ensureDirectory(_path):
   #d = os.path.dirname(_path)
   #print d
   if not os.path.exists(_path):
      os.makedirs(_path)

def selectSamples(inputfile, selList, sType = 'DATA'):
    f = open(inputfile, 'r')
    tmp_file = open('.tmp_sampleFile%s.txt' %sType, 'w')
    checkedList = []
    typeList    = []
    for line in f.readlines():
        #print "line ", line
        if '#' in line or not len(line.rstrip('\r')): continue
        for _sample in selList:
            #print "_sample", _sample
            _sample = _sample.replace('*','.*')
            if _sample == line.split()[2] or re.search(_sample, line.split()[2]):
                if not sType == 'SYNCH':
                    tmp_file.write(line)
                else:
                    tmp_splitline = line.split()
                    tmp_splitline[0] = 'synching'
                    tmp_file.write('  '.join(tmp_splitline+['\n']))
                checkedList.append(_sample)
                typeList   .append(int(line.split()[-2]))
                #typeList   .append(int(line.split()[-1])) 
                #i changed this to 2!
    for _selSample in selList:
        if not '*' in _selSample:
            if _selSample not in checkedList:
                print 'ERROR: some samples weren\'t selected, check all sample names!'
                print 'check this sample:', _selSample
                sys.exit('exiting...')
        else:
            print 'you used some wildcards in selecting the samples. be careful with that!'
    print typeList
    if not len(set(typeList)) == 1:
            print 'ERROR: you\'re mixing DATA and MC!'
            sys.exit('exiting...')
            
    return tmp_file.name

class color:
       purple = '\033[95m'
       cyan = '\033[96m'
       darkcyan = '\033[36m'
       blue = '\033[94m'
       green = '\033[92m'
       yellow = '\033[93m'
       red = '\033[91m'
       bold = '\033[1m'
       underline = '\033[4m'
       end = '\033[0m'

class bcolors: 
    HEADER = '\033[95m' 
    OKBLUE = '\033[94m' 
    OKGREEN = '\033[92m' 
    WARNING = '\033[93m' 
    FAIL = '\033[91m' 
    ENDC = '\033[0m' 
    BOLD = '\033[1m' 
    UNDERLINE = '\033[4m' 

class valErrs:
    def __init__(self, cen_val, cen_sys, cen_stat, fwd_val, fwd_sys, fwd_stat, name):
        self.cen_val  = cen_val
        self.cen_sys  = cen_sys
        self.cen_stat = cen_stat
        self.fwd_val  = fwd_val
        self.fwd_sys  = fwd_sys
        self.fwd_stat = fwd_stat
        self.vals = []
        self.name = name
    def totError(self):
        self.cen_err = math.sqrt((self.cen_val*self.cen_sys)**2 + self.cen_stat**2)
        self.fwd_err = math.sqrt((self.fwd_val*self.fwd_sys)**2 + self.fwd_stat**2)
        self.vals.append(self.cen_err)
        self.vals.append(self.fwd_err)
    def setVals(self, line):
        self.cen_val  = float(line.split()[-6])
        self.cen_sys  = float(line.split()[-4])
        self.cen_stat = float(line.split()[-5])
        self.vals.extend([self.cen_val, self.cen_sys, self.cen_stat])
        self.fwd_val  = float(line.split()[-3])
        self.fwd_sys  = float(line.split()[-1])
        self.fwd_stat = float(line.split()[-2])
        self.vals.extend([self.fwd_val, self.fwd_sys, self.fwd_stat])
        self.totError()
    def printValues(self):
        print '%s central %.3f +- %.3f (%.3f stat. %.3f syst.)' %(
                self.name, self.cen_val, self.cen_err, self.cen_stat, self.cen_sys)
        print '%s forward %.3f +- %.3f (%.3f stat. %.3f syst.)' %(
                self.name, self.fwd_val, self.fwd_err, self.fwd_stat, self.fwd_sys)

class ingredients:
    def __init__(self, infile, dType):
        self.infile = infile
        self.isData = (dType == 'DATA')
        self.dType  = dType
        self.rs = []
        ## rmue
        self.rmue_alone   = valErrs(-1., -1., -1., -1., -1., -1., 'r_mue alone'  ); self.rs.append(self.rmue_alone  )
        self.rmue_factor  = valErrs(-1., -1., -1., -1., -1., -1., 'r_mue factor' ); self.rs.append(self.rmue_factor )
        ## rsfof                                   
        self.rsfof_direct = valErrs(-1., -1., -1., -1., -1., -1., 'R_sfof direct'); self.rs.append(self.rsfof_direct)
        ## RT
        self.rt_region    = valErrs(-1., -1., -1., -1., -1., -1., 'R_T'          ); self.rs.append(self.rt_region   )
        ## rinout                                  
        self.rinout_dy_lm = valErrs(-1., -1., -1., -1., -1., -1., 'R_inout SR lm'); self.rs.append(self.rinout_dy_lm)
        self.rinout_dy_bz = valErrs(-1., -1., -1., -1., -1., -1., 'R_inout SR oz'); self.rs.append(self.rinout_dy_bz)
        self.rinout_dy_oz = valErrs(-1., -1., -1., -1., -1., -1., 'R_inout SR bz'); self.rs.append(self.rinout_dy_oz)
        self.rinout_dy_az = valErrs(-1., -1., -1., -1., -1., -1., 'R_inout SR az'); self.rs.append(self.rinout_dy_az)
        self.rinout_dy_hm = valErrs(-1., -1., -1., -1., -1., -1., 'R_inout SR hm'); self.rs.append(self.rinout_dy_hm)
        ## fill the values from the file
        self.readValues()
        ## check if none is unset. exit if one is.
        self.checkValues()
        ## caluclate final rsfof as well
        self.calculateFullRSFOF()
        print 'loaded all ingredients from %s for %s' %(self.infile, self.dType)

    def check(self, test, string):
        return all(i in string for i in test)

    def readValues(self):
        print 'reading values from %s for %s' %(self.infile, self.dType)
        f = open(self.infile, 'r')
        lines = f.read().splitlines()
        for line in lines:
            if '#' in line or not len(line.strip()): continue
            #rmue
            if self.check(['rmue'  , 'alone' , self.dType], line): self.rmue_alone  .setVals(line)
            if self.check(['rmue'  , 'factor', self.dType], line): self.rmue_factor .setVals(line)
            #rsfof
            if self.check(['rsfof' , 'direct', self.dType], line): self.rsfof_direct.setVals(line)
            #RT
            if self.check(['rt'    , 'region', self.dType], line): self.rt_region   .setVals(line)
            #rinout
            if self.check(['rinout', 'dy_lm' , self.dType], line): self.rinout_dy_lm.setVals(line)
            if self.check(['rinout', 'dy_bz' , self.dType], line): self.rinout_dy_bz.setVals(line)
            if self.check(['rinout', 'dy_oz' , self.dType], line): self.rinout_dy_oz.setVals(line)
            if self.check(['rinout', 'dy_az' , self.dType], line): self.rinout_dy_az.setVals(line)
            if self.check(['rinout', 'dy_hm' , self.dType], line): self.rinout_dy_hm.setVals(line)
        f.close()

    def calculateFullRSFOF(self):
        ## ========================================
        ## calculate everything for central first.
        ## ========================================
        self.rsfof_fac_cen   = self.rmue_factor.cen_val * self.rt_region.cen_val
        self.rsfof_fac_cen_e = math.sqrt(self.rmue_factor.cen_err**2 + self.rt_region.cen_err**2)
        self.rsfof_fac_cen_w = 1./(self.rsfof_fac_cen_e**2)

        self.rsfof_dir_cen    = self.rsfof_direct.cen_val
        self.rsfof_dir_cen_e  = self.rsfof_direct.cen_err
        self.rsfof_dir_cen_w  = 1./(self.rsfof_dir_cen_e**2)

        ## calculate the numerator and denominator for the weighted average
        self.rsfof_final_cen_num = self.rsfof_fac_cen * self.rsfof_fac_cen_w + self.rsfof_dir_cen * self.rsfof_dir_cen_w
        self.rsfof_final_cen_den = self.rsfof_fac_cen_w + self.rsfof_dir_cen_w
        ## here come the final value and its uncertainty
        self.rsfof_final_cen_val = self.rsfof_final_cen_num / self.rsfof_final_cen_den
        self.rsfof_final_cen_err = 1./math.sqrt(self.rsfof_final_cen_den)

        ## ==========================================================
        ## repeat everything for forward. will do that nicer later.
        ## ==========================================================
        self.rsfof_fac_fwd   = self.rmue_factor.fwd_val * self.rt_region.fwd_val
        self.rsfof_fac_fwd_e = math.sqrt(self.rmue_factor.fwd_err**2 + self.rt_region.fwd_err**2)
        self.rsfof_fac_fwd_w = 1./(self.rsfof_fac_fwd_e**2)

        self.rsfof_dir_fwd    = self.rsfof_direct.fwd_val
        self.rsfof_dir_fwd_e  = self.rsfof_direct.fwd_err
        self.rsfof_dir_fwd_w  = 1./(self.rsfof_dir_fwd_e**2)

        ## calculate the numerator and denominator for the weighted average
        self.rsfof_final_fwd_num = self.rsfof_fac_fwd * self.rsfof_fac_fwd_w + self.rsfof_dir_fwd * self.rsfof_dir_fwd_w
        self.rsfof_final_fwd_den = self.rsfof_fac_fwd_w + self.rsfof_dir_fwd_w
        ## here come the final value and its uncertainty
        self.rsfof_final_fwd_val = self.rsfof_final_fwd_num / self.rsfof_final_fwd_den
        self.rsfof_final_fwd_err = 1./math.sqrt(self.rsfof_final_fwd_den)

        f=open(self.infile, 'r')
        newlines = []
        for line in list(f):
            if len(line.split()) > 2 and line.split()[0] == 'rsfof'  and line.split()[1] == 'final' and line.split()[2] == self.dType:
                newlines.append('rsfof       final           %-4s        %.4f      %.4f      %.4f      %.4f      %.4f      %.4f\n'%(
                    str(self.dType), self.rsfof_final_cen_val, self.rsfof_final_cen_err, 0.0, 
                                     self.rsfof_final_fwd_val, self.rsfof_final_fwd_err, 0.0 ))
            else:
                newlines.append(line)
        f.close()
        g=open(self.infile, 'w')
        g.writelines(newlines)
        g.close()


    def checkValues(self):
        for thing in self.rs:
            if any(i < 0 for i in thing.vals):
                print 'ERROR: some of the ingredients aren\'t set properly'
                print thing.printValues()
                #sys.exit('exiting...')



##################################################################################
### These functions were created by Aachen in order to define the proper style ###
##################################################################################
def createMyColors():
    iIndex = 2000

    containerMyColors = []
    for color in defineMyColors.keys():
       	tempColor = r.TColor(iIndex,
       	float(defineMyColors[color][0]) / 255, float(defineMyColors[color][1]) / 255, float(defineMyColors[color][2]) / 255)
       	containerMyColors.append(tempColor)

       	myColors.update({ color: iIndex })
       	iIndex += 1

    return containerMyColors


defineMyColors = {
        'Black' : (0, 0, 0),
        'White' : (255, 255, 255),
        'Red' : (255, 0, 0),
        'DarkRed' : (128, 0, 0),
        'Green' : (0, 255, 0),
        'Blue' : (0, 0, 255),
        'Yellow' : (255, 255, 0),
        'Orange' : (255, 128, 0),
        'DarkOrange' : (255, 64, 0),
        'Magenta' : (255, 0, 255),
        'KDEBlue' : (64, 137, 210),
        'Grey' : (128, 128, 128),
        'DarkGreen' : (0, 128, 0),
        'DarkSlateBlue' : (72, 61, 139),
        'Brown' : (70, 35, 10),

        'MyBlue' : (36, 72, 206),
        'MyDarkBlue' : (18, 36, 103),
        'MyGreen' : (70, 164, 60),
        'AnnBlueTitle' : (29, 47, 126),
        'AnnBlue' : (55, 100, 255),
#        'W11AnnBlue' : (0, 68, 204),
#        'W11AnnBlue' : (63, 122, 240),
    }


myColors = {
            'W11ttbar':  855,
            'W11singlet':  854,
            'W11ZLightJets':  401,
            'W11ZbJets':  400,
            'W11WJets':  842,
            'W11Diboson':  920,
            'W11AnnBlue': 856,
            'W11Rare':  630,
            }


