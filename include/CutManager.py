

class CutManager:
   'This class serves as an on-demand cut server'

   def __init__(self):

      self.protection = '((run_Edge<=276811)  ||  (278820<=run_Edge && run_Edge<=279931))'

      ########################################################################
      ######Basic Lepton Cuts ################################################
      ########################################################################
      self.twoLeptons = "nPairLep_Edge > 0"
      self.filters = "Filters"
      self.tightCharge = 'Lep1_tightCharge_Edge > 0 && Lep2_tightCharge_Edge > 0'
      self.leptonPt = "Lep1_pt_Edge > 25 && Lep2_pt_Edge > 20."
      self.diLeptonPt = "lepsZPt_Edge > 25"
      self.leptonDR = "lepsDR_Edge > 0.1"       
      self.leptonDR4l = "lepsDR_Edge > 0.02"       
      self.leptonsMll = "lepsMll_Edge > 20"
      self.ee = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -121)"
      self.mm = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -169)"
      self.OF = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -143)"
      self.SF = "(" + self.ee + " || " +  self.mm + ")"
      self.AF = "(" + self.SF + " || " +  self.OF + ")"
      self.trigMM = "Trigger_2m"
      self.trigEE = "Trigger_2e"
      self.trigEM = "Trigger_em"
      self.trigger     = "((" + self.trigMM + " && " + self.mm + ") || (" + self.trigEE + " && " + self.ee + ") || (" + self.trigEM + " && " + self.OF + "))"
      self.goodLepton = "("+self.trigger+ '&&'+self.filters+"&&"+self.twoLeptons + "&&"  + self.leptonPt + "&&" + self.leptonDR + "&&"  + self.leptonsMll +  ")"
      self.goodLepton3l = "("+self.trigger+ "&&"+self.filters+"&&"+self.twoLeptons + "&&"  + self.leptonPt +   ")"
      self.goodLepton4l = "("+self.trigger+ "&&"+self.filters+"&&"+self.twoLeptons + "&&"  + self.leptonPt +   ")"
      
      #self.nSRA = "(nJetSel_Edge == 2 && nPFHad10_Edge == 0 && nPFLep5_Edge <= 2 && nJetSel_Edge <= 3 && nBJetMedium25_Edge == 0 && htJet35j_Edge > 500 && mt2_Edge>80 && lepsMll_Edge >= 86 && lepsMll_Edge < 96 && (abs(j1MetDPhi_Edge)>= 0.4)&& (abs(j2MetDPhi_Edge)>= 0.4))"
      self.fromZ = '(abs(Lep1_mcMatchId_Edge) == 23 && abs(Lep2_mcMatchId_Edge) == 23)'
      self.masscut = '(lepsMll_Edge >= 86 && lepsMll_Edge < 96)'
      self.nSRAbveto = "(nJetSel_Edge >= 2 && nJetSel_Edge <= 3 && nBJetMedium25_Edge == 0 && htJet35j_Edge > 500 && mt2_Edge>80 && (abs(j1MetDPhi_Edge)>= 0.4)&& (abs(j2MetDPhi_Edge)>= 0.4))"
      self.nSRBbveto = "(nJetSel_Edge >= 4 && nJetSel_Edge <= 5 && nBJetMedium25_Edge == 0 && htJet35j_Edge > 500 && mt2_Edge>80 && (abs(j1MetDPhi_Edge)>= 0.4)&& (abs(j2MetDPhi_Edge)>= 0.4))"
      self.nSRCbveto = "(nJetSel_Edge >= 6                      && nBJetMedium25_Edge == 0 &&                        mt2_Edge>80 && (abs(j1MetDPhi_Edge)>= 0.4)&& (abs(j2MetDPhi_Edge)>= 0.4))"
      self.nSRAb = "(nJetSel_Edge >= 2 && nJetSel_Edge <= 3 && nBJetMedium25_Edge >= 1 && htJet35j_Edge > 200 && mt2_Edge>100 && (abs(j1MetDPhi_Edge)>= 0.4)&& (abs(j2MetDPhi_Edge)>= 0.4))"
      self.nSRBb = "(nJetSel_Edge >= 4 && nJetSel_Edge <= 5 && nBJetMedium25_Edge >= 1 && htJet35j_Edge > 200 && mt2_Edge>100 && (abs(j1MetDPhi_Edge)>= 0.4)&& (abs(j2MetDPhi_Edge)>= 0.4))"
      self.nSRCb = "(nJetSel_Edge >= 6                      && nBJetMedium25_Edge >= 1 &&                        mt2_Edge>100 && (abs(j1MetDPhi_Edge)>= 0.4)&& (abs(j2MetDPhi_Edge)>= 0.4))"

      self.oldThirdLepVeto = '(nLepLoose_Edge == 2 && nPFHad10_Edge == 0 && nPFLep5_Edge <= 2 )'
      self.ThirdLeptonVeto = '(nPFHad10_Edge == 0 && nPFLep5_Edge <= 2)'
      self.ThirdLeptonVetoOLD = 'nLepLoose_Edge == 2'
      self.tightIso = 'max(Lep1_miniRelIso_Edge, Lep2_miniRelIso_Edge) < 0.05'
      self.twoTightLeptons = 'nLepTight_Edge == 2'
      self.threeTightLeptons = 'nLepTight_Edge == 3'
      self.fourTightLeptons = 'nLepTight_Edge == 4'
      self.twoLooseLeptons = 'nLepLoose_Edge == 2'
      self.fourLooseLeptons = 'nLepLoose_Edge > 3'
      ########################################################################
      ######Basic Jets and bJets Cuts###############################################
      ########################################################################
      self.nj2 = "(nJetSel_Edge >= 2)"
      self.njless2 = "(nJetSel_Edge < 2)"
      self.nbjless2 = "(nBJetMedium25_Edge < 2)"
      self.nj1 = "(nJetSel_Edge >= 1)"
      self.nj0 = "(nJetSel_Edge >= 0)"
      self.nj25eq0 = "(nJet25_Edge == 0)"
      self.nj25eq0_up = "(nJet25_jecUp_Edge == 0)"
      self.nj25eq0_dn = "(nJet25_jecDn_Edge == 0)"
      self.nj35eq0 = "(nJet35_Edge == 0)"
      self.nj25eq1 = "(nJet25_Edge == 1)"
      self.njExtEtaeq0 = "(nJet25extEta_Edge == 0)"
      self.njeq0 = "(nJetSel_Edge == 0)"
      self.njeq1 = "(nJetSel_Edge == 1)"
      self.nbj2 = "(nBJetMedium25_Edge >= 2)"
      self.nbj1 = "(nBJetMedium25_Edge >= 1)"
      self.nbj0 = "(nBJetMedium25_Edge >= 0)"
      self.bveto = "(nBJetMedium25_Edge == 0)"
      self.bvetoMedium25 = "(nBJetMedium25_Edge == 0)"
      self.bvetoLoose35 = "(nBJetLoose35_Edge == 0)"
      self.bvetoLoose25 = "(nBJetLoose25_Edge == 0)"
      self.binv = "(!(nBJetMedium25_Edge == 0))"
      self.njExact2 = '(nJetSel_Edge == 2)'
      self.nbExact2 = '(nBJetMedium25_Edge == 2)'
      self.mjj110  = 'dphiMjj_Edge < 110.'
      self.mbb150  = 'mbb_Edge < 150.'
      self.mZ2g20  = 'mZ2_Edge > 20.'
      self.mZ1inZOLD  = 'mZ1_Edge >= 86 && mZ1_Edge < 96'
      self.mZ1inZ  = 'mllBestZ_Edge >= 86 && mllBestZ_Edge < 96'
      self.mZ2inZ  = 'mZ2_Edge >= 86 && mZ2_Edge < 96'
      self.ZmassInc  = 'mZ1_Edge >= 80 && mZ1_Edge < 100'
      self.lepsFromZ  = 'abs(Lep1_mcMatchId_Edge) == 23 && abs(Lep2_mcMatchId_Edge) == 23'
      self.lepsNotFromZ  = '(abs(Lep1_mcMatchId_Edge) == 24 || abs(Lep2_mcMatchId_Edge) == 24)'
      #self.lepsNotFromZ  = '(abs(Lep1_mcMatchId_Edge) == 24 && abs(Lep2_mcMatchId_Edge) == 23) || (abs(Lep1_mcMatchId_Edge) == 23 && abs(Lep2_mcMatchId_Edge) == 24)'
      
      ########################################################################
      ######Basic MET cuts##########################################################
      ########################################################################
      self.METl50  = "(MET_pt_Edge < 50)"
      self.METl100 = "(MET_pt_Edge < 100)"
      self.METl150 = "(MET_pt_Edge < 150)"
      self.METg20 = "(MET_pt_Edge > 20)"
      self.METg30 = "(MET_pt_Edge > 30)"
      self.METg60 = "(MET_pt_Edge > 60)"
      self.METg65 = "(MET_pt_Edge > 65)"
      self.METg50 = "(MET_pt_Edge > 50)"
      self.METg80 = "(MET_pt_Edge >  80)"
      self.METg70 = "(MET_pt_Edge >  70)"
      self.METg50 = "(MET_pt_Edge >= 50)"
      self.METg100 = "(MET_pt_Edge >= 100)"
      self.METg120 = "(MET_pt_Edge >= 120)"
      self.METg150 = "(MET_pt_Edge >= 150)"
      self.METg200 = "(MET_pt_Edge >= 200)"
      self.METg250 = "(MET_pt_Edge >= 250)"
      self.JETg50 = "(JetSel_Edge_pt >= 50)"
      self.JETg100 = "(JetSel_Edge_pt >= 100)"
      self.JETg150 = "(JetSel_Edge_pt >= 150)"
      self.JETg200 = "(JetSel_Edge_pt >= 200)"
      self.JETg250 = "(JetSel_Edge_pt >= 250)"
      self.JET = "(JetSel_Edge_pt >= 250)"
      
      self.MET50_100 = "(MET_pt_Edge >= 50 && MET_pt_Edge < 100)"
      self.MET100_150 = "(MET_pt_Edge >= 100 && MET_pt_Edge < 150)"
      self.MET150_250 = "(MET_pt_Edge >= 150 && MET_pt_Edge < 250)"
      self.MET250 = "(MET_pt_Edge >= 250)"
      self.lep1W = "abs(Lep1_mcMatchId_Edge) == 24 "
      self.lep1Z = "abs(Lep1_mcMatchId_Edge) == 23 "
      self.lep2W = "abs(Lep2_mcMatchId_Edge) == 24 "
      self.lep2Z = "abs(Lep2_mcMatchId_Edge) == 23 "
      self.bothTau = "Lep1_mcMatchTau_Edge == 1 && Lep2_mcMatchTau_Edge == 1 "
      self.bothNoTau = "Lep1_mcMatchTau_Edge == 0 && Lep2_mcMatchTau_Edge == 0 "
      self.dPhiJETMET = " (abs(j1MetDPhi_Edge)>= 0.4)&& (abs(j2MetDPhi_Edge)>= 0.4)"
      self.dPhiJET1MET = "abs(j1MetDPhi_Edge)>= 0.4 && j1MetDPhi_Edge> -10 && nJetSel_Edge > 0"
      
      self.mass_700_25 = "GenSusyMScan1_Edge == 700 && GenSusyMScan2_Edge == 25"
      self.mass_600_25 = "GenSusyMScan1_Edge == 600 && GenSusyMScan2_Edge == 25"
      self.mass_650_25 = "GenSusyMScan1_Edge == 650 && GenSusyMScan2_Edge == 25"
      self.mass_550_25 = "GenSusyMScan1_Edge == 550 && GenSusyMScan2_Edge == 25"

      ########################################################################
      ######Basic NLL cut#####################################################
      ########################################################################
      self.NLL = '(nll(MET_pt_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) > 21)'

      ########################################################################
      ######Basic MT2 cuts####################################################
      ########################################################################
      self.mT2_40 = "(mt2_Edge >= 40)"
      self.mT2_50 = "(mt2_Edge >= 50)"
      self.mT2_60 = "(mt2_Edge >= 60)"
      self.mT2_80 = "(mt2_Edge >= 80)"
      self.mT2_90 = "(mt2_Edge >= 90)"
      self.mT2_100 = "(mt2_Edge >= 100)"
      self.mT2b200 = '(mt2bb_Edge >= 200.)'

      self.FSCentralJetCleaning = 'FS_central_jets_Edge > 0'
      
      ########################################################################
      ######Basic Mass cuts###################################################
      ########################################################################
      self.wideZ = "lepsMll_Edge > 76 && lepsMll_Edge < 106"
      self.mZ140 =  "mZ1_Edge < 140"
      self.WmT50 =  "WmT_Edge > 50"
      self.MT2MET =  "((nJet25_Edge == 0 ) && (mt2_Edge < 90 && mt2_Edge > 40))||(nJet25_Edge > 0)"
      self.mZ4_120 = "mZ2_Edge < 120 && mZ2_Edge >= 60"
      self.Zmass = "lepsMll_Edge >= 86 && lepsMll_Edge < 96"
      self.ZmassLoose = "lepsMll_Edge >= 76 && lepsMll_Edge < 106"
      self.belowZmass = "lepsMll_Edge < 86"
      self.aboveZmass = "lepsMll_Edge > 96"
      self.mll20_60= "lepsMll_Edge < 60. && lepsMll_Edge >= 20 "
      self.mll60_86 = "lepsMll_Edge < 86. && lepsMll_Edge >= 60 "
      self.mll60_120= "lepsMll_Edge < 120. && lepsMll_Edge >= 60 "
      self.mll81_101= "lepsMll_Edge < 101. && lepsMll_Edge >= 81 "
      self.mll86_96= "lepsMll_Edge < 96. && lepsMll_Edge >= 86 "
      self.mll70_110= "lepsMll_Edge < 110. && lepsMll_Edge >= 70 "
      self.mll96_150= "lepsMll_Edge < 150. && lepsMll_Edge >= 96 "
      self.mll150_200= "lepsMll_Edge < 200. && lepsMll_Edge >= 150 "
      self.mll200_300= "lepsMll_Edge < 300. && lepsMll_Edge >= 200 "
      self.mll300_400= "lepsMll_Edge < 400. && lepsMll_Edge >= 300 "
      self.mll400= "lepsMll_Edge >= 400 "
      self.ZmassExtended = "lepsMll_Edge >= 60 && lepsMll_Edge < 120"
      self.ZmassExtendedRsfof = "lepsMll_Edge >= 70 && lepsMll_Edge < 110"
      self.Zveto = "!(lepsMll_Edge >= 86 && lepsMll_Edge < 96)"
      self.ZvetoExt = "!(lepsMll_Edge >= 76 && lepsMll_Edge < 106)"
      self.baseline = self.AddList([self.nj2, self.mT2_80, self.dPhiJETMET,self.goodLepton])
      self.baselineNoTrigger = self.AddList([self.METg100, self.nj2, self.mT2_80, self.dPhiJETMET,self.goodLepton])
      self.baselineNoMT2 = self.AddList([self.nj2,  self.dPhiJETMET,self.goodLepton])
      self.baselineNoMT2NoTrigger = self.AddList([self.nj2,  self.dPhiJETMET,self.goodLepton])

      ########################################################################
      ######Edge regions #####################################################
      ########################################################################
      self.JETMETBaseline = self.AddList([self.METg150, self.nj2, self.mT2_80])
      self.JETMETBaselineNoMT2 = self.AddList([self.METg150, self.nj2])
      self.JETMETBaselineNoMT2SF_86_96 = self.AddList([self.METg150, self.nj2, self.SF, self.mll86_96])
      self.SignalRegionNoDPhi = self.AddList([self.goodLepton, self.JETMETBaseline]) 
      self.SignalRegion = self.AddList([self.goodLepton, self.JETMETBaseline, self.dPhiJETMET]) 
      self.SignalRegionNoMT2 = self.AddList([self.goodLepton, self.JETMETBaselineNoMT2, self.dPhiJETMET]) 
      self.SignalRegionNLL = self.AddList([self.goodLepton, self.JETMETBaseline, self.dPhiJETMET, self.NLL]) 

      self.RSFOFDirectCR = self.AddList([self.METg100, self.METl150, self.njExact2, self.donot(self.mll70_110)])
      self.RSFOFDirectCROF = self.AddList([self.METg100, self.METl150, self.njExact2, self.donot(self.mll70_110), self.OF])
      self.RSFOFDirectCRNoJet = self.AddList([self.METg100, self.METl150, self.donot(self.mll70_110)])
      self.RSFOFDirectCRNoMll = self.AddList([self.METg100, self.METl150, self.njExact2])
      self.RSFOFDirectCRNoMET = self.AddList([self.njExact2, self.donot(self.mll70_110)])
      self.RSFOFDirectSignalRegion = self.AddList([self.goodLepton, self.nj2, self.METg100, self.donot(self.mll70_110)])
      self.RSFOFDirectSignalRegionNoMET = self.AddList([self.goodLepton, self.nj2, self.donot(self.mll70_110)])
      self.RSFOFDirectSignalRegionNoJet = self.AddList([self.goodLepton, self.METg100, self.donot(self.mll70_110)])
      self.RSFOFSleptonSignalRegion = self.AddList([self.goodLepton, self.METg100, self.ZvetoExt])
      self.RSFOFSleptonSignalRegionNoMET = self.AddList([self.goodLepton,  self.ZvetoExt])                                                       

      self.RSFOFCR = self.AddList([self.METg100, self.METl150, self.njExact2, self.donot(self.mll70_110)])
      self.RSFOFCRnoMET = self.AddList([self.njExact2, self.donot(self.mll70_110)])
      self.RSFOFCROF = self.AddList([self.METg100, self.METl150, self.njExact2, self.donot(self.mll70_110), self.OF])
      self.RSFOFCRNoJet = self.AddList([self.METg100, self.METl150, self.donot(self.mll70_110)])
      self.RSFOFCRNoMll = self.AddList([self.METg100, self.METl150, self.njExact2])
      self.RSFOFCRNoMET = self.AddList([self.njExact2, self.donot(self.mll70_110)])
      self.RSFOFSR = self.AddList([self.nj2, self.METg150, self.donot(self.mll70_110)])
      self.RSFOFSRNoMET = self.AddList([self.nj2, self.donot(self.mll70_110)])
      self.RSFOFSRNoJet = self.AddList([self.METg150, self.donot(self.mll70_110)])
      self.RSFOFSleptonSR =  self.AddList([self.bveto, self.METg100, self.ZvetoExt, 'htJet25j_Edge  == 0  &&  mt2_Edge > 90'])  
      self.RSFOFSleptonSRNoMT2 =  self.AddList([self.bveto, self.METg100, self.ZvetoExt, 'htJet25j_Edge  == 0'])  
      self.RSFOFSleptonSRNoNJet =  self.AddList([self.bveto, self.METg100, self.ZvetoExt, 'mt2_Edge > 90'])  


      self.DYControlRegion = self.AddList([self.METl50, self.nj2, self.ZmassExtended])
      self.DYControlRegionNoJet = self.AddList([self.METl50, self.ZmassExtended])
      self.DYControlRegionNoMllNoNJet = self.AddList([self.METl50])
      self.DYControlRegionNoMll = self.AddList([self.METl50, self.nj2])
      self.DYControlRegionNoMET = self.AddList([self.nj2, self.ZmassExtended])
      self.DYControlRegionNoMllNoMET = self.AddList([self.nj2])

      self.Edge20Mll60   = '(lepsMll_Edge < 60)'                                              
      self.Edge60Mll86   = '(lepsMll_Edge > 60 ) && (lepsMll_Edge < 86 )'
      self.Edge96Mll150  = '(lepsMll_Edge > 96 ) && (lepsMll_Edge < 150)'
      self.Edge150Mll200 = '(lepsMll_Edge > 150) && (lepsMll_Edge < 200)'
      self.Edge200Mll300 = '(lepsMll_Edge > 200) && (lepsMll_Edge < 300)'
      self.Edge300Mll400 = '(lepsMll_Edge > 300) && (lepsMll_Edge < 400)'
      self.Edge400MllInf = '(lepsMll_Edge > 400)'
                                                                                            
      self.ttBarLike    = ' nll(MET_pt_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) < 21'
      self.NonttBarLike = ' nll(MET_pt_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) > 21'
      self.loNLL = self.ttBarLike
      self.hiNLL = self.NonttBarLike                                                        
      self.mt290 = 'mt2_Edge > 90'                                                        

      ########################################################################
      ######EWK signal regions ###############################################
      self.slepBaseline  = self.AddList([self.bvetoLoose25, self.METg100, self.ZvetoExt, self.ThirdLeptonVeto,  'htJet35j_Edge  == 0  && Lep1_pt_Edge > 50'])
      self.slep  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, self.ZvetoExt, self.goodLepton, "nJet25_Edge == 0 && mt2_Edge > 90&& Lep1_pt_Edge > 50"])
      self.slepLowPt  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, self.ZvetoExt, self.goodLepton])
      self.slepRegPt  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, self.ZvetoExt, self.goodLepton])
      self.slep0jet  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, self.ZvetoExt, ' nJet25_Edge  == 0 &&  mt2_Edge > 90 && Lep1_pt_Edge > 50'])
      self.slepInclJet  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, self.ZvetoExt, 'mt2_Edge > 90 && Lep1_pt_Edge > 50'])
      self.slepExclJet  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, self.ZvetoExt, 'nJet25_Edge > 0 && mt2_Edge > 90 && Lep1_pt_Edge > 50'])
      self.slep0jetNoMT2  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, self.ZvetoExt,  'htJet25j_Edge  == 0'])
      self.slep0jetNoMET  = self.AddList([self.bvetoLoose25, self.ThirdLeptonVeto, self.ZvetoExt,  'htJet25j_Edge  == 0 && mt2_Edge > 90'])
      self.region3lSlepton  = self.AddList([self.threeTightLeptons, self.METg70, self.bvetoLoose25, self.nj25eq0, 'WmT_Edge > 50'])
      self.region4lSlepton  = self.AddList([self.goodLepton, self.fourTightLeptons,  self.bvetoLoose25, self.nj25eq0, "mZ2_Edge > 40"])
      self.region4lSleptonNoNJetCut  = self.AddList([self.goodLepton, self.fourTightLeptons, "mZ2_Edge > 40"])
      self.region4lSleptonRelaxedCut  = self.AddList([self.goodLepton, self.fourTightLeptons])
      self.region4lSleptonIncJet  = self.AddList([self.goodLepton, self.fourTightLeptons,  self.bvetoLoose25, "mZ2_Edge > 40"])
      self.region2l2nuSlepton  = self.AddList([self.goodLepton, self.twoTightLeptons, self.bvetoLoose25, self.nj25eq0])
      self.region4l      = self.AddList([self.goodLepton4l, self.fourTightLeptons,  self.mZ2g20, self.bveto, self.dPhiJETMET, "mZ2_Edge > 40"])
      self.region3l      = self.AddList([self.goodLepton3l, self.threeTightLeptons, self.METg70, self.bveto, self.dPhiJETMET,  'WmT_Edge > 55'])
      #self.region3l      = self.AddList([self.goodLepton, self.threeTightLeptons, self.dPhiJETMET,self.METg60, self.bveto, self.nj2])
      self.region4lInc   = self.AddList([self.goodLepton4l, self.fourLooseLeptons, self.ZmassInc, self.mZ140, self.mZ4_120 ])
      self.regionttZ     = self.AddList([self.goodLepton3l, self.threeTightLeptons, self.METg30, self.nbj2, self.dPhiJETMET ])
      #self.regionttZ     = self.AddList([self.goodLepton, self.threeTightLeptons, self.dPhiJETMET, self.nbj2, self.METg30, self.nj2])

      ########################################################################
      ######EWK signal regions ###############################################
      ########################################################################
      self.Baseline = self.AddList([self.nj2,  self.METg100,self.dPhiJETMET,self.goodLepton])
      self.BaselineNoTriggerNoNJet = self.AddList([self.goodLepton,  self.dPhiJET1MET])
      self.BaselineNoTrigger = self.AddList([self.METg100,self.dPhiJETMET,self.goodLepton])
      self.EdgeBaseline = self.AddList( [self.BaselineNoTrigger,self.METg150, self.mT2_80])
      self.ewinoWZnew           = self.AddList([self.BaselineNoTriggerNoNJet,  self.METg100, self.bveto,  self.ThirdLeptonVetoOLD, self.ZmassLoose, self.mT2_80])
      self.ewinoWZ              = self.AddList([self.BaselineNoTrigger,        self.bveto,  self.ThirdLeptonVetoOLD, self.Zmass, self.mT2_80, self.mjj110])
      self.ewinoWZNoMjj         = self.AddList([self.BaselineNoTrigger,        self.bveto,  self.ThirdLeptonVetoOLD, self.Zmass, self.mT2_80])
      self.ewinoWZResolved      = self.AddList([self.BaselineNoTriggerNoNJet,   self.METg100, self.bveto,  self.ThirdLeptonVetoOLD, self.ZmassLoose, self.mT2_80, "lepsZPt_Edge <= 200"])
      self.ewinoWZBoosted       = self.AddList([self.BaselineNoTriggerNoNJet,   self.METg100, self.bveto,  self.ThirdLeptonVetoOLD, self.ZmassLoose, self.mT2_80, "lepsZPt_Edge >  200"])
      self.ewinoWZResolved1jet  = self.AddList([self.BaselineNoTriggerNoNJet,   self.METg120, self.bveto, self.mT2_100,  self.ThirdLeptonVetoOLD, self.ZmassLoose, "lepsZPt_Edge <= 200 && lepsZPt_Edge >=0&& nJetSel_Edge ==1 &&  lepsDR_Edge < 1.5 "])
      self.ewinoWZResolved1jetNoDR  = self.AddList([self.BaselineNoTriggerNoNJet,   self.METg120, self.bveto, self.mT2_100,  self.ThirdLeptonVetoOLD, self.ZmassLoose, "lepsZPt_Edge <= 200 && lepsZPt_Edge >=0&& nJetSel_Edge ==1"])
      self.ewinoWZResolved1jetNoMT2  = self.AddList([self.BaselineNoTriggerNoNJet,   self.METg120, self.mT2_100, self.bveto, self.ThirdLeptonVetoOLD, self.ZmassLoose, "lepsZPt_Edge <= 200 && lepsZPt_Edge >=0&& nJetSel_Edge ==1 &&  lepsDR_Edge < 1.5 "])
      self.ewinoWZResolved1jetNoMET  = self.AddList([self.BaselineNoTriggerNoNJet,   self.METg120, self.mT2_100, self.bveto, self.ThirdLeptonVetoOLD, self.ZmassLoose, "lepsZPt_Edge <= 200 && lepsZPt_Edge >=0&& nJetSel_Edge ==1 &&  lepsDR_Edge < 1.5 "])
      self.ewinoWZResolved2jet  = self.AddList([self.BaselineNoTriggerNoNJet,   self.METg120, self.bveto, self.mT2_100,  self.ThirdLeptonVetoOLD, self.ZmassLoose, "lepsZPt_Edge <= 200 && lepsZPt_Edge >=0&& nJetSel_Edge > 1 &&  hardMjj_Edge < 100 && hardMjj_Edge >= 0"])
      self.ewinoWZResolved2jetNoMjj  = self.AddList([self.BaselineNoTriggerNoNJet,   self.METg120, self.bveto, self.mT2_100,  self.ThirdLeptonVetoOLD, self.ZmassLoose, "lepsZPt_Edge <= 200 && lepsZPt_Edge >=0&& nJetSel_Edge > 1"])
      self.ewinoWZResolved2jetNoMT2  = self.AddList([self.BaselineNoTriggerNoNJet,   self.METg120, self.mT2_100, self.bveto, self.ThirdLeptonVetoOLD, self.ZmassLoose, "lepsZPt_Edge <= 200 && lepsZPt_Edge >=0&& nJetSel_Edge > 1 &&  hardMjj_Edge < 100 && hardMjj_Edge >= 0"])
      self.ewinoWZResolved2jetNoMET  = self.AddList([self.BaselineNoTriggerNoNJet,   self.METg120, self.mT2_100, self.bveto, self.ThirdLeptonVetoOLD, self.ZmassLoose, "lepsZPt_Edge <= 200 && lepsZPt_Edge >=0&& nJetSel_Edge > 1 &&  hardMjj_Edge < 100 && hardMjj_Edge >= 0"])
      self.ewinoWZBoostedHP     = self.AddList([self.BaselineNoTriggerNoNJet,   self.METg80,self.bveto, self.mT2_100,  self.ThirdLeptonVetoOLD, self.ZmassLoose, self.dPhiJETMET , "lepsZPt_Edge > 200 && nFatJetSel_Edge > 0"])
      #self.ewinoWZBoostedHP     = self.AddList([self.BaselineNoTriggerNoNJet,   self.METg80,self.bveto, self.mT2_100,  self.ThirdLeptonVetoOLD, self.ZmassLoose, "lepsZPt_Edge > 200 && FatJetSel_Edge_tau2[0]/FatJetSel_Edge_tau1[0] < 0.15 && (abs(j2MetDPhi_Edge)>= 0.4)"])
      self.ewinoWZBoostedHPNoTau     = self.AddList([self.BaselineNoTriggerNoNJet,   self.METg80,self.bveto, self.mT2_100,  self.ThirdLeptonVetoOLD, self.ZmassLoose, "lepsZPt_Edge > 200 && (abs(j2MetDPhi_Edge)>= 0.4)"])
      self.ewinoWZBoostedHPNoMT2     = self.AddList([self.BaselineNoTriggerNoNJet,    self.mT2_100, self.METg80, self.bveto,  self.ThirdLeptonVetoOLD, self.ZmassLoose, "lepsZPt_Edge > 200 && FatJetSel_Edge_tau2[0]/FatJetSel_Edge_tau1[0] < 0.15 && (abs(j2MetDPhi_Edge)>= 0.4)"])
      self.ewinoWZBoostedHPNoMET     = self.AddList([self.BaselineNoTriggerNoNJet,    self.mT2_100, self.METg80, self.bveto,  self.ThirdLeptonVetoOLD, self.ZmassLoose, "lepsZPt_Edge > 200 && FatJetSel_Edge_tau2[0]/FatJetSel_Edge_tau1[0] < 0.15 && (abs(j2MetDPhi_Edge)>= 0.4)"])
      self.ewinoWZBoostedLP     = self.AddList([self.BaselineNoTriggerNoNJet,    self.mT2_60, self.METg100, self.bveto, self.ThirdLeptonVetoOLD, self.ZmassLoose, "lepsZPt_Edge > 200 && FatJetSel_Edge_tau2[0]/FatJetSel_Edge_tau1[0] >= 0.15 && FatJetSel_Edge_tau2[0]/FatJetSel_Edge_tau1[0] < 0.65 && (abs(j2MetDPhi_Edge)>= 0.4)"])
      self.ewinoWZBoostedLPNoTau     = self.AddList([self.BaselineNoTriggerNoNJet,    self.mT2_60, self.METg100, self.bveto, self.ThirdLeptonVetoOLD, self.ZmassLoose, "lepsZPt_Edge > 200 && (abs(j2MetDPhi_Edge)>= 0.4)"])
      self.ewinoWZBoostedLPNoMT2     = self.AddList([self.BaselineNoTriggerNoNJet,    self.mT2_60, self.METg100, self.bveto,  self.ThirdLeptonVetoOLD, self.ZmassLoose, "lepsZPt_Edge > 200 && FatJetSel_Edge_tau2[0]/FatJetSel_Edge_tau1[0] >= 0.15 && FatJetSel_Edge_tau2[0]/FatJetSel_Edge_tau1[0] < 0.65 && (abs(j2MetDPhi_Edge)>= 0.4)"])
      self.ewinoWZBoostedLPNoMET     = self.AddList([self.BaselineNoTriggerNoNJet,    self.mT2_60, self.METg100, self.bveto,  self.ThirdLeptonVetoOLD, self.ZmassLoose, "lepsZPt_Edge > 200 && FatJetSel_Edge_tau2[0]/FatJetSel_Edge_tau1[0] >= 0.15 && FatJetSel_Edge_tau2[0]/FatJetSel_Edge_tau1[0] < 0.65 && (abs(j2MetDPhi_Edge)>= 0.4)"])
      self.ewinoWZBaseline     = self.AddList([self.BaselineNoTrigger,       self.bveto,  self.ThirdLeptonVetoOLD, self.Zmass, self.dPhiJET1MET, "(signalRegion1 ==1 || signalRegion2 ==1 || signalRegion3 ==1 || signalRegion4 ==1         || signalRegion5 ==1 || signalRegion6 ==1 || signalRegion7 ==1 || signalRegion8 ==1 || signalRegion9 ==1 || signalRegion10 ==1 || signalRegion11 ==1 || signalRegion12 ==1 )&& signalRegionXX == 0 "])
      self.ewinoWZNoTrigger    = self.AddList([self.BaselineNoTrigger,       self.bveto,  self.ThirdLeptonVetoOLD, self.Zmass, self.mT2_80, self.mjj110])
      self.ewinoWZExtMll  = self.AddList([self.baseline, self.bveto, self.mjj110, self.ThirdLeptonVeto])
      self.ewinoZH        = self.AddList([self.baselineNoMT2NoTrigger,self.nbExact2, self.mT2b200, self.mbb150, self.ThirdLeptonVeto, self.Zmass])
      self.ewinoZHExtMll  = self.AddList([self.baselineNoMT2,self.nbExact2, self.mT2b200, self.mbb150, self.ThirdLeptonVeto])
      self.strongOnZBVeto     = self.AddList([self.baseline,self.ThirdLeptonVeto,self.mT2_80,  self.bveto])
      self.strongOnZWithB     = self.AddList([self.baseline,self.ThirdLeptonVeto,self.mT2_100, self.nbj1])
      self.strongOnZBase      = self.AddList([self.baseline,self.Zmass, self.ThirdLeptonVeto, self.OR(self.strongOnZBVeto,self.strongOnZWithB)])
      ########################################################################
      ######Basic HT cuts#####################################################
      ########################################################################
      self.HT = "(htJet35j_Edge > 200) "
      self.triggerHT17 = "( HLT_HT_Edge == 1  )" 
      self.triggerHT16 = "( HLT_PFHT200_Edge > 0.5 ||HLT_PFHT250_Edge > 0.5 ||HLT_PFHT300_Edge > 0.5 ||HLT_PFHT350_Edge > 0.5 ||HLT_PFHT400_Edge > 0.5 ||HLT_PFHT475_Edge > 0.5 ||HLT_PFHT600_Edge > 0.5 ||HLT_PFHT650_Edge > 0.5 ||HLT_PFHT800_Edge > 0.5  )" 
      self.triggerMET17 = " (HLT_PFMET120_PFMHT120_IDTight_PFHT60_Edge || HLT_PFMET120_PFMHT120_IDTight_Edge)"
      self.MET = "MET_pt_Edge > 130"
      self.num16 = self.AddList([self.goodLepton,self.donot(self.JETMETBaselineNoMT2SF_86_96), self.donot(self.RSFOFDirectCROF), self.filters, self.HT, self.triggerHT16, self.trigger])
      self.den16 = self.AddList([self.goodLepton,self.donot(self.JETMETBaselineNoMT2SF_86_96), self.donot(self.RSFOFDirectCROF), self.filters, self.HT, self.triggerHT16])  
      self.num17 = self.AddList([self.goodLepton,self.donot(self.JETMETBaselineNoMT2SF_86_96), self.donot(self.RSFOFDirectCROF), self.filters, self.MET, self.triggerMET17, self.trigger])
      self.den17 = self.AddList([self.goodLepton,self.donot(self.JETMETBaselineNoMT2SF_86_96), self.donot(self.RSFOFDirectCROF), self.filters, self.MET, self.triggerMET17])  

   def donot(self, cut):
     return '(!' + self.brackets(cut) + ')'

   def brackets(self, cut):
      if cut != '':
          return '('+cut+')'

   def AddList(self, cutlist):
      returncut = ''
      for cut in cutlist:
          if cut != '':
              returncut += cut
              if not cutlist.index(cut) == len(cutlist)-1:
                  returncut += ' && '
      return self.brackets(returncut)
  
   def Add(self, cut1, cut2):
      if cut1 == '':
          return cut2 
      if cut2 == '':
          return cut1
      return self.brackets(cut1 + " && " + cut2 )
  
   def OR(self, cut1, cut2):

      return self.brackets(cut1 + " || " + cut2 )

   def MaxRun(self, run):
      
      return self.brackets("run <= %d"%(run))

   def MinRun(self, run):
      
      return self.brackets("run >  %d"%(run))

      
