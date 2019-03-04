

class CutManager:
   'This class serves as an on-demand cut server'

   def __init__(self):

      self.protection = '((run_Edge<=276811)  ||  (278820<=run_Edge && run_Edge<=279931))'

      ########################################################################
      ######Basic Lepton Cuts ################################################
      ########################################################################
      self.twoLeptons = "nPairLep_Edge > 0 &&  Flag_HBHENoiseFilter_Edge ==1 && Flag_HBHENoiseIsoFilter_Edge ==1 && Flag_EcalDeadCellTriggerPrimitiveFilter_Edge == 1 && Flag_goodVertices_Edge == 1 && Flag_globalTightHalo2016Filter_Edge ==1 && Flag_badChargedHadronFilter_Edge == 1 && Flag_badMuonMoriond2017_Edge == 1 && Flag_badCloneMuonMoriond2017_Edge == 1"
      self.tightCharge = 'Lep1_tightCharge_Edge > 0 && Lep2_tightCharge_Edge > 0'
      self.leptonPt = "Lep1_pt_Edge > 25 && Lep2_pt_Edge > 20."
      self.sleptonLeptonPt = "Lep1_pt_Edge > 20 && Lep2_pt_Edge > 20."
      self.diLeptonPt = "lepsZPt_EdgeXXX > 25"
      self.leptonDR = "lepsDR_Edge > 0.1"       
      self.leptonDR4l = "lepsDR_Edge > 0.02"       
      self.leptonsMll = "lepsMll_Edge > 20"
      self.ee = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -121)"
      self.mm = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -169)"
      self.OF = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -143)"
      self.SF = "(" + self.ee + " || " +  self.mm + ")"
      self.AF = "(" + self.SF + " || " +  self.OF + ")"
      self.trigMMNEW = "( HLT_BIT_HLT_Mu30_TkMu11_v_Edge ==1 ||HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v_Edge == 1  || HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v_Edge ==1 || HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_Edge == 1 || HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v_Edge == 1 ||HLT_BIT_HLT_Mu27_TkMu8_v_Edge==1 )"
      self.trigEENEW = "( HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge == 1 || HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge==1 || HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v_Edge==1 || HLT_BIT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v_Edge==1)"
      self.trigEMNEW = "( HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge ==1 || HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_Edge == 1 || HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v_Edge ==1 || HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v_Edge == 1 || HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v_Edge ==1 || HLT_BIT_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v_Edge ==1 ||  HLT_BIT_HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v_Edge ==1  || HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge == 1 || HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge == 1 || HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_Edge == 1)"
      self.trigger = "((" + self.trigMMNEW + " && " + self.mm + ") || (" + self.trigEENEW + " && " + self.ee + ") || (" + self.trigEMNEW + " && " + self.OF + "))"
      self.triggerForCR = "((" + self.trigMMNEW + " && " + self.mm + ") || (" + self.trigEENEW + " && " + self.ee + ") || (" + self.trigEMNEW + " && " + self.OF + "))"
      self.goodLeptonLowPt = "("+self.trigger + "&&" +self.twoLeptons + "&&" + self.leptonDR + "&&"  + self.leptonsMll +  ")"
      self.goodLepton = "("+self.trigger + "&&" +self.twoLeptons + "&&"  + self.leptonPt + "&&" + self.leptonDR + "&&"  + self.leptonsMll +  ")"
      self.goodLepton3l = "("+self.trigger + "&&" +self.twoLeptons + "&&" + self.leptonDR  +  ")"
      self.goodLepton4l = "("+self.trigger + "&&" +self.twoLeptons + ")"
      self.goodLeptonSlepton = "("+self.trigger + "&&" +self.twoLeptons  +"&&" + self.sleptonLeptonPt + "&&" + self.leptonDR + "&&"  + self.leptonsMll +  ")"
      self.goodLeptonSignal = "( nPairLep_Edge > 0 &&"  + self.leptonPt + "&&" + self.leptonDR + "&&"  + self.leptonsMll +  ")"
      self.goodLeptonNoTrigger = "(" +self.twoLeptons  +"&&" + self.leptonPt + "&&" + self.leptonDR + "&&"  + self.leptonsMll +  ")"
      self.oldThirdLepVeto = '(nLepLoose_Edge == 2 && nPFHad10_Edge == 0 && nPFLep5_Edge <= 2 )'
      self.ThirdLeptonVeto = '(nPFHad10_Edge == 0 && nPFLep5_Edge <= 2)'
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
      self.METl50  = "(met_Edge < 50)"
      self.METl100 = "(met_Edge < 100)"
      self.METl150 = "(met_Edge < 150)"
      self.METg20 = "(met_Edge > 20)"
      self.METg30 = "(met_Edge > 30)"
      self.METg60 = "(met_Edge > 60)"
      self.METg65 = "(met_Edge > 65)"
      self.METg50 = "(met_Edge > 50)"
      self.METg80 = "(met_Edge >  80)"
      self.METg70 = "(met_Edge >  70)"
      self.METg50 = "(met_Edge >= 50)"
      self.METg100 = "(met_Edge >= 100)"
      self.METg150 = "(met_Edge >= 150)"
      self.METg200 = "(met_Edge >= 200)"
      self.METg250 = "(met_Edge >= 250)"
      self.vlowMET = "(met_Edge >= 100 && met_Edge < 150)"
      self.lowMET = "(met_Edge >= 150 && met_Edge < 225)"
      self.medMET = "(met_Edge >= 225 && met_Edge < 300)"
      self.highMET = "(met_Edge >= 300)"
      self.JETg50 = "(JetSel_Edge_pt >= 50)"
      self.JETg100 = "(JetSel_Edge_pt >= 100)"
      self.JETg150 = "(JetSel_Edge_pt >= 150)"
      self.JETg200 = "(JetSel_Edge_pt >= 200)"
      self.JETg250 = "(JetSel_Edge_pt >= 250)"
      self.JET = "(JetSel_Edge_pt >= 250)"
      
      self.MET50_100 = "(met_Edge >= 50 && met_Edge < 100)"
      self.MET100_150 = "(met_Edge >= 100 && met_Edge < 150)"
      self.MET150_250 = "(met_Edge >= 150 && met_Edge < 250)"
      self.MET250 = "(met_Edge >= 250)"
      self.lep1W = "abs(Lep1_mcMatchId_Edge) == 24 "
      self.lep1Z = "abs(Lep1_mcMatchId_Edge) == 23 "
      self.lep2W = "abs(Lep2_mcMatchId_Edge) == 24 "
      self.lep2Z = "abs(Lep2_mcMatchId_Edge) == 23 "
      self.bothTau = "Lep1_mcMatchTau_Edge == 1 && Lep2_mcMatchTau_Edge == 1 "
      self.bothNoTau = "Lep1_mcMatchTau_Edge == 0 && Lep2_mcMatchTau_Edge == 0 "
      self.dPhiJETMET = " (abs(j1MetDPhi_Edge)>= 0.4)&& (abs(j2MetDPhi_Edge)>= 0.4)"
      
      ########################################################################
      ######Basic NLL cut#####################################################
      ########################################################################
      self.NLL = '(nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) > 21)'

      ########################################################################
      ######Basic MT2 cuts####################################################
      ########################################################################
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
      self.belowZmass = "lepsMll_Edge < 86"
      self.aboveZmass = "lepsMll_Edge > 96"
      self.mll20_60= "lepsMll_Edge < 60. && lepsMll_Edge >= 20 "
      self.mll60_86= "lepsMll_Edge < 86. && lepsMll_Edge >= 60 "
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
      self.OnZ = "(lepsMll_Edge >= 76 && lepsMll_Edge < 106)"
      self.baseline = self.AddList([self.nj2, self.mT2_80, self.dPhiJETMET,self.goodLepton])
      self.baselineNoTrigger = self.AddList([self.METg100, self.nj2, self.mT2_80, self.dPhiJETMET,self.goodLeptonNoTrigger])
      self.baselineNoMT2 = self.AddList([self.nj2,  self.dPhiJETMET,self.goodLepton])
      self.baselineNoMT2NoTrigger = self.AddList([self.nj2,  self.dPhiJETMET,self.goodLeptonNoTrigger])

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

      self.RSFOFDirectControlRegion = self.AddList([self.METg100, self.METl150, self.njExact2, self.donot(self.mll70_110)])
      self.RSFOFDirectControlRegionOF = self.AddList([self.METg100, self.METl150, self.njExact2, self.donot(self.mll70_110), self.OF])
      self.RSFOFDirectControlRegionNoJet = self.AddList([self.goodLepton, self.METg100, self.METl150, self.donot(self.mll70_110)])
      self.RSFOFDirectControlRegionNoMll = self.AddList([self.goodLepton,  self.METg100, self.METl150, self.njExact2])
      self.RSFOFDirectControlRegionNoMET = self.AddList([self.goodLepton,  self.njExact2, self.donot(self.mll70_110)])
      self.RSFOFDirectSignalRegion = self.AddList([self.goodLepton, self.nj2, self.METg150, self.donot(self.mll70_110)])
      self.RSFOFDirectSignalRegionNoMET = self.AddList([self.goodLepton, self.nj2, self.donot(self.mll70_110)])
      self.RSFOFDirectSignalRegionNoJet = self.AddList([self.goodLepton, self.METg150, self.donot(self.mll70_110)])
      self.RSFOFSleptonSignalRegion = self.AddList([self.goodLepton, self.METg150, self.ZvetoExt])
      self.RSFOFSleptonSignalRegionNoMET = self.AddList([self.goodLepton,  self.ZvetoExt])                                                       

      self.RSFOFCR = self.AddList([self.METg100, self.METl150, self.njExact2, self.donot(self.mll70_110)])
      self.RSFOFCRnoMET = self.AddList([self.njExact2, self.donot(self.mll70_110)])
      self.RSFOFCROF = self.AddList([self.METg100, self.METl150, self.njExact2, self.donot(self.mll70_110), self.OF])
      self.RSFOFCRNoJet = self.AddList([self.goodLepton, self.METg100, self.METl150, self.donot(self.mll70_110)])
      self.RSFOFCRNoMll = self.AddList([self.goodLepton,  self.METg100, self.METl150, self.njExact2])
      self.RSFOFCRNoMET = self.AddList([self.goodLepton,  self.njExact2, self.donot(self.mll70_110)])
      self.RSFOFSR = self.AddList([self.goodLepton, self.nj2, self.METg150, self.donot(self.mll70_110)])
      self.RSFOFSRNoMET = self.AddList([self.goodLepton, self.nj2, self.donot(self.mll70_110)])
      self.RSFOFSRNoJet = self.AddList([self.goodLepton, self.METg150, self.donot(self.mll70_110)])
      self.RSFOFSleptonSR =  self.AddList([self.goodLepton, self.bveto, self.METg100, self.ZvetoExt, 'htJet25j_Edge  == 0  &&  mt2_Edge > 90'])  
      self.RSFOFSleptonSRNoMT2 =  self.AddList([self.goodLepton, self.bveto, self.METg100, self.ZvetoExt, 'htJet25j_Edge  == 0'])  
      self.RSFOFSleptonSRNoNJet =  self.AddList([self.goodLepton, self.bveto, self.METg100, self.ZvetoExt, 'mt2_Edge > 90'])  


      self.DYControlRegion = self.AddList([self.goodLepton, self.METl50, self.nj2, self.ZmassExtended])
      self.DYControlRegionNoJet = self.AddList([self.goodLepton, self.METl50, self.ZmassExtended])
      self.DYControlRegionNoMllNoNJet = self.AddList([self.goodLepton, self.METl50])
      self.DYControlRegionNoMll = self.AddList([self.goodLepton, self.METl50, self.nj2])
      self.DYControlRegionNoMET = self.AddList([self.goodLepton, self.nj2, self.ZmassExtended])
      self.DYControlRegionNoMllNoMET = self.AddList([self.goodLepton, self.nj2])

      self.Edge20Mll60   = '(lepsMll_Edge < 60)'                                              
      self.Edge60Mll86   = '(lepsMll_Edge > 60 ) && (lepsMll_Edge < 86 )'
      self.Edge96Mll150  = '(lepsMll_Edge > 96 ) && (lepsMll_Edge < 150)'
      self.Edge150Mll200 = '(lepsMll_Edge > 150) && (lepsMll_Edge < 200)'
      self.Edge200Mll300 = '(lepsMll_Edge > 200) && (lepsMll_Edge < 300)'
      self.Edge300Mll400 = '(lepsMll_Edge > 300) && (lepsMll_Edge < 400)'
      self.Edge400MllInf = '(lepsMll_Edge > 400)'
                                                                                            
      self.ttBarLike    = ' nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) < 21'
      self.NonttBarLike = ' nll(met_Edge, lepsZPt_Edge, sum_mlb_Edge, lepsDPhi_Edge) > 21'
      self.loNLL = self.ttBarLike
      self.hiNLL = self.NonttBarLike                                                        
      self.mt290 = 'mt2_Edge > 90'                                                        

      ########################################################################
      ######EWK signal regions ###############################################
      self.slepBaseline  = self.AddList([self.bvetoLoose25, self.METg100, self.ZvetoExt, self.ThirdLeptonVeto,  'htJet35j_Edge  == 0  && Lep1_pt_Edge > 50'])
      self.incl  = self.AddList([self.METg50])
      self.slep  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, self.ZvetoExt, self.goodLeptonSignal, "nJet25_Edge == 0 && mt2_Edge > 90&& Lep1_pt_Edge > 50"])
      self.slepNoMT2  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, self.ZvetoExt, self.goodLeptonSignal, "nJet25_Edge == 0 &&  Lep1_pt_Edge > 50"])
      self.slepNoMET  = self.AddList([self.bvetoLoose25, self.ThirdLeptonVeto, self.ZvetoExt, self.goodLeptonSignal, "nJet25_Edge == 0 && mt2_Edge > 90&& Lep1_pt_Edge > 50"])
      self.slepLowPt  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, self.ZvetoExt, self.goodLeptonSignal])
      self.slepRegPt  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, self.ZvetoExt, self.goodLeptonSignal])
      self.slep0jetForMT2  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, self.ZvetoExt, ' nJet25_Edge  == 0 && Lep1_pt_Edge > 50'])
      self.slep0jet        = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, self.ZvetoExt, ' nJet25_Edge  == 0 &&  mt2_Edge > 90 && Lep1_pt_Edge > 50'])
      self.slep0jetOnZForMT2 = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto,self.OnZ,   ' nJet25_Edge  == 0 && Lep1_pt_Edge > 50'])
      self.slep0jetOnZ  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto,self.OnZ,   ' nJet25_Edge  == 0 &&  mt2_Edge > 90 && Lep1_pt_Edge > 50'])
      self.slep0jetNoMll  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, ' nJet25_Edge  == 0 &&  mt2_Edge > 90 && Lep1_pt_Edge > 50'])
      self.slep0jetNoMT2  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, self.ZvetoExt, ' nJet25_Edge  == 0 && Lep1_pt_Edge > 50'])
      self.slep0jetNoMET  = self.AddList([self.bvetoLoose25, self.ThirdLeptonVeto, self.ZvetoExt, ' nJet25_Edge  == 0 &&  mt2_Edge > 90 && Lep1_pt_Edge > 50'])
      self.slep0jetNoMll  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto,  ' nJet25_Edge  == 0 &&  mt2_Edge > 90 && Lep1_pt_Edge > 50'])
      self.slep0jetNoMllNoLepVeto  = self.AddList([self.bvetoLoose25, self.METg100, ' nJet25_Edge  == 0 &&  mt2_Edge > 90 && Lep1_pt_Edge > 50'])
      self.slep0jetOnZNoLepVeto  = self.AddList([self.bvetoLoose25, self.METg100, ' nJet25_Edge  == 0 &&  mt2_Edge > 90 && Lep1_pt_Edge > 50'])
      self.slep0jetNoLepVeto  = self.AddList([self.bvetoLoose25, self.METg100, self.ZvetoExt, ' nJet25_Edge  == 0 &&  mt2_Edge > 90 && Lep1_pt_Edge > 50'])
      self.slepInclJet  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, self.ZvetoExt, 'mt2_Edge > 90 && Lep1_pt_Edge > 50'])
      self.slepInclJetForMT2  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, self.ZvetoExt, 'Lep1_pt_Edge > 50'])
      self.slepExclJet  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, self.ZvetoExt, 'nJet25_Edge > 0 && mt2_Edge > 90 && Lep1_pt_Edge > 50'])
      self.slepExclJetForMT2  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, self.ZvetoExt, 'nJet25_Edge > 0 && Lep1_pt_Edge > 50'])
      #self.slep0jetNoMT2  = self.AddList([self.bvetoLoose25, self.METg100, self.ThirdLeptonVeto, self.ZvetoExt,  'htJet25j_Edge  == 0'])
      #self.slep0jetNoMET  = self.AddList([self.bvetoLoose25, self.ThirdLeptonVeto, self.ZvetoExt,  'htJet25j_Edge  == 0 && mt2_Edge > 90'])
      self.region3lSlepton  = self.AddList([self.threeTightLeptons, self.METg70, self.bvetoLoose25, self.nj25eq0, 'WmT_Edge > 50'])
      self.region4lSlepton  = self.AddList([self.goodLepton, self.fourTightLeptons,  self.bvetoLoose25, self.nj25eq0, "mZ2_Edge > 40"])
      self.region4lSleptonNoNJetCut  = self.AddList([self.goodLepton, self.fourTightLeptons, "mZ2_Edge > 40"])
      self.region4lSleptonRelaxedCut  = self.AddList([self.goodLepton, self.fourTightLeptons])
      self.region4lSleptonIncJet  = self.AddList([self.goodLepton, self.fourTightLeptons,  self.bvetoLoose25, "mZ2_Edge > 40"])
      self.region2l2nuSlepton  = self.AddList([self.goodLepton, self.twoTightLeptons, self.bvetoLoose25, self.nj25eq0])
      self.region4l  = self.AddList([self.goodLepton4l, self.fourTightLeptons,  self.dPhiJETMET, self.mZ2g20, self.bveto, self.nj2])
      self.region3l  = self.AddList([self.goodLepton, self.threeTightLeptons, self.dPhiJETMET,self.METg60, self.bveto, self.nj2])
      self.region4lInc  = self.AddList([self.goodLepton, self.fourLooseLeptons, self.ZmassInc, self.mZ140, self.mZ4_120 ])
      self.regionttZ = self.AddList([self.goodLepton, self.threeTightLeptons, self.dPhiJETMET, self.nbj2, self.METg30, self.nj2])
      self.regionEdge3l  = self.AddList([self.goodLepton, self.threeTightLeptons, self.njless2, self.dPhiJETMET, self.METg60, self.bveto])  
      self.regionEdge4l  = self.AddList([self.goodLepton, self.fourTightLeptons,  self.njless2, self.dPhiJETMET, self.mZ2g20, self.bveto])
      self.regionEdgettZ = self.AddList([self.goodLepton, self.threeTightLeptons, self.nbjless2, self.njless2, self.dPhiJETMET, self.METg30])

      ########################################################################
      ######EWK signal regions ###############################################
      ########################################################################
      self.Baseline = self.AddList([self.nj2, self.METg100,self.dPhiJETMET,self.goodLepton])
      self.BaselineNoTrigger = self.AddList([self.nj2, self.METg100,self.dPhiJETMET,self.goodLeptonNoTrigger])
      self.EdgeBaseline = self.AddList( [self.baselineNoTrigger,self.METg150, self.mT2_80])
      self.ewinoWZ      = self.AddList([self.baselineNoTrigger, self.bveto, self.mjj110, self.ThirdLeptonVeto, self.Zmass])
      self.ewinoWZNoTrigger = self.AddList([self.baselineNoTrigger, self.bveto, self.mjj110, self.ThirdLeptonVeto, self.Zmass])
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
      self.triggerHT = "( HLT_BIT_HLT_PFHT200_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT250_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT300_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT350_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT400_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT475_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT600_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT650_v_Edge > 0.5 ||HLT_BIT_HLT_PFHT800_v_Edge > 0.5  )" 
      self.numerator =   self.AddList([self.goodLeptonNoTrigger,  self.donot(self.JETMETBaselineNoMT2SF_86_96), self.donot(self.RSFOFDirectControlRegionOF), self.HT, self.triggerHT, self.trigger])
      self.denominator = self.AddList([self.goodLeptonNoTrigger,  self.donot(self.JETMETBaselineNoMT2SF_86_96), self.donot(self.RSFOFDirectControlRegionOF), self.HT, self.triggerHT])  

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

      
