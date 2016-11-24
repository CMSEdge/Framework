

class CutManager:
   'This class serves as an on-demand cut server'

   def __init__(self):

      #self.twoLeptons = "nPairLep_Edge > 0 && hbheFilterIso_Edge > 0 && hbheFilterNew25ns_Edge > 0 && Flag_eeBadScFilter_Edge > 0 "
      self.twoLeptons = "nPairLep_Edge > 0 && passesFilters_Edge > 0 "
      self.tightCharge = 'Lep1_tightCharge_Edge > 0 && Lep2_tightCharge_Edge > 0'
      #self.trigMMc = "(HLT_mu17mu8_dz_Edge  > 0 || HLT_mu17mu8_Edge > 0 || HLT_mu17tkmu8_Edge > 0 || HLT_mu17tkmu8_dz_Edge > 0 || HLT_mu27tkmu8_Edge > 0 || HLT_mu30tkmu11_noniso_Edge > 0 )"
      self.trigMMc = "(HLT_mu17mu8_Edge > 0 || HLT_mu17tkmu8_Edge > 0 || HLT_mu27tkmu8_Edge > 0 || HLT_mu30tkmu11_noniso_Edge > 0 )"
      self.trigEEc = "(HLT_el17el12_dz_Edge > 0 || HLT_el23el12_dz_Edge > 0 || HLT_doubleele33_noniso_Edge > 0)"
      self.trigEMc = "(HLT_mu8el17_Edge > 0 || HLT_mu8el23_Edge > 0 || HLT_mu17el12_Edge > 0 || HLT_mu23el12_Edge > 0 || HLT_mu30ele30_noniso_Edge > 0)"
      self.leptonPt = "Lep1_pt_Edge > 25 && Lep2_pt_Edge > 20."
      self.leptonDR = "lepsDR_Edge > 0.1"       
      self.ECALCrack = "abs(abs(Lep1_eta_Edge) - 1.5) > 0.1 && abs(abs(Lep2_eta_Edge) - 1.5) > 0.1"
      self.leptonsMll = "lepsMll_Edge > 20"
      self.goodLepton = self.twoLeptons + "&&" + self.leptonPt + "&&" + self.leptonDR + "&&" + self.ECALCrack + "&&" + self.leptonsMll# + '&&' + self.tightCharge
      self.ee = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -121)"
      self.mm = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -169)"
      self.OF = "(Lep1_pdgId_Edge * Lep2_pdgId_Edge == -143)"
      self.SF = "(" + self.ee + " || " +  self.mm + ")"
      self.AF = "(" + self.SF + " || " +  self.OF + ")"
      self.nj2 = "(nJetSel_Edge >= 2)"
      self.nj1 = "(nJetSel_Edge >= 1)"
      self.nj0 = "(nJetSel_Edge >= 0)"
      self.nbj2 = "(nBJetMedium25_Edge >= 2)"
      self.nbj1 = "(nBJetMedium25_Edge >= 1)"
      self.nbj0 = "(nBJetMedium25_Edge >= 0)"
      self.nbExact2 = '(nBJetMedium25_Edge == 2)'
      self.MET100 = "(met_Edge > 100)"
      self.MET150 = "(met_Edge > 150)"
      self.MET200 = "(met_Edge > 200)"
      self.mT2_80 = "(mt2_Edge > 80)"
      self.mT2_100 = "(mt2_Edge > 100)"
      self.bVeto  = "(nBJetMedium25_Edge == 0)"
      self.ThirdLeptonVeto = '(nLepLoose_Edge == 2)'
      self.JetMETPhi04 = "abs(j1MetDPhi_Edge) >  0.4 && abs(j2MetDPhi_Edge) > 0.4"
      # to do cuts #####################################
      print 'still to do this'
      self.mjj110  = '1'
      self.mbb150  = '1'
      self.mT2b200 = '1'
      self.mbb150  = '1'
######################
      self.narrowZMass = '(lepsMll_Edge > 86 && lepsMll_Edge < 96)'
      self.trigger = "((" + self.trigMMc + " && " + self.mm + ") || (" + self.trigEEc + " && " + self.ee + ") || (" + self.trigEMc + " && " + self.OF + "))"
      self.Baseline = self.AddList([self.nj2,self.MET100,self.JetMETPhi04,self.goodLepton, self.trigger])
      self.BaselineNoTrigger = self.AddList([self.nj2,self.MET100,self.JetMETPhi04,self.goodLepton])
      self.EdgeBaseline = self.AddList( [self.MET150, self.mT2_80])
      self.ewinoCharNeu = self.AddList( [self.MET150, self.bVeto,  self.mjj110, self.narrowZMass]); 
      print 'has this a third lepton veto?'
      self.ewinoNeuNeu  = self.AddList( [self.nbExact2, self.narrowZMass, self.mT2b200, self.mbb150])
      self.strongOnZBVeto    = self.AddList( [ self.mT2_80  , self.bVeto ])
      self.strongOnZWithB    = self.AddList( [ self.mT2_100 , self.nbj1   ])
      self.strongOnZBase     = self.AddList( [self.narrowZMass, self.ThirdLeptonVeto, self.OR(strongOnZBase,strongOnZWithB)])

      self.lowmass = "lepsMll_Edge > 20 && lepsMll_Edge < 81"
      self.Zmass = "lepsMll_Edge > 81 && lepsMll_Edge < 101"
      self.loMass= "lepsMll_Edge <  81."
      self.hiMass= "lepsMll_Edge > 101."
      self.ZmassExtended = "lepsMll_Edge > 61 && lepsMll_Edge < 121"
      self.Zveto = "!(lepsMll_Edge > 81 && lepsMll_Edge < 101)"
      self.highmass = "lepsMll_Edge > 101"
      self.SignalRegionBaseLine          = self.AddList([self.goodLepton, self.JetMETBaseline, self.trigger ]) 
      self.SignalRegionBaseLineNoTrigger = self.AddList([self.goodLepton, self.JetMETBaseline]) 
#      self.SignalRegionBaseLine = self.AddList([self.goodLepton, self.trigger, self.JetMETBaseline]) 
#      self.SignalRegionBaseLineNoTrigger = self.AddList([self.goodLepton, self.JetMETBaseline]) 

      self.region3l = '(nLepTight_Edge == 3 && met_Edge > 60 && nBJetMedium25_Edge == 0)'
      self.region4l = '(nLepTight_Edge == 4)'

      self.ewinoSR = '('+self.goodLepton +'&&'+ self.SF +'&&'+ self.nj2 +'&&'+ self.Zmass +'&& nBJetMedium25_Edge == 0 && abs(j1MetDPhi_Edge) > 1. && mt2_Edge > 80. && nLepLoose_Edge == 2 )'
      
      ##### Needed by RSFOF direct calculation ########
      self.RSFOFDirectControlRegion = "((met_Edge > 100 && met_Edge < 150 && nJetSel_Edge == 2) && ((lepsMll_Edge > 20 && lepsMll_Edge < 70) || lepsMll_Edge > 110))"
      self.RSFOFDirectControlRegionNoMll = "((met_Edge > 100 && met_Edge < 150 && nJetSel_Edge == 2) && (lepsMll_Edge>20))"
      self.RSFOFDirectControlRegionNoMET = "((nJetSel_Edge == 2) && ((lepsMll_Edge > 20 && lepsMll_Edge < 70) || lepsMll_Edge > 110))"
      self.RSFOFDirectSignalRegion = "((met_Edge > 150 && nJetSel_Edge >= 2) && ((lepsMll_Edge > 20 && lepsMll_Edge < 70) || lepsMll_Edge > 110))"
      self.RSFOFDirectSignalRegionNoMET = "((nJetSel_Edge >= 2) && ((lepsMll_Edge > 20 && lepsMll_Edge < 70) || lepsMll_Edge > 110))"
      self.RSFOFDirectSignalRegionNoJet = "((nJetSel_Edge >= 2) && ((lepsMll_Edge > 20 && lepsMll_Edge < 70) || lepsMll_Edge > 110) && met_Edge > 150.)"
      self.RSFOFDirectSignalRegionNLL   = "((nJetSel_Edge >= 2) && met_Edge > 150. && ((lepsMll_Edge > 20 && lepsMll_Edge < 70) || lepsMll_Edge > 110))"
      ##### Needed by rmue calculation################ 
      self.DYControlRegion = "(met_Edge < 50 && nJetSel_Edge >= 2 && lepsMll_Edge > 81 && lepsMll_Edge < 101)"
      self.DYControlRegionNoMll = "(met_Edge < 50 && nJetSel_Edge >= 2)"
      self.DYControlRegionNoMET = "(nJetSel_Edge >= 2 && lepsMll_Edge > 81 && lepsMll_Edge < 101)"
      self.DYControlRegionNoJet = "(met_Edge < 50 && lepsMll_Edge > 81 && lepsMll_Edge < 101)"
      ##### Needed by RT calculation################ 
      self.HT = "(htJet35j_Edge > 200)"
      self.triggerHT = "(HLT_htall_Edge > 0 || HLT_htmet_Edge > 0 || HLT_atall_Edge > 0)"
      self.numerator = self.AddList([self.goodLepton, self.donot(self.Baseline), self.donot(self.RSFOFDirectControlRegion), self.HT, self.triggerHT, self.trigger])
      self.denominator = self.AddList([self.goodLepton, self.donot(self.Baseline), self.donot(self.RSFOFDirectControlRegion), self.HT, self.triggerHT])

      ## for checks of the excess
      self.tightIso = 'max(Lep1_miniRelIso_Edge, Lep2_miniRelIso_Edge) < 0.05'


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

      
