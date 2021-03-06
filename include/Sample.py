import ROOT as r
from array import array
from ROOT import TTree, TFile, TCut, TH1F, TH2F, TH3F, THStack, TCanvas, TChain, SetOwnership
import include.LeptonSF
import include.FastSimSF
import include.CutManager as CutManager
import copy
import time
import math
from multiprocessing import Pool


chunkSize = 100000# chunks of 100k events (so we have the cluster active :) )

def runTH1Tasks (task):
   sample, blockName, params  = task
   return (blockName, sample.getTH1F(*params))

class Sample:
   'Common base class for all Samples'

   def __init__(self, name, friendlocation, xsection, isdata, doKfactor, isScan, isOnEOS):
      self.name = name
      self.location = friendlocation
      self.xSection = xsection
      self.doKfactor = doKfactor
      self.isData = isdata
      self.ttree = TChain('Events')
      self.fttree = TChain('sf/t')
      if isOnEOS:
         self.ttree.Add(friendlocation)
      else:
         for part in name.split('+'):
            self.ttree.Add(friendlocation+'/%s/nanoAODskim/Events.root'%part)
            self.fttree.Add(friendlocation+'/1_triggers/evVarFriend_%s.root'%part)
      self.isScan = isScan
      if not self.isData and not self.isScan:
        gw = 0.
        for i in self.ttree:
            gw = abs(i.genWeight_Edge)
            if gw: break
        self.count = 0
        for part in name.split('+'):
           ftfile = TFile(friendlocation+'/%s/nanoAODskim/Events.root'%part)
           self.count += ftfile.Get('SumWeights').GetBinContent(1)/abs(gw)
      #else:
          #self.count = self.ftfile.Get('Count').GetBinContent(1)
          #self.count = self.ftfile.Get('sf/t').GetEntries()
      self.puWeight   = '1.0'
      self.SFWeight   = '1.0'
      self.btagWeight = '1.0'
      self.triggWeight = '1.0'
      self.ISRWeight  = '1.0'

      if not self.isData and not self.isScan  > 0:
        cuts = CutManager.CutManager()
        self.lumWeight = '%s /'%self.xSection + str(self.count)
        self.puWeight    = "PileupW_Edge"
        self.btagWeight  = "weight_btagsf_Edge"
        self.SFWeight = "LepSF(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge,year)*LepSF(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge,year)"
        self.triggWeight = 'TriggerSF(Lep1_pdgId_Edge,Lep2_pdgId_Edge,year)'

        #((%s)*(%4.5f)+(%s)*(%4.5f)+(%s)*(%4.5f))'%(cuts.mm,0.878/0.947,
        #                                                               cuts.ee,0.900/0.953,
        #                                                               cuts.OF,0.844/0.905)

      if isinstance(self.isScan,tuple):                       
         self.puWeight    = "PileupW_Edge"                     
         self.btagWeight  = "weight_btagsf_Edge"               
         self.SFWeight = "LepSF(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSF(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge)*LepSFFastSim(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSFFastSim(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge)"                         
         self.ISRWeight = '1'# ISRweight_Edge'                     
         smsCount = self.ftfile.Get('CountSMS')                
         xbin = smsCount.GetXaxis().FindBin(self.isScan[0])    
         ybin = smsCount.GetYaxis().FindBin(self.isScan[1])    
         zbin = smsCount.GetZaxis().FindBin(1)                 
         bin  = smsCount.GetBin(xbin,ybin,zbin)                
         counts = smsCount.GetBinContent(bin)                  
         if self.isScan[2] == 'TChiWZ':                        
            fxsec = open('datacards/charneuXsec.txt','r')      
            xsec  = eval(fxsec.read())[self.isScan[0]][0]      
            br    = 0.102   
            fxsec.close()   
         else:              
            raise RuntimeError('No model %s avialable'%self.isScan[2])                            
         self.lumWeight = xsec*br / counts                     
                            
      elif not isinstance(self.isScan,tuple)  and self.isScan > 0:                                  
        self.lumWeight  =  1.0                                 
        self.xSection = self.isScan                            
        self.puWeight    = "PileupW_Edge"                      
        self.btagWeight  = "weight_btagsf_Edge"                
        self.SFWeight = "LepSF(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSF(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge)*LepSFFastSim(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSFFastSim(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge)"                          
        self.ISRWeight = 'ISRweight_Edge'                      
        if self.isScan == True:                                
            self.smsCount =   self.ftfile.Get('CountSMS')      
            print "using this self.smsCount ", self.smsCount, " for ", self.name                  
        else:               
            self.smsCount =  self.ftfile.Get('Events').GetEntries()                                 
            self.lumWeight = self.xSection / self.smsCount     
            print "using this self.smsCount ", self.smsCount, " for ", self.name                  


      elif self.isScan > 0:
        self.lumWeight  =  1.0
        self.xSection = self.isScan
        self.puWeight    = "PileupW_Edge"
        self.btagWeight  = "weight_btagsf_Edge"
        self.SFWeight = "LepSF(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge,year)*LepSF(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge,year)*LepSFFastSim(Lep1_pt_Edge,Lep1_eta_Edge,Lep1_pdgId_Edge)*LepSFFastSim(Lep2_pt_Edge,Lep2_eta_Edge,Lep2_pdgId_Edge)"
        self.ISRWeight = 'ISRweight_Edge'
        if self.isScan == True:
            self.smsCount =   self.ftfile.Get('CountSMS')
        else:
            #self.smsCount =  self.ftfile.Get('sf/t').GetEntries()
            #self.lumWeight = self.xSection / self.smsCount
            self.lumWeight = self.isScan

   def printSample(self):
      print "#################################"
      print "Sample Name: ", self.name
      print "Sample Location: ", self.location
      print "Sample XSection: ", self.xSection
      print "Sample IsData: ", self.isData
      print "Sample LumWeight: ", self.lumWeight
      print "#################################"


   def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel,  extraWeight, doKfactorGENVar, chunk=-1, forceDataWeight=False):
      self.ttree.AddFriend(self.fttree)
   #def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel, ofBin=True, extraWeight="1", ylabel = 'Events', doKfactorGENVar = "noKFactor"):
#      if doKfactorGENVar == 'doTGraph':
#          fileD =  TFile("kfactors.root");
#          hScale = fileD.Get("kfactorsPt");
#          hScale.Fit("pol6")
#          fit = hScale.GetFunction("pol6");
#          fileD.Close()    
      ofBin = True
      ylabel = 'Events'
      if(xmin == xmax):                                          
          _nbins = len(nbin)-1
          _arr = array('d', nbin)
          h = TH1F(name+'_noOF_%d'%chunk, "", _nbins, _arr)
          _newarr = _arr + array('d', [ 2*_arr[-1]-_arr[-2] ])
          h_of = TH1F(name, "", _nbins+1, _newarr)                                     
      else:
           h = TH1F(name+'_noOF_%d'%chunk, "", nbin, xmin, xmax)
           bw = int((xmax-xmin)/nbin)
           #if ylabel == "au":
           #    ylabel = "A.U"
           #else:
           ylabel = ylabel +"/ " + str(bw) + " GeV"
           h_of = TH1F(name + '%d'%chunk, '', nbin+1, xmin, xmax+bw)
      h.Sumw2()
      h.GetXaxis().SetTitle(xlabel)
      h.GetYaxis().SetTitle(ylabel)

      h_of.Sumw2()
      h_of.GetXaxis().SetTitle(xlabel)
      h_of.GetYaxis().SetTitle(ylabel)

      addCut = ""
      
      kf = "1"  
      if(self.isData == 0):
          if (self.doKfactor == 1 ): #this is the kfactor for ZZto4l
              if doKfactorGENVar == 'ZZmass':
                  kf = "(1.23613*(abs(GENmassZZ_Edge)>0.0&&abs(GENmassZZ_Edge)<=25.0)+1.1755*(abs(GENmassZZ_Edge)>25.0 &&abs(GENmassZZ_Edge)<=50.0 )+ 1.1704*(abs(GENmassZZ_Edge)>50.0 &&abs(GENmassZZ_Edge)<=75.0 )+ 1.0314*(abs(GENmassZZ_Edge)>75.0 &&abs(GENmassZZ_Edge)<=100.0)+1.0528*(abs(GENmassZZ_Edge)>100.0&&abs(GENmassZZ_Edge)<=125.0)+1.1128*(abs(GENmassZZ_Edge)>125.0&&abs(GENmassZZ_Edge)<=150.0)+1.1336*(abs(GENmassZZ_Edge)>150.0&&abs(GENmassZZ_Edge)<=175.0)+1.1035*(abs(GENmassZZ_Edge)>175.0&&abs(GENmassZZ_Edge)<=200.0)+1.1005*(abs(GENmassZZ_Edge)>200.0&&abs(GENmassZZ_Edge)<=225.0)+1.1097*(abs(GENmassZZ_Edge)>225.0&&abs(GENmassZZ_Edge)<=250.0)+1.1206*(abs(GENmassZZ_Edge)>250.0&&abs(GENmassZZ_Edge)<=275.0)+1.1158*(abs(GENmassZZ_Edge)>275.0&&abs(GENmassZZ_Edge)<=300.0)+1.1390*(abs(GENmassZZ_Edge)>300.0&&abs(GENmassZZ_Edge)<=325.0)+1.1485*(abs(GENmassZZ_Edge)>325.0&&abs(GENmassZZ_Edge)<=350.0)+1.1461*(abs(GENmassZZ_Edge)>350.0&&abs(GENmassZZ_Edge)<=375.0)+1.1457*(abs(GENmassZZ_Edge)>375.0&&abs(GENmassZZ_Edge)<=400.0)+1.1382*(abs(GENmassZZ_Edge)>400.0&&abs(GENmassZZ_Edge)<=425.0)+1.1552*(abs(GENmassZZ_Edge)>425.0&&abs(GENmassZZ_Edge)<=450.0)+1.1367*(abs(GENmassZZ_Edge)>450.0&&abs(GENmassZZ_Edge)<=475.0)+1.1322*(abs(GENmassZZ_Edge)>475.0))"
              if doKfactorGENVar == 'ZZpt':
                  kf =  "(0.64155*(abs(GENptZZ_Edge)>0.0&&abs(GENptZZ_Edge)<=5.0)+1.0998*(abs(GENptZZ_Edge)>5.0 &&abs(GENptZZ_Edge)<=10.0) +1.2939*(abs(GENptZZ_Edge)>10.0&&abs(GENptZZ_Edge)<=15.0)+1.3785*(abs(GENptZZ_Edge)>15.0&&abs(GENptZZ_Edge)<=20.0)+1.4243*(abs(GENptZZ_Edge)>20.0&&abs(GENptZZ_Edge)<=25.0)+1.4503*(abs(GENptZZ_Edge)>25.0&&abs(GENptZZ_Edge)<=30.0)+1.4701*(abs(GENptZZ_Edge)>30.0&&abs(GENptZZ_Edge)<=35.0)+1.4882*(abs(GENptZZ_Edge)>35.0&&abs(GENptZZ_Edge)<=40.0)+1.5057*(abs(GENptZZ_Edge)>40.0&&abs(GENptZZ_Edge)<=45.0)+1.5021*(abs(GENptZZ_Edge)>45.0&&abs(GENptZZ_Edge)<=50.0)+1.5091*(abs(GENptZZ_Edge)>50.0&&abs(GENptZZ_Edge)<=55.0)+1.5246*(abs(GENptZZ_Edge)>55.0&&abs(GENptZZ_Edge)<=60.0)+1.5240*(abs(GENptZZ_Edge)>60.0&&abs(GENptZZ_Edge)<=65.0)+1.5241*(abs(GENptZZ_Edge)>65.0&&abs(GENptZZ_Edge)<=70.0)+1.5542*(abs(GENptZZ_Edge)>70.0&&abs(GENptZZ_Edge)<=75.0)+1.5254*(abs(GENptZZ_Edge)>75.0&&abs(GENptZZ_Edge)<=80.0)+1.5789*(abs(GENptZZ_Edge)>80.0&&abs(GENptZZ_Edge)<=85.0)+1.5303*(abs(GENptZZ_Edge)>85.0&&abs(GENptZZ_Edge)<=90.0)+1.5614*(abs(GENptZZ_Edge)>90.0&&abs(GENptZZ_Edge)<=95.0)+1.5446*(abs(GENptZZ_Edge)>95.0&&abs(GENptZZ_Edge)<=100.0)+1.5722*(abs(GENptZZ_Edge)>100.0))"
              if doKfactorGENVar == 'noKFactor':
                  kf = "1"
          if (self.doKfactor == 2): #this is the kfactor for ZZto2l2nu
              if doKfactorGENVar == 'ZZmass':
                  kf = "(1.25094*(abs(GENmassZZ_Edge)>0.0&&abs(GENmassZZ_Edge)<=25.0)+1.2245*(abs(GENmassZZ_Edge)>25.0 &&abs(GENmassZZ_Edge)<=50.0 )+1.1928*(abs(GENmassZZ_Edge)>50.0 &&abs(GENmassZZ_Edge)<=75.0 )+1.0459*(abs(GENmassZZ_Edge)>75.0 &&abs(GENmassZZ_Edge)<=100.0)+1.0832*(abs(GENmassZZ_Edge)>100.0&&abs(GENmassZZ_Edge)<=125.0)+1.0999*(abs(GENmassZZ_Edge)>125.0&&abs(GENmassZZ_Edge)<=150.0)+1.1669*(abs(GENmassZZ_Edge)>150.0&&abs(GENmassZZ_Edge)<=175.0)+1.1039*(abs(GENmassZZ_Edge)>175.0&&abs(GENmassZZ_Edge)<=200.0)+1.1059*(abs(GENmassZZ_Edge)>200.0&&abs(GENmassZZ_Edge)<=225.0)+1.1069*(abs(GENmassZZ_Edge)>225.0&&abs(GENmassZZ_Edge)<=250.0)+1.1119*(abs(GENmassZZ_Edge)>250.0&&abs(GENmassZZ_Edge)<=275.0)+1.1352*(abs(GENmassZZ_Edge)>275.0&&abs(GENmassZZ_Edge)<=300.0)+1.1189*(abs(GENmassZZ_Edge)>300.0&&abs(GENmassZZ_Edge)<=325.0)+1.1389*(abs(GENmassZZ_Edge)>325.0&&abs(GENmassZZ_Edge)<=350.0)+1.1546*(abs(GENmassZZ_Edge)>350.0&&abs(GENmassZZ_Edge)<=375.0)+1.1734*(abs(GENmassZZ_Edge)>375.0&&abs(GENmassZZ_Edge)<=400.0)+1.2009*(abs(GENmassZZ_Edge)>400.0&&abs(GENmassZZ_Edge)<=425.0)+1.1891*(abs(GENmassZZ_Edge)>425.0&&abs(GENmassZZ_Edge)<=450.0)+1.1854*(abs(GENmassZZ_Edge)>450.0&&abs(GENmassZZ_Edge)<=475.0)+1.12864*(abs(GENmassZZ_Edge)>475.0))"
              if doKfactorGENVar == 'ZZpt':
                  kf  = "(0.7436*(abs(GENptZZ_Edge)>0.0&&abs(GENptZZ_Edge)<=5.0)+1.14789*(abs(GENptZZ_Edge)>5.0&&abs(GENptZZ_Edge)<=10.0)+1.33815*(abs(GENptZZ_Edge)>10.0&&abs(GENptZZ_Edge)<=15.0)+1.41420*(abs(GENptZZ_Edge)>15.0&&abs(GENptZZ_Edge)<=20.0)+1.45511*(abs(GENptZZ_Edge)>20.0&&abs(GENptZZ_Edge)<=25.0)+1.47569*(abs(GENptZZ_Edge)>25.0&&abs(GENptZZ_Edge)<=30.0)+1.49053*(abs(GENptZZ_Edge)>30.0&&abs(GENptZZ_Edge)<=35.0)+1.50622*(abs(GENptZZ_Edge)>35.0&&abs(GENptZZ_Edge)<=40.0)+1.50328*(abs(GENptZZ_Edge)>40.0&&abs(GENptZZ_Edge)<=45.0)+1.52186*(abs(GENptZZ_Edge)>45.0&&abs(GENptZZ_Edge)<=50.0)+1.52043*(abs(GENptZZ_Edge)>50.0&&abs(GENptZZ_Edge)<=55.0)+1.53977*(abs(GENptZZ_Edge)>55.0&&abs(GENptZZ_Edge)<=60.0)+1.53491*(abs(GENptZZ_Edge)>60.0&&abs(GENptZZ_Edge)<=65.0)+1.51772*(abs(GENptZZ_Edge)>65.0&&abs(GENptZZ_Edge)<=70.0)+1.54494*(abs(GENptZZ_Edge)>70.0&&abs(GENptZZ_Edge)<=75.0)+1.57762*(abs(GENptZZ_Edge)>75.0&&abs(GENptZZ_Edge)<=80.0)+1.55078*(abs(GENptZZ_Edge)>80.0&&abs(GENptZZ_Edge)<=85.0)+1.57078*(abs(GENptZZ_Edge)>85.0&&abs(GENptZZ_Edge)<=90.0)+1.56162*(abs(GENptZZ_Edge)>90.0&&abs(GENptZZ_Edge)<=95.0)+1.54183*(abs(GENptZZ_Edge)>95.0&&abs(GENptZZ_Edge)<=100.0)+1.58485*(abs(GENptZZ_Edge)>100.0))"
              if doKfactorGENVar == 'noKFactor':
                  kf = "1"
          
          if (self.doKfactor == 0): #all other MC has no NNLO/NLO reweighting
              kf = "1"  
          cut = cut + "*1000*( %s *"%self.lumWeight + str(lumi) + " * genWeight_Edge/abs(genWeight_Edge) * " + self.puWeight +" * "+self.SFWeight+" * "+self.btagWeight+" * "+extraWeight +" * "+kf + " * " + self.triggWeight + " )" 
      else: 
         addDataFilters = "&&(  Flag_eeBadScFilter_Edge == 1)"
         cut = "("+ cut + addDataFilters+ ")" + "* (" + (extraWeight if forceDataWeight else '1') +")"
      if (self.doKfactor == 1): print "doing ", doKfactorGENVar, "for ZZ4l kfactor!"
      if (self.doKfactor == 2): print "doing ", doKfactorGENVar, "for  ZZ2l kfactor!"

      if (chunk > -1):
         allEntries = self.ttree.GetEntries()
         firstEntry = chunkSize * chunk
         maxEntries = chunkSize # the last chunk may go beyond the end of the tree, but ROOT stops anyway so we don't care
      else:
         firstEntry = 0
         maxEntries = allEntries+1 # +1 just in case 
      self.ttree.Project(h.GetName(), var, cut, options, maxEntries, firstEntry)

      for _bin in range(1, h.GetNbinsX()+2):
          h_of.SetBinContent(_bin, h.GetBinContent(_bin))
          h_of.SetBinError  (_bin, h.GetBinError  (_bin))
      return (h_of if ofBin else h)

    
   def getTH2F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel, extraWeight):
   
     if(xmin == xmax) and (ymax == ymin):
        h = TH2F(name, "", len(nbinx)-1, array('d', nbinx),len(nbiny)-1, array('d', nbiny))
     elif (xmin == xmax):
        h = TH2F(name, "", len(nbinx)-1, array('d', nbinx),nbiny,ymin,ymax)
     elif (ymin == ymax):
        h = TH2F(name, "", nbinx,xmin,xmax,len(nbiny)-1, array('d', nbiny))
     else: 
        h = TH2F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax)
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     
     if(self.isData == 0):
        cut = cut + "* ( " + str(self.lumWeight*lumi) + " * genWeight_Edge/abs(genWeight_Edge) * " + self.puWeight + " * " + self.SFWeight + " * " + self.btagWeight + " * " +  self.triggWeight  + "*" + extraWeight + " )" 
     else: 
        cut = cut + "* ( " + extraWeight + ")"
     self.ttree.Project(name, var, cut, options) 
     return h

   def getTH3F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax, cut, options, xlabel, ylabel, zlabel, extraWeight):
   
     if(xmin == xmax) and (ymax == ymin) and (zmax == zmin):
        h = TH3F(name, "", len(nbinx)-1, array('d', nbinx), len(nbiny)-1, array('d', nbiny), len(nbinz)-1, array('d', nbinz))
     else: 
        h = TH3F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax)
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     h.GetZaxis().SetTitle(zlabel)
     
     if(self.isData == 0):
        cut = cut + "* ( " + str(self.lumWeight*lumi) + " * genWeight_Edge/abs(genWeight_Edge) * " + self.puWeight + " * " + self.SFWeight + " * " + self.btagWeight + " * " +  self.triggWeight  + "*" + extraWeight + " )" 
     else: 
        cut = cut + "* ( " + extraWeight + ")"
     self.ttree.Project(name, var, cut, options) 

     return h

class Block:
   'Common base class for all Sample Blocks'

   def __init__(self, name, label, color, isdata, doKfactor):
      self.name  = name
      self.color = color
      self.isData = isdata
      self.doKfactor = doKfactor
      self.label = label
      self.samples = []

   def printBlock(self):

      print "####################"
      print "Block Name: ", self.name
      print "Block Color: ", self.color
      print "Block IsData: ", self.isData
      print "####################"
      print "This block contains the following Samples"

      for l in self.samples:
        l.printSample()
     

   def addSample(self, s):
      self.samples.append(s)

   def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel, extraWeight,doKFactorGENVar,forceDataWeight):
     for _is,s in enumerate(self.samples):
       
       AuxName = "auxT1_sample" + s.name
       haux = s.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel, extraWeight, doKFactorGENVar, forceDataWeight)
       if not _is:
          h = haux.Clone(name+'_blockHisto')
       else:
          h.Add(haux)
       del haux

     h.SetLineColor(self.color)
     h.SetMarkerColor(self.color)
     h.GetYaxis().SetTitle('Events')
     h.SetTitle(self.label)

     return h

   def getTH2F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel, extraWeight):
     if(xmin == xmax) and (ymax == ymin):
        h = TH2F(name, "", len(nbinx)-1, array('d', nbinx),len(nbiny)-1, array('d', nbiny))
     elif (xmin == xmax):
        h = TH2F(name, "", len(nbinx)-1, array('d', nbinx),nbiny,ymin,ymax)
     elif (ymin == ymax):
        h = TH2F(name, "", nbinx,xmin,xmax,len(nbiny)-1, array('d', nbiny))
     else: 
        h = TH2F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax)

     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     
     for s in self.samples:
     
       AuxName = "auxT2_block" + s.name
       haux = s.getTH2F(lumi, AuxName, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel, extraWeight)
       h.Add(haux)
       del haux

     return h   

   def getTH3F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax, cut, options, xlabel, ylabel, zlabel, extraWeight):
     if(xmin == xmax) and (ymax == ymin) and (zmax == zmin):
        h = TH3F(name, "", len(nbinx)-1, array('d', nbinx),len(nbiny)-1, array('d', nbiny),len(nbinz)-1, array('d', nbinz))
     else: 
        h = TH3F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax)

     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     h.GetZaxis().SetTitle(zlabel)
     
     for s in self.samples:
     
       AuxName = "auxT3_block" + s.name
       haux = s.getTH3F(lumi, AuxName, var, nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax, cut, options, xlabel, ylabel, zlabel,extraWeight)
       h.Add(haux)
       del haux

     return h   

       

class Tree:
   'Common base class for a physics meaningful tree'

   def __init__(self, fileNameAndMap, name, isdata, isScan = False, isOnEOS = 0):
      self.name  = name
      self.isData = isdata
      self.blocks = []
      self.isScan = isScan
      self.isOnEOS = isOnEOS
      self.parseFileName(fileNameAndMap)
      self.pool = Pool(80)
      self.split = 1 

   def parseFileName(self, fileNameAndMap):

      fileName, defMap = fileNameAndMap
      f = open(fileName)

      for l in f.readlines():
        if (l[0] == "#" or len(l) < 2):
          continue
        l = l.format(**defMap)
        splitedLine = str.split(l)
        block       = splitedLine[0]
        theColor    = splitedLine[1]
        name        = splitedLine[2]
        label       = splitedLine[3]
        flocation   = splitedLine[4]
        xsection    = splitedLine[5]
        isdata      = int(splitedLine[6])
        doKfactor   = int(splitedLine[7])

        #color = 0
        # plusposition = theColor.find("+")
        # if(plusposition == -1):
        #   color = eval(theColor)
        # else:
        #   color = eval(theColor[0:plusposition])
        #   color = color + int(theColor[plusposition+1:len(theColor)])

        color = eval(theColor)

        sample = Sample(name, flocation, xsection, isdata, doKfactor, self.isScan, self.isOnEOS)
        coincidentBlock = [l for l in self.blocks if l.name == block]
        #print name
        if(coincidentBlock == []):

          newBlock = Block(block, label, color, isdata, doKfactor)
          newBlock.addSample(sample)
          self.addBlock(newBlock)

        else:

          coincidentBlock[0].addSample(sample)





   def printTree(self):

      print "######"
      print "Tree Name: ", self.name
      print "Tree IsData: ", self.isData
      print "######"
      print "This Tree contains the following Blocks"

      for l in self.blocks:
        l.printBlock()
     

   def addBlock(self, b):
      self.blocks.append(b)



   def getYields(self, lumi, var, xmin, xmax, cut, forceDataWeight=False):
  
      h = self.getTH1F(lumi, "yields", var, 1, xmin, xmax, cut, "", "",forceDataWeight)
      nbinmin = h.FindBin(xmin)
      nbinmax = h.FindBin(xmax)
      error = r.Double()
      value = h.IntegralAndError(nbinmin, nbinmax, error)
      y = [value, error]
      
      del h
      return y




      
   def getStack(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel, weight, kfactor, forceDataWeight=False):
   
      hs = THStack(name, "")
      tasks = [] 
      for b in self.blocks:
         for s in  b.samples:
            AuxName = "auxStack_block_" + name + "_" + b.name + '_' + s.name
            
            if self.split: 
               chunks = int(math.ceil(s.ttree.GetEntries() / float(chunkSize)))
               for i in range(chunks):
                  tasks.append( (s, b.name, (lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel, weight, kfactor, i,forceDataWeight)))
            else:
               tasks.append( (s, b.name, (lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel, weight, kfactor,forceDataWeight))  )
         
      temp_results = self.pool.map_async( runTH1Tasks, tasks,1 ) 
      while not temp_results.ready():
         print "Jobs left: %d out of %d"%(temp_results._number_left, len(tasks))
         #print tasks
         time.sleep(5)
      results = temp_results.get()
      for b in self.blocks:
         b.first = True
         for bname, hist in results:
            if b.name == bname: 
               if b.first:
                  mergedHist = hist.Clone( hist.GetName() + '_blockHisto')
                  b.first = False
               else:
                  mergedHist.Add( hist )
            del hist
         mergedHist.SetLineColor(b.color)
         mergedHist.SetFillColor(b.color)
         mergedHist.SetMarkerColor(b.color)
         mergedHist.GetYaxis().SetTitle('Events')
         mergedHist.SetTitle(b.label)
         for i in range(1,mergedHist.GetNbinsX()+1):
            if mergedHist.GetBinContent(i) < 0 : mergedHist.SetBinContent(i,0)
         hs.Add(mergedHist)
         del mergedHist
      can_aux = TCanvas("can_%s_%s"%(name, b.name))
      can_aux.cd()
      hs.Draw()

      del can_aux
      ylabel = "Events"
      if xmax != xmin:
         hs.GetXaxis().SetTitle(xlabel)
         b = int((xmax-xmin)/nbin)
         ylabel = "Events / " + str(b) + " GeV"
      else:     
         ylabel = "Events"
   
      hs.GetYaxis().SetTitle(ylabel)
      return hs   

   def getYieldTable(self, lumi, cut, weight, kfactor,forceDataWeight=False):
  
      h = self.getStack(lumi, 'yield_table', '1', 1, 0, 2, cut, '', '', weight, kfactor,forceDataWeight)
      table = {} 
      hists = [] 
      for hist in h.GetHists():
         hists.append(hist) # so we can erase them through the stack
         title = hist.GetTitle()
         error = r.Double()
         value = hist.IntegralAndError(0, 2, error)
         table[title]  = (value, error)
      SetOwnership(h,0)
      return table



   def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel, extraWeight, doKFactorGENVar,forceDataWeight=False):
    
      hs = self.getStack(lumi, name, var, nbin, xmin, xmax, cut, options, xlabel, extraWeight, doKFactorGENVar,forceDataWeight)
      h = hs.GetStack().Last().Clone(name + '_treeHisto')
      del hs 
      return h 

     # for ib,b in enumerate(self.blocks):
       
     #   AuxName = "auxh1_block_" + name + "_" + b.name
     #   haux = b.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel, extraWeight, doKFactorGENVar)
     #   if not ib:
     #      h = haux.Clone(name+'_treeHisto')
     #   else:
     #      h.Add(haux)
     #   del haux
       
     # return h

   def getTH2F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, extraWeight='1'):
     if(xmin == xmax) and (ymax == ymin):
        h = TH2F(name, "", len(nbinx)-1, array('d', nbinx),len(nbiny)-1, array('d', nbiny))
     elif (xmin == xmax):
        h = TH2F(name, "", len(nbinx)-1, array('d', nbinx),nbiny,ymin,ymax)
     elif (ymin == ymax):
        h = TH2F(name, "", nbinx,xmin,xmax,len(nbiny)-1, array('d', nbiny))
     else: 
        h = TH2F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax)
        
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     
     for b in self.blocks:
     
       AuxName = "aux_block" + name + "_" + b.name
       haux = b.getTH2F(lumi, AuxName, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel, extraWeight)
       h.Add(haux)
       del haux

     return h   

   def getTH3F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax, cut, options, xlabel, ylabel, zlabel, extraWeight='1'):
     if(xmin == xmax) and (ymax == ymin) and (zmax == zmin):
        h = TH3F(name, "", len(nbinx)-1, array('d', nbinx),len(nbiny)-1, array('d', nbiny), len(nbinz)-1, array('d', nbinz))
     else: 
        h = TH3F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax)
        
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     h.GetZaxis().SetTitle(zlabel)
     
     for b in self.blocks:
     
       AuxName = "aux_block" + name + "_" + b.name
       haux = b.getTH3F(lumi, AuxName, var, nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax, cut, options, xlabel, ylabel, zlabel, extraWeight)
       h.Add(haux)
       del haux

     return h   

