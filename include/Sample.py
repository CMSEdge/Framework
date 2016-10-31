import ROOT as r
from array import array
from ROOT import TTree, TFile, TCut, TH1F, TH2F, TH3F, THStack, TCanvas


class Sample:
   'Common base class for all Samples'

   def __init__(self, name, friendlocation, xsection, isdata, isScan):
      self.name = name
      self.location = friendlocation
      self.xSection = xsection
      self.isData = isdata
      ftfileloc = friendlocation+'/evVarFriend_'+self.name+'.root' 
      self.ftfile = TFile(ftfileloc)
      self.ttree = self.ftfile.Get('sf/t')
      self.isScan = isScan
      if not self.isData:
        gw = 0.
        for i in self.ttree:
            gw = abs(i.genWeight_Edge)
            if gw: break
        self.count = self.ftfile.Get('SumGenWeights').GetBinContent(1)/abs(gw)
      else:
        self.count = self.ftfile.Get('Count').GetBinContent(1)
      self.lumWeight =   1.0
      self.puWeight   = '1.0'
      self.SFWeight   = '1.0'
      self.btagWeight = '1.0'
      self.trigScale  = '1.0'
      ## no longer needed self.emSF_trig = 0.89
      ## no longer needed self.eeSF_trig = 0.97
      ## no longer needed self.mmSF_trig = 0.94
      if not self.isData:
        self.trigScale  = 'weight_trigger_Edge'
        self.lumWeight = self.xSection / self.count
        self.puWeight    = "PileupW_Edge"
        self.btagWeight  = "weight_btagsf_Edge"
        self.SFWeight    = "weight_LepSF_Edge"
      if self.isScan:
        self.trigScale  = 'weight_trigger_Edge'
        self.lumWeight  =  1.0
        self.puWeight   = "PileupW_Edge"
        self.btagWeight = "weight_btagsf_Edge"
        self.SFWeight    = "weight_LepSF_Edge"
        self.smsCount =  self.ftfile.Get('CountSMS')
   def printSample(self):
      print "#################################"
      print "Sample Name: ", self.name
      print "Sample Location: ", self.location
      print "Sample XSection: ", self.xSection
      print "Sample IsData: ", self.isData
      print "Sample LumWeight: ", self.lumWeight
      print "#################################"


   def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel, ofBin=True):
     
      if(xmin == xmax):
        _nbins = len(nbin)-1
        _arr = array('d', nbin)
        h = TH1F(name+'_noOF', "", _nbins, _arr)
        _newarr = _arr + array('d', [ 2*_arr[-1]-_arr[-2] ])
        h_of = TH1F(name, "", _nbins+1, _newarr)
        ylabel = "# events"
      else:
        h = TH1F(name+'_noOF', "", nbin, xmin, xmax)
        bw = int((xmax-xmin)/nbin)
        ylabel = "Events / " + str(bw) + " GeV"
        h_of = TH1F(name, '', nbin+1, xmin, xmax+bw)

      h.Sumw2()
      h.GetXaxis().SetTitle(xlabel)
      h.GetYaxis().SetTitle(ylabel)

      h_of.Sumw2()
      h_of.GetXaxis().SetTitle(xlabel)
      h_of.GetYaxis().SetTitle(ylabel)

      addCut = ""
      if self.isData:
        if(name.find("DoubleMuon") != -1):
          addCut = "(!((Lep1_pdgId_Edge * Lep2_pdgId_Edge == -121) || (Lep1_pdgId_Edge * Lep2_pdgId_Edge == -143)))"
          cut = cut + "* ( " + addCut + " ) "
        if(name.find("DoubleEG") != -1):
          addCut = "(!((Lep1_pdgId_Edge * Lep2_pdgId_Edge == -169) || (Lep1_pdgId_Edge * Lep2_pdgId_Edge == -143)))"
          cut = cut + "* ( " + addCut + " ) "
        if(name.find("MuonEG") != -1):
          addCut = "(!((Lep1_pdgId_Edge * Lep2_pdgId_Edge == -121) || (Lep1_pdgId_Edge * Lep2_pdgId_Edge == -169)))"
          cut = cut + "* ( " + addCut + " ) "
           
      if(self.isData == 0):
        cut = cut + "* ( " + str(self.lumWeight*lumi) + " * genWeight_Edge/abs(genWeight_Edge) * " + self.puWeight + " * " + self.SFWeight + " * " + self.btagWeight + " * " + self.trigScale + " )" 
      
      self.ttree.Project(h.GetName(), var, cut, options) 

      for _bin in range(1, h.GetNbinsX()+2):
          h_of.SetBinContent(_bin, h.GetBinContent(_bin))
          h_of.SetBinError  (_bin, h.GetBinError  (_bin))
      
      return (h_of if ofBin else h)

   def getTH2F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel, extraWeight='1'):
   
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
<<<<<<< HEAD
        cut = cut + "* ( " + str(self.lumWeight*lumi) + " * genWeight_Edge/abs(genWeight_Edge) * " + self.puWeight + " * " + self.SFWeight +  " * " + self.btagWeight + " * " + self.trigScale + " )" 
=======
        cut = cut + "* ( " + str(self.lumWeight*lumi) + " * genWeight_Edge/abs(genWeight_Edge) * " + self.puWeight + " * " + self.SFWeight +  " * " + self.btagWeight + " * " + extraWeight +" )" 
>>>>>>> Systematics partially implemented
     
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
<<<<<<< HEAD
        cut = cut + "* ( " + str(self.lumWeight*lumi) + " * genWeight_Edge/abs(genWeight_Edge) * " + self.puWeight + " * " + self.SFWeight +  " * " + self.btagWeight + " * " + self.trigScale + " )" 
=======
        cut = cut + "* ( " + str(self.lumWeight*lumi) + " * genWeight_Edge/abs(genWeight_Edge) * " + self.puWeight + " * " + self.SFWeight +  " * " + self.btagWeight + " * " + extraWeight + " )" 
>>>>>>> Systematics partially implemented
     
     self.ttree.Project(name, var, cut, options) 
     return h

class Block:
   'Common base class for all Sample Blocks'

   def __init__(self, name, label, color, isdata):
      self.name  = name
      self.color = color
      self.isData = isdata
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

   def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel, ofBin = True):

     for _is,s in enumerate(self.samples):
       AuxName = "auxT1_sample" + s.name
       haux = s.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel, ofBin)
       if not _is:
          h = haux.Clone(name+'_blockHisto')
       else:
          h.Add(haux)
       del haux

     h.SetLineColor(self.color)
     h.SetMarkerColor(self.color)
     h.SetTitle(self.label)

     return h

   def getTH2F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel):
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
       haux = s.getTH2F(lumi, AuxName, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel)
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

   def __init__(self, fileName, name, isdata, isScan = False):
      self.name  = name
      self.isData = isdata
      self.blocks = []
      self.isScan = isScan
      self.parseFileName(fileName)

   def parseFileName(self, fileName):

      f = open(fileName)

      for l in f.readlines():

        if (l[0] == "#" or len(l) < 2):
          continue

        splitedLine = str.split(l)
        block       = splitedLine[0]
        theColor    = splitedLine[1]
        name        = splitedLine[2]
        label       = splitedLine[3]
        flocation   = splitedLine[4]
        xsection    = float(splitedLine[5])
        isdata      = int(splitedLine[6])

        color = 0
        plusposition = theColor.find("+")
        if(plusposition == -1):
          color = eval(theColor)
        else:
          color = eval(theColor[0:plusposition])
          color = color + int(theColor[plusposition+1:len(theColor)])

        sample = Sample(name, flocation, xsection, isdata, self.isScan)
        coincidentBlock = [l for l in self.blocks if l.name == block]

        if(coincidentBlock == []):

          newBlock = Block(block, label, color, isdata)
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



   def getYields(self, lumi, var, xmin, xmax, cut):
  
      h = self.getTH1F(lumi, "yields", var, 1, xmin, xmax, cut, "", "")
      nbinmin = h.FindBin(xmin)
      nbinmax = h.FindBin(xmax)
      error = r.Double()
      value = h.IntegralAndError(nbinmin, nbinmax, error)
      y = [value, error]
      
      del h
      return y

   def getStack(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel):
   

     hs = THStack(name, "")

     for b in self.blocks:
     
       AuxName = "auxStack_block_" + name + "_" + b.name
       haux = b.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel)
       haux.SetFillColor(b.color)
       hs.Add(haux)
       del haux


     can_aux = TCanvas("can_%s_%s"%(name, b.name))
     can_aux.cd()
     hs.Draw()

     del can_aux

     if xmax != xmin:
       hs.GetXaxis().SetTitle(xlabel)
       b = int((xmax-xmin)/nbin)
       ylabel = "Events / " + str(b) + " GeV"
     else:     
       ylabel = "# events"
   
     hs.GetYaxis().SetTitle(ylabel)
     return hs   


   def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel, ofBin = True):
     
     for ib,b in enumerate(self.blocks):
       AuxName = "auxh1_block_" + name + "_" + b.name
       haux = b.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel, ofBin)
       if not ib:
          h = haux.Clone(name+'_treeHisto')
       else:
          h.Add(haux)
       del haux

     return h

   def getTH2F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel):
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
       haux = b.getTH2F(lumi, AuxName, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel)
       h.Add(haux)
       del haux

     return h   

   def getTH3F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax, cut, options, xlabel, ylabel, zlabel, extraWeight):
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

