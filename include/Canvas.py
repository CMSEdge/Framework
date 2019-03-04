from ROOT import TCanvas, TLegend,TPie,  TPad, TLine, TLatex, TGraphAsymmErrors, TH1F, THStack, TGraphErrors, TLine, TPaveStats, TGraph, TArrow
import ROOT as r
import os, copy, math, array
from array import array
class Canvas:
   'Common base class for all Samples'

   def __init__(self, name, _format, x1, y1, x2, y2, ww=0, hh=0):
      self.name = name
      self.format = _format
      self.plotNames    = [name + "." + i for i in _format.split(',')]
      self.plotNamesLog = [name + "_log." + i for i in _format.split(',')]
      self.myCanvas = TCanvas(name, name) if not ww else TCanvas(name, name, ww, hh)
      self.ToDraw = []
      self.orderForLegend = []
      self.histos = []
      self.lines = []
      self.arrows= []
      self.latexs= []
      self.bands = []
      self.options = []
      self.labels = []      
      self.labelsOption = []
      self.myLegend = TLegend(x1, y1, x2, y2)
      self.myLegend.SetFillColor(r.kWhite)
      self.myLegend.SetTextFont(42)
      self.myLegend.SetTextSize(0.04)
      self.myLegend.SetLineWidth(0)
      self.myLegend.SetBorderSize(0)              
      r.gStyle.SetPadRightMargin(0.05)

   def changeLabelsToNames(self):
      newlabels = []
      for il,lab in enumerate(self.labels):
         print 'changing label %s to %s'%(lab, self.histos[il].GetName())
         newlabels.append(self.histos[il].GetName())
      self.labels = newlabels

   def banner(self, isData, lumi):
     
      latex = TLatex()
      latex.SetNDC();                         
      latex.SetTextAngle(0);                  
      latex.SetTextColor(r.kBlack);           
      latex.SetTextFont(42);                  
      latex.SetTextAlign(31);                 
      latex.SetTextSize(0.06);                
      latex.DrawLatex(0.25, 0.93, "#bf{CMS}") 
               
      latexb = TLatex()                      
      latexb.SetNDC();
      latexb.SetTextAngle(0);
      latexb.SetTextColor(r.kBlack);
      latexb.SetTextFont(42);
      latexb.SetTextAlign(31);
      latexb.SetTextSize(0.04);            

      latexd = TLatex()
      latexd.SetNDC();
      latexd.SetTextAngle(90);
      latexd.SetTextColor(r.kBlack);
      latexd.SetTextFont(42);
      latexd.SetTextAlign(31);
      latexd.SetTextSize(0.06);
      #latexd.DrawLatex(0.15, 0.93, "Events")  

      #if(isData):
      #  latexb.DrawLatex(0.44, 0.93, "#it{Preliminary}")
      #else:
      #  latexb.DrawLatex(0.44, 0.93, "#it{Simulation}")

      text_lumi = str(lumi) + " fb^{-1} (13 TeV)"
      latexc = TLatex()
      latexc.SetNDC();
      latexc.SetTextAngle(0);
      latexc.SetTextColor(r.kBlack);
      latexc.SetTextFont(42);
      latexc.SetTextAlign(31);
      latexc.SetTextSize(0.04);
      latexc.DrawLatex(0.90, 0.93, text_lumi)                

   def banner2(self, isData, lumi):
    
      latex = TLatex()
      latex.SetNDC();
      latex.SetTextAngle(0);
      latex.SetTextColor(r.kBlack);
      latex.SetTextFont(42);
      latex.SetTextAlign(31);
      latex.SetTextSize(0.06);
      latex.DrawLatex(0.23, 0.93, "#bf{CMS}")

      latexb = TLatex()
      latexb.SetNDC();
      latexb.SetTextAngle(0);
      latexb.SetTextColor(r.kBlack);
      latexb.SetTextFont(42);
      latexb.SetTextAlign(31);
      latexb.SetTextSize(0.04);            

      latexd = TLatex()
      latexd.SetNDC();
      latexd.SetTextAngle(90);
      latexd.SetTextColor(r.kBlack);
      latexd.SetTextFont(42);
      latexd.SetTextAlign(31);
      latexd.SetTextSize(0.04);          
      #latexd.DrawLatex(0.07, 0.93, "Events")  

      #if(isData):
      #  latexb.DrawLatex(0.38, 0.93, "#it{Preliminary}")
      #else:
      latexb.DrawLatex(0.415, 0.93, "#it{Supplementary}")
      


      latexe = TLatex()
      latexe.SetNDC();
      latexe.SetTextColor(r.kBlack);
      latexe.SetTextFont(42);
      latexe.SetTextAlign(31);
      latexe.SetTextSize(0.03);          
      latexe.DrawLatex(0.33, 0.86, "arXiv: 1806.05264")

      text_lumi = str(lumi) + " fb^{-1} (13 TeV)"
      latexc = TLatex()
      latexc.SetNDC();
      latexc.SetTextAngle(0);
      latexc.SetTextColor(r.kBlack);
      latexc.SetTextFont(42);
      latexc.SetTextAlign(31);
      latexc.SetTextSize(0.05);
      latexc.DrawLatex(0.90, 0.93, text_lumi)

   def banner3(self, isData, lumi):
     
      latex = TLatex()
      latex.SetNDC();                         
      latex.SetTextAngle(0);                  
      latex.SetTextColor(r.kBlack);           
      latex.SetTextFont(42);                  
      latex.SetTextAlign(31);                 
      latex.SetTextSize(0.07);                
      latex.DrawLatex(0.1, 1.22, "#bf{CMS}") 
               
      latexb = TLatex()                      
      latexb.SetNDC();
      latexb.SetTextAngle(0);
      latexb.SetTextColor(r.kBlack);
      latexb.SetTextFont(42);
      latexb.SetTextAlign(31);
      latexb.SetTextSize(0.05);            
                                                             
      #if(isData):
      #  latexb.DrawLatex(0.34, 1.22, "#it{Preliminary}")
      #else:
      #  latexb.DrawLatex(0.34, 1.22, "#it{Simulation}")
                                                             
      text_lumi = str(lumi) + " fb^{-1} (13 TeV)"
      latexc = TLatex()
      latexc.SetNDC();
      latexc.SetTextAngle(0);
      latexc.SetTextColor(r.kBlack);
      latexc.SetTextFont(42);
      latexc.SetTextAlign(31);
      latexc.SetTextSize(0.05);


   def addBand(self, x1, y1, x2, y2, color, opacity):

      grshade = TGraph(4)
      grshade.SetPoint(0,x1,y1)
      grshade.SetPoint(1,x2,y1)
      grshade.SetPoint(2,x2,y2)
      grshade.SetPoint(3,x1,y2)
      #grshade.SetFillStyle(3001)
      grshade.SetFillColorAlpha(color, opacity)
      self.bands.append(grshade)

   def addLine(self, x1, y1, x2, y2, color, thickness = 0.):
      line = TLine(x1,y1,x2,y2)
      line.SetLineColor(color)
      line.SetLineStyle(2)
      if thickness:
          line.SetLineWidth(thickness)
      self.lines.append(line)

   def addArrow(self, x1, y1, x2, y2, color, option, thickness = 0.):
      arrow = TArrow(x1,y1,x2,y2, 0.05, option)
      arrow.SetLineColor(color)
      if thickness:
          arrow.SetLineWidth(thickness)
      self.arrows.append(arrow)

   def addLatex(self, x1, y1, text, font=42, size = 0.04):
      lat = [x1, y1, text, font, size]
      self.latexs.append(lat)

   def makeOFHisto(self, histo):
      nbinsx = histo.GetNbinsX()
      xmin = histo.GetXaxis().GetXmin(); xmax = histo.GetXaxis().GetXmax();
      newhisto = r.TH1F(histo.GetName() +'_withOFBin', histo.GetTitle()+'_withOFBin', nbinsx+1, xmin, xmax+(xmax-xmin)/nbinsx)
      newhisto.Sumw2()
      newhisto.SetMarkerColor(histo.GetMarkerColor())
      newhisto.SetMarkerStyle(histo.GetMarkerStyle())
      newhisto.SetMarkerSize (histo.GetMarkerSize ())
      newhisto.SetLineColor(histo.GetLineColor())
      newhisto.SetLineStyle(histo.GetLineStyle())
      newhisto.SetMaximum(histo.GetMaximum())
      for i in range(1,nbinsx+2):
         newhisto.SetBinContent(i,histo.GetBinContent(i))
         newhisto.SetBinError  (i,histo.GetBinError  (i))
      return newhisto
        
 
   def addHisto(self, h, option, label, labelOption, color, ToDraw, orderForLegend, doOF = False):

      if(color != ""):
          h.SetLineColor(color)
          h.SetMarkerColor(color)
      if(label == ""):
          label = h.GetTitle()

      self.histos.append(h if not doOF else self.makeOFHisto(h))
      self.options.append(option)
      self.labels.append(label)
      self.labelsOption.append(labelOption)
      self.ToDraw.append(ToDraw)
      self.orderForLegend.append(orderForLegend)

   def addGraph(self, h, option, label, labelOption, color, ToDraw, orderForLegend):

      if(color != ""):
          h.SetLineColor(color)
          h.SetMarkerColor(color)
      if(label == ""):
          label = h.GetTitle()

      self.histos.append(h)
      self.options.append(option)
      self.labels.append(label)
      self.labelsOption.append(labelOption)
      self.ToDraw.append(ToDraw)
      self.orderForLegend.append(orderForLegend)


   def addStack(self, h, option, ToDraw, orderForLegend):

      legendCounter = orderForLegend
      if(orderForLegend < len(self.orderForLegend)):
          legendCounter = len(self.orderForLegend)

      self.addHisto(h, option, "", "", "", ToDraw, -1)  
      for h_c in h.GetHists():
          self.addHisto(h_c, "H", h_c.GetTitle(), "F", "", 0, legendCounter)
          legendCounter = legendCounter + 1                                          
       
   def addPies(self, h, option, ToDraw, orderForLegend):
   
      legendCounter = orderForLegend
      if(orderForLegend < len(self.orderForLegend)):
          legendCounter = len(self.orderForLegend)
   
      self.addHisto(h, option, "", "", "", ToDraw, -1)  
      for h_c in h.GetHists():
          self.addHisto(h_c, "H", h_c.GetTitle(), "F", "", 0, legendCounter)
          legendCounter = legendCounter + 1                                       

   def makeLegend(self):

      for i in xrange(len(self.histos)-1, -1, -1):
      #for i in range(0, len(self.histos)):
          print "i", i
          for j in xrange(len(self.orderForLegend)-1, -1, -1):
              #for j in range(0, len(self.orderForLegend)):
              if(self.orderForLegend[j] == i):
              #if(self.orderForLegend[j] != -1 and self.orderForLegend[j] == i):
                  self.myLegend.AddEntry(self.histos[j], self.labels[j], self.labelsOption[j])
                  #self.myLegend.AddEntry(self.histos[j], self.labels[j], self.labelsOption[j])
          

   def ensurePath(self, _path):
      d = os.path.dirname(_path)
      if not os.path.exists(d):
         os.makedirs(d)                 

   def saveRatio(self, legend, isData, log, lumi, hdata, hMC, r_ymin=0, r_ymax=2, label ="Data/Prediction"):

      self.myCanvas.cd()

      pad1 = TPad("pad1", "pad1", 0, 0.25, 1, 1.0)
      pad1.SetBottomMargin(0.12)
      pad1.Draw()                                     
      pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.25)
      pad2.SetTopMargin(0.1);
      pad2.SetBottomMargin(0.3);
      pad2.Draw();                                      

      pad1.cd()
      if(log):
          pad1.SetLogy(1)

      for i in range(0, len(self.histos)):
          if(self.ToDraw[i] != 0):
              #self.histos[i].SetMinimum(0.00001)
              #self.histos[i].SetMinimum(0.1)
              #self.histos[i].SetMinimum(0.0001)
              #self.histos[i].SetMaximum(0.5)
              self.histos[i].Draw(self.options[i])

      if(legend):
          self.makeLegend()
          self.myLegend.Draw()

      for band in self.bands:
          band.Draw('f')
  
      for line in self.lines:
          line.Draw()
  
      for arrow in self.arrows:
          arrow.Draw()
  
      for latex in self.latexs:
          lat = TLatex()
          lat.SetNDC()
          lat.SetTextSize(latex[-1])
          lat.SetTextFont(latex[-2])
          lat.DrawLatex(latex[0], latex[1], latex[2])
  
      
      if type(hMC) != list:
          hMClist = [hMC]
      else: hMClist = hMC

      ratios = []

      for tmp_hMC in hMClist:
          ind = hMClist.index(tmp_hMC)
          tmp_ratio = hdata.Clone(tmp_hMC.GetName()+'_ratio')
          tmp_ratio.Divide(tmp_hMC)

          tmp_ratio.SetTitle("")
          tmp_ratio.GetYaxis().SetRangeUser(r_ymin, r_ymax);
          tmp_ratio.GetYaxis().SetTitle(label);
          tmp_ratio.GetYaxis().CenterTitle();
          tmp_ratio.GetYaxis().SetLabelSize(0.12);
          tmp_ratio.GetXaxis().SetLabelSize(0.12);
          tmp_ratio.GetYaxis().SetTitleOffset(0.3);
          tmp_ratio.GetYaxis().SetNdivisions(4);
          tmp_ratio.GetYaxis().SetTitleSize(0.12);
          tmp_ratio.GetXaxis().SetTitleSize(0.12);
          tmp_ratio.GetXaxis().SetLabelOffset(0.08);
          tmp_ratio.GetXaxis().SetTitle('');
          tmp_ratio.SetMarkerStyle(tmp_hMC.GetMarkerStyle());
          tmp_ratio.SetFillColorAlpha(r.kBlue+3,1.0)
          tmp_ratio.SetFillStyle(3017)
          tmp_ratio.SetMarkerColor(r.kBlack);
          #tmp_ratio.SetMarkerSize(0.6);
          #tmp_ratio.SetMarkerSize(0.6);
          #tmp_ratio.SetMarkerColor(r.kBlack if len(hMClist) == 1 else tmp_hMC.GetMarkerColor());
          #tmp_ratio.SetLineColor  (r.kBlack if len(hMClist) == 1 else tmp_hMC.GetLineColor  ());
          tmp_ratio.SetLineColor  (tmp_hMC.GetLineColor());
          tmp_ratio.SetLineStyle(tmp_hMC.GetLineStyle())
          ratios.append(tmp_ratio)
          xmin = tmp_ratio.GetBinLowEdge(1)
          xmax = tmp_ratio.GetBinLowEdge(tmp_ratio.GetNbinsX()+1)

      #tmp_ratio.Draw("E,SAME");
      pad2.cd();  
      for rat in ratios:
          rat.Draw('PE2,same');

      line = TLine(xmin, 1, xmax, 1)
      line.SetLineColor(r.kGray+2);
      line.Draw('');

      pad1.cd()
      self.banner2(isData, lumi)
      for i,plotName in enumerate(self.plotNames):
          pad1.cd()
          pad1.SetLogy(0)
          path    = 'plots/'+plotName
          pathlog = 'plots/'+self.plotNamesLog[i]
          self.ensurePath(path)
          self.myCanvas.SaveAs(path)
          if not '.root' in pathlog:
              pad1.cd()
              pad1.SetLogy()
              self.myCanvas.SaveAs(pathlog)

          

      pad1.IsA().Destructor(pad1) 
      pad2.IsA().Destructor(pad2) 
      self.myLegend.IsA().Destructor(self.myLegend)
      self.myCanvas.IsA().Destructor(self.myCanvas)                                                                                                                                            


   def save(self, legend, isData, log, lumi, labelx, ymin=0, ymax=0):

      self.myCanvas.cd()
      
      if(log):
          self.myCanvas.GetPad(0).SetLogy(1)
     
      for i in range(0, len(self.histos)):
          if(self.ToDraw[i] != 0):        
              if ymin and ymax:
                  self.histos[i].GetYaxis().SetRangeUser(ymin, ymax)
              self.histos[i].Draw(self.options[i])

      for band in self.bands:
          band.Draw('f')
  
      for line in self.lines:
          line.Draw()
  
      for arrow in self.arrows:
          arrow.Draw()
  
      for latex in self.latexs:
          lat = TLatex()
          lat.SetNDC()
          lat.SetTextSize(latex[-1])
          lat.SetTextFont(latex[-2])
          lat.DrawLatex(latex[0], latex[1], latex[2])
  
      if(legend):
          self.makeLegend()
          self.myLegend.Draw()

      lat = TLatex()
      lat.SetNDC()
      lat.SetTextSize(0.05)
      lat.SetTextFont(42)
      lat.DrawLatex(0.46, 0.04, labelx)
      
      
      self.banner(isData, lumi)
      for plotName in self.plotNames:
          path = 'plots/'+plotName
          self.ensurePath(path)
          self.myCanvas.SaveAs(path)

      self.myLegend.IsA().Destructor(self.myLegend)
      self.myCanvas.IsA().Destructor(self.myCanvas)

   def savePie(self, legend, lumi, labelx):

      cpie = TCanvas("cpie","TPie test",700,700)
      
      pad1 = TPad("pad1", "pad1", 0.1, 0.1, 0.75, 0.75)
      pad1.SetBottomMargin(0.12)
      pad1.Draw()   
      
      pad1.cd()
    
      colors = []
      names = []
      vals = []
      
      for i in range(0, len(self.histos)):
          vals.append(self.histos[i].Integral())
          colors.append(self.histos[i].GetLineColor())
          names.append(self.histos[i].GetName())
      v = array('d', vals)
      c = array('i', colors)
      pie4 = TPie("p4","",len(v),v,c);
      
      pie4.SetRadius(.45);
      pie4.SetLabelsOffset(.01);
      pie4.SetLabelsOffset(100);
      pie4.SetEntryLineWidth(1,2);
      pie4.SetEntryLineWidth(2,2);
      pie4.SetEntryLineWidth(3,2);
      pie4.SetEntryLineWidth(4,2);
      pie4.SetEntryLineWidth(5,2);
      pie4.SetEntryLineWidth(6,2);
      pie4.SetEntryLineWidth(7,2);
      pie4.Draw();                                      

      self.makeLegend()
      self.myLegend.Draw()              

      lat = TLatex()
      lat.SetNDC()
      lat.SetTextSize(0.06)
      lat.SetTextFont(42)
      lat.DrawLatex(0.12, -0.1, "Slepton signal region, "+labelx)

      self.banner3(True, lumi)
      for plotName in self.plotNames:
          path = 'plots/'+plotName
          #self.ensurePath(path)
          cpie.SaveAs(path)

      self.myLegend.IsA().Destructor(self.myLegend)
      #cpie.Destructor(self.myCanvas)

 #del self.myCanvas                                                                                      
