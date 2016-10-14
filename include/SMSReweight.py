import ROOT as r
from array import array
from ROOT import TTree, TFile, TCut, TH1F, TH2F, TH3F, THStack, TCanvas
import Sample, helper

class SMSReweight:
    'Provides weight for SMS scans'

    def __init__(self, smsName, sampleFile='samples.dat',smsXSection=""):
        self.smsName    = smsName
        self.sampleFile = sampleFile
        self.SetDataset()
        if smsXSection=="": self.defaultSMSXsection()
        else: self.smsXSection = smsXSection
        self.GetSMSCounter()

    def SetDataset(self):
        if self.smsName == 'T6bbllslepton':
            self.datasets = ['SMS_T6bbllslepton_mSbottom400to575_mLSP150to550',
                             'SMS_T6bbllslepton_mSbottom600to775_mLSP150to725',
                             'SMS_T6bbllslepton_mSbottom800to950_mLSP150to900']
        else:
            print '[E][SMSReWeight]: SMS not implemented'

    def defaultSMSXsection(self):
        if self.smsName == 'T6bbllslepton':
            self.xsecFile = 'datacards/sbottomXsec.txt'
        else:
            print '[E][SMSReWeight]: SMS not implemented'
        
    def GetSMSCounter(self):
        self.SetDataset()
        self.tree = Sample.Tree(helper.selectSamples(self.sampleFile, self.datasets, 'SIG'), 'SIG'  , 0, isScan = True)
        self.ngen = self.tree.blocks[0].samples[0].smsCount
        for ind, i in enumerate(self.tree.blocks[0].samples):
            if ind: self.ngen.Add(i.smsCount,1.)

    def getReWeightFactor(self, mSMS1, mSMS2):
        xsecf  = open(self.xsecFile,'r')
        xsecs  = eval(xsecf.read())
        xsec   = xsecs[mSMS1][0]
        ngen   = self.ngen.GetBinContent(self.ngen.FindBin(mSMS1,mSMS2,1))
        print xsec, ngen
        return xsec/ngen

    def reWeightHisto(self, mSMS1, mSMS2, histo):
        histoRW = histo.Clone(histo.GetName()+"_RW")
        histoRW.Scale( getReWeightFactor(mSMS1, mSMS2) )
        return histoRW

