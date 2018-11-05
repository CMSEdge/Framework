import ROOT, math, optparse, copy

##########33
## problemas e historias
##
#1) el breikindancing
#2) electrones: no hay mva tight vs miniiso 0.4 (solo loose)
#3) muones:     no hay ip3d < 8, solo ip3d < 4


# Electrons:
# GsfElectronToMVATightIDEmuTightIP2DSIP3D4
# MVAVLooseElectronToMini
# MVATightElectronToConvVetoIHit0
# Reco scale factors

# We need the "Emu" version on the first one for the trigger emulation cuts.

# Muons:
# Medium ID 
# miniIso0.2 vs Medium ID 
# dxy0.05, dz0.1 vs Medium ID 
# SIP3D<4 vs Medium ID 
# (reco scale factors when available)

def LoadLeptonSF():
    ROOT.gROOT.LoadMacro('include/LeptonSF.C+')
    # Electrons
    elFile = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/stuffForMoriond/LepSFs/ELEC/scaleFactors.root')
    ROOT.SetElecID (copy.deepcopy( elFile.Get('GsfElectronToMVATightIDEmuTightIP2DSIP3D4')))
    ROOT.SetElecIP (copy.deepcopy( elFile.Get('MVAVLooseElectronToMini'))) 
    ROOT.SetElecISO(copy.deepcopy( elFile.Get('MVATightElectronToConvVetoIHit0')))# its not iso, but who cares
    elFile.Close()
    elFile = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/stuffForMoriond/LepSFs/ELEC/egammaEffi.txt_EGM2D.root')
    ROOT.SetElecTrk(copy.deepcopy( elFile.Get('EGamma_SF2D')))

    # Muons
    muFile = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/stuffForMoriond/LepSFs/MUON/TnP_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root')
    ROOT.SetMuonID( copy.deepcopy( muFile.Get('SF') ))
    muFile.Close()
    muFile = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/stuffForMoriond/LepSFs/MUON/TnP_NUM_MiniIsoLoose_DENOM_MediumID_VAR_map_pt_eta.root')
    ROOT.SetMuonISO( copy.deepcopy( muFile.Get('SF') ))
    muFile.Close()
    muFile = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/stuffForMoriond/LepSFs/MUON/TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root')
    ROOT.SetMuonIP1( copy.deepcopy(muFile.Get('SF')))
    muFile.Close()
    muFile = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/stuffForMoriond/LepSFs/MUON/TnP_NUM_TightIP3D_DENOM_MediumID_VAR_map_pt_eta.root')
    ROOT.SetMuonIP2( copy.deepcopy(muFile.Get('SF')))
    muFile.Close()
    muFile = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/stuffForMoriond/LepSFs/MUON/Tracking_EfficienciesAndSF_BCDEFGH.root')
    ROOT.SetMuonTrk( copy.deepcopy(muFile.Get('ratio_eff_aeta_dr030e030_corr')))
    muFile.Close()

LoadLeptonSF()

if __name__ == "__main__":
    print 'El', ROOT.LepSF(35,-2.3,11 ,'ElUp')
    print 'El', ROOT.LepSF(35, 2.3,11 ,'ElUp')
    print 'El', ROOT.LepSF(500, 2.3,11,'ElUp')
    print 'Mu', ROOT.LepSF(35,-2.3,13 ,'ElUp')
    print 'Mu', ROOT.LepSF(35,-2.3,13 ,'ElDn')
    print 'Mu', ROOT.LepSF(35,-2.3,13 ,'')
    print 'Mu', ROOT.LepSF(35,-2.3,11 ,'MuUp')
    print 'Mu', ROOT.LepSF(35,-2.3,11 ,'MuDn')
    print 'Mu', ROOT.LepSF(35,-2.3,11 ,'')

