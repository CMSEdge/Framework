import ROOT, math, optparse, copy


def LoadFastSimSF():


    ROOT.gROOT.LoadMacro('include/FastSimSF.C+')
    # Electrons
    elFile = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/stuffForMoriond/LepSFs/Full_Fast/sf_el_tight2d3d.root')
    ROOT.SetFSElecID (copy.deepcopy( elFile.Get('histo2D'))); elFile.Close()
    elFile = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/stuffForMoriond/LepSFs/Full_Fast/sf_el_mini01.root')    
    ROOT.SetFSElecIP (copy.deepcopy( elFile.Get('histo2D'))) ; elFile.Close()
    elFile = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/stuffForMoriond/LepSFs/Full_Fast/sf_el_inhit_eq0.root')
    ROOT.SetFSElecISO(copy.deepcopy( elFile.Get('histo2D')))# its not iso, but who cares
    elFile.Close()

    # Muons
    muFile = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/stuffForMoriond/LepSFs/Full_Fast/sf_mu_mediumID.root')
    ROOT.SetFSMuonID( copy.deepcopy( muFile.Get('histo2D') ))
    muFile.Close()
    muFile = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/stuffForMoriond/LepSFs/Full_Fast/sf_mu_mediumID_mini02.root')
    ROOT.SetFSMuonISO( copy.deepcopy( muFile.Get('histo2D') ))
    muFile.Close()
    muFile = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/stuffForMoriond/LepSFs/Full_Fast/sf_mu_mediumID_tightIP2D.root')
    ROOT.SetFSMuonIP1( copy.deepcopy(muFile.Get('histo2D')))
    muFile.Close()
    muFile = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/stuffForMoriond/LepSFs/Full_Fast/sf_mu_mediumID_tightIP3D.root')
    ROOT.SetFSMuonIP2( copy.deepcopy(muFile.Get('histo2D')))
    muFile.Close()

LoadFastSimSF()

if __name__ == "__main__":
    print 'El', ROOT.LepSFFastSim(35,-2.3,11 ,'ElUp')
    print 'El', ROOT.LepSFFastSim(35, 2.3,11 ,'ElUp')
    print 'El', ROOT.LepSFFastSim(500, 2.3,11,'ElUp')
    print 'Mu', ROOT.LepSFFastSim(35,-2.3,13 ,'ElUp')
    print 'Mu', ROOT.LepSFFastSim(35,-2.3,13 ,'ElDn')
    print 'Mu', ROOT.LepSFFastSim(35,-2.3,13 ,'')
    print 'Mu', ROOT.LepSFFastSim(35,-2.3,11 ,'ElUp')
    print 'Mu', ROOT.LepSFFastSim(35,-2.3,11 ,'ElDn')
    print 'Mu', ROOT.LepSFFastSim(35,-2.3,11 ,'')

