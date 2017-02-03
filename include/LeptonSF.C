#include<iostream>
#include"TH1.h"
#include"TH2.h"


using namespace std;

TH2F* hElecID  = 0;
TH2F* hElecIP  = 0;
TH2F* hElecISO = 0;
TH2F* hElecTrk = 0;
TH2F* hMuonID  = 0;
TH2F* hMuonIP1  = 0;
TH2F* hMuonIP2  = 0;
TH2F* hMuonISO = 0;
TH2F* hMuonTrk = 0;


void SetElecID (TH2F* h){ hElecID  = h; }
void SetElecIP (TH2F* h){ hElecIP  = h; }
void SetElecISO(TH2F* h){ hElecISO = h; }
void SetElecTrk(TH2F* h){ hElecTrk = h; }
void SetMuonID (TH2F* h){ hMuonID  = h; }
void SetMuonIP1 (TH2F* h){ hMuonIP1  = h; }
void SetMuonIP2 (TH2F* h){ hMuonIP2  = h; }
void SetMuonISO(TH2F* h){ hMuonISO = h; }
void SetMuonTrk(TH2F* h){ hMuonTrk = h; }

Double_t _getIDoverRECO(Double_t pt, Double_t eta, Int_t pdgId, TString sys)
{

  eta = TMath::Abs(eta);
  Double_t sf   = 0.;
  Double_t sf_e = 0.;
  if (TMath::Abs(pdgId) == 13){
    if (pt > 120) pt = 119.9;
    sf   = hMuonID->GetBinContent( hMuonID->FindBin(pt,eta));
    // error hardcoded to 3 later on
    // sf_e = (sys.Contains("Mu")) ? hMuonID->GetBinError  ( hMuonID->FindBin(pt,eta)) : 0.;
    sf_e = 0;
  }
  else{
    if (pt > 200) pt = 199.9;
    sf   = hElecID->GetBinContent( hElecID->FindBin(pt,eta));
    sf_e = hElecID->GetBinError  ( hElecID->FindBin(pt,eta));
    sf_e = (sys.Contains("El")) ? hElecID->GetBinError  ( hElecID->FindBin(pt,eta)) : 0.;
  }
  if (sys.Contains("Up"))     return sf+sf_e;
  else if(sys.Contains("Dn")) return sf-sf_e;
  else return sf;
}


Double_t _getIPoverID(Double_t pt, Double_t eta, Int_t pdgId, TString sys)
{
  eta = TMath::Abs(eta);
  Double_t sf   = 0.;
  Double_t sf_e = 0.;
  if (TMath::Abs(pdgId) == 13){
    if (pt > 120) pt = 119.9;
    sf    = hMuonIP1->GetBinContent( hMuonIP1->FindBin(pt,eta));
    sf   *= hMuonIP2->GetBinContent( hMuonIP2->FindBin(pt,eta));
    //sf_e = (sys.Contains("Mu")) ? hMuonIP->GetBinError  ( hMuonIP->FindBin(pt,eta)) : 0.;
    // error hardcoded to 3 later on
    // sf_e = (sys.Contains("Mu")) ? hMuonID->GetBinError  ( hMuonID->FindBin(pt,eta)) : 0.;
    sf_e = 0;
  }
  else{
    if (pt > 200) pt = 199.9;
    sf   = hElecIP->GetBinContent( hElecIP->FindBin(pt,eta));
    sf_e = hElecIP->GetBinError  ( hElecIP->FindBin(pt,eta));
    sf_e = (sys.Contains("El")) ? hElecIP->GetBinError  ( hElecIP->FindBin(pt,eta)) : 0.;
  }
  if (sys.Contains("Up"))     return sf+sf_e;
  else if(sys.Contains("Dn")) return sf-sf_e;
  else return sf;
}


Double_t _getISOoverIP(Double_t pt, Double_t eta, Int_t pdgId, TString sys)
{
  eta = TMath::Abs(eta);
  Double_t sf   = 0.;
  Double_t sf_e = 0.;
  if (TMath::Abs(pdgId) == 13){
    if (pt > 120) pt = 119.9;
    sf   = hMuonISO->GetBinContent( hMuonISO->FindBin(pt,eta));
    // sf_e = (sys.Contains("Mu")) ? hMuonISO->GetBinError  ( hMuonISO->FindBin(pt,eta)) : 0.;
    // error hardcoded to 3 later on
    // sf_e = (sys.Contains("Mu")) ? hMuonID->GetBinError  ( hMuonID->FindBin(pt,eta)) : 0.;
    sf_e = 0;
  }
  else{
    if (pt > 200) pt = 199.9;
    sf   = hElecISO->GetBinContent( hElecISO->FindBin(pt,eta));
    sf_e = hElecISO->GetBinError  ( hElecISO->FindBin(pt,eta));
    sf_e = (sys.Contains("El")) ? hElecISO->GetBinError  ( hElecISO->FindBin(pt,eta)) : 0.;
  }
  if (sys.Contains("Up"))     return sf+sf_e;
  else if(sys.Contains("Dn")) return sf-sf_e;
  else return sf;
}


Double_t _getTracking(Double_t eta, Int_t pdgId)
{
  if (TMath::Abs(pdgId) == 13){
    return 1.;
    //  return hMuonTrk->GetBinContent( hMuonTrk->FindBin(20.,eta));
  }
  else{
    return hElecTrk->GetBinContent( hElecTrk->FindBin(20.,eta));
  }

}

Double_t LepSF(Double_t pt, Double_t eta, Int_t pdgId, TString sys="")
{
  // sys can be 
  // "": nominal
  // "ElUp": obvious
  // "ElDn": obvious
  // "MuUp": obvious
  // "MuDn": obvious

  Double_t sf = _getIDoverRECO(pt,eta,pdgId,sys)*_getIPoverID(pt,eta,pdgId,sys)*_getISOoverIP(pt,eta,pdgId,sys);
  if (TMath::Abs(pdgId) == 13){
    if (sys.Contains("MuUp")) return sf*1.03;
    else if (sys.Contains("MuDn")) return sf*0.97;
    else return sf;
  }
  return sf;

}


// Hasta aqui era fullsim/data
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////


// Double_t _getIDoverRECO(Double_t pt, Double_t eta, Int_t pdgId, TString sys)
// {
//   Double_t sf   = 0.;
//   Double_t sf_e = 0.;
//   if (TMath::Abs(pdgId) == 13){
//     sf   = hMuonID->GetBinContent( hMuonID->FindBin(pt,eta));
//     sf_e = (sys.Contains("Mu")) ? hMuonID->GetBinError  ( hMuonID->FindBin(pt,eta)) : 0.;
//   }
//   else{
//     sf   = hElecID->GetBinContent( hElecID->FindBin(pt,eta));
//     sf_e = hElecID->GetBinError  ( hElecID->FindBin(pt,eta));
//     sf_e = (sys.Contains("El")) ? hMuonID->GetBinError  ( hMuonID->FindBin(pt,eta)) : 0.;
//   }
//   if (sys.Contains("Up"))     return sf+sf_e;
//   else if(sys.Contains("Dn")) return sf-sf_e;
//   else return sf;
// }


// Double_t _getIPoverID(Double_t pt, Double_t eta, Int_t pdgId, TString sys)
// {
//   Double_t sf   = 0.;
//   Double_t sf_e = 0.;
//   if (TMath::Abs(pdgId) == 13){
//     sf   = hMuonIP->GetBinContent( hMuonIP->FindBin(pt,eta));
//     sf_e = (sys.Contains("Mu")) ? hMuonIP->GetBinError  ( hMuonIP->FindBin(pt,eta)) : 0.;
//   }
//   else{
//     sf   = hElecIP->GetBinContent( hElecIP->FindBin(pt,eta));
//     sf_e = hElecIP->GetBinError  ( hElecIP->FindBin(pt,eta));
//     sf_e = (sys.Contains("El")) ? hMuonIP->GetBinError  ( hMuonIP->FindBin(pt,eta)) : 0.;
//   }
//   if (sys.Contains("Up"))     return sf+sf_e;
//   else if(sys.Contains("Dn")) return sf-sf_e;
//   else return sf;
// }


// Double_t _getISOoverIP(Double_t pt, Double_t eta, Int_t pdgId, TString sys)
// {
//   Double_t sf   = 0.;
//   Double_t sf_e = 0.;
//   if (TMath::Abs(pdgId) == 13){
//     sf   = hMuonISO->GetBinContent( hMuonISO->FindBin(pt,eta));
//     sf_e = (sys.Contains("Mu")) ? hMuonISO->GetBinError  ( hMuonISO->FindBin(pt,eta)) : 0.;
//   }
//   else{
//     sf   = hElecISO->GetBinContent( hElecISO->FindBin(pt,eta));
//     sf_e = hElecISO->GetBinError  ( hElecISO->FindBin(pt,eta));
//     sf_e = (sys.Contains("El")) ? hMuonISO->GetBinError  ( hMuonISO->FindBin(pt,eta)) : 0.;
//   }
//   if (sys.Contains("Up"))     return sf+sf_e;
//   else if(sys.Contains("Dn")) return sf-sf_e;
//   else return sf;
// }



Double_t FastSIM(Double_t pt, Double_t eta, Int_t pdgId, TString sys)
{
  // sys can be 
  // "": nominal
  // "ElUp": obviously
  // "ElDn": obviously
  // "MuUp": obviously
  // "MuDn": obviously
  cout << "Check histo limits" << endl;
  return -1;
}