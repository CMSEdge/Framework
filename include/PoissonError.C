//#ifndef POISSONERROR_HH
//#define POISSONERROR_HH
//=============================================================================
// Name:         PoissonErrors.C
//
// Goal:         Compute Poisson errors on observed number of events 
//
// References:   This code implements a few of the examples discussed in 
//               detail in an article by Bob Cousins
//               - Why isn't every physicist a Bayesian?
//                 Robert D. Cousins, (UCLA) . UCLA-HEP-94-005, Sep 1994. 27pp. 
//                 Published in Am.J.Phys.63:398,1995.
//        
// Authors/Date: Aart Heijboer and Ivo van Vulpen  (May 2011)
//=============================================================================

#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "Math/ProbFuncMathCore.h"
#include "TLatex.h"

#include <iostream>
using namespace std;



//=======================================================================================================================
void AddText( Double_t txt_x = 0.50, Double_t txt_y = 0.50, const char * txt = "dummy", Double_t txt_size = 0.045, 
	      Double_t txt_angle = 0., const char * Alignment = "left", Int_t UseNormalizedSize = 1, Int_t txt_color = 1){
//=======================================================================================================================
// -----------------------------------------------------------------------------------------------
// Small private function to plot things on the screen (because I got irritated by Root standards)
// -----------------------------------------------------------------------------------------------
  Int_t txt_align = 12;

  if ( !strcmp(Alignment, "left"))   { txt_align = 12; }
  if ( !strcmp(Alignment, "right"))  { txt_align = 32; }
  if ( !strcmp(Alignment, "center")) { txt_align = 22; }

  TLatex* t1 = new TLatex( txt_x, txt_y, txt);
  if(UseNormalizedSize) {t1->SetNDC(kTRUE);} // <- use NDC coordinate
  t1->SetTextSize(txt_size);
  t1->SetTextAlign(txt_align);
  t1->SetTextAngle(txt_angle);
  t1->SetTextColor(txt_color);
  t1->Draw();

} // end AddText()




	      


//============================================
double  *PoissonError(int n_obs, int Ioption = 4){
//============================================
// ----------------------------------------------------------------------
// Input: o n_obs   (observed number of events)
//        o Ioption (type of error to be returned to user)
//          Ioption: 1 = Classical central
//                   2 = Likelihood ratio 
//	             3 = Bayesian  central (flat prior - left/right integral = 16% )
//                   4 = Bayesian  central (flat prior - integrate pdf (equaly prob) 68%)
//          For detailed discussion see:
//          Why isn't every physicist a Bayesian?
//          Robert D. Cousins, (UCLA) . UCLA-HEP-94-005, Sep 1994. 27pp. 
//          Published in Am.J.Phys.63:398,1995.
//
// Output: ErrorRange
// Extra:   Iplot:   1/2 for plots on each of the methods for n_obs=4 Hindeloopen
//
// Howto:  root> .L PoissonError.C++
//         root> PoissonError(n_obs,1)    // method 1
// ----------------------------------------------------------------------

  //-----------------
  //-- Standard stuff
  //-----------------
  gROOT->Clear();
  gROOT->Delete();

  //-- Define variables and prepare histograms
  double Lambda_min = 0.;
  double Lambda_max = n_obs + 6*sqrt( n_obs );
  int    Nsteps_lambda = int(1e4);	           
  TH1D *h_likelihood = new TH1D("h_likelihood", "Likelihood for lambda   ", Nsteps_lambda, Lambda_min,Lambda_max);
  TH1D *h_min2loglikelihood = new TH1D("h_min2loglikelihood", "Likelihood for lambda   ", Nsteps_lambda, Lambda_min,Lambda_max);
  TH1D *h_pdf_full   = new TH1D("h_pdf_full",   "Pdf for lambda          ", Nsteps_lambda, Lambda_min,Lambda_max);

  
  //define output
  double ErrorRange[3];
  
  //-- Define a few standard parameters
  double Integral_fraction_1sigma = ROOT::Math::gaussian_cdf(-1.,1.,0.); // roughly 15.87 %


  // ===================================================================================
  // [A] Prepare Likelihood: Loop over hypotheses of mean of the 'true' Poisson (lambda)
  // ===================================================================================
  double lambda = 0.;
  for(Int_t i_bin_lambda = 1; i_bin_lambda <= Nsteps_lambda; i_bin_lambda++){
    lambda = h_likelihood->GetBinCenter(i_bin_lambda);
    if(fabs(lambda)< 1e-9) lambda = 1e-9;
    
    //-- Compute Poisson probability == Likelihood
    Double_t PoissonProbability = TMath::Poisson(n_obs,lambda);
    Double_t LogLikelihood      = -2.*TMath::Log(PoissonProbability); 

    //-- Fill histogram for Likelihood and probability density function
    h_likelihood->Fill(lambda, LogLikelihood);
    h_min2loglikelihood->Fill(lambda, LogLikelihood);
    h_pdf_full->Fill(lambda, PoissonProbability*1);                  // Bayes prior is constant '(*1)'    
  } // end loop over lambda hypotheses

  //-- save characteristic values
  int    bin_central        = h_likelihood->GetMinimumBin();              // bin with smallest -2*Log(likelihood)
  double LogLikelihood_min  = h_likelihood->GetBinContent(bin_central);
  double Lambda_central     = h_likelihood->GetBinCenter(bin_central);


  // ========================
  // [B] Compute error region
  // ========================

  //-- Compute error region using various options
  double Lambda_1sigma_low = -1.;
  double Lambda_1sigma_up  = -1.;

  // -----------------------------------------------------
  // Option 1: Frequentist (Classical central)
  // -----------------------------------------------------
  if(Ioption == 1){
    int Nobs_max = n_obs+100;  

    for(int i_bin = 1; i_bin < h_pdf_full->GetNbinsX(); i_bin++){ // loop over lambda
      lambda = h_pdf_full->GetBinCenter(i_bin);    
      double Poisson_sum_low = 0.;
      double Poisson_sum_up  = 0.;

      // lower value
      for(int i_obs = n_obs; i_obs < Nobs_max; i_obs++){ // loop from n_obs to infinity
         Poisson_sum_low += TMath::Poisson(i_obs,lambda);
      } 
      if(Poisson_sum_low > Integral_fraction_1sigma && Lambda_1sigma_low < 0.) {
       	 Lambda_1sigma_low = lambda;
      }

      // upper value
      for(int i_obs = 0; i_obs <= n_obs; i_obs++){ // loop from 0 to n_obs
         Poisson_sum_up += TMath::Poisson(i_obs,lambda);
      }
      if(Poisson_sum_up < Integral_fraction_1sigma && Lambda_1sigma_up < 0.) {
         Lambda_1sigma_up = lambda;
       }        
    } // end loop over lambda  
  } // end Ioption == 1


  // -----------------------------------------------------------------------
  // Option 2: Likelihood ratio
  //           Find region with delta -2log(likelihood) w.r.t minimum <1.
  // -----------------------------------------------------------------------
  if(Ioption == 2){
    Lambda_1sigma_low = Lambda_central;
    Lambda_1sigma_up  = Lambda_central;

    int i_bin_lambda     = -1;
    double LogLikelihood = -1.;
    int KeepMoving = 0;
  
    //-- left boundary  
    i_bin_lambda = bin_central;
    KeepMoving = (i_bin_lambda == 1) ? 0 : 1;
    while ( KeepMoving ){
      i_bin_lambda--;
      Lambda_1sigma_low = h_likelihood->GetBinCenter(i_bin_lambda);
      LogLikelihood     = h_likelihood->GetBinContent(i_bin_lambda);
      // decide when to stop
      if(LogLikelihood > LogLikelihood_min+1. || i_bin_lambda==1){
        KeepMoving = 0;
      }
    }

    //-- right boundary  
    i_bin_lambda = bin_central;
    KeepMoving = (i_bin_lambda == h_likelihood->GetNbinsX()) ? 0 : 1;
    while ( KeepMoving ){
      i_bin_lambda++;
      Lambda_1sigma_up = h_likelihood->GetBinCenter(i_bin_lambda);
      LogLikelihood    = h_likelihood->GetBinContent(i_bin_lambda);
      // decide when to stop
      if(LogLikelihood > LogLikelihood_min+1. || i_bin_lambda==h_likelihood->GetNbinsX()){
        KeepMoving = 0;
      }
    }
    
  } // end Ioption == 2



  // --------------------------------------------------------------
  // Option 3: Bayes central: flat prior
  //           Find left-error (16%) and right-error (16%) from PDF
  // --------------------------------------------------------------
  if(Ioption == 3){
    double integral_full     = h_pdf_full->Integral();
    double integral_fraction = 0.;

    for(int i_bin = 1; i_bin < h_pdf_full->GetNbinsX(); i_bin++){ // start loop over bins
       lambda = h_pdf_full->GetBinCenter(i_bin);
       integral_fraction = h_pdf_full->Integral(1,i_bin)/integral_full;  // fraction of PDF < i_bin
       if(integral_fraction > Integral_fraction_1sigma && Lambda_1sigma_low < 0.) {
       	  Lambda_1sigma_low = lambda;
       }
       if(integral_fraction > (1.-Integral_fraction_1sigma) && Lambda_1sigma_up < 0.) {
       	  Lambda_1sigma_up = lambda;
       }
    }  // end loop over bins

  } // end Ioption == 3


  // --------------------------------------------------------------
  // Option 4: Bayes central: flat prior
  //           Equal probability
  // --------------------------------------------------------------
  if(Ioption == 4){
  
    int Nbins = h_pdf_full->GetNbinsX();
    double Integral_tot  = h_pdf_full->Integral(1,Nbins);

    //-- start values
    int    i_bin_mostlikely = h_pdf_full->GetMaximumBin();  
    double i_bin_left       = i_bin_mostlikely;
    double i_bin_right      = i_bin_mostlikely;
    double prob_left         = h_pdf_full->GetBinContent(i_bin_mostlikely);
    double prob_right        = h_pdf_full->GetBinContent(i_bin_mostlikely);
  
    int    KeepMoving  = 1;
    int    Direction   = ( i_bin_mostlikely != Nbins ) ? 1 : -1;  // start right unless you are in the right most bin
    double integral    = 0.;

    while( KeepMoving ){      

      //-- Move right
      if( Direction == 1 ){      
	      i_bin_right += 1;
 	      prob_right = h_pdf_full->GetBinContent(i_bin_right); 
      }
      
      //-- Move left
      if( Direction == -1 ){      
	      i_bin_left -= 1;
 	      prob_left = h_pdf_full->GetBinContent(i_bin_left); 
      }
      
      //-- Decide which way to go    
      Direction = (prob_right < prob_left ) ? -1 : 1;
      if( Direction ==  1 && i_bin_right == Nbins ){ Direction = -1; } 
      if( Direction == -1 && i_bin_left  ==  1    ){ Direction =  1; } 
      
      integral = h_pdf_full->Integral(i_bin_left, i_bin_right) / Integral_tot;

      //-- Check if you can stop
      if( integral > (1.-2.*Integral_fraction_1sigma) ){
     	KeepMoving = 0;
      }

    } // end KeepMoving

    Lambda_1sigma_low = h_pdf_full->GetBinCenter(i_bin_left);
    Lambda_1sigma_up  = h_pdf_full->GetBinCenter(i_bin_right);

  } // end Ioption == 4






  //===============
  //-- Plots
  //===============

  TCanvas *canvas1;
  canvas1 = new TCanvas("canvas1","Standard Canvas",600,400);	
  canvas1->SetLeftMargin(0.125);
  canvas1->SetBottomMargin(0.125); 

  //---------------------------------  
  //-- Plot option 1: classic central
  //---------------------------------
  if(Ioption == 1){  
     Int_t Nobs_max = n_obs + 50;
     canvas1->Divide(1,2);

     // [A] Poisson Error frequentist upper limit
     TH1D *h_1   = new TH1D("h_1",   "poisson probability (sec)      ", Nobs_max+1, -0.5,Nobs_max+0.5);
     TH1D *h_up  = new TH1D("h_up",  "poisson probability (n >= nobs)", Nobs_max+1, -0.5,Nobs_max+0.5);
     for(Int_t i_obs = 0; i_obs <= Nobs_max; i_obs++){
       Double_t PoissonProb = TMath::Poisson(i_obs, Lambda_1sigma_low);
       h_1->Fill(Double_t(i_obs), PoissonProb);
       if(i_obs >= n_obs){
          h_up->Fill(Double_t(i_obs), PoissonProb);
       }
     }
     canvas1->cd(1);
     h_1->SetAxisRange(-0.5,n_obs*2,"X");
     h_1->Draw();
     AddText( 0.900, 0.035, "Number of events",0.060,0.,"right");             // X-axis
     AddText( 0.050, 0.900, "Poisson probability",0.060,90.,"right");         // Y-axis
     AddText( 0.800, 0.750, Form("#lambda = %5.2f",Lambda_1sigma_low),0.080,0.,"right");      
     h_up->SetFillColor(5);   
     h_up->Draw("same");
     h_1->Draw("axis same");
 
     // [B] Poisson Error frequentist upper limit
     TH1D *h_2   = new TH1D("h_2",   "poisson probability (sec)      ", Nobs_max+1, -0.5,Nobs_max+0.5);
     TH1D *h_low = new TH1D("h_low", "poisson probability (n <= nobs)", Nobs_max+1, -0.5,Nobs_max+0.5);  
     for(Int_t i_obs = 0; i_obs <= Nobs_max; i_obs++){
       Double_t PoissonProb = TMath::Poisson(i_obs, Lambda_1sigma_up);
       h_2->Fill(Double_t(i_obs), PoissonProb);
       if(i_obs <= n_obs){
          h_low->Fill(Double_t(i_obs), PoissonProb);
       }
     }
     canvas1->cd(2);
     h_2->SetAxisRange(-0.5,n_obs*2,"X");
     h_2->Draw();
     AddText( 0.900, 0.035, "Number of events",0.060,0.,"right");             // X-axis
     AddText( 0.050, 0.900, "Poisson probability",0.060,90.,"right");         // Y-axis
     AddText( 0.300, 0.750, Form("#lambda = %5.2f", Lambda_1sigma_up),0.080,0.,"right");      
     h_low->SetFillColor(5);   
     h_low->Draw("same");
     h_2->Draw("axis same");

  } // end plot with option 1 (classic central)



  //----------------------------------
  //-- Plot option 2: likelihood ratio
  //----------------------------------
   if(Ioption == 2){
      h_min2loglikelihood->SetAxisRange(0.,Lambda_max,"X");
      h_min2loglikelihood->SetAxisRange(0.,40.,"Y"); 
      h_min2loglikelihood->Draw();         
      AddText( 0.900, 0.035, "#lambda",0.060,0.,"right");                // X-axis
      AddText( 0.050, 0.900, "-2*Log(likelihood)",0.060,90.,"right");    // Y-axis
      //canvas1->Print("./Plots/PoissonError_likelihood.gif");
   }   // end plot 2 (likelihood ratio)


  //---------------------------------------------
  //-- Plot option 3 and 4: Bayes with flat prior
  //---------------------------------------------
  TH1D *h_pdf_1sigma = (TH1D*)h_pdf_full->Clone("h_pdf_1sigma"); h_pdf_1sigma->Reset();
  if(Ioption == 3 || Ioption == 4){

      //-- prepare histogram for 1 sigma region
      lambda = 0.;
      double prob   = 0.;
      double integral_left  = 0.;
      double integral_tot   = 0.;
      double integral_right = 0.;
      for(int i_bin = 1; i_bin < h_pdf_full->GetNbinsX(); i_bin++){
        lambda = h_pdf_full->GetBinCenter(i_bin);      
        prob   = h_pdf_full->GetBinContent(i_bin);       
        if(lambda <= Lambda_1sigma_low || lambda >= Lambda_1sigma_up){
          h_pdf_1sigma->SetBinContent(i_bin,prob);
  	    }
        integral_tot += prob;
        if(lambda <= Lambda_1sigma_low){ integral_left  += prob; }
        if(lambda >= Lambda_1sigma_up) { integral_right += prob; }
      }
      integral_left  = integral_left  / integral_tot;
      integral_right = integral_right / integral_tot;
      
      //-- draw plot
      h_pdf_full->Draw();
      h_pdf_1sigma->SetFillColor(5);
      h_pdf_1sigma->Draw("same");
      h_pdf_full->Draw("axis same");
      AddText( 0.900, 0.035, "#lambda",0.060,0.,"right");        // X-axis
      AddText( 0.050, 0.900, "Probability",0.060,90.,"right");   // Y-axis
      AddText( 0.65, 0.800, Form("N(obs) = %d",n_obs),0.050,0.,"left");  
   }   // end plot 3 or 4 (Bayes)




  //----------------------------
  //-- print to screen for users
  //----------------------------
  printf("\n\n");
  printf(" Observed number of events: Nobs = %d\n",n_obs);
  printf("\n");
//if(Ioption == 1){ printf(" Error Range using method %d: Classical central     \n", Ioption);}
//  if(Ioption == 2){ printf(" Error Range using method %d: Method: Likelihood ratio      \n", Ioption);}
//  if(Ioption == 3){ printf(" Error Range using method %d: Method: Bayesian (central)    \n", Ioption);}
//  if(Ioption == 4){ printf(" Error Range using method %d: Method: Bayesian (equal prob) \n", Ioption);}
  printf("  --> Confidence level interval = (%4.2f,%4.2f)\n",Lambda_1sigma_low,Lambda_1sigma_up);
  printf("\n\n");
  ErrorRange[0] = Lambda_1sigma_low;
  ErrorRange[1] = Lambda_central;
  ErrorRange[2] = Lambda_1sigma_up;
    	   
  return ErrorRange;
} // end PoissonError()








