#include <TH1F.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <iostream>
#include <TStyle.h>
#include <TF1.h>
#include <TLegend.h>

#include "TDRStyle.h"

std::string ftoa0(double i) 
{
  char res[10];
  sprintf(res, "%2.0f", i);
  std::string ret(res);
  return ret;
}

std::string ftoa2(double i) 
{
  char res[10];
  sprintf(res, "%2.2f", i);
  std::string ret(res);
  return ret;
}

std::string ftoa3(double i) 
{
  char res[10];
  sprintf(res, "%2.3f", i);
  std::string ret(res);
  return ret;
}

Double_t crystalBall(Double_t *x, Double_t *par) // 5 Parameters
{
  // Double shouldered 
  Double_t std=(x[0]-par[0])/par[1];
  Double_t A=pow(par[3]/par[2], par[3])*exp(-0.5*pow(par[2], 2));
  Double_t B=par[3]/par[2]-par[2];
  Double_t result=0;
  if (std<=par[2]) // Gaussian Region
  {
    result=exp(-0.5*pow(std, 2));
  }
  else if (std>par[2]) // Power Law Region
  {
    result=A/pow(B+std, par[3]);
  }
  result=result*par[4];
  
  return result;
}

Double_t GaussExp(Double_t *x, Double_t *par)
{
  // Single shouldered
  Double_t std=(x[0]-par[0])/par[1];
  Double_t result=0;
  
  if (par[2]<std) result=exp(par[2]*par[2]/2.-par[2]*std); // high
  else if (std<=par[2]) result=exp(-0.5*pow(std, 2));
  
  result=result*par[3];
  
  return result;
}

void ThrowAndCatch()
{
  TF1 *f_parent=new TF1("f_parent", GaussExp, 0, 30, 4);
  f_parent->SetParameter(0, 5);
  f_parent->SetParameter(1, 1);
  f_parent->SetParameter(2, 0.5);
  f_parent->SetParameter(3, 1);
  
  TH1F *h_distribution=new TH1F("h_distribution", "; x; N", 30, 0, 30); h_distribution->Sumw2();
  h_distribution->GetYaxis()->SetTitleOffset(1.4);
  h_distribution->SetLineColor(kBlack);
  h_distribution->SetMarkerStyle(20);
  double nThrows=1e3;
  for (unsigned int i=0; i<nThrows; ++i)
     h_distribution->Fill(f_parent->GetRandom());
     
  TF1 *f_crystalBall=new TF1("f_crystalBall", crystalBall, 0, 30, 5);
  f_crystalBall->SetParLimits(0, 3.1, 7.1);   // mean
  f_crystalBall->SetParLimits(1, 0.1, 5.1);    // stddev
  f_crystalBall->SetParLimits(2, 0.1, 5.1);  // switch
  f_crystalBall->SetParLimits(3, 0.1, 200.1);  // exponent
  f_crystalBall->SetParLimits(4, nThrows*0.01, nThrows*2.);
  f_crystalBall->SetLineColor(kBlue);
  h_distribution->Fit(f_crystalBall, "LMEVR");
  double cb_mean=f_crystalBall->GetParameter(0);   double cb_mean_error=f_crystalBall->GetParError(0);
  double cb_std=f_crystalBall->GetParameter(1);    double cb_std_error=f_crystalBall->GetParError(1);
  double cb_alpha1=f_crystalBall->GetParameter(2); double cb_alpha1_error=f_crystalBall->GetParError(2);
  double cb_n1=f_crystalBall->GetParameter(3);     double cb_n1_error=f_crystalBall->GetParError(3);
  double cb_chi2=f_crystalBall->GetChisquare();
  double cb_NDF=f_crystalBall->GetNDF();
  std::cout<<"Crystal Ball fit_chi2/NDF = "<<cb_chi2/cb_NDF<<std::endl;
  
  TF1 *f_gaussExp=new TF1("f_gaussExp", GaussExp, 0, 30, 4);
  f_gaussExp->SetParLimits(0, 3.1, 7.1);
  f_gaussExp->SetParLimits(1, 0.1, 5.1);
  f_gaussExp->SetParLimits(2, 0.1, 5.1);
  f_gaussExp->SetParLimits(3, nThrows*0.01, nThrows*2.);
  f_gaussExp->SetLineColor(kRed);
  h_distribution->Fit(f_gaussExp, "LMER+");
  double fit_mean=f_gaussExp->GetParameter(0);    double fit_mean_error=f_gaussExp->GetParError(0);
  double fit_std=f_gaussExp->GetParameter(1);     double fit_std_error=f_gaussExp->GetParError(1);
  double fit_kH=f_gaussExp->GetParameter(2);      double fit_kH_error=f_gaussExp->GetParError(2);
  double fit_chi2=f_gaussExp->GetChisquare();
  double fit_NDF=f_gaussExp->GetNDF();
  std::cout<<"GaussExp fit_chi2/NDF = "<<fit_chi2/fit_NDF<<std::endl;
  
  TLegend *leg=new TLegend(0.89, 0.89, 0.60, 0.65);
  leg->SetFillStyle(1);
  leg->SetFillColor(kWhite); 
  leg->SetLineColor(kWhite);
  leg->AddEntry(f_gaussExp, "GaussExp Fit", "L");
  leg->AddEntry((TObject*)0, ("#bar{x} = "+ftoa2(fit_mean)+" #pm "+ftoa2(fit_mean_error)+" GeV").c_str(), "");
  leg->AddEntry((TObject*)0, ("#sigma = "+ftoa2(fit_std)+" #pm "+ftoa2(fit_std_error)+" GeV").c_str(), "");
  leg->AddEntry((TObject*)0, ("k = "+ftoa2(fit_kH)+" #pm "+ftoa2(fit_kH_error)).c_str(), "");
  leg->AddEntry((TObject*)0, ("#chi^{2}/NDF = "+ftoa2(fit_chi2/fit_NDF)).c_str(), "");
  
  TLegend *leg2=new TLegend(0.26, 0.89, 0.55, 0.65);
  leg2->SetFillStyle(1); 
  leg2->SetFillColor(kWhite); 
  leg2->SetLineColor(kWhite);
  leg2->AddEntry(f_crystalBall, "Crystal Ball Fit", "L");
  leg2->AddEntry((TObject*)0, ("#bar{x} = "+ftoa2(cb_mean)+" #pm "+ftoa2(cb_mean_error)+" GeV").c_str(), "");
  leg2->AddEntry((TObject*)0, ("#sigma = "+ftoa2(cb_std)+" #pm "+ftoa2(cb_std_error)+" GeV").c_str(), "");
  leg2->AddEntry((TObject*)0, ("#alpha = "+ftoa2(cb_alpha1)+" #pm "+ftoa2(cb_alpha1_error)).c_str(), "");
  leg2->AddEntry((TObject*)0, ("n = "+ftoa0(cb_n1)+" #pm "+ftoa0(cb_n1_error)).c_str(), "");
  leg2->AddEntry((TObject*)0, ("#chi^{2}/NDF = "+ftoa2(cb_chi2/cb_NDF)).c_str(), "");
  
  TLegend *leg3=new TLegend(0.50, 0.60, 0.89, 0.45);
  leg3->SetFillStyle(1); 
  // leg3->SetFillColor(kWhite); 
  leg3->SetLineColor(kWhite);
  leg3->AddEntry((TObject*)0, ("            (a)"), "");
  leg3->AddEntry((TObject*)0, ((ftoa0(nThrows)+" data points").c_str()), "");
  
  gROOT->SetStyle("Plain");
  TStyle *tdrStyle=setTDRStyle();
  tdrStyle->cd();
  tdrStyle->SetOptStat(0);
  
  TCanvas *c_parent=new TCanvas("c_parent", "c_parent", 700, 700);
  c_parent->SetLogy();
  f_parent->Draw();
  leg->Draw();
  leg2->Draw();
  c_parent->SaveAs("c_parent.png");
  
  TCanvas *c_distribution=new TCanvas("c_distribution", "c_distribution", 700, 700);
  c_distribution->SetLogy();
  h_distribution->SetMinimum(0.1);
  h_distribution->SetMaximum(h_distribution->GetMaximum()*5e2);
  h_distribution->Draw();
  leg->Draw();
  leg2->Draw();
  leg3->Draw();
  c_distribution->SaveAs("c_distribution.png");
  c_distribution->SaveAs("c_distribution.pdf");
}
  
  
