#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <iostream>
#include <TStyle.h>

#include "TDRStyle.h"

double SR_lo=320.;
double SR_hi=2000.;

std::string ftoa2(double i) 
{
  char res[10];
  sprintf(res, "%2.2f", i);
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

void DisplayFitBackground()
{

  // TFile *data_8TeVData2012=new TFile("/Users/souvik/HbbHbb/Analysis/StandardSelection_bFirst_CMVA0p71_dR1p5_SR17p5/MMMM_nominal/a_KinFit/Histograms_BJetPlusX_Run2012BCD_Skim.root");
  TFile *data_8TeVData2012=new TFile("/Users/souvik/HbbHbb/Analysis/StandardSelection_bFirst_CMVA0p71_dR1p5_SR17p5_MMMMbar/MMMMbar/a_KinFit/Histograms_BJetPlusX_Run2012BCD_Skim.root");
  TH1F *h_mX_SR=(TH1F*)data_8TeVData2012->Get("h_mX_SR");
  
  gROOT->SetStyle("Plain");
  TStyle *tdrStyle=setTDRStyle();
  tdrStyle->cd();
  tdrStyle->SetOptStat(0);
  
  // h_mX_SR->Rebin(2);
  h_mX_SR->SetMarkerStyle(20);
  h_mX_SR->SetLineColor(kBlack);
  h_mX_SR->GetXaxis()->SetRangeUser(SR_lo, SR_hi);
  
  TF1 *f_crystalBall=new TF1("f_crystalBall", crystalBall, SR_lo, SR_hi, 5);
  f_crystalBall->SetParLimits(0, 450.1, 550.1);   // mean
  f_crystalBall->SetParLimits(1, 20.1, 70.1);    // stddev
  f_crystalBall->SetParLimits(2, 0.1, 10.1);  // switch
  f_crystalBall->SetParLimits(3, 0.1, 200.1);  // exponent
  f_crystalBall->SetLineColor(kGreen);
  h_mX_SR->Fit(f_crystalBall, "R");
  double cb_mean=f_crystalBall->GetParameter(0);
  double cb_std=f_crystalBall->GetParameter(1);
  double cb_alpha1=f_crystalBall->GetParameter(2);
  double cb_n1=f_crystalBall->GetParameter(3);
  double cb_chi2=f_crystalBall->GetChisquare();
  double cb_NDF=f_crystalBall->GetNDF();
  
  TF1 *f_gaussExp=new TF1("f_gaussExp", GaussExp, SR_lo, SR_hi, 4);
  f_gaussExp->SetParLimits(0, 450.1, 520.1);
  f_gaussExp->SetParLimits(1, 20.1, 70.1);
  f_gaussExp->SetParLimits(2, 0.1, 10.1);
  f_gaussExp->SetLineColor(kRed);
  h_mX_SR->Fit(f_gaussExp, "R+");
  double fit_mean=f_gaussExp->GetParameter(0);
  double fit_std=f_gaussExp->GetParameter(1);
  double fit_kH=f_gaussExp->GetParameter(2);
  double fit_chi2=f_gaussExp->GetChisquare();
  double fit_NDF=f_gaussExp->GetNDF();
  
  h_mX_SR->SetTitle("; m_{HH#rightarrow(b#bar{b})(b#bar{b})} (GeV); Events / 10 GeV");
  h_mX_SR->GetYaxis()->SetTitleOffset(1.7);
  
  TLine *l_mean=new TLine(fit_mean, 0, fit_mean, h_mX_SR->GetMaximum());
  TLine *l_hi=new TLine(fit_mean+fit_std*fit_kH, 0, fit_mean+fit_std*fit_kH, h_mX_SR->GetMaximum());
  l_hi->SetLineStyle(2);
  
  TLegend *leg=new TLegend(0.89, 0.89, 0.60, 0.60);
  leg->SetFillStyle(1); 
  leg->SetFillColor(kWhite); 
  leg->SetLineColor(kWhite);
  leg->AddEntry(h_mX_SR, "CMS pp collision data", "EP");
  leg->AddEntry((TObject*)0, "17.9 fb^{-1} (8 TeV)", "");
  leg->AddEntry(f_gaussExp, "GaussExp Fit", "L");
  leg->AddEntry((TObject*)0, ("#bar{x} = "+ftoa2(fit_mean)+" GeV").c_str(), "");
  leg->AddEntry((TObject*)0, ("#sigma = "+ftoa2(fit_std)+" GeV").c_str(), "");
  leg->AddEntry((TObject*)0, ("k = "+ftoa2(fit_kH)).c_str(), "");
  leg->AddEntry((TObject*)0, ("#chi^{2}/NDF = "+ftoa2(fit_chi2/fit_NDF)).c_str(), "");
  
  TLegend *leg2=new TLegend(0.18, 0.89, 0.4, 0.65);
  leg2->SetFillStyle(1); 
  leg2->SetFillColor(kWhite); 
  leg2->SetLineColor(kWhite);
  leg2->AddEntry(f_crystalBall, "Crystal Ball Fit", "L");
  leg2->AddEntry((TObject*)0, ("#bar{x} = "+ftoa2(cb_mean)+" GeV").c_str(), "");
  leg2->AddEntry((TObject*)0, ("#sigma = "+ftoa2(cb_std)+" GeV").c_str(), "");
  leg2->AddEntry((TObject*)0, ("#alpha = "+ftoa2(cb_alpha1)).c_str(), "");
  leg2->AddEntry((TObject*)0, ("n = "+ftoa2(cb_n1)).c_str(), "");
  leg2->AddEntry((TObject*)0, ("#chi^{2}/NDF = "+ftoa2(cb_chi2/cb_NDF)).c_str(), "");
  
  TCanvas *c_mX_SR=new TCanvas("c_mX_SR", "c_mX_SR", 700, 700);
  h_mX_SR->SetMaximum(h_mX_SR->GetMaximum()*1.5);
  h_mX_SR->Draw("Ep9");
  l_mean->Draw("same");
  l_hi->Draw("same");
  leg->Draw();
  leg2->Draw();
  c_mX_SR->SaveAs("c_mX_SR.png");
  c_mX_SR->SaveAs("c_mX_SR.pdf");
  
  TCanvas *c_mX_SR_log=new TCanvas("c_mX_SR_log", "c_mX_SR_log", 700, 700);
  c_mX_SR_log->SetLogy();
  h_mX_SR->SetMinimum(0.1);
  h_mX_SR->SetMaximum(h_mX_SR->GetMaximum()*100.);
  h_mX_SR->Draw("Ep9");
  l_mean->Draw("same");
  l_hi->Draw("same");
  leg->Draw();
  leg2->Draw();
  c_mX_SR_log->SaveAs("c_mX_SR_log.png");
  c_mX_SR_log->SaveAs("c_mX_SR_log.pdf");
  
}
