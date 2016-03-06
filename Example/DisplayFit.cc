#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLine.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TROOT.h"

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

Double_t crystalBall(Double_t *x, Double_t *par) // 7 Parameters
{
  // Double shouldered 
  Double_t std=(x[0]-par[0])/par[1];
  Double_t A=pow(par[3]/par[2], par[3])*exp(-0.5*pow(par[2], 2));
  Double_t B=par[3]/par[2]-par[2];
  Double_t C=pow(par[5]/fabs(par[4]), par[5])*exp(-0.5*pow(par[4], 2));
  Double_t D=par[5]/fabs(par[4])-fabs(par[4]);
  Double_t result=0;
  if (std>=par[4] && std<=par[2]) // Gaussian Region
  {
    result=exp(-0.5*pow(std, 2));
  }
  else if (std>par[2]) // Power Law Region
  {
    result=A/pow(B+std, par[3]);
  }
  else if (std<par[4]) // Power Law Region
  {
    result=C/pow(D-std, par[5]);
  }
  result=result*par[6];
  
  // std::cout<<"par[1] = "<<par[1]<<std::endl;
  
  return result;
}

Double_t ExpGaussExp(Double_t *x, Double_t *par)
{
  // Double shouldered
  Double_t std=(x[0]-par[0])/par[1];
  Double_t result=0;
  
  if (par[2]<std) result=exp(par[2]*par[2]/2.-par[2]*std); // high
  else if (-par[3]<std && std<=par[2]) result=exp(-0.5*pow(std, 2));
  else if (std<=-par[3]) result=exp(par[3]*par[3]/2.+par[3]*std); // low
  
  result=result*par[4];
  
  return result;
}

void DisplayFit()
{
  // TFile *file=new TFile("DelphesFiles/Histograms/Output_LowPtSUSY_Tree_HZZ.root");
  TFile *file=new TFile("DelphesFiles/Histograms/Output_LowPtSUSY_Tree_HGG.root");
  
  // TH1F *h=(TH1F*)file->Get("h_InvariantMass_HZZ");
  TH1F *h=(TH1F*)file->Get("h_InvariantMass_PhPh");
  
  double binSize=0.2; //GeV
  double range_min=90; // GeV
  double range_max=160; // GeV
  
  h->Rebin(binSize*10.);
  h->GetXaxis()->SetRangeUser(range_min, range_max);
  h->SetLineColor(kBlack);
  h->SetMarkerStyle(20);
  
  // Fits
  TF1 *f_crystalBall=new TF1("f_crystalBall", crystalBall, range_min+10, range_max-10, 7);
  f_crystalBall->SetParLimits(0, 115.1, 135.1);   // mean
  f_crystalBall->SetParLimits(1, 1.1, 10.1);    // stddev
  f_crystalBall->SetParLimits(2, 0.1, 10.1);  // switch
  f_crystalBall->SetParLimits(3, 0.1, 200.1);  // exponent
  f_crystalBall->SetParLimits(4, -10.1, -0.1);  // switch
  f_crystalBall->SetParLimits(5, 0.1, 200.1);  // exponent
  f_crystalBall->SetParLimits(6, 800.1, 2e3); // normalization
  f_crystalBall->SetLineColor(kBlue);
  h->Fit(f_crystalBall, "LMER");
  double cb_mean=f_crystalBall->GetParameter(0);    double cb_mean_error=f_crystalBall->GetParError(0);
  double cb_std=f_crystalBall->GetParameter(1);     double cb_std_error=f_crystalBall->GetParError(1);
  double cb_alpha1=f_crystalBall->GetParameter(2);  double cb_alpha1_error=f_crystalBall->GetParError(2);
  double cb_n1=f_crystalBall->GetParameter(3);      double cb_n1_error=f_crystalBall->GetParError(3);
  double cb_alpha2=f_crystalBall->GetParameter(4);  double cb_alpha2_error=f_crystalBall->GetParError(4);
  double cb_n2=f_crystalBall->GetParameter(5);      double cb_n2_error=f_crystalBall->GetParError(5);
  double cb_chi2=f_crystalBall->GetChisquare();
  double cb_NDF=f_crystalBall->GetNDF();
  
  TF1 *f_expGaussExp=new TF1("f_expGaussExp", ExpGaussExp, range_min+10, range_max-10, 5);
  f_expGaussExp->SetParLimits(0, 115.1, 135.1); // mean
  f_expGaussExp->SetParLimits(1, 1.1, 10.1); // stddev
  f_expGaussExp->SetParLimits(2, 0.1, 10.1);  // high exponent
  f_expGaussExp->SetParLimits(3, 0.1, 10.1);  // low exponent
  f_expGaussExp->SetParLimits(4, 800.1, 2e3); // normalization
  f_expGaussExp->SetLineColor(kRed);
  f_expGaussExp->SetLineWidth(2);
  h->Fit(f_expGaussExp, "LMER+");
  double fit_mean=f_expGaussExp->GetParameter(0);   double fit_mean_error=f_expGaussExp->GetParError(0);
  double fit_std=f_expGaussExp->GetParameter(1);    double fit_std_error=f_expGaussExp->GetParError(1);
  double fit_kH=f_expGaussExp->GetParameter(2);     double fit_kH_error=f_expGaussExp->GetParError(2);
  double fit_kL=f_expGaussExp->GetParameter(3);     double fit_kL_error=f_expGaussExp->GetParError(3);
  double fit_chi2=f_expGaussExp->GetChisquare();
  double fit_NDF=f_expGaussExp->GetNDF();
  
  h->SetTitle("; m_{H #rightarrow #gamma#gamma} (GeV); Events / 200 MeV");
  h->GetYaxis()->SetTitleOffset(1.7);
  
  TLine *l_mean=new TLine(fit_mean, 0, fit_mean, h->GetMaximum());
  TLine *l_lo=new TLine(fit_mean-fit_std*fit_kL, 0, fit_mean-fit_std*fit_kL, h->GetMaximum());
  TLine *l_hi=new TLine(fit_mean+fit_std*fit_kH, 0, fit_mean+fit_std*fit_kH, h->GetMaximum());
  l_lo->SetLineStyle(2);
  l_hi->SetLineStyle(2);
  
  TLegend *leg=new TLegend(0.91, 0.89, 0.65, 0.55);
  leg->SetFillStyle(1); 
  leg->SetFillColor(kWhite); 
  leg->SetLineColor(kWhite);
  leg->AddEntry(h, "Signal Sample", "EP");
  leg->AddEntry(f_expGaussExp, "ExpGaussExp Fit", "L");
  leg->AddEntry((TObject*)0, ("#bar{x} = "+ftoa2(fit_mean)+" #pm "+ftoa2(fit_mean_error)+" GeV").c_str(), "");
  leg->AddEntry((TObject*)0, ("#sigma = "+ftoa2(fit_std)+" #pm "+ftoa2(fit_std_error)+" GeV").c_str(), "");
  leg->AddEntry((TObject*)0, ("k_{L} = "+ftoa2(fit_kL)+" #pm "+ftoa2(fit_kL_error)).c_str(), "");
  leg->AddEntry((TObject*)0, ("k_{H} = "+ftoa2(fit_kH)+" #pm "+ftoa2(fit_kH_error)).c_str(), "");
  leg->AddEntry((TObject*)0, ("#chi^{2}/NDF = "+ftoa2(fit_chi2/fit_NDF)).c_str(), "");
  
  TLegend *leg2=new TLegend(0.17, 0.89, 0.42, 0.55);
  leg2->SetFillStyle(1); 
  leg2->SetFillColor(kWhite); 
  leg2->SetLineColor(kWhite);
  leg2->AddEntry(f_crystalBall, "Crystal Ball Fit", "L");
  leg2->AddEntry((TObject*)0, ("#bar{x} = "+ftoa2(cb_mean)+" #pm "+ftoa2(cb_mean_error)+" GeV").c_str(), "");
  leg2->AddEntry((TObject*)0, ("#sigma = "+ftoa2(cb_std)+" #pm "+ftoa2(cb_std_error)+" GeV").c_str(), "");
  leg2->AddEntry((TObject*)0, ("#alpha_{L} = "+ftoa2(cb_alpha1)+" #pm "+ftoa2(cb_alpha1_error)).c_str(), "");
  leg2->AddEntry((TObject*)0, ("n_{L} = "+ftoa2(cb_n1)+" #pm "+ftoa0(cb_n1_error)).c_str(), "");
  leg2->AddEntry((TObject*)0, ("#alpha_{H} = "+ftoa2(cb_alpha2)+" #pm "+ftoa2(cb_alpha2_error)).c_str(), "");
  leg2->AddEntry((TObject*)0, ("n_{H} = "+ftoa2(cb_n2)+" #pm "+ftoa0(cb_n2_error)).c_str(), "");
  leg2->AddEntry((TObject*)0, ("#chi^{2}/NDF = "+ftoa2(cb_chi2/cb_NDF)).c_str(), "");
  
  gROOT->SetStyle("Plain");
  TStyle *tdrStyle=setTDRStyle();
  tdrStyle->cd();
  tdrStyle->SetOptStat(0);
  
  TPave *box1=new TPave(0.18, 0.72, 0.4, 0.68, 4, "NDC");
  box1->SetLineColor(kRed);
  box1->SetLineWidth(3);
  box1->SetFillStyle(0);
  
  TPave *box2=new TPave(0.18, 0.635, 0.4, 0.595, 4, "NDC");
  box2->SetLineColor(kRed);
  box2->SetLineWidth(3);
  box2->SetFillStyle(0);
  
  TCanvas *c=new TCanvas("c", "c", 700, 700);
  h->SetMaximum(h->GetMaximum()*1.2);
  h->Draw("Ep9");
  f_expGaussExp->Draw("same");
  l_mean->Draw("same");
  l_lo->Draw("same");
  l_hi->Draw("same");
  leg->Draw();
  leg2->Draw();
  box1->Draw();
  box2->Draw();
  c->SaveAs("c.png");
  c->SaveAs("c.pdf");
  
  TCanvas *c_log=new TCanvas("c_log", "c_log", 700, 700);
  c_log->SetLogy();
  h->SetMinimum(0.1);
  h->Draw("Ep9");
  f_expGaussExp->Draw("same");
  l_mean->Draw("same");
  l_lo->Draw("same");
  l_hi->Draw("same");
  leg->Draw();
  leg2->Draw();
  box1->Draw();
  box2->Draw();
  h->Draw("Ep9 same");
  f_expGaussExp->Draw("same");
  c_log->SaveAs("c_log.png");
  c_log->SaveAs("c_log.pdf");
}
