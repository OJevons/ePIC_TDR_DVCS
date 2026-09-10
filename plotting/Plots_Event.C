using namespace std;

#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include "./ePIC_style.C"

// INITIAL SETTINGS
// Choose to print plots
bool kSAVE = true;
// PATH TO MAIN REPO DIRECTORY
TString sMainDir = "/home/ojj2x/eic/ePIC_TDR_DVCS/";
TString sFilePath = sMainDir + "rootfiles/";
TString sFigsPath = sMainDir + "figs/";

//---------------------------------------------------------------------
// Helper functions for printing strings
//---------------------------------------------------------------------

// 1. Decompose a double into scientific notation (mantissa, exponent)
void frexp10(double x, double &mantissa, int &exp10){
  if(x == 0.){
    mantissa = 0.;
    exp10 = 0;
    return;
  }
  
  exp10 = static_cast<int>(std::floor(std::log10(std::fabs(x))));
  mantissa = x / std::pow(10., exp10);
}

// 2. Print (a +/- b) such that both terms have the same exponent
TString combineScientific(double val, double err){
  double val_m, err_m;
  int val_e, err_e;
  
  // Decompose value and error
  frexp10(val, val_m, val_e);
  frexp10(err, err_m, err_e);
  
  // Use exponent of value as the common exponent
  int exp_common = val_e;
  
  // Rescale error to common exponent
  double scaled_val = val_m;
  double scaled_err = err_m * std::pow(10., err_e-exp_common);
  
  double rounded_val = std::round(scaled_val * 100.)/100.;
  double rounded_err = std::round(scaled_err * 100.)/100.;

  TString outstring = Form("(%.2f #pm %.2f)x10^{%i}",scaled_val,scaled_err,exp_common);
  
  return outstring;
}

//---------------------------------------------------------------------
// Calculation of FT from exponential (and error)
//---------------------------------------------------------------------
Double_t FTFunc(Double_t *x, Double_t *par){
  Float_t xx = x[0];
  Float_t A = par[0];
  Float_t B = par[1];
  Float_t tmin = par[2];
  Float_t tmax = par[3];

  Double_t factor = A / (B*B + xx*xx);
  Double_t expomin = TMath::Exp(-B*tmin);
  Double_t expomax = TMath::Exp(-B*tmax);
  
  Double_t trigmin = B*expomin*TMath::Cos(xx*tmin) - xx*expomin*TMath::Sin(xx*tmin);
  Double_t trigmax = xx*expomax*TMath::Sin(xx*tmax) - B*expomax*TMath::Cos(xx*tmax);

  Double_t f = factor*(trigmin+trigmax);

  return f;
}

Double_t FTErrFunc(Double_t *x, Double_t *par){
  Float_t b = x[0];
  Float_t A = par[0];
  Float_t dA = par[1];
  Float_t B = par[2];
  Float_t dB = par[3];
  Float_t tmin = par[4];
  Float_t tmax = par[5];

  Double_t expmin = TMath::Exp(-B*tmin);
  Double_t expmax = TMath::Exp(-B*tmax);
  Double_t den = 1/(B*B + b*b);
  Double_t den2 = den*den;

  Double_t dfda = den * ( B*expmin*TMath::Cos(b*tmin) - b*expmin*TMath::Sin(b*tmin) + b*expmax*TMath::Sin(b*tmax) - B*expmax*TMath::Cos(b*tmax) );

  Double_t dfdb_1 = (expmin*den2)*TMath::Cos(b*tmin)*( b*b*(1 - B*tmin) - B*B*(B*tmin + 1) );
  Double_t dfdb_2 = (expmax*den2)*TMath::Cos(b*tmax)*( b*b*(1 - B*tmax) - B*B*(B*tmax + 1) );
  Double_t dfdb_3 = (b*expmin*den2)*TMath::Sin(b*tmin)*( 2*B + tmin*(B*B + b*b) );
  Double_t dfdb_4 = (b*expmax*den2)*TMath::Sin(b*tmax)*( 2*B + tmax*(B*B + b*b) );
  Double_t dfdb = A * (dfdb_1 - dfdb_2 + dfdb_3 - dfdb_4);
  
  Double_t df2 = (dfda*dfda*dA*dA) + (dfdb*dfdb*dB*dB);
  
  return TMath::Sqrt(df2);
}

//---------------------------------------------------------------------
// Extract sliced RMS of 2D histo
//---------------------------------------------------------------------

TH1D* extractRMSSlice(TString outputHistoName, TH2D* twoDHisto){

	int num_bins  = twoDHisto->GetNbinsX();
	double xBinWidth = twoDHisto->GetXaxis()->GetBinWidth(1); 
	double xMin = twoDHisto->GetXaxis()->GetBinCenter(1) - xBinWidth*0.5;
	double xMax = twoDHisto->GetXaxis()->GetBinCenter(num_bins) + xBinWidth*0.5;

	TH1D * finalResoHisto = new TH1D(outputHistoName, outputHistoName, num_bins, xMin, xMax);

	TH1D* tmp;
	double rmsReso = 0.0;
	double rmsErr = 0.0;
	for(int bin = 1; bin < num_bins+1; bin++){
	  rmsReso = 0.0;
	  tmp = (TH1D*)twoDHisto->ProjectionY("NEIN", bin, bin);
	  
	  rmsReso = tmp->GetRMS();
	  rmsErr  = tmp->GetRMSError();
	  finalResoHisto->SetBinContent(bin, rmsReso);
	  finalResoHisto->SetBinError(bin, rmsErr);
	} 

	return finalResoHisto;
}

//---------------------------------------------------------------------
// Calculate histogram scaling factors to get to desired luminosity
// Give luminosity in full scientific (expect 1e15 for 1fb-1)
//---------------------------------------------------------------------
double calcScaling(TString energy, TString hel, float lumi){
  // Holding variables - No. of events generated, integrated cross-section
  // These vary by beam settings
  Double_t fXSint{0}, NEv{0};

  if(energy == "9x130"){
    NEv = 1e6;
    if(hel.Contains("T")){
      if(hel == "emhTm") fXSint = 6.91502943051105e-9;
      else if(hel == "emhTp") fXSint = 7.00501745732437e-9;
      else if(hel == "ephTm") fXSint = 6.92070657974331e-9;
      else if(hel == "ephTp") fXSint = 7.00302210771398e-9;
    }
    else if(hel.Contains("L")){
      if(hel == "emhTm") fXSint = 6.94645727655662e-9;
      else if(hel == "emhTp") fXSint = 6.97831632951577e-9;
      else if(hel == "ephTm") fXSint = 6.97649277129601e-9;
      else if(hel == "ephTp") fXSint = 6.95324310420876e-9;
    }
    else fXSint = 1;
  }
  else if(energy == "9x275"){
    NEv = 2.5e6;
    if(hel.Contains("T")){
      if(hel == "emhTm") fXSint = 7.65782690151261e-9;
      else if(hel == "emhTp") fXSint = 7.73125664480267e-9;
      else if(hel == "ephTm") fXSint = 7.65495841373398e-9;
      else if(hel == "ephTp") fXSint = 7.73434664504977e-9;
    }
    else if(hel.Contains("L")){
      if(hel == "emhTm") fXSint = 7.6845279388898e-9;
      else if(hel == "emhTp") fXSint = 7.70176410396224e-9;
      else if(hel == "ephTm") fXSint = 7.70861721413337e-9;
      else if(hel == "ephTp") fXSint = 7.68892983520446e-9;
    }
    else fXSint = 1;
  }

  Double_t genlumi = NEv/fXSint;
  Double_t scale = lumi/genlumi;

  return scale;
}

//---------------------------------------------------------------------
// MAIN
// 
// ePIC TDR plots - single energy
//---------------------------------------------------------------------
void Plots_Event(TString campaign = "26.07.1", TString energy = "9x130", TString hel = "emhTm"){
  // Print plot settings
  cout<<"---------------------------------"<<endl;
  cout<<"Processing plots - single energy"<<endl;
  cout<<"\tCampaign: "<<campaign<<endl;
  cout<<"\tEnergy: "<<energy<<endl;
  cout<<"\tBeam helicity: "<<hel<<endl;
  cout<<"---------------------------------"<<endl;

  // Set beam energies
  Float_t fEBeam{0}, fPBeam{0};
  if(energy == "9x130" || energy == "9x275"){
    fEBeam = 9.;
    if(energy == "9x130") fPBeam = 130.;
    else if(energy == "9x275") fPBeam = 275.;
  }
  else{
    cout<<"Invalid beam energy."<<endl;
    return;
  }

  // Load chosen input file
  TString sIn = sFilePath + "ePIC_DVCS_"+campaign+"_"+energy+"_"+hel+".root";
  TFile* fIn = TFile::Open(sIn);

  //--------------------------------------------------------------------
  // Load histograms from file
  //--------------------------------------------------------------------
  // 1) Inclusive, Q2 - distribution, efficiency, 2D response, resolution
  TH1D* h_Q2_MC  = (TH1D*)fIn->Get("q2_mc");
  TH1D* h_Q2_Acc = (TH1D*)fIn->Get("q2_acc");
  TH1D* h_Q2_Rec = (TH1D*)fIn->Get("q2_reco");
  TH2D* h_Q2_2D  = (TH2D*)fIn->Get("q2_2d");
  TH2D* h_Q2_Res = (TH2D*)fIn->Get("q2_pctres");
  
  // 2) Inclusive, xB - distribution, efficiency, 2D response, resolution
  TH1D* h_xB_MC  = (TH1D*)fIn->Get("xb_mc");
  TH1D* h_xB_Acc = (TH1D*)fIn->Get("xb_acc");
  TH1D* h_xB_Rec = (TH1D*)fIn->Get("xb_reco");
  TH2D* h_xB_2D  = (TH2D*)fIn->Get("xb_2d");
  TH2D* h_xB_Res = (TH2D*)fIn->Get("xb_pctres");
  
  // 3) Inclusive, y - distribution, efficiency, 2D response, resolution
  TH1D* h_y_MC  = (TH1D*)fIn->Get("y_mc");
  TH1D* h_y_Acc = (TH1D*)fIn->Get("y_acc");
  TH1D* h_y_Rec = (TH1D*)fIn->Get("y_reco");
  TH2D* h_y_2D  = (TH2D*)fIn->Get("y_2d");
  TH2D* h_y_Res = (TH2D*)fIn->Get("y_pctres");
  
  // 4) Exclusive, t - distribution, efficiency, resolution
  TH1D* h_t_MC    = (TH1D*)fIn->Get("t_truth");
  TH1D* h_t_B0Acc = (TH1D*)fIn->Get("t_b0acc");
  TH1D* h_t_RPAcc = (TH1D*)fIn->Get("t_rpacc");
  TH1D* h_t_LCAcc = (TH1D*)fIn->Get("t_lcacc");
  TH1D* h_t_B0Rec = (TH1D*)fIn->Get("t_b0reco");
  TH1D* h_t_RPRec = (TH1D*)fIn->Get("t_rpreco");
  TH1D* h_t_LCRec = (TH1D*)fIn->Get("t_lcreco");
  // Resolution - absolute
  TH2D* h_t_B0Res = (TH2D*)fIn->Get("tresb0_2d");
  TH2D* h_t_RPRes = (TH2D*)fIn->Get("tresrp_2d");
  TH2D* h_t_LCRes = (TH2D*)fIn->Get("treslc_2d");
  // Resolution - percentage
  TH2D* h_t_B0ResPct = (TH2D*)fIn->Get("tresb0pct_2d");
  TH2D* h_t_RPResPct = (TH2D*)fIn->Get("tresrppct_2d");
  TH2D* h_t_LCResPct = (TH2D*)fIn->Get("treslcpct_2d");

  // 5) Exclusive, pT_miss
  TH1D* h_pTmiss_MC = (TH1D*)fIn->Get("ptmiss3_mc");
  TH1D* h_pTmiss_RP = (TH1D*)fIn->Get("ptmiss3_rp");

  // 6) Full event - (E-pz)
  TH1D* h_EmPz_MC = (TH1D*)fIn->Get("empz3_mc");
  TH1D* h_EmPz_RP = (TH1D*)fIn->Get("empz3_rp");  

  // 7) 2D inclusive correlations - x/Q2, x/t
  TH2D* h_xBvQ2_Rec = (TH2D*)fIn->Get("2d_xvq2_rp");
  TH2D* h_xBvt_Rec = (TH2D*)fIn->Get("2d_xvt_rp");


  // CALCULATIONS
  // 1. Scaling factor to EIC lumi
  double EIClumi{1.};
  if(energy == "9x130") EIClumi = 1e15;
  else if(energy == "9x275") EIClumi = 2.5e15;
  else EIClumi = 2.5e15;
  Double_t scaleToEIC = calcScaling(energy, hel, EIClumi);

  // Scale histograms
  h_Q2_MC->Scale(scaleToEIC);
  h_Q2_Acc->Scale(scaleToEIC);
  h_Q2_Rec->Scale(scaleToEIC);
  //
  h_xB_MC->Scale(scaleToEIC);
  h_xB_Acc->Scale(scaleToEIC);
  h_xB_Rec->Scale(scaleToEIC);
  //
  h_y_MC->Scale(scaleToEIC);
  h_y_Acc->Scale(scaleToEIC);
  h_y_Rec->Scale(scaleToEIC);
  //
  h_t_MC->Scale(scaleToEIC);
  h_t_B0Acc->Scale(scaleToEIC);
  h_t_RPAcc->Scale(scaleToEIC);
  h_t_LCAcc->Scale(scaleToEIC);
  h_t_B0Rec->Scale(scaleToEIC);
  h_t_RPRec->Scale(scaleToEIC);
  h_t_LCRec->Scale(scaleToEIC);
  //
  h_EmPz_MC->Scale(scaleToEIC);
  h_EmPz_RP->Scale(scaleToEIC);
  //
  h_pTmiss_RP->Scale(scaleToEIC);
    
  // 2. Efficiency and corrected reco
  TH1D* h_Q2_Eff = (TH1D*)h_Q2_Acc->Clone("q2_eff");
  h_Q2_Eff->Divide(h_Q2_MC);
  TH1D* h_Q2_Corr = (TH1D*)h_Q2_Rec->Clone("q2_corr");
  h_Q2_Corr->Divide(h_Q2_Eff);
  //
  TH1D* h_xB_Eff = (TH1D*)h_xB_Acc->Clone("xB_eff");
  h_xB_Eff->Divide(h_xB_MC);
  TH1D* h_xB_Corr = (TH1D*)h_xB_Rec->Clone("xB_corr");
  h_xB_Corr->Divide(h_xB_Eff);
  //
  TH1D* h_y_Eff = (TH1D*)h_y_Acc->Clone("y_eff");
  h_y_Eff->Divide(h_y_MC);
  TH1D* h_y_Corr = (TH1D*)h_y_Rec->Clone("y_corr");
  h_y_Corr->Divide(h_y_Eff);
  //
  TH1D* h_t_B0Eff = (TH1D*)h_t_B0Acc->Clone("t_b0eff");
  h_t_B0Eff->Divide(h_t_MC);
  TH1D* h_t_B0Corr = (TH1D*)h_t_B0Rec->Clone("t_b0corr");
  h_t_B0Corr->Divide(h_t_B0Eff);
  TH1D* h_t_RPEff = (TH1D*)h_t_RPAcc->Clone("t_rpeff");
  h_t_RPEff->Divide(h_t_MC);
  TH1D* h_t_RPCorr = (TH1D*)h_t_RPRec->Clone("t_rpcorr");
  h_t_RPCorr->Divide(h_t_RPEff);
  TH1D* h_t_LCEff = (TH1D*)h_t_LCAcc->Clone("t_lceff");
  h_t_LCEff->Divide(h_t_MC);
  TH1D* h_t_LCCorr = (TH1D*)h_t_LCRec->Clone("t_lccorr");
  h_t_LCCorr->Divide(h_t_LCEff);

  //--------------------------------------------------------------------------------------
  // Draw histograms
  //--------------------------------------------------------------------------------------
  // Global style
  gROOT->ProcessLine("set_ePIC_style()");
  gStyle->SetCanvasPreferGL(kTRUE);

  // CANVAS: Q2 distribution
  TCanvas* cQ2 = new TCanvas("cq2","",1000,1000);
  TPad* upperq2 = new TPad("upperq2","",0.05,0.3,0.95,0.95);
  upperq2->SetBottomMargin(0);
  upperq2->Draw();
  TPad* lowerq2 = new TPad("lowerq2","",0.05,0.05,0.95,0.3);
  lowerq2->SetBottomMargin(0.22);
  lowerq2->SetTopMargin(0);
  lowerq2->Draw();
  // UPPER - Distribution
  upperq2->cd();
  gPad->SetLogy();
  // Markers and lines
  h_Q2_MC->SetLineColor(kBlack);
  h_Q2_MC->SetLineWidth(2);
  h_Q2_Rec->SetMarkerColor(kP6Blue);
  h_Q2_Rec->SetMarkerStyle(kOpenCircle);
  h_Q2_Rec->SetMarkerSize(1.5);
  h_Q2_Rec->SetLineColor(kP6Blue);
  h_Q2_Rec->SetLineWidth(2);
  h_Q2_Corr->SetMarkerColor(kP6Blue);
  h_Q2_Corr->SetMarkerStyle(kFullCircle);
  h_Q2_Corr->SetMarkerSize(1.5);
  h_Q2_Corr->SetLineColor(kP6Blue);
  h_Q2_Corr->SetLineWidth(2);
  // Axes
  h_Q2_MC->GetYaxis()->SetTitle("Counts / 0.02 GeV^{2}");
  h_Q2_MC->GetYaxis()->SetTitleSize(0.05);
  h_Q2_MC->GetYaxis()->SetLabelSize(0.05);
  // Drawing
  h_Q2_MC->Draw("hist");
  h_Q2_Rec->Draw("pe same");
  h_Q2_Corr->Draw("pe same");
  // Labels
  TLatex* tePICLabel = new TLatex(0.17, 0.86, 
				  Form("#splitline{#bf{ePIC} Performance %s, %s GeV}{ep #rightarrow e'p'#gamma, L_{proj} = %.1f fb^{-1}}", campaign.Data(), energy.Data(), EIClumi/1e15));
  tePICLabel->SetNDC();
  tePICLabel->SetTextSize(0.05);
  tePICLabel->Draw("same");
  // Legend
  TLegend* lQ2 = new TLegend(0.65, 0.6, 0.99, 0.81);
  lQ2->SetLineWidth(0);
  lQ2->SetFillStyle(0);
  lQ2->AddEntry(h_Q2_MC, "MC gen.", "l");
  lQ2->AddEntry(h_Q2_Rec, "Raw reco.", "lp");
  lQ2->AddEntry(h_Q2_Corr, "Corr. reco.", "lp");
  lQ2->Draw();
  
  // LOWER - Ratio/efficiency
  lowerq2->cd();
  // Markers and lines
  h_Q2_Eff->SetLineColor(kBlack);
  h_Q2_Eff->SetLineWidth(2);
  // Axes
  h_Q2_Eff->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
  h_Q2_Eff->GetXaxis()->SetTitleSize(0.11);
  h_Q2_Eff->GetXaxis()->SetTitleOffset(0.93);
  h_Q2_Eff->GetXaxis()->SetLabelSize(0.11);
  h_Q2_Eff->GetYaxis()->SetTitle("Acc/MC");
  h_Q2_Eff->GetYaxis()->SetTitleSize(0.13);
  h_Q2_Eff->GetYaxis()->SetTitleOffset(0.4);
  h_Q2_Eff->GetYaxis()->SetLabelSize(0.13);
  // Draw
  h_Q2_Eff->Draw();
  // Save figure
  if(kSAVE) cQ2->SaveAs(sFigsPath + "TDR_" + energy +"_Q2.png");
  cQ2->Close();

  // CANVAS: Q2 2D
  TCanvas* cQ2_2D = new TCanvas("cq2_2d","",1000,1000);
  gPad->SetLogz();
  gPad->SetRightMargin(0.12);
  h_Q2_2D->GetXaxis()->SetTitle("Q^{2}_{MC} [GeV^2]");
  h_Q2_2D->GetYaxis()->SetTitle("Q^{2}_{reco} [GeV^2]");
  h_Q2_2D->Draw("colz");
  // Save figure
  if(kSAVE) cQ2_2D->SaveAs(sFigsPath + "TDR_" + energy +"_Q2_2D.png");
  cQ2_2D->Close();

  // CANVAS: Q2 resolution
  TCanvas* cQ2_Res = new TCanvas("cq2_res","",1000,1000);
  gPad->SetLogz();
  gPad->SetRightMargin(0.12);
  h_Q2_Res->GetYaxis()->SetTitle("(Q^{2}_{reco}-Q^{2}_{MC})/Q^{2}_{MC}");
  h_Q2_Res->GetXaxis()->SetTitle("Q^{2}_{MC} [GeV^{2}]");
  h_Q2_Res->Draw("colz");
  // Save figure
  if(kSAVE) cQ2_Res->SaveAs(sFigsPath + "TDR_" + energy +"_Q2_Res.png");
  cQ2_Res->Close();

  // CANVAS: xB distribution
  TCanvas* cxB = new TCanvas("cxb","",1000,1000);
  TPad* upperxb = new TPad("upperxb","",0.05,0.3,0.95,0.95);
  upperxb->SetBottomMargin(0);
  upperxb->Draw();
  TPad* lowerxb = new TPad("lowerxb","",0.05,0.05,0.95,0.3);
  lowerxb->SetBottomMargin(0.22);
  lowerxb->SetTopMargin(0);
  lowerxb->Draw();
  // UPPER - Distribution
  upperxb->cd();
  gPad->SetLogy();
  gPad->SetLogx();
  // Markers and lines
  h_xB_MC->SetLineColor(kBlack);
  h_xB_MC->SetLineWidth(2);
  h_xB_Rec->SetMarkerColor(kP6Blue);
  h_xB_Rec->SetMarkerStyle(kOpenCircle);
  h_xB_Rec->SetMarkerSize(1.5);
  h_xB_Rec->SetLineColor(kP6Blue);
  h_xB_Rec->SetLineWidth(2);
  h_xB_Corr->SetMarkerColor(kP6Blue);
  h_xB_Corr->SetMarkerStyle(kFullCircle);
  h_xB_Corr->SetMarkerSize(1.5);
  h_xB_Corr->SetLineColor(kP6Blue);
  h_xB_Corr->SetLineWidth(2);
  // Axes
  h_xB_MC->GetYaxis()->SetTitle("Counts");
  h_xB_MC->GetYaxis()->SetTitleSize(0.05);
  h_xB_MC->GetYaxis()->SetLabelSize(0.05);
  // Drawing
  h_xB_MC->Draw("hist");
  h_xB_Rec->Draw("pe same");
  h_xB_Corr->Draw("pe same");
  // Labels
  TLatex* tePICLabel_R = new TLatex(0.39, 0.87, 
				    Form("#splitline{#bf{ePIC} Performance %s, %s GeV}{ep #rightarrow e'p'#gamma, L_{proj} = %.1f fb^{-1}}", campaign.Data(), energy.Data(), EIClumi/1e15));
  tePICLabel_R->SetNDC();
  tePICLabel_R->SetTextSize(0.05);
  tePICLabel_R->Draw("same");
  // Legend
  TLegend* lxB = new TLegend(0.65, 0.6, 0.99, 0.81);
  lxB->SetLineWidth(0);
  lxB->SetFillStyle(0);
  lxB->AddEntry(h_xB_MC, "MC gen.", "l");
  lxB->AddEntry(h_xB_Rec, "Raw reco.", "lp");
  lxB->AddEntry(h_xB_Corr, "Corr. reco.", "lp");
  lxB->Draw();
  
  // LOWER - Ratio/efficiency
  lowerxb->cd();
  gPad->SetLogx();
  // Markers and lines
  h_xB_Eff->SetLineColor(kBlack);
  h_xB_Eff->SetLineWidth(2);
  // Axes
  h_xB_Eff->GetXaxis()->SetTitle("x_{B}");
  h_xB_Eff->GetXaxis()->SetTitleSize(0.11);
  h_xB_Eff->GetXaxis()->SetTitleOffset(0.93);
  h_xB_Eff->GetXaxis()->SetLabelSize(0.11);
  h_xB_Eff->GetYaxis()->SetTitle("Acc/MC");
  h_xB_Eff->GetYaxis()->SetTitleSize(0.13);
  h_xB_Eff->GetYaxis()->SetTitleOffset(0.4);
  h_xB_Eff->GetYaxis()->SetLabelSize(0.13);
  // Draw
  h_xB_Eff->Draw();
  // Save figure
  if(kSAVE) cxB->SaveAs(sFigsPath + "TDR_" + energy +"_xB.png");
  cxB->Close();

  // CANVAS: xB 2D
  TCanvas* cxB_2D = new TCanvas("cxb_2d","",1000,1000);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetLogz();
  gPad->SetRightMargin(0.12);
  h_xB_2D->GetXaxis()->SetTitle("x_{B,MC}");
  h_xB_2D->GetYaxis()->SetTitle("x_{B,reco}");
  h_xB_2D->Draw("colz");
  // Save figure
  if(kSAVE) cxB_2D->SaveAs(sFigsPath + "TDR_" + energy +"_xB_2D.png");
  cxB_2D->Close();

  // CANVAS: xB resolution
  TCanvas* cxB_Res = new TCanvas("cxb_res","",1000,1000);
  gPad->SetLogx();
  gPad->SetLogz();
  gPad->SetRightMargin(0.12);
  h_xB_Res->GetYaxis()->SetTitle("(x_{B,reco}-x_{B,MC})/x_{B,MC}");
  h_xB_Res->GetXaxis()->SetTitle("x_{B,MC}");
  h_xB_Res->Draw("colz");
  // Save figure
  if(kSAVE) cxB_Res->SaveAs(sFigsPath + "TDR_" + energy +"_xB_Res.png");
  cxB_Res->Close();

  // CANVAS: y distribution
  TCanvas* cy = new TCanvas("cy","",1000,1000);
  TPad* uppery = new TPad("uppery","",0.05,0.3,0.95,0.95);
  uppery->SetBottomMargin(0);
  uppery->Draw();
  TPad* lowery = new TPad("lowery","",0.05,0.05,0.95,0.3);
  lowery->SetBottomMargin(0.22);
  lowery->SetTopMargin(0);
  lowery->Draw();
  // UPPER - Distribution
  uppery->cd();
  gPad->SetLogy();
  // Markers and lines
  h_y_MC->SetLineColor(kBlack);
  h_y_MC->SetLineWidth(2);
  h_y_Rec->SetMarkerColor(kP6Blue);
  h_y_Rec->SetMarkerStyle(kOpenCircle);
  h_y_Rec->SetMarkerSize(1.5);
  h_y_Rec->SetLineColor(kP6Blue);
  h_y_Rec->SetLineWidth(2);
  h_y_Corr->SetMarkerColor(kP6Blue);
  h_y_Corr->SetMarkerStyle(kFullCircle);
  h_y_Corr->SetMarkerSize(1.5);
  h_y_Corr->SetLineColor(kP6Blue);
  h_y_Corr->SetLineWidth(2);
  // Axes
  h_y_MC->GetYaxis()->SetTitle("Counts");
  h_y_MC->GetYaxis()->SetTitleSize(0.05);
  h_y_MC->GetYaxis()->SetLabelSize(0.05);
  // Drawing
  h_y_MC->Draw("hist");
  h_y_Rec->Draw("pe same");
  h_y_Corr->Draw("pe same");
  // Labels
  tePICLabel->Draw("same");
  // Legend
  TLegend* lY = new TLegend(0.17, 0.07, 0.51, 0.28);
  lY->SetLineWidth(0);
  lY->SetFillStyle(0);
  lY->AddEntry(h_y_MC, "MC gen.", "l");
  lY->AddEntry(h_y_Rec, "Raw reco.", "lp");
  lY->AddEntry(h_y_Corr, "Corr. reco.", "lp");
  lY->Draw();
  
  // LOWER - Ratio/efficiency
  lowery->cd();
  // Markers and lines
  h_y_Eff->SetLineColor(kBlack);
  h_y_Eff->SetLineWidth(2);
  // Axes
  h_y_Eff->GetXaxis()->SetTitle("y");
  h_y_Eff->GetXaxis()->SetTitleSize(0.11);
  h_y_Eff->GetXaxis()->SetTitleOffset(0.93);
  h_y_Eff->GetXaxis()->SetLabelSize(0.11);
  h_y_Eff->GetYaxis()->SetTitle("Acc/MC");
  h_y_Eff->GetYaxis()->SetTitleSize(0.13);
  h_y_Eff->GetYaxis()->SetTitleOffset(0.4);
  h_y_Eff->GetYaxis()->SetLabelSize(0.13);
  // Draw
  h_y_Eff->Draw();
  // Save figure
  if(kSAVE) cy->SaveAs(sFigsPath + "TDR_" + energy +"_y.png");
  cy->Close();

  // CANVAS: y 2D
  TCanvas* cy_2D = new TCanvas("cy_2d","",1000,1000);
  gPad->SetLogz();
  gPad->SetRightMargin(0.12);
  h_y_2D->GetXaxis()->SetTitle("y_{MC}");
  h_y_2D->GetYaxis()->SetTitle("y_{reco}");
  h_y_2D->Draw("colz");
  // Save figure
  if(kSAVE) cy_2D->SaveAs(sFigsPath + "TDR_" + energy +"_y_2D.png");
  cy_2D->Close();

  // CANVAS: y resolution
  TCanvas* cy_Res = new TCanvas("cy_res","",1000,1000);
  gPad->SetLogz();
  gPad->SetRightMargin(0.12);
  h_y_Res->GetYaxis()->SetTitle("(y_{reco}-y_{MC})/y_{MC}");
  h_y_Res->GetYaxis()->SetTitle("y_{MC}");
  h_y_Res->Draw("colz");
  // Save figure
  if(kSAVE) cy_Res->SaveAs(sFigsPath + "TDR_" + energy +"_y_Res.png");
  cy_Res->Close();

  
  // CANVAS: t distribution
  TCanvas* ct = new TCanvas("ct","",1000,1000);
  TPad* uppert = new TPad("uppery","",0.05,0.3,0.95,0.95);
  uppert->SetBottomMargin(0);
  uppert->Draw();
  TPad* lowert = new TPad("lowery","",0.05,0.05,0.95,0.3);
  lowert->SetBottomMargin(0.22);
  lowert->SetTopMargin(0);
  lowert->Draw();
  // UPPER - Distribution
  uppert->cd();
  gPad->SetLogy();
  // Markers and lines
  h_t_MC->SetLineColor(kBlack);
  h_t_MC->SetLineWidth(2);
  h_t_B0Rec->SetMarkerColor(kP6Blue);
  h_t_B0Rec->SetMarkerStyle(kOpenCircle);
  h_t_B0Rec->SetMarkerSize(1.5);
  h_t_B0Rec->SetLineColor(kP6Blue);
  h_t_B0Rec->SetLineWidth(2);
  h_t_B0Corr->SetMarkerColor(kP6Blue);
  h_t_B0Corr->SetMarkerStyle(kFullCircle);
  h_t_B0Corr->SetMarkerSize(1.5);
  h_t_B0Corr->SetLineColor(kP6Blue);
  h_t_B0Corr->SetLineWidth(2);
  h_t_RPRec->SetMarkerColor(kP6Grape);
  h_t_RPRec->SetMarkerStyle(kOpenCircle);
  h_t_RPRec->SetMarkerSize(1.5);
  h_t_RPRec->SetLineColor(kP6Grape);
  h_t_RPRec->SetLineWidth(2);
  h_t_RPCorr->SetMarkerColor(kP6Grape);
  h_t_RPCorr->SetMarkerStyle(kFullCircle);
  h_t_RPCorr->SetMarkerSize(1.5);
  h_t_RPCorr->SetLineColor(kP6Grape);
  h_t_RPCorr->SetLineWidth(2);
  h_t_LCRec->SetMarkerColor(kP6Yellow);
  h_t_LCRec->SetMarkerStyle(kOpenCircle);
  h_t_LCRec->SetMarkerSize(1.5);
  h_t_LCRec->SetLineColor(kP6Yellow);
  h_t_LCRec->SetLineWidth(2);
  h_t_LCCorr->SetMarkerColor(kP6Yellow);
  h_t_LCCorr->SetMarkerStyle(kFullCircle);
  h_t_LCCorr->SetMarkerSize(1.5);
  h_t_LCCorr->SetLineColor(kP6Yellow);
  h_t_LCCorr->SetLineWidth(2);
  // Axes
  h_t_MC->GetYaxis()->SetTitle("Counts / 0.1 GeV^{2}");
  h_t_MC->GetYaxis()->SetTitleSize(0.05);
  h_t_MC->GetYaxis()->SetLabelSize(0.05);
  h_t_MC->GetYaxis()->SetRangeUser(1,10*h_t_MC->GetMaximum());
  // Drawing
  h_t_MC->Draw("hist");
  h_t_B0Rec->Draw("pe same");
  h_t_RPRec->Draw("pe same");
  h_t_LCRec->Draw("pe same");
  h_t_B0Corr->Draw("pe same");
  h_t_RPCorr->Draw("pe same");
  h_t_LCCorr->Draw("pe same");
  // Labels
  tePICLabel_R->Draw("same");
  // Legend
  TLegend* lt = new TLegend(0.17, 0.05, 0.71, 0.25);
  lt->SetLineWidth(0);
  lt->SetFillStyle(0);
  lt->SetNColumns(2);
  lt->AddEntry(h_t_MC, "MC gen.", "l");
  lt->AddEntry((TObject*)0, "", "");
  lt->AddEntry(h_t_B0Rec, "Raw reco. (B0)", "lp");
  lt->AddEntry(h_t_B0Corr, "Corr. reco. (B0)", "lp");
  lt->AddEntry(h_t_RPRec, "Raw reco. (RP)", "lp");
  lt->AddEntry(h_t_RPCorr, "Corr. reco. (RP)", "lp");
  lt->AddEntry(h_t_LCRec, "Raw reco. (e'#gamma)", "lp");
  lt->AddEntry(h_t_LCCorr, "Corr. reco. (e'#gamma)", "lp");
  lt->Draw();
  
  // LOWER - Ratio/efficiency
  lowert->cd();
  // Markers and lines
  h_t_B0Eff->SetLineColor(kP6Blue);
  h_t_B0Eff->SetLineWidth(2);
  h_t_RPEff->SetLineColor(kP6Grape);
  h_t_RPEff->SetLineWidth(2);
  h_t_LCEff->SetLineColor(kP6Yellow);
  h_t_LCEff->SetLineWidth(2);
  // Axes
  h_t_LCEff->GetXaxis()->SetTitle("|t| [GeV^{2}]");
  h_t_LCEff->GetXaxis()->SetTitleSize(0.11);
  h_t_LCEff->GetXaxis()->SetTitleOffset(0.93);
  h_t_LCEff->GetXaxis()->SetLabelSize(0.11);
  h_t_LCEff->GetYaxis()->SetTitle("Acc/MC");
  h_t_LCEff->GetYaxis()->SetTitleSize(0.13);
  h_t_LCEff->GetYaxis()->SetTitleOffset(0.4);
  h_t_LCEff->GetYaxis()->SetLabelSize(0.13);
  h_t_LCEff->GetYaxis()->SetNdivisions(505);
  h_t_LCEff->GetYaxis()->SetRangeUser(0.,1.05);
  // Draw
  h_t_LCEff->Draw();
  h_t_B0Eff->Draw("same");
  h_t_RPEff->Draw("same");
  // Save figure
  if(kSAVE) ct->SaveAs(sFigsPath + "TDR_" + energy +"_t.png");
  ct->Close();


  // CANVAS: Overlaid t-resolutions (absolute)
  TH1D* h_tResB0_Proj = (TH1D*)h_t_B0Res->ProjectionY("tresb0_py");
  TH1D* h_tResRP_Proj = (TH1D*)h_t_RPRes->ProjectionY("tresrp_py");
  TH1D* h_tResLC_Proj = (TH1D*)h_t_LCRes->ProjectionY("treslc_py");

  TCanvas* ctRes_All = new TCanvas("ctres_all","",1000,1000);
  // Set draw options - markers and lines
  h_tResB0_Proj->SetLineColor(kP6Blue);
  h_tResB0_Proj->SetFillColor(kP6Blue);
  h_tResB0_Proj->SetFillStyle(3001);
  h_tResRP_Proj->SetLineColor(kP6Grape);
  h_tResRP_Proj->SetFillColor(kP6Grape);
  h_tResRP_Proj->SetFillStyle(3001);
  h_tResLC_Proj->SetLineColor(kP6Yellow);
  h_tResLC_Proj->SetFillColor(kP6Yellow);
  h_tResLC_Proj->SetFillStyle(3001);
  // Set draw options - axes
  // Using Method L plot for baseline
  gPad->SetLogy();
  h_tResLC_Proj->GetXaxis()->SetTitle("#Deltat [GeV^{2}]");
  h_tResLC_Proj->GetXaxis()->SetRangeUser(-2.,2.);
  h_tResLC_Proj->GetYaxis()->SetTitle("Counts/0.02 GeV^{ 2}");
  h_tResLC_Proj->GetYaxis()->SetRangeUser(1,1000*h_tResLC_Proj->GetMaximum());
  // Draw
  // Method L first, others depending on which has more entries
  h_tResLC_Proj->Draw("hist");
  if(h_tResB0_Proj->GetEntries() > h_tResRP_Proj->GetEntries()){
    h_tResB0_Proj->Draw("hist same");
    h_tResRP_Proj->Draw("hist same");
  }
  else{
    h_tResRP_Proj->Draw("hist same");
    h_tResB0_Proj->Draw("hist same");
  }
  // Add text
  TLatex* tePICLabel_tres = new TLatex(-1.9, 0.25*h_tResLC_Proj->GetMaximum(),
				       Form("#splitline{#bf{ePIC} Performance %s, %s GeV}{ep #rightarrow e'p'#gamma}", campaign.Data(), energy.Data()));
  tePICLabel_tres->SetTextSize(0.037);
  tePICLabel_tres->Draw("same");
  // Add legend
  TLegend* ltRes = new TLegend(0.16,0.70,0.45,0.84);
  ltRes->SetLineWidth(0);
  ltRes->SetFillStyle(0);
  ltRes->SetTextSize(0.037);
  TString sB0Res_leg = "BABE (B0), RMS = " + combineScientific(h_tResB0_Proj->GetRMS(),h_tResB0_Proj->GetRMSError()) + " GeV^{2}";
  TString sRPRes_leg = "BABE (RP), RMS = " + combineScientific(h_tResRP_Proj->GetRMS(),h_tResRP_Proj->GetRMSError()) + " GeV^{2}";
  TString sLCRes_leg = "eXBE, RMS = " + combineScientific(h_tResLC_Proj->GetRMS(),h_tResLC_Proj->GetRMSError()) + " GeV^{2}";
  ltRes->AddEntry(h_tResB0_Proj,sB0Res_leg,"lp");
  ltRes->AddEntry(h_tResRP_Proj,sRPRes_leg,"lp");
  ltRes->AddEntry(h_tResLC_Proj,sLCRes_leg,"lp");
  ltRes->Draw();
  // Save figure
  if(kSAVE) ctRes_All->SaveAs(sFigsPath + "TDR_" + energy +"_tRes_All.png");
  ctRes_All->Close();

  // CANVAS - Overlaid t-resolutions (relative)
  TH1D* h_tResB0Pct_Proj = (TH1D*)h_t_B0ResPct->ProjectionY("tresb0pct_py");
  TH1D* h_tResRPPct_Proj = (TH1D*)h_t_RPResPct->ProjectionY("tresrppct_py");
  TH1D* h_tResLCPct_Proj = (TH1D*)h_t_LCResPct->ProjectionY("treslcpct_py");

  TCanvas* ctResPct_All = new TCanvas("ctrespct_all","",1000,1000);
  // Set draw options - markers and lines
  h_tResB0Pct_Proj->SetLineColor(kP6Blue);
  h_tResB0Pct_Proj->SetFillColor(kP6Blue);
  h_tResB0Pct_Proj->SetFillStyle(3001);
  h_tResRPPct_Proj->SetLineColor(kP6Grape);
  h_tResRPPct_Proj->SetFillColor(kP6Grape);
  h_tResRPPct_Proj->SetFillStyle(3001);
  h_tResLCPct_Proj->SetLineColor(kP6Yellow);
  h_tResLCPct_Proj->SetFillColor(kP6Yellow);
  h_tResLCPct_Proj->SetFillStyle(3001);
  // Set draw options - axes
  // Using Method L plot for baseline
  h_tResLCPct_Proj->GetXaxis()->SetTitle("#Deltat/|t_{MC}|");
  h_tResLCPct_Proj->GetXaxis()->SetRangeUser(-0.05,1.);
  h_tResLCPct_Proj->GetYaxis()->SetRangeUser(1,1.25*h_tResLCPct_Proj->GetMaximum());
  // Draw
  // Method L first, others depending on which has more entries
  h_tResLCPct_Proj->Draw("hist");
  if(h_tResB0Pct_Proj->GetEntries() > h_tResRPPct_Proj->GetEntries()){
    h_tResB0Pct_Proj->Draw("hist same");
    h_tResRPPct_Proj->Draw("hist same");
  }
  else{
    h_tResRPPct_Proj->Draw("hist same");
    h_tResB0Pct_Proj->Draw("hist same");
  }
  // Add text
  TLatex* tePICLabel_trespct = new TLatex(-0.01, 0.91*h_tResLCPct_Proj->GetMaximum(),
					  Form("#splitline{#bf{ePIC} Performance %s, %s GeV}{ep #rightarrow e'p'#gamma}", campaign.Data(), energy.Data()));
  tePICLabel_trespct->SetTextSize(0.037);
  tePICLabel_trespct->Draw("same");
  // Add legend
  TLegend* ltResPct = new TLegend(0.30,0.64,0.64,0.78);
  ltResPct->SetLineWidth(0);
  ltResPct->SetFillStyle(0);
  ltResPct->SetTextSize(0.032);
  TString sB0ResPct_leg = "BABE (B0), RMS = " + combineScientific(h_tResB0Pct_Proj->GetRMS(),h_tResB0Pct_Proj->GetRMSError());
  TString sRPResPct_leg = "BABE (RP), RMS = " + combineScientific(h_tResRPPct_Proj->GetRMS(),h_tResRPPct_Proj->GetRMSError());
  TString sLCResPct_leg = "eXBE, RMS = " + combineScientific(h_tResLCPct_Proj->GetRMS(),h_tResLCPct_Proj->GetRMSError());
  ltResPct->AddEntry(h_tResB0Pct_Proj,sB0ResPct_leg,"lf");
  ltResPct->AddEntry(h_tResRPPct_Proj,sRPResPct_leg,"lf");
  ltResPct->AddEntry(h_tResLCPct_Proj,sLCResPct_leg,"lf");
  ltResPct->Draw();
  // Save figure
  if(kSAVE) ctResPct_All->SaveAs(sFigsPath + "TDR_" + energy +"_tResPct_All.png");
  ctResPct_All->Close();


  // CANVAS: Event (E-pz)
  TCanvas* cEmPz = new TCanvas("cempz","",1000,1000);
  gPad->SetLogy();
  // Markers and lines
  h_EmPz_MC->SetLineColor(kBlack);
  h_EmPz_MC->SetLineWidth(2);
  h_EmPz_RP->SetLineColor(kP6Blue);
  h_EmPz_RP->SetLineWidth(2);
  // Axes
  h_EmPz_MC->GetXaxis()->SetTitle("(E-p_{z})_{e'p'#gamma} [GeV]");
  h_EmPz_MC->GetXaxis()->SetRangeUser(0.,30.);
  h_EmPz_MC->GetYaxis()->SetTitle("Counts / 0.25 GeV");
  float empz_high = 0.8*h_EmPz_MC->GetMaximum();
  h_EmPz_MC->GetYaxis()->SetRangeUser(0.5,50*empz_high);
  // Draw
  h_EmPz_MC->Draw("hist");
  h_EmPz_RP->Draw("hist same");
  // Text
  TLatex* tePICLabel_empz = new TLatex(0.17, 0.86, 
				  Form("#splitline{#bf{ePIC} Performance %s, %s GeV}{ep #rightarrow e'p'#gamma, L_{proj} = %.1f fb^{-1}}", campaign.Data(), energy.Data(), EIClumi/1e15));
  tePICLabel_empz->SetNDC();
  tePICLabel_empz->SetTextSize(0.04);
  tePICLabel_empz->Draw("same");
  // Lines
  TLine* lLower = new TLine(15.,0.,15.,empz_high);
  lLower->SetLineColor(kP6Red);
  lLower->SetLineWidth(2);
  lLower->Draw();
  TLine* lUpper = new TLine(25.,0.,25.,empz_high);
  lUpper->SetLineColor(kP6Red);
  lUpper->SetLineWidth(2);
  lUpper->Draw();
  // Legend
  TLegend* lEmPz = new TLegend(0.18,0.58,0.53,0.72);
  lEmPz->SetLineWidth(0);
  lEmPz->SetFillStyle(0);
  lEmPz->SetTextSize(0.032);
  lEmPz->AddEntry(h_EmPz_MC,"MC gen.","l");
  lEmPz->AddEntry(h_EmPz_RP,"Raw reco.","l");
  lEmPz->Draw();
  // Save figure
  if(kSAVE) cEmPz->SaveAs(sFigsPath + "TDR_" + energy +"_ct_All.png");
  cEmPz->Close();
  

  // CANVAS: Missing pT
  TCanvas* cpTmiss3 = new TCanvas("cempz","",1000,1000);
  gPad->SetLogy();
  // Markers and lines
  h_pTmiss_RP->SetLineColor(kBlack);
  h_pTmiss_RP->SetLineWidth(2);
  // Axes
  h_pTmiss_RP->GetXaxis()->SetTitle("p_{T,miss} [GeV]");
  h_pTmiss_RP->GetXaxis()->SetTitleOffset(1.10);
  h_pTmiss_RP->GetYaxis()->SetTitle("Counts / 0.01 GeV");
  h_pTmiss_RP->GetYaxis()->SetRangeUser(0.5,5*h_pTmiss_RP->GetMaximum());
  // Draw
  h_pTmiss_RP->Draw("hist");
  // Text
  TLatex* tePICLabel_ptm3 = new TLatex(0.3, 0.87, 
				       Form("#splitline{#bf{ePIC} Performance %s, %s GeV}{ep #rightarrow e'p'#gamma, L_{proj} = %.1f fb^{-1}}", campaign.Data(), energy.Data(), EIClumi/1e15));
  tePICLabel_ptm3->SetNDC();
  tePICLabel_ptm3->SetTextSize(0.04);
  tePICLabel_ptm3->Draw("same");
  // Lines
  TLine* lptm3 = new TLine(0.5,0.,0.5,0.75*h_pTmiss_RP->GetMaximum());
  lptm3->SetLineColor(kP6Red);
  lptm3->SetLineWidth(2);
  lptm3->Draw();
  // Save figure
  if(kSAVE) cpTmiss3->SaveAs(sFigsPath + "TDR_" + energy +"_pTmiss3.png");
  cpTmiss3->Close();


  // CANVAS: 2d event x/Q2
  TCanvas* cxBvQ2 = new TCanvas("cxbvq2","",1000,1000);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetLogz();
  gPad->SetRightMargin(0.12);
  h_xBvQ2_Rec->GetXaxis()->SetTitle("x_{B}");
  h_xBvQ2_Rec->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
  h_xBvQ2_Rec->Draw("colz");
  // Save figure
  if(kSAVE) cxBvQ2->SaveAs(sFigsPath + "TDR_" + energy +"_xBvQ2_Rec.png");
  cxBvQ2->Close();


  // CANVAS: 2d event x/Q2
  TCanvas* cxBvt = new TCanvas("cxbvt","",1000,1000);
  gPad->SetLogx();
  gPad->SetLogz();
  gPad->SetRightMargin(0.12);
  h_xBvt_Rec->GetXaxis()->SetTitle("x_{B}");
  h_xBvt_Rec->GetYaxis()->SetTitle("|t| [GeV^{2}]");
  h_xBvt_Rec->Draw("colz");
  // Save figure
  if(kSAVE) cxBvt->SaveAs(sFigsPath + "TDR_" + energy +"_xBvt_Rec.png");
  cxBvQ2->Close();

  
  return;
}
