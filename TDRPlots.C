using namespace std;

#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include "../ePIC_style.C"


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
// MAIN
//---------------------------------------------------------------------
// Parameter `energy` declares which energy setting to use for most plots (all except t-resolution)
void TDRPlots(TString energy = "10x100"){
  TString campaign="25.10.2";

  cout<<"Processing TDR plots"<<endl;
  cout<<"Campaigin used: "<<campaign<<endl;
  
  //---------------------------------------------------------------------
  // Get data file from DVCSAnalysis script
  //---------------------------------------------------------------------
  // Need plots for all three 'standard' EIC energy settings
  TFile* f5x41   = TFile::Open("./ePIC_DVCS_TDR_5x41.root");
  TFile* f10x100 = TFile::Open("./ePIC_DVCS_TDR_10x100.root");
  TFile* f18x275 = TFile::Open("./ePIC_DVCS_TDR_18x275.root");
  
  //---------------------------------------------------------------------
  // Prepare histograms for single energy plots
  // Will be pulled from files separately depending on energy
  //---------------------------------------------------------------------
  // Plot 1: Detector occupancies (eta) for all species
  TH1D* h_eta_MCp;
  TH1D* h_eta_MCe;
  TH1D* h_eta_MCg;
  TH1D* h_eta_RPp;
  TH1D* h_eta_RPe;
  TH1D* h_eta_RPg;
  TH1D* h_eta_RPPp;

  // Plot 2: Photon angular resolution
  TH1D* h_PhotRes_theta;
  TH2D* h_PhotRes2D_theta;

  // Plot 3: Mandelstam t-distribution
  TH1D* h_t_Truth;
  TH1D* h_t_B0Acc;
  TH1D* h_t_RPAcc;
  TH1D* h_t_B0Reco;
  TH1D* h_t_RPReco;
  
  //---------------------------------------------------------------------
  // Extract histograms
  //---------------------------------------------------------------------
  // Case 1: 5x41 GeV default
  if(energy == "5x41"){
    h_eta_MCp   = (TH1D*)f5x41->Get("eta_MCp")->Clone("eta_MC_prot");
    h_eta_MCe   = (TH1D*)f5x41->Get("eta_MCe")->Clone("eta_MC_elec");
    h_eta_MCg   = (TH1D*)f5x41->Get("eta_MCg")->Clone("eta_MC_phot");
    h_eta_RPp   = (TH1D*)f5x41->Get("eta_RPp")->Clone("eta_RE_B0prot");
    h_eta_RPPp  = (TH1D*)f5x41->Get("eta_RPPp")->Clone("eta_RE_RPprot");
    h_eta_RPe   = (TH1D*)f5x41->Get("eta_RPe")->Clone("eta_RE_elec");
    h_eta_RPg   = (TH1D*)f5x41->Get("eta_RPg")->Clone("eta_RE_phot");
    
    h_PhotRes_theta   = (TH1D*)f5x41->Get("photres_theta")->Clone("thetares_phot");
    h_PhotRes2D_theta = (TH2D*)f5x41->Get("photres2d_theta")->Clone("thetares_phot_2d");
    
    h_t_Truth  = (TH1D*)f5x41->Get("t_truth");
    h_t_B0Acc  = (TH1D*)f5x41->Get("t_b0acc");
    h_t_RPAcc  = (TH1D*)f5x41->Get("t_rpacc");
    h_t_B0Reco = (TH1D*)f5x41->Get("t_b0reco");
    h_t_RPReco = (TH1D*)f5x41->Get("t_rpreco");
  }
  // Case 2: 10x100 GeV default
  else if(energy == "10x100"){
    h_eta_MCp   = (TH1D*)f10x100->Get("eta_MCp")->Clone("eta_MC_prot");
    h_eta_MCe   = (TH1D*)f10x100->Get("eta_MCe")->Clone("eta_MC_elec");
    h_eta_MCg   = (TH1D*)f10x100->Get("eta_MCg")->Clone("eta_MC_phot");
    h_eta_RPp   = (TH1D*)f10x100->Get("eta_RPp")->Clone("eta_RE_B0prot");
    h_eta_RPPp  = (TH1D*)f10x100->Get("eta_RPPp")->Clone("eta_RE_RPprot");
    h_eta_RPe   = (TH1D*)f10x100->Get("eta_RPe")->Clone("eta_RE_elec");
    h_eta_RPg   = (TH1D*)f10x100->Get("eta_RPg")->Clone("eta_RE_phot");
    
    h_PhotRes_theta   = (TH1D*)f10x100->Get("photres_theta")->Clone("thetares_phot");
    h_PhotRes2D_theta = (TH2D*)f10x100->Get("photres2d_theta")->Clone("thetares_phot_2d");

    h_t_Truth  = (TH1D*)f10x100->Get("t_truth");
    h_t_B0Acc  = (TH1D*)f10x100->Get("t_b0acc");
    h_t_RPAcc  = (TH1D*)f10x100->Get("t_rpacc");
    h_t_B0Reco = (TH1D*)f10x100->Get("t_b0reco");
    h_t_RPReco = (TH1D*)f10x100->Get("t_rpreco");
  }
  // Case 3: 18x275 GeV default
  else if(energy == "18x275"){
    h_eta_MCp   = (TH1D*)f18x275->Get("eta_MCp")->Clone("eta_MC_prot");
    h_eta_MCe   = (TH1D*)f18x275->Get("eta_MCe")->Clone("eta_MC_elec");
    h_eta_MCg   = (TH1D*)f18x275->Get("eta_MCg")->Clone("eta_MC_phot");
    h_eta_RPp   = (TH1D*)f18x275->Get("eta_RPp")->Clone("eta_RE_B0prot");
    h_eta_RPPp  = (TH1D*)f18x275->Get("eta_RPPp")->Clone("eta_RE_RPprot");
    h_eta_RPe   = (TH1D*)f18x275->Get("eta_RPe")->Clone("eta_RE_elec");
    h_eta_RPg   = (TH1D*)f18x275->Get("eta_RPg")->Clone("eta_RE_phot");
    
    h_PhotRes_theta = (TH1D*)f18x275->Get("photres_theta")->Clone("thetares_phot");
    h_PhotRes2D_theta = (TH2D*)f18x275->Get("photres2d_theta")->Clone("thetares_phot_2d");

    h_t_Truth  = (TH1D*)f18x275->Get("t_truth");
    h_t_B0Acc  = (TH1D*)f18x275->Get("t_b0acc");
    h_t_RPAcc  = (TH1D*)f18x275->Get("t_rpacc");
    h_t_B0Reco = (TH1D*)f18x275->Get("t_b0reco");
    h_t_RPReco = (TH1D*)f18x275->Get("t_rpreco");
  }
  // Call errors
  h_eta_MCp->Sumw2();
  h_eta_MCe->Sumw2();
  h_eta_MCg->Sumw2();
  h_eta_RPp->Sumw2();
  h_eta_RPe->Sumw2();
  h_eta_RPg->Sumw2();
  h_eta_RPPp->Sumw2();
  h_PhotRes_theta->Sumw2();
  h_PhotRes2D_theta->Sumw2();
  h_t_Truth->Sumw2();
  h_t_B0Acc->Sumw2();
  h_t_RPAcc->Sumw2();
  h_t_B0Reco->Sumw2();
  h_t_RPReco->Sumw2();
  
  // Plot 4: t-resolution
  // Need from all 3 energy settings anyway
  TH2D* h_tresb0_2d_5x41 = (TH2D*)f5x41->Get("h_tresb0pct_2d")->Clone("tresb0pct_2d_5x41");
  TH2D* h_tresrp_2d_5x41 = (TH2D*)f5x41->Get("h_tresrppct_2d")->Clone("tresrppct_2d_5x41");
  TH2D* h_tresb0_2d_10x100 = (TH2D*)f10x100->Get("h_tresb0pct_2d")->Clone("tresb0pct_2d_10x100");
  TH2D* h_tresrp_2d_10x100 = (TH2D*)f10x100->Get("h_tresrppct_2d")->Clone("tresrppct_2d_10x100");
  TH2D* h_tresb0_2d_18x275 = (TH2D*)f18x275->Get("h_tresb0pct_2d")->Clone("tresb0pct_2d_18x275");
  TH2D* h_tresrp_2d_18x275 = (TH2D*)f18x275->Get("h_tresrppct_2d")->Clone("tresrppct_2d_18x275");
  h_tresb0_2d_5x41->Sumw2();
  h_tresrp_2d_5x41->Sumw2();
  h_tresb0_2d_10x100->Sumw2();
  h_tresrp_2d_10x100->Sumw2();
  h_tresb0_2d_18x275->Sumw2();
  h_tresrp_2d_18x275->Sumw2();
   
  TH1D* h_tresB0_5x41 = extractRMSSlice("extracted_tb0_rms_5x41",h_tresb0_2d_5x41);
  TH1D* h_tresB0_10x100 = extractRMSSlice("extracted_tb0_rms_10x100",h_tresb0_2d_10x100);
  TH1D* h_tresB0_18x275 = extractRMSSlice("extracted_tb0_rms_18x275",h_tresb0_2d_18x275);
  TH1D* h_tresRP_5x41 = extractRMSSlice("extracted_trp_rms_5x41",h_tresrp_2d_5x41);
  TH1D* h_tresRP_10x100 = extractRMSSlice("extracted_trp_rms_10x100",h_tresrp_2d_10x100);
  TH1D* h_tresRP_18x275 = extractRMSSlice("extracted_trp_rms_18x275",h_tresrp_2d_18x275);

  //---------------------------------------------------------------------
  // Calculations on histograms
  // Scaling to 5fb-1
  //---------------------------------------------------------------------
  // Calculations - SCALING FACTOR TO 5FB-1
  Double_t fXSint{0};
  if(energy == "5x41") fXSint = 0.471255145271248e-9;
  else if(energy == "10x100") fXSint = 0.510611375008106e-9;
  else if(energy == "18x275") fXSint = 0.556957853417171e-9;

  Double_t lumi = 1e6/fXSint;
  Double_t scaleTo5 = 5e15/lumi;

  // Apply scalings
  h_eta_MCp->Scale(scaleTo5);
  h_eta_MCe->Scale(scaleTo5);
  h_eta_MCg->Scale(scaleTo5);
  h_eta_RPp->Scale(scaleTo5);
  h_eta_RPe->Scale(scaleTo5);
  h_eta_RPg->Scale(scaleTo5);
  h_eta_RPPp->Scale(scaleTo5);

  //---------------------------------------------------------------------
  // Detector efficiency corrections
  //---------------------------------------------------------------------
  
  // Calculate detector acceptances
  // Duplicate MCA histograms
  TH1D* h_t_B0Eff = (TH1D*)h_t_B0Acc->Clone("t_b0eff");
  TH1D* h_t_RPEff = (TH1D*)h_t_RPAcc->Clone("t_rpeff");
  // Divide by MC truth for acceptance ratios
  h_t_B0Eff->Divide(h_t_Truth);
  h_t_RPEff->Divide(h_t_Truth);

  // Correct for detector acceptances (need separate for B0 and RP)
  // Duplicate reco. histograms
  TH1D* h_t_B0Corr = (TH1D*)h_t_B0Reco->Clone("t_b0corr");
  TH1D* h_t_RPCorr = (TH1D*)h_t_RPReco->Clone("t_rpcorr");
  // Divide by efficiencies
  h_t_B0Corr->Divide(h_t_B0Eff);
  h_t_RPCorr->Divide(h_t_RPEff);

  //---------------------------------------------------------------------
  // Combine corrected points
  //---------------------------------------------------------------------
  TH1D* h_t_GoodCorr = (TH1D*)h_t_Truth->Clone("t_goodcorr");
  h_t_GoodCorr->Reset();
  TH1D* h_t_GoodReco = (TH1D*)h_t_Truth->Clone("t_goodreco"); // FOR PLOTTING LATER
  h_t_GoodReco->Reset(); 

  // Combine points into "good" histograms
  // Use points where acceptance > 5%
  for(int bin{1}; bin<h_t_Truth->GetNbinsX()+1; ++bin){    
    // Check B0 histograms
    if(h_t_B0Eff->GetBinContent(bin) >= 0.05){
      h_t_GoodCorr->SetBinContent(bin,h_t_B0Corr->GetBinContent(bin));
      h_t_GoodCorr->SetBinError(bin,h_t_B0Corr->GetBinError(bin));
      h_t_GoodReco->SetBinContent(bin,h_t_B0Reco->GetBinContent(bin));
      h_t_GoodReco->SetBinError(bin,h_t_B0Reco->GetBinError(bin));
    }
    // Then check RP
    // RP will overwrite B0 in cases of overlap
    if(h_t_RPEff->GetBinContent(bin) >= 0.05){
      h_t_GoodCorr->SetBinContent(bin,h_t_RPCorr->GetBinContent(bin));
      h_t_GoodCorr->SetBinError(bin,h_t_RPCorr->GetBinError(bin));
      h_t_GoodReco->SetBinContent(bin,h_t_RPReco->GetBinContent(bin));
      h_t_GoodReco->SetBinError(bin,h_t_RPReco->GetBinError(bin));
    }
  }
  
  // Fit combined points and extract result
  TFitResultPtr t_fitres = h_t_GoodCorr->Fit("expo","S0","",0,1.6);
  TF1* fitfunc = h_t_GoodCorr->GetFunction("expo");
  
  // Calculate errors on points
  // Need to extract points from temporary histogram to calculate errors
  vector<double> pts;
  for(int bin{1}; bin<h_t_Truth->GetNbinsX()+1;bin++){
    // If "good" histogram has data, use that
    if(h_t_GoodCorr->GetBinContent(bin) > 1e-6) pts.push_back(h_t_GoodCorr->GetBinContent(bin));
    // If not, evaluate fit at bin centre
    else pts.push_back(fitfunc->Eval(h_t_GoodCorr->GetBinCenter(bin)));
  }
  // Convert variable size array (vector pts) to fixed size (double x[size])
  const int Npoints = (int)pts.size();
  double x[Npoints];
  copy(pts.begin(), pts.end(), x);
  double err[Npoints];
  // Calculate 68.3% CI (1 sigma) from fit
  t_fitres->GetConfidenceIntervals(1,1,1,x,err,0.683,false);

  // Combine fit result with original corrected points
  TH1D* h_t_CorrFit = (TH1D*)h_t_Truth->Clone("t_corrfit");
  h_t_CorrFit->Reset();
  TH1D* h_t_CorrDraw = (TH1D*)h_t_Truth->Clone("t_corrdraw"); // For drawing purposes - isolate reco. points used
  h_t_CorrDraw->Reset();

  for(int bin{1}; bin<h_t_Truth->GetNbinsX()+1;bin++){
    // Use corrected reco. if used in fit
    if(h_t_GoodReco->GetBinContent(bin) > 1e-6){
      h_t_CorrFit->SetBinContent(bin,h_t_GoodCorr->GetBinContent(bin));
      h_t_CorrFit->SetBinError(bin,h_t_GoodCorr->GetBinError(bin));
      // And add to 'drawing' histogram
      h_t_CorrDraw->SetBinContent(bin,h_t_GoodCorr->GetBinContent(bin));
      h_t_CorrDraw->SetBinError(bin,h_t_GoodCorr->GetBinError(bin));
    }
    // If not used for fit, use fit result (and ignore drawing histogram)
    // Only use if truth histogram has entries (to ignore tails above generated limit)
    else if(h_t_Truth->GetBinContent(bin) > 1e-6){
      h_t_CorrFit->SetBinContent(bin,x[bin-1]);
      h_t_CorrFit->SetBinError(bin,err[bin-1]);
    }
  }

  
  //--------------------------------------------------------------------------------------
  // Draw histograms
  //--------------------------------------------------------------------------------------
  // Global style
  gROOT->ProcessLine("set_ePIC_style()");
  gStyle->SetCanvasPreferGL(kTRUE);
 
  // CANVAS: DETECTOR OCCUPANCIES (SEPARATED BY SPECIES)
  TCanvas* cEta = new TCanvas("cEta","",1500,500);
  cEta->Divide(3,1,-1);
  
  // Pad 1 - electrons
  cEta->cd(1);
  gPad->SetLogy();
  // Set Drawing options
  // Markers and Lines
  h_eta_MCe->SetLineColor(kBlack);
  h_eta_MCe->SetLineWidth(2);
  h_eta_RPe->SetLineColor(kBlack);
  h_eta_RPe->SetMarkerColor(kBlack);
  h_eta_RPe->SetMarkerStyle(20);
  h_eta_RPe->SetFillColorAlpha(kP6Gray,0.5);
  // Axes
  h_eta_MCe->SetMaximum(15*h_eta_MCe->GetMaximum());
  h_eta_MCe->GetYaxis()->SetRangeUser(0.9, h_eta_MCe->GetMaximum());
  h_eta_MCe->GetXaxis()->SetTitle("#eta");
  h_eta_MCe->GetYaxis()->SetTitleOffset(1.25);
  h_eta_MCe->GetYaxis()->SetTitle("Counts");
  // Draw
  h_eta_MCe->Draw("HIST");
  h_eta_RPe->Draw("ep same");
  // Add legend
  TLegend* lcEta_p1 = new TLegend(0.67,0.60,1.,0.73);
  lcEta_p1->AddEntry(h_eta_MCe, "#bf{EpIC} MC", "l");
  lcEta_p1->AddEntry(h_eta_RPe, "Raw reco.", "pl");
  lcEta_p1->Draw();
  // Add text
  TLatex* tePICLabel_eta = new TLatex(0.17, 0.92, "#splitline{#bf{ePIC} Performance " + campaign + ", 10x100 GeV}{ep #rightarrow e'p'#gamma, L_{proj} = 5 fb^{-1}}");
  tePICLabel_eta->SetNDC();
  tePICLabel_eta->Draw("same");
  TLatex* tLabel1_eta = new TLatex(0.68, 0.79, "#splitline{Single e^{-}}{Q^{2} #geq 1 GeV^{2}}");
  tLabel1_eta->SetNDC();
  tLabel1_eta->Draw("same");
  
  // Pad 2 - photons
  cEta->cd(2);
  gPad->SetLogy();
  // Markers and Lines
  h_eta_MCg->SetLineColor(kP6Red);
  h_eta_MCg->SetLineWidth(2);
  h_eta_RPg->SetLineColor(kP6Red);
  h_eta_RPg->SetMarkerColor(kP6Red);
  h_eta_RPg->SetMarkerStyle(20);
  // Axes
  h_eta_MCg->SetMaximum(h_eta_MCe->GetMaximum()); // Match maxima across species
  h_eta_MCg->GetYaxis()->SetRangeUser(0.9, h_eta_MCg->GetMaximum());
  h_eta_MCg->GetXaxis()->SetTitle("#eta");
  // Draw
  h_eta_MCg->Draw("HIST");
  h_eta_RPg->Draw("ep same");
  // Add text
  TLatex* tLabel2_eta = new TLatex(0.07, 0.92, "Single #gamma");
  tLabel2_eta->SetNDC();
  tLabel2_eta->Draw("same");
  // Add legend
  TLegend* lcEta_p2 = new TLegend(0.06,0.77,0.71,0.90);
  lcEta_p2->AddEntry(h_eta_MCg, "#bf{EpIC} MC", "l");
  lcEta_p2->AddEntry(h_eta_RPg, "Raw reco.", "pl");
  lcEta_p2->Draw();

  // Pad 3 - protons
  cEta->cd(3);
  gPad->SetLogy();
  // Markers and Lines
  h_eta_MCp->SetLineColor(kP6Blue);
  h_eta_MCp->SetLineWidth(2);
  h_eta_RPp->SetLineColor(kP6Violet);
  h_eta_RPp->SetMarkerColor(kP6Violet);
  h_eta_RPp->SetMarkerStyle(20);
  h_eta_RPPp->SetLineColor(kP6Blue);
  h_eta_RPPp->SetMarkerColor(kP6Blue);
  h_eta_RPPp->SetMarkerStyle(20);
  // Axes
  h_eta_MCp->SetMaximum(h_eta_MCe->GetMaximum()); // Match maxima across species
  h_eta_MCp->GetYaxis()->SetRangeUser(0.9, h_eta_MCp->GetMaximum());
  h_eta_MCp->GetXaxis()->SetTitle("#eta");
  h_eta_MCp->GetYaxis()->SetTitle("Counts");
  // Draw
  h_eta_MCp->Draw("HIST");
  h_eta_RPp->Draw("ep same");
  h_eta_RPPp->Draw("ep same");
  // Add legend
  TLegend* lcEta_p3 = new TLegend(0.05,0.74,0.72,0.9);
  lcEta_p3->AddEntry(h_eta_MCp, "#bf{EpIC} MC", "l");
  lcEta_p3->AddEntry(h_eta_RPp, "Reco. B0", "pl");
  lcEta_p3->AddEntry(h_eta_RPPp, "Reco. RP", "pl");
  lcEta_p3->Draw();
  // Add text
  TLatex* tLabel3_eta = new TLatex(0.07, 0.92, "Single p");
  tLabel3_eta->SetNDC();
  tLabel3_eta->Draw("same");
  
  // Save figure
  cEta->SaveAs("TDR_" + energy + "_Eta.png");

  
  // CANVAS: Photon angular resolution - side-by-side
  TCanvas* cPhotRes = new TCanvas("cPhotRes","",1500,750);
  cPhotRes->Divide(2,1);
  
  // Pad 1 - 1D resolution
  cPhotRes->cd(1);
  // Markers and Lines
  h_PhotRes_theta->SetLineColor(kBlack);
  h_PhotRes_theta->SetMarkerColor(kBlack);
  h_PhotRes_theta->SetMarkerStyle(20);
  // Axes
  h_PhotRes_theta->GetXaxis()->SetTitle("#Delta#theta_{#gamma} [deg]");
  h_PhotRes_theta->GetXaxis()->SetRangeUser(-10,10);
  h_PhotRes_theta->GetYaxis()->SetTitle("Counts / 0.5 deg");
  h_PhotRes_theta->GetYaxis()->SetTitleOffset(1.45);
  h_PhotRes_theta->GetYaxis()->SetRangeUser(0.1,1.5*h_PhotRes_theta->GetMaximum());
  // Draw
  h_PhotRes_theta->Draw("ep");
  // Add text
  TLatex* tePICLabel_photres = new TLatex(-8.8, 0.9*h_PhotRes_theta->GetMaximum(), "#splitline{#bf{ePIC} Performance " + campaign + ", " + energy + " GeV}{ep #rightarrow e'p'#gamma, all reco.#gamma}");
  tePICLabel_photres->SetTextSize(0.04);
  tePICLabel_photres->Draw("same");
  // Fit plot and extract fitted function
  TFitResultPtr fitres = h_PhotRes_theta->Fit("gaus","S0","",-10,10);
  // Add text: fitted parameters
  TString sFitRes = "#sigma = " + combineScientific(fitres->Parameter(2),fitres->ParError(2)) + " deg";
  TLatex* tFitRes = new TLatex(-8.8, 0.75*h_PhotRes_theta->GetMaximum(), sFitRes);
  tFitRes->Draw("same");
  
  // Pad 2 - 2D resolution against MC theta
  cPhotRes->cd(2);
  gPad->SetLogz();
  h_PhotRes2D_theta->GetYaxis()->SetTitle("#Delta#theta_{#gamma} [deg]");
  h_PhotRes2D_theta->GetYaxis()->SetRangeUser(-60,60);
  h_PhotRes2D_theta->GetZaxis()->SetMaxDigits(2);
  h_PhotRes2D_theta->Draw("colz");
  gPad->Update();
  auto palette = (TPaletteAxis*)h_PhotRes2D_theta->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.88);
  palette->SetY1NDC(0.12);
  palette->SetX2NDC(0.92);
  palette->SetY2NDC(0.95);
  
  // NOTE: WILL NEED TO MANUALLY RESIZE TPAD FOR 2D HISTOGRAM

  cPhotRes->SaveAs("TDR_" + energy +"_PhotRes.png");

  // BACKUP HISTOGRAMS - PROJECTION Y OF 2D PLOT @ 90 AND 170-175 DEGREES
  // Projection in centre of barrel
  TH1D* hProj_Cent = h_PhotRes2D_theta->ProjectionY("proj_cent",h_PhotRes2D_theta->GetXaxis()->FindBin(89.9),h_PhotRes2D_theta->GetXaxis()->FindBin(90.1));
  // Projection in backward direction (where most photons go)
  TH1D* hProj_FBwd = h_PhotRes2D_theta->ProjectionY("proj_fbwd",h_PhotRes2D_theta->GetXaxis()->FindBin(170),h_PhotRes2D_theta->GetXaxis()->FindBin(175));

  TCanvas* cPhotRes_a = new TCanvas("cPhotRes_a","",1500,750);
  cPhotRes_a->Divide(2,1);
  cPhotRes_a->cd(1);
  hProj_Cent->SetLineColor(kBlack);
  hProj_Cent->SetMarkerColor(kBlack);
  hProj_Cent->SetMarkerStyle(20);
  hProj_Cent->GetXaxis()->SetTitle("#Delta#theta_{#gamma} [deg]");
  hProj_Cent->GetXaxis()->SetRangeUser(-10,10);
  hProj_Cent->GetYaxis()->SetTitle("Counts / 0.5 deg");
  hProj_Cent->GetYaxis()->SetTitleOffset(1.45);
  hProj_Cent->GetYaxis()->SetRangeUser(0.1,1.5*hProj_Cent->GetMaximum());
  hProj_Cent->Draw("ep");
  TLatex* tePICLabel_photres_a = new TLatex(-8.8, 0.9*hProj_Cent->GetMaximum(), "#splitline{#bf{ePIC} Performance " + campaign + ", " + energy + " GeV}{ep #rightarrow e'p'#gamma, all reco.#gamma (#theta_{#gamma}: [89.5, 90.5] deg)}");
  tePICLabel_photres_a->SetTextSize(0.04);
  tePICLabel_photres_a->Draw("same");
  TFitResultPtr fitres_a = hProj_Cent->Fit("gaus","S0","",-10,10);
  TString sFitRes_a = "#sigma = " + combineScientific(fitres_a->Parameter(2),fitres_a->ParError(2)) + " deg";
  TLatex* tFitRes_a = new TLatex(-8.8, 0.75*hProj_Cent->GetMaximum(), sFitRes_a);
  tFitRes_a->Draw("same");
  
  cPhotRes_a->cd(2);
  gPad->SetLogz();
  h_PhotRes2D_theta->Draw("colz");
  gPad->Update();
  cPhotRes_a->SaveAs("TDR_" + energy +"_PhotRes_Cent.png");

  TCanvas* cPhotRes_b = new TCanvas("cPhotRes_b","",1500,750);
  cPhotRes_b->Divide(2,1);
  cPhotRes_b->cd(1);
  hProj_FBwd->SetLineColor(kBlack);
  hProj_FBwd->SetMarkerColor(kBlack);
  hProj_FBwd->SetMarkerStyle(20);
  hProj_FBwd->GetXaxis()->SetTitle("#Delta#theta_{#gamma} [deg]");
  hProj_FBwd->GetXaxis()->SetRangeUser(-10,10);
  hProj_FBwd->GetYaxis()->SetTitle("Counts / 0.5 deg");
  hProj_FBwd->GetYaxis()->SetTitleOffset(1.45);
  hProj_FBwd->GetYaxis()->SetRangeUser(0.1,1.5*hProj_FBwd->GetMaximum());
  hProj_FBwd->Draw("ep");
  TLatex* tePICLabel_photres_b = new TLatex(-8.8, 0.9*hProj_FBwd->GetMaximum(), "#splitline{#bf{ePIC} Performance " + campaign + ", " + energy + " GeV}{ep #rightarrow e'p'#gamma, all reco.#gamma (#theta_{#gamma}: [169.5, 175.5] deg)}");
  tePICLabel_photres_b->SetTextSize(0.04);
  tePICLabel_photres_b->Draw("same");
  TFitResultPtr fitres_b = hProj_FBwd->Fit("gaus","S0","",-10,10);
  TString sFitRes_b = "#sigma = " + combineScientific(fitres_b->Parameter(2),fitres_b->ParError(2)) + " deg";
  TLatex* tFitRes_b = new TLatex(-8.8, 0.75*hProj_FBwd->GetMaximum(), sFitRes_b);
  tFitRes_b->Draw("same");
  
  cPhotRes_b->cd(2);
  gPad->SetLogz();
  h_PhotRes2D_theta->Draw("colz");
  gPad->Update();
  cPhotRes_b->SaveAs("TDR_" + energy +"_PhotRes_FBwd.png");

  // CANVAS: t-distribution
  // DUPLICATE DRAWN HISTOGRAMS FOR SCALING
  // MC truth
  TH1D* h_t_TruthScaled  = (TH1D*)h_t_Truth->Clone("t_truthscaled");
  // Reco - B0 raw, RP raw and "good"
  TH1D* h_t_B0RecoScaled = (TH1D*)h_t_B0Reco->Clone("t_b0recoscaled");
  TH1D* h_t_RPRecoScaled = (TH1D*)h_t_RPReco->Clone("t_rprecoscaled");
  TH1D* h_t_GoodRecoScaled = (TH1D*)h_t_GoodReco->Clone("t_goodrecoscaled");
  // Corrected - With and without fit points
  TH1D* h_t_CorrFitScaled = (TH1D*)h_t_CorrFit->Clone("t_corrfitscaled");
  TH1D* h_t_CorrDrawScaled = (TH1D*)h_t_CorrDraw->Clone("t_corrdrawscaled");

  // APPLY SCALINGS
  h_t_TruthScaled->Scale(scaleTo5);
  h_t_B0RecoScaled->Scale(scaleTo5);
  h_t_RPRecoScaled->Scale(scaleTo5);
  h_t_GoodRecoScaled->Scale(scaleTo5);
  h_t_CorrFitScaled->Scale(scaleTo5);
  h_t_CorrDrawScaled->Scale(scaleTo5);
  
  TCanvas* ctDist = new TCanvas("ctDist","",1200,800);
  gPad->SetLogy();

  // Drawing options - Markers and lines
  h_t_TruthScaled->SetLineColor(kBlack); // TO BE DRAW AS BARS
  // Reco - blue
  // Open circles - all reco.; closed circles - used in fit calculation
  h_t_GoodRecoScaled->SetLineColor(kP6Blue);
  h_t_GoodRecoScaled->SetMarkerColor(kP6Blue);
  h_t_GoodRecoScaled->SetMarkerStyle(20);
  h_t_GoodRecoScaled->SetMarkerSize(1.5);
  h_t_B0RecoScaled->SetLineColor(kP6Blue);
  h_t_B0RecoScaled->SetMarkerColor(kP6Blue);
  h_t_B0RecoScaled->SetMarkerStyle(24);
  h_t_B0RecoScaled->SetMarkerSize(1.5);
  h_t_RPRecoScaled->SetLineColor(kP6Blue);
  h_t_RPRecoScaled->SetMarkerColor(kP6Blue);
  h_t_RPRecoScaled->SetMarkerStyle(24);
  h_t_RPRecoScaled->SetMarkerSize(1.5);
  // Fitted - filled diamonds
  // Yellow - full (with fit); red - corr. reco. only
  h_t_CorrFitScaled->SetLineColor(kP6Yellow);
  h_t_CorrFitScaled->SetMarkerColor(kP6Yellow);
  h_t_CorrFitScaled->SetMarkerStyle(33);
  h_t_CorrFitScaled->SetMarkerSize(2.5);
  h_t_CorrDrawScaled->SetLineColor(kP6Red);
  h_t_CorrDrawScaled->SetMarkerColor(kP6Red);
  h_t_CorrDrawScaled->SetMarkerStyle(33);
  h_t_CorrDrawScaled->SetMarkerSize(2.5);
  // Axes
  h_t_TruthScaled->GetXaxis()->SetTitle("|t| [GeV^{2}]");
  h_t_TruthScaled->GetYaxis()->SetTitle("Counts / 0.1 GeV^{2}");
  h_t_TruthScaled->GetYaxis()->SetRangeUser(0.5, 10*h_t_TruthScaled->GetMaximum());
  // Draw 
  h_t_TruthScaled->Draw("HIST");
  h_t_B0RecoScaled->Draw("ep same");
  h_t_RPRecoScaled->Draw("ep same");
  h_t_GoodRecoScaled->Draw("ep same"); // Draw ON TOP of open points
  h_t_CorrFitScaled->Draw("ep same");
  h_t_CorrDrawScaled->Draw("ep same"); // Draw ON TOP of yellow points
  // Add legend
  TLegend* lctDist = new TLegend(0.14,0.21,0.38,0.51);
  lctDist->AddEntry(h_t_TruthScaled, "#bf{EpIC} MC", "l");
  lctDist->AddEntry(h_t_B0RecoScaled, "Raw reco.", "pl");
  lctDist->AddEntry(h_t_GoodRecoScaled, "Raw reco. (#varepsilon #geq 5%)", "pl");
  lctDist->AddEntry(h_t_CorrDrawScaled, "Corrected reco.", "pl");
  lctDist->AddEntry(h_t_CorrFitScaled, "Points from fit", "pl");
  lctDist->Draw();
  // Add text
  TLatex* tePICLabel_t = new TLatex(0.42, 0.87, "#splitline{#bf{ePIC} Performance " + campaign + ", 10x100 GeV}{ep #rightarrow e'p'#gamma, L_{proj} = 5 fb^{-1}}");
  tePICLabel_t->SetNDC();
  tePICLabel_t->Draw("same");
  TLatex* tLabel1_t = new TLatex(0.42, 0.73, "#splitline{Full final state, Q^{2} #geq 1 GeV^{2}, M^{2}_{miss} #leq 1 GeV^{2}}{Expo. fit to reco. where #varepsilon #geq 5%}");
  tLabel1_t->SetNDC();
  tLabel1_t->Draw("same");

  ctDist->SaveAs("TDR_" + energy +"_tdist.png");

  // CANVAS: t-resolution
  TCanvas* ctRes = new TCanvas("ctRes","",1500,500);
  ctRes->Divide(3,1,-1);
  Float_t tRes_common_ystart{1e-3}; 
  Float_t tRes_common_yend{3.5e-1};
  Float_t tRes_common_xstart{-0.01}; 
  Float_t tRes_common_xend{1.9};

  // Pad 1 - 5x41
  ctRes->cd(1);
  //gPad->SetLogy();
  // Markers and lines
  h_tresB0_5x41->SetMarkerStyle(20);
  h_tresB0_5x41->SetMarkerSize(1.5);
  h_tresB0_5x41->SetMarkerColor(kP6Red);
  h_tresB0_5x41->SetLineColor(kP6Red);
  h_tresRP_5x41->SetMarkerStyle(24);
  h_tresRP_5x41->SetMarkerSize(1.5);
  h_tresRP_5x41->SetMarkerColor(kP6Red);
  h_tresRP_5x41->SetLineColor(kP6Red);
  // Axes
  h_tresB0_5x41->SetTitle("");
  h_tresB0_5x41->GetYaxis()->SetTitle("RMS(#Deltat/|t|_{MC})");
  h_tresB0_5x41->GetYaxis()->SetRangeUser(tRes_common_ystart,tRes_common_yend);
  h_tresB0_5x41->GetXaxis()->SetRangeUser(tRes_common_xstart,tRes_common_xend);
  h_tresB0_5x41->GetYaxis()->SetNdivisions(505);
  h_tresB0_5x41->GetXaxis()->SetTitle("|t|_{MC} [GeV^{2}]");
  // Draw
  h_tresB0_5x41->Draw();
  h_tresRP_5x41->Draw("same");
  // Add text
  TLatex* tePICLabel_tres = new TLatex(0.17, 0.91, "#splitline{#bf{ePIC} Performance " + campaign + "}{ep #rightarrow e'p'#gamma, L_{proj} = 5 fb^{-1}}");
  tePICLabel_tres->SetNDC();
  tePICLabel_tres->Draw("same");
  // Add legend
  TLegend* ltres_p1 = new TLegend(0.46,0.19,0.98,0.36);
  ltres_p1->SetHeader("-t = (p' - p)^{2} resolution, 5x41 GeV");
  ltres_p1->AddEntry((TObject*)0, "", "");
  ltres_p1->AddEntry(h_tresB0_5x41,"p' in B0","pl");
  ltres_p1->AddEntry(h_tresRP_5x41,"p' in RP","pl");
  ltres_p1->Draw();
  

  // Pad 2 - 10x100
  ctRes->cd(2);
  //gPad->SetLogy();
  // Markers and lines
  h_tresB0_10x100->SetMarkerStyle(21);
  h_tresB0_10x100->SetMarkerSize(1.3);
  h_tresB0_10x100->SetMarkerColor(kBlue);
  h_tresB0_10x100->SetLineColor(kBlue);
  h_tresRP_10x100->SetMarkerStyle(25);
  h_tresRP_10x100->SetMarkerSize(1.3);
  h_tresRP_10x100->SetMarkerColor(kBlue);
  h_tresRP_10x100->SetLineColor(kBlue);
  // Axes
  h_tresB0_10x100->SetTitle("");
  h_tresB0_10x100->GetYaxis()->SetRangeUser(tRes_common_ystart,tRes_common_yend);
  h_tresB0_10x100->GetYaxis()->SetNdivisions(505);
  h_tresB0_10x100->GetXaxis()->SetRangeUser(tRes_common_xstart,tRes_common_xend);
  h_tresB0_10x100->GetXaxis()->SetTitle("|t|_{MC} [GeV^{2}]");
  // Draw
  h_tresB0_10x100->Draw();
  h_tresRP_10x100->Draw("same");
  // Add legend
  TLegend* ltres_p2 = new TLegend(0.05,0.78,0.75,0.95);
  ltres_p2->SetHeader("10x100 GeV");
  ltres_p2->AddEntry((TObject*)0, "", "");
  ltres_p2->AddEntry(h_tresB0_10x100,"p' in B0","pl");
  ltres_p2->AddEntry(h_tresRP_10x100,"p' in RP","pl");
  ltres_p2->Draw();


  // Pad 3 - 18x275
  ctRes->cd(3);
  //gPad->SetLogy();
  // Markers and lines
  h_tresB0_18x275->SetMarkerStyle(22);
  h_tresB0_18x275->SetMarkerSize(1.4);
  h_tresB0_18x275->SetMarkerColor(kP6Grape);
  h_tresB0_18x275->SetLineColor(kP6Grape);
  h_tresRP_18x275->SetMarkerStyle(26);
  h_tresRP_18x275->SetMarkerSize(1.4);
  h_tresRP_18x275->SetMarkerColor(kP6Grape);
  h_tresRP_18x275->SetLineColor(kP6Grape);
  // Axes
  h_tresB0_18x275->SetTitle("");
  h_tresB0_18x275->GetYaxis()->SetRangeUser(tRes_common_ystart,tRes_common_yend);
  h_tresB0_18x275->GetYaxis()->SetNdivisions(505);
  h_tresB0_18x275->GetXaxis()->SetRangeUser(tRes_common_xstart,tRes_common_xend);
  h_tresB0_18x275->GetXaxis()->SetTitle("|t|_{MC} [GeV^{2}]");
  // Draw
  h_tresB0_18x275->Draw();
  h_tresRP_18x275->Draw("same");
  // Add legend
  TLegend* ltres_p3 = new TLegend(0.05,0.78,0.75,0.95);
  ltres_p3->SetHeader("18x275 GeV");
  ltres_p3->AddEntry((TObject*)0, "", "");
  ltres_p3->AddEntry(h_tresB0_18x275,"Proton in B0","pl");
  ltres_p3->AddEntry(h_tresRP_18x275,"Proton in RP","pl");
  ltres_p3->Draw();
    
  ctRes->SaveAs("TDR_tres.png");

  // Save histograms in ROOT file for potential use later
  TFile* fOut = new TFile("./hists_TDR.root","RECREATE");
  fOut->cd();
  h_eta_MCe->Write();
  h_eta_RPe->Write();
  h_eta_MCg->Write();
  h_eta_RPg->Write();
  h_eta_MCp->Write();
  h_eta_RPp->Write();
  h_eta_RPPp->Write();
  cEta->Write();
  h_PhotRes_theta->Write();
  h_PhotRes2D_theta->Write();
  cPhotRes->Write();
  hProj_Cent->Write();
  cPhotRes_a->Write();
  hProj_FBwd->Write();
  cPhotRes_b->Write();
  h_t_TruthScaled->Write();
  h_t_B0RecoScaled->Write();
  h_t_RPRecoScaled->Write();
  h_t_GoodRecoScaled->Write();
  h_t_CorrFitScaled->Write();
  h_t_CorrDrawScaled->Write();
  ctDist->Write();
  h_tresB0_5x41->Write();
  h_tresB0_10x100->Write();
  h_tresB0_18x275->Write();
  h_tresRP_5x41->Write();
  h_tresRP_10x100->Write();
  h_tresRP_18x275->Write();
  ctRes->Write();

  // Close files
  f5x41->Close();
  f10x100->Close();
  f18x275->Close();
  fOut->Close();

  return;
}
