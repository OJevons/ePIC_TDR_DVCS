using namespace std;

#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include "../plotting/ePIC_style.C"

bool kSAVE = false;

double calcScalingDVCS(TString energy, TString hel, float lumi){
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

double calcScalingDIS(TString energy, TString q2range, int n_gen, float lumi){
  // Holding variables - No. of events generated, integrated cross-section
  // These vary by beam settings
  Double_t fXSint{0}

  // Cross-sections taken from pythia8 samples
  if(energy == "9x130"){
    if(q2range == "lo") fXSint = 5.29241e-7;
    else if(q2range == "hi") fXSint = 4.06442e-8;
    else fXSint = 1;
  }
  else if(energy == "9x275"){
    if(q2range == "lo") fXSint = 6.00707e-7;
    else if(q2range == "hi") fXSint = 5.24191e-8;
    else fXSint = 1;
  }

  Double_t genlumi = (float)n_gen/fXSint;
  Double_t scale = lumi/genlumi;

  return scale;
}

//---------------------------------------------------------------------
// MAIN
// 
// DIS physics background plots to DVCS analysis
//---------------------------------------------------------------------
void Plots_DISBkg(TString campaign = "26.07.1", TString energy = "9x130", TString hel = "emhTm"){
  // Print plot settings
  cout<<"\n--------------------------------------"<<endl;
  cout<<"Processing DIS background plots"<<endl;
  cout<<"\tCampaign: "<<campaign<<endl;
  cout<<"\tBeam energy combination: "<<energy<<endl;
  cout<<"\tPolarization state: "<<hel<<endl;
  cout<<"--------------------------------------\n"<<endl;

  // Set beam energies
  Float_t fEBeam{0}, fPBeam{0};
  if(energy != "9x130" && energy != "9x275"){
    cout<<"Invalid beam energy for early science!"<<endl;
    return;
  }
  else{
    fEBeam = 9.0;

    if(energy == "9x130") fPBeam = 130.;
    else if(energy == "9x275") fPBeam = 275.;
  }
  
  //--------------------------------------------------------------------
  // Load histograms from file - DVCS baseline
  //--------------------------------------------------------------------
  TString sDVCS = "../rootfiles/ePIC_DVCS_"+campaign+"_"+energy+"_"+hel+"_simple.root";
  TFile* fDVCS = TFile::Open(sDVCS);
  // Mandelstam t distributions
  // MC truth
  TH1D* h_t_Truth = (TH1D*)fDVCS->Get("t_truth");
  // BABE method (using recoil proton)
  TH1D* h_t_B0Acc = (TH1D*)fDVCS->Get("t_b0acc");
  TH1D* h_t_B0Rec = (TH1D*)fDVCS->Get("t_b0reco");
  TH1D* h_t_RPAcc = (TH1D*)fDVCS->Get("t_rpacc");
  TH1D* h_t_RPRec = (TH1D*)fDVCS->Get("t_rpreco");
  // Light-cone semi-inclusive method
  TH1D* h_t_LCAcc = (TH1D*)fDVCS->Get("t_lcacc");
  TH1D* h_t_LCRec = (TH1D*)fDVCS->Get("t_lcreco");

  //--------------------------------------------------------------------
  // Load histograms from file - DIS comparison
  // "Low Q2" -> minQ2=1
  // "High Q2" -> minQ2=10
  //--------------------------------------------------------------------
  TString sDISLQ2 = "../rootfiles/ePIC_DIS_"+campaign+"_"+energy+"_minQ2=1.root";
  TFile* fDISLQ2 = TFile::Open(sDISLQ2);
  TH1D* h_tDISLQ2_B0Rec = (TH1D*)fDISLQ2->Get("t_b0reco");
  TH1D* h_tDISLQ2_RPRec = (TH1D*)fDISLQ2->Get("t_rpreco");
  TH1D* h_tDISLQ2_LCRec = (TH1D*)fDISLQ2->Get("t_lcreco");
  TString sDISHQ2 = "../rootfiles/ePIC_DIS_"+campaign+"_"+energy+"_minQ2=10.root";
  TFile* fDISHQ2 = TFile::Open(sDISLQ2);
  TH1D* h_tDISHQ2_B0Rec = (TH1D*)fDISHQ2->Get("t_b0reco");
  TH1D* h_tDISHQ2_RPRec = (TH1D*)fDISHQ2->Get("t_rpreco");
  TH1D* h_tDISHQ2_LCRec = (TH1D*)fDISHQ2->Get("t_lcreco");
 
  //---------------------------------------------------------------------
  // Calculations on histograms
  // Scaling to 5fb-1
  //---------------------------------------------------------------------
  // 1. Scaling factor to EIC lumi
  // -> SET LUMINOSITY TO SCALE TO
  double EIClumi{1.};
  if(energy == "9x130") EIClumi = 1e15;
  else if(energy == "9x275") EIClumi = 2.5e15;
  else EIClumi = 2.5e15;
  
  // Calculations - SCALING FACTOR TO EIC luminosities
  // -> DIS: NEED TO KNOW NO. OF GENERATED EVENTS
  scale_DVCS = calcScalingDVCS(energy, hel, EIClumi);
  int n_gen_lo = (TH1D*)fDISLQ2->Get("count")->GetEntries();
  scale_DISLo = calcScalingDIS(energy, "lo", n_gen_lo, EIClumi);
  int n_gen_hi = (TH1D*)fDISHQ2->Get("count")->GetEntries();
  scale_DISHi = calcScalingDIS(energy, "hi", n_gen_hi, EIClumi);
  
  //cout<<"DVCS data represents "<<lumi/1e15<<" fb-1\n\tDIS low Q2 = "<<lumi_DISlo/1e15<<" fb-1\n\tDIS high Q2 = "<<lumi_DIShi/1e15<<" fb-1\n"<<endl;

  //---------------------------------------------------------------------
  // Detector efficiency corrections
  //
  // Calculate efficiency: MC accepted/MC truth
  // 1. Duplicate MCA histograms
  // 2. Divide by MC truth
  // 3. Duplicate reco. histograms
  // 4. Divide reco. by efficiency
  //---------------------------------------------------------------------
  // Mandelstam t - for DVCS baseline
  TH1D* h_t_B0Eff = (TH1D*)h_t_B0Acc->Clone("t_b0eff");
  TH1D* h_t_RPEff = (TH1D*)h_t_RPAcc->Clone("t_rpeff");
  TH1D* h_t_LCEff = (TH1D*)h_t_LCAcc->Clone("t_lceff");
  h_t_B0Eff->Divide(h_t_Truth);
  h_t_RPEff->Divide(h_t_Truth);
  h_t_LCEff->Divide(h_t_Truth);
  TH1D* h_t_B0Corr = (TH1D*)h_t_B0Rec->Clone("t_b0corr");
  TH1D* h_t_RPCorr = (TH1D*)h_t_RPRec->Clone("t_rpcorr");
  TH1D* h_t_LCCorr = (TH1D*)h_t_LCRec->Clone("t_lccorr");
  h_t_B0Corr->Divide(h_t_B0Eff);
  h_t_RPCorr->Divide(h_t_RPEff);
  h_t_LCCorr->Divide(h_t_LCEff);

  //---------------------------------------------------------------------
  // Mandelstam t: combine corrected points
  //---------------------------------------------------------------------
  TH1D* h_t_GoodCorr = (TH1D*)h_t_Truth->Clone("t_goodcorr");
  h_t_GoodCorr->Reset();
  TH1D* h_t_GoodReco = (TH1D*)h_t_Truth->Clone("t_goodreco"); // FOR PLOTTING PURPOSES LATER
  h_t_GoodReco->Reset(); 

  // Combine reconstructed into "good" histograms
  // If proton methods have efficiency > 5%, use those (if both B0 and RP, use B0)
  // If not, use light-cone e'gamma method
  for(int bin{1}; bin<h_t_Truth->GetNbinsX()+1; ++bin){
    // 1. Check RP
    if(h_t_RPEff->GetBinContent(bin) >= 0.05){
      h_t_GoodCorr->SetBinContent(bin,h_t_RPCorr->GetBinContent(bin));
      h_t_GoodCorr->SetBinError(bin,h_t_RPCorr->GetBinError(bin));
      h_t_GoodReco->SetBinContent(bin,h_t_RPRec->GetBinContent(bin));
      h_t_GoodReco->SetBinError(bin,h_t_RPRec->GetBinError(bin));
    }
    // 2. Then check B0 (will override RP in case of overlap)
    if(h_t_B0Eff->GetBinContent(bin) >= 0.05){
      h_t_GoodCorr->SetBinContent(bin,h_t_B0Corr->GetBinContent(bin));
      h_t_GoodCorr->SetBinError(bin,h_t_B0Corr->GetBinError(bin));
      h_t_GoodReco->SetBinContent(bin,h_t_B0Rec->GetBinContent(bin));
      h_t_GoodReco->SetBinError(bin,h_t_B0Rec->GetBinError(bin));
    }
    // 3. If neither proton reconstruction method has good reconstruction, use e' and gamma
    else{
      h_t_GoodCorr->SetBinContent(bin,h_t_LCCorr->GetBinContent(bin));
      h_t_GoodCorr->SetBinError(bin,h_t_LCCorr->GetBinError(bin));
      h_t_GoodReco->SetBinContent(bin,h_t_LCRec->GetBinContent(bin));
      h_t_GoodReco->SetBinError(bin,h_t_LCRec->GetBinError(bin));
    }
  }

  //--------------------------------------------------------------------------------------
  // Draw histograms
  //--------------------------------------------------------------------------------------
  // Global style
  gROOT->ProcessLine("set_ePIC_style()");
  gStyle->SetCanvasPreferGL(kTRUE);
  
  // Draw options and scaling
  h_t_Truth->Scale(scale_DVCS);
  h_t_GoodCorr->Scale(scale_DVCS);
  h_t_GoodReco->Scale(scale_DVCS);
  h_tDISLQ2_B0Rec->Scale(scale_DISLo);
  h_tDISLQ2_RPRec->Scale(scale_DISLo);
  h_tDISHQ2_B0Rec->Scale(scale_DISHi);
  h_tDISHQ2_RPRec->Scale(scale_DISHi);
  
  h_t_B0Rec->Scale(scaleTo5);
  h_t_B0Rec->SetLineColor(kP6Blue);
  h_t_B0Rec->SetLineWidth(2);
  h_t_B0Rec->SetMarkerColor(kP6Blue);
  h_t_B0Rec->SetMarkerStyle(24);
  h_t_B0Rec->SetMarkerSize(2);

  h_t_RPRec->Scale(scaleTo5);
  h_t_RPRec->SetLineColor(kP6Blue);
  h_t_RPRec->SetLineWidth(2);
  h_t_RPRec->SetMarkerColor(kP6Blue);
  h_t_RPRec->SetMarkerStyle(25);
  h_t_RPRec->SetMarkerSize(2);

  h_t_LCRec->Scale(scaleTo5);
  h_t_LCRec->SetLineColor(kP6Blue);
  h_t_LCRec->SetLineWidth(2);
  h_t_LCRec->SetMarkerColor(kP6Blue);
  h_t_LCRec->SetMarkerStyle(26);
  h_t_LCRec->SetMarkerSize(2);
  
  // Set draw options - markers and Lines
  h_t_Truth->SetLineColor(kBlack);
  h_t_Truth->SetLineWidth(2);
  h_t_GoodCorr->SetLineColor(kP6Blue);
  h_t_GoodCorr->SetLineWidth(2);
  h_t_GoodCorr->SetMarkerColor(kP6Blue);
  h_t_GoodCorr->SetMarkerStyle(20);
  h_t_GoodCorr->SetMarkerSize(2);
  h_tDISLQ2_B0Rec->SetLineColor(kP6Yellow);
  h_tDISLQ2_B0Rec->SetLineWidth(2);
  h_tDISLQ2_RPRec->SetLineColor(kP6Yellow);
  h_tDISLQ2_RPRec->SetLineWidth(2);
  h_tDISHQ2_B0Rec->SetLineColor(kP6Red);
  h_tDISHQ2_B0Rec->SetLineWidth(2);
  h_tDISHQ2_RPRec->SetLineColor(kP6Red);
  h_tDISHQ2_RPRec->SetLineWidth(2);
  // Set draw options - axes
  h_t_Truth->GetYaxis()->SetRangeUser(1,500*h_t_Truth->GetMaximum());
  h_t_Truth->GetYaxis()->SetTitle("Counts/0.1 GeV^{2}");
  h_t_Truth->GetXaxis()->SetTitle("|t| [GeV^{2}]");


  // CANVAS - DVCS vs DIS low Q2 ONLY
  TCanvas* ct_DISlo = new TCanvas("ct_dislo","",1000,1000);
  gPad->SetLogy();
  // Draw
  h_t_Truth->Draw("hist");
  h_t_GoodCorr->Draw("same");
  h_tDISLQ2_B0Rec->Draw("ep same");
  h_tDISLQ2_RPRec->Draw("ep same");

  h_t_B0Rec->Draw("same");
  h_t_RPRec->Draw("same");
  h_t_LCRec->Draw("same");  
  // Add text
  TLatex* tePICLabel_tDIS = new TLatex(0.1, 0.18*h_t_Truth->GetMaximum(), "#splitline{#bf{ePIC} Performance " + campaign + ", " + energy + " GeV}{ep, 10x130 GeV , L_{proj} = 5 fb^{-1}}");
  tePICLabel_tDIS->SetTextSize(0.037);
  tePICLabel_tDIS->Draw("same");
  // Add legend
  TLegend* lt_DISlo = new TLegend(0.46,0.61,0.85,0.79);
  lt_DISlo->SetLineWidth(0);
  lt_DISlo->SetFillStyle(0);
  lt_DISlo->SetTextSize(0.037);
  lt_DISlo->SetHeader("Comparing to low Q^{2} DIS");
  lt_DISlo->AddEntry(h_t_Truth,"EpIC MC truth","l");
  lt_DISlo->AddEntry(h_t_GoodCorr,"Corrected reco. DVCS","lp");
  lt_DISlo->AddEntry(h_tDISLQ2_B0Rec,"Raw reco. DIS","lp");
  lt_DISlo->Draw();
  // Save figure
  if(kSAVE) ct_DISlo->SaveAs("../figs/DISComp_" + energy +"_t_Q2lo.png");
  
  
  // CANVAS - DVCS vs DIS high Q2 ONLY
  TCanvas* ct_DIShi = new TCanvas("ct_dishi","",1000,1000);
  gPad->SetLogy();
  // Draw
  h_t_Truth->Draw("hist");
  h_t_GoodCorr->Draw("same");
  h_tDISHQ2_B0Rec->Draw("e same");
  h_tDISHQ2_RPRec->Draw("e same");

  h_t_B0Rec->Draw("same");
  h_t_RPRec->Draw("same");
  h_t_LCRec->Draw("same");  
  // Add text
  tePICLabel_tDIS->Draw("same");
  // Add legend
  TLegend* lt_DIShi = new TLegend(0.46,0.61,0.85,0.79);
  lt_DIShi->SetLineWidth(0);
  lt_DIShi->SetFillStyle(0);
  lt_DIShi->SetTextSize(0.037);
  lt_DIShi->SetHeader("Comparing to high Q^{2} DIS");
  lt_DIShi->AddEntry(h_t_Truth,"EpIC MC truth","l");
  lt_DIShi->AddEntry(h_t_GoodCorr,"Corrected reco. DVCS","lp");
  lt_DIShi->AddEntry(h_tDISHQ2_B0Rec,"Raw reco. DIS","lp");
  lt_DIShi->Draw();
  // Save figure
  if(kSAVE) ct_DIShi->SaveAs("../figs/DISComp_" + energy +"_t_Q2hi.png");
  
  
  // CANVAS - DVCS vs DIS both energies
  TCanvas* ct_DISall = new TCanvas("ct_disall","",1000,1000);
  gPad->SetLogy();
  // Draw
  h_t_Truth->Draw("hist");
  h_t_GoodCorr->Draw("same");
  h_tDISLQ2_B0Rec->Draw("e same");
  h_tDISLQ2_RPRec->Draw("e same");
  h_tDISHQ2_B0Rec->Draw("e same");
  h_tDISHQ2_RPRec->Draw("e same");

  h_t_B0Rec->Draw("same");
  h_t_RPRec->Draw("same");
  h_t_LCRec->Draw("same");  
  // Add text
  tePICLabel_tDIS->Draw("same");
  // Add legend
  TLegend* lt_DISall = new TLegend(0.41,0.62,0.80,0.80);
  lt_DISall->SetLineWidth(0);
  lt_DISall->SetFillStyle(0);
  lt_DISall->SetTextSize(0.037);
  lt_DISall->SetHeader("DIS - BABE only");
  lt_DISall->AddEntry(h_t_Truth,"EpIC MC truth","l");
  lt_DISall->AddEntry(h_t_GoodCorr,"Reco. DVCS (open: raw.; closed: corr.)","lp");
  lt_DISall->AddEntry(h_tDISLQ2_B0Rec,"Raw reco. DIS - Q^{2} < 1 GeV^{2}","lp");
  lt_DISall->AddEntry(h_tDISHQ2_B0Rec,"Raw reco. DIS - Q^{2} < 10 GeV^{2}","lp");
  lt_DISall->Draw();
  // Save figure
  if(kSAVE) ct_DISall->SaveAs("../figs/DISComp_" + energy +"_t_Q2hi.png");
  

  // CANVAS - eXBE calculation
  h_tDISLQ2_LCRec->Scale(scaleTo5_DISlo);
  h_tDISHQ2_LCRec->Scale(scaleTo5_DIShi);
  
  h_tDISLQ2_LCRec->SetLineColor(kP6Grape);
  h_tDISLQ2_LCRec->SetLineWidth(2);
  h_tDISHQ2_LCRec->SetLineColor(kP6Gray);
  h_tDISHQ2_LCRec->SetLineWidth(2);

  TCanvas* ct_DISeg = new TCanvas("ct_diseg","",1000,1000);
  gPad->SetLogy();
  h_t_Truth->Draw("hist");
  h_t_GoodCorr->Draw("same");
  h_tDISLQ2_LCRec->Draw("e same");
  h_tDISHQ2_LCRec->Draw("e same");

  h_t_B0Rec->Draw("same");
  h_t_RPRec->Draw("same");
  h_t_LCRec->Draw("same");  
  // Add text
  tePICLabel_tDIS->Draw("same");
  // Add legend
  TLegend* lt_DISeg = new TLegend(0.15,0.16,0.55,0.35);
  lt_DISeg->SetLineWidth(0);
  lt_DISeg->SetFillStyle(0);
  lt_DISeg->SetTextSize(0.037);
  lt_DISeg->SetHeader("DIS - eXBE only");
  lt_DISeg->AddEntry(h_t_Truth,"EpIC MC truth","l");
  lt_DISeg->AddEntry(h_t_GoodCorr,"Reco. DVCS (open: raw.; closed: corr.)","lp");
  lt_DISeg->AddEntry(h_tDISLQ2_LCRec,"Raw reco. DIS - Q^{2} < 1 GeV^{2}","lp");
  lt_DISeg->AddEntry(h_tDISHQ2_LCRec,"Raw reco. DIS - Q^{2} < 10 GeV^{2}","lp");
  lt_DISeg->Draw();
  // Save figure
  if(kSAVE) ct_DISeg->SaveAs("../figs/DISComp_" + energy +"_t_Q2hi.png");


  // CANVAS - Ratios DIS/DVCS
  // Sum DIS histograms
  TH1D* h_tDIS_AllBABE = (TH1D*)h_tDISLQ2_B0Rec->Clone("tdis_allbabe");
  h_tDIS_AllBABE->Add(h_tDISLQ2_RPRec);
  h_tDIS_AllBABE->Add(h_tDISHQ2_B0Rec);
  h_tDIS_AllBABE->Add(h_tDISHQ2_RPRec);
  TH1D* h_tDIS_AllRec = (TH1D*)h_tDIS_AllBABE->Clone("tdis_allrec");
  h_tDIS_AllRec->Add(h_tDISLQ2_LCRec);
  h_tDIS_AllRec->Add(h_tDISHQ2_LCRec);
  h_tDIS_AllRec->SetLineColor(kP6Red);

  cout<<"\n------------- Ratios - t-integrated -------------"
      <<"\nDIS, BABE only:"
      <<"\n\t DIS/DVCS corrected = "<<h_tDIS_AllBABE->Integral()/h_t_GoodCorr->Integral()
      <<"\n\t DIS/DVCS raw (RP only) = "<<h_tDIS_AllBABE->Integral()/h_t_RPRec->Integral()
      <<"\n\t DIS/DVCS raw (B0 only) = "<<h_tDIS_AllBABE->Integral()/h_t_B0Rec->Integral()
      <<"\n\t DIS/DVCS raw (eXBE only) = "<<h_tDIS_AllBABE->Integral()/h_t_LCRec->Integral()
      <<"\nDIS, all reco.:"
      <<"\n\t DIS/DVCS corrected = "<<h_tDIS_AllRec->Integral()/h_t_GoodCorr->Integral()
      <<"\n\t DIS/DVCS raw (RP only) = "<<h_tDIS_AllRec->Integral()/h_t_RPRec->Integral()
      <<"\n\t DIS/DVCS raw (B0 only) = "<<h_tDIS_AllRec->Integral()/h_t_B0Rec->Integral()
      <<"\n\t DIS/DVCS raw (eXBE only) = "<<h_tDIS_AllRec->Integral()/h_t_LCRec->Integral()<<endl;

  TCanvas* ct_Ratio = new TCanvas("ct_ratio","",1800,900);
  // Split into pads - 1 left, 4 right
  TPad* pLeft = new TPad("left_pad","",0.025,0.025,0.475,0.975);
  TPad* pRight1 = new TPad("right_pad1","",0.500,0.025,0.975,0.263);
  pRight1->SetTopMargin(0);
  pRight1->SetBottomMargin(0.3);
  TPad* pRight2 = new TPad("right_pad2","",0.500,0.263,0.975,0.501);
  pRight2->SetTopMargin(0);
  pRight2->SetBottomMargin(0);
  TPad* pRight3 = new TPad("right_pad3","",0.500,0.501,0.975,0.739);
  pRight3->SetTopMargin(0);
  pRight3->SetBottomMargin(0);
  TPad* pRight4 = new TPad("right_pad4","",0.500,0.739,0.975,0.975);
  pRight4->SetBottomMargin(0);
  pLeft->Draw();
  pRight1->Draw();
  pRight2->Draw();
  pRight3->Draw();
  pRight4->Draw();
   
  // Left pad - distributions
  pLeft->cd();
  gPad->SetLogy();
  // Draw histograms
  h_t_Truth->Draw("hist");
  h_t_GoodCorr->Draw("same");
  h_tDIS_AllRec->Draw("same");
  h_t_B0Rec->Draw("same");
  h_t_RPRec->Draw("same");
  h_t_LCRec->Draw("same");
  // Add ePIC text label
  tePICLabel_tDIS->Draw("same");
  // Add legend
  TLegend* lt_Ratio = new TLegend(0.41,0.62,0.78,0.80);
  lt_Ratio->SetLineWidth(0);
  lt_Ratio->SetFillStyle(0);
  lt_Ratio->SetTextSize(0.032);
  lt_Ratio->AddEntry(h_t_Truth,"EpIC MC truth","l");
  lt_Ratio->AddEntry(h_t_GoodCorr,"Reco. DVCS (corrected)","lp");
  lt_Ratio->AddEntry(h_t_B0Rec,"Reco. DVCS (raw BABE - B0)","lp");
  lt_Ratio->AddEntry(h_t_RPRec,"Reco. DVCS (raw BABE - RP)","lp");
  lt_Ratio->AddEntry(h_t_LCRec,"Reco. DVCS (raw eXBE)","lp");
  lt_Ratio->AddEntry(h_tDIS_AllRec,"Reco. DIS (10x100)","l");
  lt_Ratio->Draw();
  // Right pad - ratios
  // Starting from top: DIS/DVCS corrected
  pRight4->cd();
  TH1D* hRatio_Corr = (TH1D*)h_tDIS_AllRec->Clone("hratio_corr");
  hRatio_Corr->Divide(h_t_GoodCorr);
  hRatio_Corr->SetLineColor(kBlack);
  hRatio_Corr->GetYaxis()->SetTitle("Ratio");
  hRatio_Corr->GetYaxis()->SetTitleSize(0.13);
  hRatio_Corr->GetYaxis()->SetTitleOffset(0.35);
  hRatio_Corr->GetYaxis()->SetLabelSize(0.13);
  hRatio_Corr->GetYaxis()->SetNdivisions(505);
  // Pad axes
  gPad->SetLogy();
  hRatio_Corr->GetYaxis()->SetRangeUser(2e-3,9);
  // Draw
  hRatio_Corr->Draw();
  // Add text
  TLatex* tePICLabel_RatCorr = new TLatex(0.12, 1.1, "All DIS/DVCS corrected");
  tePICLabel_RatCorr->SetTextSize(0.14);
  tePICLabel_RatCorr->Draw("same");
  // Next: DIS/DVCS (raw BABE B0)
  pRight3->cd();
  TH1D* hRatio_B0 = (TH1D*)h_tDIS_AllBABE->Clone("hratio_b0");
  hRatio_B0->Divide(h_t_B0Rec);
  hRatio_B0->SetLineColor(kBlack);
  hRatio_B0->GetYaxis()->SetTitle("Ratio");
  hRatio_B0->GetYaxis()->SetTitleSize(0.13);
  hRatio_B0->GetYaxis()->SetTitleOffset(0.35);
  hRatio_B0->GetYaxis()->SetLabelSize(0.13);
  // Pad axes
  gPad->SetLogy();
  hRatio_B0->GetYaxis()->SetRangeUser(2e-3,9);
  // Draw
  hRatio_B0->Draw();
  // Add text
  TLatex* tePICLabel_RatB0 = new TLatex(0.12, 1.1, "DIS - BABE/DVCS - BABE (B0)");
  tePICLabel_RatB0->SetTextSize(0.14);
  tePICLabel_RatB0->Draw("same");
  // Next: DIS/DVCS (raw BABE RP)
  pRight2->cd();
  TH1D* hRatio_RP = (TH1D*)h_tDIS_AllBABE->Clone("hratio_rp");
  hRatio_RP->Divide(h_t_RPRec);
  hRatio_RP->SetLineColor(kBlack);
  hRatio_RP->GetYaxis()->SetTitle("Ratio");
  hRatio_RP->GetYaxis()->SetTitleSize(0.13);
  hRatio_RP->GetYaxis()->SetTitleOffset(0.35);
  hRatio_RP->GetYaxis()->SetLabelSize(0.13);
  // Pad axes
  gPad->SetLogy();
  hRatio_RP->GetYaxis()->SetRangeUser(2e-3,9);
  // Draw
  hRatio_RP->Draw();
  // Add text
  TLatex* tePICLabel_RatRP = new TLatex(0.12, 1.1, "DIS - BABE/DVCS - BABE (RP)");
  tePICLabel_RatRP->SetTextSize(0.14);
  tePICLabel_RatRP->Draw("same");
  // Last: DIS/DVCS (raw eXBE)
  pRight1->cd();
  TH1D* hRatio_LC = (TH1D*)h_tDIS_AllRec->Clone("hratio_lc");
  hRatio_LC->Divide(h_t_LCRec);
  hRatio_LC->SetLineColor(kBlack);
  hRatio_LC->GetYaxis()->SetTitle("Ratio");
  hRatio_LC->GetYaxis()->SetTitleSize(0.13);
  hRatio_LC->GetYaxis()->SetTitleOffset(0.35);
  hRatio_LC->GetYaxis()->SetLabelSize(0.13);
  hRatio_LC->GetYaxis()->SetNdivisions(505);
  hRatio_LC->GetXaxis()->SetTitle("|t| [GeV^{2}]");
  hRatio_LC->GetXaxis()->SetTitleSize(0.13);
  hRatio_LC->GetXaxis()->SetTitleOffset(0.95);
  hRatio_LC->GetXaxis()->SetLabelSize(0.14);
  // Pad axes
  gPad->SetLogy();
  hRatio_LC->GetYaxis()->SetRangeUser(2e-3,9);
  // Draw
  hRatio_LC->Draw();
  // Add text
  TLatex* tePICLabel_RatLC = new TLatex(0.12, 1.1, "DIS - all/DVCS - eXBE");
  tePICLabel_RatLC->SetTextSize(0.14);
  tePICLabel_RatLC->Draw("same");
  if(kSAVE) ct_Ratio->SaveAs("../figs/DISComp_" + energy +"_Ratio.png");

  
  // Print ratios
  cout<<"|t| [GeV^2]\tDIS(BABE)/DVCS(B0)\tDIS(BABE)/DVCS(RP)\tDIS(All)/DVCS(Corr.)"<<endl;
  for(int bin{1}; bin<= h_t_Truth->GetNbinsX(); bin++){
    cout<<"["<<h_t_Truth->GetBinLowEdge(bin)<<" - "<<h_t_Truth->GetBinLowEdge(bin+1)<<"]\t"<<hRatio_B0->GetBinContent(bin)<<"\t"<<hRatio_RP->GetBinContent(bin)<<"\t"<<hRatio_Corr->GetBinContent(bin)<<endl;
  }

  // Close canvases
  /*ct_DISlo->Close();
  ct_DIShi->Close();
  ct_DISall->Close();
  ct_DISeg->Close();
  ct_Ratio->Close();*/

  return;
}
