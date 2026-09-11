using namespace std;

#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include "./ePIC_style.C"

bool kSAVE = false;

//---------------------------------------------------------------------
// MAIN
// 
// DIS physics background plots to DVCS analysis
//---------------------------------------------------------------------
void Plots_DISBkg(TString campaign = "26.02.0", TString energy = "10x130", TString hel = "minus"){
  // Print plot settings
  cout<<"\n--------------------------------------"<<endl;
  cout<<"Processing DIS background plots"<<endl;
  cout<<"\tCampaign: "<<campaign<<endl;
  cout<<"\tDVCS sample energy: "<<energy<<endl;
  cout<<"\te- DVCS sample beam helicity: "<<hel<<endl;
  cout<<"--------------------------------------\n"<<endl;

  // Set beam energies
  Float_t fEBeam{0}, fPBeam{0};
  if(energy == "5x41"){
    fEBeam = 5.;
    fPBeam = 41.;
  }
  else if(energy == "10x100" || energy == "10x130" || energy == "10x250"){
    fEBeam = 10.;
    if(energy == "10x100") fPBeam = 100.;
    if(energy == "10x130") fPBeam = 130.;
    if(energy == "10x250") fPBeam = 250.;
  }
  else if(energy == "18x275"){
    fEBeam = 18.;
    fPBeam = 275;
  }
  else{
    cout<<"Invalid beam energy."<<endl;
    return;
  }

  //--------------------------------------------------------------------
  // Load histograms from file - DVCS baseline
  //--------------------------------------------------------------------
  TString sDVCS = "$EIC_WORK_DIR/DVCS_Analysis/RootFiles/ePIC_DVCS_"+campaign+"_"+energy+"_"+hel+".root";
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
  TFile* fDISLQ2 = TFile::Open("$EIC_WORK_DIR/DVCS_Analysis/RootFiles/ePIC_DIS_26.02.0_10x100_minQ2=1_MoreVetoes.root");
  TH1D* h_tDISLQ2_B0Rec = (TH1D*)fDISLQ2->Get("t_b0reco");
  TH1D* h_tDISLQ2_RPRec = (TH1D*)fDISLQ2->Get("t_rpreco");
  TH1D* h_tDISLQ2_LCRec = (TH1D*)fDISLQ2->Get("t_lcreco");
  TFile* fDISHQ2 = TFile::Open("$EIC_WORK_DIR/DVCS_Analysis/RootFiles/ePIC_DIS_26.02.0_10x100_minQ2=10_MoreVetoes.root");
  TH1D* h_tDISHQ2_B0Rec = (TH1D*)fDISHQ2->Get("t_b0reco");
  TH1D* h_tDISHQ2_RPRec = (TH1D*)fDISHQ2->Get("t_rpreco");
  TH1D* h_tDISHQ2_LCRec = (TH1D*)fDISHQ2->Get("t_lcreco");
 
  //---------------------------------------------------------------------
  // Calculations on histograms
  // Scaling to 5fb-1
  //---------------------------------------------------------------------
  // Calculations - SCALING FACTOR TO 5FB-1
  Double_t fXSint{0}, NEv{0};
  if(energy == "5x41"){
    NEv = 1e6;
    if(hel == "minus") fXSint = 5.02849207976343e-9;
    if(hel == "plus") fXSint = 5.0297340302691e-9;
  }
  else if(energy == "10x100"){
    NEv = 1.3e6;
    if(hel == "minus") fXSint = 6.49238407587587e-9;
    if(hel == "plus") fXSint = 6.49388763015155e-9;
  }
  else if(energy == "10x130"){
    NEv = 1.3e6;
    if(hel == "minus") fXSint = 6.73508639641424e-9;
    if(hel == "plus") fXSint = 6.74022528459919e-9;
  }
  else if(energy == "10x250"){
    NEv = 1e6;
    if(hel == "minus") fXSint = 7.34712332835456e-9;
    if(hel == "plus") fXSint = 7.34179803172751e-9;
  }
  else if(energy == "18x275"){
    NEv = 1e6;
    if(hel == "minus") fXSint = 7.96572820985694e-9;
    if(hel == "plus") fXSint = 7.95421883076817e-9;
  }
  Double_t lumi = NEv/fXSint;
  Double_t scaleTo5 = 5e15/lumi;

  // Need 5fb-1 scaling factor for DIS data too
  Double_t fXS_DISlo{5.55970643e-7};
  Double_t fXS_DIShi{3.99634216e-8};
  Double_t lumi_DISlo = 1e6/fXS_DISlo;
  Double_t lumi_DIShi = 1e6/fXS_DIShi;
  Double_t scaleTo5_DISlo = 5e15/lumi_DISlo;
  Double_t scaleTo5_DIShi = 5e15/lumi_DIShi;

  cout<<"DVCS data represents "<<lumi/1e15<<" fb-1\n\tDIS low Q2 = "<<lumi_DISlo/1e15<<" fb-1\n\tDIS high Q2 = "<<lumi_DIShi/1e15<<" fb-1\n"<<endl;

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
  h_t_Truth->Scale(scaleTo5);
  h_t_GoodCorr->Scale(scaleTo5);
  h_t_GoodReco->Scale(scaleTo5);
  h_tDISLQ2_B0Rec->Scale(scaleTo5_DISlo);
  h_tDISLQ2_RPRec->Scale(scaleTo5_DISlo);
  h_tDISHQ2_B0Rec->Scale(scaleTo5_DIShi);
  h_tDISHQ2_RPRec->Scale(scaleTo5_DIShi);
  
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
  if(kSAVE) ct_DISlo->SaveAs("figs/DISComp_" + energy +"_t_Q2lo.png");
  
  
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
  if(kSAVE) ct_DIShi->SaveAs("figs/DISComp_" + energy +"_t_Q2hi.png");
  
  
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
  if(kSAVE) ct_DISall->SaveAs("figs/DISComp_" + energy +"_t_Q2hi.png");
  

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
  if(kSAVE) ct_DISeg->SaveAs("figs/DISComp_" + energy +"_t_Q2hi.png");


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
  if(kSAVE) ct_Ratio->SaveAs("figs/DISComp_" + energy +"_Ratio.png");

  
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
