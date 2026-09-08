using namespace std;

#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

// Ragged Q2/xB/|t| binning reader - the SAME one used by the analysis.
#include "DVCSBinning.hh"
// ePIC plotting style
#include "../DVCS_ep/ePIC_style.C"

// DECIDE ON OUTPUT BEHAVIOUR
Bool_t kPRINT{kFALSE};
Bool_t kFIGS{kFALSE};
Bool_t kDEBUG{kFALSE};

// Calculate histogram scaling factors to get to desired luminosity
// Give luminosity in full scientific (expect 1e15 for 1fb-1)
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
// Trento Phi and Asymmetries - 26.07.1 campaign and onwards
//---------------------------------------------------------------------

void DVCSAsymPlots(TString campaign = "26.07.1", TString energy = "9x130", TString haddir = "T"){
  if(haddir != "T" && haddir != "L"){
    cout<<"Invalid hadron polarisation direction \""<<haddir<<"\" (valid - L/T)"<<endl;
    return;
  }

  // Print plot settings
  cout<<"---------------------------------------------"<<endl;
  cout<<"Processing TDR plots - Asymmetries (";
  if(haddir == "T") cout<<"transverse hadron)"<<endl;
  if(haddir == "L") cout<<"longitudinal hadron)"<<endl;
  cout<<"\tCampaign: "<<campaign<<endl;
  cout<<"\tEnergy: "<<energy<<endl;
  if(energy == "9x130") cout<<"\tAsymmetries for 1 fb-1 total ep"<<endl;
  else if(energy == "9x275") cout<<"\tAsymmetries for 2.5 fb-1 total ep"<<endl;
  cout<<"---------------------------------------------"<<endl;
  
  // Using ePIC global style
  //gROOT->ProcessLine("set_ePIC_style()");
  //gStyle->SetCanvasPreferGL(kTRUE);

  // Set beam energies
  Float_t fEBeam{0}, fPBeam{0};
  fEBeam = 9.;
  
  if(energy == "9x130") fPBeam = 130.;
  else if(energy == "9x275") fPBeam = 275;
  else{
    cout<<"Invalid beam energy."<<endl;
    return;
  }

  TString path = "$EIC_WORK_DIR/DVCS_Analysis/RootFiles/";
  //path = "/scratch/oliver/";
  TString sInMM = path + "ePIC_DVCS_"+campaign+"_"+energy+"_emh"+haddir+"m.root";
  gSystem->ExpandPathName(sInMM);
  TString sInMP = path + "ePIC_DVCS_"+campaign+"_"+energy+"_emh"+haddir+"p.root";
  gSystem->ExpandPathName(sInMP);
  TString sInPM = path + "ePIC_DVCS_"+campaign+"_"+energy+"_eph"+haddir+"m.root";
  gSystem->ExpandPathName(sInPM);
  TString sInPP = path + "ePIC_DVCS_"+campaign+"_"+energy+"_eph"+haddir+"p.root";
  gSystem->ExpandPathName(sInPP);

  // Check files exist  
  TFile* fInMM;
  TFile* fInMP;
  TFile* fInPM;
  TFile* fInPP;
  bool kFileCheck{true};
  if(gSystem->AccessPathName(sInMM)){
    cout<<"ERROR: File "+sInMM+" does not exist"<<endl;
    kFileCheck = false;
  }
  if(gSystem->AccessPathName(sInMP)){
    cout<<"File "+sInMP+" does not exist"<<endl;
    kFileCheck = false;
  }
  if(gSystem->AccessPathName(sInPM)){
    cout<<"File "+sInPM+" does not exist"<<endl;
    kFileCheck = false;
  }
  if(gSystem->AccessPathName(sInPP)){
    cout<<"File "+sInPP+" does not exist"<<endl;
    kFileCheck = false;
  }
  
  if(!kFileCheck) return;
  else{
    fInMM = TFile::Open(sInMM);
    fInMP = TFile::Open(sInMP);
    fInPM = TFile::Open(sInPM);
    fInPP = TFile::Open(sInPP);
  }
  
  //---------------------------------------------------------------------
  // Ragged binning (same file the analysis used to book the histograms)
  //---------------------------------------------------------------------
  // NOTE: keep bins.txt reachable from the run directory, exactly like the
  //       analysis. Change the name here if it is configuration-specific.
  DVCSBinning binning;
  TString sBinFile = "bins_"+energy+".txt";
  if(!binning.load(sBinFile.Data())){
    cout<<"[DVCSAsymPlots] FATAL: could not load bins.txt"<<endl;
    return;
  }
  //binning.print();

  gStyle->SetOptStat(00000000);

  int page = 0; // global page counter -> zero-padded temp PDFs stay ordered

  // Small helper: colour a histogram and return it (nullptr-safe).
  auto style = [](TH1D* h, int col)->TH1D*{
    if(!h) return nullptr;
    h->SetLineColor(col);
    h->SetLineWidth(2);
    return h;
  };

  // Calculate scaling factors for EIC luminosities
  double EIClumi{1.};
  if(energy == "9x130") EIClumi = 1e15;
  else if(energy == "9x275") EIClumi = 2.5e15;

  // Scaling/division based on expected polarisation of beams
  EIClumi /= 2;

  double scaleMM = calcScaling(energy, "emhTm", EIClumi/8);
  double scaleMP = calcScaling(energy, "emhTp", EIClumi/8);
  double scalePM = calcScaling(energy, "ephTm", EIClumi/8);
  double scalePP = calcScaling(energy, "ephTp", EIClumi/8);

  // Print phi_h/Q2/x/t binning
  // Load single histo for printing
  TH1D* hPhiBins = (TH1D*)fInPP->Get("tphi_b0mc[0][0][0]")->Clone("tphibins");
  // Check if phi bins in degrees or rad (want to print in rad)
  bool kPhiRad{true};
  if(hPhiBins->GetBinLowEdge(hPhiBins->GetNbinsX()+1) == 360) kPhiRad = false;
  
  cout<<"\n---------------------------------------------------------------------------------------------------"<<endl;
  cout<<"Printing binning:"<<endl;
  
  cout<<"Bin centres (phi_h): [";
  for(int phibin{0}; phibin < hPhiBins->GetNbinsX(); phibin++){
    if(!kPhiRad) cout<<hPhiBins->GetBinCenter(phibin+1)*TMath::DegToRad();
    else cout<<hPhiBins->GetBinCenter(phibin+1);
    if(phibin < hPhiBins->GetNbinsX()-1) cout<<",";
  }
  cout<<"]"<<endl;
  cout<<"---------------------------------------------------------------------------------------------------"<<endl;


  for(int q=0; q<binning.nQ2(); q++){
    for(int region=DVCSBinning::kB0; region<=DVCSBinning::kRP; region++){

      const char* preg  = (region==DVCSBinning::kB0) ? "b0" : "rp";
      const char* rname = (region==DVCSBinning::kB0) ? "B0" : "RP";

      for(int x=0; x<binning.nXB(region,q); x++){
        const int nt = binning.nT(region,q,x);
	for(int t=0; t<nt; t++){
	  // Set up canvas
	  TCanvas* c = new TCanvas(Form("cDiff_q%d_%s_x%d_t%d",q,preg,x,t),"",1200,800);
	  c->Divide(2,2);
	  
	  // Extract and scale histograms
	  // Gen.
	  TH1D* hmc_mm = style((TH1D*)fInMM->Get(Form("tphi_%smc[%d][%d][%d]",preg,q,x,t))->Clone("tphi_recmm"), kBlue);
	  TH1D* hmc_mp = style((TH1D*)fInMP->Get(Form("tphi_%smc[%d][%d][%d]",preg,q,x,t))->Clone("tphi_recmp"), kRed);
	  TH1D* hmc_pm = style((TH1D*)fInPM->Get(Form("tphi_%smc[%d][%d][%d]",preg,q,x,t))->Clone("tphi_recpm"), kBlack);
	  TH1D* hmc_pp = style((TH1D*)fInPP->Get(Form("tphi_%smc[%d][%d][%d]",preg,q,x,t))->Clone("tphi_recpp"), kGreen+2);
	  hmc_mm->Scale(scaleMM);
	  hmc_mp->Scale(scaleMP);
	  hmc_pm->Scale(scalePM);
	  hmc_pp->Scale(scalePP);
	  // Reco.
	  TH1D* hrec_mm = style((TH1D*)fInMM->Get(Form("tphi_%sreco[%d][%d][%d]",preg,q,x,t))->Clone("tphi_recmm"), kBlue);
	  TH1D* hrec_mp = style((TH1D*)fInMP->Get(Form("tphi_%sreco[%d][%d][%d]",preg,q,x,t))->Clone("tphi_recmp"), kRed);
	  TH1D* hrec_pm = style((TH1D*)fInPM->Get(Form("tphi_%sreco[%d][%d][%d]",preg,q,x,t))->Clone("tphi_recpm"), kBlack);
	  TH1D* hrec_pp = style((TH1D*)fInPP->Get(Form("tphi_%sreco[%d][%d][%d]",preg,q,x,t))->Clone("tphi_recpp"), kGreen+2);
	  hrec_mm->Scale(scaleMM);
	  hrec_mp->Scale(scaleMP);
	  hrec_pm->Scale(scalePM);
	  hrec_pp->Scale(scalePP);

	  float Pe = 0.7; // electron beam polarisation - 70%
	  float Ph = 0.7; // hadron beam polarisation - 70%

	  // Create subhistograms for asymmetries
	  // 0) Full data sum - needed for all asymmatries
	  TH1D* hsum = (TH1D*)hrec_pp->Clone("phi_sum");
	  hsum->Add(hrec_pm);
	  hsum->Add(hrec_mp);
	  hsum->Add(hrec_mm);

	  // 1) Beam Spin Asymmetry - sum over hadron helicity states
	  TH1D* hbsa_plus = (TH1D*)hrec_pp->Clone("bsa_plus");
	  hbsa_plus->Add(hrec_pm);
	  TH1D* hbsa_minus = (TH1D*)hrec_mm->Clone("bsa_minus");
	  hbsa_minus->Add(hrec_mp);	  
	  // Scale to beam polarisations; NEED TO OVERWRITE ERRORS
	  hbsa_plus->Scale(Pe);
	  hbsa_minus->Scale(Pe);
	  for(int bin{1}; bin<hbsa_plus->GetNbinsX()+1;bin++){
	    hbsa_plus->SetBinError(bin,TMath::Sqrt(hbsa_plus->GetBinContent(bin)));
	    hbsa_minus->SetBinError(bin,TMath::Sqrt(hbsa_minus->GetBinContent(bin)));
	  }
	  
	  TH1D* hbsa_diff = (TH1D*)hbsa_plus->Clone("bsa_diff");
	  hbsa_diff->Add(hbsa_minus,-1);
	  TH1D* hbsa_sum = (TH1D*)hbsa_plus->Clone("bsa_sum");
	  hbsa_sum->Add(hbsa_minus);
	  TH1D* hbsa = (TH1D*)hbsa_diff->Clone("hbsa");
	  hbsa->Divide(hbsa_sum);

	  TH1D* hbsamc_plus = (TH1D*)hmc_pp->Clone("bsamc_plus");
	  hbsamc_plus->Add(hmc_pm);
	  TH1D* hbsamc_minus = (TH1D*)hmc_mm->Clone("bsamc_minus");
	  hbsamc_minus->Add(hmc_mp);
	  TH1D* hbsamc = (TH1D*)hbsamc_plus->GetAsymmetry(hbsamc_minus);

	  // 2) Target Spin Asymmetry - sum over electron helicity states
	  TH1D* htsa_plus = (TH1D*)hrec_pp->Clone("tsa_plus");
	  htsa_plus->Add(hrec_mp);
	  TH1D* htsa_minus = (TH1D*)hrec_mm->Clone("tsa_minus");
	  htsa_minus->Add(hrec_pm);
	  // Scale to beam polarisations; NEED TO OVERWRITE ERRORS
	  htsa_plus->Scale(Ph);
	  htsa_minus->Scale(Ph);
	  for(int bin{1}; bin<htsa_plus->GetNbinsX()+1;bin++){
	    htsa_plus->SetBinError(bin,TMath::Sqrt(htsa_plus->GetBinContent(bin)));
	    htsa_minus->SetBinError(bin,TMath::Sqrt(htsa_minus->GetBinContent(bin)));
	  }
	  
	  TH1D* htsa_diff = (TH1D*)htsa_plus->Clone("tsa_diff");
	  htsa_diff->Add(htsa_minus,-1);
	  TH1D* htsa_sum = (TH1D*)htsa_plus->Clone("tsa_sum");
	  htsa_sum->Add(htsa_minus);
	  TH1D* htsa = (TH1D*)htsa_diff->Clone("htsa");
	  htsa->Divide(htsa_sum);
	  
	  TH1D* htsamc_plus = (TH1D*)hmc_pp->Clone("tsamc_plus");
	  htsamc_plus->Add(hmc_mp);
	  TH1D* htsamc_minus = (TH1D*)hmc_mm->Clone("tsamc_minus");
	  htsamc_minus->Add(hmc_pm);
	  TH1D* htsamc = (TH1D*)htsamc_plus->GetAsymmetry(htsamc_minus);
	  
	  // 3) Double Spin Asymmetry - same spin minus different spin
	  TH1D* hdsa_plus = (TH1D*)hrec_pp->Clone("dsa_plus");
	  hdsa_plus->Add(hrec_mm);
	  TH1D* hdsa_minus = (TH1D*)hrec_mp->Clone("dsa_minus");
	  hdsa_minus->Add(hrec_pm);
	  // Scale to beam polarisations; NEED TO OVERWRITE ERRORS
	  hdsa_plus->Scale(Pe);
	  hdsa_plus->Scale(Ph);
	  hdsa_minus->Scale(Pe);
	  hdsa_minus->Scale(Ph);
	  for(int bin{1}; bin<hdsa_plus->GetNbinsX()+1;bin++){
	    hdsa_plus->SetBinError(bin,TMath::Sqrt(hdsa_plus->GetBinContent(bin)));
	    hdsa_minus->SetBinError(bin,TMath::Sqrt(hdsa_minus->GetBinContent(bin)));
	  }
	  
	  TH1D* hdsa_diff = (TH1D*)hdsa_plus->Clone("dsa_diff");
	  hdsa_diff->Add(hdsa_minus,-1);
	  hdsa_diff->Scale(Pe);
	  hdsa_diff->Scale(Ph);
	  TH1D* hdsa_sum = (TH1D*)hdsa_plus->Clone("dsa_sum");
	  hdsa_sum->Add(hdsa_minus);
	  TH1D* hdsa = (TH1D*)hdsa_diff->Clone("hdsa");
	  hdsa->Divide(hdsa_sum);

	  TH1D* hdsamc_plus = (TH1D*)hmc_pp->Clone("dsamc_plus");
	  hdsamc_plus->Add(hmc_mm);
	  TH1D* hdsamc_minus = (TH1D*)hmc_mp->Clone("dsamc_minus");
	  hdsamc_minus->Add(hmc_pm);
	  TH1D* hdsamc = (TH1D*)hdsamc_plus->GetAsymmetry(hdsamc_minus);

	  // Scale asymmetries for polarisations
	  hbsa->Scale(1./Pe);
	  htsa->Scale(1./Ph);
	  hdsa->Scale(1./(Pe*Ph));
	  /*hbsamc->Scale(1./Pe);
	  htsamc->Scale(1./Ph);
	  hdsamc->Scale(1./(Pe*Ph));*/

	  // Pad 1 - reconstructed phi
	  c->cd(1);
	  // Histogram with maximum bin defines the frame; scale y to fit all.
          double ymax=0.;
	  for(TH1D* h : {hrec_mm,hrec_mp,hrec_pm,hrec_pp}) if(h){ if(h->GetMaximum()>ymax) ymax=std::max(ymax,h->GetMaximum()); }
	  
	  // Start from first polarisation state
	  hrec_mm->SetTitle(Form("%s: Q^{2}[%.2f,%.2f]  x_{B}[%.4g,%.4g]  |t|[%.3g,%.3g];#phi_{h} [deg];Counts/%.3g deg",
				 rname,
				 binning.q2Low(q), binning.q2High(q),
				 binning.xBLow(region,q,x), binning.xBHigh(region,q,x),
				 binning.tLow(region,q,x,t), binning.tHigh(region,q,x,t),
				 hrec_mm->GetBinLowEdge(hrec_mm->GetNbinsX()+1)/hrec_mm->GetNbinsX()));
          hrec_mm->SetMinimum(0.);
          hrec_mm->SetMaximum(1.3*ymax + 1.);
	  hrec_mm->GetYaxis()->SetTitleSize(0.05);
	  hrec_mm->GetYaxis()->SetTitleOffset(1.);
          hrec_mm->GetYaxis()->SetLabelSize(0.05);
	  hrec_mm->GetXaxis()->SetTitleSize(0.05);
	  hrec_mm->GetXaxis()->SetTitleOffset(0.90);
          hrec_mm->GetXaxis()->SetLabelSize(0.05);
	  hrec_mm->Draw("hist");
          for(TH1D* h : {hrec_mp,hrec_pm,hrec_pp}) if(h) h->Draw("hist same");
	  TLegend* leg = new TLegend(0.30,0.65,0.75,0.89);
	  leg->SetLineWidth(0);
	  leg->SetFillStyle(1);	 
	  leg->SetNColumns(2);
	  leg->AddEntry(hrec_mm,"e-, p-","l");
	  leg->AddEntry(hrec_mp,"e-, p+","l");
	  leg->AddEntry(hrec_pm,"e+, p-","l");
	  leg->AddEntry(hrec_pp,"e+, p+","l");
	  leg->Draw();
	  
	  // Pad 2 - BSA
	  c->cd(2);
	  ymax = 0;
	  for(TH1D* h : {hbsa,hbsamc}) if(h){ if(h->GetMaximum()>ymax) ymax=std::max(ymax,h->GetMaximum()); }
	  style(hbsa, kBlue);
	  hbsa->SetMarkerStyle(kFullDiamond);
	  hbsa->SetMarkerColor(kBlue);
	  hbsa->SetMarkerSize(2);
	  style(hbsamc, kBlack);
	  hbsamc->SetMarkerStyle(kOpenCircle);
	  hbsamc->SetMarkerColor(kBlack);
	  hbsamc->SetMarkerSize(2);
	  hbsa->SetTitle(";#phi_{h} [deg];A_{LU}");
	  hbsa->GetYaxis()->SetTitleSize(0.05);
	  hbsa->GetYaxis()->SetTitleOffset(1.);
          hbsa->GetYaxis()->SetLabelSize(0.05);
	  hbsa->GetXaxis()->SetTitleSize(0.05);
	  hbsa->GetXaxis()->SetTitleOffset(0.90);
          hbsa->GetXaxis()->SetLabelSize(0.05);
	  hbsa->SetMinimum(-2*ymax);
          hbsa->SetMaximum(2*ymax);
	  hbsa->Draw();
	  hbsamc->Draw("same");
	  TLegend* leg2 = new TLegend(0.57,0.68,0.96,0.87);
	  leg2->SetLineWidth(0);
	  leg2->SetFillStyle(1);
	  leg2->AddEntry(hbsamc,"MC gen.","lp");
	  leg2->AddEntry(hbsa,"Raw reco.","lp");
	  leg2->Draw();

	  // Pad 3 - TSA
	  c->cd(3);
	  ymax = 0;
	  for(TH1D* h : {htsa,htsamc}) if(h){ if(h->GetMaximum()>ymax) ymax=std::max(ymax,h->GetMaximum()); }
	  style(htsa, kBlue);
	  htsa->SetMarkerStyle(kFullDiamond);
	  htsa->SetMarkerColor(kBlue);
	  htsa->SetMarkerSize(2);
	  style(htsamc, kBlack);
	  htsamc->SetMarkerStyle(kOpenCircle);
	  htsamc->SetMarkerColor(kBlack);
	  htsamc->SetMarkerSize(1.5);
	  if(haddir == "T") htsa->SetTitle(";#phi_{h} [deg];A_{UT}");
	  else htsa->SetTitle(";#phi_{h} [deg];A_{UL}");
	  htsa->GetYaxis()->SetTitleSize(0.05);
	  htsa->GetYaxis()->SetTitleOffset(1.);
          htsa->GetYaxis()->SetLabelSize(0.05);
	  htsa->GetXaxis()->SetTitleSize(0.05);
	  htsa->GetXaxis()->SetTitleOffset(0.90);
          htsa->GetXaxis()->SetLabelSize(0.05);
	  htsa->SetMinimum(-2*ymax);
          htsa->SetMaximum(2*ymax);
	  htsa->Draw();
	  htsamc->Draw("same");
	  
	  // Pad 4 - DSA
	  c->cd(4);
	  ymax = 0;
	  for(TH1D* h : {hdsa,hdsamc}) if(h){ if(h->GetMaximum()>ymax) ymax=std::max(ymax,h->GetMaximum()); }
	  style(hdsa, kBlue);
	  hdsa->SetMarkerStyle(kFullDiamond);
	  hdsa->SetMarkerColor(kBlue);
	  hdsa->SetMarkerSize(2);
	  style(hdsamc, kBlack);
	  hdsamc->SetMarkerStyle(kOpenCircle);
	  hdsamc->SetMarkerColor(kBlack);
	  hdsamc->SetMarkerSize(2);
	  if(haddir == "T") hdsa->SetTitle(";#phi_{h} [deg];A_{LT}");
	  else hdsa->SetTitle(";#phi_{h} [deg];A_{LL}");
	  hdsa->GetYaxis()->SetTitleSize(0.05);
	  hdsa->GetYaxis()->SetTitleOffset(1.);
          hdsa->GetYaxis()->SetLabelSize(0.05);
	  hdsa->GetXaxis()->SetTitleSize(0.05);
	  hdsa->GetXaxis()->SetTitleOffset(0.90);
          hdsa->GetXaxis()->SetLabelSize(0.05);
	  hdsa->SetMinimum(-2*ymax);
          hdsa->SetMaximum(2*ymax);
	  hdsa->Draw();
	  hdsamc->Draw("same");
	  
	  c->Print(Form("DVCSasym_temp%03d.pdf", ++page));
	  c->Close();

	  // Print Asymmetry values
	  float counts =  hrec_mm->Integral()+hrec_mp->Integral()+hrec_pm->Integral()+hrec_pp->Integral();
	  TH1D* h_Q2Diff_Tot = (TH1D*)fInMM->Get(Form("q2diff_%s[%d][%d][%d]",preg,q,x,t))->Clone("q2tot");
	  h_Q2Diff_Tot->Add((TH1D*)fInMP->Get(Form("q2diff_%s[%d][%d][%d]",preg,q,x,t)));
	  h_Q2Diff_Tot->Add((TH1D*)fInPM->Get(Form("q2diff_%s[%d][%d][%d]",preg,q,x,t)));
	  h_Q2Diff_Tot->Add((TH1D*)fInPP->Get(Form("q2diff_%s[%d][%d][%d]",preg,q,x,t)));
	  float q2mean = h_Q2Diff_Tot->GetMean();
	  TH1D* h_XBDiff_Tot = (TH1D*)fInMM->Get(Form("xbdiff_%s[%d][%d][%d]",preg,q,x,t))->Clone("xbtot");
	  h_XBDiff_Tot->Add((TH1D*)fInMP->Get(Form("xbdiff_%s[%d][%d][%d]",preg,q,x,t)));
	  h_XBDiff_Tot->Add((TH1D*)fInPM->Get(Form("xbdiff_%s[%d][%d][%d]",preg,q,x,t)));
	  h_XBDiff_Tot->Add((TH1D*)fInPP->Get(Form("xbdiff_%s[%d][%d][%d]",preg,q,x,t)));
	  float xbmean = h_XBDiff_Tot->GetMean();
	  TH1D* h_TDiff_Tot = (TH1D*)fInMM->Get(Form("tdiff_%s[%d][%d][%d]",preg,q,x,t))->Clone("ttot");
	  h_TDiff_Tot->Add((TH1D*)fInMP->Get(Form("tdiff_%s[%d][%d][%d]",preg,q,x,t)));
	  h_TDiff_Tot->Add((TH1D*)fInPM->Get(Form("tdiff_%s[%d][%d][%d]",preg,q,x,t)));
	  h_TDiff_Tot->Add((TH1D*)fInPP->Get(Form("tdiff_%s[%d][%d][%d]",preg,q,x,t)));
	  float tmean = h_TDiff_Tot->GetMean();
	  
	  if(region == DVCSBinning::kB0) cout<<"\nHigh-t: ["<<q<<"]["<<x<<"]["<<t<<"] = ["<<q2mean<<"]["<<xbmean<<"]["<<tmean<<"]"<<endl;
	  else cout<<"\nLow-t: ["<<q<<"]["<<x<<"]["<<t<<"] = ["<<q2mean<<"]["<<xbmean<<"]["<<tmean<<"]"<<endl;
	  cout<<"  - Total raw events (scaled to "<<EIClumi/1e15<<" fb-1) = "<<counts<<endl;
	  
	  cout<<"  - A_LU (raw rec.): [";
	  for(int phi{0}; phi < hPhiBins->GetNbinsX(); phi++){
	    cout<<hbsa->GetBinContent(phi+1);
	    if(phi < hPhiBins->GetNbinsX()-1) cout<<",";
	  }
	  cout<<"]"<<endl;
	  
	  cout<<"  - dA_LU (raw rec.): [";
	  for(int phi{0}; phi < hPhiBins->GetNbinsX(); phi++){
	    cout<<hbsa->GetBinError(phi+1);
	    if(phi < hPhiBins->GetNbinsX()-1) cout<<",";
	  }
	  cout<<"]"<<endl;

	  if(haddir == "T") cout<<"  - A_UT (raw rec.): [";
	  else if(haddir == "L") cout<<"  - A_UL (raw rec.): [";
	  for(int phi{0}; phi < hPhiBins->GetNbinsX(); phi++){
	    cout<<htsa->GetBinContent(phi+1);
	    if(phi < hPhiBins->GetNbinsX()-1) cout<<",";
	  }
	  cout<<"]"<<endl;
	  
	  if(haddir == "T") cout<<"  - dA_UT (raw rec.): [";
	  else if(haddir == "L") cout<<"  - dA_UL (raw rec.): [";
	  for(int phi{0}; phi < hPhiBins->GetNbinsX(); phi++){
	    cout<<htsa->GetBinError(phi+1);
	    if(phi < hPhiBins->GetNbinsX()-1) cout<<",";
	  }
	  cout<<"]"<<endl;

	  if(haddir == "T") cout<<"  - A_LT (raw rec.): [";
	  else if(haddir == "L") cout<<"  - A_LL (raw rec.): [";
	  for(int phi{0}; phi < hPhiBins->GetNbinsX(); phi++){
	    cout<<hdsa->GetBinContent(phi+1);
	    if(phi < hPhiBins->GetNbinsX()-1) cout<<",";
	  }
	  cout<<"]"<<endl;
	  
	  if(haddir == "T") cout<<"  - dA_LT (raw rec.): [";
	  else if(haddir == "L") cout<<"  - dA_LL (raw rec.): [";
	  for(int phi{0}; phi < hPhiBins->GetNbinsX(); phi++){
	    cout<<hdsa->GetBinError(phi+1);
	    if(phi < hPhiBins->GetNbinsX()-1) cout<<",";
	  }
	  cout<<"]"<<endl;

	  cout<<"  - Sigma(N) (raw scaled): [";
	  for(int phi{0}; phi < hPhiBins->GetNbinsX(); phi++){
	    float Sigma = hrec_mm->GetBinContent(phi+1) + hrec_mp->GetBinContent(phi+1) + hrec_pm->GetBinContent(phi+1) + hrec_pp->GetBinContent(phi+1);
	    cout<<Sigma;
	    if(phi < hPhiBins->GetNbinsX()-1) cout<<",";
	  }
	  cout<<"]"<<endl;

	} // rof (t bins)
      } // rof (xB bins)
    } // rof (detector regions)
  } // rof (Q2 bins)

  //std::cout<<"...Cleaning up files..."<<std::endl;
  TString filePlots = "$EIC_WORK_DIR/DVCS_Analysis/Plots/DVCSAsymPlots_" + campaign + "_" + energy + "_" + haddir + ".pdf";
  //std::cout<<"Moving plots to "<<filePlots<<std::endl;
  gSystem->Exec("pdfunite DVCSasym_temp*.pdf " + filePlots);
  gSystem->Exec("rm DVCSasym_temp*.pdf");

  return;
}
