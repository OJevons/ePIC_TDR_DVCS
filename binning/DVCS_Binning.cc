using namespace std;

#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include "./ePIC_ExcKinUtils.h"

// Calculate histogram scaling factors to get to desired luminosity
// Give luminosity in full scientific (expect 1e15 for 1fb-1)
double calcScaling(TString energy, TString hel, float lumi){
  // Holding variables - No. of events generated, integrated cross-section
  // These vary by beam settings
  Double_t fXSint{0}, NEv{0};
  if(energy == "9x130"){
    NEv = 1e6;
    if(hel == "emhTm") fXSint = 6.91502943051105e-9;
    else if(hel == "emhTp") fXSint = 7.00501745732437e-9;
    else if(hel == "ephTm") fXSint = 6.92070657974331e-9;
    else if(hel == "ephTp") fXSint = 7.00302210771398e-9;
  }
  else if(energy == "9x275"){
    NEv = 2.5e6;
    if(hel == "emhTm") fXSint = 7.65782690151261e-9;
    else if(hel == "emhTp") fXSint = 7.73125664480267e-9;
    else if(hel == "ephTm") fXSint = 7.65495841373398e-9;
    else if(hel == "ephTp") fXSint = 7.73434664504977e-9;
  }

  Double_t genlumi = NEv/fXSint;
  Double_t scale = lumi/genlumi;

  return scale;
}

bool checkPol(TString ehel, TString pdir, TString phel){
  // Check polarization flags match naming convention
  bool goodPol{true};
  if(ehel != "m" && ehel != "p"){
    cout<<"Electron helicity state "<<ehel<<" invalid. Use \"p/m\"."<<endl;
    goodPol = false;
  }
  if(phel != "m" && phel != "p"){
    cout<<"Proton helicity state \""<<ehel<<"\" invalid. Use \"p/m\"."<<endl;
    goodPol = false;
  }
  if(pdir != "L" && pdir != "T"){
    cout<<"Proton helicity state \""<<ehel<<"\" invalid. Use \"L/T\"."<<endl;
    goodPol = false;
  }

  return goodPol;
}

void DVCS_Binning(TString campaign = "26.07.1", TString energy = "9x130", TString ehel = "m", TString pdir = "T", TString phel = "m"){
  // Input validation on polarisation states
  if(!checkPol(ehel, pdir, phel)) return;
  
  TString sPolState = "e"+ehel+"h"+pdir+phel;

  // Print plot settings
  cout<<"---------------------------------"<<endl;
  cout<<"Creating binning scheme"<<endl;
  cout<<"\tCampaign: "<<campaign<<endl;
  cout<<"\tEnergy: "<<energy<<endl;
  cout<<"\tBeam polarisation state: "<<sPolState<<endl;
  cout<<"---------------------------------"<<endl;
  
  TString sIn = "$EIC_WORK_DIR/DVCS_Analysis/RootFiles/ePIC_DVCS_"+campaign+"_"+energy+"_"+sPolState+".root";
  TFile* fIn = TFile::Open(sIn);

  // Load histograms - expecting 3 Q2 bins
  double q2edges[4]  = {0.99, 1.4, 2.5, 100.};
  const int nQ2bins = sizeof(q2edges)/sizeof(q2edges[0]) - 1;

  TH2D* h_xvt_q2diff[3];

  // Extract bins for B0 t-range
  float t_gap{-1.};
  if(energy == "9x130") t_gap = 0.45;
  else if(energy == "9x275") t_gap = 1.9;
  float t_lo{0.01}; // V. low edge for BH singularity
  
  // Using TCutG to select ranges
  TCutG* cutIsB0 = new TCutG("isB0",4);
  cutIsB0->SetPoint(0, 0., t_gap);
  cutIsB0->SetPoint(1, 1., t_gap);
  cutIsB0->SetPoint(2, 1., 2.);
  cutIsB0->SetPoint(3, 0., 2.);

  // Calculate scaling factor to reach early science lumi
  double scaleToEIC{1.};
  if(energy == "9x130") scaleToEIC = calcScaling(energy,sPolState,1e15);
  else if(energy == "9x275") scaleToEIC = calcScaling(energy,sPolState,2.5e15);

  // Loop over number of Q2 bins
  for(int q{0}; q<3; q++){
    cout<<"--------------------------------------------------------"<<endl;
    cout<<"Q2 bin "<<q+1<<" - ["<<q2edges[q]<<" - "<<q2edges[q+1]<<"] GeV2\n"<<endl;
    
    h_xvt_q2diff[q] = (TH2D*)fIn->Get(Form("xvtdiff_rp[%i]",q));
    //h_xvt_q2diff[q]->Scale(scaleToEIC);

    //---------------------------------------------------------------------------------------------------
    // Plot to show bin edges
    TCanvas* c = new TCanvas("c","",900,900);
    gStyle->SetOptStat(00000000);
    // Plot 2D x/t histogram
    h_xvt_q2diff[q]->SetTitle(Form("x_{B}:|t|, %.2f<Q^{2}<%.2f GeV^{2}",q2edges[q],q2edges[q+1]));
    gPad->SetLogx();
    gPad->SetLogz();
    h_xvt_q2diff[q]->Draw("col");
    //---------------------------------------------------------------------------------------------------
    
    // Step 1 - Create temporary histograms for B0 and RP t-ranges
    TH2D* tempB0 = (TH2D*)h_xvt_q2diff[q]->Clone("tempb0");
    tempB0->Reset();
    TH2D* tempRP = (TH2D*)h_xvt_q2diff[q]->Clone("temprp");
    tempRP->Reset();

    for(int binx{1}; binx<tempB0->GetNbinsX()+1;binx++){
      for(int biny{1}; biny<tempB0->GetNbinsY()+1;biny++){
	float x = tempB0->GetXaxis()->GetBinCenter(binx);
	float y = tempB0->GetYaxis()->GetBinCenter(biny);

	if(cutIsB0->IsInside(x,y)) tempB0->Fill(x,y,h_xvt_q2diff[q]->GetBinContent(binx,biny)); //fi (IsInside TCutG - for B0)
	else tempRP->Fill(x,y,h_xvt_q2diff[q]->GetBinContent(binx,biny)); //fi (else - for RP)
      } //rof (y bins)
    }   //rof (x bins)
    

    // Step 2 - B0 bins (t-integrated)
    // Split into xB bins of integral at least 1000 (minimum 1)
    
    // Case 1 - all in 1 bin (< 1900 entries)
    if(tempB0->Integral() <= 1900){
      cout<<"(1) High-t (few bins) - \n|t| = {"<<t_gap<<", 2.}\nx = {0., 1.}"<<endl;
      
      //---------------------------------------------------------------------------------------------------
      // Add boxes to plots
      TBox* boxB0 = new TBox(0,t_gap,1,2.);
      boxB0->SetFillStyle(0);
      boxB0->SetLineColor(kRed);
      boxB0->SetLineWidth(2);
      boxB0->Draw();
      //---------------------------------------------------------------------------------------------------

    } //fi (integral of B0 histogram)
    
    // Case 2 - able to split
    else{
      // Project down along xB and split by entries
      TH1D* projx = tempB0->ProjectionX();
      vector<float> x_edges;
      splitByEntries_array(projx,1000,x_edges);
      // Print bins
      cout<<"(1) B0 region - \n|t| = {"<<t_gap<<", 2.}\nxB = {";
      for(int i{0}; i<x_edges.size();i++){
	if(i==x_edges.size()-1) cout<<x_edges[i]<<"}"<<endl;
	else cout<<x_edges[i]<<", ";

	//---------------------------------------------------------------------------------------------------
	// Add boxes to plots
	TBox* boxB0 = new TBox(x_edges[i],t_gap,x_edges[i+1],2.);
	boxB0->SetFillStyle(0);
	boxB0->SetLineColor(kRed);
	boxB0->SetLineWidth(2);
	boxB0->Draw();
	//---------------------------------------------------------------------------------------------------
	
      } //rof (No. of xB edges)
    }   //fi (else)


    // Step 3 - RP bins
    // Project along xB and split by bins of integral ~4k
    TH1D* projx = tempRP->ProjectionX();
    vector<float> x_edges;
    splitByEntries_array(projx,4000,x_edges);
    //splitByEntries_print(projx,4000);

    cout<<"(2) Low-t (fine binning) - \nxB = {";
    for(int j{0}; j<x_edges.size();j++){
      if(j==x_edges.size()-1) cout<<x_edges[j]<<"}"<<endl;
      else cout<<x_edges[j]<<", ";
    } //rof (No. of xB edges)
    
    for(int xb{0}; xb<x_edges.size()-1;xb++){
      // For each xB bin, ProjectionY along |t| and split into 4 t-bins
      TH1D* projy;
      // If last set of bins, just run to end of histogram
      if(xb==x_edges.size()-1) projy = tempRP->ProjectionY("",projx->FindBin(x_edges[xb]), projx->GetNbinsX());
      else projy = tempRP->ProjectionY("",projx->FindBin(x_edges[xb]), projx->FindBin(x_edges[xb+1])-1);

      vector<float> t_edges;
      splitByBins_array(projy,4,t_edges);
      //splitByBins_print(projy,4);
      
      // Overwrite elements - first bin just beyond BH singularity
      t_edges[0] = t_lo;
      // Overwrite elements - last bin goes to B0/RP gap
      t_edges[t_edges.size()-1] = t_gap;

      cout<<"xB: {"<<x_edges[xb]<<", "<<x_edges[xb+1]<<"}\t|t| = {";
      for(int k{0}; k<t_edges.size(); k++){
	if(k==t_edges.size()-1) cout<<t_edges[k]<<"}"<<endl;
	else cout<<t_edges[k]<<", ";
	
	//---------------------------------------------------------------------------------------------------
	// Add boxes to plots
	TBox* boxRP = new TBox(x_edges[xb],t_edges[k],x_edges[xb+1],t_edges[k+1]);
	boxRP->SetFillStyle(0);
	boxRP->SetLineColor(kRed);
	boxRP->SetLineWidth(2);
	boxRP->Draw();
	//---------------------------------------------------------------------------------------------------

      } // rof (No. of t edges)
    } //rof (No. of xB edges)
    
  } //rof (loop over Q2 bins)

  return;
}
