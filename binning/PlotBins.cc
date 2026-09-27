///----------------------------------------------------------------------------------------
///
/// Macro for plotting 3D DVCS bins
/// NEEDS: input DVCS analysis file (with xB:t distributions), .txt file with bin edges
///
///----------------------------------------------------------------------------------------

using namespace std;

#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

// Ragged Q2/xB/|t| binning reader - the SAME one used by the analysis.
#include "../include/DVCSBinning.hh"
// ePIC plotting style
#include "../plotting/ePIC_style.C"

// ePIC simulation campaign tag, beam energy and polarization state string need to match existing simulation file
void PlotBins(TString campaign = "26.07.1", TString energy = "9x130", TString pol = "emhTm"){

  // Load analysis ROOT file
  TString sIn = "../rootfiles/ePIC_DVCS_"+campaign+"_"+energy+"_"+pol+"_diff.root";
  TFile* fIn = TFile::Open(sIn);

  // Declare x:t histograms
  TH2D* h_xvt_q2diff[3];

  // Declare binning scheme object
  DVCSBinning binning;

  // Set .txt file for binning
  TString sBinFile = "./bins_"+energy+".txt";
  if(!binning.load(sBinFile.Data())){
    std::cerr << "[ePIC_DVCS] FATAL: could not load binning from " << sBinFile << std::endl;
    return;
  }
  binning.print();
  const int nQ2bins = binning.nQ2();

  // Set ePIC plotting style
  //gROOT->ProcessLine("set_ePIC_style()");
  //gStyle->SetCanvasPreferGL(kTRUE);
  
  // Run over Q2 bins - 1 canvas per bin
  for(int q{0}; q<nQ2bins; q++){

    // Load 2D histogram
    h_xvt_q2diff[q] = (TH2D*)fIn->Get(Form("xvtdiff_rp[%i]",q));

    //---------------------------------------------------------------------------------------------------
    // Create canvas and draw initial histogram
    TCanvas* c = new TCanvas("c","",900,900);
    gStyle->SetOptStat(00000000);
    // Plot 2D x/t histogram
    h_xvt_q2diff[q]->SetTitle(Form("x_{B}:|t|, %.2f<Q^{2}<%.2f GeV^{2}",binning.q2Low(q), binning.q2High(q)));
    gPad->SetLogx();
    gPad->SetLogz();
    gPad->SetRightMargin(0.12);
    h_xvt_q2diff[q]->GetYaxis()->SetTitle("|t| [GeV^{2}]");
    h_xvt_q2diff[q]->GetYaxis()->SetTitleOffset(1.3);
    h_xvt_q2diff[q]->GetXaxis()->SetTitle("x_{B}");
    h_xvt_q2diff[q]->Draw("colz");
    //---------------------------------------------------------------------------------------------------

    // Step 1: Draw bins for B0 region
    int reg{DVCSBinning::kB0};
    // Find number of x bins for region
    int nx = binning.nXB(reg,q);
    
    // Loop over x bins
    for(int x{0}; x<nx; x++){
      int nt = binning.nT(reg,q,x);
      // Loop over t bins
      for(int t{0}; t<nt; t++){
	// Create TBox based on edges of x/t bins
	TBox* boxB0 = new TBox(binning.xBLow(reg,q,x),binning.tLow(reg,q,x,t),binning.xBHigh(reg,q,x),binning.tHigh(reg,q,x,t));
	boxB0->SetFillStyle(0);
	boxB0->SetLineColor(kRed);
	boxB0->SetLineWidth(2);
	boxB0->Draw();
      } // EOL - t
    } // EOL - x
    
    
    // Step 2: bins for RP region
    reg = DVCSBinning::kRP;
    nx = binning.nXB(reg,q);
    // Loop over x bins
    for(int x{0}; x<nx; x++){
      int nt = binning.nT(reg,q,x);
      // Loop over t bins
      for(int t{0}; t<nt; t++){
	// Create TBox based on edges of x/t bins
	TBox* boxRP = new TBox(binning.xBLow(reg,q,x),binning.tLow(reg,q,x,t),binning.xBHigh(reg,q,x),binning.tHigh(reg,q,x,t));
	boxRP->SetFillStyle(0);
	boxRP->SetLineColor(kRed);
	boxRP->SetLineWidth(2);
	boxRP->Draw();
      } // EOL - t
    } // EOL - x

    // Add text to binning plot
    TLatex* tePICLabel = new TLatex(0.15, 0.85,
				    Form("ePIC #bf{Simulation %s, %s GeV}", campaign.Data(), energy.Data()));
    tePICLabel->SetNDC();
    tePICLabel->SetTextSize(0.04);
    tePICLabel->Draw("same");
    // Save plot
    c->SaveAs(Form("3Dbinning_%s_%i.png",energy.Data(),q));
    c->Close();
    
  } // EOL - Q2
  
  return;
}

