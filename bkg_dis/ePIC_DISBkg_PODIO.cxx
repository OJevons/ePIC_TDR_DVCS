// ePIC DVCS analysis class definition
// USE CASE - DIS PHYSICS BACKGROUND

#include "preLoadLib.hh"

// Data model headers
#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/utils/vector_utils.h"
#include "edm4hep/utils/kinematics.h"
#include "edm4eic/ClusterCollection.h"
#include "edm4eic/MCRecoParticleAssociationCollection.h"
#include "podio/Frame.h"
#include "podio/ROOTReader.h"

// ROOT Includes
#include <TSystem.h>
#include <TMath.h>
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <Math/LorentzVector.h>
#include <TTree.h>
#include <TString.h>

#include <fstream>

// Class header include
#include "../include/ePIC_DVCS_TASK.h"

//----------------------------------------------------
//----------------------------------------------------
//                    CONSTRUCTORS
//----------------------------------------------------
//----------------------------------------------------

// Default constructor
ePIC_DVCS_TASK::ePIC_DVCS_TASK(){
}

// Specific constructor
ePIC_DVCS_TASK::ePIC_DVCS_TASK(TString camp, TString energy, TString sett){
  // Set stored campaign attributes
  setDate(camp);
  setEnergy(energy);
  setSetting(sett);

  setBeamMomenta();
}


//----------------------------------------------------
//----------------------------------------------------
//                    SETTERS
//----------------------------------------------------
//----------------------------------------------------

// Set input file list
void ePIC_DVCS_TASK::setInFileList(TString name){
  sInList = name;
  std::cout<<"Input file list used: "<<name<<std::endl;
}

// Set output file name and create new
void ePIC_DVCS_TASK::setOutFile(TString name){
  std::cout<<"Output ROOT file: "<<name<<std::endl;
  fOutFile = new TFile(name,"RECREATE");
}

void ePIC_DVCS_TASK::setBeamMomenta(){
  if(sEnergy == "5x41"){
    fPBeam_p=41.0;
    fPBeam_e=5.0;
  }
  else if(sEnergy == "10x100"){
    fPBeam_p=100.0;
    fPBeam_e=10.0;
  }
  else if(sEnergy == "10x130"){
    fPBeam_p=130.0;
    fPBeam_e=10.0;
  }
  else if(sEnergy == "10x250"){
    fPBeam_p=250.0;
    fPBeam_e=10.0;
  }
  else if(sEnergy == "10x275"){
    fPBeam_p=275.0;
    fPBeam_e=10.0;
  }
  else if(sEnergy == "18x275"){
    fPBeam_p=275.0;
    fPBeam_e=18.0;
  }
  else{
    fPBeam_p=100.0;
    fPBeam_e=10.0;
  }
}

void ePIC_DVCS_TASK::setMomCutFactors(Float_t factore = 1.,Float_t factorp = 1. ){
  fPMaxFactor_e = factore;
  fPMaxFactor_p = factorp;
}

//----------------------------------------------------
//----------------------------------------------------
//                    APPLY CUTS
//----------------------------------------------------
//----------------------------------------------------

// Single particle cuts - electron
Bool_t ePIC_DVCS_TASK::applyCuts_Electron(P3EVector beame, std::vector<P3EVector> scate){
  Bool_t passCuts{kTRUE};
   
  // EVENT CUTS
  // Require single particle in final state
  if(scate.size() != 1) passCuts = kFALSE;
  //if(scate.size() == 0) passCuts = kFALSE;
  // Return out of function if array is not filled
  if(!passCuts) return passCuts;

  // KINEMATIC CUTS
  // 1. Momentum
  if(scate[0].P() > (fPBeam_e*fPMaxFactor_e)) passCuts = kFALSE;
  // 2. Q2
  fQ2 = calcQ2_Elec(beame, scate[0]);
  if(fQ2 < fMinQ2) passCuts = kFALSE;

  return passCuts;
}

// Single particle cuts - photon
Bool_t ePIC_DVCS_TASK::applyCuts_Photon(std::vector<P3EVector> scatg){
   Bool_t passCuts{kTRUE};
   
   // EVENT CUTS
   // Require single particle in final state
   if(scatg.size() != 1) passCuts = kFALSE;
   // Return out of function if array is not filled
   if(!passCuts) return passCuts;

   // KINEMATIC CUTS
   // Energy cut - beam bkg photons at low energy
   if(scatg[0].E() < 1) passCuts = kFALSE;

   return passCuts;
}

// Single particle cuts - proton
Bool_t ePIC_DVCS_TASK::applyCuts_Proton(std::vector<P3EVector> scatp, TString sProtonDet="all"){
  Bool_t passCuts{kTRUE};
  
  // EVENT CUTS
  // Require single particle in final state
  if(scatp.size() != 1) passCuts = kFALSE;
  // Return out of function if array is not filled
  if(!passCuts) return passCuts;
  
  // KINEMATIC CUTS
  // 1. Momentum
  if(scatp[0].P() > (fPBeam_p*fPMaxFactor_p)) passCuts = kFALSE;
  
  // 2. Scattered proton theta (ensure within B0, Roman Pots or 'all')
  // If invalid detector name used, consider all
  if(sProtonDet != "B0" && sProtonDet != "RP" && sProtonDet != "all") sProtonDet="all";
  Float_t fMinPTheta{0.};
  Float_t fMaxPTheta{0.};

  // B0 angular acceptance: 5.5 mrad - 20 mrad
  if(sProtonDet == "B0"){
    fMinPTheta = 0.0055;
    fMaxPTheta = 0.02;
  }
  // RP angluar acceptance: < 5.0 mrad
  else if(sProtonDet == "RP"){
    fMinPTheta = 0.;
    fMaxPTheta = 0.005;
  }
  // Full FF proton acceptance: < 20 mrad
  else if(sProtonDet == "all"){
    fMinPTheta = 0.;
    fMaxPTheta = 0.02;
  } 
  if(scatp[0].Theta()<fMinPTheta || scatp[0].Theta()>fMaxPTheta) passCuts = kFALSE;

  return passCuts;
}

// Event-level cuts (DVCS kinematics)
Bool_t ePIC_DVCS_TASK::applyCuts_DVCS(TString sProtonDet="all"){
  Bool_t passCuts{kTRUE};

  if(sProtonDet != "B0" && sProtonDet != "RP" && sProtonDet != "all") sProtonDet="all";
  
  // E-pz cut
  if(fEmPz < fMin_EmPz || fEmPz > fMax_EmPz) passCuts = kFALSE;
  // Missing pT cut
  if(fpTmiss > fMax_pTmiss) passCuts = kFALSE;

  return passCuts;
}

// Combination of all cuts
Bool_t ePIC_DVCS_TASK::applyCuts_All(P3EVector beame, P3EVector beamp, vector<P3EVector> scate, vector<P3EVector> scatp, vector<P3EVector> scatg, TString sProtonDet="all"){
  Bool_t passCuts{kTRUE};
  
  // 1. Electron cuts
  // Need to calculate Q2 first - set to zero if no detected electron
  if(scate.size() == 0) fQ2 = 0;
  else fQ2 = calcQ2_Elec(beame, scate[0]);
  passCuts = applyCuts_Electron(beame,scate);
  // Exit from function if failure
  if(!passCuts) return passCuts;

  // 2. Photon cuts
  passCuts = applyCuts_Photon(scatg);
  // Exit from function if failure
  if(!passCuts) return passCuts;

  // 3. Proton cuts (if not semi-inclusive)
  if(sProtonDet != "semi"){
    passCuts = applyCuts_Proton(scatp, sProtonDet);
    // Exit from function if failure
    if(!passCuts) return passCuts;
  }

  // 4. Event cuts
  if(sProtonDet == "semi"){
    // E-pz cut (general to fully exclusive and semi-inclusive reco.)
    fEmPz = (scate[0]+scatg[0]).E() - (scate[0]+scatg[0]).Pz();
    //if(fEmPz < fMin_EmPz || fEmPz > fMax_EmPz) passCuts = false;
    // Exit from function if failure
    //if(!passCuts) return passCuts;
    
    fpTmiss = calcPtMiss_2Body(beame, beamp, scate[0], scatg[0]);
    passCuts = applyCuts_DVCS(sProtonDet);
  }
  else{
    // E-pz cut (general to fully exclusive and semi-inclusive reco.)
    fEmPz = (scate[0]+scatg[0]+scatp[0]).E() - (scate[0]+scatg[0]+scatp[0]).Pz();
    //if(fEmPz < fMin_EmPz || fEmPz > fMax_EmPz) passCuts = false;
    // Exit from function if failure
    //if(!passCuts) return passCuts;

    // Last cuts - fully exclusive DVCS event kinematics
    fpTmiss = calcPtMiss_3Body(beame, beamp, scate[0], scatp[0], scatg[0]);
    passCuts = applyCuts_DVCS(sProtonDet);
  }

  return passCuts;
}

//----------------------------------------------------
//----------------------------------------------------
//            UNDO AFTERBURNER PROCEDURE
//----------------------------------------------------
//----------------------------------------------------

// Undo AB and calculate boost vectors - DO THIS FIRST FOR EACH EVENT
// USE BEAM VECTORS
void ePIC_DVCS_TASK::undoAfterburnAndCalc(P3EVector& p, P3EVector& k){
  // Correction Feb. '26 - Use nomimal beam momenta for ab calculation
  P3EVector p_beam(fXAngle*fPBeam_p, 0., fPBeam_p, fPBeam_p);
  P3EVector e_beam(0., 0., -fPBeam_e, fPBeam_e);
    
  // Define boost vector to CoM frame
  P3EVector CoM_boost = p_beam+e_beam;
  vBoostToCoM.SetXYZ(-CoM_boost.X()/CoM_boost.E(), -CoM_boost.Y()/CoM_boost.E(), -CoM_boost.Z()/CoM_boost.E());
  
  // Apply boost to beam vectors
  p_beam = boost(p_beam, vBoostToCoM);
  e_beam = boost(e_beam, vBoostToCoM);
  
  // Calculate rotation angles and create rotation objects
  fRotY = -1.0*TMath::ATan2(p_beam.X(), p_beam.Z());
  fRotX = 1.0*TMath::ATan2(p_beam.Y(), p_beam.Z());

  rotAboutY = RotationY(fRotY);
  rotAboutX = RotationX(fRotX);

  // Apply rotation to beam vectors
  p_beam = rotAboutY(p_beam);
  p_beam = rotAboutX(p_beam);
  e_beam = rotAboutY(e_beam);
  e_beam = rotAboutX(e_beam);

  // Define boost vector back to head-on frame
  P3EVector HoF_boost(0., 0., CoM_boost.Z(), CoM_boost.E());
  vBoostToHoF.SetXYZ(HoF_boost.X()/HoF_boost.E(), HoF_boost.Y()/HoF_boost.E(), HoF_boost.Z()/HoF_boost.E());

  // Apply boost back to head on frame to beam vectors
  p_beam = boost(p_beam, vBoostToHoF);
  e_beam = boost(e_beam, vBoostToHoF);

  // Make changes to input vectors
  p.SetPxPyPzE(p_beam.X(), p_beam.Y(), p_beam.Z(), calcE(p_beam.Vect(),fMass_proton));
  k.SetPxPyPzE(e_beam.X(), e_beam.Y(), e_beam.Z(), calcE(e_beam.Vect(),fMass_electron));
}

// Undo afterburn procedure only
void ePIC_DVCS_TASK::undoAfterburn(P3EVector& a){
  Float_t mass = a.M();
  
  // Undo AB procedure for single vector, a^{mu}
  a = boost(a, vBoostToCoM); // BOOST TO COM FRAME
  a = rotAboutY(a);          // ROTATE TO Z-AXIS
  a = rotAboutX(a);          // ROTATE TO Z-AXIS
  a = boost(a, vBoostToHoF); // BOOST BACK TO HEAD ON FRAME

  a.SetPxPyPzE(a.X(), a.Y(), a.Z(), calcE(a.Vect(),mass));
}

//----------------------------------------------------
//----------------------------------------------------
//              KINEMATIC CALCULATIONS
//----------------------------------------------------
//----------------------------------------------------

// Alternate Trento Phi calculation using particle 4-vectors
// Source: Bachetta, A. et al; JHEP 02 (2007); eq. 2.3
Double_t ePIC_DVCS_TASK::calcTrentoPhi_4Vec(P3EVector k, P3EVector p, P3EVector kprime, P3EVector pprime){
  // First, extract 4-vector components into arrays
  Float_t ppcomp[4]{0.,0.,0.,0.};
  Float_t kcomp[4]{0.,0.,0.,0.};
  pprime.GetCoordinates(ppcomp);
  k.GetCoordinates(kcomp);
  
  // Need q vector for gT calculations
  P3EVector q = k-kprime;
  
  // Need kT^2 and p'T^2 for denominator
  Float_t kTcomp[4]{0.,0.,0.,0.};
  Float_t ppTcomp[4]{0.,0.,0.,0.};
  
  // Loop over indices of 4-vectors
  Float_t numCos{0};
  Float_t numSin{0};
  
  for(int i{0}; i<4; i++){
    for(int j{0}; j<4; j++){
      kTcomp[i] += calcgT_ij(q,p,i,j)*kcomp[j];
      ppTcomp[i] += calcgT_ij(q,p,i,j)*ppcomp[j];

      numCos += kcomp[i]*ppcomp[j]*calcgT_ij(q,p,i,j);
      numSin += kcomp[i]*ppcomp[j]*calcepsT_ij(q,p,i,j);
    } // Loop over index j
  } // Loop over index i

  P3EVector kT(kTcomp[0],kTcomp[1],kTcomp[2],kTcomp[3]);
  P3EVector ppT(ppTcomp[0],ppTcomp[1],ppTcomp[2],ppTcomp[3]);

  Float_t kTSq = (kTcomp[3]*kTcomp[3])-(kTcomp[0]*kTcomp[0])-(kTcomp[1]*kTcomp[1])-(kTcomp[2]*kTcomp[2]);
  Float_t ppTSq = (ppTcomp[3]*ppTcomp[3])-(ppTcomp[0]*ppTcomp[0])-(ppTcomp[1]*ppTcomp[1])-(ppTcomp[2]*ppTcomp[2]);

  Float_t den = TMath::Sqrt(kTSq*ppTSq);
  Float_t cosphi = -numCos/den;
  Float_t sinphi = -numSin/den;

  return TMath::ACos(cosphi)*TMath::Sign(1.,TMath::ASin(sinphi));
}

// Calculate value of transverse metric tensor for 4-vector phi calculation
// Source: Bachetta, A. et al; JHEP 02 (2007); eq. 2.4
Double_t ePIC_DVCS_TASK::calcgT_ij(P3EVector q, P3EVector p, Int_t i, Int_t j){
  // Extract 4-vector components into arrays
  Float_t pcomp[4]{0.,0.,0.,0.};
  Float_t qcomp[4]{0.,0.,0.,0.};
  p.GetCoordinates(pcomp);
  q.GetCoordinates(qcomp);

  // First term comes from standard metric tensor
  Int_t gij{0};
  
  // Metric tensor only contains diagonal components
  if(i==j && i==3) gij = 1;
  else if(i==j) gij = -1;
  
  // Second term comes from cross-product of 4-vector terms
  Float_t fQ2 = -q.M2();
  Float_t fpq = p.Dot(q);
  Float_t fGamma2 = (TMath::Power(fMass_proton,2)*fQ2)/TMath::Power(fpq,2);
  Float_t fPQTerm = qcomp[i]*pcomp[j] + pcomp[i]*qcomp[j];
  fPQTerm /= (1+fGamma2)*fpq;
  
  // Final term comes from products of terms in the same vector
  Float_t fQQTerm = (qcomp[i]*qcomp[j])/fQ2;
  fQQTerm -= (pcomp[i]*pcomp[j])/p.M2();
  fQQTerm *= fGamma2/(1+fGamma2);

  Float_t gT = gij - fPQTerm + fQQTerm;
  return gT;
}

// Calculate value of transverse Levi-Civita tensor for 4-vector phi calculation
// Source: Bachetta, A. et al; JHEP 02 (2007); eq. 2.4
Double_t ePIC_DVCS_TASK::calcepsT_ij(P3EVector q, P3EVector p, Int_t i, Int_t j){
  // Extract 4-vector components into arrays
  Float_t pcomp[4]{0.,0.,0.,0.};
  Float_t qcomp[4]{0.,0.,0.,0.};
  p.GetCoordinates(pcomp);
  q.GetCoordinates(qcomp);

  // Need gamma2 and p.q for calculation
  Float_t fQ2 = -q.M2();
  Float_t fpq = p.Dot(q);
  Float_t fGamma2 = (TMath::Power(fMass_proton,2)*fQ2)/TMath::Power(fpq,2);

  Float_t epsT{0};
  
  // Sum over all combinations of k and l
  for(int k{0}; k<4; k++){
    for(int l{0}; l<4; l++){
      epsT += LeviCivita(i,j,k,l)*(1./fpq)*(1./TMath::Sqrt(1+fGamma2))*pcomp[k]*qcomp[l];
    }
  }

  return epsT;
}

// Need Levi-Civita symbol for epsT calculation
// Using form of epsilon using Pi notation
Int_t ePIC_DVCS_TASK::LeviCivita(int i, int j, int k, int l){
  int ai[4]{i+1,j+1,k+1,l+1};

  int eps{1};

  for(int ind1{3}; ind1>0; ind1--){
    for(int ind2{0}; ind2<ind1; ind2++){
      // epsilon = 0 if any 2 indices are the same
      if(ai[ind1] == ai[ind2]) eps*=0;
      else eps*=TMath::Sign(1.,ai[ind1]-ai[ind2]);
    }
  }
  return eps;
}

// Calculate angle between planes of qp and qg
Double_t ePIC_DVCS_TASK::calcPhiQPQG(P3EVector k, P3EVector p, P3EVector kprime, P3EVector gprime){
  MomVector p3 = p.Vect();
  MomVector k3 = k.Vect();
  MomVector kp3 = kprime.Vect();
  MomVector gp3 = gprime.Vect();
  MomVector q3 = k3-kp3;

  MomVector qp = q3.Cross(p3);
  MomVector qg = q3.Cross(gp3);

  return Angle(qp,qg);
}

// Calculation of cone angle (angle between measured photon and expected photon)
Double_t ePIC_DVCS_TASK::calcConeAngle(P3EVector k, P3EVector p, P3EVector kprime, P3EVector pprime, P3EVector gprime){
  // Initial state vectors
  MomVector p3 = p.Vect(); // Proton beam
  MomVector k3 = k.Vect(); // Electron beam
  MomVector pi = p3+k3;    // Total initial momenta
  // Final state vectors
  MomVector pp3 = pprime.Vect(); // Scattered proton
  MomVector kp3 = kprime.Vect(); // Scattered electron
  MomVector gp3 = gprime.Vect(); // Real photon
  MomVector pfkp = pp3+kp3;      // Sum of scattered proton and scattered electron momenta

  // Return angle between real photon and expected photon (calc. from missing momentum of ep system)
  MomVector gExpected = pi-pfkp;

  return Angle(gp3,gExpected);
}

//----------------------------------------------------
//----------------------------------------------------
//                     DO ANALYSIS
//----------------------------------------------------
//----------------------------------------------------

void ePIC_DVCS_TASK::doAnalysis(){

  //---------------------------------------------------------
  // Setup: Load input file list
  //---------------------------------------------------------
  ifstream fileListStream;
  fileListStream.open(sInList);
  string fileName;
  TFile* inputRootFile;
  
  //---------------------------------------------------------
  // Setup: Declare histograms
  //---------------------------------------------------------
  // Starting with TDR histograms
  // 1a) Eta - MC particles
  TH1D* h_eta_MCp   = new TH1D("eta_MCp",";#eta_{p'}(MC)", 275, -11.0, 11.0);
  TH1D* h_eta_MCe   = new TH1D("eta_MCe",";#eta_{e'}(MC)", 275, -11.0, 11.0);
  TH1D* h_eta_MCg   = new TH1D("eta_MCg",";#eta_{#gamma}(MC)", 275, -11.0, 11.0);
  // 1c) Eta - reco. particles
  TH1D* h_eta_RPp   = new TH1D("eta_RPp",";#eta_{p'}(Reco)", 275, -11.0, 11.0);
  TH1D* h_eta_RPPp   = new TH1D("eta_RPPp",";#eta_{p'}(Reco)", 275, -11.0, 11.0);
  TH1D* h_eta_RPe   = new TH1D("eta_RPe",";#eta_{e'}(Reco)", 275, -11.0, 11.0);
  TH1D* h_eta_RPg   = new TH1D("eta_RPg",";#eta_{#gamma}(Reco)", 275, -11.0, 11.0);
  // 1c) E/p for electron candidates (before cuts applied)
  TH1D* h_EoverP_elec = new TH1D("eoverp_elec",";E/p",100,0.5,1.5);
  // 1d) E/eta coverage
  TH2D* h_2D_EvEta_g = new TH2D("2d_eveta_g",";#eta_{#gamma};E_{#gamma} [GeV]",200,-4.,4.,100,0.,50.);
  TH2D* h_2D_EvEta_e = new TH2D("2d_eveta_e",";#eta_{e'};E_{e'} [GeV]",200,-4.,4.,(Int_t)4*fPBeam_e, 0., 2.*fPBeam_e);
  TH2D* h_2D_EvEta_p = new TH2D("2d_eveta_p",";#eta_{p'};E_{p'} [GeV]",150,4.,10.,(Int_t)4.*fPBeam_p, 0., 2.*fPBeam_p);

  // 2) t distribution
  TH1D* h_t_Truth  = new TH1D("t_truth" ,";|t|(MC) [(GeV/#it{c}^{2})^{2}]"           , 20, 0., 2.);
  TH1D* h_t_B0Acc  = new TH1D("t_b0acc" ,";|t|(MC|Reco. - B0) [(GeV/#it{c}^{2})^{2}]", 20, 0., 2.);
  TH1D* h_t_RPAcc  = new TH1D("t_rpacc" ,";|t|(MC|Reco. - RP) [(GeV/#it{c}^{2})^{2}]", 20, 0., 2.);
  TH1D* h_t_B0Reco = new TH1D("t_b0reco",";|t|(Reco. - B0) [(GeV/#it{c}^{2})^{2}]"   , 20, 0., 2.);
  TH1D* h_t_RPReco = new TH1D("t_rpreco",";|t|(Reco. - RP) [(GeV/#it{c}^{2})^{2}]"   , 20, 0., 2.);
  TH1D* h_t_LCAcc  = new TH1D("t_lcacc" ,";|t_{e'#gamma}|(MC|Reco.) [(GeV/#it{c}^{2})^{2}]", 20, 0., 2.);
  TH1D* h_t_LCReco = new TH1D("t_lcreco",";|t_{e''gamma}|(Reco.) [(GeV/#it{c}^{2})^{2}]"   , 20, 0., 2.);

  // 3) Counters
  TH1D* hCount = new TH1D("count",";;",2,0,2);
  hCount->GetXaxis()->SetBinLabel(1,"Not_epg");
  hCount->GetXaxis()->SetBinLabel(2,"Is_epg");
  const int nCuts{17};
  TString cutname[17] = {"inc-sing_ele","e'-Q^{2}",
                         "inc-sing_pho","#gamma-E_{#gamma}",
                         "inc-sing_pro","p'-#theta_{p}",
                         "e'p'#gamma-mult","e'p'#gamma-NegTrackVeto","e'p'#gamma-FFVeto","e'p'#gamma-(E-pz)","e'p'#gamma-p_{T,miss}"
                         "e'#gamma-mult","e'#gamma-NegTrackVeto","e'#gamma-FFVeto","e'#gamma-BarrelPosVeto","e'#gamma-(E-pz)","e'#gamma-p_{T,miss}"};
  TH1D* hPassCuts_MC = new TH1D("passcuts_mc",";;",nCuts,0,nCuts);
  TH1D* hPassCuts_Rec = new TH1D("passcuts_rec",";;",nCuts,0,nCuts);
  
  for(int bin{1}; bin<=nCuts; bin++){
    hPassCuts_MC->GetXaxis()->SetBinLabel(bin,cutname[bin-1]);
    hPassCuts_Rec->GetXaxis()->SetBinLabel(bin,cutname[bin-1]);
  }

  // 4) 2-particle exclusivity variables
  // Missing mass, momentum and energy
  Float_t protlower = (Float_t)0.8*fPBeam_p;
  Float_t protupper = (Float_t)1.2*fPBeam_p;
  Int_t nbins = protupper-protlower;
  
  TH1D* h_MM2eg_Rec = new TH1D("mm2eg_rec",";M^{2}_{miss, e'#gamma}(Reco.) [GeV^{2}]",400,-100,100);
  TH1D* h_PMeg_Rec  = new TH1D("pmeg_rec",";P_{miss, e'#gamma}(Reco.) [GeV]",5*nbins,protlower,protupper);
  TH1D* h_PTMeg_Rec  = new TH1D("ptmeg_rec",";P_{T, miss, e'#gamma}(Reco.) [GeV]",300,-5.,10.);
  TH1D* h_EMeg_Rec  = new TH1D("emeg_rec",";E_{miss, e'#gamma}(Reco.) [GeV]",800,-50,350);
  
  TH1D* h_PTMepg_Rec  = new TH1D("ptmepg_rec",";P_{T, miss, e'p'#gamma}(Reco.) [GeV]",300,-5.,10.);

  //---------------------------------------------------------
  // Loop over files in list
  //---------------------------------------------------------
  int fileCounter{0};

  // 4-vectors for beam particles - need these defined outside of file loop
  P3EVector beame4(0,0,0,-1);     // Beam electron (generated)
  P3EVector beamp4(0,0,0,-1);     // Beam proton (generated)
  
  // Start file loop
  while(getline(fileListStream,fileName)){
    std::cout<<"Input file "<<fileCounter<<" : "<<fileName<<std::endl;
    
    // Open podio reader
    // New reader for each file
    auto reader = podio::ROOTReader();
    reader.openFile(fileName);
    Int_t nEntries = reader.getEntries("events");
    std::cout<<"File has "<<nEntries<<" events..."<<std::endl;
    
    // Case of taking average beams from file
    if(!kUseEventBeams){
      // MUST DO THIS FIRST
      // Full run over tree in first file before anything else
      // Calculate beams from average of individual event beam particles
      if(fileCounter==0){
	// Accumulator variables
	P3EVector beame4_acc(0,0,0,-1);
	P3EVector beamp4_acc(0,0,0,-1);
	
	for(size_t ev = 0; ev < reader.getEntries("events"); ev++){
	  const auto event = podio::Frame(reader.readNextEntry("events"));
	  
	  // Define holding particles
	  edm4hep::MCParticle beame_evt;
	  edm4hep::MCParticle beamp_evt;

	  // LOOP AND FIND PARTICLE OBJECTS
	  // MCParticles
	  auto& mcparts = event.get<edm4hep::MCParticleCollection>("MCParticlesHeadOnFrameNoBeamFX");
	  for(const auto& mcp : mcparts){
	    if(mcp.getPDG() == 11 && mcp.getGeneratorStatus() == 4) beame_evt = mcp;
	    if(mcp.getPDG() == 2212 && mcp.getGeneratorStatus() == 4) beamp_evt = mcp;
	  } // BEAM PARTICLES FOR EVENT FOUND
	  
	  // Add to accumulator
	  XYZVector e3vec_temp(beame_evt.getMomentum().x,beame_evt.getMomentum().y,beame_evt.getMomentum().z);
	  P3EVector e4vec_temp(e3vec_temp.X(), e3vec_temp.Y(), e3vec_temp.Z(), calcE(e3vec_temp,fMass_electron));
	  XYZVector p3vec_temp(beamp_evt.getMomentum().x,beamp_evt.getMomentum().y,beamp_evt.getMomentum().z);
	  P3EVector p4vec_temp(p3vec_temp.X(), p3vec_temp.Y(), p3vec_temp.Z(), calcE(p3vec_temp,fMass_proton));

	  beame4_acc += e4vec_temp;
	  beamp4_acc += p4vec_temp;	  
	} // End of event loop - FOR AVERAGED BEAMS
	
	// Divide by number of events in file
	beame4.SetCoordinates(beame4_acc.X()/nEntries, beame4_acc.Y()/nEntries, beame4_acc.Z()/nEntries, beame4_acc.E()/nEntries);
	beamp4.SetCoordinates(beamp4_acc.X()/nEntries, beamp4_acc.Y()/nEntries, beamp4_acc.Z()/nEntries, beamp4_acc.E()/nEntries);

	undoAfterburnAndCalc(beamp4,beame4);

	std::cout<<"First file - beams\n\te:"<<beame4<<"\n\tp:"<<beamp4<<std::endl;
      } // fi (fileCounter == 0)
      else std::cout<<"Using beams from first file."<<std::endl;
    } // fi (!kUseEventBeams)

    // (Re)Run reader for main events
    for(size_t ev = 0; ev < reader.getEntries("events"); ev++){
      // Load next event
      //const auto event = podio::Frame(reader.readNextEntry("events"));
      const auto event = podio::Frame(reader.readEntry("events",ev));
    
      // Define booleans for vetoes
      bool kBarrelPos_MC{false};      // Positive track outside FF region?
      bool kOtherFFTrack_MC{false};   // Track in OMD/ZDC?
      bool kNonElecNeg_MC{false};     // Non-electron -ve tracks in barrel
      bool kBarrelPos_Rec{false};
      bool kOtherFFTrack_Rec{false};
      bool kNonElecNeg_Rec{false};
      
      // 4-vectors for MC raw particles
      vector<P3EVector> scate4_gen;   // Scattered electron (generated)
      vector<P3EVector> scatp4_gen;   // Scattered proton (generated)
      vector<P3EVector> scatg4_gen;   // Scattered photon (generated)
      // 4-vectors for associated MC particles (ONLY SCATTERED)
      vector<P3EVector> scate4_aso;   // Scattered electron (associated MC)
      vector<P3EVector> scatp4_aso;   // Scattered proton (associated MC)
      vector<P3EVector> scatg4_aso;   // Scattered photon (associated MC)
      // 4-vectors for reconstructed particles (SEPARATE PROTONS FOR B0 AND ROMAN POTS)
      vector<P3EVector> scate4_rec;   // Scattered electron (reconstructed)
      vector<P3EVector> scatp4_rec;   // Scattered proton (B0 reconstructed)
      vector<P3EVector> scatp4_rom;   // Scattered proton (Roman Pots reconstructed)
      vector<P3EVector> scatg4_rec;   // Scattered photon (reconstructed)     
      
      // MC truth
      auto& mcparts = event.get<edm4hep::MCParticleCollection>("MCParticlesHeadOnFrameNoBeamFX");
      for(const auto& mcp : mcparts){
	// If using beams per event, look for generatorStatus 4
	if(kUseEventBeams){
	  if(mcp.getPDG() == 11 && mcp.getGeneratorStatus() == 4){
	    beame4.SetCoordinates(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z, 
				  calcE(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z, fMass_electron));
	  }
	  if(mcp.getPDG() == 2212 && mcp.getGeneratorStatus() == 4){
	    beamp4.SetCoordinates(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z, 
				  calcE(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z, fMass_proton));
	  }
	  undoAfterburnAndCalc(beamp4,beame4);
	} // fi (kUseEventBeams)
	
	// Then look for rest of MC particles
	if(mcp.getGeneratorStatus() == 1){
	  P3EVector temp(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z, 
			 calcE(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z, mcp.getMass()));
	  
	  if(mcp.getPDG() == 11)   scate4_gen.push_back(temp);
	  if(mcp.getPDG() == 22)   scatg4_gen.push_back(temp);
	  if(mcp.getPDG() == 2212) scatp4_gen.push_back(temp);

	  // VETOES
	  // Check for positive track in barrel
	  if(mcp.getCharge() == 1 && TMath::Abs(temp.Eta()) < 4) kBarrelPos_MC = true;  // Positive track at -4 < eta < 4
	  // Check for particles that would hit OMD/ZDC
	  if(mcp.getCharge() == 0 && temp.Eta() > 4) kOtherFFTrack_MC = true;           // Neutral at high eta - ZDC
	  if(temp.Eta() > 4 && temp.Pz()/temp.E() < 0.65) kOtherFFTrack_MC = true;      // Particle at eta > 4, xL < 0.6 - OMD
	  // Check for non-electron -ve particles
	  if(mcp.getCharge() == -1 && mcp.getPDG() != 11) kNonElecNeg_MC = true;        // Negative particle, not electron

	} // fi (mcp.getGeneratorStatus() == 1)
      } // END OF MCPARTICLES LOOP

      // Reconstructed and associated particles (electrons/photons)
      const auto& assocReco = event.get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedParticleAssociations");
      for(const auto& mcreco : assocReco){
	// Declare holding vector variables
	P3EVector temp_mca(0,0,0,0);
	P3EVector temp_rec(0,0,0,0);
	
	// CASE 1: Using explicit MC matching
	// Only fill arrays if associated MC exists
	if(kUseExplicitMatch){
	  // Electrons - from sim.getPDG() flag
	  if(mcreco.getSim().getGeneratorStatus() == 1 && mcreco.getSim().getPDG() == 11){
	    // Simulated mass comes from MC truth - can trust sim.getMass()
	    temp_mca.SetCoordinates(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z,
				    calcE(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z, mcreco.getSim().getMass()));
	    // Reconstructed mass comes from PID hypothesis - cannot trust rec.getMass()
	    // Set mass by hand
	    temp_rec.SetCoordinates(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z,
				    calcE(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z, fMass_electron));
	  
	    // Undo afterburner
	    undoAfterburn(temp_mca);
	    undoAfterburn(temp_rec);
	    
	    // Add to particle arrays
	    scate4_aso.push_back(temp_mca); 
	    scate4_rec.push_back(temp_rec); 
	  }
	  // Photons - from sim.getPDG() flag
	  else if(mcreco.getSim().getGeneratorStatus() == 1 && mcreco.getSim().getPDG() == 22){
	    // Set LorentzVector coordinates by hand - as for electron case (but for massless reco.)
	    temp_mca.SetCoordinates(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z,
				    calcE(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z, mcreco.getSim().getMass()));
	    temp_rec.SetCoordinates(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z,
				    calcE(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z, 0.));
	  
	    // Undo afterburner
	    undoAfterburn(temp_mca);
	    undoAfterburn(temp_rec);
	 
	    // Add to particle arrays
	    scatg4_aso.push_back(temp_mca); 
	    scatg4_rec.push_back(temp_rec); 
	  }
	} // fi (kUseExplicitMatch)
	
	// Not using explicit matching?
	// Use reconstructed particle properties for PID
	else{
	  // CASE 2: Using ePIC PID
	  if(kUsePID){
	    // PID CODE HERE
	    // LOOKING FOR ELECTRONS AND PHOTONS
	  } // fi (kUsePID)
	  
	  // CASE 3: Using other properties for PID (charge, E/p, etc.)
	  else{
	    // Look for electrons - start from Q = -1
	    if(mcreco.getRec().getCharge() == -1){
	      // Skip if particle is missing ECAL clusters
	      if(!mcreco.getRec().getClusters()) continue;
	      
	      // Choose electrons from E/p 0.8 - 1.2
	      float clus_e{0};
	      for(auto& clust : mcreco.getRec().getClusters()) clus_e += clust.getEnergy();
	      float eoverp = clus_e/edm4hep::utils::magnitude(mcreco.getRec().getMomentum());
	      
	      h_EoverP_elec->Fill(eoverp);

	      // Electron finder cuts
	      if(eoverp < 0.8 || eoverp > 1.2){
		// VETO - any negative tracks that aren't electrons
		kNonElecNeg_Rec = true;
		continue;
	      }
	      
	      // IF WANTING TO ADD CLUSTER ISOLATION - DO HERE
	      // OTHER ELECTRON FINDER CUTS
	      
	      // Create Lorentz vectors
	      // Simulated mass comes from MC truth - can trust sim.getMass()
	      temp_mca.SetCoordinates(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z,
				      calcE(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z, mcreco.getSim().getMass()));
	      // Reconstructed mass comes from PID hypothesis - cannot trust rec.getMass()
	      // Set mass by hand
	      temp_rec.SetCoordinates(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z,
				      calcE(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z, fMass_electron));
	      
	      // Undo afterburner
	      undoAfterburn(temp_mca);
	      undoAfterburn(temp_rec);
	      
	      // Add to particle arrays
	      scate4_aso.push_back(temp_mca); 
	      scate4_rec.push_back(temp_rec);
	    } // End of electron finding
	    
	    // Look for photons - Q = 0 
	    // Will also catch cases where MC electron is missing track, but not clusters (clusters w/o track reconstructs as neutral)
	    if(mcreco.getRec().getCharge() == 0){
	      // Set LorentzVector coordinates by hand - as for electron case (but for massless reco.)
	      temp_mca.SetCoordinates(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z,
				      calcE(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z, mcreco.getSim().getMass()));
	      temp_rec.SetCoordinates(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z,
				      calcE(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z, 0.));
	      
	      // Undo afterburner
	      undoAfterburn(temp_mca);
	      undoAfterburn(temp_rec);
	      
	      // Add to particle arrays
	      scatg4_aso.push_back(temp_mca); 
	      scatg4_rec.push_back(temp_rec);
	    } // End of neutral finding
	  } // fi (!kUsePID)
	  
	} // fi (!kUseExplicitMatch)
	
	// VETOES
	// Positive tracks in barrel (B0 tracks ONLY in `TruthSeededCharged`)
	if(mcreco.getRec().getCharge() == 1) kBarrelPos_Rec = true;
      } // End of ReconstructedParticleAssociations
      

      // Now look for protons
      // Using ReconstructedTruthSeededChargedParticles
      const auto& assocTSReco = event.get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedTruthSeededChargedParticleAssociations");
      for(const auto& mcreco : assocTSReco){
	// Declare holding vector variables
	P3EVector temp_mca(0,0,0,0);
	P3EVector temp_rec(0,0,0,0);
      
	// CASE 1: Using explicit MC matching
	// Only fill arrays if associated MC exists
	if(kUseExplicitMatch){
	  // Select protons from sim.getPDG() flag
	  if(mcreco.getSim().getGeneratorStatus() == 1 && mcreco.getSim().getPDG() == 2212){
	    // Simulated mass comes from MC truth - can trust sim.getMass()
	    temp_mca.SetCoordinates(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z,
				    calcE(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z, mcreco.getSim().getMass()));
	    // Reconstructed mass comes from PID hypothesis - cannot trust rec.getMass()
	    // Set mass by hand
	    temp_rec.SetCoordinates(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z,
				    calcE(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z, fMass_proton));
	  
	    // Undo afterburner
	    undoAfterburn(temp_mca);
	    undoAfterburn(temp_rec);
	    
	    // Add to particle arrays
	    scatp4_aso.push_back(temp_mca); 
	    scatp4_rec.push_back(temp_rec); 
	  } // fi (sim.PDG() == 2212)
	} // fi (kUseExplicitMatch)

	// Not using explicit matching?
	// Use reconstructed particle properties for PID
	else{
	  // CASE 2: Using ePIC PID
	  if(kUsePID){
	    // PID CODE HERE
	    // LOOKING FOR ELECTRONS AND PHOTONS
	  } // fi (kUsePID)

	  // CASE 3: Using charge for PID
	  // Will ensure particle has associated track
	  else{
	    if(mcreco.getRec().getCharge() == 1){
	      // Simulated mass comes from MC truth - can trust sim.getMass()
	      temp_mca.SetCoordinates(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z,
				      calcE(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z, mcreco.getSim().getMass()));
	      // Reconstructed mass comes from PID hypothesis - cannot trust rec.getMass()
	      // Set mass by hand
	      temp_rec.SetCoordinates(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z,
				      calcE(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z, fMass_proton));
	  
	      // Undo afterburner
	      undoAfterburn(temp_mca);
	      undoAfterburn(temp_rec);

	      // Add to particle arrays
	      scatp4_aso.push_back(temp_mca); 
	      scatp4_rec.push_back(temp_rec); 
	    } // End of Q = 1
	  } // fi (!kUsePID)

	} // fi (!kUseExplicitMatch)
      } // End of ReconstructedTruthSeededChargedParticlesAssociations

      // Now get Roman Pot tracks
      // Branches from edm4eic::ReconstructedParticle without MC associations
      const auto& RPparts = event.get<edm4eic::ReconstructedParticleCollection>("ForwardRomanPotRecParticles");
      for(const auto& rpreco : RPparts){
	// Declare holding vector
	P3EVector temp_rec(0,0,0,0);
	
	// Assume all Roman Pot tracks are from reconstructed protons
	temp_rec.SetCoordinates(rpreco.getMomentum().x, rpreco.getMomentum().y, rpreco.getMomentum().z,
				calcE(rpreco.getMomentum().x, rpreco.getMomentum().y, rpreco.getMomentum().z, fMass_proton));

	// No need to undo afterburner - this is built in to the track reconstruction
	// Add to particle arrays
	scatp4_rom.push_back(temp_rec);
      } // End of ForwardRomanPotRecParticles

      // VETOES
      // OMD
      const auto& OMDparts = event.get<edm4eic::ReconstructedParticleCollection>("ForwardOffMRecParticles");
      for(const auto& omdreco : OMDparts){
	// Check if track exists
	if(omdreco.getTracks()) kOtherFFTrack_Rec = true;
      }
      // ZDC
      const auto& ZDCparts = event.get<edm4eic::ReconstructedParticleCollection>("ReconstructedFarForwardZDCNeutrals");
      for(const auto& zdcreco : ZDCparts){
	// Check for associated clusters
	if(zdcreco.getClusters()) kOtherFFTrack_Rec = true;
      }
      

      //---------------------------------------------------------
      // Fill histograms
      //---------------------------------------------------------

      // Look at all MC
      for(int ele_ind{0}; ele_ind<scate4_gen.size(); ele_ind++){
	h_eta_MCe->Fill(scate4_gen[ele_ind].Eta());
	h_2D_EvEta_e->Fill(scate4_gen[ele_ind].Eta(), scate4_gen[ele_ind].E());
      }
      for(int pho_ind{0}; pho_ind<scatg4_gen.size(); pho_ind++){
	h_eta_MCg->Fill(scatg4_gen[pho_ind].Eta());
	h_2D_EvEta_g->Fill(scatg4_gen[pho_ind].Eta(), scatg4_gen[pho_ind].E());
      }
      for(int pro_ind{0}; pro_ind<scatp4_gen.size(); pro_ind++){
	h_eta_MCp->Fill(scatp4_gen[pro_ind].Eta());
	h_2D_EvEta_p->Fill(scatp4_gen[pro_ind].Eta(), scatp4_gen[pro_ind].E());
      }
      
      //cout<<"[DEBUG]: ALL MC FILLED"<<endl;

      // Count number of generated e'p'gamma events from event cuts
      // Flag event if so (only want to look at fake DVCS events)
      bool kGenDVCS{false};
      if(applyCuts_All(beame4, beamp4, scate4_gen, scatp4_gen, scatg4_gen, "all") && kNonElecNeg_MC){
	hCount->Fill(1);
	kGenDVCS = true;
      }
      else hCount->Fill(0);
      
      //cout<<"[DEBUG]: EVENT COUNTER FILLED"<<endl;
      
      // Look at all reco. for eta
      for(int ele_ind{0}; ele_ind<scate4_rec.size(); ele_ind++)	h_eta_RPe->Fill(scate4_rec[ele_ind].Eta());  // Electrons
      for(int pho_ind{0}; pho_ind<scatg4_rec.size(); pho_ind++)	h_eta_RPg->Fill(scatg4_rec[pho_ind].Eta());  // Photons
      for(int pro_ind{0}; pro_ind<scatp4_rec.size(); pro_ind++)	h_eta_RPp->Fill(scatp4_rec[pro_ind].Eta());  // B0 protons
      for(int pro_ind{0}; pro_ind<scatp4_rom.size(); pro_ind++)	h_eta_RPPp->Fill(scatp4_rom[pro_ind].Eta()); // RP tracks - assume proton

      //cout<<"[DEBUG]: ALL RECO FILLED"<<endl;

      // Count no. of events which pass cuts - MC only
      // ---------------------Single particle---------------------
      // Electron
      if(scate4_gen.size() == 1 && !kGenDVCS){
        hPassCuts_MC->Fill(0);

        // ...and Q2 > 1 GeV2
	fQ2 = calcQ2_Elec(beame4, scate4_gen[0]);
        if(fQ2 >= 1.) hPassCuts_MC->Fill(1); //fi (Q2 cut)
      } //fi (electrons inclusive)
      // Photon
      if(scatg4_gen.size() == 1 && !kGenDVCS){
        hPassCuts_MC->Fill(2);

        // ...and E_{gamma} > 1 GeV
	if(scatg4_gen[0].E() >= 1.) hPassCuts_MC->Fill(3); //fi (E_gamma cut)
      }//fi (photons inclusive)
      // Proton
      if(scatp4_gen.size() == 1 && !kGenDVCS){
        hPassCuts_MC->Fill(4);

        // ...and proton track theta cut (FF region)
	if(scatp4_gen[0].Theta() <= 0.02) hPassCuts_MC->Fill(5); //fi (p'theta cut)
      } //fi (protons inclusive)
      // ------------------Event---------------------
      // Full DVCS multiplicity...
      if(scate4_gen.size() == 1 && scatg4_gen.size() == 1 && scatp4_gen.size() == 1 && !kGenDVCS){
        hPassCuts_MC->Fill(6);
        // ...and no extra -ve tracks
	if(!kNonElecNeg_MC){
          hPassCuts_MC->Fill(7);

          // ...and no ZDC/OMD particles
	  if(!kOtherFFTrack_MC){
            hPassCuts_MC->Fill(8);

            // ...and (E-pz) cut
	    float EmPz = (scate4_gen[0]+scatp4_gen[0]+scatg4_gen[0]).E() - (scate4_gen[0]+scatp4_gen[0]+scatg4_gen[0]).Pz();
            if(EmPz >= 15. && EmPz <= 25.){
              hPassCuts_MC->Fill(9);

              // ...and missing pT cut
	      float pTmiss = calcPtMiss_3Body(beame4, beamp4, scate4_gen[0], scatp4_gen[0], scatg4_gen[0]);
              if(pTmiss <= 0.5) hPassCuts_MC->Fill(10); //fi (pTmiss cut)
	    } //fi (E-pz cut)
	  }   //fi (ZDC/OMD veto)
	}     //fi (non-electron -ve particles)
      }       //fi (e'p'gamma)
      // e'gamma final state (fakes for eXBE reco.
      if(scate4_gen.size() == 1 && scatg4_gen.size() == 1 && !kGenDVCS){
        hPassCuts_MC->Fill(11);

        // ...and non-elec -ve veto
	if(!kNonElecNeg_MC){
          hPassCuts_MC->Fill(12);
          // ...and ZDC/OMD veto
	  if(!kOtherFFTrack_MC){
            hPassCuts_MC->Fill(13);

            // ...and barrel +ve track veto
	    if(!kBarrelPos_MC){
              hPassCuts_MC->Fill(14);

              // ...and (E-pz cut)
	      float EmPz = (scate4_gen[0]+scatg4_gen[0]).E() - (scate4_gen[0]+scatg4_gen[0]).Pz();
              if(EmPz >= 15. && EmPz <= 25.){
                hPassCuts_MC->Fill(15);

                // ...and missing pT cut
		float pTmiss = calcPtMiss_2Body(beame4, beamp4, scate4_gen[0], scatg4_gen[0]);
                if(pTmiss <= 0.5) hPassCuts_MC->Fill(16); //fi (pTmiss cut)
	      } //fi (E-pz cut)
	    }   //fi (barrel +ve tracks)
	  }     //fi (ZDC/OMD veto)
	}       //fi (non-elec -ve veto)
      }         //fi (e'gamma)
      
      //cout<<"[DEBUG]: MC PASSCUTS FILLED"<<endl;

      // Count no. of events which pass cuts - reconstructed only
      // ---------------------Single particle---------------------
      // Electron
      if(scate4_rec.size() == 1 && !kGenDVCS){
        hPassCuts_Rec->Fill(0);

        // ...and Q2 > 1 GeV2
	fQ2 = calcQ2_Elec(beame4, scate4_rec[0]);
        if(fQ2 >= 1.) hPassCuts_Rec->Fill(1); //fi (Q2 cut)
      } //fi (electrons inclusive)
      // Photon
      if(scatg4_rec.size() == 1 && !kGenDVCS){
        hPassCuts_Rec->Fill(2);

        // ...and E_{gamma} > 1 GeV
	if(scatg4_rec[0].E() >= 1.) hPassCuts_Rec->Fill(3); //fi (E_gamma cut)
      }//fi (photons inclusive)
      // Proton - B0
      if(scatp4_rec.size() == 1 && scatp4_rom.size() == 0 && !kGenDVCS){
        hPassCuts_Rec->Fill(4);

        // ...and proton track theta cut (B0 region)
	if(scatp4_rec[0].Theta() >= 0.0055 && scatp4_rec[0].Theta() <= 0.02) hPassCuts_Rec->Fill(5); //fi (p'theta cut)
      } //fi (B0 proton)
      // Proton - RP
      if(scatp4_rom.size() == 1 && scatp4_rec.size() == 0 && !kGenDVCS){
        hPassCuts_Rec->Fill(4);

        // ...and proton track theta cut (RP region)
	if(scatp4_rom[0].Theta() > 0 && scatp4_rom[0].Theta() <= 0.005) hPassCuts_Rec->Fill(5); //fi (p'theta cut)
      } //fi (RP proton)
      // ------------------Event---------------------
      // Full DVCS multiplicity (using B0 protons)...
      if(scate4_rec.size() == 1 && scatg4_rec.size() == 1 && scatp4_rec.size() == 1 && scatp4_rom.size() == 0 && !kGenDVCS){
        hPassCuts_Rec->Fill(6);
        // ...and no extra -ve tracks
	if(!kNonElecNeg_Rec){
          hPassCuts_Rec->Fill(7);

          // ...and no ZDC/OMD particles
	  if(!kOtherFFTrack_Rec){
            hPassCuts_Rec->Fill(8);

            // ...and (E-pz) cut
	    float EmPz = (scate4_rec[0]+scatp4_rec[0]+scatg4_rec[0]).E() - (scate4_rec[0]+scatp4_rec[0]+scatg4_rec[0]).Pz();
            if(EmPz >= 15. && EmPz <= 25.){
              hPassCuts_Rec->Fill(9);

              // ...and missing pT cut
	      float pTmiss = calcPtMiss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rec[0], scatg4_rec[0]);
              if(pTmiss <= 0.5) hPassCuts_Rec->Fill(10); //fi (pTmiss cut)
	    } //fi (E-pz cut)
	  }   //fi (ZDC/OMD veto)
	}     //fi (non-electron -ve particles)
      }       //fi (e'p'gamma)
      // Full DVCS multiplicity (using RP tracks)...
      if(scate4_rec.size() == 1 && scatg4_rec.size() == 1 && scatp4_rom.size() == 1 && scatp4_rec.size() == 0 && !kGenDVCS){
        hPassCuts_Rec->Fill(6);
        // ...and no extra -ve tracks
	if(!kNonElecNeg_Rec){
          hPassCuts_Rec->Fill(7);

          // ...and no ZDC/OMD particles
	  if(!kOtherFFTrack_Rec){
            hPassCuts_Rec->Fill(8);

            // ...and (E-pz) cut
	    float EmPz = (scate4_rec[0]+scatp4_rom[0]+scatg4_rec[0]).E() - (scate4_rec[0]+scatp4_rom[0]+scatg4_rec[0]).Pz();
            if(EmPz >= 15. && EmPz <= 25.){
              hPassCuts_Rec->Fill(9);

              // ...and missing pT cut
	      float pTmiss = calcPtMiss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rom[0], scatg4_rec[0]);
              if(pTmiss <= 0.5) hPassCuts_Rec->Fill(10); //fi (pTmiss cut)
	    } //fi (E-pz cut)
	  }   //fi (ZDC/OMD veto)
	}     //fi (non-electron -ve particles)
      }       //fi (e'p'gamma)
      // e'gamma final state (fakes for eXBE reco.)
      if(scate4_rec.size() == 1 && scatg4_rec.size() == 1 && !kGenDVCS){
        hPassCuts_Rec->Fill(11);

        // ...and non-elec -ve veto
	if(!kNonElecNeg_Rec){
          hPassCuts_Rec->Fill(12);
          // ...and ZDC/OMD veto
	  if(!kOtherFFTrack_Rec){
            hPassCuts_Rec->Fill(13);

            // ...and barrel +ve track veto
	    if(!kBarrelPos_Rec){
              hPassCuts_Rec->Fill(14);

              // ...and (E-pz cut)
	      float EmPz = (scate4_rec[0]+scatg4_rec[0]).E() - (scate4_rec[0]+scatg4_rec[0]).Pz();
              if(EmPz >= 15. && EmPz <= 25.){
                hPassCuts_Rec->Fill(15);

                // ...and missing pT cut
		float pTmiss = calcPtMiss_2Body(beame4, beamp4, scate4_rec[0], scatg4_rec[0]);
                if(pTmiss <= 0.5) hPassCuts_Rec->Fill(16); //fi (pTmiss cut)
	      } //fi (E-pz cut)
	    }   //fi (barrel +ve tracks)
	  }     //fi (ZDC/OMD veto)
	}       //fi (non-elec -ve veto)
      }         //fi (e'gamma)

      //cout<<"[DEBUG]: RECO PASSCUTS FILLED"<<endl;

      // Mandelstam t distributions
      // MC truth
      if(kGenDVCS) h_t_Truth->Fill(calcT_BABE(beamp4,scatp4_gen[0])); 
      // Reconstructed and MC accepted - B0 only
      if(!kGenDVCS                                                                  // Not generated DVCS
	 && applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rec, scatg4_rec, "B0") // Single e', single photon, single p', proton w/in B0
	 && scatp4_rom.size() == 0                                                  // No RP tracks (assume all RP are protons)
	 && !kNonElecNeg_Rec                                                        // Non-electron -ve track veto
	 && !kOtherFFTrack_Rec){                                                    // OMD/ZDC veto
	
	// Calculations
	Float_t t_acc = calcT_BABE(beamp4,scatp4_aso[0]);
	Float_t t_rec = calcT_BABE(beamp4,scatp4_rec[0]);
	// Distributions
       	h_t_B0Acc->Fill(t_acc);
	h_t_B0Reco->Fill(t_rec);
	//cout<<"[DEBUG]: RECO T (B0) FILLED"<<endl;
      }
      // Reconstructed and accepted - RP only
      if(!kGenDVCS                                                                  // Not generated DVCS
	 && applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rom, scatg4_rec, "RP") // Single e', single photon, single p', proton w/in RP (assume RP = p')
	 && scatp4_rec.size() == 0                                                  // No B0 tracks
	 && !kNonElecNeg_Rec                                                        // Non-electron -ve track veto
	 && !kOtherFFTrack_Rec){                                                    // OMD/ZDC veto
	
	//cout<<"[DEBUG]: TRYING FULL EPG (RP)"<<endl;
	//cout<<"[DEBUG]: N_REC = "<<scatp4_rom.size()<<"\tN_GEN = "<<scatp4_gen.size()<<endl;
	
	if(scatp4_gen.size() > 0){
	  Float_t t_acc = calcT_BABE(beamp4,scatp4_gen[0]);
	  h_t_RPAcc->Fill(t_acc);
	}
	Float_t t_rec = calcT_BABE(beamp4,scatp4_rom[0]);
	h_t_RPReco->Fill(t_rec);
	//cout<<"[DEBUG]: RECO T (RP) FILLED"<<endl;
      }
      // Semi-inclusive calculation - ignore if proton is detected or not
      // Also plot missing kinematics
      if(!kGenDVCS                                                                // Not generated DVCS
	 && applyCuts_Electron(beame4,scate4_rec) && applyCuts_Photon(scatg4_rec) // Single e', single photon, Q2 > 1 GeV
	 && !kBarrelPos_Rec                                                       // Barrel +ve track veto
	 && !kNonElecNeg_Rec                                                      // Non-electron -ve track veto
	 && !kOtherFFTrack_Rec){                                                  // OMD/ZDC veto
	
	//cout<<"[DEBUG]: TRYING SEMI-INCLUSIVE"<<endl;
	// Calculations
	Float_t t_acc = calcT_MethodL(beame4,beamp4,scate4_aso[0],fMass_proton,scatg4_aso[0]);
	Float_t t_rec = calcT_MethodL(beame4,beamp4,scate4_rec[0],fMass_proton,scatg4_rec[0]);
	// Distributions
	h_t_LCAcc->Fill(t_acc);
	h_t_LCReco->Fill(t_rec);

	h_MM2eg_Rec->Fill(calcM2Miss_2Body(beame4, beamp4, scate4_rec[0], scatg4_rec[0]));
	h_PMeg_Rec->Fill(calcPMiss_2Body(beame4, beamp4, scate4_rec[0], scatg4_rec[0]));
	h_PTMeg_Rec->Fill(calcPtMiss_2Body(beame4, beamp4, scate4_rec[0], scatg4_rec[0]));
	h_EMeg_Rec->Fill(calcEMiss_2Body(beame4, beamp4, scate4_rec[0], scatg4_rec[0]));

	//cout<<"[DEBUG]: RECO T (SEMI-INC) FILLED"<<endl;
      }

      // Plotting full e'p'gamma missing pT (w/out pT cut applied)
      if(!kGenDVCS
	 && applyCuts_Electron(beame4,scate4_rec) && applyCuts_Photon(scatg4_rec) && applyCuts_Proton(scatp4_rec, "B0") && scatp4_rom.size() == 0
	 && (scate4_rec[0]+scatg4_rec[0]).E() - (scate4_rec[0]+scatg4_rec[0]).Pz() >= fMin_EmPz && (scate4_rec[0]+scatg4_rec[0]).E() - (scate4_rec[0]+scatg4_rec[0]).Pz() <= fMax_EmPz
	 && !kBarrelPos_Rec                                                    
	 && !kNonElecNeg_Rec                                                  
	 && !kOtherFFTrack_Rec) h_PTMepg_Rec->Fill(calcPtMiss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rec[0], scatg4_rec[0]));
      if(!kGenDVCS
	 && applyCuts_Electron(beame4,scate4_rec) && applyCuts_Photon(scatg4_rec) && applyCuts_Proton(scatp4_rom, "RP") && scatp4_rec.size() == 0
	 && (scate4_rec[0]+scatg4_rec[0]).E() - (scate4_rec[0]+scatg4_rec[0]).Pz() >= fMin_EmPz && (scate4_rec[0]+scatg4_rec[0]).E() - (scate4_rec[0]+scatg4_rec[0]).Pz() <= fMax_EmPz
	 && !kBarrelPos_Rec                                                    
	 && !kNonElecNeg_Rec                                                  
	 && !kOtherFFTrack_Rec) h_PTMepg_Rec->Fill(calcPtMiss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rom[0], scatg4_rec[0]));
      
      
      //cout<<"[DEBUG]: MANDELSTAM T FILLED"<<endl;

    } // END OF EVENT LOOP - MAIN LOOP

    fileCounter++;
  } // END OF FILE LIST
  
  //------------------------------------------------------------
  // Write to output file
  //------------------------------------------------------------
  fOutFile->cd();
  
  // Eta - MC
  h_eta_MCp->Write();
  h_eta_MCe->Write();
  h_eta_MCg->Write();
  // Eta - reco.
  h_eta_RPp->Write();
  h_eta_RPPp->Write();
  h_eta_RPe->Write();
  h_eta_RPg->Write();
  // Other single-particle plots
  h_EoverP_elec->Write();
  h_2D_EvEta_e->Write();
  h_2D_EvEta_g->Write();
  h_2D_EvEta_p->Write();
  // t distributions
  h_t_Truth->Write();
  h_t_B0Acc->Write();
  h_t_RPAcc->Write();
  h_t_B0Reco->Write();
  h_t_RPReco->Write();
  h_t_LCAcc->Write();
  h_t_LCReco->Write();
  hCount->Write();
  hPassCuts_MC->Write();
  hPassCuts_Rec->Write();
  // 2-body exclusivity variables
  h_MM2eg_Rec->Write();
  h_PMeg_Rec->Write();
  h_PTMeg_Rec->Write();
  h_EMeg_Rec->Write();
  h_PTMepg_Rec->Write();

  return;
}
