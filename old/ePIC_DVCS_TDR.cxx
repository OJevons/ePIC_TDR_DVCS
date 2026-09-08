// ePIC DVCS analysis class definition

// ROOT Includes
#include <TSystem.h>
#include <TMath.h>
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <Math/LorentzVector.h>

// Class header include
#include "ePIC_DVCS_TASK.h"

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

   //----------------------------------
   // INSERT ANY OTHER PHOTON CUTS HERE
   //----------------------------------

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
  // Need to know beam proton energy for minimum theta in RP
  Float_t beamP{0};
  if(sEnergy == "5x41"){
    beamP=41.0;
  }
  else if(sEnergy == "10x100"){
    beamP=100.0;
  }
  else if(sEnergy == "10x130"){
    beamP=130.0;
  }
  else if(sEnergy == "18x275"){
    beamP=275.0;
  }
  else{
    beamP=100.0;
  }
  // RP momentum acceptance < 200 MeV
  
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

  // 1. MAXIMUM T CUT FOR ROMAN POTS
  if(sProtonDet != "B0" && sProtonDet != "RP" && sProtonDet != "all") sProtonDet="all";
  //if(sProtonDet == "RP" && ft > fMaxt_RP) passCuts = kFALSE;

  // 2. BJORKEN X CUT (removing tail from reconstructed histogram)
  //if(TMath::Log10(fxB) < fxB_Tail) passCuts = kFALSE;
  
  // 3. MAXIMUM MISSING MASS^2
  if(TMath::Abs(fM2miss) > fMax_M2miss) passCuts = kFALSE;

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

  // 3. Proton cuts
  passCuts = applyCuts_Proton(scatp, sProtonDet);
  // Exit from function if failure
  if(!passCuts) return passCuts;

  // 4. Event cuts
  // Need to calculate t, xB and MM2 first (e'p'y final state guaranteed by this point)
  fxB = calcX_Elec(beame, scate[0], beamp);
  fM2miss = calcM2Miss_3Body(beame, beamp, scate[0], scatp[0], scatg[0]);
  ft = calcT_BABE(beamp, scatp[0]);
  passCuts = applyCuts_DVCS(sProtonDet);

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
  // Holding vectors for beam - undoing crossing angle ONLY
  P3EVector p_beam(fXAngle*p.E(), 0., p.E(), p.E());
  P3EVector e_beam(0., 0., -k.E(), k.E());
  
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

// Calculate angle between hadronic and leptonic planes (Trento phi)
// Using planes defined by [k, q] and [q, p']
// Source: Bachetta, A. et al; Phys. Rev. D (2004); eq. 16
Double_t ePIC_DVCS_TASK::calcTrentoPhi_qp(P3EVector k, P3EVector kprime, P3EVector pprime){  
  MomVector k3 = k.Vect();
  MomVector kp3 = kprime.Vect();
  MomVector pp3 = pprime.Vect();
  MomVector qhat3 = (k3-kp3).Unit();

  // Define leptonic plane using beam and scattered electron
  MomVector lNorm = qhat3.Cross(k3);
  lNorm /= lNorm.R();
  // Define hadronic plane using q vector and scattered proton
  MomVector hNorm = qhat3.Cross(pp3);
  hNorm /= hNorm.R();

  // Angle() function just returns magnitude of angle
  // If p' vector has a component parallel to the leptonic normal, should be positive. If opposite, negative.
  return TMath::Sign(1.,pp3.Dot(lNorm))*Angle(lNorm,hNorm);
}

// Alternate Trento Phi calculation (using final state photon and proton)
// Angle between planes defined by [k, k'] and [g, p']
Double_t ePIC_DVCS_TASK::calcTrentoPhi_pg(P3EVector k, P3EVector kprime, P3EVector pprime, P3EVector gprime){  
  MomVector k3 = k.Vect();
  MomVector kp3 = kprime.Vect();
  MomVector pp3 = pprime.Vect();
  MomVector gp3 = gprime.Vect();

  // Define leptonic plane using beam and scattered electron
  MomVector lNorm = kp3.Cross(k3);
  lNorm /= lNorm.R();
  // Define hadronic plane using scattered proton and final state photon
  MomVector hNorm = gp3.Cross(pp3);
  hNorm /= hNorm.R();

  // Angle() function just returns magnitude of angle
  // If p' vector has a component parallel to the leptonic normal, should be positive. If opposite, negative.
  return TMath::Sign(1.,pp3.Dot(lNorm))*Angle(lNorm,hNorm);
}

// Alternate Trento Phi calculation (using final state photon)
// Angle between planes defined by [k, k'] and [q, g]
Double_t ePIC_DVCS_TASK::calcTrentoPhi_qg(P3EVector k, P3EVector kprime, P3EVector gprime){  
  MomVector k3 = k.Vect();
  MomVector kp3 = kprime.Vect();
  MomVector gp3 = gprime.Vect();
  MomVector q3 = k3-kp3;

  // Define leptonic plane using beam and scattered electron
  MomVector lNorm = kp3.Cross(k3);
  lNorm /= lNorm.R();
  // Define hadronic plane using q vector and final state photon
  MomVector hNorm = q3.Cross(gp3);
  hNorm /= hNorm.R();

  // Angle() function just returns magnitude of angle
  // If g' vector has a component parallel to the leptonic normal, should be positive. If opposite, negative.
  return TMath::Sign(1.,gp3.Dot(lNorm))*Angle(lNorm,hNorm);
}

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
  // IF TESTING, LOAD TEST FILE LIST
  if(sSett != "hiAcc" && sSett != "hiDiv") sInList="./inputFileList_ePIC_25.10.2_10x100_hiAcc.list";
  
  ifstream fileListStream;
  fileListStream.open(sInList);
  string fileName;
  TFile* inputRootFile;

  //---------------------------------------------------------
  // Setup: Declare histograms
  //---------------------------------------------------------
  // Used for TDR:
  // 1) Single particle pseudorapidity distributions
  //    --> Scattered electron, real photon and scattered proton (separating B0 from Roman Pot tracks)
  // 2) Photon single particle theta resolution
  //    --> As 1D distribution and as 2D as function of generated photon theta
  // 3) Mandelstam t distribution
  //    --> [MC/Reco/acc.]
  // 4) Mandelstam t resolution
  //    --> Need 2D plot as function of generated t so that plotting macro can do resolutions

  // 1a) Eta - MC particles
  TH1D* h_eta_MCp   = new TH1D("eta_MCp",";#eta_{p'}(MC)", 275, -11.0, 11.0);
  TH1D* h_eta_MCe   = new TH1D("eta_MCe",";#eta_{e'}(MC)", 275, -11.0, 11.0);
  TH1D* h_eta_MCg   = new TH1D("eta_MCg",";#eta_{#gamma}(MC)", 275, -11.0, 11.0);
  // 1b) Eta - reco. particles
  TH1D* h_eta_RPp   = new TH1D("eta_RPp",";#eta_{p'}(Reco)", 275, -11.0, 11.0);
  TH1D* h_eta_RPPp   = new TH1D("eta_RPPp",";#eta_{p'}(Reco)", 275, -11.0, 11.0);
  TH1D* h_eta_RPe   = new TH1D("eta_RPe",";#eta_{e'}(Reco)", 275, -11.0, 11.0);
  TH1D* h_eta_RPg   = new TH1D("eta_RPg",";#eta_{#gamma}(Reco)", 275, -11.0, 11.0);

  // 2) Photon theta resolution
  TH1D* h_PhotRes_theta = new TH1D("photres_theta",";#theta_{#gamma}(Reco)-#theta_{#gamma}(MC) [deg]",360,-90,90);
  TH2D* h_PhotRes2D_theta = new TH2D("photres2d_theta",";#theta_{#gamma, MC} [deg]; #Delta#theta_{#gamma} [deg]",370,0,185,360,-90,90);

  // 3) t distribution
  TH1D* h_t_Truth  = new TH1D("t_truth" ,";|t|(MC) [(GeV/#it{c}^{2})^{2}]"           , 20, 0., 2.);
  TH1D* h_t_B0Acc  = new TH1D("t_b0acc" ,";|t|(MC|Reco. - B0) [(GeV/#it{c}^{2})^{2}]", 20, 0., 2.);
  TH1D* h_t_RPAcc  = new TH1D("t_rpacc" ,";|t|(MC|Reco. - RP) [(GeV/#it{c}^{2})^{2}]", 20, 0., 2.);
  TH1D* h_t_B0Reco = new TH1D("t_b0reco",";|t|(Reco. - B0) [(GeV/#it{c}^{2})^{2}]"   , 20, 0., 2.);
  TH1D* h_t_RPReco = new TH1D("t_rpreco",";|t|(Reco. - RP) [(GeV/#it{c}^{2})^{2}]"   , 20, 0., 2.);

  // 4) t resolution - as absolute or as percentage (plot as preferred)
  TH2D* h_tResB0_2d = new TH2D("tresb0_2d",";|t|_{MC} [(GeV/#it{c})^{2}];#delta t [(GeV/#it{c})^{2}]", 20, 0., 2., 100, -5., 5.);
  TH2D* h_tResRP_2d = new TH2D("tresrp_2d",";|t|_{MC} [(GeV/#it{c})^{2}];#delta t [(GeV/#it{c})^{2}]", 20, 0., 2., 100, -5., 5.);
  TH2D* h_tResB0Pct_2d = new TH2D("tresb0pct_2d",";|t|_{MC} [(GeV/#it{c})^{2}];#Deltat/t_{MC}", 20, 0., 2., 200, -1., 1.);
  TH2D* h_tResRPPct_2d = new TH2D("tresrppct_2d",";|t|_{MC} [(GeV/#it{c})^{2}];#Deltat/t_{MC}", 20, 0., 2., 200, -1., 1.);


  //---------------------------------------------------------
  // Loop over files in list
  //---------------------------------------------------------
  int fileCounter{0};

  // 4-vectors for beam particles - need these defined outside of file loop
  P3EVector beame4(0,0,0,-1);     // Beam electron (generated)
  P3EVector beamp4(0,0,0,-1);     // Beam proton (generated)
  
  // Start file loop
  while(getline(fileListStream,fileName)){
    //---------------------------------------------------------
    // Open 1 file from list at a time, then get event TTree
    //---------------------------------------------------------
    // Get file
    TString tmp{fileName};
    std::cout<<"Input file "<<fileCounter<<" : "<<tmp<<std::endl;
    auto inputRootFile = TFile::Open(tmp);
    if(!inputRootFile){ std::cout<<"MISSING_ROOT_FILE"<<tmp<<endl; continue;}
    fileCounter++;

    // Get TTree
    TTree * evtTree = (TTree*)inputRootFile->Get("events");
    if (!(inputRootFile->GetListOfKeys()->Contains("events"))) continue;
    int numEvents = evtTree->GetEntries();
    std::cout<<"File has "<<numEvents<<" events..."<<std::endl;

    //---------------------------------------------------------
    // Declare TTreeReader, and choose appropriate branches
    //---------------------------------------------------------
    TTreeReader tree_reader(evtTree);
    // MC particles
    TTreeReaderArray<double> mc_px_array        = {tree_reader, "MCParticles.momentum.x"};
    TTreeReaderArray<double> mc_py_array        = {tree_reader, "MCParticles.momentum.y"};
    TTreeReaderArray<double> mc_pz_array        = {tree_reader, "MCParticles.momentum.z"};
    TTreeReaderArray<double> mc_mass_array      = {tree_reader, "MCParticles.mass"};
    TTreeReaderArray<int>    mc_genStatus_array = {tree_reader, "MCParticles.generatorStatus"};
    TTreeReaderArray<int>    mc_pdg_array       = {tree_reader, "MCParticles.PDG"};
    // Reconstructed/MC particle associations - BARREL (using ReconstructedParticles branch)
    TTreeReaderArray<unsigned int> assoc_rec_id = {tree_reader, "ReconstructedParticleAssociations.recID"};
    TTreeReaderArray<unsigned int> assoc_sim_id = {tree_reader, "ReconstructedParticleAssociations.simID"};
    // Reconstructed particles - BARREL (using ReconstructedParticles branch)
    TTreeReaderArray<float>  re_px_array     = {tree_reader, "ReconstructedParticles.momentum.x"};
    TTreeReaderArray<float>  re_py_array     = {tree_reader, "ReconstructedParticles.momentum.y"};
    TTreeReaderArray<float>  re_pz_array     = {tree_reader, "ReconstructedParticles.momentum.z"};
    TTreeReaderArray<float>  re_e_array      = {tree_reader, "ReconstructedParticles.energy"};
    TTreeReaderArray<float>  re_charge_array = {tree_reader, "ReconstructedParticles.charge"};
    TTreeReaderArray<float>  re_mass_array   = {tree_reader, "ReconstructedParticles.mass"};
    TTreeReaderArray<int>    re_pdg_array    = {tree_reader, "ReconstructedParticles.PDG"};
    // Reconstructed/MC particle associations - B0 (using ReconstructedTruthSeededChargedParticles branch)
    TTreeReaderArray<unsigned int> tsassoc_rec_id = {tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.recID"};
    TTreeReaderArray<unsigned int> tsassoc_sim_id = {tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.simID"};
    // Reconstructed particles - B0 (using ReconstructedTruthSeededChargedParticles branch)
    TTreeReaderArray<float>  tsre_px_array     = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x"};
    TTreeReaderArray<float>  tsre_py_array     = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y"};
    TTreeReaderArray<float>  tsre_pz_array     = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z"};
    TTreeReaderArray<float>  tsre_e_array      = {tree_reader, "ReconstructedTruthSeededChargedParticles.energy"};
    TTreeReaderArray<float>  tsre_charge_array = {tree_reader, "ReconstructedTruthSeededChargedParticles.charge"};
    TTreeReaderArray<float>  tsre_mass_array   = {tree_reader, "ReconstructedTruthSeededChargedParticles.mass"};
    // RP hits
    TTreeReaderArray<float> global_hit_RP_x = {tree_reader, "ForwardRomanPotRecParticles.referencePoint.x"};
    TTreeReaderArray<float> global_hit_RP_y = {tree_reader, "ForwardRomanPotRecParticles.referencePoint.y"};
    TTreeReaderArray<float> global_hit_RP_z = {tree_reader, "ForwardRomanPotRecParticles.referencePoint.z"};
    TTreeReaderArray<float> rp_px_array     = {tree_reader, "ForwardRomanPotRecParticles.momentum.x"};
    TTreeReaderArray<float> rp_py_array     = {tree_reader, "ForwardRomanPotRecParticles.momentum.y"};
    TTreeReaderArray<float> rp_pz_array     = {tree_reader, "ForwardRomanPotRecParticles.momentum.z"};
    TTreeReaderArray<float> rp_mass_array   = {tree_reader, "ForwardRomanPotRecParticles.mass"};
    TTreeReaderArray<int>   rp_pdg_array    = {tree_reader, "ForwardRomanPotRecParticles.PDG"};

    tree_reader.SetEntriesRange(0, evtTree->GetEntries());
    
    // If averaging beams from file
    if(!kUseEventBeams){
      // MUST DO THIS FIRST
      // Full run over tree in first file before anything else
      // Calculate beams from average of individual event beam particles
      if(fileCounter == 1){
	// Accumulator variables
	P3EVector beame4_acc(0,0,0,-1);
	P3EVector beamp4_acc(0,0,0,-1);
	
	while (tree_reader.Next()){
	  // Beams for each event
	  P3EVector beame4_evt(0,0,0,-1);
	  P3EVector beamp4_evt(0,0,0,-1);
	  TVector3 mctrk;
	  for(int imc=0;imc<mc_px_array.GetSize();imc++){
	    mctrk.SetXYZ(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);

	    // Beam particles ==> Generator Status 4
	    if(mc_genStatus_array[imc] == 4){
	      // Proton
	      if(mc_pdg_array[imc] == 2212){ 
		beamp4_evt.SetCoordinates(mctrk.X(), mctrk.Y(), mctrk.Z(), calcE(mctrk,mc_mass_array[imc]));
	      }
	      // Electron
	      else if(mc_pdg_array[imc] == 11){ 
		beame4_evt.SetCoordinates(mctrk.X(), mctrk.Y(), mctrk.Z(), calcE(mctrk,mc_mass_array[imc]));
	      }		
	    } // Found beam particles for event
	  } // End MCParticles loop
	  
	  // Add found beams to accumulator
	  beame4_acc += beame4_evt;
	  beamp4_acc += beamp4_evt;
	} // DONE while loop
	
	// Divide by number of events in file
	beame4.SetCoordinates(beame4_acc.X()/evtTree->GetEntries(), beame4_acc.Y()/evtTree->GetEntries(), beame4_acc.Z()/evtTree->GetEntries(), beame4_acc.E()/evtTree->GetEntries());
	beamp4.SetCoordinates(beamp4_acc.X()/evtTree->GetEntries(), beamp4_acc.Y()/evtTree->GetEntries(), beamp4_acc.Z()/evtTree->GetEntries(), beamp4_acc.E()/evtTree->GetEntries());
      
	
	// Undo afterburn on beam particles and calculate "postburn" variables
	undoAfterburnAndCalc(beamp4,beame4);
	
	// Restart TTreeReader for first file
	tree_reader.Restart();

	std::cout<<"File 1 - beams\n\te:"<<beame4<<"\n\tp:"<<beamp4<<std::endl;
      } // fi (fileCounter == 1)
      else std::cout<<"Using beams from first file."<<std::endl;
      
    } // fi (kUseEventBeams)

    //---------------------------------------------------------
    // Loop over all entries in event TTree
    //---------------------------------------------------------
    while (tree_reader.Next()){ 
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

      // Holding 3-vectors
      TVector3 mctrk, assoctrk, recotrk;
      
      // Get beams for each event
      if(kUseEventBeams){
	for(int imc=0;imc<mc_px_array.GetSize();imc++){
	  mctrk.SetXYZ(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);
	  
	  // Beam particles ==> Generator Status 4
	  if(mc_genStatus_array[imc] == 4){
	    // Proton
	    if(mc_pdg_array[imc] == 2212){ 
	      beamp4.SetCoordinates(mctrk.X(), mctrk.Y(), mctrk.Z(), calcE(mctrk,mc_mass_array[imc]));
	    }
	    // Electron
	    else if(mc_pdg_array[imc] == 11){ 
	      beame4.SetCoordinates(mctrk.X(), mctrk.Y(), mctrk.Z(), calcE(mctrk,mc_mass_array[imc]));
	    }		
	  } // Found beam particles
	}
	// Undo afterburn on beam particles and calculate "postburn" variables
	undoAfterburnAndCalc(beamp4,beame4);
      }

      //---------------------------------------------------------
      // Fill particle holding arrays
      //---------------------------------------------------------
      // 1. MC generated
      for(int imc=0;imc<mc_px_array.GetSize();imc++){
	mctrk.SetXYZ(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);
	P3EVector q_scat(mctrk.X(),mctrk.Y(),mctrk.Z(),calcE(mctrk,mc_mass_array[imc]));
	
	// Undo afterburner
	undoAfterburn(q_scat);

	// Look for scattered particles ==> Generator status 1
	if(mc_genStatus_array[imc] == 1){
	  if(mc_pdg_array[imc] == 2212){
	    scatp4_gen.push_back(q_scat);
	  }
	  if(mc_pdg_array[imc] == 11){
	    scate4_gen.push_back(q_scat);
	  }
	  if(mc_pdg_array[imc] == 22){
	    scatg4_gen.push_back(q_scat);
	  }
	} // Found scattered particles
      }// End of generated particles loop
      
      // 2. and 3. Associated MC tracks and their matched reconstucted tracks
      // ONLY IF USING EXPLICIT TRACK MATCHING
      if(kUseExplicitMatch){
	unsigned int mc_assoc_index = -1;
	// LOOK FOR ELECTRONS AND PHOTONS (using ReconstructedParticleAssociations)
	for(unsigned int iAssoc{0};iAssoc<assoc_rec_id.GetSize();iAssoc++){
	  mc_assoc_index = assoc_sim_id[iAssoc];

	  // If reco track isn't associated to an MC track, then skip
	  if(mc_assoc_index == -1) continue;
	
	  assoctrk.SetXYZ(mc_px_array[mc_assoc_index], mc_py_array[mc_assoc_index], mc_pz_array[mc_assoc_index]); 
	  P3EVector q_assoc(assoctrk.X(),assoctrk.Y(),assoctrk.Z(),calcE(assoctrk,mc_mass_array[mc_assoc_index]));
	  undoAfterburn(q_assoc);
	  recotrk.SetXYZ(re_px_array[iAssoc], re_py_array[iAssoc], re_pz_array[iAssoc]);
	  P3EVector q_reco(recotrk.X(),recotrk.Y(),recotrk.Z(),calcE(recotrk,mc_mass_array[mc_assoc_index]));
	  undoAfterburn(q_reco);

	  // Fill track vectors based on associated PID
	  // Electrons
	  if(mc_genStatus_array[mc_assoc_index] == 1 && mc_pdg_array[mc_assoc_index] == 11){ 
	    scate4_aso.push_back(q_assoc); 
	    scate4_rec.push_back(q_reco); 
	  }
	  // Photons
	  if(mc_genStatus_array[mc_assoc_index] == 1 && mc_pdg_array[mc_assoc_index] == 22){ 
	    scatg4_aso.push_back(q_assoc); 
	    scatg4_rec.push_back(q_reco); 
	  }
	} // End of associations loop

	mc_assoc_index=-1; // Reset association index
	// THEN LOOK FOR PROTONS (using ReconstructedTruthSeededChargedParticleAssociations)
	for(unsigned int iTSAssoc{0};iTSAssoc<tsassoc_rec_id.GetSize();iTSAssoc++){
	  mc_assoc_index = tsassoc_sim_id[iTSAssoc];

	  // Only care about protons here (PID 2212)
	  if(mc_assoc_index != -1 && mc_genStatus_array[mc_assoc_index] == 1 && mc_pdg_array[mc_assoc_index] == 2212){
	    assoctrk.SetXYZ(mc_px_array[mc_assoc_index], mc_py_array[mc_assoc_index], mc_pz_array[mc_assoc_index]);
	    P3EVector q_assoc(assoctrk.X(),assoctrk.Y(),assoctrk.Z(),calcE(assoctrk,mc_mass_array[mc_assoc_index]));
	    undoAfterburn(q_assoc);
	    recotrk.SetXYZ(tsre_px_array[iTSAssoc], tsre_py_array[iTSAssoc], tsre_pz_array[iTSAssoc]);
	    P3EVector q_reco(recotrk.X(),recotrk.Y(),recotrk.Z(),calcE(recotrk,mc_mass_array[mc_assoc_index]));
	    undoAfterburn(q_reco);

	    scatp4_aso.push_back(q_assoc); 
	    scatp4_rec.push_back(q_reco); 
	  }
	} // End of truth-seeded association loop
      }

      // 2. and 3. Associated MC and reconstructed tracks
      // IF NOT USING EXPLICIT MATCHING (default behaviour)
      else if(!kUseExplicitMatch){
	// 2. Associated MC (scattered)
	for(unsigned int iAssoc{0};iAssoc<assoc_rec_id.GetSize();iAssoc++){
	  unsigned int mc_assoc_index = assoc_sim_id[iAssoc];
	  assoctrk.SetXYZ(mc_px_array[mc_assoc_index], mc_py_array[mc_assoc_index], mc_pz_array[mc_assoc_index]);
	  P3EVector q_assoc(assoctrk.X(),assoctrk.Y(),assoctrk.Z(),calcE(assoctrk,mc_mass_array[mc_assoc_index]));
	  // Undo afterburner
	  undoAfterburn(q_assoc);
	  // Look for scattered particles ==> Generator status 1
	  if(mc_genStatus_array[mc_assoc_index] == 1){
	    if(mc_pdg_array[mc_assoc_index] == 11){ scate4_aso.push_back(q_assoc); }
	    if(mc_pdg_array[mc_assoc_index] == 22){ scatg4_aso.push_back(q_assoc); }
	  } // Found associated particles
	}// End of associated particles loop
	
	// 2a. Associated MC protons (found using different association branch)
	for(unsigned int iTSAssoc{0};iTSAssoc<tsassoc_rec_id.GetSize();iTSAssoc++){
	  unsigned int mc_assoc_index = tsassoc_sim_id[iTSAssoc];
	  assoctrk.SetXYZ(mc_px_array[mc_assoc_index], mc_py_array[mc_assoc_index], mc_pz_array[mc_assoc_index]); 
	  P3EVector q_assoc(assoctrk.X(),assoctrk.Y(),assoctrk.Z(),calcE(assoctrk,mc_mass_array[mc_assoc_index]));
	  undoAfterburn(q_assoc);
	  
	  if(mc_genStatus_array[mc_assoc_index] == 1 && mc_pdg_array[mc_assoc_index] == 2212){ scatp4_aso.push_back(q_assoc); }
	  
	} // End of truth-seeded association loop
	
	// 3. Reconstructed particles
	// Start with ACTS reconstructed particles (barrel and B0)
	if(re_px_array.GetSize() == 0) continue;
	else{
	  for(int ireco{0}; ireco<re_px_array.GetSize(); ireco++){
	    recotrk.SetXYZ(re_px_array[ireco], re_py_array[ireco], re_pz_array[ireco]);
	    // If not using ePIC PID, assume particles based on charge of track
	    // Look for electrons and photons from ReconstructedParticles branch
	    if(!kUsePID){
	      // Negative => ELECTRON
	      if(re_charge_array[ireco] == -1){
		P3EVector q_reco(recotrk.X(),recotrk.Y(),recotrk.Z(),calcE(recotrk,fMass_electron));
		undoAfterburn(q_reco);
		scate4_rec.push_back(q_reco);
		//std::cout<<"[DEBUG] Reco. e' (-ve) FOUND"<<std::endl;
	      }
	      // Neutral => REAL PHOTON
	      if(re_charge_array[ireco] == 0){
		P3EVector q_reco(recotrk.X(),recotrk.Y(),recotrk.Z(),recotrk.Mag());
		undoAfterburn(q_reco);
		scatg4_rec.push_back(q_reco);
		//std::cout<<"[DEBUG] Reco. g (neutral) FOUND"<<std::endl;
	      }
	    }
	    // Using ePIC PID
	    else{
	      P3EVector q_reco(recotrk.X(),recotrk.Y(),recotrk.Z(),calcE(recotrk,re_mass_array[ireco]));
	      // Undo afterburner
	      undoAfterburn(q_reco);
	      if(re_pdg_array[ireco] == 11){ scate4_rec.push_back(q_reco); }
	      if(re_pdg_array[ireco] == 22){ scatg4_rec.push_back(q_reco); }
	    }
	  }// End of ACTS reconstructed particles loop
	}
	// 3a. Reconstructed B0 protons (ReconstructedTruthSeededChargedParticles branch)
	for(int ireco{0}; ireco<tsre_px_array.GetSize(); ireco++){
	  recotrk.SetXYZ(tsre_px_array[ireco], tsre_py_array[ireco], tsre_pz_array[ireco]);
	  
	  // Using track charge
	  if(!kUsePID){
	    // Positive => PROTON
	    if(re_charge_array[ireco] == 1){
	      P3EVector q_reco(recotrk.X(),recotrk.Y(),recotrk.Z(),calcE(recotrk,fMass_proton));
	      undoAfterburn(q_reco);
	      scatp4_rec.push_back(q_reco);
	      //std::cout<<"[DEBUG] Reco. B0 p' (+ve) FOUND"<<std::endl;
	    }
	  }
	  
	  // Not using track charge
	  else{
	    P3EVector q_reco(recotrk.X(),recotrk.Y(),recotrk.Z(),calcE(recotrk,tsre_mass_array[ireco]));
	    undoAfterburn(q_reco);
	    // Select on protons with eta cut
	    if(q_reco.Eta() > 4.2){ scatp4_rec.push_back(q_reco); }
	  }
	}// End of truth seeded charged particles
      }

      // Add in RP hits - only looking at protons
      // NO NEED TO UNDO AFTERBURNER FOR FF DETECTORS - NOT APPLIED IN FIRST PLACE
      for(int irpreco{0}; irpreco<rp_px_array.GetSize(); irpreco++){
 	recotrk.SetXYZ(rp_px_array[irpreco], rp_py_array[irpreco], rp_pz_array[irpreco]);	
	P3EVector q_rpreco(recotrk.X(),recotrk.Y(),recotrk.Z(),calcE(recotrk,rp_mass_array[irpreco]));
	if(rp_pdg_array[irpreco] == 2212){
	  scatp4_rom.push_back(q_rpreco);
	}
      }// End of RP reconstructed particles loop
      // NEED TO COMBINE B0 and RP tracks to avoid double counting
      vector<P3EVector> scatp4_all;   
      for(auto i:scatp4_rom) scatp4_all.push_back(i);
      for(auto j:scatp4_rec) scatp4_all.push_back(j);
      	

      //---------------------------------------------------------
      // Fill histograms
      //---------------------------------------------------------
      // Eta - generated particles
      // Need Q2 for electron cuts
      if(scate4_gen.size() == 0) fQ2 = 0;
      else fQ2 = calcQ2_Elec(beame4, scate4_gen[0]);
      
      if(applyCuts_Electron(beame4,scate4_gen)) h_eta_MCe->Fill(scate4_gen[0].Eta());
      if(applyCuts_Photon(scatg4_gen)) 	        h_eta_MCg->Fill(scatg4_gen[0].Eta());
      if(applyCuts_Proton(scatp4_gen, "all"))   h_eta_MCp->Fill(scatp4_gen[0].Eta());
      
      // Eta - reco. particles
      // Need Q2 for electron cuts
      if(scate4_rec.size() == 0) fQ2 = 0;
      else fQ2 = calcQ2_Elec(beame4, scate4_rec[0]);

      if(applyCuts_Electron(beame4,scate4_rec))	h_eta_RPe->Fill(scate4_rec[0].Eta());
      if(applyCuts_Photon(scatg4_rec))    	h_eta_RPg->Fill(scatg4_rec[0].Eta());
      if(applyCuts_Proton(scatp4_rec, "B0"))	h_eta_RPp->Fill(scatp4_rec[0].Eta());
      if(applyCuts_Proton(scatp4_rom, "RP"))    h_eta_RPPp->Fill(scatp4_rom[0].Eta());
      
      // Photon theta resolution
      Float_t th_rec{0}, th_gen{0};
      if(applyCuts_Photon(scatg4_rec) && applyCuts_Photon(scatg4_aso)){
	th_gen = scatg4_aso[0].Theta()*TMath::RadToDeg();
	th_rec = scatg4_rec[0].Theta()*TMath::RadToDeg();
	
	h_PhotRes_theta->Fill(th_rec-th_gen);
	h_PhotRes2D_theta->Fill(th_gen, th_rec-th_gen);
      }

      // Mandelstam t distributions
      // MC truth
      if(applyCuts_All(beame4, beamp4, scate4_gen, scatp4_gen, scatg4_gen, "all")) h_t_Truth->Fill(calcT_BABE(beamp4,scatp4_gen[0]));
      // Reconstructed and MC accepted - B0 only
      if(applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rec, scatg4_rec, "B0") && scatp4_rom.size() == 0){
	// Calculations
	Float_t t_acc = calcT_BABE(beamp4,scatp4_aso[0]);
	Float_t t_rec = calcT_BABE(beamp4,scatp4_rec[0]);
	// Distributions
       	h_t_B0Acc->Fill(t_acc);
	h_t_B0Reco->Fill(t_rec);
      }
      // Reconstructed and accepted - RP only
      if(applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rom, scatg4_rec, "RP") && scatp4_rec.size() == 0){
	// Calculations
	Float_t t_acc = calcT_BABE(beamp4,scatp4_gen[0]);
	Float_t t_rec = calcT_BABE(beamp4,scatp4_rom[0]);
	// Distributions
	h_t_RPAcc->Fill(t_acc);
	h_t_RPReco->Fill(t_rec);
      }
    

      // Mandelstam t-resolution
      // ASSUME THAT GENERATED POSITIVE TRACK MATCHES RECONSTRUCTED POSITIVE TRACK
      Float_t t_rec{0}, t_gen{0};
      if(applyCuts_All(beame4, beamp4, scate4_gen, scatp4_gen, scatg4_gen, "all") && applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rec, scatg4_rec, "B0") && scatp4_rom.size()==0){
	t_gen = calcT_BABE(beamp4, scatp4_gen[0]);
	t_rec = calcT_BABE(beamp4, scatp4_rec[0]);
	h_tResB0_2d->Fill(t_gen, TMath::Abs(t_rec-t_gen));
	h_tResB0Pct_2d->Fill(t_gen, TMath::Abs(t_rec-t_gen)/t_gen);
      }
      if(applyCuts_All(beame4, beamp4, scate4_gen, scatp4_gen, scatg4_gen, "all") && applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rom, scatg4_rec, "RP") && scatp4_rec.size()==0){
	t_gen = calcT_BABE(beamp4, scatp4_gen[0]);
	t_rec = calcT_BABE(beamp4, scatp4_rom[0]);
	h_tResRP_2d->Fill(t_gen, TMath::Abs(t_rec-t_gen));
	h_tResRPPct_2d->Fill(t_gen, TMath::Abs(t_rec-t_gen)/t_gen);
      }
      
    } // END EVENT/TTREEREADER LOOP

    inputRootFile->Close();
  } // END OF FILE LOOP


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
  // Photon theta resolution
  h_PhotRes_theta->Write();
  h_PhotRes2D_theta->Write();
  // t distributions
  h_t_Truth->Write();
  h_t_B0Acc->Write();
  h_t_RPAcc->Write();
  h_t_B0Reco->Write();
  h_t_RPReco->Write();
  // 2D t resolution
  h_tResB0_2d->Write();
  h_tResRP_2d->Write();
  h_tResB0Pct_2d->Write();
  h_tResRPPct_2d->Write();
  
  fOutFile->Close();
  return;
}
