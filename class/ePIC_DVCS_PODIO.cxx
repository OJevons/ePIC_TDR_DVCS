// ePIC DVCS analysis class definition
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
// Flexible (ragged) Q2/xB/|t| binning read from a text file
#include "../include/DVCSBinning.hh"

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
void ePIC_DVCS_TASK::setOutFileName(TString name){
  std::cout<<"Output ROOT files: "<<name<<std::endl;
  sOutFileName = name;
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
  else if(sEnergy == "9x130"){
    fPBeam_p=130.0;
    fPBeam_e=9.0;
  }
  else if(sEnergy == "10x130"){
    fPBeam_p=130.0;
    fPBeam_e=10.0;
  }
  else if(sEnergy == "10x250"){
    fPBeam_p=250.0;
    fPBeam_e=10.0;
  }
  else if(sEnergy == "9x275"){
    fPBeam_p=275.0;
    fPBeam_e=9.0;
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
   // 3. Force energy loss by the electron
   if(scate[0].E() > beame.E()) passCuts = kFALSE;

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
    
    // Last cuts - fully exclusive DVCS event kinematics
    fpTmiss = calcPtMiss_2Body(beame, beamp, scate[0], scatg[0]);
    passCuts = applyCuts_DVCS(sProtonDet);
  }
  else{
    // E-pz cut (general to fully exclusive and semi-inclusive reco.)
    fEmPz = (scate[0]+scatg[0]+scatp[0]).E() - (scate[0]+scatg[0]+scatp[0]).Pz();
    
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
  // Holding vectors for beam - undoing crossing angle ONLY
  //P3EVector p_beam(fXAngle*p.E(), 0., p.E(), p.E());
  //P3EVector e_beam(0., 0., -k.E(), k.E());
  
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

// Calculate angle between hadronic and leptonic planes (Trento phi)
// Using planes defined by [k, q] and [q, g']
// Source: Bachetta, A. et al; Phys. Rev. D (2004); eq. 16
Double_t ePIC_DVCS_TASK::calcTrentoPhi_qg(P3EVector k, P3EVector p, P3EVector kprime, P3EVector gprime){  
  // Before calculating angle, boost into target rest frame
  //MomVector vTgtRest = p.BoostToCM();

  // Before calculating angle, boost into gamma*-p rest frame
  // Calculate q in lab frame
  P3EVector q = (k-kprime);
  // Boost vector
  MomVector vTgtRest = (p+q).BoostToCM();

  k = boost(k,vTgtRest);
  kprime = boost(kprime,vTgtRest);
  gprime = boost(gprime,vTgtRest);

  MomVector k3 = k.Vect();
  MomVector kp3 = kprime.Vect();
  MomVector gp3 = gprime.Vect();
  MomVector qhat3 = (k3-kp3).Unit();

  // Define leptonic plane using virtual photon and scattered electron
  MomVector lNorm = qhat3.Cross(kp3);
  lNorm /= lNorm.R();
  // Define hadronic plane using q vector and scattered photon
  MomVector hNorm = qhat3.Cross(gp3);
  hNorm /= hNorm.R();

  // Angle() function just returns magnitude of angle
  // If photon vector has a component parallel to the leptonic normal, should be positive. If opposite, negative.
  float phi = TMath::Sign(1.,gp3.Dot(lNorm))*Angle(lNorm,hNorm);

  if (phi < 0) return phi+2*TMath::Pi();
  else return phi;
}

// DIFF CALC. G.HILL PHD THESIS 2008
// EQ. 5.2 (p102)

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
  TH1I* h_mult_rec_neg = new TH1I("mult_rec_neg",";N_{-1} (rec)",10,0,10);
  TH1I* h_mult_rec_neu = new TH1I("mult_rec_neu",";N_{0} (rec)",10,0,10);

  // Starting with TDR histograms
  // 1a) Eta - MC particles
  TH1D* h_eta_MCp   = new TH1D("eta_MCp",";#eta_{p'}(MC)", 275, -11.0, 11.0);
  TH1D* h_eta_MCe   = new TH1D("eta_MCe",";#eta_{e'}(MC)", 275, -11.0, 11.0);
  TH1D* h_eta_MCg   = new TH1D("eta_MCg",";#eta_{#gamma}(MC)", 275, -11.0, 11.0);
  // 1a) Eta - MC particles
  TH1D* h_eta_MCAp   = new TH1D("eta_MCAp",";#eta_{p'}(MCA)", 275, -11.0, 11.0);
  TH1D* h_eta_MCAe   = new TH1D("eta_MCAe",";#eta_{e'}(MCA)", 275, -11.0, 11.0);
  TH1D* h_eta_MCAg   = new TH1D("eta_MCAg",";#eta_{#gamma}(MCA)", 275, -11.0, 11.0);
  // 1b) Eta - reco. particles
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
  // DEBUGGING
  TH1D* h_eta_FakePhot = new TH1D("eta_fakephot",";#eta_{e^{-} #rightarrow #gamma}(MC)", 275, -11.0, 11.0);
  TH2D* h_eta_ElecGamma_MC = new TH2D("eta_elecgamma_mc",";#eta_{e'}(MC);#eta_{#gamma}(MC)", 275, -11.0, 11.0, 275, -11.0, 11.0);
  TH2D* h_eta_ElecGamma_RP = new TH2D("eta_elecgamma_rp",";#eta_{e'}(Reco);#eta_{#gamma}(Reco)", 275, -11.0, 11.0, 275, -11.0, 11.0);

  // 2) Photon theta resolution
  TH1D* h_PhotRes_theta = new TH1D("photres_theta",";#theta_{#gamma}(Reco)-#theta_{#gamma}(MC) [deg]",360,-90,90);
  TH2D* h_PhotRes2D_theta = new TH2D("photres2d_theta",";#theta_{#gamma, MC} [deg]; #Delta#theta_{#gamma} [deg]",370,0,185,360,-90,90);

  // 3) t distribution
  TH1D* h_t_Truth  = new TH1D("t_truth" ,";|t|(MC) [(GeV/#it{c}^{2})^{2}]"           , 20, 0., 2.);
  TH1D* h_t_B0Acc  = new TH1D("t_b0acc" ,";|t|(MC|Reco. - B0) [(GeV/#it{c}^{2})^{2}]", 20, 0., 2.);
  TH1D* h_t_RPAcc  = new TH1D("t_rpacc" ,";|t|(MC|Reco. - RP) [(GeV/#it{c}^{2})^{2}]", 20, 0., 2.);
  TH1D* h_t_LCAcc  = new TH1D("t_lcacc" ,";|t_{e'#gamma}|(MC|Reco.) [(GeV/#it{c}^{2})^{2}]", 20, 0., 2.);
  TH1D* h_t_B0Reco = new TH1D("t_b0reco",";|t|(Reco. - B0) [(GeV/#it{c}^{2})^{2}]"   , 20, 0., 2.);
  TH1D* h_t_RPReco = new TH1D("t_rpreco",";|t|(Reco. - RP) [(GeV/#it{c}^{2})^{2}]"   , 20, 0., 2.);
  TH1D* h_t_LCReco = new TH1D("t_lcreco",";|t_{e''gamma}|(Reco.) [(GeV/#it{c}^{2})^{2}]"   , 20, 0., 2.);
  TH1D* h_t_Truth_Fine = new TH1D("t_truth_fine" ,"Counts/0.005 GeV^{2};|t|(MC) [(GeV/#it{c}^{2})^{2}]", 400, 0., 2.);
  
  // 4) t resolution - as absolute or as percentage (plot as preferred)
  TH2D* h_tResB0_2d = new TH2D("tresb0_2d",";|t|_{MC} [(GeV/#it{c})^{2}];#Deltat [(GeV/#it{c})^{2}]", 20, 0., 2., 500, -5., 5.);
  TH2D* h_tResRP_2d = new TH2D("tresrp_2d",";|t|_{MC} [(GeV/#it{c})^{2}];#Deltat [(GeV/#it{c})^{2}]", 20, 0., 2., 500, -5., 5.);
  TH2D* h_tResB0Pct_2d = new TH2D("tresb0pct_2d",";|t|_{MC} [(GeV/#it{c})^{2}];#Deltat/t_{MC}", 20, 0., 2., 200, -1., 1.);
  TH2D* h_tResRPPct_2d = new TH2D("tresrppct_2d",";|t|_{MC} [(GeV/#it{c})^{2}];#Deltat/t_{MC}", 20, 0., 2., 200, -1., 1.);
  TH2D* h_tResLC_2d    = new TH2D("treslc_2d",";|t_{e'#gamma}|_{MC} [(GeV/#it{c})^{2}];#Deltat [(GeV/#it{c})^{2}]", 20, 0., 2., 500, -5., 5.);
  TH2D* h_tResLCPct_2d = new TH2D("treslcpct_2d",";|t_{e'#gamma}|_{MC} [(GeV/#it{c})^{2}];#Deltat/t_{MC}", 20, 0., 2., 200, -1., 1.);

  // 5) Inclusive event kinematics - distributions, 2D response and resolution
  // 5a) 1D distributions
  TH1D* h_Q2_MC   = new TH1D("q2_mc"  , "Counts/0.2 GeV^{2};Q^{2}(MC) [GeV^{2}]"     , 550, 0., 110.);
  TH1D* h_Q2_Acc  = new TH1D("q2_acc" , "Counts/0.2 GeV^{2};Q^{2}(MC|Reco) [GeV^{2}]", 550, 0., 110.);
  TH1D* h_Q2_Reco = new TH1D("q2_reco", "Counts/0.2 GeV^{2};Q^{2}(Reco) [GeV^{2}]"   , 550, 0., 110.);
  TH1D* h_Q2_ExcReco = new TH1D("q2_excreco", "Counts/0.2 GeV^{2};Q^{2}(Reco) [GeV^{2}]"   , 550, 0., 110.);
  TH1D* h_xB_MC   = new TH1D("xb_mc"  , ";x_{B}(MC)"     , 2500, 0., 1.);
  TH1D* h_xB_Acc  = new TH1D("xb_acc" , ";x_{B}(MC|Reco)", 2500, 0., 1.);
  TH1D* h_xB_Reco = new TH1D("xb_reco", ";x_{B}(Reco)"   , 2500, 0., 1.);
  TH1D* h_y_MC   = new TH1D("y_mc"  , ";y(MC)"     , 100, 0., 1.);
  TH1D* h_y_Acc  = new TH1D("y_acc" , ";y(MC|Reco)", 100, 0., 1.);
  TH1D* h_y_Reco = new TH1D("y_reco", ";y(Reco)"   , 100, 0., 1.);
  // 5b) 2D distributions
  TH2D* h_Q2_2d = new TH2D("q2_2d",";Q^{2}(MC) [GeV^{2}];Q^{2}(Reco.) [GeV^{2}]", 500, 0., 100., 500, 0., 100.);
  TH2D* h_xB_2d = new TH2D("xb_2d",";x_{B}(MC);x_{B}(Reco.)", 2500, 0., 1., 2500, 0., 1.);
  TH2D* h_y_2d  = new TH2D("y_2d" ,";y(MC);y(Reco.)", 100, 0.,  1.,100, 0., 1.);
  // 5c) Resolutions
  TH2D* h_PctResQ2 = new TH2D("q2_pctres",";Q^{2} [GeV^{2}];#DeltaQ^{2}/Q^{2}", 500, 0., 100., 200, -1., 1.);
  TH2D* h_PctResxB = new TH2D("xb_pctres",";x_{B};#Deltax_{B}/x_{B}", 1e4, 0., 1., 1000, -5., 5.);
  TH2D* h_PctResy  = new TH2D("y_pctres" ,";y;#Deltay/y", 100, 0.,  1., 200, -1., 1.);
  // 5d) Cross-variable coverages
  TH2D* h_2D_xVQ2_MC = new TH2D("2d_xvq2_mc",";x_{B,MC};Q^{2}_{MC} [GeV^{2}]",1e4,0.,1.,200,0.,100.);
  TH2D* h_2D_xVQ2_RP = new TH2D("2d_xvq2_rp",";x_{B,Reco};Q^{2}_{Reco} [GeV^{2}]",1e4,0.,1.,200,0.,100.);
  TH2D* h_2D_xVt_MC = new TH2D("2d_xvt_mc",";x_{B,MC};|t|_{MC} [GeV^{2}]",5e4,0.,1.,100,0.,2.);
  TH2D* h_2D_xVt_RP = new TH2D("2d_xvt_rp",";x_{B,MC};|t|_{RP} [GeV^{2}]",5e4,0.,1.,100,0.,2.);
  TH2D* h_2D_xVtSin2_MC = new TH2D("2d_xvtsin2_mc",";x_{B,MC};|t|_{MC} [GeV^{2}]",5e4,0.,1.,100,0.,2.);
  TH2D* h_2D_xVtSin2_RP = new TH2D("2d_xvtsin2_rp",";x_{B,MC};|t|_{RP} [GeV^{2}]",5e4,0.,1.,100,0.,2.);
  TH2D* h_2D_xVtCos2_MC = new TH2D("2d_xvtcos2_mc",";x_{B,MC};|t|_{MC} [GeV^{2}]",5e4,0.,1.,100,0.,2.);
  TH2D* h_2D_xVtCos2_RP = new TH2D("2d_xvtcos2_rp",";x_{B,MC};|t|_{RP} [GeV^{2}]",5e4,0.,1.,100,0.,2.);

  // 6) Cut variables
  // 6a) Proton track theta
  TH1D* h_theta_MCp = new TH1D("theta_mcp", ";#theta_{p'}(MC) [mrad]", 200, 0., 50.);
  TH1D* h_theta_B0p = new TH1D("theta_b0p", ";#theta_{p'}(Reco. - B0) [mrad]", 200, 0., 50.);
  TH1D* h_theta_RPp = new TH1D("theta_RPp", ";#theta_{p'}(Reco. - RP) [mrad]", 200, 0., 50.);
  // 6b) Full event missing mass (squared)
  TH1D* h_M2miss3_MC  = new TH1D("m2miss3_mc" ,";M^{2}_{miss}(MC) [(GeV/#it{c}^{2})^{2}]",600,-150,150);
  TH1D* h_M2miss3_MCA = new TH1D("m2miss3_mca",";M^{2}_{miss}(MC|Reco) [(GeV/#it{c}^{2})^{2}]",600,-150,150);
  TH1D* h_M2miss3_RP  = new TH1D("m2miss3_rp" ,";M^{2}_{miss}(Reco) [(GeV/#it{c}^{2})^{2}]",600,-150,150);
  // 6c) E-pz (to reject radiative events)
  TH1D* h_EmPz3_MC = new TH1D("empz3_mc",";E-p_{z}(e'p'#gamma - MC) [GeV/#it{c}]",600,0.,150.);
  TH1D* h_EmPz3_RP = new TH1D("empz3_rp",";E-p_{z}(e'p'#gamma - Reco.) [GeV/#it{c}]",600,0.,150.);
  TH1D* h_EmPz2_MC = new TH1D("empz2_mc",";E-p_{z}(e'#gamma - MC) [GeV/#it{c}]",600,0.,150.);
  TH1D* h_EmPz2_RP = new TH1D("empz2_rp",";E-p_{z}(e'#gamma - MC) [GeV/#it{c}]",600,0.,150.);
  // 6d) Missing momentum (full final state)
  TH1D* h_Pmiss3_MC  = new TH1D("pmiss3_mc",";p_{miss}(MC) [GeV/#it{c}]",(Int_t)240*fPBeam_e,0.,(Int_t)1.2*fPBeam_e);
  TH1D* h_Pmiss3_MCA = new TH1D("pmiss3_mca",";p_{miss}(MCA) [GeV/#it{c}]",(Int_t)240*fPBeam_e,0.,(Int_t)1.2*fPBeam_e);
  TH1D* h_Pmiss3_RP  = new TH1D("pmiss3_rp",";p_{miss}(RP) [GeV/#it{c}]",(Int_t)240*fPBeam_e,0.,(Int_t)1.2*fPBeam_e);
  TH1D* h_Ptmiss3_MC  = new TH1D("ptmiss3_mc",";p_{T,miss}(MC) [GeV/#it{c}]",300,0.,3.);
  TH1D* h_Ptmiss3_MCA = new TH1D("ptmiss3_mca",";p_{T,miss}(MCA) [GeV/#it{c}]",300,0.,3.);
  TH1D* h_Ptmiss3_RP  = new TH1D("ptmiss3_rp",";p_{T,miss}(RP) [GeV/#it{c}]",300,0.,3.);
  // 6e) Missing momentum (semi-inclusive final state)
  Float_t protlower = (Float_t)0.8*fPBeam_p;
  Float_t protupper = (Float_t)1.2*fPBeam_p;
  Int_t nbins = protupper-protlower;
  TH1D* h_Pmiss2eg_MC  = new TH1D("pmiss2eg_mc",";p_{miss,e'#gamma}(MC) [GeV/#it{c}]",5*nbins,protlower,protupper);
  TH1D* h_Pmiss2eg_MCA = new TH1D("pmiss2eg_mca",";p_{miss,e'#gamma}(MCA) [GeV/#it{c}]",5*nbins,protlower,protupper);
  TH1D* h_Pmiss2eg_RP  = new TH1D("pmiss2eg_rp",";p_{miss,e'#gamma}(RP) [GeV/#it{c}]",5*nbins,protlower,protupper);
  TH1D* h_Ptmiss2eg_MC  = new TH1D("ptmiss2eg_mc",";p_{T,miss,e'#gamma}(MC) [GeV/#it{c}]",300,0.,3.);
  TH1D* h_Ptmiss2eg_MCA = new TH1D("ptmiss2eg_mca",";p_{T,miss,e'#gamma}(MCA) [GeV/#it{c}]",300,0.,3.);
  TH1D* h_Ptmiss2eg_RP  = new TH1D("ptmiss2eg_rp",";p_{T,miss,e'#gamma}(RP) [GeV/#it{c}]",300,0.,3.);  
  // 6e) Counting effect of cuts
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
  
  // 7) Trento phi
  TH1D* h_TPhi_MC    = new TH1D("tphi_mc"    ,";#phi_{h}(MC) [rad]"     , 18, 0, 360);
  TH1D* h_TPhi_B0Acc = new TH1D("tphi_b0acc" ,";#phi_{h}(MC|Reco) [rad]", 18, 0, 360);
  TH1D* h_TPhi_RPAcc = new TH1D("tphi_rpacc" ,";#phi_{h}(MC|Reco) [rad]", 18, 0, 360);
  TH1D* h_TPhi_B0Reco = new TH1D("tphi_b0reco",";#phi_{h}(Reco) [rad]"  , 18, 0, 360);
  TH1D* h_TPhi_RPReco = new TH1D("tphi_rpreco",";#phi_{h}(Reco) [rad]"  , 18, 0, 360);

  TH1D* h_TPhiRes = new TH1D("tphires",";#delta#phi_{h}", 160, -40, 40);

  // 8) 3D differential distributions: Trento phi  ->  [Q2][xB][|t|]
  // -------------------------------------------------------------------
  //  Binning is now read from an external text file (see DVCSBinning.hh /
  //  bins.txt) so that it can vary from one configuration to another. The
  //  number of Q2/xB/|t| bins is defined entirely by that file. The binning
  //  is region-dependent (B0 vs RP) and ragged (the |t| edges may depend on
  //  the xB slice).
  // -------------------------------------------------------------------
  DVCSBinning binning;
  // NOTE: change this if your binning file lives elsewhere or is named per
  //       configuration, e.g. "bins_"+camp+"_"+energy+".txt".
  TString sBinFile = "bins_"+sEnergy+".txt";
  if(!binning.load(sBinFile.Data())){
    std::cerr << "[ePIC_DVCS] FATAL: could not load binning from " << sBinFile << std::endl;
    return;
  }
  binning.print();
  const int nQ2bins = binning.nQ2();

  // Ragged histogram containers: [Q2][xB][|t|]
  typedef std::vector<std::vector<std::vector<TH1D*>>> H1Grid;

  // MC-truth Trento phi
  H1Grid h_TPhiDiff_B0MC, h_TPhiDiff_RPMC;
  // B0 region (B0 binning)
  H1Grid h_TPhiDiff_B0Acc, h_TPhiDiff_B0Reco, h_TPhiResDiff_B0;
  H1Grid h_Q2Diff_B0, h_xBDiff_B0, h_tDiff_B0;
  // RP region (RP binning)
  H1Grid h_TPhiDiff_RPAcc, h_TPhiDiff_RPReco, h_TPhiResDiff_RP;
  H1Grid h_Q2Diff_RP, h_xBDiff_RP, h_tDiff_RP;

  // xB:t 2D distribution - 1 per Q2 bin (continuous, unchanged)
  std::vector<TH2D*> h_2D_xVtDiff_RP(nQ2bins, nullptr);

  // Helper: book a ragged [Q2][xB][|t|] grid of TH1D for a given region.
  auto bookGrid = [&](H1Grid& grid, int reg, const char* nameBase, const char* yAxis,
                      int nYbins, double yLo, double yHi){
    grid.resize(nQ2bins);
    for(int q{0}; q<nQ2bins; q++){
      const int nx = binning.nXB(reg,q);
      grid[q].resize(nx);
      for(int x{0}; x<nx; x++){
        const int nt = binning.nT(reg,q,x);
        grid[q][x].resize(nt, nullptr);
        for(int t{0}; t<nt; t++){
          TString hname  = Form("%s[%i][%i][%i]", nameBase, q, x, t);
          TString htitle = Form("%.2e<Q^{2}<%.2e GeV^{2}, %.2e<x_{B}<%.2e, %.2f<|t|<%.2f;%s;",
                                binning.q2Low(q), binning.q2High(q),
                                binning.xBLow(reg,q,x), binning.xBHigh(reg,q,x),
                                binning.tLow(reg,q,x,t), binning.tHigh(reg,q,x,t), yAxis);
          grid[q][x][t] = new TH1D(hname, htitle, nYbins, yLo, yHi);
        }
      }
    }
  };

  // Book all differential families (y-axis binning kept as in the original code)
  bookGrid(h_TPhiDiff_B0MC,   DVCSBinning::kB0, "tphi_b0mc",   "#phi_{h}(MC - high-#it{t}) [deg]",10,  0.,  360.);
  bookGrid(h_TPhiDiff_B0Acc,  DVCSBinning::kB0, "tphi_b0acc",  "#phi_{h}(MC|Reco) [deg]",         10,  0.,  360.);
  bookGrid(h_TPhiDiff_B0Reco, DVCSBinning::kB0, "tphi_b0reco", "#phi_{h}(Reco) [deg]",            10,  0.,  360.);
  bookGrid(h_TPhiResDiff_B0,  DVCSBinning::kB0, "tphi_b0res",  "#delta#phi_{h}(Reco - B0) [deg]", 450, -45., 45.);
  bookGrid(h_TPhiDiff_RPMC,   DVCSBinning::kRP, "tphi_rpmc",   "#phi_{h}(MC - low-#it{t}) [deg]", 10,  0.,  360.);
  bookGrid(h_TPhiDiff_RPAcc,  DVCSBinning::kRP, "tphi_rpacc",  "#phi_{h}(MC|Reco) [deg]",         10,  0.,  360.);
  bookGrid(h_TPhiDiff_RPReco, DVCSBinning::kRP, "tphi_rpreco", "#phi_{h}(Reco) [deg]",            10,  0.,  360.);
  bookGrid(h_TPhiResDiff_RP,  DVCSBinning::kRP, "tphi_rpres",  "#delta#phi_{h}(Reco - RP) [deg]", 450, -45., 45.);

  bookGrid(h_Q2Diff_B0, DVCSBinning::kB0, "q2diff_b0", "Q^{2}(MCA) [GeV^{2}]", 550,   0., 110.);
  bookGrid(h_xBDiff_B0, DVCSBinning::kB0, "xbdiff_b0", "x_{B}",              10000,   0.,   1.);
  bookGrid(h_tDiff_B0,  DVCSBinning::kB0, "tdiff_b0",  "|t|(MCA) [GeV^{2}]",   20,   0.,   2.);
  bookGrid(h_Q2Diff_RP, DVCSBinning::kRP, "q2diff_rp", "Q^{2}(MCA) [GeV^{2}]", 550,   0., 110.);
  bookGrid(h_xBDiff_RP, DVCSBinning::kRP, "xbdiff_rp", "x_{B}",              10000,   0.,   1.);
  bookGrid(h_tDiff_RP,  DVCSBinning::kRP, "tdiff_rp",  "|t|(MCA) [GeV^{2}]",   20,   0.,   2.);

  // 2D xB:|t| per Q2 bin
  for(int q{0}; q<nQ2bins; q++)
    h_2D_xVtDiff_RP[q] = new TH2D(Form("xvtdiff_rp[%i]",q),";x_{B,MC};|t|_{RP} [GeV^{2}]",5e4,0.,1.,100,0.,2.);
  
  
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
    
      // Defin booleans for vetoes
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
      
      Int_t nRecNeg{0}, nRecNeu{0};

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
	      
	      nRecNeg++;

	      // Choose electrons from E/p 0.8 - 1.2
	      float clus_e{0};
	      for(auto& clust : mcreco.getRec().getClusters()) clus_e += clust.getEnergy();
	      float eoverp = clus_e/edm4hep::utils::magnitude(mcreco.getRec().getMomentum());
	      
	      h_EoverP_elec->Fill(eoverp);

	      // E-finder: veto on "non-electron" -ves
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
	      nRecNeu++;
	      
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
      
      if(nRecNeg > 9) nRecNeg = 9;
      if(nRecNeu > 9) nRecNeu = 9;
      h_mult_rec_neg->Fill(nRecNeg);
      h_mult_rec_neu->Fill(nRecNeu);

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

      //cout<<"[DEBUG] Particles found"<<endl;
      
      // VETOES
      // OMD
      const auto& OMDparts = event.get<edm4eic::ReconstructedParticleCollection>("ForwardOffMRecParticles");
      for(const auto& omdreco : OMDparts){
	// Check if track exists
	if(omdreco.getTracks()) kOtherFFTrack_Rec = true;
      }
      // ZDC
      const auto& ZDCparts = event.get<edm4eic::ReconstructedParticleCollection>("ReconstructedHcalFarForwardZDCNeutrals");
      for(const auto& zdcreco : ZDCparts){
	// Check for associated clusters
	if(zdcreco.getClusters()) kOtherFFTrack_Rec = true;
      }
      
      //cout<<"[DEBUG] Vetoes found"<<endl;

      //---------------------------------------------------------
      // Fill histograms
      //---------------------------------------------------------
      // Single-particle kinematics
      // Eta and inclusive Q2/xB/y

      // MC particles
      // Need Q2 for electron cuts
      if(scate4_gen.size() == 0) fQ2 = 0;
      else fQ2 = calcQ2_Elec(beame4, scate4_gen[0]);
      
      if(applyCuts_Electron(beame4,scate4_gen)){
	h_eta_MCe->Fill(scate4_gen[0].Eta());
	h_Q2_MC->Fill(fQ2);
	h_xB_MC->Fill(calcX_Elec(beame4, beamp4, scate4_gen[0]));
	h_y_MC->Fill(calcY_Elec(beame4, beamp4, scate4_gen[0]));

	h_2D_EvEta_e->Fill(scate4_gen[0].Eta(), scate4_gen[0].E());
	h_2D_xVQ2_MC->Fill(calcX_Elec(beame4, beamp4, scate4_gen[0]), fQ2);
      }
      if(applyCuts_Photon(scatg4_gen)){
	h_eta_MCg->Fill(scatg4_gen[0].Eta());
	h_2D_EvEta_g->Fill(scatg4_gen[0].Eta(), scatg4_gen[0].E());
      }
      if(applyCuts_Proton(scatp4_gen, "all")){
	h_eta_MCp->Fill(scatp4_gen[0].Eta());
	h_2D_EvEta_p->Fill(scatp4_gen[0].Eta(), scatp4_gen[0].E());
      }

      // For cut histogram - don't apply proton theta cut
      if(scatp4_gen.size() == 1){
	h_theta_MCp->Fill(scatp4_gen[0].Theta()*1000);
      }
      if(applyCuts_Electron(beame4,scate4_gen) && applyCuts_Photon(scatg4_gen)){
	h_eta_ElecGamma_MC->Fill(scate4_gen[0].Eta(),scatg4_gen[0].Eta());
	h_EmPz2_MC->Fill((scate4_gen[0]+scatg4_gen[0]).E() - (scate4_gen[0]+scatg4_gen[0]).Pz());
      }
	
      if(applyCuts_Electron(beame4,scate4_aso)) h_eta_MCAe->Fill(scate4_aso[0].Eta());
      if(applyCuts_Photon(scatg4_aso))	        h_eta_MCAg->Fill(scatg4_gen[0].Eta());
      if(applyCuts_Proton(scatp4_aso, "all"))	h_eta_MCAp->Fill(scatp4_aso[0].Eta());
      
      // Reconstructed particles and associated MC
      // Need Q2 for electron cuts
      if(scate4_rec.size() == 0) fQ2 = 0;
      else fQ2 = calcQ2_Elec(beame4, scate4_rec[0]);

      if(scate4_rec.size() == 1) h_Q2_Reco->Fill(fQ2);
      if(applyCuts_Electron(beame4,scate4_rec) && applyCuts_Electron(beame4,scate4_aso)){
	h_eta_RPe->Fill(scate4_rec[0].Eta());
	h_Q2_Acc->Fill(calcQ2_Elec(beame4, scate4_aso[0]));
	h_xB_Acc->Fill(calcX_Elec(beame4, beamp4, scate4_aso[0]));
	h_xB_Reco->Fill(calcX_Elec(beame4, beamp4, scate4_rec[0]));
	h_y_Acc->Fill(calcY_Elec(beame4, beamp4, scate4_aso[0]));
	h_y_Reco->Fill(calcY_Elec(beame4, beamp4, scate4_rec[0]));

	h_2D_xVQ2_RP->Fill(calcX_Elec(beame4, beamp4, scate4_rec[0]), fQ2);
      }
      if(applyCuts_Photon(scatg4_rec))       h_eta_RPg->Fill(scatg4_rec[0].Eta());
      if(applyCuts_Proton(scatp4_rec, "B0")) h_eta_RPp->Fill(scatp4_rec[0].Eta());
      if(applyCuts_Proton(scatp4_rom, "RP")) h_eta_RPPp->Fill(scatp4_rom[0].Eta());
      // For cut histogram - don't apply proton theta cut
      if(scatp4_rec.size() == 1 && scatp4_rom.size() == 0) h_theta_B0p->Fill(scatp4_rec[0].Theta()*1000);
      if(scatp4_rom.size() == 1 && scatp4_rec.size() == 0) h_theta_RPp->Fill(scatp4_rom[0].Theta()*1000);
      if(applyCuts_Electron(beame4,scate4_rec) && applyCuts_Photon(scatg4_rec)) h_eta_ElecGamma_RP->Fill(scate4_rec[0].Eta(),scatg4_rec[0].Eta());
      
      // 2D inclusive distributions
      if(applyCuts_Electron(beame4,scate4_rec) && applyCuts_Electron(beame4,scate4_aso)){
	double q2_mc = calcQ2_Elec(beame4, scate4_aso[0]);
	double q2_rec = calcQ2_Elec(beame4, scate4_rec[0]);
	double q2_frac = (q2_rec-q2_mc)/q2_mc;
	h_Q2_2d->Fill(q2_mc,q2_rec);
	h_PctResQ2->Fill(q2_mc, q2_frac);

	double xB_mc = calcX_Elec(beame4, beamp4, scate4_aso[0]);
	double xB_rec = calcX_Elec(beame4, beamp4, scate4_rec[0]);
	double xB_frac = (xB_rec-xB_mc)/xB_mc;
	h_xB_2d->Fill(xB_mc,xB_rec);
	h_PctResxB->Fill(xB_mc, xB_frac);

	double y_mc = calcY_Elec(beame4, beamp4, scate4_aso[0]);
	double y_rec = calcY_Elec(beame4, beamp4, scate4_rec[0]);
	double y_frac = (y_rec-y_mc)/y_mc;
	h_y_2d->Fill(y_mc,y_rec);
	h_PctResy->Fill(y_mc, y_frac);
      }

      // Photon theta resolution
      Float_t th_rec{0}, th_gen{0};
      if(applyCuts_Photon(scatg4_rec) && applyCuts_Photon(scatg4_aso)){
	th_gen = scatg4_aso[0].Theta()*TMath::RadToDeg();
	th_rec = scatg4_rec[0].Theta()*TMath::RadToDeg();
	
	h_PhotRes_theta->Fill(th_rec-th_gen);
	h_PhotRes2D_theta->Fill(th_gen, th_rec-th_gen);
      }

      //cout<<"[DEBUG] Inclusive histos filled"<<endl;

      // Full DVCS event distributions - Mandelstam t and Trento phi
      // MC truth
      if(applyCuts_Electron(beame4,scate4_gen) 
	 && applyCuts_Photon(scatg4_gen) 
	 && applyCuts_Proton(scatp4_gen, "all")) h_EmPz3_MC->Fill((scate4_gen[0]+scatp4_gen[0]+scatg4_gen[0]).E() - (scate4_gen[0]+scatp4_gen[0]+scatg4_gen[0]).Pz());
	 
      //cout<<"[DEBUG] MC E-Pz (full evt.) filled"<<endl;

      if(applyCuts_All(beame4, beamp4, scate4_gen, scatp4_gen, scatg4_gen, "all")      // Single-particle and exclusive event kinematic cuts
	 && !kNonElecNeg_MC && !kOtherFFTrack_MC && !kBarrelPos_MC){                   // Vetoes
	// Calculations
	Float_t t_gen = calcT_BABE(beamp4,scatp4_gen[0]);
	Float_t tphi_gen = calcTrentoPhi_qg(beame4, beamp4, scate4_gen[0], scatg4_gen[0]);
	// Distributions - Q2/xB integrated
	h_t_Truth->Fill(t_gen);
	h_t_Truth_Fine->Fill(t_gen);
	
	h_TPhi_MC->Fill(tphi_gen*TMath::RadToDeg());

	// Q2/xB diff.
	// Need Q2 and xB for events
	Float_t q2_gen = calcQ2_Elec(beame4, scate4_gen[0]);
	Float_t xB_gen = calcX_Elec(beame4, beamp4, scate4_gen[0]);
	
	// Find bin numbers for RP scheme
	int binq2 = binning.findQ2(q2_gen);
	int binxB = binning.findXB(DVCSBinning::kRP, binq2, xB_gen);
	int bint  = binning.findT (DVCSBinning::kRP, binq2, binxB, t_gen);
	if(binq2>=0 && binxB>=0 && bint>=0) h_TPhiDiff_RPMC[binq2][binxB][bint]->Fill(tphi_gen*TMath::RadToDeg());

	// Find bin numbers for B0 scheme
	binxB = binning.findXB(DVCSBinning::kB0, binq2, xB_gen);
	bint  = binning.findT (DVCSBinning::kB0, binq2, binxB, t_gen);
	if(binq2>=0 && binxB>=0 && bint>=0) h_TPhiDiff_B0MC[binq2][binxB][bint]->Fill(tphi_gen*TMath::RadToDeg());
	
	h_2D_xVt_MC->Fill(xB_gen,t_gen);
	h_2D_xVtSin2_MC->Fill(xB_gen,t_gen,TMath::Power(TMath::Sin(tphi_gen),2));
	h_2D_xVtCos2_MC->Fill(xB_gen,t_gen,TMath::Power(TMath::Cos(tphi_gen),2));
       
	h_Pmiss3_MC->Fill(calcPMiss_3Body(beame4, beamp4, scate4_gen[0], scatp4_gen[0], scatg4_gen[0]));
	h_Ptmiss3_MC->Fill(calcPtMiss_3Body(beame4, beamp4, scate4_gen[0], scatp4_gen[0], scatg4_gen[0]));
      	
	Float_t pmiss = calcPMiss_2Body(beame4, beamp4, scate4_gen[0], scatg4_gen[0]);
	Float_t ptmiss = calcPtMiss_2Body(beame4, beamp4, scate4_gen[0], scatg4_gen[0]);
	h_Pmiss2eg_MC->Fill(pmiss);
	h_Ptmiss2eg_MC->Fill(ptmiss);
      }
      
      //cout<<"[DEBUG] MC TPhi filled"<<endl;
 
      // Reconstructed and MC accepted - B0 only
      if(applyCuts_Electron(beame4,scate4_rec) && applyCuts_Photon(scatg4_rec) && applyCuts_Proton(scatp4_rec, "B0") 
	 && scatp4_rom.size() == 0)
	 h_EmPz3_RP->Fill((scate4_rec[0]+scatp4_rec[0]+scatg4_rec[0]).E() - (scate4_rec[0]+scatp4_rec[0]+scatg4_rec[0]).Pz());

      //cout<<"[DEBUG] Reco. E-Pz (full evt., B0 proton) filled"<<endl;

      if(applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rec, scatg4_rec, "B0") && scatp4_rom.size() == 0     // Single-particle and exclusive event kinematic cuts
	 && !kNonElecNeg_Rec && !kOtherFFTrack_Rec && !kBarrelPos_Rec){                                        // Vetoes
	// Calculations
	Float_t t_acc = calcT_BABE(beamp4,scatp4_aso[0]);
	Float_t t_rec = calcT_BABE(beamp4,scatp4_rec[0]);
	Float_t tphi_acc = calcTrentoPhi_qg(beame4, beamp4, scate4_aso[0], scatg4_aso[0]);
	Float_t tphi_rec = calcTrentoPhi_qg(beame4, beamp4, scate4_rec[0], scatg4_rec[0]);
	// Distributions - Q2/xB integrated
	h_t_B0Acc->Fill(t_acc);
	h_t_B0Reco->Fill(t_rec);
	h_TPhi_B0Acc->Fill(tphi_acc*TMath::RadToDeg());
	h_TPhi_B0Reco->Fill(tphi_rec*TMath::RadToDeg());

	h_TPhiRes->Fill((tphi_rec-tphi_acc)*TMath::RadToDeg());

	// Q2/xB diff. - accepted
	// Need Q2 and xB for events
	Float_t q2_acc = calcQ2_Elec(beame4, scate4_aso[0]);
	Float_t xB_acc = calcX_Elec(beame4, beamp4, scate4_aso[0]);

	// Find bin numbers (B0 region binning)
	int binq2 = binning.findQ2(q2_acc);
	int binxB = binning.findXB(DVCSBinning::kB0, binq2, xB_acc);
	int bint  = binning.findT (DVCSBinning::kB0, binq2, binxB, t_acc);
	if(binq2>=0 && binxB>=0 && bint>=0) h_TPhiDiff_B0Acc[binq2][binxB][bint]->Fill(tphi_acc*TMath::RadToDeg());
	
	// Q2/xB diff. - reconstructed
	// Need Q2 and xB for events
	Float_t q2_rec = calcQ2_Elec(beame4, scate4_rec[0]);
	Float_t xB_rec = calcX_Elec(beame4, beamp4, scate4_rec[0]);

	if(binq2>=0 && binxB>=0 && bint>=0){
	  h_TPhiDiff_B0Reco[binq2][binxB][bint]->Fill(tphi_rec*TMath::RadToDeg());
	  h_TPhiResDiff_B0[binq2][binxB][bint]->Fill((tphi_rec-tphi_acc)*TMath::RadToDeg());
	  // Fill MCA histograms based on reco. kinematics
	  h_Q2Diff_B0[binq2][binxB][bint]->Fill(q2_acc);
	  h_xBDiff_B0[binq2][binxB][bint]->Fill(xB_acc);
	  h_tDiff_B0[binq2][binxB][bint]->Fill(t_acc);
	}

	h_Q2_ExcReco->Fill(calcQ2_Elec(beame4, scate4_rec[0]), TMath::Power(TMath::Sin(tphi_rec),2));
	
	h_2D_xVt_RP->Fill(xB_rec,t_rec);
	if(binq2>=0) h_2D_xVtDiff_RP[binq2]->Fill(xB_rec,t_rec,TMath::Power(TMath::Sin(tphi_rec),2));
	h_2D_xVtSin2_RP->Fill(xB_rec,t_rec,TMath::Power(TMath::Sin(tphi_rec),2));
	h_2D_xVtCos2_RP->Fill(xB_rec,t_rec,TMath::Power(TMath::Cos(tphi_rec),2));

	h_Pmiss3_MCA->Fill(calcPMiss_3Body(beame4, beamp4, scate4_aso[0], scatp4_aso[0], scatg4_aso[0]));
	h_Ptmiss3_MCA->Fill(calcPtMiss_3Body(beame4, beamp4, scate4_aso[0], scatp4_aso[0], scatg4_aso[0]));
	h_Pmiss3_RP->Fill(calcPMiss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rec[0], scatg4_rec[0]));
	h_Ptmiss3_RP->Fill(calcPtMiss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rec[0], scatg4_rec[0]));
      }
      
      //cout<<"[DEBUG] TPhi and differential histos. (B0 proton) filled"<<endl;

      // Reconstructed and accepted - RP only
      if(applyCuts_Electron(beame4,scate4_rec) && applyCuts_Photon(scatg4_rec) && applyCuts_Proton(scatp4_rom, "RP") 
	 && scatp4_rec.size() == 0)
	 h_EmPz3_RP->Fill((scate4_rec[0]+scatp4_rom[0]+scatg4_rec[0]).E() - (scate4_rec[0]+scatp4_rom[0]+scatg4_rec[0]).Pz());

      //cout<<"[DEBUG] Reco. E-Pz (full evt., RP proton) filled"<<endl;

      if(applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rom, scatg4_rec, "RP") && scatp4_rec.size() == 0      // Single-particle and exclusive event kinematic cuts
	 && !kNonElecNeg_Rec && !kOtherFFTrack_Rec && !kBarrelPos_Rec                                           // Vetoes
	 && scatp4_gen.size() > 0 && scate4_aso.size() > 0 && scatg4_aso.size() > 0){                           // Ensure associations exist
       	// Calculations
	Float_t t_acc = calcT_BABE(beamp4,scatp4_gen[0]);
	Float_t t_rec = calcT_BABE(beamp4,scatp4_rom[0]);
	Float_t tphi_acc = calcTrentoPhi_qg(beame4, beamp4, scate4_aso[0], scatg4_aso[0]);
	Float_t tphi_rec = calcTrentoPhi_qg(beame4, beamp4, scate4_rec[0], scatg4_rec[0]);
	// Distributions - Q2/xB integrated
       	h_t_RPAcc->Fill(t_acc);
	h_t_RPReco->Fill(t_rec);
	h_TPhi_RPAcc->Fill(tphi_acc*TMath::RadToDeg());
	h_TPhi_RPReco->Fill(tphi_rec*TMath::RadToDeg());	
	
	h_TPhiRes->Fill((tphi_rec-tphi_acc)*TMath::RadToDeg());

	//cout<<"[DEBUG] 3D integrated histos filled"<<endl;

	// Q2/xB diff. - accepted
	// Need Q2 and xB for events
	Float_t q2_acc = calcQ2_Elec(beame4, scate4_aso[0]);
	Float_t xB_acc = calcX_Elec(beame4, beamp4, scate4_aso[0]);

	// Find bin numbers (RP region binning)
	int binq2 = binning.findQ2(q2_acc);
	int binxB = binning.findXB(DVCSBinning::kRP, binq2, xB_acc);
	int bint  = binning.findT (DVCSBinning::kRP, binq2, binxB, t_acc);
	if(binq2>=0 && binxB>=0 && bint>=0) h_TPhiDiff_RPAcc[binq2][binxB][bint]->Fill(tphi_acc*TMath::RadToDeg());

	// Q2/xB diff. - reconstructed
	Float_t q2_rec = calcQ2_Elec(beame4, scate4_rec[0]);
	Float_t xB_rec = calcX_Elec(beame4, beamp4, scate4_rec[0]);
	
	if(binq2>=0 && binxB>=0 && bint>=0){
	  h_TPhiDiff_RPReco[binq2][binxB][bint]->Fill(tphi_rec*TMath::RadToDeg());
	  h_TPhiResDiff_RP[binq2][binxB][bint]->Fill((tphi_rec-tphi_acc)*TMath::RadToDeg());
	  // Fill MCA histograms based on reco. kinematics
	  h_Q2Diff_RP[binq2][binxB][bint]->Fill(q2_acc);
	  h_xBDiff_RP[binq2][binxB][bint]->Fill(xB_acc);
	  h_tDiff_RP[binq2][binxB][bint]->Fill(t_acc);
	}
	
	//cout<<"[DEBUG] 3D diff. histos filled"<<endl;
 
	//cout<<xB_rec<<"\t"<<t_rec<<"\t"<<binq2<<"\t"<<q2_acc<<endl;
	
	h_Q2_ExcReco->Fill(calcQ2_Elec(beame4, scate4_rec[0]), TMath::Power(TMath::Sin(tphi_rec),2));

	h_2D_xVt_RP->Fill(xB_rec,t_rec);
	if(binq2>=0) h_2D_xVtDiff_RP[binq2]->Fill(xB_rec,t_rec,TMath::Power(TMath::Sin(tphi_rec),2));
	h_2D_xVtSin2_RP->Fill(xB_rec,t_rec,TMath::Power(TMath::Sin(tphi_rec),2));
	h_2D_xVtCos2_RP->Fill(xB_rec,t_rec,TMath::Power(TMath::Cos(tphi_rec),2));

	//cout<<"[DEBUG] 2D x:t distribution filled"<<endl;

	h_Pmiss3_MCA->Fill(calcPMiss_3Body(beame4, beamp4, scate4_aso[0], scatp4_gen[0], scatg4_aso[0]));
	h_Ptmiss3_MCA->Fill(calcPtMiss_3Body(beame4, beamp4, scate4_aso[0], scatp4_gen[0], scatg4_aso[0]));
	h_Pmiss3_RP->Fill(calcPMiss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rom[0], scatg4_rec[0]));
	h_Ptmiss3_RP->Fill(calcPtMiss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rom[0], scatg4_rec[0]));
      
	//cout<<"[DEBUG] Missing kinematics filled"<<endl;
      }
      
      //cout<<"[DEBUG] TPhi and differential histos. (RP proton) filled"<<endl;

      // Semi-inclusive calculation - ignore if proton is detected or not
      if(applyCuts_Electron(beame4,scate4_rec) && applyCuts_Photon(scatg4_rec))	h_EmPz2_RP->Fill((scate4_rec[0]+scatg4_rec[0]).E() - (scate4_rec[0]+scatg4_rec[0]).Pz());

      //cout<<"[DEBUG] Reco. E-Pz (semi-inclusive) filled"<<endl;

      if(applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rec, scatg4_rec, "semi")     // Single-particle and semi-inclusive event kinematic cuts
	 && !kNonElecNeg_Rec && !kOtherFFTrack_Rec && !kBarrelPos_Rec){                // Vetoes
	// Calculations
	Float_t t_acc = calcT_MethodL(beame4,beamp4,scate4_aso[0],fMass_proton,scatg4_aso[0]);
	Float_t t_rec = calcT_MethodL(beame4,beamp4,scate4_rec[0],fMass_proton,scatg4_rec[0]);
	
	// Distributions - Q2/xB integrated
	h_t_LCAcc->Fill(t_acc);
	h_t_LCReco->Fill(t_rec);
	
	// Q2/xB diff. - accepted
	Float_t q2_acc = calcQ2_Elec(beame4, scate4_aso[0]);
	Float_t xB_acc = calcX_Elec(beame4, beamp4, scate4_aso[0]);

	// Q2/xB diff. - reconstructed
	Float_t q2_rec = calcQ2_Elec(beame4, scate4_rec[0]);
	Float_t xB_rec = calcX_Elec(beame4, beamp4, scate4_rec[0]);
	
	h_Pmiss2eg_MCA->Fill(calcPMiss_2Body(beame4, beamp4, scate4_aso[0], scatg4_aso[0]));
	h_Ptmiss2eg_MCA->Fill(calcPtMiss_2Body(beame4, beamp4, scate4_aso[0], scatg4_aso[0]));
	h_Pmiss2eg_RP->Fill(calcPMiss_2Body(beame4, beamp4, scate4_rec[0], scatg4_rec[0]));
	h_Ptmiss2eg_RP->Fill(calcPtMiss_2Body(beame4, beamp4, scate4_rec[0], scatg4_rec[0]));
      }
      
      //cout<<"[DEBUG] Semi-inclusive kinematics filled"<<endl;

      // Mandelstam t-resolution
      // ASSUME THAT GENERATED POSITIVE TRACK MATCHES RECONSTRUCTED POSITIVE TRACK
      Float_t t_rec{0}, t_gen{0};
      if(applyCuts_All(beame4, beamp4, scate4_gen, scatp4_gen, scatg4_gen, "all") 
	 && applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rec, scatg4_rec, "B0") && scatp4_rom.size()==0){
	t_gen = calcT_BABE(beamp4, scatp4_gen[0]);
	t_rec = calcT_BABE(beamp4, scatp4_rec[0]);
	h_tResB0_2d->Fill(t_gen, TMath::Abs(t_rec-t_gen));
	h_tResB0Pct_2d->Fill(t_gen, TMath::Abs(t_rec-t_gen)/t_gen);
      }
      if(applyCuts_All(beame4, beamp4, scate4_gen, scatp4_gen, scatg4_gen, "all") 
	 && applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rom, scatg4_rec, "RP") && scatp4_rec.size()==0){
	t_gen = calcT_BABE(beamp4, scatp4_gen[0]);
	t_rec = calcT_BABE(beamp4, scatp4_rom[0]);
	h_tResRP_2d->Fill(t_gen, TMath::Abs(t_rec-t_gen));
	h_tResRPPct_2d->Fill(t_gen, TMath::Abs(t_rec-t_gen)/t_gen);
      }
      if(applyCuts_All(beame4, beamp4, scate4_gen, scatp4_gen, scatg4_gen, "semi")
	 && applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rec, scatg4_rec, "semi")){
	 // applyCuts_Electron(beame4,scate4_rec) && applyCuts_Photon(scatg4_rec)
	 //&& applyCuts_Electron(beame4,scate4_gen) && applyCuts_Photon(scatg4_gen)){
	t_gen = calcT_MethodL(beame4,beamp4,scate4_gen[0],fMass_proton,scatg4_gen[0]);
	t_rec = calcT_MethodL(beame4,beamp4,scate4_rec[0],fMass_proton,scatg4_rec[0]);
	h_tResLC_2d->Fill(t_gen, TMath::Abs(t_rec-t_gen));
	h_tResLCPct_2d->Fill(t_gen, TMath::Abs(t_rec-t_gen)/t_gen);
      }

      //cout<<"[DEBUG] Mandelstam t filled"<<endl;

      // Missing mass histos - only apply multiplicity and track quality cuts
      if(applyCuts_Electron(beame4,scate4_gen) && 
	 applyCuts_Photon(scatg4_gen)          &&
	 applyCuts_Proton(scatp4_gen, "all")     ) h_M2miss3_MC->Fill(calcM2Miss_3Body(beame4, beamp4, scate4_gen[0], scatp4_gen[0], scatg4_gen[0]));
      if(applyCuts_Electron(beame4,scate4_aso) && 
	 applyCuts_Photon(scatg4_aso)          &&
	 applyCuts_Proton(scatp4_aso, "all")     ) h_M2miss3_MCA->Fill(calcM2Miss_3Body(beame4, beamp4, scate4_aso[0], scatp4_aso[0], scatg4_aso[0]));
      if(applyCuts_Electron(beame4,scate4_rec) && 
	 applyCuts_Photon(scatg4_rec)          &&
	 applyCuts_Proton(scatp4_rec, "B0")      ) h_M2miss3_RP->Fill(calcM2Miss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rec[0], scatg4_rec[0]));
      if(applyCuts_Electron(beame4,scate4_rec) && 
	 applyCuts_Photon(scatg4_rec)          &&
	 applyCuts_Proton(scatp4_rom, "RP")      ) h_M2miss3_RP->Fill(calcM2Miss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rom[0], scatg4_rec[0]));

      //cout<<"[DEBUG] Missing kinematics kinematics filled"<<endl;

      // Count no. of events which pass cuts - MC only
      // ---------------------Single particle---------------------
      // Electron
      if(scate4_gen.size() == 1){
        hPassCuts_MC->Fill(0);

        // ...and Q2 > 1 GeV2
	fQ2 = calcQ2_Elec(beame4, scate4_gen[0]);
        if(fQ2 >= 1.) hPassCuts_MC->Fill(1); //fi (Q2 cut)
      } //fi (electrons inclusive)
      // Photon
      if(scatg4_gen.size() == 1){
        hPassCuts_MC->Fill(2);

        // ...and E_{gamma} > 1 GeV
	if(scatg4_gen[0].E() >= 1.) hPassCuts_MC->Fill(3); //fi (E_gamma cut)
      }//fi (photons inclusive)
      // Proton
      if(scatp4_gen.size() == 1){
        hPassCuts_MC->Fill(4);

        // ...and proton track theta cut (FF region)
	if(scatp4_gen[0].Theta() <= 0.02) hPassCuts_MC->Fill(5); //fi (p'theta cut)
      } //fi (protons inclusive)
      // ------------------Event---------------------
      // Full DVCS multiplicity...
      if(scate4_gen.size() == 1 && scatg4_gen.size() == 1 && scatp4_gen.size() == 1){
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
      if(scate4_gen.size() == 1 && scatg4_gen.size() == 1){
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
      
      //cout<<"[DEBUG] Cut counter (MC) filled"<<endl;

      // Count no. of events which pass cuts - reconstructed only
      // ---------------------Single particle---------------------
      // Electron
      if(scate4_rec.size() == 1){
        hPassCuts_Rec->Fill(0);

        // ...and Q2 > 1 GeV2
	fQ2 = calcQ2_Elec(beame4, scate4_rec[0]);
        if(fQ2 >= 1.) hPassCuts_Rec->Fill(1); //fi (Q2 cut)
      } //fi (electrons inclusive)
      // Photon
      if(scatg4_rec.size() == 1){
        hPassCuts_Rec->Fill(2);

        // ...and E_{gamma} > 1 GeV
	if(scatg4_rec[0].E() >= 1.) hPassCuts_Rec->Fill(3); //fi (E_gamma cut)
      }//fi (photons inclusive)
      // Proton - B0
      if(scatp4_rec.size() == 1 && scatp4_rom.size() == 0){
        hPassCuts_Rec->Fill(4);

        // ...and proton track theta cut (B0 region)
	if(scatp4_rec[0].Theta() >= 0.0055 && scatp4_rec[0].Theta() <= 0.02) hPassCuts_Rec->Fill(5); //fi (p'theta cut)
      } //fi (B0 proton)
      // Proton - RP
      if(scatp4_rom.size() == 1 && scatp4_rec.size() == 0){
        hPassCuts_Rec->Fill(4);

        // ...and proton track theta cut (RP region)
	if(scatp4_rom[0].Theta() > 0 && scatp4_rom[0].Theta() <= 0.005) hPassCuts_Rec->Fill(5); //fi (p'theta cut)
      } //fi (RP proton)
      // ------------------Event---------------------
      // Full DVCS multiplicity (using B0 protons)...
      if(scate4_rec.size() == 1 && scatg4_rec.size() == 1 && scatp4_rec.size() == 1 && scatp4_rom.size() == 0){
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
      if(scate4_rec.size() == 1 && scatg4_rec.size() == 1 && scatp4_rom.size() == 1 && scatp4_rec.size() == 0){
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
      if(scate4_rec.size() == 1 && scatg4_rec.size() == 1){
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
      
      //cout<<"[DEBUG] Cut counter (rec.) filled"<<endl;


    }// END OF EVENT LOOP

    fileCounter++;
  } // END OF FILE LIST
  
  
  //------------------------------------------------------------
  // Write to output file
  //------------------------------------------------------------
  TString sSimple = sOutFileName+"_simple.root";
  TFile* fOutSimple = new TFile(sSimple,"RECREATE");
  fOutSimple->cd();
  
  h_mult_rec_neg->Write();
  h_mult_rec_neu->Write();
  // Eta - MC
  h_eta_MCp->Write();
  h_eta_MCe->Write();
  h_eta_MCg->Write();
  h_eta_MCAp->Write();
  h_eta_MCAe->Write();
  h_eta_MCAg->Write();
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
  // Photon theta resolution
  h_PhotRes_theta->Write();
  h_PhotRes2D_theta->Write();
  // t distributions
  h_t_Truth->Write();
  h_t_Truth_Fine->Write();
  h_t_B0Acc->Write();
  h_t_RPAcc->Write();
  h_t_B0Reco->Write();
  h_t_RPReco->Write();
  h_t_LCAcc->Write();
  h_t_LCReco->Write();
  // 2D t resolution
  h_tResB0_2d->Write();
  h_tResRP_2d->Write();
  h_tResB0Pct_2d->Write();
  h_tResRPPct_2d->Write();
  h_tResLC_2d->Write();
  h_tResLCPct_2d->Write();
  // Trento phi
  h_TPhi_MC->Write();
  h_TPhi_B0Acc->Write();
  h_TPhi_RPAcc->Write();
  h_TPhi_B0Reco->Write();
  h_TPhi_RPReco->Write();
  h_TPhiRes->Write();
  // Inclusive kinematic distributions - 1D
  h_Q2_MC->Write();
  h_Q2_Acc->Write();
  h_Q2_Reco->Write();
  h_Q2_ExcReco->Write();
  h_xB_MC->Write();
  h_xB_Acc->Write();
  h_xB_Reco->Write();
  h_y_MC->Write();
  h_y_Acc->Write();
  h_y_Reco->Write();
  // Inclusive kinematic distributions - 2D
  h_Q2_2d->Write();
  h_xB_2d->Write();
  h_y_2d->Write();
  // Inclusive kinematic resolutions
  h_PctResQ2->Write();
  h_PctResxB->Write();
  h_PctResy->Write();
  // Cross-variable coverages
  h_2D_xVQ2_MC->Write();
  h_2D_xVQ2_RP->Write();
  h_2D_xVt_MC->Write();
  h_2D_xVt_RP->Write();
  // Cuts plots
  h_theta_MCp->Write();
  h_theta_B0p->Write();
  h_theta_RPp->Write();
  h_M2miss3_MC->Write();
  h_M2miss3_MCA->Write();
  h_M2miss3_RP->Write();
  h_Pmiss3_MC->Write();
  h_Pmiss3_MCA->Write();
  h_Pmiss3_RP->Write();
  h_Ptmiss3_MC->Write();
  h_Ptmiss3_MCA->Write();
  h_Ptmiss3_RP->Write();
  h_Pmiss2eg_MC->Write();
  h_Ptmiss2eg_MC->Write();
  h_Pmiss2eg_MCA->Write();
  h_Ptmiss2eg_MCA->Write();
  h_Pmiss2eg_RP->Write();
  h_Ptmiss2eg_RP->Write();
  h_EmPz3_MC->Write();
  h_EmPz3_RP->Write();
  h_EmPz2_MC->Write();
  h_EmPz2_RP->Write();

  h_eta_FakePhot->Write();
  h_eta_ElecGamma_MC->Write();
  h_eta_ElecGamma_RP->Write();

  hPassCuts_MC->Write();
  hPassCuts_Rec->Write();
  
  fOutSimple->Close();

  TString sDiff = sOutFileName+"_diff.root";
  TFile* fOutDiff = new TFile(sDiff,"RECREATE");
  fOutDiff->cd();

  for(int q{0}; q<nQ2bins; q++){
    // B0 region grids
    for(int x{0}; x<binning.nXB(DVCSBinning::kB0,q); x++){
      for(int t{0}; t<binning.nT(DVCSBinning::kB0,q,x); t++){
	h_Q2Diff_B0[q][x][t]->Write();
	h_xBDiff_B0[q][x][t]->Write();
	h_tDiff_B0[q][x][t]->Write();
	h_TPhiDiff_B0MC[q][x][t]->Write();
	h_TPhiDiff_B0Acc[q][x][t]->Write();
	h_TPhiDiff_B0Reco[q][x][t]->Write();
	h_TPhiResDiff_B0[q][x][t]->Write();
      }
    }
    // RP region grids (and MC truth, which uses the RP scheme)
    for(int x{0}; x<binning.nXB(DVCSBinning::kRP,q); x++){
      for(int t{0}; t<binning.nT(DVCSBinning::kRP,q,x); t++){
	h_Q2Diff_RP[q][x][t]->Write();
	h_xBDiff_RP[q][x][t]->Write();
	h_tDiff_RP[q][x][t]->Write();
	h_TPhiDiff_RPMC[q][x][t]->Write();
	h_TPhiDiff_RPAcc[q][x][t]->Write();
	h_TPhiDiff_RPReco[q][x][t]->Write();
	h_TPhiResDiff_RP[q][x][t]->Write();
      }
    }
  }
  
  h_2D_xVtSin2_MC->Write();
  h_2D_xVtCos2_MC->Write();
  h_2D_xVtSin2_RP->Write();
  h_2D_xVtCos2_RP->Write();
  for(int q{0}; q<nQ2bins; q++) h_2D_xVtDiff_RP[q]->Write();
  
  fOutDiff->Close();

  return;
}
