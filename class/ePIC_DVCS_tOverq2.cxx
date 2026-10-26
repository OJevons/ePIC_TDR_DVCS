// ePIC DVCS analysis class definition
#include "../include/preLoadLib.hh"

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
    // Photon-proton coplanarity cut
    float phigp = calcAzimGP(beame, beamp, scate[0], scatp[0], scatg[0]); 
    float phidiff = TMath::Abs(phigp)-TMath::Pi();
    if(TMath::Abs(phidiff) > 0.18) passCuts = false;
    
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

// Calculation of difference in azimuthal angle of real photon and scattered proton, in the qp CoM frame
Double_t ePIC_DVCS_TASK::calcAzimGP(P3EVector k, P3EVector p, P3EVector kprime, P3EVector pprime, P3EVector gprime){
  // Before calculating angle, boost into gamma*-p rest frame
  // Calculate q in lab frame
  P3EVector q = (k-kprime);
  // Boost vector
  MomVector vTgtRest = (p+q).BoostToCM();

  gprime = boost(gprime,vTgtRest);
  pprime = boost(pprime,vTgtRest);

  // Extract azimuthal angles and calculate difference
  double phi_gp = gprime.Phi();
  double phi_pp = pprime.Phi();

  for(float phi : {phi_gp,phi_pp}){ if(phi < 0) phi += 2*TMath::Pi(); }

  return phi_gp-phi_pp;
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
  
  //|t|/Q2 vs x
  TH2D* h_tq2Vx_MC = new TH2D("tq2vx_mc",";x_{B} (MC);|t|/Q^{2} (MC)",1e5,0.,1.,1e5,0.,2.);
  TH2D* h_tq2Vx_RP = new TH2D("tq2vx_rp",";x_{B} (Rec);|t|/Q^{2} (Rec)",1e5,0.,1.,1e5,0.,2.);
  
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
      //const auto& ZDCparts = event.get<edm4eic::ReconstructedParticleCollection>("ReconstructedFarForwardZDCNeutrals");
      for(const auto& zdcreco : ZDCparts){
	// Check for associated clusters
	if(zdcreco.getClusters()) kOtherFFTrack_Rec = true;
      }
      
      //cout<<"[DEBUG] Vetoes found"<<endl;

      //---------------------------------------------------------
      // Fill histograms
      //---------------------------------------------------------
      
      bool kAllVetoes_MC = !kNonElecNeg_MC && !kOtherFFTrack_MC && !kBarrelPos_MC;
      bool kAllVetoes_Rec = !kNonElecNeg_Rec && !kOtherFFTrack_Rec && !kBarrelPos_Rec;
      
      if(applyCuts_All(beame4, beamp4, scate4_gen, scatp4_gen, scatg4_gen, "all") && kAllVetoes_MC){
	Float_t t_gen = calcT_BABE(beamp4,scatp4_gen[0]);
	Float_t q2_gen = calcQ2_Elec(beame4, scate4_gen[0]);
	Float_t xB_gen = calcX_Elec(beame4, beamp4, scate4_gen[0]);
	
	h_tq2Vx_MC->Fill(xB_gen,t_gen/q2_gen);
      }
      
      if(applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rec, scatg4_rec, "B0") && scatp4_rom.size() == 0  && kAllVetoes_Rec){
	Float_t t_rec = calcT_BABE(beamp4,scatp4_rec[0]);
	Float_t q2_rec = calcQ2_Elec(beame4, scate4_rec[0]);
	Float_t xB_rec = calcX_Elec(beame4, beamp4, scate4_rec[0]);
	
	h_tq2Vx_RP->Fill(xB_rec,t_rec/q2_rec);
      }
      
      
      if(applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rom, scatg4_rec, "RP") && scatp4_rec.size() == 0 && kAllVetoes_Rec
	 && scatp4_gen.size() > 0 && scate4_aso.size() > 0 && scatg4_aso.size() > 0){
	Float_t t_rec = calcT_BABE(beamp4,scatp4_rom[0]);
	Float_t q2_rec = calcQ2_Elec(beame4, scate4_rec[0]);
	Float_t xB_rec = calcX_Elec(beame4, beamp4, scate4_rec[0]);
	
	h_tq2Vx_RP->Fill(xB_rec,t_rec/q2_rec);
      }
      
    }// END OF EVENT LOOP
    
    fileCounter++;
  } // END OF FILE LIST
  
  
  //------------------------------------------------------------
  // Write to output file
  //------------------------------------------------------------
  TString sOut = sOutFileName+"_kin.root";
  TFile* fOut = new TFile(sOut,"RECREATE");
  fOut->cd();
  h_tq2Vx_MC->Write();
  h_tq2Vx_RP->Write();
  fOut->Close();

  return;
}
