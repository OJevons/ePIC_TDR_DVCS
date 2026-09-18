//---------------------------------------------------------------------------------------
//
// Utility functions; miscellaneous calculations
//
// Author: O. Jevons, 18/09/26
//
//---------------------------------------------------------------------------------------

#include "TMath.h"

// Aliases for common 3/4-vector types
using P3EVector=ROOT::Math::PxPyPzEVector;
using P3MVector=ROOT::Math::PxPyPzMVector;
using MomVector=ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>,ROOT::Math::DefaultCoordinateSystemTag>;
using ROOT::Math::XYZVector;

// Other ROOT::Math aliases
using ROOT::Math::VectorUtil::Angle;
using ROOT::Math::VectorUtil::boost;

//-----------------------------------------------------------------------------------------------------------------------------
// FUNCTION DEFINITIONS
//
// NOTE: Templating applied for brevity
//       4-vector functions are valid for any type which contains the operators E(), P(), Pt() and M2()
//       This includes TLorentzVector (legacy) and ALL variants of ROOT::Math::LorentzVector class
//
//-----------------------------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------------------------
// FUNCTION DEFINITIONS
//-----------------------------------------------------------------------------------------------------------------------------

// Calculate energy from momentum and mass
// 1. Using vector structures for momentum
// Works for ANY structure which contains Mag2() operator
template<typename P>
Double_t calcE(const P& mom, const Float_t& M){ 
  return TMath::Sqrt(mom.Mag2() + TMath::Power(M,2)); 
}
// 2. Using separate floats for momentum components
Double_t calcE(const Float_t& px, const Float_t& py, const Float_t& pz, const Float_t& M){ 
  return TMath::Sqrt(TMath::Power(px,2) + TMath::Power(py,2) + TMath::Power(pz,2) + TMath::Power(M,2)); 
}

