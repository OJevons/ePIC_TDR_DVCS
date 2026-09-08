//---------------------------------------------------------------------------------------
//
// Utility functions; calculation of kinematic quantities for exclusive ePIC analyses
//
// Author: O. Jevons, 27/02/25
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

// Extra utilities - calculate 10ths of histogram contents
template<typename h>
void splitByTenths(const h& hist){
  int iTenthBin{1};

  // Count over all bins in histogram
  for(int iHistBin{1}; iHistBin<hist->GetNbinsX(); iHistBin++){
    // Break out of loop after 90% mark (last 10% up to last bin)
    if(iTenthBin == 10) break;

    // Else calculate fraction of events up to ith bin
    float frac = (float)hist->Integral(1,iHistBin)/hist->Integral();
    
    // First time fraction of events is greater than the n*10% mark, print previous bin and increment n
    if(frac > (float)iTenthBin/10){
      cout<<(int)iTenthBin*10<<"% \t bin no. "<<iHistBin-1<<"\t bin centre "<<hist->GetBinCenter(iHistBin-1)<<endl;
      iTenthBin++;
    }
  }

  return;
}

// Extra utilities - calculate fractions of histogram contents
// Use bins for internal bookkeeping but prints edges for analysis use
// 1. Print to screen
template<typename h>
void splitByBins_print(const h& hist, int nbins){
  int iFracBin{1}; // index of new split bin
  int binlo{1};    // low split bin - only used for Integral()
  float edgelo{0.}, edgehi{0.}; // low and high edges of split bins
  float thres{0.}; // threshold of subdivided integral
  
  // Exit early if trying to divde by more than the original number of bins
  if(nbins >= hist->GetNbinsX()){
    cout<<"Trying to divide into too many bins. Original histogram has "<<hist->GetNbinsX()<<" bins."<<endl;
    return;
  }

  // If valid, start subdivision process
  // Start from first bin in histogram
  edgelo = hist->GetBinLowEdge(1);
  thres = hist->Integral()/nbins;
  
  cout<<"Splitting "<<hist->GetName()<<" into "<<nbins<<" bins:     ";
  // Count over all bins in histogram
  for(int bin{1}; bin<hist->GetNbinsX(); bin++){
    // Break out of loop if in last new bin
    if(iFracBin == nbins){
      edgehi = hist->GetBinLowEdge(hist->GetNbinsX()+1);
      cout<<"["<<edgelo<<" - "<<edgehi<<"]"<<endl;
      break;
    }

    // Else calculate fraction of events up to ith bin
    float frac{0.};
    frac = hist->Integral(binlo,bin);

    // If above threshold, set rebin high edge, print bin limits, and set new rebin low edge
    if(frac > thres){
      edgehi = hist->GetBinLowEdge(bin+1);
      cout<<"["<<edgelo<<" - "<<edgehi<<"]     "<<frac<<"/"<<thres<<endl;

      binlo = bin+1;
      edgelo = edgehi;
      iFracBin++;
    }
  }

  return;
}

// 2. Save in array
template<typename h>
void splitByBins_array(const h& hist, int nbins, vector<float>& edges){
  int iFracBin{1}; // index of new split bin
  float edgelo{0.}, edgehi{0.}; // low and high edges of split bins
  int binlo{1}; // low and split bin - only used for Integral()
  float thres{0.}; // threshold of subdivided integral

  edges.clear();

  // Exit early if trying to divde by more than the original number of bins
  if(nbins >= hist->GetNbinsX()){
    cout<<"Trying to divide into too many bins. Original histogram has "<<hist->GetNbinsX()<<" bins."<<endl;
    return;
  }

  // If valid, start subdivision process
  // Start from first bin in histogram
  edgelo = hist->GetBinLowEdge(1);
  edges.push_back(edgelo);
  thres = hist->Integral()/nbins;
  
  // Count over all bins in histogram
  for(int bin{1}; bin<hist->GetNbinsX(); bin++){
    // Break out of loop if in last new bin
    if(iFracBin == nbins){
      edgehi = hist->GetBinLowEdge(hist->GetNbinsX()+1);
      edges.push_back(edgehi);
      break;
    }

    // Else calculate fraction of events up to ith bin
    float frac{0.};
    frac = hist->Integral(binlo,bin);

    // If above threshold, set rebin high edge, print bin limits, and set new rebin low edge
    if(frac > thres){
      edgehi = hist->GetBinLowEdge(bin+1);
      edges.push_back(edgehi);

      binlo = bin+1;
      iFracBin++;
    }
  }

  return;
}


// Extra utilities - split histogram into bins of a given no. of entries
// Use bins for internal bookkeeping, but print edges for analysis use
// 1. Print to screen
template<typename h>
void splitByEntries_print(const h& hist, int nentries){
  // Exit early if second parameter is zero or negative
  if(nentries <= 0){
    cout<<"Enter a valid number of entries (positive integer)"<<endl;
    return;
  }
  // Exit early if 'nentries' is greater than the no. of events in the histogram
  if(nentries > hist->GetEntries()){
    cout<<"Attempted no. of entries too high. Original histogram contains "<<hist->GetEntries()<<" events"<<endl;
    return;
  }

  cout<<"Splitting histogram "<<hist->GetName()<<" ("<<hist->Integral()<<" entries) into bins of roughly "<<nentries<<" (minimum "<<nentries*0.95<<")"<<endl;

  int inewbin{1};
  int binlo{1};
  float edgelo{0.}, edgehi{0.};
  
  // Start from lowest bin edge in histogram
  edgelo = hist->GetBinLowEdge(1);

  for(int bin{1}; bin<hist->GetNbinsX(); bin++){
    // Last bin - just run to end of histogram (& break out of loop)
    if(inewbin == (int)hist->Integral()/nentries){
      edgehi = hist->GetBinLowEdge(hist->GetNbinsX()+1);
      cout<<"["<<edgelo<<" - "<<edgehi<<"] ("<<hist->Integral(binlo,hist->GetNbinsX())<<")"<<endl;
      break;
    }

    // If count above threshold, print new bin edges
     if(hist->Integral(binlo,bin) > nentries*0.95){
      edgehi = hist->GetBinLowEdge(bin+1);
      cout<<"["<<edgelo<<" - "<<edgehi<<"] ("<<hist->Integral(binlo,bin)<<")     ";
      
      binlo = bin+1;
      edgelo = edgehi;
      inewbin++;
    }
  }

  return;
}

// 2. Save to array
template<typename h>
void splitByEntries_array(const h& hist, int nentries, vector<float>& edges){
  // Exit early if second parameter is zero or negative
  if(nentries <= 0){
    cout<<"Enter a valid number of entries (positive integer)"<<endl;
    return;
  }
  
  // Exit early if 'nentries' is greater than the no. of events in the histogram
  if(nentries > hist->GetEntries()){
    cout<<"Attempted no. of entries too high. Original histogram contains "<<hist->GetEntries()<<" events"<<endl;
    return;
  }

  int count{0};
  int inewbin{1};
  int binlo{1};
  float edgelo{0.}, edgehi{0.};
  edges.clear();
  
  // Start from lowest bin edge in histogram
  edgelo = hist->GetBinLowEdge(1);
  edges.push_back(edgelo);

  for(int bin{1}; bin<hist->GetNbinsX(); bin++){
    // Last bin - just run to end of histogram (& break out of loop)
    if(inewbin == (int)hist->Integral()/nentries){
      edgehi = hist->GetBinLowEdge(hist->GetNbinsX()+1);
      edges.push_back(edgehi);
      break;
    }

    // Else - calculate integral of entries from new bin low edge to current bin
    // If count above threshold, print new bin edges and reset count
    if(hist->Integral(binlo,bin) > nentries*0.95){
      edgehi = hist->GetBinLowEdge(bin+1);
      edges.push_back(edgehi);

      binlo = bin+1;
      inewbin++;
    }
  }

  return;
}

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


// Calculate Mandelstam t - BABE method using tRECO conventions
// Uses incoming proton BEam and scattered BAryon 4-vectors
// Another way of saying t = -(p' - p)^2
//--------------------------------------------------------------
// NEEDS: p, p'
//--------------------------------------------------------------
template <typename V>
Double_t calcT_BABE(const V& be, const V& ba){
  double t = (ba - be).M2();
  
  return TMath::Abs(t);
}

// Calculate Mandelstam t - eX method using tRECO conventions
// Uses difference between the beam and scattered electron and all of the (non-scattered) final state
// e.g. for DVCS, X is a photon; for electroproduction, it is the species of interest; etc...
//--------------------------------------------------------------
// NEEDS: e, e', X   [q, X]
//--------------------------------------------------------------
// 1. Separate vectors for beam and scattered electrons
template <typename V>
Double_t calcT_eX(const V& e, const V& ep, const V& X){
  double t = (e - ep - X).M2();
  return TMath::Abs(t);
}
// 2. Giving virtual photon vector directly
template <typename V>
Double_t calcT_eX(const V& q, const V& X){
  double t = (q - X).M2();
  return TMath::Abs(t);
}

// Calculate Mandelstam t - eXBA method using tRECO conventions
// Include final state baryon information into eX method
//--------------------------------------------------------------
// NEEDS: e, e', p', X 
// CANNOT JUST GIVE q-VECTOR; NEED INFO. FROM e'
//--------------------------------------------------------------
template <typename V>
Double_t calcT_eXBA(const V& e, const V& ep, const V& pp, const V& X){
  // Extract info. from vectors - e' energy and theta
  double E_ep = ep.E();
  double theta_ep = ep.Theta();

  // Intermediate calculations - q vector and HFS sigma
  P3EVector q(e.X()-ep.X(), e.Y()-ep.Y(), e.Z()-ep.Z(), e.E()-ep.E());
  double sigma_h = (pp+X).E() - (pp+X).Z();
  double sigterm = sigma_h/2;
  double eterm = (E_ep*(1+TMath::Cos(theta_ep)))/2;
  
  P3EVector pcorr(q.X(), q.Y(), -sigterm-eterm, sigterm-eterm);
  
  double t = (pcorr - X).M2();
  return TMath::Abs(t);
}

// Calculate Mandelstam t - eXBE method using tRECO conventions
// Include initial beam proton information into eX method
//--------------------------------------------------------------
// NEEDS: e, p, ep, pp, X    [q, p, pp, X]
//--------------------------------------------------------------
// 1. Separate vectors for beam and scattered electrons
template<typename V>
Double_t calcT_eXBE(const V& e, const V& p, const V& ep, const V& pp, const V& X){
  // Calculate 'missing' momentum, ignoring scattered baryon vector
  P3EVector p4miss((e+p-ep-X).X(),(e+p-ep-X).Y(),(e+p-ep-X).Z(),(e+p-ep-X).E());
    
  // Define corrected momentum vector using missing momentum and scattered baryon mass
  Float_t pmiss_mag = p4miss.Vect().R();
  Float_t pcorr_mag = TMath::Sqrt(TMath::Power(pmiss_mag,2) + TMath::Power(pp.M(),2));
  P3EVector pcorr(p4miss.Vect().X(), p4miss.Vect().Y(), p4miss.Vect().Z(), pcorr_mag);

  double t = (pcorr-p).M2();
  return TMath::Abs(t);
}
// 2. Giving virtual photon vector directly
template<typename V>
Double_t calcT_eXBE(const V& p, const V& q, const V& pp, const V& X){
  // Calculate 'missing' momentum, ignoring scattered baryon vector
  P3EVector p4miss((p+q-X).X(),(p+q-X).Y(),(p+q-X).Z(),(p+q-X).E());
    
  // Define corrected momentum vector using missing momentum and scattered baryon mass
  Float_t pmiss_mag = p4miss.Vect().R();
  Float_t pcorr_mag = TMath::Sqrt(TMath::Power(pmiss_mag,2) + TMath::Power(pp.M(),2));
  P3EVector pcorr(p4miss.Vect().X(), p4miss.Vect().Y(), p4miss.Vect().Z(), pcorr_mag);

  double t = (pcorr-p).M2();
  return TMath::Abs(t);
}
// 3. Using scalar mass of scattered baryon (with e/e')
template<typename V>
Double_t calcT_eXBE(const V& e, const V& p, const V& ep, const Float_t& mb, const V& X){
  // Calculate 'missing' momentum, ignoring scattered baryon vector
  P3EVector p4miss((e+p-ep-X).X(),(e+p-ep-X).Y(),(e+p-ep-X).Z(),(e+p-ep-X).E());
    
  // Define corrected momentum vector using missing momentum and scattered baryon mass
  Float_t pmiss_mag = p4miss.Vect().R();
  Float_t pcorr_mag = TMath::Sqrt(TMath::Power(pmiss_mag,2) + TMath::Power(mb,2));
  P3EVector pcorr(p4miss.Vect().X(), p4miss.Vect().Y(), p4miss.Vect().Z(), pcorr_mag);

  double t = (pcorr-p).M2();
  return TMath::Abs(t);
}
// 4. Using scalar mass of scattered baryon (with q)
template<typename V>
Double_t calcT_eXBE(const V& p, const V& q, const Float_t& mb, const V& X){
  // Calculate 'missing' momentum, ignoring scattered baryon vector
  P3EVector p4miss((p+q-X).X(),(p+q-X).Y(),(p+q-X).Z(),(p+q-X).E());
    
  // Define corrected momentum vector using missing momentum and scattered baryon mass
  Float_t pmiss_mag = p4miss.Vect().R();
  Float_t pcorr_mag = TMath::Sqrt(TMath::Power(pmiss_mag,2) + TMath::Power(mb,2));
  P3EVector pcorr(p4miss.Vect().X(), p4miss.Vect().Y(), p4miss.Vect().Z(), pcorr_mag);

  double t = (pcorr-p).M2();
  return TMath::Abs(t);
}

// Calculate Mandelstam t - eBABE method using tRECO conventions
// Include electron (beam and scattered) information into BABE method
//--------------------------------------------------------------
// NEEDS: e, p, ep, pp    [q, p, pp]
//--------------------------------------------------------------
// 1. Separate vectors for beam and scattered electrons
template<typename V>
Double_t calcT_eBABE(const V& e, const V& p, const V& ep, const V& pp){
  // Calculate needed vectors
  P3EVector q(e.X()-ep.X(), e.Y()-ep.Y(), e.Z()-ep.Z(), e.E()-ep.E());
  P3EVector pcorr(-q.X(), -q.Y(), pp.Z(), pp.E());

  double t = (pcorr - p).M2();
  return TMath::Abs(t);
}
// 2. Giving virtual photon vector directly
template <typename V>
Double_t calcT_eBABE(const V& p, const V& q, const V& pp){
  P3EVector pcorr(-q.X(), -q.Y(), pp.Z(), pp.E());

  double t = (pcorr - p).M2();
  return TMath::Abs(t);
}

// Calculate Mandelstam t - XBABE method using tRECO conventions
// Includes further final state information into BABE method
//--------------------------------------------------------------
// NEEDS: p, pp, X
// CANNOT PROVIDE HFS BY ITSELF
//--------------------------------------------------------------
template<typename V>
Double_t calcT_XBABE(const V& p, const V& pp, const V& X){
  P3EVector pcorr(-X.X(), -X.Y(), pp.Z(), pp.E());

  double t = (pcorr - p).M2();
  return TMath::Abs(t);
}

// Calculate Mandelstam t - eXBABE method using tRECO conventions
// Uses full event information
//--------------------------------------------------------------
// NEEDS: e, p, ep, pp, X    [q, p, pp, X]
//--------------------------------------------------------------
// 1. Separate vectors for beam and scattered electrons
template<typename V>
Double_t calcT_eXBABE(const V& e, const V& p, const V& ep, const V& pp, const V& X){
  // Calculate 'missing' momentum, ignoring scattered baryon vector
  P3EVector p4miss((e+p-ep-X).X(),(e+p-ep-X).Y(),(e+p-ep-X).Z(),(e+p-ep-X).E());
    
  // Define corrected momentum vector using missing momentum and scattered baryon mass
  Float_t pmiss_mag = p4miss.Vect().R();
  ROOT::Math::Polar3DVector pcorr_vect(pmiss_mag, pp.Theta(), pp.Phi());
  Float_t pcorr_mag = TMath::Sqrt(TMath::Power(pmiss_mag,2) + TMath::Power(pp.M(),2));
  
  P3EVector pcorr(pcorr_vect.X(), pcorr_vect.Y(), pcorr_vect.Z(), pcorr_mag);

  double t = (pcorr-p).M2();
  return TMath::Abs(t);
}
// 2. Giving virtual photon vector directly
template<typename V>
Double_t calcT_eXBABE(const V& p, const V& q, const V& pp, const V& X){
  // Calculate 'missing' momentum, ignoring scattered baryon vector
  P3EVector p4miss((p+q-X).X(),(p+q-X).Y(),(p+q-X).Z(),(p+q-X).E());
    
  // Define corrected momentum vector using missing momentum and scattered baryon mass
  Float_t pmiss_mag = p4miss.Vect().R();
  ROOT::Math::Polar3DVector pcorr_vect(pmiss_mag, pp.Theta(), pp.Phi());
  Float_t pcorr_mag = TMath::Sqrt(TMath::Power(pmiss_mag,2) + TMath::Power(pp.M(),2));
  
  P3EVector pcorr(pcorr_vect.X(), pcorr_vect.Y(), pcorr_vect.Z(), pcorr_mag);

  double t = (pcorr-p).M2();
  return TMath::Abs(t);
}

// Calculate Mandelstam t - eHe method from ECCE EDT paper
// A. Bylinkin et al., Nucl. Instrum. Meth. A 1052, 168238 (2023); eq. 9
// Needs full event information - give vectors in lab frame
template<typename V>
Double_t calcT_eHe(const V& e, const V& p, const V& ep, const V& pp, const V& X){
  // Extract intermediate quantities from 4-vectors
  double M = p.M();
  double nu = e.E() - ep.E();
  double Q2 = (e-ep).M2();

  // Calculate cos(theta between virtual and real photon)
  P3EVector q((e-ep).X(), (e-ep).Y(), (e-ep).Z(), (e-ep).E());

  double cTheta = TMath::Cos( (q-X).Theta() );
  
  double cosTerm = TMath::Sqrt((nu*nu)+Q2)*cTheta;
  double num = M*Q2 + (2*M*nu)*(nu-cosTerm);
  double den = M + nu - cosTerm;

  double t = num/den;
  return TMath::Abs(t);
}

// Calculate Mandelstam t - "Method L"
// Described in presentation from Jihee Kim 
// ePIC Exclusive, Diffractive and Tagging PWG meeting; 10th November 2025
template<typename V>
Double_t calcT_MethodL(const V& e, const V& p, const V& ep, const Float_t& mb, const V& X){
  // Calculate missing hadron vector from rest of final state
  P3EVector pmiss = p - (ep + X - e);
  
  // Express scattered hadron in terms of lightcone variables
  Double_t pplus = pmiss.E() + pmiss.Pz();
  Double_t pT2 = TMath::Power(pmiss.Px(), 2) + TMath::Power(pmiss.Py(), 2);
  
  // Correct missing momentum using known mass of scattered hadron
  Double_t num = TMath::Power(mb,2) + pT2;
  Double_t pminus = num/pplus;
  P3EVector pcorr(pmiss.Px(), pmiss.Py(), (pplus-pminus)/2, (pplus+pminus)/2);
  
  //double t = (pcorr-p).M2();
  double t = -(p-pcorr).M2();
  return TMath::Abs(t);
}

// Calculate missing kinematics (mass/energy/momentum)
// 3-body final state: ab->cdf
// Missing momentum
template <typename V>
Double_t calcPMiss_3Body(const V& a, const V& b, const V& c, const V& d, const V& f){ 
  return (a+b-c-d-f).P(); 
}
// Missing transverse momentum
template <typename V>
Double_t calcPtMiss_3Body(const V& a, const V& b, const V& c, const V& d, const V& f){
  return (a+b-c-d-f).Pt(); 
}
// Missing energy
template <typename V>
Double_t calcEMiss_3Body(const V& a, const V& b, const V& c, const V& d, const V& f){
  return (a+b-c-d-f).E(); 
}
// Missing mass (squared)
template <typename V>
Double_t calcM2Miss_3Body(const V& a, const V& b, const V& c, const V& d, const V& f){
  Float_t fEMiss = (a+b-c-d-f).E();
  Float_t fPMiss = (a+b-c-d-f).P();

  Float_t fM2Miss = TMath::Power(fEMiss,2) - TMath::Power(fPMiss,2);
  return fM2Miss;
}

// 2-body final state: ab->cd
// Missing momentum
template <typename V>
Double_t calcPMiss_2Body(const V& a, const V& b, const V& c, const V& d){
  return (a+b-c-d).P();
}
// Missing transverse momentum
template <typename V>
Double_t calcPtMiss_2Body(const V& a, const V& b, const V& c, const V& d){
  return (a+b-c-d).Pt();
}
// Missing energy
template <typename V>
Double_t calcEMiss_2Body(const V& a, const V& b, const V& c, const V& d){
  return (a+b-c-d).E();
}
// Missing mass (squared)
template <typename V>
Double_t calcM2Miss_2Body(const V& a, const V& b, const V& c, const V& d){
  Float_t fEMiss = (a+b-c-d).E();
  Float_t fPMiss = (a+b-c-d).P();

  Float_t fM2Miss = TMath::Power(fEMiss,2) - TMath::Power(fPMiss,2);
  return fM2Miss;
}

// Calculate Trento Phi
// Calculate angle between hadronic and leptonic planes (Trento phi)
// Using planes defined by [k, q] and [q, g']
// Source: Bachetta, A. et al; Phys. Rev. D (2004); eq. 16
template <typename V>
Double_t calcTrentoPhi_qg(const V& k, const V& p, const V& kprime, const V& gprime){  
  // Before calculating angle, boost into gamma*-p rest frame
  // Calculate q in lab frame
  V q = (k-kprime);
  // Boost vector
  MomVector vTgtRest = (p+q).BoostToCM();

  V kB = boost(k,vTgtRest);
  V kpB = boost(kprime,vTgtRest);
  V gpB = boost(gprime,vTgtRest);

  MomVector k3 = kB.Vect();
  MomVector kp3 = kpB.Vect();
  MomVector gp3 = gpB.Vect();
  MomVector qhat3 = (k3-kp3).Unit();

  // Define leptonic plane using virtual photon and scattered electron
  MomVector lNorm = qhat3.Cross(kp3);
  lNorm /= lNorm.R();
  // Define hadronic plane using q vector and scattered photon
  MomVector hNorm = qhat3.Cross(gp3);
  hNorm /= hNorm.R();

  // Angle() function just returns magnitude of angle
  // If photon vector has a component parallel to the leptonic normal, should be positive. If opposite, negative.
  return TMath::Sign(1.,gp3.Dot(lNorm))*Angle(lNorm,hNorm);
}
