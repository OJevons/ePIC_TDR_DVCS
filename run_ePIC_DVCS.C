// Running macro for ePIC DVCS analysis
#include "ePIC_DVCS_TDR.cxx"
const float fMass_proton{0.938272};
const Float_t fMass_electron{0.000511};

void run_ePIC_DVCS(TString camp="Camp", TString energy="10x100", TString sett="test", TString comment="X"){
  
  // Write to screen
  std::cout<<"----------------------------------------------------"<<std::endl;
  std::cout<<"                 ePIC DVCS Analysis                 "<<std::endl;
  std::cout<<"----------------------------------------------------"<<std::endl;
  std::cout<<std::endl;
  
  std::cout<<"Settings:"<<std::endl;
  std::cout<<"\tCampaign - "<<camp<<std::endl;
  std::cout<<"\tBeam energy - "<<energy<<" GeV"<<std::endl;
  std::cout<<"\tBeam setting - "<<sett<<std::endl;
  
  // Initialise DVCS task object
  ePIC_DVCS_TASK *objDVCS = new ePIC_DVCS_TASK(camp,energy,sett);
  
  // Set input file list and output file
  TString sInFileList = "./inputFileList_ePIC_"+camp+"_"+energy+"_"+sett+".list";
  objDVCS->setInFileList(sInFileList);
  TString sOutFileName = "./ePIC_DVCS_TDR_"+energy+".root";
  objDVCS->setOutFile(sOutFileName);
  
  // Set DVCS cut values
  objDVCS->setMomCutFactors(10.,10.);  // Expected maximum scattered particle momenta (as multiple of beam momenta)
  objDVCS->setMin_Q2(1);               // Minimum Q2 - GeV^2
  objDVCS->setMax_tRP(2);              // Maximum t for Roman Pot proton tracks - GeV^2
  objDVCS->setMax_M2miss(1);           // Maximum MM2 for fully-exclusive DVCS events - GeV^2
  
  // Set other behaviours
  objDVCS->setUsePID(kFALSE);          // Don't use ePIC PID response
  objDVCS->setUseExplicitMatch(kTRUE); // Use explicit MC matching for PID
  objDVCS->setUseEventBeams(kFALSE);   // Calculate beams as average of first file in list

  // Run analysis
  objDVCS->doAnalysis();
}

