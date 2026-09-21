// Running macro for ePIC DVCS analysis
//#include "class/ePIC_DVCS_TASK.cxx"
//#include "class/ePIC_DVCS_Res.cxx"
#include "class/ePIC_DVCS_PODIO.cxx"
//#include "class/ePIC_DVCS_BeamGas.cxx"

const float fMass_proton{0.938272};
const float fMass_electron{0.000511};

void run_ePIC_DVCS(TString camp="Camp", TString energy="10x100", TString sett="test", TString comment="X"){

  std::cout<<"----------------------------------------------------"<<std::endl;
  std::cout<<"                 ePIC DVCS Analysis                 "<<std::endl;
  std::cout<<"----------------------------------------------------"<<std::endl;
  std::cout<<std::endl;

  // Check pathing
  char* DVCS_ep_Path_char;
  DVCS_ep_Path_char = getenv("DVCS_ep");
  if (DVCS_ep_Path_char == nullptr) {
      cerr << "!!!!! ERROR !!!!! DVCS_ep environment variable not set !!!!! ERROR !!!!!" << endl;
      cerr << "!!!!! ERROR !!!!! Source the setup.sh (or .csh) script and rerun !!!!! ERROR !!!!!" << endl;
      exit(0);
  }
  TString DVCS_ep_Path(DVCS_ep_Path_char);
  
  // Initialize DVCS analysis object
  std::cout<<"Settings:"<<std::endl;
  std::cout<<"\tCampaign - "<<camp<<std::endl;
  std::cout<<"\tBeam energy - "<<energy<<" GeV"<<std::endl;
  std::cout<<"\tBeam setting - "<<sett<<std::endl;
  
  ePIC_DVCS_TASK *objDVCS = new ePIC_DVCS_TASK(camp,energy,sett);
  
  TString sInFileList = DVCS_ep_Path+"/filelists/inputFileList_ePIC_"+camp+"_"+energy+"_"+sett+".list";
  objDVCS->setInFileList(sInFileList);

  // Output file name
  // DO NOT INCLUDE ".root" EXTENSION
  TString sOutFileName;
  if(comment=="X") sOutFileName = DVCS_ep_Path+"rootfiles/ePIC_DVCS_"+camp+"_"+energy+"_"+sett;
  else sOutFileName = DVCS_ep_Path+"rootfiles/ePIC_DVCS_"+camp+"_"+energy+"_"+comment;
   
  objDVCS->setOutFileName(sOutFileName);
  
  // Set DVCS cut values
  objDVCS->setMomCutFactors(10.,10.);
  objDVCS->setMin_Q2(1);         // GeV^2
  objDVCS->setMax_pTmiss(0.5);   // GeV
  objDVCS->setEmPzCuts(15.,25.); // GeV

  // Set other behaviours
  objDVCS->setUsePID(kFALSE);
  objDVCS->setUseExplicitMatch(kFALSE);
  objDVCS->setUseEventBeams(kFALSE);

  objDVCS->doAnalysis();
}

