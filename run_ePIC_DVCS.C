// Running macro for ePIC DVCS analysis
//#include "ePIC_DVCS_TASK.cxx"
//#include "ePIC_DVCS_Res.cxx"
#include "ePIC_DVCS_PODIO.cxx"
//#include "ePIC_DVCS_NoDiff.cxx"

//#include "ePIC_DISBkg_PODIO.cxx"
//#include "ePIC_DVCS_BeamGas.cxx"
//#include "ePIC_DVCS_Test.cxx"

const float fMass_proton{0.938272};
const float fMass_electron{0.000511};

void run_ePIC_DVCS(TString camp="Camp", TString energy="10x100", TString sett="test", TString comment="X"){

  std::cout<<"----------------------------------------------------"<<std::endl;
  std::cout<<"                 ePIC DVCS Analysis                 "<<std::endl;
  std::cout<<"----------------------------------------------------"<<std::endl;
  std::cout<<std::endl;
  
  // Initialize DVCS analysis object
  std::cout<<"Settings:"<<std::endl;
  std::cout<<"\tCampaign - "<<camp<<std::endl;
  std::cout<<"\tBeam energy - "<<energy<<" GeV"<<std::endl;
  std::cout<<"\tBeam setting - "<<sett<<std::endl;
  
  ePIC_DVCS_TASK *objDVCS = new ePIC_DVCS_TASK(camp,energy,sett);
  
  TString sInFileList = "./filelists/inputFileList_ePIC_"+camp+"_"+energy+"_"+sett+".list";
  //sInFileList = "./filelists/inputFileList_test.list";
  objDVCS->setInFileList(sInFileList);

  TString sOutFileName;
  if(comment=="X") sOutFileName = "rootfiles/ePIC_DVCS_"+camp+"_"+energy+"_"+sett;
  else sOutFileName = "rootfiles/ePIC_DVCS_"+camp+"_"+energy+"_"+comment;
  
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

