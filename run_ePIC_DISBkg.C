// Running macro for ePIC DVCS analysis
// ----> VARIANT FOR DIS PHYSICS BACKGROUND
#include "./bkg_dis/ePIC_DISBkg_PODIO.cxx"

const float fMass_proton{0.938272};
const Float_t fMass_electron{0.000511};

void run_ePIC_DISBkg(TString camp="Camp", TString energy="9x130", TString sett="test", TString comment="X"){

  std::cout<<"----------------------------------------------------"<<std::endl;
  std::cout<<"            ePIC DVCS Analysis - DIS bkg            "<<std::endl;
  std::cout<<"----------------------------------------------------"<<std::endl;
  std::cout<<std::endl;
  
  // Initialize DVCS analysis object
  std::cout<<"Settings:"<<std::endl;
  std::cout<<"\tCampaign - "<<camp<<std::endl;
  std::cout<<"\tBeam energy - "<<energy<<" GeV"<<std::endl;
  std::cout<<"\tBeam setting - "<<sett<<std::endl;
  
  ePIC_DVCS_TASK *objDIS = new ePIC_DVCS_TASK(camp,energy,sett);
  
  TString sInFileList = "./bkg_dis/inputFileList_DIS_"+camp+"_"+energy+"_minQ2="+sett+".list";
  sInFileList="./bkg_dis/inputFileList_DIS_test.list";
  objDIS->setInFileList(sInFileList);

  // Output file name
  // DO NOT INCLUDE ".root" EXTENSION
  TString sOutFileName;
  if(comment=="X") sOutFileName = "./rootfiles/ePIC_DIS_"+camp+"_"+energy+"_minQ2="+sett;
  else sOutFileName = "./rootfiles/ePIC_DIS_"+camp+"_"+energy+"_"+comment;
  sOutFileName = "./bkg_dis/ePIC_DIS_test";

  objDIS->setOutFileName(sOutFileName);

  // Set DVCS cut values
  objDIS->setMomCutFactors(10.,10.);
  objDIS->setMin_Q2(1);         // GeV^2
  objDIS->setMax_pTmiss(0.5);   // GeV
  objDIS->setEmPzCuts(15.,25.); // GeV

  // Set other behaviours
  objDIS->setUsePID(kFALSE);
  objDIS->setUseExplicitMatch(kFALSE);
  objDIS->setUseEventBeams(kFALSE);

  objDIS->doAnalysis();
}

