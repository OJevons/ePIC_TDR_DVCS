// Running macro for ePIC DVCS analysis
// ----> VARIANT FOR DIS PHYSICS BACKGROUND
#include "./ePIC_DISBkg_PODIO.cxx"

const float fMass_proton{0.938272};
const Float_t fMass_electron{0.000511};

void run_ePIC_DISBkg(TString camp="Camp", TString energy="10x100", TString sett="test", TString comment="X"){

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
  
  TString sInFileList = "./filelists/inputFileList_ePIC_"+camp+"_"+energy+"_"+sett+".list";
  sInFileList="./filelists/inputFileList_ePIC_DIS_26.02.0_10x100_minQ2=1.list";
  objDIS->setInFileList(sInFileList);

  TString sOutFileName;
  if(comment=="X") sOutFileName = "$EIC_WORK_DIR/DVCS_Ana/rootfiles/ePIC_DVCS_"+camp+"_"+energy+"_"+sett+".root";
  else sOutFileName = "$EIC_WORK_DIR/DVCS_Ana/rootfiles/ePIC_DVCS_"+camp+"_"+energy+"_"+comment+".root";
  sOutFileName = "$EIC_WORK_DIR/DVCS_Analysis/RootFiles/ePIC_DIS_26.02.0_10x100_minQ2=1_NewCuts-.root";

  objDIS->setOutFile(sOutFileName);

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

