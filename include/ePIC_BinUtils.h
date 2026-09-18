//---------------------------------------------------------------------------------------
//
// Utility functions for splitting histograms
// - By bins or by fixed number of entries
// Author: O. Jevons, 18/09/26
//
//---------------------------------------------------------------------------------------

#include "TMath.h"

//-----------------------------------------------------------------------------------------------------------------------------
// FUNCTION DEFINITIONS
//
// NOTE: Templating applied for brevity
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
      cout<<"["<<edgelo<<" - "<<edgehi<<"] ("<<hist->Integral(binlo,bin)<<")\n";
      
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
