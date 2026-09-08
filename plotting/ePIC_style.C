//
// ePIC Style
// version 0.1
//

#include <iostream>
#include "TROOT.h"

//===========================
TStyle* ePIC_style() 
{
  TStyle *ePIC_style = new TStyle("ePIC","ePIC style");

  // use plain black on white colors
  Int_t icol=0; // WHITE
  ePIC_style->SetFrameBorderMode(icol);
  ePIC_style->SetFrameFillColor(icol);
  ePIC_style->SetCanvasBorderMode(icol);
  ePIC_style->SetCanvasColor(icol);
  ePIC_style->SetPadBorderMode(icol);
  ePIC_style->SetPadColor(icol);
  ePIC_style->SetStatColor(icol);

  // set margin sizes
  ePIC_style->SetPadBottomMargin(0.12);
  ePIC_style->SetPadTopMargin(0.05);
  ePIC_style->SetPadLeftMargin(0.12);
  ePIC_style->SetPadRightMargin(0.05);

  // set title offsets (for axis label)
  ePIC_style->SetTitleXOffset(1.1);
  ePIC_style->SetTitleYOffset(1.1);
  ePIC_style->SetTitleOffset(1.1,"z"); // Set the offset for Z axis titles expliticly to avoid it being cut off

  // use large fonts
  //Int_t font=72; // Helvetica italics
  const Int_t g_font = 42; // Helvetica
  const Double_t g_tsize = 0.045;
  const Double_t g_tsize_label = 0.042;
  ePIC_style->SetTextFont(g_font);
  ePIC_style->SetTextSize(g_tsize);
  
  ePIC_style->SetLabelFont(g_font,"x");
  ePIC_style->SetTitleFont(g_font,"x");
  ePIC_style->SetLabelFont(g_font,"y");
  ePIC_style->SetTitleFont(g_font,"y");
  ePIC_style->SetLabelFont(g_font,"z");
  ePIC_style->SetTitleFont(g_font,"z");
  
  ePIC_style->SetLabelSize(g_tsize_label,"x");
  ePIC_style->SetTitleSize(g_tsize,"x");
  ePIC_style->SetLabelSize(g_tsize_label,"y");
  ePIC_style->SetTitleSize(g_tsize,"y");
  ePIC_style->SetLabelSize(g_tsize_label,"z");
  ePIC_style->SetTitleSize(g_tsize,"z");

  // use bold lines and markers
  ePIC_style->SetMarkerStyle(20);
  ePIC_style->SetMarkerSize(1.2);
  ePIC_style->SetHistLineWidth(2.);
  ePIC_style->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of error bar caps
  ePIC_style->SetEndErrorSize(0.);

  // statistics
  ePIC_style->SetOptStat(0);
  //ePIC_style->SetOptFit(1);
  ePIC_style->SetOptFit(0);
  ePIC_style->SetOptDate(0);
  ePIC_style->SetStatFontSize(g_tsize);

  // get rid of grid
  ePIC_style->SetPadGridX(0);
  ePIC_style->SetPadGridY(0);

  // legend modification
  ePIC_style->SetLegendBorderSize(0);
  ePIC_style->SetLegendFillColor(0);
  ePIC_style->SetLegendFont(g_font);

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
  std::cout << "ePIC_style: ROOT6 mode" << std::endl;
  ePIC_style->SetLegendTextSize(g_tsize);
  ePIC_style->SetPalette(kBird);
#else
  std::cout << "ePIC_style: ROOT5 mode" << std::endl;
  // color palette - manually define 'kBird' palette only available in ROOT 6
  Int_t alpha = 0;
  Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
  Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
  TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
#endif

  ePIC_style->SetNumberContours(80);

  return ePIC_style;

}

//===========================
void set_ePIC_style()
{
  static TStyle* epic_style = 0;
  std::cout << "ePIC style: Applying nominal settings." << std::endl ;
  if ( epic_style==0 ) epic_style = ePIC_style();
  gROOT->SetStyle("ePIC");
  gROOT->ForceStyle();
}
