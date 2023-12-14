#ifndef canvas_margin_h
#define canvas_margin_h

#include "TStyle.h"

// fixOverlay: Redraws the axis

void fixOverlay() {
  gPad->RedrawAxis();
}

void setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);
  
// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  // tdrStyle->SetErrorMarker(20);
  //tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);
  
//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->cd();

}

void canvas_margin(TCanvas *c1){
  c1->SetTopMargin( 0.05 );
  c1->SetBottomMargin( 0.13 );
  c1->SetRightMargin( 0.05 );
  c1->SetLeftMargin( 0.16 );

}

void canvas_margin_limit(TCanvas *c1){
  c1->SetTopMargin( 0.05 );
  c1->SetBottomMargin( 0.13 );
  c1->SetRightMargin( 0.05 );
  c1->SetLeftMargin( 0.20 );

}

void canvas_margin(TCanvas *c1, TPad *c1_up, TPad *c1_down){
  c1_up->SetTopMargin( 0.07 );
  c1_up->SetBottomMargin( 0.02 );
  c1_up->SetLeftMargin( 0.15 );
  c1_up->SetRightMargin( 0.03 );

  c1_down->SetTopMargin( 0.03 );
  c1_down->SetBottomMargin( 0.4 );
  c1_down->SetLeftMargin( 0.15 );
  c1_down->SetRightMargin( 0.03 );
  //c1_down->SetGridx();
  //c1_down->SetGridy();
  
  c1->SetTopMargin( 0.05 );
  c1->SetBottomMargin( 0.13 );
  c1->SetRightMargin( 0.05 );
  c1->SetLeftMargin( 0.16 );

}

void canvas_margin_twoplots(TCanvas *c1, TPad *c1_up, TPad *c1_down){
  c1_up->SetTopMargin( 0.07 );
  c1_up->SetBottomMargin( 0.017 );
  c1_up->SetLeftMargin( 0.15 );
  c1_up->SetRightMargin( 0.03 );

  c1_down->SetTopMargin( 0.017 );
  c1_down->SetBottomMargin( 0.15 );
  c1_down->SetLeftMargin( 0.15 );
  c1_down->SetRightMargin( 0.03 );

  c1->SetTopMargin( 0.05 );
  c1->SetBottomMargin( 0.13 );
  c1->SetRightMargin( 0.05 );
  c1->SetLeftMargin( 0.16 );

}

void hist_axis(TH1D *hist){

  hist->SetTitle("");

  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleOffset(1.10);

  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetXaxis()->SetTitleSize(0.05);

}

void hist_axis_limit(TH1D *hist){

  hist->SetTitle("");

  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleOffset(1.30);

  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetXaxis()->SetTitleSize(0.055);
  hist->GetXaxis()->SetTitleOffset(1.00);

}

void hist_axis(THStack *hist){

  hist->SetTitle("");

  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleOffset(1.10);

  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetXaxis()->SetTitleSize(0.05);

}

void hist_axis(TGraph *hist){

  hist->SetTitle("");

  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleOffset(1.10);

  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetXaxis()->SetTitleSize(0.05);

}

void hist_axis(TGraphAsymmErrors *hist){

  hist->SetTitle("");

  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleOffset(1.10);

  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetXaxis()->SetTitleSize(0.05);

}

void hist_axis(TH2D *hist){

  hist->SetTitle("");

  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleOffset(1.10);

  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetXaxis()->SetTitleSize(0.05);

}

void hist_axis(TH2F *hist){

  hist->SetTitle("");

  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleOffset(1.10);

  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetXaxis()->SetTitleSize(0.05);

}

void hist_axis(THStack *hist, TH1D *hist_compare){

  hist->SetTitle("");

  //==== top plot
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.07);
  hist->GetYaxis()->SetTitleOffset(1.02);
  //==== hide x-axis for top plot
  hist->GetXaxis()->SetLabelSize(0);

  //==== bottom plot
  hist_compare->SetTitle("");
  hist_compare->GetXaxis()->SetLabelSize(0.10);
  hist_compare->GetXaxis()->SetTitleSize(0.15);
  hist_compare->GetYaxis()->SetLabelSize(0.08);
  hist_compare->GetYaxis()->SetTitleSize(0.12);
  hist_compare->GetYaxis()->SetTitleOffset(0.5);
  hist_compare->SetFillColorAlpha(45,0.35);

}

void hist_axis(TH1D *hist, TH1D *hist_compare){

  hist->SetTitle("");

  //==== top plot
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.07);
  hist->GetYaxis()->SetTitleOffset(1.02);
  //==== hide x-axis for top plot
  hist->GetXaxis()->SetLabelSize(0);

  //==== bottom plot
  hist_compare->SetTitle("");

  hist_compare->GetXaxis()->SetLabelSize(0.13);
  hist_compare->GetXaxis()->SetTitleSize(0.18);

  hist_compare->GetYaxis()->SetLabelSize(0.13);
  hist_compare->GetYaxis()->SetTitleSize(0.16);
  hist_compare->GetYaxis()->SetTitleOffset(0.4);
  hist_compare->SetFillColorAlpha(45,0.35);

}

void hist_axis_twoplots(TH1D *hist, TH1D *hist_compare){

  hist->SetTitle("");

  //==== top plot
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.07);
  hist->GetYaxis()->SetTitleOffset(1.02);
  //==== hide x-axis for top plot
  hist->GetXaxis()->SetLabelSize(0);

  //==== bottom plot
  hist_compare->SetTitle("");

  hist_compare->GetXaxis()->SetLabelSize(0.03);
  hist_compare->GetXaxis()->SetTitleSize(0.06);
  hist_compare->GetYaxis()->SetTitleOffset(1.2);

  hist_compare->GetYaxis()->SetLabelSize(0.05);
  hist_compare->GetYaxis()->SetTitleSize(0.05);
  hist_compare->GetYaxis()->SetTitleOffset(0.4);
  hist_compare->SetFillColorAlpha(45,0.35);

}

void hist_axis(TGraph *hist, TGraph *hist_compare){

  hist->SetTitle("");

  //==== top plot
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.07);
  hist->GetYaxis()->SetTitleOffset(1.02);
  //==== hide x-axis for top plot
  hist->GetXaxis()->SetLabelSize(0);

  //==== bottom plot
  hist_compare->SetTitle("");
  hist_compare->GetXaxis()->SetLabelSize(0.10);
  hist_compare->GetXaxis()->SetTitleSize(0.15);
  hist_compare->GetYaxis()->SetLabelSize(0.08);
  hist_compare->GetYaxis()->SetTitleSize(0.12);
  hist_compare->GetYaxis()->SetTitleOffset(0.5);
  hist_compare->SetFillColorAlpha(45,0.35);

}

#endif
