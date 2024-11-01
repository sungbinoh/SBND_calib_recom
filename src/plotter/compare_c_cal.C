#include "canvas_margin.h"

void compare_c_cal(){
  setTDRStyle();
  gStyle->SetLineWidth(3);
  double ccal[6] = {1.909, 1.925, 1.929, 1.926, 1.930, 1.929};
  double ccal_err[6] = {0.004, 0.003, 0.003, 0.003, 0.003, 0.003};

  double y[6] = {1., 2., 3., 4., 5., 6.};
  double y_err[6] = {0., 0., 0., 0., 0., 0.};

  TCanvas *c = new TCanvas("", "", 600, 800);
  TPad *pad1 = new TPad("", "", 0.25, 0, 0.90, 0.95);
  pad1 -> SetLeftMargin( 0.15 );
  pad1 -> SetRightMargin( 0.03 );
  pad1 -> Draw();
  pad1 -> cd();
  TH1F * pad1_template = new TH1F("", "", 1, 1.899, 1.941);
  TString label[3] = {"Plane0", "Plane1", "Plane2"};
  pad1_template -> GetYaxis() -> SetRangeUser(0., 7.);
  pad1_template -> GetYaxis() -> SetLabelSize(0);
  pad1_template -> GetYaxis() -> SetNdivisions(0000);
  pad1_template -> GetXaxis() -> SetTitle("Collectoin Plane C_{cal.} [ #times 10^{-2} ]");
  pad1_template -> GetXaxis() -> CenterTitle(true);
  pad1_template -> GetXaxis() -> SetNdivisions(2105);
  pad1_template -> GetXaxis() -> SetLabelSize(0.05);
  pad1_template -> SetStats(0);
  pad1_template -> Draw();

  TGraphErrors *err_gr = new TGraphErrors(6, ccal, y, ccal_err, y_err);
  err_gr -> SetLineColor(kBlue);
  err_gr -> Draw("ezpsame");
  
  gPad->RedrawAxis();

  c -> cd();
  TLatex latex0, latex1, latex2, latex3, latex4, latex5, latex6;
  double text_size = 0.035;
  latex0.SetNDC();
  latex1.SetNDC();
  latex2.SetNDC();
  latex3.SetNDC();
  latex4.SetNDC();
  latex5.SetNDC();
  latex6.SetNDC();
  latex0.SetTextSize(text_size);
  latex1.SetTextSize(text_size);
  latex2.SetTextSize(text_size);
  latex3.SetTextSize(text_size);  
  latex4.SetTextSize(text_size);
  latex5.SetTextSize(text_size);
  latex6.SetTextSize(text_size);
  latex0.SetTextAlign(31);
  latex1.SetTextAlign(31);
  latex2.SetTextAlign(31);
  latex3.SetTextAlign(31);
  latex4.SetTextAlign(31);
  latex5.SetTextAlign(31);
  latex6.SetTextAlign(31);

  latex0.DrawLatex(0.30, 0.226, "Run 14480");
  latex1.DrawLatex(0.30, 0.337, "Run 14608");
  latex2.DrawLatex(0.30, 0.448, "Run 14685");
  latex3.DrawLatex(0.30, 0.559, "Run 14724");
  latex4.DrawLatex(0.30, 0.670, "Run 14833");
  latex5.DrawLatex(0.30, 0.781, "Run 14860");
  
  /*
  latex0.DrawLatex(0.08, 0.74, "Previous BDT");

  latex1.DrawLatex(0.10, 0.60, "No GENIE");
  latex2.DrawLatex(0.05, 0.56, "for Ext. and Signal");

  latex3.DrawLatex(0.10, 0.45, "No GENIE");
  latex4.DrawLatex(0.11, 0.42, "for Ext.");

  latex5.DrawLatex(0.10, 0.29, "GENIE for");
  latex6.DrawLatex(0.12, 0.25, "all MC");
  */
  
  TLatex latex_SBND;
  latex_SBND.SetNDC();
  latex_SBND.SetTextSize(0.035);
  latex_SBND.SetTextAlign(31);
  latex_SBND.DrawLatex(0.885, 0.91, "#font[62]{SBND} #font[42]{#it{#scale[0.8]{Preliminary}}}");

  TString pdfname;
  TString WORKING_DIR = getenv("PLOT_PATH");
  WORKING_DIR = WORKING_DIR + "/comparison/";
  pdfname = WORKING_DIR + "ccal_comparison.pdf";
  c -> SaveAs(pdfname);

}
