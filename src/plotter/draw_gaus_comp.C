#include "canvas_margin.h"
#include "mylib.h"
#include "LanGausFit.h"
#include "TRandom3.h"
#include "BetheBloch.h"
#include <iostream>

void draw(TString y_var, double x_min, double x_max, double y_min, double y_max, TString y_title){

  TString this_id = "id";
  TString gr_name = "rr_vs_" + y_var;

  TString input_file_dir = getenv("OUTPUTROOT_PATH");
  TFile *f = new TFile(input_file_dir + "/gaus_comp_fit.root");
  TGraphErrors *gr_plane0 = (TGraphErrors*)gDirectory -> Get("plane0_" + gr_name);
  TGraphErrors *gr_plane1 = (TGraphErrors*)gDirectory -> Get("plane1_" + gr_name);
  TGraphErrors *gr_plane2 = (TGraphErrors*)gDirectory -> Get("plane2_" + gr_name);

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D * template_h = new TH1D("", "", 1., x_min, x_max);
  template_h -> SetStats(0);
  template_h -> GetYaxis() -> SetRangeUser(y_min, y_max);
  template_h -> GetXaxis() -> SetTitle("Residual range [cm]");
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(y_title);
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> Draw();

  gr_plane0 -> SetLineColor(kBlue);
  gr_plane0 -> SetMarkerColor(kBlue);
  gr_plane0 -> Draw("epsame");

  gr_plane1 -> SetLineColor(kRed);
  gr_plane1 -> SetMarkerColor(kRed);
  gr_plane1 -> Draw("epsame");

  gr_plane2 -> SetLineColor(kGreen);
  gr_plane2 -> SetMarkerColor(kGreen);
  gr_plane2 -> Draw("epsame");

  TLegend *l = new TLegend(0.6, 0.6, 0.85, 0.85);
  l -> AddEntry(gr_plane0, "Plane 0", "lp");
  l -> AddEntry(gr_plane1, "Plane 1", "lp");
  l -> AddEntry(gr_plane2, "Plane 2", "lp");
  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_particle;
  latex_ProtoDUNE.SetNDC();
  latex_particle.SetNDC();
  latex_particle.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_particle.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND Data} Run 14860");
  latex_particle.DrawLatex(0.95, 0.96, "Cathode Passing Stopping Tracks");

  TString output_plot_dir = getenv("PLOT_PATH");
  TString outfile_str = output_plot_dir + "/comparison/plane_comparison_" + y_var + ".pdf";
  c -> SaveAs(outfile_str);
  c -> Close();
  
}

void draw_gaus_comp(){

  setTDRStyle();
  draw("MPV", 5., 100., 1., 4., "dE/dx MPV");
  draw("sigmaG", 5., 100., 0., 1., "#sigma_{G}");
  draw("sigmaL", 5., 100., 0., 0.2, "#sigma_{L}");

}
