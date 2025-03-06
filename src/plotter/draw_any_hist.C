#include "canvas_margin.h"
#include "mylib.h"
#include <iostream>

void draw_2D_hist(TString filename, TString histname, double x_min, double x_max, double y_min, double y_max, double rebin_x, double rebin_y, TString title_x, TString title_y, TString top_left, TString top_right){

  TString input_file_dir = getenv("OUTPUTROOT_PATH");
  TFile *f = new TFile(input_file_dir + "/" + filename);
  TH2D *hist_2D = (TH2D*)gDirectory -> Get(histname);
  hist_2D -> RebinX(rebin_x);
  hist_2D -> RebinY(rebin_y);

  TCanvas *c = new TCanvas("", "", 1600, 1200);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  gStyle->SetLineWidth(2);
  c -> SetRightMargin(0.15);

  double z_max = hist_2D -> GetMaximum();
  TH2D * template_h = new TH2D("", "", 1., x_min, x_max, 1., y_min, y_max);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(title_x);
  template_h -> GetXaxis() -> SetTitleSize(0.045);
  template_h -> GetXaxis() -> SetTitleOffset(1.0);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(title_y);
  template_h -> GetYaxis() -> SetTitleOffset(1.1);
  template_h -> GetYaxis() -> SetTitleSize(0.045);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetZaxis() -> SetRangeUser(0., z_max * 1.1);
  template_h -> GetZaxis() -> SetTitle("Entries");
  template_h -> GetZaxis() -> SetTitleSize(0.045);
  template_h -> GetZaxis() -> SetLabelSize(0.035);
  template_h -> Draw("colz");

  hist_2D -> Draw("colzsame");

  gPad->RedrawAxis();

  TLatex latex_top_left, latex_top_right;
  latex_top_left.SetNDC();
  latex_top_left.SetTextSize(0.03);
  latex_top_right.SetNDC();
  latex_top_right.SetTextAlign(31);
  latex_top_right.SetTextSize(0.03);
  latex_top_left.DrawLatex(0.16, 0.96, top_left);
  latex_top_right.DrawLatex(0.90, 0.96, top_right);

  TString output_plot_dir = getenv("PLOT_PATH");
  TString outfile_str = output_plot_dir + "/any_plot/" + histname + ".png";
  c -> SaveAs(outfile_str);

  c -> Close();
}

void draw_any_hist(){

  setTDRStyle();

  TString filename = "output_dqdx_uniform_data.root";

  TString planes[] = {"plane0", "plane1", "plane2"};
  for(int i = 0; i < 3; i++){
    TString this_plane = planes[i];
    draw_2D_hist(filename, "x_vs_y_" + this_plane + "_passing_cathode_trklen_400", -250., 250., -250., 250., 1., 1., "x [cm]", "y [cm]", "#font[62]{SBND} Run (14608, 14685, 14724, 14833, 14860)", this_plane);
    draw_2D_hist(filename, "x_vs_z_" + this_plane + "_passing_cathode_trklen_400", -250., 250., -50., 550., 1., 1., "x [cm]", "z [cm]", "#font[62]{SBND} Run (14608, 14685, 14724, 14833, 14860)", this_plane);
    draw_2D_hist(filename, "z_vs_y_" + this_plane + "_passing_cathode_trklen_400", -50., 550., -250., 250., 1., 1., "z [cm]", "y [cm]", "#font[62]{SBND} Run (14608, 14685, 14724, 14833, 14860)", this_plane);


  }
    
}
