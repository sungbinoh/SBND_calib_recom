#include "canvas_margin.h"
#include "mylib.h"
#include <iostream>

void draw_bragg_2D(int run_num = 0, int rebin_x = 1, int rebin_y = 1){

  bool isdata = false;
  TString run_str = "";

  if(run_num != 0){
    isdata = true;
    run_str = TString::Format("%d", run_num);
  }

  setTDRStyle();

  TString input_file_dir = getenv("OUTPUTROOT_PATH");
  TFile *f = new TFile(input_file_dir + "/output_recom_" + run_str + ".root");
  TH2D *hist_2D = (TH2D*)gDirectory -> Get("rr_vs_dqdx_trklen_60cm_passing_cathode_coszx");
  hist_2D -> RebinX(rebin_x);
  hist_2D -> RebinY(rebin_y);
  
  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  c -> SetRightMargin(0.15);

  double z_max = hist_2D -> GetMaximum();
  TH1D * template_h = new TH1D("", "", 1., 0., 200.);
  template_h -> SetBinContent(1, -999999999999.);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("Residual range [cm]");
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("dQ/dx [ADC/cm]");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(-100., 3000.);
  template_h -> Draw("colz");
 
  hist_2D -> Draw("colzsame");

  TString particle_label_str = "Cathode Passing Stopping Tracks";
  TLatex latex_SBND, latex_method;
  latex_SBND.SetNDC();
  latex_method.SetNDC();
  latex_method.SetTextAlign(31);
  latex_SBND.SetTextSize(0.03);
  latex_method.SetTextSize(0.03);
  if(isdata) latex_SBND.DrawLatex(0.16, 0.96, "#font[62]{SBND Data} Run " + run_str);
  else latex_SBND.DrawLatex(0.16, 0.96, "#font[62]{SBND Simulation} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_method.DrawLatex(0.90, 0.96, particle_label_str);

  TString output_plot_dir = getenv("PLOT_PATH");
  TString outfile_str = output_plot_dir + "/MC/" + "/recom_fit/2D/rr_vs_dqdx.pdf";
  if(isdata) outfile_str = output_plot_dir + "/Run" + run_str + "/rr_vs_dqdx.pdf";
  c -> SaveAs(outfile_str);

  c -> Close();  
}
