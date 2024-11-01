#include "canvas_margin.h"
#include "mylib.h"

void draw(TString filename, TString evt, TString id, TString plane, TString wire){

  TString input_file_dir = getenv("OUTPUTROOT_PATH");
  TFile *f = new TFile(input_file_dir + "/" + filename);
  TString this_hist_str = "evt" + evt + "_id" + id + "_plane" + plane + "_wire" + wire;
  TH1D* this_h = (TH1D*)gDirectory -> Get(this_hist_str);

  TCanvas *c = new TCanvas("", "", 900, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  gStyle -> SetLineWidth(2);
  
  double y_min = this_h -> GetMinimum();
  double y_max = this_h -> GetMaximum();
  double y_up = y_max + 0.2 * (y_max - y_min);
  double y_down = y_min - 0.2 * (y_max - y_min);
  
  this_h -> SetStats(0);
  this_h -> GetYaxis() -> SetRangeUser(y_down, y_up);
  this_h -> GetXaxis() -> SetTitle("Time [tick]");
  this_h -> GetXaxis() -> SetTitleSize(0.037);
  this_h -> GetXaxis() -> SetTitleOffset(1.4);
  this_h -> GetXaxis() -> SetLabelSize(0.035);
  this_h -> GetYaxis() -> SetTitle("ADC");
  this_h -> GetYaxis() -> SetTitleSize(0.05);
  this_h -> GetYaxis() -> SetLabelSize(0.035);
  this_h -> SetLineColor(kBlack);
  this_h -> SetLineWidth(2);
  this_h -> Draw("hist");

  TLatex latex_ProtoDUNE, latex_particle;
  latex_ProtoDUNE.SetNDC();
  latex_particle.SetNDC();
  latex_particle.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_particle.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND Data} Run 14860");
  latex_particle.DrawLatex(0.95, 0.96, "Evt: " + evt + ", Track ID: " + id + ", Plane" + plane + ", Wire: " + wire);
  
  TString output_plot_dir = getenv("PLOT_PATH");
  TString outfile_str = output_plot_dir + "/waveform/" + this_hist_str + ".pdf";
  c -> SaveAs(outfile_str);
  c -> Close();
}

void draw_waveform(){

  setTDRStyle();
  TString input_file = "output_waveform_14860_evt8473_id2.root";
  draw(input_file, "8473", "2", "0", "1462");
  draw(input_file, "8473", "2", "1", "643");
  draw(input_file, "8473", "2", "2", "1161");
}
