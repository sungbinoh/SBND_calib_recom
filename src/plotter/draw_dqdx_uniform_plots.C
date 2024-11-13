#include "canvas_margin.h"
#include "mylib.h"
#include "LanGausFit.h"
#include "TRandom3.h"
#include "BetheBloch.h"
#include <iostream>

bool isdata = false;
TString run_str = "";

double get_y_deg40to45(TGraphErrors *gr){

  double y_40to45 = -1.;
  int nPoints = gr -> GetN();
  const double* x_ptr = gr -> GetX();
  const double* y_ptr = gr -> GetY();
  std::vector<double> x_plane0(x_ptr, x_ptr + nPoints);
  std::vector<double> y_plane0(y_ptr, y_ptr + nPoints);
  for (int i = 0; i < nPoints; ++i) {
    double this_x = x_plane0[i];
    if(this_x > 40. && this_x < 45.) y_40to45 = y_plane0[i];
  }

  return y_40to45;
}

void draw(double x_min, double x_max, double y_min, double y_max, TString suffix, TString y_title){

  
  TString input_file_dir = getenv("OUTPUTROOT_PATH");
  TFile *f = new TFile(input_file_dir + "/dqdx_unif_" + run_str + "_selected1or2.root");
  TGraphErrors *gr_plane0_east = (TGraphErrors*)gDirectory -> Get("plane0_east_" + suffix);
  TGraphErrors *gr_plane1_east = (TGraphErrors*)gDirectory -> Get("plane1_east_" + suffix);
  TGraphErrors *gr_plane2_east = (TGraphErrors*)gDirectory -> Get("plane2_east_" + suffix);
  
  TGraphErrors *gr_plane0_west = (TGraphErrors*)gDirectory -> Get("plane0_west_" + suffix);
  TGraphErrors *gr_plane1_west = (TGraphErrors*)gDirectory -> Get("plane1_west_" + suffix);
  TGraphErrors *gr_plane2_west = (TGraphErrors*)gDirectory -> Get("plane2_west_" + suffix);
  
  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  gStyle->SetLineWidth(2);
  
  TH1D * template_h = new TH1D("", "", 1., x_min, x_max);
  template_h -> SetStats(0);
  template_h -> GetYaxis() -> SetRangeUser(y_min, y_max);
  template_h -> GetXaxis() -> SetTitle("Zenith Angle [#circ]");
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(y_title);
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> Draw();

  double plane0_40to45_y = get_y_deg40to45(gr_plane0_east);
  double plane1_40to45_y = get_y_deg40to45(gr_plane1_east);
  double plane2_40to45_y = get_y_deg40to45(gr_plane2_east);

  double percentage = 0.02;
  TBox * box_1p_plane0_east = new TBox(x_min, plane0_40to45_y - plane0_40to45_y * percentage, x_max, plane0_40to45_y + plane0_40to45_y * percentage);
  box_1p_plane0_east -> SetFillColorAlpha(kAzure+8, 0.8);
  box_1p_plane0_east -> Draw("same");

  TBox * box_1p_plane1_east = new TBox(x_min, plane1_40to45_y - plane1_40to45_y * percentage, x_max, plane1_40to45_y + plane1_40to45_y * percentage);
  box_1p_plane1_east -> SetFillColorAlpha(kRed-10, 0.8);
  box_1p_plane1_east -> Draw("same");

  TBox * box_1p_plane2_east = new TBox(x_min, plane2_40to45_y - plane2_40to45_y * percentage, x_max, plane2_40to45_y + plane2_40to45_y * percentage);
  box_1p_plane2_east -> SetFillColorAlpha(kGreen-10, 0.8);
  box_1p_plane2_east -> Draw("same");
  
  gr_plane0_east -> SetLineColor(kBlue);
  gr_plane0_east -> SetMarkerColor(kBlue);
  gr_plane0_east -> SetMarkerStyle(22);
  gr_plane0_east -> Draw("epsame");

  gr_plane1_east -> SetLineColor(kRed);
  gr_plane1_east -> SetMarkerColor(kRed);
  gr_plane1_east -> SetMarkerStyle(22);
  gr_plane1_east -> Draw("epsame");

  gr_plane2_east -> SetLineColor(kGreen+2);
  gr_plane2_east -> SetMarkerColor(kGreen+2);
  gr_plane2_east -> SetMarkerStyle(22);
  gr_plane2_east -> Draw("epsame");

  gr_plane0_west -> SetLineColor(kBlue);
  gr_plane0_west -> SetMarkerColor(kBlue);
  gr_plane0_west -> SetMarkerStyle(32);
  gr_plane0_west -> SetLineStyle(2);
  gr_plane0_west -> Draw("epsame");

  gr_plane1_west -> SetLineColor(kRed);
  gr_plane1_west -> SetMarkerColor(kRed);
  gr_plane1_west -> SetMarkerStyle(32);
  gr_plane1_west -> SetLineStyle(2);
  gr_plane1_west -> Draw("epsame");

  gr_plane2_west -> SetLineColor(kGreen+2);
  gr_plane2_west -> SetMarkerColor(kGreen+2);
  gr_plane2_west -> SetMarkerStyle(32);
  gr_plane2_west -> SetLineStyle(2);
  gr_plane2_west -> Draw("epsame");
  
  TH1D *east_temp = new TH1D("", "", 1., x_min, x_max);
  east_temp -> SetMarkerStyle(22);
  east_temp -> SetFillColorAlpha(kBlack, 0.8);
  TH1D *west_temp = new TH1D("", "", 1., x_min, x_max);
  west_temp -> SetMarkerStyle(32);
  west_temp -> SetLineStyle(2);
  TLegend *l = new TLegend(0.20, 0.70, 0.90, 0.90);
  l -> SetNColumns(2);
  l -> AddEntry(gr_plane0_east, "#color[13]{Plane 0}", "lp");
  l -> AddEntry(east_temp, "#color[13]{Track coming from East}", "lp");
  l -> AddEntry(gr_plane1_east, "#color[13]{Plane 1}", "lp");
  l -> AddEntry(west_temp, "#color[13]{Track coming from West}", "lp");
  l -> AddEntry(gr_plane2_east, "#color[13]{Plane 2}", "lp");
  l -> AddEntry(east_temp, "#color[13]{<dQ/dx>_{Trunc.}^{From East.} #pm 2% of Zenith Angle = [40#circ, 45#circ)}", "f");
		
  l -> Draw("same");

  gPad->RedrawAxis();
  
  TString sample_str = "";
  if(isdata){
    sample_str = "#font[62]{SBND Data} Run " + run_str;
    if(sample_str.Contains("data")) sample_str = "#font[62]{SBND Data} Run 14608, 14685, 14724, 14833 and 14860";
  }
  else sample_str = "#font[62]{SBND MC} 2023B CV";
  
  TLatex latex_ProtoDUNE, latex_particle, latex_particle2;
  latex_ProtoDUNE.SetNDC();
  latex_particle.SetNDC();
  latex_particle.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_particle.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, sample_str);
  latex_particle.DrawLatex(0.93, 0.90, "Cathode-Crossing Through-Going Tracks");

  latex_particle2.SetNDC();
  latex_particle2.SetTextAlign(31);
  latex_particle2.SetTextSize(0.03);
  //latex_particle2.DrawLatex(0.93, 0.87, "But Not Crossing the Anode");
  
  TString output_plot_dir = getenv("PLOT_PATH");
  TString outfile_str = output_plot_dir + "/comparison/dqdx_unif/plane_n_zenith_comparison_dqdx_trun_mean_" + run_str + "_selected1or2_" + suffix + ".pdf";
  if(!isdata) outfile_str = output_plot_dir + "/comparison/dqdx_unif/plane_n_zenith_comparison_dqdx_trun_mean_" + run_str + "_selected1or2_MC_"  + suffix + ".pdf";
  c -> SaveAs(outfile_str);
  c -> Close();
  
}

void draw_each_plane(double x_min, double x_max, double y_min, double y_max, TString plane, TString y_title){

  TString input_file_dir = getenv("OUTPUTROOT_PATH");
  TFile *f = new TFile(input_file_dir + "/dqdx_unif_" + run_str + "_selected1or2.root");
  TGraphErrors *gr_east_NW = (TGraphErrors*)gDirectory -> Get(plane + "_east_NW");
  TGraphErrors *gr_east_NE = (TGraphErrors*)gDirectory -> Get(plane + "_east_NE");
  TGraphErrors *gr_east_Core = (TGraphErrors*)gDirectory -> Get(plane + "_east_TPCcore");

  TGraphErrors *gr_west_NW = (TGraphErrors*)gDirectory -> Get(plane + "_west_NW");
  TGraphErrors *gr_west_NE = (TGraphErrors*)gDirectory -> Get(plane + "_west_NE");
  TGraphErrors *gr_west_Core = (TGraphErrors*)gDirectory -> Get(plane + "_west_TPCcore");

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  gStyle->SetLineWidth(2);

  TH1D * template_h = new TH1D("", "", 1., x_min, x_max);
  template_h -> SetStats(0);
  template_h -> GetYaxis() -> SetRangeUser(y_min, y_max);
  template_h -> GetXaxis() -> SetTitle("Zenith Angle [#circ]");
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(y_title);
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> Draw();

  gr_east_NW -> SetLineColor(kBlue);
  gr_east_NW -> SetMarkerColor(kBlue);
  gr_east_NW -> SetMarkerStyle(22);
  gr_east_NW -> Draw("epsame");

  gr_east_NE -> SetLineColor(kRed);
  gr_east_NE -> SetMarkerColor(kRed);
  gr_east_NE -> SetMarkerStyle(22);
  gr_east_NE -> Draw("epsame");

  gr_east_Core -> SetLineColor(kGreen+2);
  gr_east_Core -> SetMarkerColor(kGreen+2);
  gr_east_Core -> SetMarkerStyle(22);
  gr_east_Core -> Draw("epsame");

  gr_west_NW -> SetLineColor(kBlue);
  gr_west_NW -> SetMarkerColor(kBlue);
  gr_west_NW -> SetMarkerStyle(32);
  gr_west_NW -> SetLineStyle(2);
  gr_west_NW -> Draw("epsame");

  gr_west_NE -> SetLineColor(kRed);
  gr_west_NE -> SetMarkerColor(kRed);
  gr_west_NE -> SetMarkerStyle(32);
  gr_west_NE -> SetLineStyle(2);
  gr_west_NE -> Draw("epsame");

  gr_west_Core -> SetLineColor(kGreen+2);
  gr_west_Core -> SetMarkerColor(kGreen+2);
  gr_west_Core -> SetMarkerStyle(32);
  gr_west_Core -> SetLineStyle(2);
  gr_west_Core -> Draw("epsame");

  TH1D *east_temp = new TH1D("", "", 1., x_min, x_max);
  east_temp -> SetMarkerStyle(22);
  east_temp -> SetFillColorAlpha(kBlack, 0.8);
  TH1D *west_temp = new TH1D("", "", 1., x_min, x_max);
  west_temp -> SetMarkerStyle(32);
  west_temp -> SetLineStyle(2);
  TLegend *l = new TLegend(0.20, 0.70, 0.90, 0.90);
  l -> SetNColumns(2);
  l -> AddEntry(gr_east_NW, "#color[13]{Hits in TPC NW}", "lp");
  l -> AddEntry(east_temp, "#color[13]{Track coming from East}", "lp");
  l -> AddEntry(gr_east_NE, "#color[13]{Hits in TPC NE}", "lp");
  l -> AddEntry(west_temp, "#color[13]{Track coming from West}", "lp");
  l -> AddEntry(gr_east_Core, "#color[13]{Hits in TPC Core}", "lp");
  l -> AddEntry(east_temp, "#color[13]{<dQ/dx>_{Trunc.}^{From East.} #pm 2% of Zenith Angle = [40#circ, 45#circ)}", "f");

  l -> Draw("same");

  gPad->RedrawAxis();

  TString sample_str = "";
  if(isdata){
    sample_str = "#font[62]{SBND Data} Run " + run_str;
    if(sample_str.Contains("data")) sample_str = "#font[62]{SBND Data} Run 14608, 14685, 14724, 14833 and 14860";
  }
  else sample_str = "#font[62]{SBND MC} 2023B CV";

  TLatex latex_ProtoDUNE, latex_particle, latex_particle2;
  latex_ProtoDUNE.SetNDC();
  latex_particle.SetNDC();
  latex_particle.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_particle.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, sample_str);
  latex_particle.DrawLatex(0.93, 0.90, "Cathode-Crossing Through-Going Tracks");

  latex_particle2.SetNDC();
  latex_particle2.SetTextAlign(31);
  latex_particle2.SetTextSize(0.03);

  TString output_plot_dir = getenv("PLOT_PATH");
  TString outfile_str = output_plot_dir + "/comparison/dqdx_unif/tpc_part_n_zenith_comparison_dqdx_trun_mean_" + run_str + "_selected1or2_" + plane + ".pdf";
  if(!isdata) outfile_str = output_plot_dir + "/comparison/dqdx_unif/tpc_part_n_zenith_comparison_dqdx_trun_mean_" + run_str + "_selected1or2_MC_"  + plane + ".pdf";
  c -> SaveAs(outfile_str);
  c -> Close();
}

void draw_dqdx_uniform_plots(int run_num = 0){

  if(run_num != 0){
    isdata = true;
    run_str = TString::Format("%d", run_num);
    if(run_num == -1) run_str = "data";
  }
  else run_str = "MC";
  
  setTDRStyle();

  
  TString suffixes[] = {"TPCcore", "NE", "NW"};
  for(int i = 0; i < 3; i++){
    TString this_suffix = suffixes[i];
    draw(15., 70., 1000., 1350., this_suffix, "Truncated <dQ/dx> [ADC#times tick/cm]");
  }

  TString planes[] = {"plane0", "plane1", "plane2"};
  for(int i = 0; i < 3; i++){
    TString this_plane = planes[i];
    draw_each_plane(15., 70., 1000., 1350., this_plane, "Truncated <dQ/dx> [ADC#times tick/cm]");
  }
}
