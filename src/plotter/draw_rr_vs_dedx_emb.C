#include "canvas_margin.h"
#include "BetheBloch.h"
#include "mylib.h"
#include <iostream>

TString image_type = "pdf";

int run_num, rebin_x, rebin_y, plane;

BetheBloch *muon_BB = new BetheBloch(13);
//TSpline3 * muon_sp_range_to_KE = Get_sp_range_KE(mass_muon);

void draw_rr_vs_dedx_each_angle(TString filename, TString suffix, TString angle_bin, TString angle_str){

  TString plane_str = "";
  plane_str = TString::Format("%d", plane);

  TString plane_latex_str = "";
  if(plane == 0){
    plane_latex_str = "1st Induction Plane";
  }
  else if(plane == 1){
    plane_latex_str = "2nd Induction Plane";
  }
  else if(plane == 2){
    plane_latex_str = "Collection Plane";
  }

  setTDRStyle();

  TString input_file_dir = getenv("OUTPUTROOT_PATH");
  filename = input_file_dir + "/" + filename;
  TFile *f = new TFile(filename);
  TH2D *hist_2D = (TH2D*)gDirectory -> Get("rr_vs_dedx_plane" + plane_str + "_trklen_60cm_passing_cathode_coszx" + angle_bin);
  hist_2D -> RebinX(rebin_x);
  hist_2D -> RebinY(rebin_y);

  TH2D *hits_pitch = (TH2D*)gDirectory -> Get("rr_vs_pitch_plane" + plane_str + "_trklen_60cm_passing_cathode_coszx" + angle_bin);
  double mean_pitch = hits_pitch ->GetMean(2);
  cout << "mean_pitch : " << mean_pitch << endl;

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  gStyle->SetLineWidth(2);
  c -> SetRightMargin(0.15);

  double z_max = hist_2D -> GetMaximum();
  TH2D * template_h = new TH2D("", "", 1., 0., 200., 1., 0., 10.);
  template_h -> SetBinContent(1, -999999999999.);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("Residual Range [cm]");
  template_h -> GetXaxis() -> SetTitleSize(0.045);
  template_h -> GetXaxis() -> SetTitleOffset(1.0);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(plane_latex_str + "#color[0]{a}#frac{dE}{dx} [MeV/cm]");
  template_h -> GetYaxis() -> SetTitleOffset(0.85);
  template_h -> GetYaxis() -> SetTitleSize(0.045);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetZaxis() -> SetRangeUser(0., z_max * 1.1);
  template_h -> GetZaxis() -> SetTitle("Entries");
  template_h -> GetZaxis() -> SetTitleSize(0.045);
  template_h -> GetZaxis() -> SetLabelSize(0.035);
  template_h -> Draw("colz");

  int N_binsX = hist_2D -> GetNbinsX();
  int N_binsY = hist_2D -> GetNbinsY();
  for(int i = 1; i < N_binsX + 1; i++){
    for(int j = 1; j < N_binsY + 1; j++){
      double this_content = hist_2D -> GetBinContent(i, j);
      if(this_content < 1e-9) hist_2D -> SetBinContent(i, j, 1e-9);
      if(i == 1) hist_2D -> SetBinContent(i, j, 1e-9);
    }
  }
  hist_2D -> Draw("colsame");

  int N_BB_points = 100.;
  double rr_min = 2.;
  double rr_max = 200.;
  double rr_step = (rr_max - rr_min) / (N_BB_points + 0.);
  vector<double> BB_rr;
  vector<double> BB_dedx_MPV;
  for(int i = 0; i < N_BB_points; i++){
    double this_rr = rr_min + rr_step * (i + 0.);
    double this_KE = muon_BB -> KEFromRangeSpline(this_rr);
    double gamma = (this_KE/mass_muon)+1.0;
    double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
    double this_xi = muon_BB -> Landau_xi(this_KE, mean_pitch);
    double this_Wmax = muon_BB -> Get_Wmax(this_KE);
    double this_kappa = this_xi / this_Wmax;
    double this_dEdx_BB = muon_BB -> meandEdx(this_KE);
    double par[5] = {this_kappa, beta * beta, this_xi, this_dEdx_BB, mean_pitch};
    TF1 * this_dEdx_PDF = muon_BB -> dEdx_PDF(this_KE, mean_pitch);
    double this_dEdx_MPV = this_dEdx_PDF -> GetMaximumX();
    delete this_dEdx_PDF;

    BB_rr.push_back(this_rr);
    BB_dedx_MPV.push_back(this_dEdx_MPV);
  }
  
  TGraph *gr_BB_dedx = new TGraph(N_BB_points, &BB_rr[0], &BB_dedx_MPV[0]);
  gr_BB_dedx -> SetLineColor(kRed);
  gr_BB_dedx -> SetLineStyle(7);
  gr_BB_dedx -> SetLineWidth(2);
  gr_BB_dedx -> Draw("lsame");

  TLine *line_legd = new TLine(1, 1, 1, 1);
  line_legd -> SetLineColor(kRed-4);
  line_legd -> SetLineStyle(7);
  line_legd -> SetLineWidth(4);

  //TLegend *l = new TLegend(0.45, 0.80, 0.80, 0.92);                                                                                                                                                                                                                           
  TLegend *l = new TLegend(0.16, 0.80, 0.85, 0.95);
  l -> SetNColumns(2);
  l -> AddEntry(gr_BB_dedx, "                   ", "");
  l -> AddEntry(line_legd, "#color[0]{Muon Expectation }", "l");
  l -> SetBorderSize(0);
  l -> SetFillStyle(0);
  l -> Draw("same");

  //TString particle_label_str = "Cathode-Crossing Stopping Tracks, All#color[0]{a}#phi(E)";
  TString particle_label_str = "Cathode-Crossing Stopping Tracks";
  if(angle_bin != ""){
    particle_label_str = "Cathode-Crossing Stopping Tracks,#color[0]{a}#phi(E) = " + angle_str + "#circ";
  }
  
  TLatex latex_SBND, latex_method, latex_plane;
  latex_SBND.SetNDC();
  latex_method.SetNDC();
  latex_plane.SetNDC();
  latex_method.SetTextAlign(31);
  latex_plane.SetTextAlign(31);
  latex_SBND.SetTextSize(0.03);
  latex_method.SetTextSize(0.03);
  latex_plane.SetTextSize(0.05);
  TString sbnd_sample_str = "Data";
  if(suffix.Contains("MC")) sbnd_sample_str = "MC Simulation";
  latex_SBND.DrawLatex(0.16, 0.96, sbnd_sample_str);
  latex_method.DrawLatex(0.90, 0.96, particle_label_str);

  TLatex latex_SBND_prel;
  latex_SBND_prel.SetNDC();
  latex_SBND_prel.SetTextSize(0.04);
  latex_SBND_prel.DrawLatex(0.20, 0.865, "#color[0]{SBND Preliminary}");
  
  TString output_plot_dir = getenv("PLOT_PATH");
  TString outfile_str = output_plot_dir + "/rr_vs_dedx/" + suffix + "_rr_vs_dedx_plane" + plane_str + angle_bin + "." + image_type;

  TString outfile_dir = gSystem->DirName(outfile_str);
  if (gSystem->AccessPathName(outfile_dir)) {
    std::cout << "Directory does not exist, creating: " << outfile_dir << std::endl;
    gSystem->mkdir(outfile_dir, kTRUE);
  }
  c -> SaveAs(outfile_str);

  c -> Close();
}

void draw_rr_vs_dedx_emb(TString filename, TString suffix, int this_rebin_x = 1, int this_rebin_y = 1, int this_plane = 2){
  //void draw_rr_vs_dedx_emb(int this_run_num = 0, int this_rebin_x = 1, int this_rebin_y = 1, int this_plane = 2){

  //run_num = this_run_num;
  rebin_x = this_rebin_x;
  rebin_y = this_rebin_y;
  plane = this_plane;

  setTDRStyle();

  int N_angle_bins = 7;
  TString angle_bins[] = {"", "_phi40to50", "_phi50to60", "_phi60to70", "_phi70to80", "_phi80to85", "_phi85to90"};
  TString angle_strs[] = {"Inclusive", "40 - 50", "50 - 60", "60 - 70", "70 - 80", "80 - 85", "85 - 90"}; 
  for(int i = 0; i < N_angle_bins; i++){
    draw_rr_vs_dedx_each_angle(filename, suffix, angle_bins[i], angle_strs[i]);
  }
}
