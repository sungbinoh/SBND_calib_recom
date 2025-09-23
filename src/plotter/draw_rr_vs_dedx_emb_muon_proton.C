#include "canvas_margin.h"
#include "BetheBloch.h"
#include "mylib.h"
#include <iostream>

TString image_type = "pdf";

int run_num, rebin_x, rebin_y, plane;

BetheBloch *muon_BB = new BetheBloch(13);
BetheBloch *proton_BB = new BetheBloch(2212);

//TSpline3 * muon_sp_range_to_KE = Get_sp_range_KE(mass_muon);

void draw_rr_vs_dedx_each_angle(TString filename, TString suffix, TString angle_bin, TString angle_str, TString angle_bin_proton){

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

  TString proton_file_path = "/exp/sbnd/data/users/munjung/calibration/recomb_output/2025A/";
  TString proton_file_prefix = "Dev_twoprong_";
  if(filename.Contains("MC")){
    proton_file_prefix = "MCP2025A_twoprong_";
  }

  if(plane == 0 || plane == 1){
    proton_file_prefix = proton_file_prefix + "plane" + plane_str + "_cleaned";
  }
  else{
    proton_file_prefix = proton_file_prefix + "cleaned";
  }
  
  TFile *f_proton = new TFile(proton_file_path + proton_file_prefix + angle_bin_proton + "_hist.root");
  TH2D *hist_2D_proton = (TH2D*)gDirectory -> Get("rr_vs_dedx_proton");
  //hist_2D_proton -> RebinX(rebin_x);
  hist_2D_proton -> RebinY(25);

  TH2D *hits_pitch_proton = (TH2D*)gDirectory -> Get("rr_vs_pitch_proton");
  double mean_pitch_proton = hits_pitch_proton ->GetMean(2);
  
  double z_max = hist_2D -> GetMaximum();
  double z_max_proton = hist_2D_proton -> GetMaximum();
  double scale_proton = z_max / z_max_proton;
  cout << "z_max: " <<z_max << ", z_max_proton: " << z_max_proton << endl;
  
  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  gStyle->SetLineWidth(2);
  c -> SetRightMargin(0.15);


  TH2D * template_h = new TH2D("", "", 1., 0., 200., 1., 0., 20.);
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
  template_h -> GetZaxis() -> SetTitle("Hits");
  template_h -> GetZaxis() -> SetTitleSize(0.045);
  template_h -> GetZaxis() -> SetLabelSize(0.035);
  template_h -> Draw("colz");

  int N_binsX_proton = hist_2D_proton -> GetNbinsX();
  
  int N_binsX = hist_2D -> GetNbinsX();
  int N_binsY = hist_2D -> GetNbinsY();
  for(int i = 1; i < N_binsX + 1; i++){
    for(int j = 1; j < N_binsY + 1; j++){
      double this_content = hist_2D -> GetBinContent(i, j);
      if(i < N_binsX_proton + 1){
	double this_proton_content = hist_2D_proton -> GetBinContent(i, j);
	this_content = this_content + this_proton_content * scale_proton;
	
      }

      if(this_content < 1e-9) hist_2D -> SetBinContent(i, j, 1e-9);
      else hist_2D -> SetBinContent(i, j, this_content);

      if(i == 1) hist_2D -> SetBinContent(i, j, 1e-9);
    }
  }

  z_max = hist_2D -> GetMaximum();
  template_h -> GetZaxis() -> SetRangeUser(0., z_max * 1.1);


  hist_2D -> Draw("colsame");

  //hist_2D_proton -> Draw("colsame");

  int N_BB_points = 100.;
  double rr_min = 2.;
  double rr_max = 200.;
  double rr_step = (rr_max - rr_min) / (N_BB_points + 0.);
  vector<double> BB_rr;
  vector<double> BB_dedx_MPV;
  vector<double> BB_dedx_MPV_proton;

  for(int i = 0; i < N_BB_points; i++){
    double this_rr = rr_min + rr_step * (i + 0.);
    double this_KE = muon_BB -> KEFromRangeSpline(this_rr);
    TF1 * this_dEdx_PDF = muon_BB -> dEdx_PDF(this_KE, mean_pitch);
    double this_dEdx_MPV = this_dEdx_PDF -> GetMaximumX();
    delete this_dEdx_PDF;

    double this_KE_proton =  proton_BB -> KEFromRangeSpline(this_rr);
    TF1 * this_dEdx_PDF_proton = proton_BB -> dEdx_PDF(this_KE_proton, mean_pitch_proton);
    double this_dEdx_MPV_proton = this_dEdx_PDF_proton -> GetMaximumX();
    delete this_dEdx_PDF_proton;
    
    BB_rr.push_back(this_rr);
    BB_dedx_MPV.push_back(this_dEdx_MPV);
    BB_dedx_MPV_proton.push_back(this_dEdx_MPV_proton);
  }
  
  TGraph *gr_BB_dedx = new TGraph(N_BB_points, &BB_rr[0], &BB_dedx_MPV[0]);
  gr_BB_dedx -> SetLineColor(kRed);
  gr_BB_dedx -> SetLineStyle(7);
  gr_BB_dedx -> SetLineWidth(2);
  gr_BB_dedx -> Draw("lsame");


  TGraph *gr_BB_dedx_proton = new TGraph(N_BB_points, &BB_rr[0], &BB_dedx_MPV_proton[0]);
  gr_BB_dedx_proton -> SetLineColor(kRed);
  gr_BB_dedx_proton -> SetLineStyle(9);
  gr_BB_dedx_proton -> SetLineWidth(2);
  gr_BB_dedx_proton -> Draw("lsame");
  
  TLine *line_legd = new TLine(1, 1, 1, 1);
  line_legd -> SetLineColor(kRed-4);
  line_legd -> SetLineStyle(7);
  line_legd -> SetLineWidth(4);

  TLine *line_legd_proton = new TLine(1, 1, 1, 1);
  line_legd_proton -> SetLineColor(kRed-4);
  line_legd_proton -> SetLineStyle(9);
  line_legd_proton -> SetLineWidth(4);
  
  TLegend *l = new TLegend(0.16, 0.80, 0.85, 0.95);
  l -> SetNColumns(2);
  l -> AddEntry(gr_BB_dedx, "                   ", "");
  l -> AddEntry(line_legd, "#color[0]{Muon Expectation }", "l");
  l -> AddEntry(gr_BB_dedx, "                   ", "");
  l -> AddEntry(line_legd_proton,  "#color[0]{Proton Expectation }", "l");
  l -> SetBorderSize(0);
  l -> SetFillStyle(0);
  l -> Draw("same");

  TString particle_label_str = "Muon and Proton Candidates";
  if(angle_bin != ""){
    particle_label_str = "Cathode-Crossing Stopping Tracks,#color[0]{a}#phi(E) = " + angle_str + "#circ";
  }

  gPad->RedrawAxis();
  
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

void draw_rr_vs_dedx_emb_muon_proton(TString filename, TString suffix, int this_rebin_x = 1, int this_rebin_y = 1, int this_plane = 2){
  //void draw_rr_vs_dedx_emb(int this_run_num = 0, int this_rebin_x = 1, int this_rebin_y = 1, int this_plane = 2){

  //run_num = this_run_num;
  rebin_x = this_rebin_x;
  rebin_y = this_rebin_y;
  plane = this_plane;

  setTDRStyle();

  int N_angle_bins = 7;
  TString angle_bins[] = {"", "_phi40to50", "_phi50to60", "_phi60to70", "_phi70to80", "_phi80to85", "_phi85to90"};
  TString angle_bins_proton[] = {"", "_40to50", "_50to60", "_60to70", "_70to80", "_80to85", "_85to90"};
  TString angle_strs[] = {"Inclusive", "40 - 50", "50 - 60", "60 - 70", "70 - 80", "80 - 85", "85 - 90"}; 
  for(int i = 0; i < N_angle_bins; i++){
    draw_rr_vs_dedx_each_angle(filename, suffix, angle_bins[i], angle_strs[i], angle_bins_proton[i]);
  }
}
