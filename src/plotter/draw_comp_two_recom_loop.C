#include "canvas_margin.h"
#include "mylib.h"
#include "LanGausFit.h"
#include "TRandom3.h"
#include "BetheBloch.h"
#include <iostream>

TString comp_dir = "2024B_vs_2025A_Spring";
TString image_type = "pdf";

TRandom3 gRan(1800);
map<TString, vector<double>> fitting_results;

map<TString, std::vector<double> > MPV_vec_map;
map<TString, std::vector<double> > MPV_err_vec_map;
map<TString, std::vector<double> > dEdx_vec_map;
map<TString, std::vector<double> > dEdx_err_vec_map;

vector<double> MPV_vec;
vector<double> MPV_err_vec;
vector<double> dEdx_vec;
vector<double> dEdx_err_vec;

vector<double> vec_c_cal;
vector<double> vec_c_cal_err;

double alpha_ebm = 0.904;
double beta_90 = 0.204;
double R_ebm = 1.25;

double B_phi(double phi = 90.){
  double phi_rad = phi * TMath::Pi() / 180.;
  double out = beta_90 / (sqrt(pow(sin(phi_rad), 2.) + pow(cos(phi_rad) / R_ebm, 2.)));
  return out;
}

void Write_1D_hist(TH1D *h1, TH1D *h2, TString lgd_str1, TString lgd_str2, TString latex_str, TString title_x, TString title_y, double x_min, double x_max, int rebin, TString outname, bool do_langau_fit = false){

  TCanvas *c = new TCanvas("", "", 1600, 1200);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  h1 -> Rebin(rebin);
  h2 -> Rebin(rebin);
  double max_y = h1 -> GetMaximum();
  TH1D * template_h = new TH1D("", "", 1., x_min, x_max);
  template_h -> SetStats(0);
  template_h -> GetYaxis() -> SetRangeUser(0., max_y * 1.5);
  template_h -> GetXaxis() -> SetTitle(title_x);
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(title_y);
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> Draw();

  h1 -> SetMarkerColor(kRed);
  h1 -> SetMarkerStyle(32);
  h1 -> SetMarkerSize(0.7);
  h1 -> SetLineColor(kRed);
  h1 -> SetLineWidth(1);
  h1 -> Draw("epsame");

  h2 -> SetMarkerColor(kBlue);
  h2 -> SetMarkerStyle(32);
  h2 -> SetMarkerSize(0.7);
  h2 -> SetLineColor(kBlue);
  h2 -> SetLineWidth(1);
  h2 -> Draw("epsame");

  TLegend *l = new TLegend(0.70, 0.70, 0.850, 0.85);
  l -> AddEntry(h1, lgd_str1, "ep"); 
  l -> AddEntry(h2, lgd_str2, "ep");
  l -> Draw("same");

  TLatex latex_sbnd, latex_rr;
  latex_sbnd.SetNDC();
  latex_rr.SetNDC();
  latex_rr.SetTextAlign(31);
  latex_sbnd.SetTextSize(0.03);
  latex_rr.SetTextSize(0.03);
  latex_sbnd.DrawLatex(0.16, 0.96, "#font[62]{SBND Simulation 2024B and 2025A} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_rr.DrawLatex(0.95, 0.96, latex_str);

  TString output_plot_dir = getenv("PLOT_PATH");
  TString outfile_str = output_plot_dir + "/" + comp_dir + "/rr_vs_dqdx/" + outname + "." + image_type;
  TString outfile_dir = gSystem->DirName(outfile_str);
  if (gSystem->AccessPathName(outfile_dir)) {
    std::cout << "Directory does not exist, creating: " << outfile_dir << std::endl;
    gSystem->mkdir(outfile_dir, kTRUE);
  }
  c -> SaveAs(outfile_str);
  c -> Close();

}

void Fit_and_compare(TString file1, TString file2, TString comp_2d_hist, int rebin_x, int rebin_y, double res_range_cut){

  TString input_file_dir = getenv("OUTPUTROOT_PATH");

  TFile *f1 = new TFile(input_file_dir + "/" + file1);
  TH2D *hist2d_1 = (TH2D*)gDirectory -> Get(comp_2d_hist);

  TFile *f2 = new TFile(input_file_dir + "/" + file2);
  TH2D *hist2d_2 = (TH2D*)gDirectory -> Get(comp_2d_hist);

  hist2d_1 -> RebinX(rebin_x);
  hist2d_1 -> RebinY(rebin_y);
  hist2d_2 -> RebinX(rebin_x);
  hist2d_2 -> RebinY(rebin_y);

  int N_binsX = hist2d_1 -> GetNbinsX();
  int N_binsY = hist2d_1 -> GetNbinsY();

  for(int i = 1; i < N_binsX + 1; i++){

    TString i_str = Form("%d", i);
    double this_ResRange = hist2d_1 -> GetXaxis() -> GetBinCenter(i);
    if(this_ResRange > res_range_cut) break;
    double this_ResRange_err = 0.5 * hist2d_1 -> GetXaxis() -> GetBinWidth(i);
    TString ResRange_range_str = Form("ResRange%.1fto%.1fcm", this_ResRange - this_ResRange_err, this_ResRange + this_ResRange_err);
    TString ResRange_range_latex = Form("Residual range : %.1f - %.1f cm", this_ResRange -this_ResRange_err, this_ResRange + this_ResRange_err);
    TString this_hist_name = ResRange_range_str;

    TH1D * this_1_1D = new TH1D("h1_" + this_hist_name, this_hist_name, N_binsY, 0., 3000.);
    TH1D * this_2_1D =new TH1D("h2_" + this_hist_name, this_hist_name, N_binsY, 0., 3000.);

    for(int j = 1; j < N_binsY + 1; j++){
      double this_MC_content = hist2d_1 -> GetBinContent(i, j);
      double this_MC_error = hist2d_1 -> GetBinError(i, j);
      this_1_1D -> SetBinContent(j, this_MC_content);
      this_1_1D -> SetBinError(j, this_MC_error);

      double this_Data_content = hist2d_2 -> GetBinContent(i, j);
      double this_Data_error = hist2d_2 -> GetBinError(i, j);
      this_2_1D -> SetBinContent(j, this_Data_content);
      this_2_1D -> SetBinError(j, this_Data_error);
    }

    double max_y = this_2_1D -> GetMaximum();
    this_1_1D -> Scale(max_y / this_1_1D -> GetMaximum());

    double x_range_low = 500.;
    double x_range_high = 2000.;

    if(this_ResRange < 15.) x_range_high = 3000.;

    // Write_1D_hist(TH1D *h1, TH1D *h2, TString lgd_str1, TString lgd_str2, TString latex_str, TString title_x, TString title_y, double x_min, double x_max, int rebin, TString outname, bool do_langau_fit)
    Write_1D_hist(this_1_1D, this_2_1D, "MC 2024B", "MC 2025A Spring", ResRange_range_latex, "dQ/dx [ADC* #times tick / cm]", "A.U. (peak norm.)", x_range_low, x_range_high, 1., "dQdx_" + ResRange_range_str, false);
  }

}

void draw_comp_two_recom_loop(){

  setTDRStyle();
  TString file1 = "output_recom_loop_emb_mc_2024b.root";
  TString file2 = "output_recom_loop_emb_mc_2025a_spring.root";
  TString comp_2d_hist = "rr_vs_corr_dqdx_plane2_trklen_60cm_passing_cathode_coszx";
  Fit_and_compare(file1, file2, comp_2d_hist, 2, 40, 100.);
}
