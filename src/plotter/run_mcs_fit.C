#include "canvas_margin.h"
#include "mylib.h"
#include "LanGausFit.h"
#include "TRandom3.h"
#include "BetheBloch.h"
#include <iostream>

bool isdata = false;
TString run_str = "";
//TString suffix = "";

TRandom3 gRan(1800);
map<TString, vector<double>> fitting_results;

map<TString, std::vector<double> > mu_vec_map;
map<TString, std::vector<double> > mu_err_vec_map;
map<TString, std::vector<double> > sigma_vec_map;
map<TString, std::vector<double> > sigma_err_vec_map;

Double_t HL_function(Double_t *x, Double_t *par){
  //x[0] = x[0] / 1000.;
  Double_t kappa_a = par[0];
  Double_t kappa_c = par[1];
  Double_t sigma_res = par[2];
  Double_t epsilon = par[3];
  Double_t mass = par[4] / 1000.;
  Double_t segment_size = par[5];

  Double_t kappa = ( (kappa_a * 1000000. / (x[0]*x[0])) + kappa_c );
  Double_t one_over_pbeta = pow(x[0]*x[0] + mass * mass, 0.5) / (x[0]*x[0]);
  Double_t root_term = pow(segment_size / 14., 0.5);

  Double_t func = kappa * one_over_pbeta * root_term * (1 + epsilon * log(segment_size / 14.)) * 1000.;
  func = pow(func * func + sigma_res * sigma_res, 0.5);
  return func;
}

void Write_1D_hist(TH1D *in, TString outname, TString var_str, TString suffix, TString particle, TString latex_str, TString title_x, TString title_y, double x_min, double x_max, int rebin, bool do_gaus_fit){

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  in -> Rebin(rebin);
  double max_y = in -> GetMaximum();
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

  in -> SetMarkerColor(kBlack);
  in -> SetMarkerStyle(32);
  in -> SetMarkerSize(0.7);
  in -> SetLineColor(kBlack);
  in -> SetLineWidth(1);
  in -> Draw("epsame");

  if(do_gaus_fit){

    TH1D *in_clone = (TH1D*)in -> Clone();
    double this_mean = in -> GetMean();
    double this_stddev = in -> GetStdDev();
    double bin_width = in -> GetBinWidth(1);
    double init_area = in -> Integral() * 0.05 * bin_width;
    double fit_x_min = this_mean - 1.8 * this_stddev;
    double fit_x_max = this_mean + 1.8 * this_stddev;
    //cout << "fit_x_min: " << fit_x_min << ", fit_x_max: " << fit_x_max << endl;
    
    TF1 *this_gaus_fit = new TF1("fit_gaus", "gaus", fit_x_min, fit_x_max);
    this_gaus_fit -> SetParameters(init_area, this_mean, this_stddev);
    in_clone -> Fit(this_gaus_fit, "RBOSQN", "", fit_x_min, fit_x_max);
    
    double this_const = this_gaus_fit -> GetParameter(0);
    double this_const_err = this_gaus_fit -> GetParError(0);
    double this_mu = this_gaus_fit -> GetParameter(1);
    double this_mu_err = this_gaus_fit -> GetParError(1);
    double this_std = this_gaus_fit -> GetParameter(2);
    double this_std_err = this_gaus_fit -> GetParError(2);
    double this_chi2 = this_gaus_fit -> GetChisquare() / this_gaus_fit -> GetNDF();

    TF1 *this_gaus = new TF1("fit_gaus", "gaus", fit_x_min, fit_x_max);
    this_gaus -> SetParameters(this_const, this_mu, this_std);
    this_gaus -> SetNpx(1000);
    this_gaus -> SetLineWidth(2);
    this_gaus -> SetLineStyle(7);
    this_gaus -> Draw("lsame");
    
    in -> Draw("epsame");

    TLegend *l = new TLegend(0.70, 0.40, 0.90, 0.80);
    l -> AddEntry(in, Form("Angle (%.0f)", in -> Integral()), "pl");
    l -> AddEntry(this_gaus, Form("#sigma : %.2f #pm %.2f", this_std, this_std_err), "l");
    l -> AddEntry(in, Form("#mu : %.2f #pm %.2f", this_mu, this_mu_err), "");
    l -> AddEntry(in, Form("Const : %.2f #pm %.2f", this_const, this_const_err), "");
    l -> AddEntry(in, Form("#chi^{2} / ndf : %.2f", this_chi2), "");
    l -> Draw("same");

    mu_vec_map[particle + var_str + suffix].push_back(this_mu);
    mu_err_vec_map[particle + var_str + suffix].push_back(this_mu_err);
    sigma_vec_map[particle + var_str + suffix].push_back(this_std);
    sigma_err_vec_map[particle + var_str + suffix].push_back(this_std_err);
  }
  else{
    double this_mean = in -> GetMean();
    double this_stddev = in -> GetStdDev();

    TLegend *l = new TLegend(0.70, 0.70, 0.850, 0.85);
    l -> AddEntry(in, Form("Mean : %.2f", this_mean), "pl");
    l -> AddEntry(in, Form("StdDev : %.2f", this_stddev), "");
    l -> Draw("same");
  
    mu_vec_map[particle + var_str + suffix].push_back(this_mean);
    sigma_vec_map[particle + var_str + suffix].push_back(this_stddev);
  }

  TString particle_label_str = "";
  if(particle == "muon") particle_label_str = "Cathode Passing Stopping Tracks";
  if(particle == "proton") particle_label_str = "Stopping Proton Candidates";

  TLatex latex_ProtoDUNE, latex_particle, latex_Nhits, latex_method;
  latex_ProtoDUNE.SetNDC();
  latex_particle.SetNDC();
  latex_Nhits.SetNDC();
  latex_method.SetNDC();
  latex_particle.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_particle.SetTextSize(0.03);
  latex_Nhits.SetTextSize(0.06);
  latex_method.SetTextSize(0.06);
  if(isdata) latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND Data} Run " + run_str);
  else latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND Simulation} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_particle.DrawLatex(0.95, 0.96, particle_label_str);
  latex_method.DrawLatex(0.18, 0.87, latex_str);

  TString output_plot_dir = getenv("PLOT_PATH");
  TString outfile_str = output_plot_dir + "/MC_2025A_CV/mcs/1d/" + outname + ".pdf";
  if(isdata) outfile_str = output_plot_dir + "/run_" + run_str + "/mcs/1d/" + outname + ".pdf";

  TString outfile_dir = gSystem->DirName(outfile_str);
  if (gSystem->AccessPathName(outfile_dir)) {
    std::cout << "Directory does not exist, creating: " << outfile_dir << std::endl;
    gSystem->mkdir(outfile_dir, kTRUE);
  }
  
  c -> SaveAs(outfile_str);
  c -> Close();
}

void Fit_p_vs_mcs_plots(TString input_file_name, TString suffix, TString suffix_leg, TString particle, double angle_low, double angle_high, double p_binning[], int num_NbinsX, int rebin_angle[]){

  //const Int_t num_NbinsX = sizeof(p_binning) / sizeof(p_binning[0]) - 1;
  
  TString this_id = "id";
  TString histname = "muon_p_vs" + suffix;
  
  TString input_file_dir = getenv("OUTPUTROOT_PATH");
  TFile *f = new TFile(input_file_dir + "/" + input_file_name);
  TH2D *hist_2D = (TH2D*)gDirectory -> Get(histname);

  int N_binsX = hist_2D -> GetNbinsX();
  int N_binsY = hist_2D -> GetNbinsY();

  TH1D *hist_1D_for_X = hist_2D -> ProjectionX("hist_1D_for_X");

  double hist_2D_x_low = 0.;
  double hist_2D_x_high = 3000.;
  if(particle == "proton") hist_2D_x_high = 15.;
  double hist_2D_x_step = (hist_2D_x_high - hist_2D_x_low) / N_binsX;

  double angle_min = -500.;
  double angle_max = 500.;

  cout << "num_NbinsX : " << num_NbinsX << ", N_binsX : " << N_binsX << ", hist_2D_x_step : " << hist_2D_x_step << endl;

  for(int i = 1; i < num_NbinsX - 1; i++){
      
    double this_p_low = p_binning[i];
    double this_p_high = p_binning[i + 1];

    int x_index_low = 1 + (this_p_low - hist_2D_x_low) / hist_2D_x_step;
    int x_index_high = 1 + (this_p_high - hist_2D_x_low) / hist_2D_x_step;
    //cout << "this_p_low : " << this_p_low << ", this_p_high : " << this_p_high << ", x_index_low : " << x_index_low << ", x_index_high : " << x_index_high << endl;
    //cout << << endl;

    TString p_str  = Form("p_%.2fto%.2f_MeVc", this_p_low, this_p_high);
    TString p_latex_str = Form("Momentum %.2f to %.2f MeVc", this_p_low, this_p_high); 
    TString this_1D_X_hist_name = "p_" + p_str + suffix;
    TString this_1D_Y_hist_name = "angle_" + p_str + suffix;

    TH1D * this_1D_X = new TH1D(this_1D_X_hist_name, this_1D_X_hist_name, x_index_high - x_index_low, this_p_low, this_p_high);
    TH1D * this_1D_Y = new TH1D(this_1D_Y_hist_name, this_1D_Y_hist_name, N_binsY, angle_min, angle_max);

    for(int j = 0; j < x_index_high - x_index_low; j++){
      double this_1D_X_hist_content = hist_1D_for_X -> GetBinContent(x_index_low + j);
      double this_1D_X_hist_err = hist_1D_for_X -> GetBinError(x_index_low + j);
      //cout << "this_1D_X_hist_content : " << this_1D_X_hist_content << ", this_1D_X_hist_err : " << this_1D_X_hist_err << ", sqrt(this_1D_X_hist_content) : " << sqrt(this_1D_X_hist_content) << endl;
      this_1D_X -> SetBinContent(j + 1, this_1D_X_hist_content);
      this_1D_X -> SetBinError(j + 1, this_1D_X_hist_err);
    }

    for(int k = 1; k < N_binsY + 1; k++){
      double this_1D_Y_content = 0.;
      double this_1D_Y_err = 0.;
      for(int j = 0; j < x_index_high - x_index_low; j++){
	double this_content = hist_2D -> GetBinContent(x_index_low + j, k);
	double this_err = hist_2D -> GetBinError(x_index_low + j, k);

	this_1D_Y_content = this_1D_Y_content + this_content;
	this_1D_Y_err = this_1D_Y_err + pow(this_err, 2.);
      }

      this_1D_Y_err = pow(this_1D_Y_err, 0.5);
      this_1D_Y -> SetBinContent(k, this_1D_Y_content);
      this_1D_Y -> SetBinError(k, this_1D_Y_err);
      
    }

    Write_1D_hist(this_1D_X, "/" + this_1D_X_hist_name, "momentum", suffix, particle, p_latex_str, "Momentum [MeV/c]", "Events", this_p_low, this_p_high, 1., false);
    Write_1D_hist(this_1D_Y, "/" + this_1D_Y_hist_name, "angle", suffix, particle, p_latex_str, "Angle [mrad.]", "Events", -300., 300., rebin_angle[i - 1], true);
  }

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  
  double p_low = 0.;
  double p_high = 2000.;
  double z_max = hist_2D -> GetMaximum();
  TH1D * template_h = new TH1D("", "", 1., p_low, p_high);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("Momentum [MeV/c]");
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(suffix_leg + " [mrad.]");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., 70.);
  template_h -> Draw("colz");

  //hist_2D -> Draw("colzsame");
  
  TGraphErrors * sigma_gr = new TGraphErrors(mu_vec_map[particle+"momentum"+suffix].size(),
					    &mu_vec_map[particle+"momentum"+suffix][0], &sigma_vec_map[particle+"angle"+suffix][0],
					    &mu_err_vec_map[particle+"momentum"+suffix][0], &sigma_err_vec_map[particle+"angle"+suffix][0]);
  sigma_gr -> SetLineColor(kBlack);
  sigma_gr -> SetLineWidth(2);
  sigma_gr -> SetMarkerColor(kBlack);
  sigma_gr -> SetMarkerSize(0.8);
  sigma_gr -> Draw("ezpsame");

  double particle_mass = mass_muon;
  double segment_size = 14.;
  TF1 *HL_four_params = new TF1("HL_four_params", HL_function, 0., 2500., 6);
  HL_four_params -> SetParameters(0.022, 9.078, 5., 0.038, particle_mass, segment_size + 0.);
  /*
  HL_four_params -> SetParLimits(0, 0., 0.3);
  HL_four_params -> SetParLimits(1, 0., 20.);
  HL_four_params -> SetParLimits(2, 0.0005, 0.05);
  */
  HL_four_params -> SetParLimits(3, 0., 1.);
  HL_four_params -> FixParameter(4, particle_mass);
  HL_four_params -> FixParameter(5, segment_size + 0.);

  TF1 *HL_three_params = new TF1("HL_three_params", HL_function, 0., 2500., 6);
  HL_three_params -> SetParameters(0.022, 9.078, 5., 0.038, particle_mass, segment_size + 0.);
  //HL_three_params -> SetParLimits(0, 0., 0.3);
  HL_three_params -> SetParLimits(1, 0., 20.);
  HL_three_params -> SetParLimits(2, 0.0005,10.);
  HL_three_params -> FixParameter(3, 0.038);
  HL_three_params -> FixParameter(4, particle_mass);
  HL_three_params -> FixParameter(5, segment_size + 0.);

  HL_four_params -> SetParNames("kappa_a", "kappa_c", "sigma_res", "epsilon");
  HL_three_params -> SetParNames("kappa_a", "kappa_c", "sigma_res");

  sigma_gr -> Fit(HL_four_params, "RBOSQN", "", 0., 2500.);
  sigma_gr -> Fit(HL_three_params, "RBOSQN", "", 0., 2500.);
  double par_HL_four[6], par_HL_three[6], par_err_HL_four[6], par_err_HL_three[6], par_HL_mb[6], par_HL_s2[6];
  HL_four_params -> GetParameters(par_HL_four);
  HL_three_params -> GetParameters(par_HL_three);
  for(int i_par = 0; i_par < 6; i_par++){
    par_err_HL_four[i_par] = HL_four_params -> GetParError(i_par);
    par_err_HL_three[i_par] = HL_three_params -> GetParError(i_par);
  }
  HL_four_params -> SetLineColor(kRed);
  HL_four_params -> SetLineStyle(7);
  HL_four_params -> SetLineWidth(3);
  //HL_four_params -> Draw("lsame");
  HL_three_params -> SetLineColor(kBlue);
  HL_three_params -> SetLineStyle(4);
  HL_three_params -> SetLineWidth(2);
  HL_three_params -> Draw("lsame");

  TF1 *HL_mb = new TF1("HL_mb", HL_function, 0., 2500., 6);
  HL_mb -> SetParameters(0.105, 11.004, par_HL_three[2], 0.038, particle_mass, segment_size + 0.);
  HL_mb -> GetParameters(par_HL_mb);
  HL_mb -> SetLineColor(kGreen);
  HL_mb -> Draw("lsame");

  /*
  TF1 *HL_mb_1 = new TF1("HL_mb_1", HL_function, 0., 2500., 6);
  HL_mb_1 -> SetParameters(10., 11.004, 0.005, 0.038, particle_mass, segment_size + 0.);
  //HL_mb_1 -> GetParameters(par_HL_mb);
  HL_mb_1 -> SetLineColor(kPink);
  HL_mb_1 -> SetLineStyle(7);
  HL_mb_1 -> Draw("lsame");
  */
  
  TF1 *HL_s2 = new TF1("HL_s2", HL_function, 0., 2500., 6);
  HL_s2 -> SetParameters(0., 13.6, par_HL_three[2], 0.038, particle_mass, segment_size + 0.);
  HL_s2 -> GetParameters(par_HL_s2);
  HL_s2 -> SetLineColor(kOrange);
  HL_s2 -> SetLineStyle(9);
  HL_s2 -> Draw("lsame");

  //TString HL_four_param_str =  Form("#kappa_{a} = %.3e GeV^{2}#timesMeV, #kappa_{c} = %.3e MeV, #sigma_{res} = %.3e rad., #varepsilon = %.1e", par_HL_four[0], par_HL_four[1], par_HL_four[2], par_HL_four[3]);
  TString HL_mb_str = Form("MicroBooNE's fit: #kappa_{a} = %.3f [GeV^{2}#timesMeV], #kappa_{c} = %.3f [MeV], #sigma_{res} = %.3f [mrad.]", par_HL_mb[0], par_HL_mb[1], par_HL_mb[2]);
  TString HL_s2_str = Form("S_{2} = 13.6 MeV: #kappa_{a} = %.1f [GeV^{2}#timesMeV], #kappa_{c} = %.1f [MeV], #sigma_{res} = %.3f [mrad.]", par_HL_s2[0], par_HL_s2[1], par_HL_s2[2]);
  
  TLegend *l = new TLegend(0.30, 0.65, 0.90, 0.90);
  l -> AddEntry(sigma_gr, "Data points", "lp");
  //l -> AddEntry(HL_four_params, HL_four_param_str, "l");
  l -> AddEntry(HL_three_params, "This fit", "l");
  l -> AddEntry(HL_three_params, Form("#kappa_{a} = %.3f #pm %.3f [GeV^{2}#timesMeV]", HL_three_params -> GetParameter(0), HL_three_params -> GetParError(0)), "");
  l -> AddEntry(HL_three_params, Form("#kappa_{c} = %.3f #pm %.3f [MeV]", HL_three_params -> GetParameter(1), HL_three_params -> GetParError(1)), "");
  l -> AddEntry(HL_three_params, Form("#sigma_{res} = %.3f #pm %.3f [mrad.]", HL_three_params -> GetParameter(2), HL_three_params -> GetParError(2)), "");
  //l -> AddEntry(HL_three_params, HL_three_param_str, "l");
  l -> AddEntry(HL_mb, HL_mb_str, "l");
  l -> AddEntry(HL_s2, HL_s2_str, "l");
  l -> Draw("same");
  
  TString particle_label_str = "";
  if(particle == "muon") particle_label_str = "Cathode Passing Stopping Tracks";
  if(particle == "proton") particle_label_str = "Stopping Proton Candidates";
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
  TString outfile_str = output_plot_dir + "/MC_2025A_CV/mcs/HL_fit/p_vs_mcs_sigma_" + particle + "_" + suffix + ".pdf";
  if(isdata) outfile_str = output_plot_dir + "/run_" + run_str + "/mcs/HL_fit/p_vs_mcs_sigma_" + particle + "_" + suffix + ".pdf";

  TString outfile_dir = gSystem->DirName(outfile_str);
  if (gSystem->AccessPathName(outfile_dir)) {
    std::cout << "Directory does not exist, creating: " << outfile_dir << std::endl;
    gSystem->mkdir(outfile_dir, kTRUE);
  }

  c -> SaveAs(outfile_str);
  c -> Close();
}

void run_mcs_fit(int run_num = 0){

  if(run_num != 0){
    isdata = true;
    run_str = TString::Format("%d", run_num);
  }
  
  setTDRStyle();
  double p_binning_muon[] = {200.,
    240., 280., 320., 360., 400.,
    440., 480., 520., 560., 600.,
    640., 680., 720., 760., 800.,
    840., 880., 920., 960., 1000.,
    1040., 1080.
  };

  int rebin_2d_angle_muon[] = {
    20, 10, 10, 10, 10,
    10,	10, 10,	10, 10,
    10, 10, 10, 10, 10,
    10, 10, 10, 10, 10,
    10, 10, 
  };

  TString filename = "output_mcs_loop_run_" + run_str + ".root";
  if(!isdata){
    filename = "output_mcs_loop_mc.root";
    //filename = "output_mcs_loop_mc_p_after.root";
    run_str = "MC";
  }

  TString suffixes[] = {"_theta_xz", "_theta_yz", "_theta_2d"};
  TString suffix_leg[] = {"#theta_{XZ}", "#theta_{YZ}", "#theta_{2D}"};
  for(int i = 0; i < 3; i++){
    TString this_suffix = suffixes[i];
    TString this_suffix_leg = suffix_leg[i];
    Fit_p_vs_mcs_plots(filename, this_suffix, this_suffix_leg, "muon", -300., 300., p_binning_muon, 23, rebin_2d_angle_muon);
  }

}
