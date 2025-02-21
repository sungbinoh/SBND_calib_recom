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

// == EBM par
double alpha_ebm = 0.904;
double beta_90 = 0.204;
double R_ebm = 1.25;

double B_phi(double phi = 90.){
  double phi_rad = phi * TMath::Pi() / 180.;
  double out = beta_90 / (sqrt(pow(sin(phi_rad), 2.) + pow(cos(phi_rad) / R_ebm, 2.)));
  return out;
}

void Write_1D_hist(TH1D *in, TString outname, TString suffix, TString particle, TString latex_str, TString plane, TString title_x, TString title_y, double x_min, double x_max, int rebin, bool do_langau_fit){

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

  if(do_langau_fit){

    TH1D *in_clone = (TH1D*)in -> Clone();
    double max_x = in -> GetBinCenter(in -> GetMaximumBin());
    double bin_width = in -> GetBinWidth(1);
    Double_t fitting_range[2];
    fitting_range[0] = 800.;
    fitting_range[1] = 3000.;
    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
    sv[0] = 40.;
    sv[1] = max_x;
    sv[2] = in -> Integral() * 0.05 * bin_width;
    sv[3] = 200.;
    for(int j=0; j<4; ++j){
      pllo[j] = 0.01*sv[j];
      plhi[j] = 100*sv[j];
    }

    Double_t chisqr;
    Int_t    ndf;
    Int_t    status;
    TF1 *this_Langau_fit = langaufit(in_clone, fitting_range, sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status, "Langau_this_Langau" + outname);

    TF1 *this_Langau =  new TF1("this_Langau", langaufun, fitting_range[0], fitting_range[1], 4);
    this_Langau -> SetParameters(this_Langau_fit -> GetParameters());
    this_Langau -> SetNpx(1000);
    this_Langau -> SetLineColor(kRed);
    this_Langau -> SetLineWidth(2);
    this_Langau -> Draw("lsame");

    double this_Landau_sigma = this_Langau_fit ->GetParameter(0);
    double this_Landau_sigma_err =this_Langau_fit -> GetParError(0);
    double this_MPV = this_Langau_fit -> GetParameter(1);
    double this_MPV_err = this_Langau_fit -> GetParError(1);
    double this_par2 = this_Langau_fit -> GetParameter(2);
    double this_par2_err =this_Langau_fit -> GetParError(2);
    double this_Gaus_sigma= this_Langau_fit -> GetParameter(3);
    double this_Gaus_sigma_err= this_Langau_fit -> GetParError(3);

    in -> Draw("epsame");

    TLegend *l = new TLegend(0.60, 0.40, 0.80, 0.80);
    l -> AddEntry(in, Form("dQ/dx (%.0f)", in -> Integral()), "pl");
    l -> AddEntry(this_Langau, Form("#sigma_{Landau} : %.2f #pm %.2f", this_Landau_sigma, this_Landau_sigma_err), "l");
    l -> AddEntry(in, Form("MPV : %.2f #pm %.2f", this_MPV, this_MPV_err), "");
    l -> AddEntry(in, Form("#sigma_{Gaus} : %.2f #pm %.2f", this_Gaus_sigma, this_Gaus_sigma_err), "");
    l -> AddEntry(in, Form("Par2 : %.2f #pm %.2f", this_par2, this_par2_err), "");
    l -> AddEntry(in, Form("#chi^{2} / ndf : %.2f", chisqr / ndf), "");
    l -> Draw("same");

    MPV_vec_map[particle+plane+suffix].push_back(this_MPV);
    MPV_err_vec_map[particle+plane+suffix].push_back(this_MPV_err);
    MPV_vec.push_back(this_MPV);
    MPV_err_vec.push_back(this_MPV_err);
  }
  else{

    double this_mean = in -> GetMean();
    double this_stddev = in -> GetStdDev();

    TLegend *l = new TLegend(0.70, 0.70, 0.850, 0.85);
    l -> AddEntry(in, Form("Mean : %.2f", this_mean), "pl");
    l -> AddEntry(in, Form("StdDev : %.2f", this_stddev), "");
    l -> Draw("same");
  
    dEdx_vec_map[particle+plane+suffix].push_back(this_mean);
    dEdx_err_vec_map[particle+plane+suffix].push_back(this_stddev);
    dEdx_vec.push_back(this_mean);
    dEdx_err_vec.push_back(this_stddev);
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
  TString outfile_str = output_plot_dir + "/MC_2024B_CV/recom/" + outname + ".pdf";
  if(isdata) outfile_str = output_plot_dir + "/Run" + run_str + "/" + outname + ".pdf";
  c -> SaveAs(outfile_str);
  c -> Close();

}

void Fit_dEdx_MPV_vs_dqdx_plots(TString input_file_name, TString suffix, double phi, TString plane, TString particle, int rebin_y, double dEdx_low, double dEdx_high, double dEdx_MPV_binning[], int num_NbinsX, int rebin_dqdx[]){

  //const Int_t num_NbinsX = sizeof(dEdx_MPV_binning) / sizeof(dEdx_MPV_binning[0]) - 1;
  
  TString this_id = "id";
  TString histname = "dEdx_MPV_vs_corr_dqdx_" + plane + "_trklen_60cm_passing_cathode_coszx" + suffix;
  
  TString input_file_dir = getenv("OUTPUTROOT_PATH");
  TFile *f = new TFile(input_file_dir + "/" + input_file_name);
  TH2D *hist_2D = (TH2D*)gDirectory -> Get(histname);

  //hist_2D -> Rebin(num_NbinsX, "RebinX", dEdx_MPV_binning);
  hist_2D -> RebinY(rebin_y);

  int N_binsX = hist_2D -> GetNbinsX();
  int N_binsY = hist_2D -> GetNbinsY();

  TH1D *hist_1D_for_X = hist_2D -> ProjectionX("hist_1D_for_X");


  double hist_2D_x_low = 0.;
  double hist_2D_x_high = 30.;
  if(particle == "proton") hist_2D_x_high = 15.;
  double hist_2D_x_step = (hist_2D_x_high - hist_2D_x_low) / N_binsX;

  double dqdx_max = 3000.;
  if(particle == "proton") dqdx_max = 5000.;

  cout << "num_NbinsX : " << num_NbinsX << ", N_binsX : " << N_binsX << ", hist_2D_x_step : " << hist_2D_x_step << endl;

  for(int i = 1; i < num_NbinsX - 1; i++){
      
    double this_dEdx_MPV_low = dEdx_MPV_binning[i];
    double this_dEdx_MPV_high = dEdx_MPV_binning[i + 1];

    int x_index_low = 1 + (this_dEdx_MPV_low - hist_2D_x_low) / hist_2D_x_step;
    int x_index_high = 1 + (this_dEdx_MPV_high - hist_2D_x_low) / hist_2D_x_step;
    //cout << "this_dEdx_MPV_low : " << this_dEdx_MPV_low << ", this_dEdx_MPV_high : " << this_dEdx_MPV_high << ", x_index_low : " << x_index_low << ", x_index_high : " << x_index_high << endl;
    //cout << << endl;

    TString dEdx_str  = Form("dEdx_MPV_%.2fto%.2f_MeVcm", this_dEdx_MPV_low, this_dEdx_MPV_high);
    TString dEdx_latex_str = Form("dE/dx MPV %.2f to %.2f MeV/cm", this_dEdx_MPV_low, this_dEdx_MPV_high); 
    TString this_1D_X_hist_name = "dEdx_MPV_" + dEdx_str + suffix;
    TString this_1D_Y_hist_name = "dQdx_" + dEdx_str + suffix;

    TH1D * this_1D_X = new TH1D(this_1D_X_hist_name, this_1D_X_hist_name, x_index_high - x_index_low, this_dEdx_MPV_low, this_dEdx_MPV_high);
    TH1D * this_1D_Y = new TH1D(this_1D_Y_hist_name, this_1D_Y_hist_name, N_binsY, 0., dqdx_max);

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

    Write_1D_hist(this_1D_X, "/" + plane + "_" + this_1D_X_hist_name, suffix, particle, dEdx_latex_str, plane, "dE/dx MPV [MeV/cm]", "Events", this_dEdx_MPV_low, this_dEdx_MPV_high, 1., false);
    Write_1D_hist(this_1D_Y, "/" + plane + "_" + this_1D_Y_hist_name, suffix, particle, dEdx_latex_str, plane, "dQ/dx [ADC/cm]", "Events", 0., 5000., rebin_dqdx[i - 1], true);
  }

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  c -> SetLogz();
  c -> SetRightMargin(0.15);  
  
  double z_max = hist_2D -> GetMaximum();
  TH1D * template_h = new TH1D("", "", 1., dEdx_low, dEdx_high);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("dE/dx [MeV/cm]");
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("dQ/dx [ADC/cm]");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., dqdx_max);
  template_h -> Draw("colz");

  hist_2D -> Draw("colzsame");
  
  TGraphErrors * this_gr = new TGraphErrors(MPV_vec_map[particle+plane+suffix].size(), &dEdx_vec_map[particle+plane+suffix][0], &MPV_vec_map[particle+plane+suffix][0], &dEdx_err_vec_map[particle+plane+suffix][0], &MPV_err_vec_map[particle+plane+suffix][0]);
  this_gr -> SetLineColor(kBlack);
  this_gr -> SetLineWidth(2);
  this_gr -> SetMarkerColor(kBlack);
  this_gr -> SetMarkerSize(0.8);
  this_gr -> Draw("ezpsame");

  for(unsigned int i = 0; i < MPV_vec_map[particle+plane+suffix].size(); i++){
    cout << i << ", MPV_vec_map[particle] " << plane << " : " << MPV_vec_map[particle+plane+suffix].at(i) << ", dEdx_vec : " << dEdx_vec_map[particle+plane+suffix].at(i) << endl;
  }

  double fit_x_max = 2.0;
  if(particle == "proton") fit_x_max = 9.0;  
  TF1 * f_mod_box = new TF1("f_mod_box", "(294.49153 * [0] / [1]) * log(1.4388489 * [1] * x + [2])", 1.6, fit_x_max);
  double this_B_phi = B_phi(phi);
  f_mod_box -> SetParameters(2.00, this_B_phi, alpha_ebm);
  f_mod_box -> FixParameter(1, this_B_phi);
  f_mod_box -> FixParameter(2, alpha_ebm);
  //f_mod_box -> SetLineColor(kBlack);
  f_mod_box -> SetLineColor(kRed);
  f_mod_box -> SetLineWidth(3);
  f_mod_box -> SetLineStyle(7);
  this_gr -> Fit("f_mod_box", "RNS");
  f_mod_box -> Draw("lsame");

  this_gr -> Draw("ezpsame");
  
  TLegend *l = new TLegend(0.45, 0.65, 0.80, 0.90);
  l -> AddEntry(this_gr, "Data points", "lp");
  l -> AddEntry(f_mod_box, "Fitting result", "l");
  l -> AddEntry(f_mod_box, Form("C_{cal.} = %.3f #pm %.3f #times 10^{-2} [ADC/electrons]", f_mod_box -> GetParameter(0), f_mod_box -> GetParError(0)), "");
  l -> AddEntry(f_mod_box, Form("#beta' = %.3f #pm %.3f [(kV/cm)(g/cm^{3})/MeV]", f_mod_box -> GetParameter(1), f_mod_box -> GetParError(1)), "");
  l -> AddEntry(f_mod_box, Form("#alpha = %.2f #pm %.3f",f_mod_box -> GetParameter(2), f_mod_box -> GetParError(2)), "");
  l -> Draw("same");

  vec_c_cal.push_back(f_mod_box -> GetParameter(0));
  vec_c_cal_err.push_back(f_mod_box -> GetParError(0));
  
  TString particle_label_str = "";
  if(particle == "muon") particle_label_str = "Cathode Passing Stopping Tracks";
  if(particle == "proton") particle_label_str = "Stopping Proton Candidates";
  TLatex latex_SBND, latex_method;
  latex_SBND.SetNDC();
  latex_method.SetNDC();
  latex_method.SetTextAlign(31);
  latex_SBND.SetTextSize(0.03);
  latex_method.SetTextSize(0.03);
  if(isdata) latex_SBND.DrawLatex(0.16, 0.96, "#font[62]{SBND Data} Run " + run_str + ", " + plane);
  else latex_SBND.DrawLatex(0.16, 0.96, "#font[62]{SBND Simulation} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_method.DrawLatex(0.90, 0.96, particle_label_str);

  TString output_plot_dir = getenv("PLOT_PATH");
  TString outfile_str = output_plot_dir + "/MC_2024B_CV/recom/dEdx_MPV_vs_corr_dqdx_" + particle + "_" + plane + suffix + ".pdf";
  if(isdata) outfile_str = output_plot_dir + "/Run" + run_str + "/dEdx_MPV_vs_corr_dqdx_" + particle + "_" + plane + suffix + ".pdf";
  c -> SaveAs(outfile_str);
  
  c -> Close();

}

void fit_ebm(TString particle, TString plane, double dEdx_low, double dEdx_high){
  
  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D * template_h = new TH1D("", "", 1., dEdx_low, dEdx_high);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("dE/dx [MeV/cm]");
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("dQ/dx [ADC/cm]");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., 3000.);
  template_h -> Draw("colz");

  TString suffixes[] = {"_phi40to50", "_phi50to60", "_phi60to70", "_phi70to80", "_phi80to85", "_phi85to90"};
  int phi_color_arr[] = {632, 800, 401, 418, 600, 880};
  map<TString, TGraphErrors*> grs;
  
  for(int i = 0; i < 6; i++){
    TString suffix = suffixes[i];
    grs[suffixes[i]] = new TGraphErrors(MPV_vec_map[particle+plane+suffix].size(), &dEdx_vec_map[particle+plane+suffix][0], &MPV_vec_map[particle+plane+suffix][0], &dEdx_err_vec_map[particle+plane+suffix][0], &MPV_err_vec_map[particle+plane+suffix][0]);
    grs[suffixes[i]] -> SetLineColor(phi_color_arr[i]);
    grs[suffixes[i]] -> SetLineWidth(2);
    grs[suffixes[i]] -> SetMarkerColor(phi_color_arr[i]);
    grs[suffixes[i]] -> SetMarkerSize(0.8);
    grs[suffixes[i]] -> Draw("ezpsame");
    grs[suffixes[i]] -> SetName(plane + suffix);
    grs[suffixes[i]] -> Write();
  }

  TString particle_label_str = "";
  if(particle == "muon") particle_label_str = "Cathode Passing Stopping Tracks";
  if(particle == "proton") particle_label_str = "Stopping Proton Candidates";
  TLatex latex_SBND, latex_method;
  latex_SBND.SetNDC();
  latex_method.SetNDC();
  latex_method.SetTextAlign(31);
  latex_SBND.SetTextSize(0.03);
  latex_method.SetTextSize(0.03);
  if(isdata) latex_SBND.DrawLatex(0.16, 0.96, "#font[62]{SBND Data} Run " + run_str + ", " + plane);
  else latex_SBND.DrawLatex(0.16, 0.96, "#font[62]{SBND Simulation 2024B CV} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_method.DrawLatex(0.90, 0.96, particle_label_str);

  TString output_plot_dir = getenv("PLOT_PATH");
  TString outfile_str = output_plot_dir + "/MC_2024B_CV/recom/dEdx_MPV_vs_corr_dqdx_" + particle + "_" + plane + "_allphi.pdf";
  if(isdata) outfile_str = output_plot_dir + "/Run" + run_str + "/dEdx_MPV_vs_corr_dqdx_" + particle + "_" + plane + "_allphi.pdf";
  c -> SaveAs(outfile_str);

  c -> Close();  
}

void run_recom_ebm_fit(int run_num = 0){

  if(run_num != 0){
    isdata = true;
    run_str = TString::Format("%d", run_num);
  }
  
  setTDRStyle();
  double dEdx_MPV_binning_muon[] = {0.,
				    1.6, 1.7, 1.8, 1.9, 2.0,
				    2.1, 2.3, 2.6, 3.1, 4.0
				    };
  int rebin_dqdx_muon[] = {0,
			   2, 2, 5, 5, 5,
			   5, 5, 10, 10, 10,
			   5, 5, 5, 5, 5,
  };


  double dEdx_MPV_binning_proton[] = {0.,
				      2.0, 2.4, 2.8, 3.0, 3.2, 3.4,
				      3.6, 3.8, 4.0, 4.2, 4.4, 4.8,
				      5.2, 5.6, 6.0, 6.5, 7.0, 8.0,
				      9.0, 15.};
  int rebin_dqdx_proton[] = {4, 4, 4, 4, 4,
			     2, 2, 2, 2, 2,
			     2, 2, 2, 2, 2,
			     4, 4, 4, 4, 4,
			     4, 4, 4, 4, 4,
			     4, 4, 4, 4, 4};

  TString filename = "output_recom_" + run_str + ".root";
  if(!isdata){
    filename = "output_recom_2024B_GENIE_CV.root";
    run_str = "MC";
  }

  TString suffixes[] = {"", "_phi40to50", "_phi50to60", "_phi60to70", "_phi70to80", "_phi80to85", "_phi85to90"};
  double phis[] = {65., 45., 55., 65., 75., 82.5, 87.5};
  for(int i = 0; i < 7; i++){
    TString this_suffix = suffixes[i];
    Fit_dEdx_MPV_vs_dqdx_plots(filename, this_suffix, phis[i], "plane0", "muon", 5, 1.5, 4.5, dEdx_MPV_binning_muon, 11, rebin_dqdx_muon);
    Fit_dEdx_MPV_vs_dqdx_plots(filename, this_suffix, phis[i], "plane1", "muon", 5, 1.5, 4.5, dEdx_MPV_binning_muon, 11, rebin_dqdx_muon);
    Fit_dEdx_MPV_vs_dqdx_plots(filename, this_suffix, phis[i], "plane2", "muon", 5, 1.5, 4.5, dEdx_MPV_binning_muon, 11, rebin_dqdx_muon);
  }

  fit_ebm("muon", "plane0", 1.5, 4.5);

  TString output_file_dir = getenv("OUTPUTROOT_PATH");
  TString out_root_name = output_file_dir + "/recom_fit_ebm_" + run_str + ".root";
  TFile *outfile = new TFile(out_root_name, "RECREATE");
  outfile -> cd();

  fit_ebm("muon", "plane0", 1.5, 4.5);
  fit_ebm("muon", "plane1", 1.5, 4.5);
  fit_ebm("muon", "plane2", 1.5, 4.5);
}
