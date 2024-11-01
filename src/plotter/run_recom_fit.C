#include "canvas_margin.h"
#include "mylib.h"
#include "LanGausFit.h"
#include "TRandom3.h"
#include "BetheBloch.h"
#include <iostream>

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

TString run_str = "Run14480";

void Write_1D_hist(TH1D *in, TString outname, TString particle, TString latex_str, TString title_x, TString title_y, double x_min, double x_max, int rebin, bool do_langau_fit){

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
    fitting_range[0] = 0.;
    fitting_range[1] = 5000.;
    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
    sv[0] = 20.;
    sv[1] = max_x;
    sv[2] = in -> Integral() * 0.05 * bin_width;
    sv[3] = 70.;
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

    TLegend *l = new TLegend(0.20, 0.40, 0.40, 0.80);
    l -> AddEntry(in, Form("dQ/dx (%.0f)", in -> Integral()), "pl");
    l -> AddEntry(this_Langau, Form("#sigma_{Landau} : %.2f #pm %.2f", this_Landau_sigma, this_Landau_sigma_err), "l");
    l -> AddEntry(in, Form("MPV : %.2f #pm %.2f", this_MPV, this_MPV_err), "");
    l -> AddEntry(in, Form("#sigma_{Gaus} : %.2f #pm %.2f", this_Gaus_sigma, this_Gaus_sigma_err), "");
    l -> AddEntry(in, Form("Par2 : %.2f #pm %.2f", this_par2, this_par2_err), "");
    l -> AddEntry(in, Form("#chi^{2} / ndf : %.2f", chisqr / ndf), "");
    l -> Draw("same");

    MPV_vec_map[particle].push_back(this_MPV);
    MPV_err_vec_map[particle].push_back(this_MPV_err);
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
  
    dEdx_vec_map[particle].push_back(this_mean);
    dEdx_err_vec_map[particle].push_back(this_stddev);
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
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND Simulation} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_particle.DrawLatex(0.95, 0.96, particle_label_str);
  latex_method.DrawLatex(0.18, 0.87, latex_str);

  TString output_plot_dir = getenv("PLOT_PATH");
  c -> SaveAs(output_plot_dir + "/" + run_str + "/" + outname + ".pdf");
  c -> Close();

}

void Fit_rr_vs_pitch_plots(TString input_file_name, int rebin_x, int rebin_y, double rr_low, double rr_high){

  TString this_id = "id";

  TString input_file_dir = getenv("OUTPUTROOT_PATH");
  TFile *f = new TFile(input_file_dir + "/" + input_file_name);
  TH2D *hist_2D = (TH2D*)gDirectory -> Get("rr_vs_pitch_trklen_60cm_passing_cathode_coszx");

  hist_2D -> RebinX(rebin_x);
  hist_2D -> RebinY(rebin_y);

  int N_binsX = hist_2D -> GetNbinsX();
  int N_binsY = hist_2D -> GetNbinsY();

  for(int i = 1; i < N_binsX + 1; i++){
    TString i_str = Form("%d", i);
    double this_rr = hist_2D -> GetXaxis() -> GetBinCenter(i);
    if(this_rr > rr_high || this_rr < rr_low) continue;
    double this_rr_err = 0.5 * hist_2D -> GetXaxis() -> GetBinWidth(i);
    TString rr_str = Form("rr%.0fto%.0fcm", this_rr - this_rr_err, this_rr + this_rr_err);
    if(this_rr + this_rr_err < 10.){
      rr_str = Form("rr0%.0fto0%.0fcm", this_rr - this_rr_err, this_rr + this_rr_err);
    }
    else if(this_rr - this_rr_err < 10.){
      rr_str = Form("rr0%.0fto%.0fcm", this_rr - this_rr_err, this_rr + this_rr_err);
    }
    TString rr_latex = Form("Resdual range : %.2f - %.2f cm", this_rr -this_rr_err, this_rr + this_rr_err);
    TString this_hist_name = "Pitch_" + rr_str;

    TH1D * this_hist_1D = new TH1D(this_hist_name, this_hist_name, N_binsY, 0., 2.);

    for(int j = 1; j < N_binsY + 1; j++){
      double this_content = hist_2D -> GetBinContent(i, j);
      double this_content_err = hist_2D -> GetBinError(i, j);
      this_hist_1D -> SetBinContent(j, this_content);
      this_hist_1D -> SetBinError(j, this_content_err);
    }

    double max_y = this_hist_1D -> GetMaximum();

    TCanvas *c = new TCanvas("", "", 800, 600);
    canvas_margin(c);
    gStyle -> SetOptStat(1111);

    TH1D * template_h = new TH1D("", "", 1., 0., 2.);
    template_h -> SetStats(0);
    template_h -> GetYaxis() -> SetRangeUser(0., max_y * 1.5);
    template_h -> GetXaxis() -> SetTitle("Pitch [cm]");
    template_h -> GetXaxis() -> SetTitleSize(0.037);
    template_h -> GetXaxis() -> SetTitleOffset(1.4);
    template_h -> GetXaxis() -> SetLabelSize(0.035);
    template_h -> GetYaxis() -> SetTitle("Events");
    template_h -> GetYaxis() -> SetTitleSize(0.05);
    template_h -> GetYaxis() -> SetLabelSize(0.035);
    template_h -> Draw();

    this_hist_1D -> SetMarkerColor(kBlack);
    this_hist_1D -> SetMarkerStyle(32);
    this_hist_1D -> SetMarkerSize(0.7);
    this_hist_1D -> SetLineColor(kBlack);
    this_hist_1D -> SetLineWidth(1);
    this_hist_1D -> Draw("epsame");

    TH1D *this_hist_1D_clone = (TH1D*)this_hist_1D -> Clone();
    double this_1D_mean = this_hist_1D -> GetMean();
    TLine *l_mean = new TLine(this_1D_mean, 0., this_1D_mean, max_y * 1.5);
    l_mean -> SetLineStyle(7);
    l_mean -> SetLineColor(kRed);
    l_mean -> Draw("lsame");

    TLegend *l = new TLegend(0.50, 0.40, 0.92, 0.85);
    l -> AddEntry(this_hist_1D, "Pitch [cm]", "pl");
    l -> AddEntry(l_mean, Form("Mean : %.2f", this_1D_mean), "l");
    l -> Draw("same");

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
    latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND Simulation} #font[42]{#it{#scale[0.8]{Preliminary}}}");
    latex_particle.DrawLatex(0.95, 0.96, "Cathode Passing Stopping Tracks");
    latex_method.DrawLatex(0.18, 0.87, rr_latex);

    TString output_plot_dir = getenv("PLOT_PATH");
    output_plot_dir = output_plot_dir + "/recom_fit/1D/pitch/";
    c -> SaveAs(output_plot_dir + "/" + run_str + "/" + this_hist_name + ".pdf");

    c -> Close();

    fitting_results[this_id + "_pitch_mean"].push_back(this_1D_mean);
  }
}


void Fit_rr_vs_dqdx_plots(TString input_file_name, int rebin_x, int rebin_y, double rr_low, double rr_high){

  TString this_id = "id";

  TString input_file_dir = getenv("OUTPUTROOT_PATH");
  TFile *f = new TFile(input_file_dir + "/" + input_file_name);
  TH2D *hist_2D = (TH2D*)gDirectory -> Get("rr_vs_corr_dqdx_trklen_60cm_passing_cathode_coszx");

  hist_2D -> RebinX(rebin_x);
  hist_2D -> RebinY(rebin_y);

  int N_binsX = hist_2D -> GetNbinsX();
  int N_binsY = hist_2D -> GetNbinsY();

  for(int i = 1; i < N_binsX + 1; i++){
    TString i_str = Form("%d", i);
    double this_rr = hist_2D -> GetXaxis() -> GetBinCenter(i);
    if(this_rr > rr_high || this_rr < rr_low) continue;
    double this_rr_err = 0.5 * hist_2D -> GetXaxis() -> GetBinWidth(i);
    TString rr_str = Form("rr%.0fto%.0fcm", this_rr - this_rr_err, this_rr + this_rr_err);
    if(this_rr + this_rr_err < 10.){
      rr_str = Form("rr0%.0fto0%.0fcm", this_rr - this_rr_err, this_rr + this_rr_err);
    }
    else if(this_rr - this_rr_err < 10.){
      rr_str = Form("rr0%.0fto%.0fcm", this_rr - this_rr_err, this_rr + this_rr_err);
    }
    TString rr_latex = Form("Resdual range : %.2f - %.2f cm", this_rr -this_rr_err, this_rr + this_rr_err);
    TString this_hist_name = "Corr_dQdx_" + rr_str;

    TH1D * this_hist_1D = new TH1D(this_hist_name, this_hist_name, N_binsY, 0., 3000.);

    for(int j = 1; j < N_binsY + 1; j++){
      double this_content = hist_2D -> GetBinContent(i, j);
      double this_content_err = hist_2D -> GetBinError(i, j);
      //cout << "[Fit_rr_vs_dqdx_plots] this_content : " << this_content << endl;
      this_hist_1D -> SetBinContent(j, this_content);
      this_hist_1D -> SetBinError(j, this_content_err);
    }

    double max_y = this_hist_1D -> GetMaximum();

    TCanvas *c = new TCanvas("", "", 800, 600);
    canvas_margin(c);
    gStyle -> SetOptStat(1111);

    TH1D * template_h = new TH1D("", "", 1., 0., 3000.);
    template_h -> SetStats(0);
    template_h -> GetYaxis() -> SetRangeUser(0., max_y * 1.5);
    template_h -> GetXaxis() -> SetTitle("dQ/dx [ADC/cm]");
    template_h -> GetXaxis() -> SetTitleSize(0.037);
    template_h -> GetXaxis() -> SetTitleOffset(1.4);
    template_h -> GetXaxis() -> SetLabelSize(0.035);
    template_h -> GetYaxis() -> SetTitle("Events");
    template_h -> GetYaxis() -> SetTitleSize(0.05);
    template_h -> GetYaxis() -> SetLabelSize(0.035);
    template_h -> Draw();

    this_hist_1D -> SetMarkerColor(kBlack);
    this_hist_1D -> SetMarkerStyle(32);
    this_hist_1D -> SetMarkerSize(0.7);
    this_hist_1D -> SetLineColor(kBlack);
    this_hist_1D -> SetLineWidth(1);
    this_hist_1D -> Draw("epsame");

    TH1D *this_hist_1D_clone = (TH1D*)this_hist_1D -> Clone();

    double max_x = this_hist_1D -> GetBinCenter(this_hist_1D -> GetMaximumBin());
    Double_t fitting_range[2];
    fitting_range[0] = 0.;
    fitting_range[1] = 3000.;
    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
    sv[0] = 20.;
    sv[1] = max_x;
    //sv[2] = this_hist_1D -> Integral() * 0.05;
    sv[2] = this_hist_1D -> Integral();
    sv[3] = 70.;
    for(int j=0; j<4; ++j){
      pllo[j] = 0.01*sv[j];
      plhi[j] = 100*sv[j];
    }

    Double_t chisqr;
    Int_t    ndf;
    Int_t    status;
    TF1 *this_Langau_fit = langaufit(this_hist_1D_clone, fitting_range, sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status, "Langau_this_Langau" + this_hist_name);

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

    this_hist_1D -> Draw("epsame");

    /*
    TLine *l_MPV = new TLine(this_MPV, 0., this_MPV, max_y * 1.5);
    l_mean -> SetLineStyle(7);
    l_mean -> SetLineColor(kRed);
    l_mean -> Draw("lsame");
    */
    TLegend *l = new TLegend(0.20, 0.40, 0.40, 0.80);
    l -> AddEntry(this_hist_1D, Form("dQ/dx (%.0f)", this_hist_1D -> Integral()), "pl");
    l -> AddEntry(this_Langau, Form("#sigma_{Landau} : %.2f #pm %.2f", this_Landau_sigma, this_Landau_sigma_err), "l");
    l -> AddEntry(this_hist_1D, Form("MPV : %.2f #pm %.2f", this_MPV, this_MPV_err), "");
    l -> AddEntry(this_hist_1D, Form("#sigma_{Gaus} : %.2f #pm %.2f", this_Gaus_sigma, this_Gaus_sigma_err), "");
    l -> AddEntry(this_hist_1D, Form("Par2 : %.2f #pm %.2f", this_par2, this_par2_err), "");
    l -> AddEntry(this_hist_1D, Form("#chi^{2} / ndf : %.2f", chisqr / ndf), "");
    //l -> AddEntry(l_mean, Form("Mean : %.2f", this_1D_mean), "l");
    l -> Draw("same");

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
    latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND Simulation} #font[42]{#it{#scale[0.8]{Preliminary}}}");
    latex_particle.DrawLatex(0.95, 0.96, "Cathode Passing Stopping Tracks");
    latex_method.DrawLatex(0.18, 0.87, rr_latex);

    gPad->RedrawAxis();

    TString output_plot_dir = getenv("PLOT_PATH");
    output_plot_dir = output_plot_dir + "/recom_fit/1D/dqdx/";
    c -> SaveAs(output_plot_dir + "/" + run_str + "/" + this_hist_name + ".pdf");

    c -> Close();

    fitting_results[this_id + "_rr"].push_back(this_rr);
    fitting_results[this_id + "_rr_err"].push_back(this_rr_err);
    fitting_results[this_id + "_MPV"].push_back(this_MPV);
    fitting_results[this_id + "_MPV_err"].push_back(this_MPV_err);
    fitting_results[this_id + "_sigma_gaus"].push_back(this_Gaus_sigma);
    fitting_results[this_id + "_sigma_gaus_err"].push_back(this_Gaus_sigma_err);
    fitting_results[this_id + "_sigma_Landau"].push_back(this_Landau_sigma);
    fitting_results[this_id + "_sigma_Landau_err"].push_back(this_Landau_sigma_err);
  }
}

void Fit_dEdx_MPV_vs_dqdx_plots(TString input_file_name, TString histname, TString particle, int rebin_y, double dEdx_low, double dEdx_high, double dEdx_MPV_binning[], int num_NbinsX, int rebin_dqdx[]){

  //const Int_t num_NbinsX = sizeof(dEdx_MPV_binning) / sizeof(dEdx_MPV_binning[0]) - 1;
  
  TString this_id = "id";

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

    TString dEdx_str  = Form("dEdx_MPV_%.1fto%.1f_MeVcm", this_dEdx_MPV_low, this_dEdx_MPV_high);
    TString dEdx_latex_str = Form("dE/dx MPV %.1f to %.1f MeV/cm", this_dEdx_MPV_low, this_dEdx_MPV_high); 
    TString this_1D_X_hist_name = "dEdx_MPV_" + dEdx_str;
    TString this_1D_Y_hist_name = "dQdx_" + dEdx_str;

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

    Write_1D_hist(this_1D_X, "/recom_fit/1D/dedx/" + particle + "/" + this_1D_X_hist_name, particle, dEdx_latex_str, "dE/dx MPV [MeV/cm]", "Events", this_dEdx_MPV_low, this_dEdx_MPV_high, 1., false);
    Write_1D_hist(this_1D_Y, "/recom_fit/1D/dqdx/" + particle + "/" + this_1D_Y_hist_name, particle, dEdx_latex_str, "dQ/dx [ADC/cm]", "Events", 0., 5000., rebin_dqdx[i - 1], true);
  }

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  c -> SetLogz();

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
  
  TGraphErrors * this_gr = new TGraphErrors(MPV_vec_map[particle].size(), &dEdx_vec_map[particle][0], &MPV_vec_map[particle][0], &dEdx_err_vec_map[particle][0], &MPV_err_vec_map[particle][0]);
  this_gr -> SetLineColor(kRed);
  this_gr -> SetLineWidth(2);
  this_gr -> SetMarkerColor(kRed);
  this_gr -> SetMarkerStyle(32);
  this_gr -> SetMarkerSize(0.7);
  this_gr -> Draw("epsame");

  for(unsigned int i = 0; i < MPV_vec_map[particle].size(); i++){
    cout << i << ", MPV_vec_map[particle] : " << MPV_vec_map[particle].at(i) << ", dEdx_vec : " << dEdx_vec_map[particle].at(i) << endl;
  }

  double fit_x_max = 4.5;
  if(particle == "proton") fit_x_max = 9.0;  
  TF1 * f_mod_box = new TF1("f_mod_box", "(294.49153 * [0] / [1]) * log(1.4388489 * [1] * x + [2])", 1.6, fit_x_max);
  f_mod_box -> SetParameters(2.00, 0.212, 0.93);
  f_mod_box -> FixParameter(0, 2.00);
  f_mod_box -> SetLineColor(kBlack);
  f_mod_box -> SetLineWidth(3);
  f_mod_box -> SetLineStyle(7);
  this_gr -> Fit("f_mod_box", "RNS");
  f_mod_box -> Draw("lsame");

  TLegend *l = new TLegend(0.5, 0.25, 0.85, 0.50);
  l -> AddEntry(this_gr, "Data points", "lp");
  l -> AddEntry(f_mod_box, "Fitting result", "l");
  l -> AddEntry(f_mod_box, Form("C_{cal.} = %.2f  #pm %.3f #times 10^{-2} [ADC/electrons]", f_mod_box -> GetParameter(0), f_mod_box -> GetParError(0)), "");
  l -> AddEntry(f_mod_box, Form("#beta' = %.3f #pm %.3f [(kV/cm)(g/cm^{2})/MeV]", f_mod_box -> GetParameter(1), f_mod_box -> GetParError(1)), "");
  l -> AddEntry(f_mod_box, Form("#alpha = %.2f #pm %.3f",f_mod_box -> GetParameter(2), f_mod_box -> GetParError(2)), "");
  l -> Draw("same");

  TString output_plot_dir = getenv("PLOT_PATH");
  output_plot_dir = output_plot_dir + "/" + run_str + "/" + "/recom_fit/2D/";
  c -> SaveAs(output_plot_dir + "dEdx_MPV_vs_corr_dqdx_" + particle + ".pdf");
  
  c -> Close();

}

void fit_modified_box(TString id, double x_min, double x_max, double y_min, double y_max){

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D * template_h = new TH1D("", "", 1., x_min, x_max);
  template_h -> SetStats(0);
  template_h -> GetYaxis() -> SetRangeUser(y_min, y_max);
  template_h -> GetXaxis() -> SetTitle("dE/dx [MeV/cm]");
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("dQ/dx [ADC/cm]");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> Draw();

  TGraphErrors * muon_gr = new TGraphErrors(MPV_vec_map["muon"].size(), &dEdx_vec_map["muon"][0], &MPV_vec_map["muon"][0], &dEdx_err_vec_map["muon"][0], &MPV_err_vec_map["muon"][0]);
  muon_gr -> SetLineColor(kRed);
  muon_gr -> SetLineWidth(2);
  muon_gr -> SetMarkerColor(kRed);
  muon_gr -> SetMarkerStyle(32);
  muon_gr -> SetMarkerSize(0.7);
  muon_gr -> Draw("epsame");

  TGraphErrors * proton_gr = new TGraphErrors(MPV_vec_map["proton"].size(), &dEdx_vec_map["proton"][0], &MPV_vec_map["proton"][0], &dEdx_err_vec_map["proton"][0], &MPV_err_vec_map["proton"][0]);
  proton_gr -> SetLineColor(kBlue);
  proton_gr -> SetLineWidth(2);
  proton_gr -> SetMarkerColor(kBlue);
  proton_gr -> SetMarkerStyle(32);
  proton_gr -> SetMarkerSize(0.7);
  proton_gr -> Draw("epsame");

  TGraphErrors * combined_gr = new TGraphErrors(MPV_vec.size(), &dEdx_vec[0], &MPV_vec[0], &dEdx_err_vec[0], &MPV_err_vec[0]);
  TF1 * f_mod_box = new TF1("f_mod_box", "(294.49153 * [0] / [1]) * log(1.4388489 * [1] * x + [2])", 1.6, 9.);
  f_mod_box -> SetParameters(2.00, 0.212, 0.93);
  //f_mod_box -> FixParameter(0, 2.00);
  f_mod_box -> SetLineColor(kBlack);
  f_mod_box -> SetLineWidth(2);
  f_mod_box -> SetLineStyle(7);
  combined_gr -> Fit("f_mod_box", "RNS");
  f_mod_box -> Draw("lsame");

  TLegend *l = new TLegend(0.5, 0.25, 0.85, 0.50);
  l -> AddEntry(muon_gr, "Muon points", "lp");
  l -> AddEntry(proton_gr, "Proton points", "lp");
  l -> AddEntry(f_mod_box, "Fitting result", "l");
  l -> AddEntry(f_mod_box, Form("C_{cal.} = %.2f  #pm %.3f #times 10^{-2} [ADC/electrons]", f_mod_box -> GetParameter(0), f_mod_box -> GetParError(0)), "");
  l -> AddEntry(f_mod_box, Form("#beta' = %.3f #pm %.3f [(kV/cm)(g/cm^{2})/MeV]", f_mod_box -> GetParameter(1), f_mod_box -> GetParError(1)), "");
  l -> AddEntry(f_mod_box, Form("#alpha = %.2f #pm %.3f",f_mod_box -> GetParameter(2), f_mod_box -> GetParError(2)), "");
  l -> Draw("same");

  TString output_plot_dir = getenv("PLOT_PATH");
  output_plot_dir = output_plot_dir + "/" + run_str + "/" + "/recom_fit/2D/";
  c -> SaveAs(output_plot_dir + "dEdx_MPV_vs_corr_dqdx_combined.pdf");
}



void run_recom_fit(){

  setTDRStyle();
  //Fit_rr_vs_pitch_plots("output_recom.root", 2, 4, 2., 80.);
  //Fit_rr_vs_dqdx_plots("output_recom.root", 2, 20, 2., 80.);
  double dEdx_MPV_binning_muon[] = {0.,
				    1.6, 1.7, 1.8, 1.9, 2.0, 2.1,
				    2.2, 2.3, 2.4, 2.5, 2.7, 3.0,
				    3.5, 4.5, 30.};
  int rebin_dqdx_muon[] = {1, 1, 2, 2, 2,
			   4, 4, 4, 4, 4,
			   4, 4, 4};


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

  Fit_dEdx_MPV_vs_dqdx_plots("output_recom_run14480.root", "dEdx_MPV_vs_corr_dqdx_trklen_60cm_passing_cathode_coszx", "muon", 5, 1.5, 4.7, dEdx_MPV_binning_muon, 15, rebin_dqdx_muon);
  //Fit_dEdx_MPV_vs_dqdx_plots("output_recom_muscore40.root", "dEdx_MPV_vs_corr_dqdx_proton", "proton", 5, 1.5, 10.0, dEdx_MPV_binning_proton, 20, rebin_dqdx_proton);

  fit_modified_box("", 1.5, 10.0, 0., 5000.);
}
