#include "canvas_margin.h"
#include "mylib.h"
#include "LanGausFit.h"
#include "TRandom3.h"

TRandom3 gRan(1800);
map<TString, vector<double>> fitting_results;

void Fit_1D_plots(TString input_file_name, int rebin_x, int rebin_y, double tdrift_low, double tdrift_high){

  TString this_id = "id";

  TString input_file_dir = getenv("OUTPUTROOT_PATH");
  TFile *f = new TFile(input_file_dir + "/" + input_file_name);
  TH2D *hist_2D = (TH2D*)gDirectory -> Get("tdrift_vs_dqdx");
  //TH2D *hist_2D = (TH2D*)gDirectory -> Get("tdrift_vs_corr_dqdx");

  hist_2D -> RebinX(rebin_x);
  hist_2D -> RebinY(rebin_y);

  int N_binsX = hist_2D -> GetNbinsX();
  int N_binsY = hist_2D -> GetNbinsY();

  for(int i = 1; i < N_binsX + 1; i++){
    TString i_str = Form("%d", i);
    double this_tdrift = hist_2D -> GetXaxis() -> GetBinCenter(i);
    if(this_tdrift > tdrift_high || this_tdrift < tdrift_low) continue;
    double this_tdrift_err = 0.5 * hist_2D -> GetXaxis() -> GetBinWidth(i);

    TString tdrift_str = Form("tdrift%.2fto%.2fms", this_tdrift - this_tdrift_err, this_tdrift + this_tdrift_err);
    TString tdrift_latex = Form("t_{drift} : %.2f - %.2f ms", this_tdrift -this_tdrift_err, this_tdrift + this_tdrift_err);
    TString this_hist_name = tdrift_str;

    TH1D * this_hist_1D = new TH1D(this_hist_name, this_hist_name, N_binsY, 0., 3000.);

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
    fitting_range[0] = 550.;
    fitting_range[1] = 2000.;
    Double_t sv[6], pllo[4], plhi[4], fp[4], fpe[4];
    sv[0] = 30.;
    sv[1] = 1000.; //max_x;
    sv[2] = this_hist_1D -> Integral() * 0.05;
    sv[3] = 60.;
    sv[4] = 90.;
    sv[5] = -0.03;
    for(int j=0; j<4; ++j){
      pllo[j] = 0.01*sv[j];
      plhi[j] = 100*sv[j];
    }

    Double_t chisqr;
    Int_t    ndf;
    Int_t    status;

    //TF1 *this_Langau_fit = langaufit(this_hist_1D_clone, fitting_range, sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status, "Langau_this_Langau" + this_hist_name);
    TF1 *this_Langau_fit = langaubkgfit(this_hist_1D_clone, fitting_range, sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status, "Langau_this_Langau" + this_hist_name);

    TF1 *this_Langaubkg =  new TF1("this_Langaubkg", langaubkgfun, fitting_range[0], fitting_range[1], 6);
    this_Langaubkg -> SetParameters(this_Langau_fit -> GetParameters());
    this_Langaubkg -> SetNpx(1000);
    this_Langaubkg -> SetLineColor(kMagenta);
    this_Langaubkg -> SetLineWidth(2);
    this_Langaubkg -> Draw("lsame");

    TF1 *this_bkg = new TF1("this_bkg", bkgfun, fitting_range[0], fitting_range[1], 2);
    this_bkg -> SetParameters(this_Langau_fit -> GetParameter(4), this_Langau_fit -> GetParameter(5));
    this_bkg -> SetNpx(1000);
    this_bkg -> SetLineColor(kRed - 7);
    this_bkg -> SetLineWidth(2);
    this_bkg -> Draw("lsame");

    TF1 *this_Langau = new TF1("Langau_this_Langau", langaufun, fitting_range[0], fitting_range[1], 4);
    this_Langau -> SetParameters(this_Langau_fit -> GetParameter(0), this_Langau_fit -> GetParameter(1), this_Langau_fit -> GetParameter(2), this_Langau_fit -> GetParameter(3));
    this_Langau -> SetNpx(1000);
    this_Langau -> SetLineColor(kBlue);
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
    //this_Langau -> Draw("lsame");

    TLegend *l = new TLegend(0.50, 0.40, 0.92, 0.85);
    l -> AddEntry(this_hist_1D, "dQ/dx", "pl");
    l -> AddEntry(this_Langau, Form("#sigma_{Landau} : %.2f #pm %.2f", this_Landau_sigma, this_Landau_sigma_err), "l");
    l -> AddEntry(this_hist_1D, Form("MPV : %.2f #pm %.2f", this_MPV, this_MPV_err), "");
    //l -> AddEntry(this_hist_1D, Form("Par2 : %.2f #pm %.2f", this_par2, this_par2_err), "");
    l -> AddEntry(this_hist_1D, Form("#sigma_{Gaus} : %.2f #pm %.2f", this_Gaus_sigma, this_Gaus_sigma_err), "");
    l -> AddEntry(this_hist_1D, Form("#chi^{2} / ndf : %.2f", chisqr / ndf), "");
    l -> AddEntry(this_bkg, "Background", "l");
    l -> AddEntry(this_Langaubkg, "Background + LanGau", "l");

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
    latex_particle.DrawLatex(0.95, 0.96, "Anode-Cathode Passing Tracks");
    latex_method.DrawLatex(0.18, 0.87, tdrift_latex);

    TString output_plot_dir = getenv("PLOT_PATH");
    output_plot_dir = output_plot_dir + "/lifetime/1D/";
    //output_plot_dir = output_plot_dir + "/lifetime/1D/Corr/";
    c -> SaveAs(output_plot_dir + this_hist_name + ".pdf");

    c -> Close();
    
    fitting_results[this_id + "_tdrift"].push_back(this_tdrift);
    fitting_results[this_id + "_tdrift_err"].push_back(this_tdrift_err);
    fitting_results[this_id + "_MPV"].push_back(this_MPV);
    fitting_results[this_id + "_MPV_err"].push_back(this_MPV_err);
    fitting_results[this_id + "_sigma_gaus"].push_back(this_Gaus_sigma);
    fitting_results[this_id + "_sigma_gaus_err"].push_back(this_Gaus_sigma_err);
    fitting_results[this_id + "_sigma_Landau"].push_back(this_Landau_sigma);
    fitting_results[this_id + "_sigma_Landau_err"].push_back(this_Landau_sigma_err);
  }
}

void Fit_lifetime(TString id, double x_range_down, double x_range_up, double y_range_down, double y_range_up){

  double fit_x_low = 0.3;
  double fit_x_high = 1.2;

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D * template_h = new TH1D("", "", 1., x_range_down, x_range_up);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("t_{drift} [ms]");
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("dQ/dx MPV [ADC/cm]");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(y_range_down, y_range_up);
  template_h -> Draw();

  TGraphErrors *gr_1 = new TGraphErrors(fitting_results[id + "_tdrift"].size(), &fitting_results[id + "_tdrift"][0], &fitting_results[id + "_MPV"][0],
				       &fitting_results[id + "_tdrift_err"][0], &fitting_results[id + "_MPV_err"][0]);
  gr_1 -> SetMarkerColor(kBlack);
  gr_1 -> SetMarkerStyle(32);
  gr_1 -> SetMarkerSize(0.7);
  gr_1 -> SetLineColor(kBlack);
  gr_1 -> SetLineWidth(2);
  gr_1 -> Draw("epsame");

  TF1 *gr_fit = new TF1("gr_fit", "[0] * exp(-1. * x/([1]))", fit_x_low, fit_x_high);
  gr_fit -> SetParameters(1050., 10.);
  gr_1 -> Fit(gr_fit, "RN", "", fit_x_low, fit_x_high);
  gr_fit -> SetLineColor(kRed);
  gr_fit -> SetLineWidth(3);
  gr_fit -> SetLineStyle(7);
  gr_fit -> Draw("lsame");

  TLegend *l = new TLegend(0.55, 0.65, 0.92, 0.92);
  l -> AddEntry(gr_1, "MPV from Landau*Gaussian fit", "lp");
  l -> AddEntry(gr_fit, Form("dQ/dx at APA : %.2f #pm %.2f [ADC/cm]", gr_fit -> GetParameter(0), gr_fit -> GetParError(0)), "l");
  l -> AddEntry(gr_fit, Form("Lifetime : %.2f #pm %.2f [ms]", gr_fit -> GetParameter(1), gr_fit -> GetParError(1)), "");
  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_particle;
  latex_ProtoDUNE.SetNDC();
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND Simulation} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_particle.SetNDC();
  latex_particle.SetTextSize(0.03);
  latex_particle.SetTextAlign(31);
  latex_particle.DrawLatex(0.95, 0.96, "Anode-Cathode Passing Tracks");

  TString output_plot_dir = getenv("PLOT_PATH");
  output_plot_dir = output_plot_dir + "/lifetime/";
  c -> SaveAs(output_plot_dir + "life_time.pdf");
  //c -> SaveAs(output_plot_dir + "life_time_corr.pdf");

  c -> Close();
}

void run_lifetime_fit(){

  setTDRStyle();
  Fit_1D_plots("output_lifetime.root", 2, 10, 0.3, 1.2);
  Fit_lifetime("id", 0.2, 1.3, 1010., 1100.);
}
