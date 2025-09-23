#include "canvas_margin.h"
#include "mylib.h"
#include "LanGausFit.h"
#include "TRandom3.h"

TRandom3 gRan(1800);
map<TString, vector<double>> fitting_results;
TString sample_str = "run_18255_and_18259";
void Fit_1D_plots(TString input_file_name, int rebin_x, int rebin_y, double tdrift_low, double tdrift_high, TString side, TString ngroupwires, TString plane, TString suffix, bool fit_bkg = false){

  TString input_file_dir = getenv("OUTPUTROOT_PATH");
  TFile *f = new TFile(input_file_dir + "/lifetime/whicht00/" + input_file_name);
  TString histname = "h_dQdx_tDrift" + side + "_" + ngroupwires + "wires_" + plane;

  cout << "[Fit_1D_plots] histname + suffix: " << histname + suffix << endl;
  TH2D *hist_2D = (TH2D*)gDirectory -> Get(histname + suffix);

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

    TLegend *l = new TLegend(0.55, 0.40, 0.92, 0.85);

    double max_x = this_hist_1D -> GetBinCenter(this_hist_1D -> GetMaximumBin());
    double width_x = this_hist_1D -> GetBinWidth(1);
    Double_t fitting_range[2];
    fitting_range[0] = 300.;
    fitting_range[1] = 1800.;
    Double_t sv[6], pllo[4], plhi[4], fp[4], fpe[4];
    sv[0] = 30.;
    //sv[1] = 1000.; //max_x;
    sv[1] = 1100.; //max_x; after angle cut
    sv[2] = this_hist_1D -> Integral() * 0.05 * width_x;
    sv[3] = 60.;
    sv[4] = 30.;
    sv[5] = -0.01;
    for(int j=0; j<4; ++j){
      pllo[j] = 0.01*sv[j];
      plhi[j] = 100*sv[j];
    }

    Double_t chisqr;
    Int_t    ndf;
    Int_t    status;

    double this_Landau_sigma = -1.;
    double this_Landau_sigma_err = -1.;
    double this_MPV = -1.;
    double this_MPV_err = -1.;
    double this_par2 = -1.;
    double this_par2_err = -1.;
    double this_Gaus_sigma = -1.;
    double this_Gaus_sigma_err = -1.;

    if(fit_bkg){

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

      this_Landau_sigma = this_Langau_fit ->GetParameter(0);
      this_Landau_sigma_err =this_Langau_fit -> GetParError(0);
      this_MPV = this_Langau_fit -> GetParameter(1);
      this_MPV_err = this_Langau_fit -> GetParError(1);
      this_par2 = this_Langau_fit -> GetParameter(2);
      this_par2_err =this_Langau_fit -> GetParError(2);
      this_Gaus_sigma= this_Langau_fit -> GetParameter(3);
      this_Gaus_sigma_err= this_Langau_fit -> GetParError(3);

      this_hist_1D -> Draw("epsame");
      //this_Langau -> Draw("lsame");

      l -> AddEntry(this_hist_1D, "dQ/dx", "pl");
      l -> AddEntry(this_Langau, Form("#sigma_{Landau} : %.2f #pm %.2f", this_Landau_sigma, this_Landau_sigma_err), "l");
      l -> AddEntry(this_hist_1D, Form("MPV : %.2f #pm %.2f", this_MPV, this_MPV_err), "");
      l -> AddEntry(this_hist_1D, Form("#sigma_{Gaus} : %.2f #pm %.2f", this_Gaus_sigma, this_Gaus_sigma_err), "");
      l -> AddEntry(this_hist_1D, Form("#chi^{2} / ndf : %.2f", chisqr / ndf), "");
      l -> AddEntry(this_bkg, "Background", "l");
      l -> AddEntry(this_bkg, Form("P_{0}: %.2e #pm %.2e", this_Langau_fit -> GetParameter(4), this_Langau_fit -> GetParError(4)), "");
      l -> AddEntry(this_bkg, Form("P_{1}: %.2e #pm %.2e", this_Langau_fit -> GetParameter(5), this_Langau_fit -> GetParError(5)), "");
      l -> AddEntry(this_Langaubkg, "Background + LanGau", "l");

    }
    else{
      TF1 *this_Langau_fit = langaufit(this_hist_1D_clone, fitting_range, sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status, "Langau_this_Langau" + this_hist_name);
    
      TF1 *this_Langau =  new TF1("this_Langau", langaufun, fitting_range[0], fitting_range[1], 4);
      this_Langau -> SetParameters(this_Langau_fit -> GetParameters());
      this_Langau -> SetNpx(1000);
      this_Langau -> SetLineColor(kMagenta);
      this_Langau -> SetLineWidth(2);
      this_Langau -> Draw("lsame");

      this_Landau_sigma = this_Langau_fit ->GetParameter(0);
      this_Landau_sigma_err =this_Langau_fit -> GetParError(0);
      this_MPV = this_Langau_fit -> GetParameter(1);
      this_MPV_err = this_Langau_fit -> GetParError(1);
      this_par2 = this_Langau_fit -> GetParameter(2);
      this_par2_err =this_Langau_fit -> GetParError(2);
      this_Gaus_sigma= this_Langau_fit -> GetParameter(3);
      this_Gaus_sigma_err= this_Langau_fit -> GetParError(3);

      this_hist_1D -> Draw("epsame");
      //this_Langau -> Draw("lsame");

      l -> AddEntry(this_hist_1D, "dQ/dx", "pl");
      l -> AddEntry(this_Langau, Form("#sigma_{Landau} : %.2f #pm %.2f", this_Landau_sigma, this_Landau_sigma_err), "l");
      l -> AddEntry(this_hist_1D, Form("MPV : %.2f #pm %.2f", this_MPV, this_MPV_err), "");
      //l -> AddEntry(this_hist_1D, Form("Par2 : %.2f #pm %.2f", this_par2, this_par2_err), "");
      l -> AddEntry(this_hist_1D, Form("#sigma_{Gaus} : %.2f #pm %.2f", this_Gaus_sigma, this_Gaus_sigma_err), "");
      l -> AddEntry(this_hist_1D, Form("#chi^{2} / ndf : %.2f", chisqr / ndf), "");
    }
    
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
    TString outfile_str = output_plot_dir + "/lifetime/whicht00/" + sample_str + "/1D/wire" + ngroupwires + "/" + plane + "_" + side + "/" + tdrift_str + suffix + ".pdf";
    TString outfile_dir = gSystem->DirName(outfile_str);
    if (gSystem->AccessPathName(outfile_dir)) {
      std::cout << "Directory does not exist, creating: " << outfile_dir << std::endl;
      gSystem->mkdir(outfile_dir, kTRUE);
    }
    c -> SaveAs(outfile_str);

    c -> Close();
    TString this_id = plane + side + ngroupwires;
    fitting_results[this_id + "_tdrift" + suffix].push_back(this_tdrift);
    fitting_results[this_id + "_tdrift_err" + suffix].push_back(this_tdrift_err);
    fitting_results[this_id + "_MPV" + suffix].push_back(this_MPV);
    fitting_results[this_id + "_MPV_err" + suffix].push_back(this_MPV_err);
    fitting_results[this_id + "_sigma_gaus" + suffix].push_back(this_Gaus_sigma);
    fitting_results[this_id + "_sigma_gaus_err" + suffix].push_back(this_Gaus_sigma_err);
    fitting_results[this_id + "_sigma_Landau" + suffix].push_back(this_Landau_sigma);
    fitting_results[this_id + "_sigma_Landau_err" + suffix].push_back(this_Landau_sigma_err);
  }
}

void Fit_lifetime(TString side, TString ngroupwires, TString plane, TString suffix, double x_range_down, double x_range_up, double y_range_down, double y_range_up){

  TString id = plane + side + ngroupwires;
  double fit_x_low = 0.1;
  double fit_x_high = 1.24;

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

  TGraphErrors *gr_1 = new TGraphErrors(fitting_results[id + "_tdrift" + suffix].size(), &fitting_results[id + "_tdrift" + suffix][0], &fitting_results[id + "_MPV" + suffix][0],
				       &fitting_results[id + "_tdrift_err" + suffix][0], &fitting_results[id + "_MPV_err" + suffix][0]);
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

  TLegend *l = new TLegend(0.45, 0.65, 0.92, 0.92);
  l -> AddEntry(gr_1, "MPV from Landau*Gaussian fit", "lp");
  l -> AddEntry(gr_fit, Form("dQ/dx at APA : %.2f #pm %.2f [ADC/cm]", gr_fit -> GetParameter(0), gr_fit -> GetParError(0)), "l");
  l -> AddEntry(gr_fit, Form("Lifetime : %.2f #pm %.3f [ms]", gr_fit -> GetParameter(1), gr_fit -> GetParError(1)), "");
  l -> SetFillColor(0);
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
  output_plot_dir = output_plot_dir + "/lifetime/whicht00/";
  TString outfile_str = output_plot_dir + sample_str + "/lifetime_wire" + ngroupwires + "_" + plane + "_" + side + suffix + ".pdf";
  TString outfile_dir = gSystem->DirName(outfile_str);
  if (gSystem->AccessPathName(outfile_dir)) {
    std::cout << "Directory does not exist, creating: " << outfile_dir << std::endl;
    gSystem->mkdir(outfile_dir, kTRUE);
  }
  c -> SaveAs(outfile_str);
  gr_1 -> SetName(ngroupwires + "_" + plane + "_" + side + suffix);
  gr_1 -> Write();
  c -> Close();
}

void run_lifetime_fit_sce_efield_corr(TString filename, TString outfilename, TString sample_str_){

  sample_str = sample_str_;

  setTDRStyle();

  //TString filename = "output_lifetime_2025A_Sprint25Dev_data.root";
  //TString outfilename = "fit_lifetime_2025A_SpringDev_data_bnbcosmics.root";

  //TString filename = "output_lifetime_2025SpringFinalValid_data_bnb_light.root";
  //TString outfilename = "fit_lifetime_2025B_final_valid_data_bnb_light.root";
  
  TString sides[] = {"E", "W"};
  TString planes[] = {"plane0", "plane1", "plane2"};
  TString dedx_assume_str[] = {"1p7", "1p8", "1p9", "2p0", "2p1", "2p2", "2p3"};
  int N_dedx_assume = 7;
  
  for(int i = 0; i < 3; i++){
    TString this_plane = planes[i];

    if(i != 2) continue; // FIXME: running only collection plane
    for(int j = 0; j < 2; j++){
      TString this_side = sides[j];
      Fit_1D_plots(filename, 2, 2, 0.1, 1.24, this_side, "10", this_plane, "", false);
      Fit_1D_plots(filename, 2, 2, 0.1, 1.24, this_side, "10", this_plane, "_sce_corr", false);

      for(int k = 0; k < N_dedx_assume; k++){
	TString this_dedx_assume_str = dedx_assume_str[k];
	//Fit_1D_plots(filename, 2, 2, 0.1, 1.24, this_side, "10", this_plane, "_sce_corr_e_distor_mb_dedx" + this_dedx_assume_str, false);
	Fit_1D_plots(filename, 2, 2, 0.1, 1.24, this_side, "10", this_plane, "_sce_corr_e_distor_emb_dedx" + this_dedx_assume_str, false);
      }
    }
  }

  //Fit_1D_plots("output_lifetime_test.root", 2, 20, 0.1, 1.24, "space", "angle_passing_cathode");
  //Fit_1D_plots("output_lifetime_test.root", 2, 20, 0.1, 1.24, "time", "angle_passing_cathode");
  //Fit_1D_plots("output_lifetime_test.root", 2, 20, 0.1, 1.24, "space", "trk_len", true);
  //Fit_1D_plots("output_lifetime_test.root", 2, 20, 0.1, 1.24, "time", "trk_len", true);
  //Fit_1D_plots("output_lifetime_test.root", 2, 20, 0.1, 1.24, "space", "corr_trk_len", true);
  //Fit_1D_plots("output_lifetime_test.root", 2, 20, 0.1, 1.24, "time", "corr_trk_len", true);

  //Fit_lifetime("space", "angle_passing_cathode", 0., 1.3, 950., 1050.);
  //Fit_lifetime("time", "angle_passing_cathode", 0., 1.3, 950., 1050.);
  //Fit_lifetime("space", "trk_len", 0., 1.3, 950., 1080.);
  //Fit_lifetime("time", "trk_len", 0., 1.3, 950., 1080.);
  //Fit_lifetime("space", "corr_trk_len", 0., 1.3, 1010., 1100.);
  //Fit_lifetime("time", "corr_trk_len", 0., 1.3, 1010., 1100.);
  TString output_file_dir = getenv("OUTPUTROOT_PATH");
  TString out_root_name = output_file_dir + "/lifetime/whicht00/fit_lifetime_" + outfilename + ".root";
  TFile *outfile = new TFile(out_root_name, "RECREATE");
  outfile -> cd();
  for(int i = 0; i < 3; i++){
    TString this_plane = planes[i];
    if(i != 2) continue; // FIXME: running only collection plane 

    for(int j = 0; j < 2; j++){
      TString this_side = sides[j];
      TString this_id = this_plane + this_side + "10";
      Fit_lifetime(this_side, "10", this_plane, "", 0., 1.3, 1000., 1200.);
      Fit_lifetime(this_side, "10", this_plane, "_sce_corr", 0., 1.3, 1000., 1200.);

      for(int k = 0; k < N_dedx_assume; k++){
	TString this_dedx_assume_str = dedx_assume_str[k];
	//Fit_lifetime(this_side, "10", this_plane, "_sce_corr_e_distor_mb_dedx" + this_dedx_assume_str, 0., 1.3, 1000., 1200.);
	Fit_lifetime(this_side, "10", this_plane, "_sce_corr_e_distor_emb_dedx" + this_dedx_assume_str, 0., 1.3, 1000., 1200.);
      }
    }
  }

}
