#include "canvas_margin.h"
#include "mylib.h"
#include "LanGausFit.h"
#include "TRandom3.h"
#include "BetheBloch.h"
#include <iostream>

bool isdata = false;
TString run_str = "";

TRandom3 gRan(1800);
map<TString, vector<double>> fitting_results;

void Write_1D_hist(TH1D *in, TString outname, TString particle, TString latex_str, TString plane, TString suffix, TString title_x, TString title_y, double x_min, double x_max, bool do_langau_fit, double this_ResRange){

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  //in -> Rebin(rebin);
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
    fitting_range[1] = 15.0;
    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
    sv[0] = 0.1;
    sv[1] = max_x;
    sv[2] = in -> Integral() * 0.05 * bin_width;
    sv[3] = 0.2;
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
    l -> AddEntry(this_Langau, Form("#sigma_{Landau} : %.2f #pm %.3f", this_Landau_sigma, this_Landau_sigma_err), "l");
    l -> AddEntry(in, Form("MPV : %.2f #pm %.3f", this_MPV, this_MPV_err), "");
    l -> AddEntry(in, Form("#sigma_{Gaus} : %.2f #pm %.3f", this_Gaus_sigma, this_Gaus_sigma_err), "");
    l -> AddEntry(in, Form("Par2 : %.2f #pm %.2f", this_par2, this_par2_err), "");
    l -> AddEntry(in, Form("#chi^{2} / ndf : %.2f", chisqr / ndf), "");
    l -> Draw("same");

    fitting_results[plane + suffix + "res_range"].push_back(this_ResRange);
    fitting_results[plane + suffix + "MPV"].push_back(this_MPV);
    fitting_results[plane + suffix + "MPV_err"].push_back(this_MPV_err);
    fitting_results[plane + suffix + "sigmaL"].push_back(this_Landau_sigma);
    fitting_results[plane + suffix + "sigmaL_err"].push_back(this_Landau_sigma_err);
    fitting_results[plane + suffix + "sigmaG"].push_back(this_Gaus_sigma);
    fitting_results[plane + suffix + "sigmaG_err"].push_back(this_Gaus_sigma_err);
  }
  else{
    return;
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
  if(isdata) latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND Data} Run " + run_str + ", " + plane);
  else latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND Simulation} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_particle.DrawLatex(0.95, 0.96, particle_label_str);
  latex_method.DrawLatex(0.18, 0.87, latex_str);

  TString output_plot_dir = getenv("PLOT_PATH");
  TString outfile_str = output_plot_dir + "/MC/dEdx_1d/" + outname + ".pdf";
  if(isdata) outfile_str = output_plot_dir + "/Run" + run_str + "/dEdx_1d/" + outname + ".pdf";
  c -> SaveAs(outfile_str);
  c -> Close();

}

void Fit_rr_vs_dedx(TString input_file_name, TString plane, TString suffix, double rebin_y = 1, double res_range_cut = 100.){

  //const Int_t num_NbinsX = sizeof(dEdx_MPV_binning) / sizeof(dEdx_MPV_binning[0]) - 1;
  
  TString this_id = "id";
  TString histname = "rr_vs_dedx_" + plane + "_trklen_60cm_passing_cathode_coszx" + suffix;
  
  TString input_file_dir = getenv("OUTPUTROOT_PATH");
  TFile *f = new TFile(input_file_dir + "/" + input_file_name);
  TH2D *hist_2D = (TH2D*)gDirectory -> Get(histname);

  hist_2D -> RebinX(40);
  hist_2D -> RebinY(rebin_y);

  int N_binsX = hist_2D -> GetNbinsX();
  int N_binsY = hist_2D -> GetNbinsY();
  cout << "N_binsX : " << N_binsX << endl;
  cout << "N_binsY : " << N_binsY << endl;

  for(int i = 1; i < N_binsX + 1; i++){
    TString i_str = Form("%d", i);
    double this_ResRange = hist_2D -> GetXaxis() -> GetBinCenter(i);
    if(this_ResRange > res_range_cut) break;
    double this_ResRange_err = 0.5 * hist_2D -> GetXaxis() -> GetBinWidth(i);
    TString ResRange_range_str = Form("ResRange%.1fto%.1fcm", this_ResRange - this_ResRange_err, this_ResRange + this_ResRange_err);
    TString ResRange_range_latex = Form("Residual range : %.1f - %.1f cm", this_ResRange -this_ResRange_err, this_ResRange + this_ResRange_err);
    TString this_hist_name = ResRange_range_str;
    
    TH1D * this_1D = new TH1D(this_hist_name, this_hist_name, N_binsY, 0., 20.);
    for(int j = 1; j < N_binsY + 1; j++){
      double this_content = hist_2D -> GetBinContent(i, j);
      double this_error = hist_2D -> GetBinError(i, j);
      this_1D -> SetBinContent(j, this_content);
      this_1D -> SetBinError(j, this_error);
    }
    
    Write_1D_hist(this_1D, "/" + plane + "_dEdx_" + this_hist_name + suffix, "muon", ResRange_range_latex, plane, suffix, "dE/dx [MeV/cm]", "Events", 0.0, 10.0, true, this_ResRange);
  }
}

void make_gr_and_save(){

  TString output_file_dir = getenv("OUTPUTROOT_PATH");
  TString out_root_name = output_file_dir + "/gaus_comp_fit_" + run_str + ".root";
  TFile *outfile = new TFile(out_root_name, "RECREATE");
  outfile -> cd();
  
  TString planes_str[] = {"plane0", "plane1", "plane2"};
  TString suffix[] = {"", "_NE", "_NW", "_SE", "_SW", "_cafv", "_meddqdx"};
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 7; j++){
      cout << planes_str[i] << endl;
      TString this_plane = planes_str[i];
      TString this_suffix = suffix[j];
      cout << fitting_results[this_plane + "res_range"].size() << endl;
      TGraphErrors *gr_rr_vs_MPV = new TGraphErrors(fitting_results[this_plane + this_suffix + "res_range"].size(),
						    &fitting_results[this_plane + this_suffix + "res_range"][0], &fitting_results[this_plane + this_suffix + "MPV"][0],
						    &fitting_results[this_plane + this_suffix + "res_range_err"][0], &fitting_results[this_plane + this_suffix + "MPV_err"][0]);
      gr_rr_vs_MPV -> SetName(this_plane + this_suffix + "_rr_vs_MPV");
      gr_rr_vs_MPV -> Write();
      
      TGraphErrors *gr_rr_vs_sigmaG = new TGraphErrors(fitting_results[this_plane + this_suffix + "res_range"].size(),
						       &fitting_results[this_plane + this_suffix + "res_range"][0], &fitting_results[this_plane + this_suffix + "sigmaG"][0],
						       &fitting_results[this_plane + this_suffix + "res_range_err"][0], &fitting_results[this_plane + this_suffix + "sigmaG_err"][0]);
      gr_rr_vs_sigmaG -> SetName(this_plane + this_suffix + "_rr_vs_sigmaG");
      gr_rr_vs_sigmaG -> Write();
      
      TGraphErrors *gr_rr_vs_sigmaL = new TGraphErrors(fitting_results[this_plane + this_suffix + "res_range"].size(),
						       &fitting_results[this_plane + this_suffix + "res_range"][0], &fitting_results[this_plane + this_suffix + "sigmaL"][0],
						       &fitting_results[this_plane + this_suffix + "res_range_err"][0], &fitting_results[this_plane + this_suffix + "sigmaL_err"][0]);
      gr_rr_vs_sigmaL -> SetName(this_plane + this_suffix + "_rr_vs_sigmaL");
      gr_rr_vs_sigmaL -> Write();
    }
  }
  
  outfile -> Close();
}


void run_gaus_comp_fit(int run_num = 0){

  if(run_num != 0){
    isdata = true;
    run_str = TString::Format("%d", run_num);
  }
  
  setTDRStyle();

  TString filename = "output_dedx_" + run_str + "_syst.root";
  if(!isdata){
    filename = "output_dedx_2023B_GENIE_CV.root";
    run_str = "MC";
  }

  TString suffix[] = {"", "_NE", "_NW", "_SE", "_SW", "_cafv", "_meddqdx"};
  for(int i = 0; i < 7; i++){
    Fit_rr_vs_dedx(filename, "plane0", suffix[i], 10., 200.);
    Fit_rr_vs_dedx(filename, "plane1", suffix[i], 10., 200.);
    Fit_rr_vs_dedx(filename, "plane2", suffix[i], 10., 200.);
  }
  
  make_gr_and_save();
}
