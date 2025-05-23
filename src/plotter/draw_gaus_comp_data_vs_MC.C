#include "canvas_margin.h"
#include "mylib.h"
#include "LanGausFit.h"
#include "TRandom3.h"
#include "BetheBloch.h"
#include <iostream>

BetheBloch *muon_BB = new BetheBloch(13);
//TSpline3 * muon_sp_range_to_KE = Get_sp_range_KE(mass_muon);

bool isdata = false;
TString run_str = "";

void draw(TString y_var, double x_min, double x_max, double y_min, double y_max, TString y_title){

  TString this_id = "id";
  TString gr_name = "rr_vs_" + y_var;

  TString input_file_dir = getenv("OUTPUTROOT_PATH");
  TFile *f_data = new TFile(input_file_dir + "/gaus_comp_fit_" + run_str + ".root");
  TGraphErrors *gr_plane0_data = (TGraphErrors*)gDirectory -> Get("plane0_" + gr_name);
  TGraphErrors *gr_plane1_data = (TGraphErrors*)gDirectory -> Get("plane1_" + gr_name);
  TGraphErrors *gr_plane2_data = (TGraphErrors*)gDirectory -> Get("plane2_" + gr_name);

  TFile *f_mc = new TFile(input_file_dir + "/gaus_comp_fit_MC.root");
  TGraphErrors *gr_plane0_mc = (TGraphErrors*)gDirectory -> Get("plane0_" + gr_name);
  TGraphErrors *gr_plane1_mc = (TGraphErrors*)gDirectory -> Get("plane1_" + gr_name);
  TGraphErrors *gr_plane2_mc = (TGraphErrors*)gDirectory -> Get("plane2_" + gr_name);
  
  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D * template_h = new TH1D("", "", 1., x_min, x_max);
  template_h -> SetStats(0);
  template_h -> GetYaxis() -> SetRangeUser(y_min, y_max);
  template_h -> GetXaxis() -> SetTitle("Residual range [cm]");
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(y_title);
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> Draw();

  gr_plane0_data -> SetLineColor(kBlue);
  gr_plane0_data -> SetMarkerColor(kBlue);
  gr_plane0_data -> Draw("epsame");
  gr_plane0_data -> SetMarkerStyle(22);
  
  gr_plane1_data -> SetLineColor(kRed);
  gr_plane1_data -> SetMarkerColor(kRed);
  gr_plane1_data -> Draw("epsame");
  gr_plane1_data -> SetMarkerStyle(22);

  gr_plane2_data -> SetLineColor(kGreen);
  gr_plane2_data -> SetMarkerColor(kGreen);
  gr_plane2_data -> Draw("epsame");
  gr_plane2_data -> SetMarkerStyle(22);

  gr_plane0_mc -> SetLineColor(kBlue);
  gr_plane0_mc -> SetLineStyle(2);
  gr_plane0_mc -> SetMarkerColor(kBlue);
  gr_plane0_mc -> Draw("epsame");
  gr_plane0_mc -> SetMarkerStyle(32);

  gr_plane1_mc -> SetLineColor(kRed);
  gr_plane1_mc -> SetLineStyle(2);
  gr_plane1_mc -> SetMarkerColor(kRed);
  gr_plane1_mc -> Draw("epsame");
  gr_plane1_mc -> SetMarkerStyle(32);

  gr_plane2_mc -> SetLineColor(kGreen);
  gr_plane2_mc -> SetLineStyle(2);
  gr_plane2_mc -> SetMarkerColor(kGreen);
  gr_plane2_mc -> Draw("epsame");
  gr_plane2_mc -> SetMarkerStyle(32);

  TH1D *data_temp = new TH1D("", "", 1., x_min, x_max);
  data_temp -> SetMarkerStyle(22);
  TH1D *mc_temp = new TH1D("", "", 1., x_min, x_max);
  mc_temp -> SetMarkerStyle(32);
  mc_temp -> SetLineStyle(2);

  TLegend *l = new TLegend(0.4, 0.6, 0.85, 0.90);
  l -> AddEntry(data_temp, "Data Run " + run_str, "lp");
  l -> AddEntry(mc_temp, "MC", "lp");
  l -> AddEntry(gr_plane0_data, "Plane 0", "lp");
  l -> AddEntry(gr_plane1_data, "Plane 1", "lp");
  l -> AddEntry(gr_plane2_data, "Plane 2", "lp");

  if(y_var == "MPV"){
    double mean_pitch = 0.46;
    
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
    gr_BB_dedx -> SetLineColor(kBlack);
    gr_BB_dedx -> SetLineStyle(7);
    gr_BB_dedx -> SetLineWidth(2);
    gr_BB_dedx -> Draw("lsame");
    l -> AddEntry(gr_BB_dedx, Form("Muon expectation, pitch = %.2f cm", mean_pitch), "l"); 
  }

  if(y_var == "MPV"){
    double mean_pitch = 0.33;

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
    gr_BB_dedx -> SetLineColor(kBlack);
    gr_BB_dedx -> SetLineStyle(2);
    gr_BB_dedx -> SetLineWidth(3);
    gr_BB_dedx -> Draw("lsame");
    l -> AddEntry(gr_BB_dedx, Form("Muon expectation, pitch = %.2f cm", mean_pitch), "l");
  }

  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_particle;
  latex_ProtoDUNE.SetNDC();
  latex_particle.SetNDC();
  latex_particle.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_particle.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND: Data} Run " + run_str);
  latex_particle.DrawLatex(0.95, 0.96, "Cathode Passing Stopping Tracks");

  TString output_plot_dir = getenv("PLOT_PATH");
  TString outfile_str = output_plot_dir + "/comparison/plane_comparison_" + y_var + "_" + run_str + ".pdf";
  c -> SaveAs(outfile_str);
  c -> Close();
  
}

void draw_gaus_comp_data_vs_MC(){
   
  setTDRStyle();
  draw("MPV", 5., 200., 1., 4., "dE/dx MPV");
  draw("sigmaG", 5., 200., 0., 1., "#sigma_{G}");
  draw("sigmaL", 5., 200., 0., 0.2, "#sigma_{L}");

}
