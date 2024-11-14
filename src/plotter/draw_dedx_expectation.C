#include "canvas_margin.h"
#include "BetheBloch.h"
#include "mylib.h"
#include <iostream>

TSpline3 * muon_sp_range_to_KE = Get_sp_range_KE(mass_muon);
TSpline3 * proton_sp_range_to_KE = Get_sp_range_KE(mass_proton);

void draw_dedx_expectation(){

  double mean_pitch = 0.41;

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
  template_h -> GetYaxis() -> SetTitle("#frac{dE}{dx} [MeV/cm]");
  template_h -> GetYaxis() -> SetTitleOffset(0.85);
  template_h -> GetYaxis() -> SetTitleSize(0.045);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetZaxis() -> SetRangeUser(0., 1. * 1.1);
  template_h -> GetZaxis() -> SetTitle("Entries");
  template_h -> GetZaxis() -> SetTitleSize(0.045);
  template_h -> GetZaxis() -> SetLabelSize(0.035);
  template_h -> Draw();
  
  int N_BB_points = 400.;
  double rr_min = 0.5;
  double rr_max = 200.;
  double rr_step = (rr_max - rr_min) / (N_BB_points + 0.);
  vector<double> BB_rr;
  vector<double> BB_dedx_MPV_muon;
  vector<double> BB_dedx_MPV_proton;
  for(int i = 0; i < N_BB_points; i++){
    double this_rr = rr_min + rr_step * (i + 0.);
    BB_rr.push_back(this_rr);

    double this_KE = muon_sp_range_to_KE -> Eval(this_rr);
    double gamma = (this_KE/mass_muon)+1.0;
    double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
    double this_xi = Landau_xi(this_KE, mean_pitch, mass_muon);
    double this_Wmax = Get_Wmax(this_KE, mass_muon);
    double this_kappa = this_xi / this_Wmax;
    double this_dEdx_BB = meandEdx(this_KE, mass_muon);
    double par[5] = {this_kappa, beta * beta, this_xi, this_dEdx_BB, mean_pitch};
    TF1 * this_dEdx_PDF = dEdx_PDF(par);
    double this_dEdx_MPV = this_dEdx_PDF -> GetMaximumX();
    BB_dedx_MPV_muon.push_back(this_dEdx_MPV);
    
    this_KE = proton_sp_range_to_KE -> Eval(this_rr);
    gamma = (this_KE/mass_proton)+1.0;
    beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
    this_xi = Landau_xi(this_KE, mean_pitch, mass_proton);
    this_Wmax = Get_Wmax(this_KE, mass_proton);
    this_kappa = this_xi / this_Wmax;
    this_dEdx_BB = meandEdx(this_KE, mass_proton);
    double par_2[5] = {this_kappa, beta * beta, this_xi, this_dEdx_BB, mean_pitch};
    this_dEdx_PDF = dEdx_PDF(par_2);
    this_dEdx_MPV = this_dEdx_PDF -> GetMaximumX();
    BB_dedx_MPV_proton.push_back(this_dEdx_MPV);
  }

  TGraph *gr_BB_dedx_muon = new TGraph(N_BB_points, &BB_rr[0], &BB_dedx_MPV_muon[0]);
  gr_BB_dedx_muon -> SetLineColor(kRed);
  gr_BB_dedx_muon -> SetLineStyle(7);
  gr_BB_dedx_muon -> SetLineWidth(2);
  gr_BB_dedx_muon -> Draw("lsame");

  TGraph *gr_BB_dedx_proton = new TGraph(N_BB_points, &BB_rr[0], &BB_dedx_MPV_proton[0]);
  gr_BB_dedx_proton -> SetLineColor(kBlue);
  gr_BB_dedx_proton -> SetLineStyle(7);
  gr_BB_dedx_proton -> SetLineWidth(2);
  gr_BB_dedx_proton -> Draw("lsame");
  
  TLegend *l = new TLegend(0.2 , 0.80, 0.85, 0.95);
  l -> AddEntry(gr_BB_dedx_proton, "Proton Expectation", "l");
  l -> AddEntry(gr_BB_dedx_muon, "Muon Expectation", "l");
  l -> SetBorderSize(0);
  l -> SetFillStyle(0);
  l -> Draw("same");

  
  gPad->RedrawAxis();
  
  TString output_plot_dir = getenv("PLOT_PATH");
  TString outfile_str = output_plot_dir + "/dedx_profile/rr_vs_dedx_expectation.pdf";
  c -> SaveAs(outfile_str);

  c -> Close();  


  cout << "rr\tmuon dE/dx\tproton dE/dx" << endl;
  for(unsigned int i = 0; i < N_BB_points; i++){
    cout << BB_rr.at(i) << "\t" << BB_dedx_MPV_muon.at(i) << "\t" << BB_dedx_MPV_proton.at(i) << endl;
  }

}
