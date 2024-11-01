#include "canvas_margin.h"
#include "mylib.h"
#include "LanGausFit.h"
#include "TRandom3.h"
#include "BetheBloch.h"
#include <iostream>

TSpline3 * muon_sp_range_to_KE = Get_sp_range_KE(mass_muon);

bool isdata = false;
TString run_str = "";

map<TString, double> map_fitpar;

void draw(TString plane, TString y_var, double x_min, double x_max, double y_min, double y_max, TString y_title){

  TString this_id = "id";
  TString gr_name = "rr_vs_" + y_var;

  TString suffix[] = {"", "_NE", "_NW", "_SE", "_SW", "_cafv", "_meddqdx"};
  
  TString input_file_dir = getenv("OUTPUTROOT_PATH");
  TFile *f_data = new TFile(input_file_dir + "/gaus_comp_fit_" + run_str + ".root");
  vector<TGraphErrors*> gr_vec;
  for(int i = 0; i < 7; i++){
    TString this_suffix = suffix[i];
    TString this_gr_name = plane + this_suffix + "_" + gr_name;
    TGraphErrors *this_gr = (TGraphErrors*)gDirectory -> Get(this_gr_name) -> Clone();
    gr_vec.push_back(this_gr);
  }

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

  TLegend *l = new TLegend(0.4, 0.6, 0.85, 0.90);
  int colors[] = {1, 632, 800, 400, 416, 600, 880};
  TString legend_str[] = {"Central", "NE", "NW", "SE", "SW", "CA-FV", "Med. dQ/dx"};
  for(int i = 0; i < 7; i++){
    TGraphErrors *this_gr = gr_vec.at(i);
    int this_color = colors[i];
    if(i == 0){
      this_gr -> SetLineColor(kBlack);
      this_gr -> SetFillColorAlpha(kGray, 0.7);
      //this_gr -> SetFillStyle(3005);
      this_gr -> Draw("3same");
      l -> AddEntry(this_gr, legend_str[i], "fl");
    }
    else{
      this_gr -> SetLineColor(this_color);
      this_gr -> SetLineWidth(2);
      this_gr -> SetMarkerColor(this_color);
      this_gr -> Draw("epsame");
      l -> AddEntry(this_gr, legend_str[i], "ep");
    }
  }

  gr_vec.at(0) -> Draw("3same");
  TGraphErrors *central_clone = (TGraphErrors*)gr_vec.at(0) -> Clone();
  central_clone -> SetLineColor(kBlack);
  central_clone -> SetLineWidth(2);
  central_clone -> Draw("lsame");
  /*
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
  */
  
  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_particle;
  latex_ProtoDUNE.SetNDC();
  latex_particle.SetNDC();
  latex_particle.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_particle.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND: Data} Run " + run_str + ", " + plane + ", #font[42]{#it{#scale[1.0]{Preliminary}}}");
  latex_particle.DrawLatex(0.95, 0.96, "Cathode Passing Stopping Tracks");

  TString output_plot_dir = getenv("PLOT_PATH");
  TString outfile_str = output_plot_dir + "/comparison/syst_" + plane + "_comparison_" + y_var + "_" + run_str + ".pdf";
  c -> SaveAs(outfile_str);
  c -> Close();
  
}

void draw_gaus_comp_syst(int run_num = 0){
  
  if(run_num != 0){
    isdata = true;
    run_str = TString::Format("%d", run_num);
  }
  
  setTDRStyle();
  TString planes[] = {"plane0", "plane1", "plane2"};
  for(int i = 0; i < 3; i++){
    TString this_plane = planes[i];
    draw(this_plane, "MPV", 7., 200., 1., 4., "dE/dx MPV");
    draw(this_plane, "sigmaG", 7., 200., 0., 1., "#sigma_{G}");
    draw(this_plane, "sigmaL", 7., 200., 0., 0.3, "#sigma_{L}");
  }
}
