#include "canvas_margin.h"
#include "mylib.h"
#include <iostream>

bool isdata = false;
TString run_str = "";

// == EBM par
double alpha_ebm = 0.904;
//double beta_90 = 0.204;
double R_ebm = 1.25;

double B_phi(double beta_90 = 0.204, double phi = 90.){
  double phi_rad = phi * TMath::Pi() / 180.;
  double out = beta_90 / (sqrt(pow(sin(phi_rad), 2.) + pow(cos(phi_rad) / R_ebm, 2.)));
  return out;
}

void fit_plane(TString plane, TString particle, double dEdx_low, double dEdx_high, double dQdx_low, double dQdx_high){

  TString suffixes[] = {"_phi40to50", "_phi50to60", "_phi60to70", "_phi70to80", "_phi80to85", "_phi85to90"};
  
  double phis[] = {45., 55., 65., 75., 82.5, 87.5};
  double phis_err[] = {5., 5., 5., 5., 2.5, 2.5};
  vector<TString> y_labels = {"#phi = [40, 50)#circ",  "#phi = [50, 60)#circ",  "#phi = [60, 70)#circ",  "#phi = [70, 80)#circ",  "#phi = [80, 85)#circ",  "#phi = [85, 90)#circ"};
  int phi_color_arr[] = {632, 800, 401, 418, 600, 880};
  TString input_file_dir = getenv("OUTPUTROOT_PATH");
  TFile *f = new TFile(input_file_dir + "/recom_fit_ebm_MC.root");
  vector<TGraphErrors*> grs;
  for(int i = 0; i < 6; i++){
    TString suffix = suffixes[i];
    TString this_gr_str = plane + suffix;
    cout << "Get " << this_gr_str << endl;
    TGraphErrors *this_gr = (TGraphErrors*)gDirectory -> Get(this_gr_str);
    grs.push_back(this_gr);
  }

  // == Draw plots with all grs
  cout << "Draw plots with all grs" << endl;
  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  gStyle->SetLineWidth(2);
  
  TH1D * template_h = new TH1D("", "", 1., dEdx_low, dEdx_high);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("dE/dx [MeV/cm]");
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("dQ/dx [ADC/cm]");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(dQdx_low, dQdx_high);
  template_h -> Draw();

  TLegend *l = new TLegend(0.2, 0.5, 0.40, 0.9);

  vector<double> phi_vec;
  vector<double> phi_err_vec;
  vector<double> beta_vec;
  vector<double> beta_err_vec;
  for(int i = 1; i < 5; i++){
    TGraphErrors* this_gr = grs.at(i);
    this_gr -> Draw("ezpsame");
    l -> AddEntry(this_gr, y_labels.at(i), "lp");

    TF1 * f_mod_box = new TF1("f_mod_box", "(294.49153 * [0] / [1]) * log(1.4388489 * [1] * x + [2])", 1.6, 4.0);
    double this_B_phi = B_phi(0.204, phis[i]);
    f_mod_box -> SetParameters(2.00, this_B_phi, alpha_ebm);
    f_mod_box -> FixParameter(0, 2.00);
    f_mod_box -> FixParameter(2, alpha_ebm);
    f_mod_box -> SetLineColor(phi_color_arr[i]);
    f_mod_box -> SetLineWidth(3);
    f_mod_box -> SetLineStyle(7);
    this_gr -> Fit("f_mod_box", "RNS");
    f_mod_box -> Draw("lsame");
    l -> AddEntry(this_gr, Form("Fit #beta' = %.3f #pm %.3f", f_mod_box -> GetParameter(1), f_mod_box -> GetParError(1)), "");

    phi_vec.push_back(phis[i]);
    phi_err_vec.push_back(phis_err[i]);
    beta_vec.push_back(f_mod_box -> GetParameter(1));
    beta_err_vec.push_back(f_mod_box -> GetParError(1));
  }

  l -> Draw("same");

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
  else latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND Simulation 2024B CV} " + plane + " #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_particle.DrawLatex(0.95, 0.96, particle_label_str);
  //latex_method.DrawLatex(0.18, 0.87, latex_str);
  
  TString output_plot_dir = getenv("PLOT_PATH");
  TString outfile_str = output_plot_dir + "/MC_2024B_CV/recom/dEdx_MPV_vs_corr_dqdx_" + particle + "_" + plane + "_allphi.pdf";
  if(isdata) outfile_str = output_plot_dir + "/Run" + run_str + "/dEdx_MPV_vs_corr_dqdx_" + particle + "_" + plane + "_allphi.pdf";
  c -> SaveAs(outfile_str);

  TH1D * template_h2 = new TH1D("", "", 1., 40., 90.);
  template_h2 -> SetStats(0);
  template_h2 -> GetXaxis() -> SetTitle("#phi [#circ]");
  template_h2 -> GetXaxis() -> SetTitleSize(0.037);
  template_h2 -> GetXaxis() -> SetTitleOffset(1.4);
  template_h2 -> GetXaxis() -> SetLabelSize(0.035);
  template_h2 -> GetYaxis() -> SetTitle("#beta' [(kV/cm)(g/cm^{3})/MeV]");
  template_h2 -> GetYaxis() -> SetTitleSize(0.05);
  template_h2 -> GetYaxis() -> SetLabelSize(0.035);
  template_h2 -> GetYaxis() -> SetRangeUser(0.19, 0.25);
  template_h2 -> Draw();

  TF1 *tf_betap = new TF1("ebm_betap", "0.204 / sqrt(pow(sin(x * TMath::Pi() / 180.), 2.) + pow(cos(x * TMath::Pi() / 180.) / 1.25, 2.))", 40., 90.);
  tf_betap -> SetLineColor(kBlue);
  tf_betap -> SetLineStyle(kBlue);
  tf_betap -> Draw("lsame");

  cout << "tf_betap -> Eval(75.) : " << tf_betap -> Eval(75.) << endl;
  
  TGraphErrors* beta_gr = new TGraphErrors(phi_vec.size(), &phi_vec[0], &beta_vec[0], &phi_err_vec[0], &beta_err_vec[0]);
  beta_gr -> SetLineColor(kRed);
  beta_gr -> SetMarkerColor(kRed);
  beta_gr -> Draw("ezpsame");
  
  TLegend *l2 = new TLegend(0.2, 0.7, 0.7, 0.9);
  l2 -> AddEntry(beta_gr, "Fit #beta' with fixed #alpha = 0.904 and C_{cal.} = 0.02");
  l2 -> AddEntry(tf_betap, "Ellipsoidal Box Model with #beta_{90} = 0.204 and R = 1.25", "l");
  l2 -> Draw("same");
  
  if(isdata) latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND Data} Run " + run_str);
  else latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND Simulation 2024B CV} " + plane + " #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_particle.DrawLatex(0.95, 0.96, particle_label_str);
  outfile_str = output_plot_dir + "/MC_2024B_CV/recom/EBM_beta_" + particle + "_" + plane + "_allphi.pdf";
  if(isdata) outfile_str = output_plot_dir + "/Run" + run_str + "/EBM_beta_" + particle + "_" + plane + "_allphi.pdf";
  c -> SaveAs(outfile_str);
  
  c -> Close();

}

void fit_recom_ebm(int run_num = 0){

  if(run_num != 0){
    isdata = true;
    run_str = TString::Format("%d", run_num);
  }
 if(!isdata)  run_str = "MC";

  setTDRStyle();
  
  fit_plane("plane0", "muon", 1.5, 4.0, 800., 2200.); 
  fit_plane("plane1", "muon", 1.5, 4.0, 800., 2200.);
  fit_plane("plane2", "muon", 1.5, 4.0, 800., 2200.);
 
}
