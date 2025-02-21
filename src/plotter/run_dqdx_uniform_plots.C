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

double trun_mean_min = 800.;
double trun_mean_max = 1400.;

void Write_1D_hist(TH1D *in, TString outname, TString suffix, TString particle, TString latex_str, TString plane, TString title_x, TString title_y, double x_min, double x_max, bool do_langau_fit, double this_zenith, double this_zenith_err){

  TString side_str = "east";
  if(outname.Contains("east")) side_str = "east";
  if(outname.Contains("west")) side_str = "west";
  
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

  double mean = Get_Mean_in_range(in, trun_mean_min, trun_mean_max);
  double mean_err = Get_Mean_Err_in_range(in, trun_mean_min, trun_mean_max);
  
  TLegend *l = new TLegend(0.60, 0.40, 0.80, 0.80);
  l -> AddEntry(in, Form("Tot. Entreis: %.0f, Mean: %.1f #pm %.1f", in -> Integral(), mean, mean_err), "lp");
  l -> AddEntry(in, Form("Mean: %.1f #pm %.1f", in -> Integral(), mean, mean_err), "");
  //if(in -> Integral() > 1000.) do_langau_fit = true;

  if(do_langau_fit){

    TH1D *in_clone = (TH1D*)in -> Clone();
    double max_x = in -> GetBinCenter(in -> GetMaximumBin());
    double bin_width = in -> GetBinWidth(1);
    Double_t fitting_range[2];
    fitting_range[0] = 500.;
    fitting_range[1] = 1500.0;
    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
    sv[0] = 40.;
    sv[1] = 1050.;
    sv[2] = in -> Integral() * 0.05 * bin_width;
    sv[3] = 150.;
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

    //TLegend *l = new TLegend(0.60, 0.40, 0.80, 0.80);
    l -> AddEntry(this_Langau, Form("#sigma_{Landau} : %.2f #pm %.3f", this_Landau_sigma, this_Landau_sigma_err), "l");
    l -> AddEntry(in, Form("MPV : %.2f #pm %.3f", this_MPV, this_MPV_err), "");
    l -> AddEntry(in, Form("#sigma_{Gaus} : %.2f #pm %.3f", this_Gaus_sigma, this_Gaus_sigma_err), "");
    l -> AddEntry(in, Form("Par2 : %.2f #pm %.2f", this_par2, this_par2_err), "");
    l -> AddEntry(in, Form("#chi^{2} / ndf : %.2f", chisqr / ndf), "");
  
    fitting_results[plane + "res_range"].push_back(this_zenith);
    fitting_results[plane + "MPV"].push_back(this_MPV);
    fitting_results[plane + "MPV_err"].push_back(this_MPV_err);
    fitting_results[plane + "sigmaL"].push_back(this_Landau_sigma);
    fitting_results[plane + "sigmaL_err"].push_back(this_Landau_sigma_err);
    fitting_results[plane + "sigmaG"].push_back(this_Gaus_sigma);
    fitting_results[plane + "sigmaG_err"].push_back(this_Gaus_sigma_err);
  }
  else{
    //return;
  }

  if(!std::isnan(mean)){
    fitting_results[plane + suffix + side_str + "zenith"].push_back(this_zenith);
    fitting_results[plane + suffix + side_str + "zenith_err"].push_back(this_zenith_err);
    fitting_results[plane + suffix + side_str + "trunmean"].push_back(mean);
    fitting_results[plane + suffix + side_str + "trunmean_err"].push_back(mean_err);
  }
  
  TLine *L_trun_min = new TLine(trun_mean_min, 0., trun_mean_min, max_y * 1.5);
  L_trun_min -> SetLineColor(kRed);
  L_trun_min -> SetLineStyle(7);
  L_trun_min -> Draw("lsame");

  TLine *L_trun_max = new TLine(trun_mean_max, 0., trun_mean_max, max_y * 1.5);
  L_trun_max -> SetLineColor(kRed);
  L_trun_max -> SetLineStyle(7);
  L_trun_max ->	Draw("lsame");

  l -> AddEntry(L_trun_min, "Truncated Mean Range", "l");
  l -> Draw("same");

  
  TString particle_label_str = "";
  if(particle == "muon") particle_label_str = "Cathode-Crossing Through-Going Tracks";
  if(particle == "proton") particle_label_str = "Stopping Proton Candidates";

  TString from_dir = "";
  if(outname.Contains("east")) from_dir = "East";
  if(outname.Contains("west")) from_dir = "West";
      
  TLatex latex_ProtoDUNE, latex_particle, latex_Nhits, latex_method;
  latex_ProtoDUNE.SetNDC();
  latex_particle.SetNDC();
  latex_Nhits.SetNDC();
  latex_method.SetNDC();
  latex_particle.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_particle.SetTextSize(0.03);
  latex_Nhits.SetTextSize(0.05);
  latex_method.SetTextSize(0.05);
  if(isdata) latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND Data} Run " + run_str + ", " + plane);
  else latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND Simulation} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_particle.DrawLatex(0.95, 0.96, particle_label_str);
  latex_method.DrawLatex(0.18, 0.87, latex_str);
  latex_Nhits.DrawLatex(0.18, 0.82, "From " + from_dir);
  
  TString output_plot_dir = getenv("PLOT_PATH");
  TString outfile_str = output_plot_dir + "/MC/dqdx_unif/1d/" + outname + "_" + suffix + ".pdf";
  if(isdata) outfile_str = output_plot_dir + "/Run" + run_str + "/dqdx_unif/1d/"+ outname + "_" + suffix + ".pdf";
  c -> SaveAs(outfile_str);
  c -> Close();

}

void Draw_dqdx_comp(TString input_file_name, TString plane, TString suffix, double rebin_y = 1, double trun_min = 800., double trun_max = 1400.){

  trun_mean_min = trun_min;
  trun_mean_max = trun_max;

  TString histname_east = "zenith_vs_dqdx_from_east_" + suffix + "_" + plane + "_passing_cathode_trklen_400";
  TString histname_west = "zenith_vs_dqdx_from_west_" + suffix + "_" + plane + "_passing_cathode_trklen_400";

  cout << "calling\n" << histname_east << "\n" << histname_west << endl;
  
  TString input_file_dir = getenv("OUTPUTROOT_PATH");
  TFile *f = new TFile(input_file_dir + "/" + input_file_name);
  TH2D *hist_2D_east = (TH2D*)gDirectory -> Get(histname_east);
  TH2D *hist_2D_west = (TH2D*)gDirectory -> Get(histname_west);

  hist_2D_east -> RebinX(5);
  hist_2D_west -> RebinX(5);
  
  hist_2D_east -> RebinY(rebin_y);
  hist_2D_west -> RebinY(rebin_y);

  cout << "Rebinned" << endl;
  
  int N_binsX = hist_2D_east -> GetNbinsX();
  int N_binsY = hist_2D_west -> GetNbinsY();
  cout << "N_binsX : " << N_binsX << endl;
  cout << "N_binsY : " << N_binsY << endl;

  for(int i = 1; i < N_binsX + 1; i++){
    TString i_str = Form("%d", i);
    double this_zenith = hist_2D_east -> GetXaxis() -> GetBinCenter(i);
    double this_zenith_err = 0.5 * hist_2D_east -> GetXaxis() -> GetBinWidth(i);
    TString zenith_range_str = Form("zenith%.1fto%.1fdeg", this_zenith - this_zenith_err, this_zenith + this_zenith_err);
    TString zenith_range_latex = Form("Zenith angle : %.1f - %.1f deg", this_zenith -this_zenith_err, this_zenith + this_zenith_err);
    TString this_hist_name_east = zenith_range_str + "_east";
    TString this_hist_name_west = zenith_range_str + "_west";
    
    TH1D * this_1D_east = new TH1D(this_hist_name_east, this_hist_name_east, N_binsY, 0., 3000.);
    TH1D * this_1D_west = new TH1D(this_hist_name_west, this_hist_name_west, N_binsY, 0., 3000.);
    for(int j = 1; j < N_binsY + 1; j++){
      double this_content_east = hist_2D_east -> GetBinContent(i, j);
      double this_error_east = hist_2D_east -> GetBinError(i, j);
      this_1D_east -> SetBinContent(j, this_content_east);
      this_1D_east -> SetBinError(j, this_error_east);

      double this_content_west = hist_2D_west -> GetBinContent(i, j);
      double this_error_west = hist_2D_west -> GetBinError(i, j);
      this_1D_west -> SetBinContent(j, this_content_west);
      this_1D_west -> SetBinError(j, this_error_west);
    }

    Write_1D_hist(this_1D_east, "/" + plane + "_dqdx_" + this_hist_name_east, suffix, "muon", zenith_range_latex, plane, "dQ/dx [ADC/cm]", "Events", 0.0, 3000.0, false, this_zenith, this_zenith_err);
    Write_1D_hist(this_1D_west, "/" + plane + "_dqdx_" + this_hist_name_west, suffix, "muon", zenith_range_latex, plane, "dQ/dx [ADC/cm]", "Events", 0.0, 3000.0, false, this_zenith, this_zenith_err);
  }
}

void make_gr_and_save(){

  TString output_file_dir = getenv("OUTPUTROOT_PATH");
  TString out_root_name = output_file_dir + "/dqdx_unif_" + run_str + "_selected1or2.root";
  TFile *outfile = new TFile(out_root_name, "RECREATE");
  outfile -> cd();
  
  TString planes_str[] = {"plane0", "plane1", "plane2"};
  TString side_str[] = {"east", "west"};
  TString suffixes[] = {"TPCcore", "NE", "NW"};
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 2; j++){
      for(int k = 0; k < 3; k++){
	cout << planes_str[i] << endl;
	TString this_plane = planes_str[i];
	TString this_side = side_str[j];
	TString this_suffix = suffixes[k];
	cout << fitting_results[this_plane + this_side + "zenith"].size() << endl;
	TGraphErrors *gr_zenith_mean_dqdx = new TGraphErrors(fitting_results[this_plane + this_suffix + this_side + "zenith"].size(),
							     &fitting_results[this_plane + this_suffix + this_side + "zenith"][0], &fitting_results[this_plane + this_suffix + this_side + "trunmean"][0],
							     &fitting_results[this_plane + this_suffix + this_side + "zenith_err"][0], &fitting_results[this_plane + this_suffix + this_side + "trunmean_err"][0]);
	gr_zenith_mean_dqdx -> SetName(this_plane + "_" + this_side + "_" + this_suffix);
	gr_zenith_mean_dqdx -> Write();
      }
    }
  }
  
  outfile -> Close();
}


void run_dqdx_uniform_plots(int run_num = 0){

  if(run_num != 0){
    isdata = true;
    run_str = TString::Format("%d", run_num);
    if(run_num == -1) run_str = "data";
  }
  
  setTDRStyle();

  TString filename = "output_dqdx_uniform_" + run_str + ".root";
  if(!isdata){
    filename = "output_dqdx_uniform_2023B_GENIE_CV.root";
    run_str = "MC";
  }

  TString suffixes[] = {"TPCcore", "NE", "NW"};
  for(int i = 0; i < 3; i++){
    TString this_suffix = suffixes[i];
    Draw_dqdx_comp(filename, "plane0", this_suffix, 10, 800., 1400.);
    Draw_dqdx_comp(filename, "plane1", this_suffix, 10, 800., 1600.);
    Draw_dqdx_comp(filename, "plane2", this_suffix, 10, 900., 1300.);
  }
  make_gr_and_save();
}
