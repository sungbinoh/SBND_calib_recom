#include "canvas_margin.h"
#include "mylib.h"
#include "TRandom3.h"

//TString input_file_name = "fit_lifetime_2025A_SpringDev_data_bnbcosmics.root";
TString input_file_name = "fit_lifetime_2025B_final_valid_data_bnb_light.root";
TString sample_str = "Data (Run 18255 and 18259)";
TString output_str = "run_18255_and_18259";
void Fit_lifetime_overlay(TString side, TString ngroupwires, TString plane, TString recom_model, double x_range_down, double x_range_up, double y_range_down, double y_range_up){

  TString input_file_dir = getenv("OUTPUTROOT_PATH");
  TFile *f = new TFile(input_file_dir + "/lifetime/whicht00/" + input_file_name);

  TString side_latex = "East TPC";
  if(side == "W") side_latex = "West TPC";

  
  TString grname = ngroupwires + "_" + plane + "_" + side;
  double fit_x_low = 0.1;
  double fit_x_high = 1.24;

  vector<TGraphErrors*> gr_vec;
  TGraphErrors* default_gr = (TGraphErrors*)gDirectory -> Get(grname);
  TGraphErrors* sce_space_gr = (TGraphErrors*)gDirectory -> Get(grname + "_sce_corr");
  gr_vec.push_back(default_gr);
  gr_vec.push_back(sce_space_gr);

  TString dedx_assume_str[] = {"1p7", "1p8", "1p9", "2p0", "2p1", "2p2", "2p3"};
  int N_dedx_assume = 7;

  for(int i = 0; i < N_dedx_assume; i++){
    TString this_dedx_assume = dedx_assume_str[i];
    cout << "[Fit_lifetime_overlay] adding " << grname + "_sce_corr_e_distor_" + recom_model + "_dedx" + this_dedx_assume << endl;
    TGraphErrors* this_gr = (TGraphErrors*)gDirectory -> Get(grname + "_sce_corr_e_distor_" + recom_model + "_dedx" + this_dedx_assume) -> Clone();
    cout << "[Fit_lifetime_overlay] added " << this_dedx_assume << endl;
    gr_vec.push_back(this_gr);
  }

  cout << "[Fit_lifetime_overlay] gr_vec.size(): " << gr_vec.size() << endl;
  
  TString lgd_strs[] = {"No SCE Corr.", "Spatial SCE Corr. Only",
                        "dE/dx = 1.7 MeV/cm",
			"dE/dx = 1.8 MeV/cm",
			"dE/dx = 1.9 MeV/cm",
			"dE/dx = 2.0 MeV/cm",
			"dE/dx = 2.1 MeV/cm",
			"dE/dx = 2.2 MeV/cm",
			"dE/dx = 2.3 MeV/cm"
  };

  TString outtxt_strs[] = {"UnCorr", "Spatial",
			   "1.7",
			   "1.8",
			   "1.9",
			   "2.0",
			   "2.1",
			   "2.2",
			   "2.3"
  };

  //int color_arr[] = {1, 920, 632, 800, 401, 418, 432, 600, 880};
  int color_arr[] = {920, 1, 632, 800, 401, 418, 432, 600, 880};

  TString recom_model_lgd_str = "";
  if(recom_model == "emb") recom_model_lgd_str = "#splitline{Ellipsoidal Modified Box Model}{(ICARUS parameters)}";
  else recom_model_lgd_str = "#splitline{Modified Box Model}{(ArgoNeuT parameters)}";
  
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

  vector<TF1*> fit_grs;

  TLegend *l = new TLegend(0.45, 0.82, 0.92, 0.92);
  TLegend *l2 = new TLegend(0.45, 0.55, 0.92, 0.80);

  TString output_file_dir = getenv("OUTPUTROOT_PATH");
  TString out_txt_name = output_file_dir + "/lifetime/whicht00/text/lifetime_run_" + output_str + "_" + side + ".txt";
  ofstream outtxtfile(out_txt_name);
  
  vector<double> one_over_tau_vec;
  for(int i = 0; i < N_dedx_assume + 2; i++){

    cout << "[Fit_lifetime_overlay] collect gr.vec.at(i)" << endl;
    TGraphErrors *this_gr = gr_vec.at(i);
    int this_color = color_arr[i];
    this_gr -> SetLineColor(this_color);
    this_gr -> SetMarkerColor(this_color);
    this_gr -> Draw("epsame");

    cout << "[Fit_lifetime_overlay] set colors for this_gr" << endl;
    TF1 *this_gr_fit = new TF1("gr_fit_" + lgd_strs[i], "[0] * exp(-1.  * [1] * x)", fit_x_low, fit_x_high);
    this_gr_fit -> SetParameters(1150., 0.02);
    this_gr -> Fit(this_gr_fit, "RN", "", fit_x_low, fit_x_high);
    this_gr_fit -> SetLineColor(this_color);
    this_gr_fit -> SetMarkerColor(this_color);
    this_gr_fit -> SetMarkerStyle(32);
    this_gr_fit -> SetLineWidth(2);
    this_gr_fit -> SetLineStyle(7);
    this_gr_fit -> Draw("lsame");

    double this_one_over_tau = this_gr_fit -> GetParameter(1);
    double this_one_over_tau_err = this_gr_fit -> GetParError(1);
    double tau_central = 1. / this_one_over_tau;
    double tau_minus_1sig = 1. / (this_one_over_tau + this_one_over_tau_err);
    double tau_plus_1sig = 1. / (this_one_over_tau - this_one_over_tau_err);
    double tau_error_minus = tau_central - tau_minus_1sig;
    double tau_error_plus = tau_plus_1sig - tau_central;

    TString this_lifetime = Form(" (#tau = %.1f^{+%.1f}_{-%1.f} ms)", tau_central, tau_error_plus, tau_error_minus);

    outtxtfile << outtxt_strs[i] << "\t" << tau_central << "\t" << tau_error_plus << "\t" << tau_error_minus << "\n";

    
    if(i < 2) l -> AddEntry(this_gr_fit, lgd_strs[i] + this_lifetime, "pl");
    if(i == 2) l2 ->  AddEntry(this_gr_fit, "Spatial + |E|-distor. SCE corr.", "");
    if(i > 1){
      l2 -> AddEntry(this_gr_fit, lgd_strs[i] + this_lifetime, "pl");
    }
  }

  /*
  gr_vec.at(1) -> Draw("epsame");
  fit_grs.at(1) -> Draw("lsame");
  */

  sce_space_gr -> Draw("epsame");
  
  l -> SetFillColor(0);
  l -> Draw("same");

  l2 -> SetFillColor(0);
  l2 -> Draw("same");


  TLatex latex_ProtoDUNE, latex_particle, latex_recom;
  latex_ProtoDUNE.SetNDC();
  latex_ProtoDUNE.SetTextSize(0.03);
  if(input_file_name.Contains("data") || input_file_name.Contains("run")){
    //latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND " + sample_str + " " + side_latex + "} #font[42]{#it{#scale[0.8]{Preliminary}}}");
    latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND " + sample_str + " " + side_latex + "}");
    TLatex not_for_public;
    not_for_public.SetNDC();
    not_for_public.SetTextSize(0.10);
    not_for_public.DrawLatex(0.20, 0.25, "#color[632]{SBND internal use only}");
  }
  else{
    latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{SBND Simulation} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  }
  latex_particle.SetNDC();
  latex_particle.SetTextSize(0.03);
  latex_particle.SetTextAlign(31);
  latex_particle.DrawLatex(0.95, 0.96, "Anode-Cathode Passing Tracks");

  latex_recom.SetNDC();
  latex_recom.SetTextSize(0.025);
  latex_recom.DrawLatex(0.20, 0.85, recom_model_lgd_str);
  
  TString output_plot_dir = getenv("PLOT_PATH");
  TString outfile_str = output_plot_dir + "/lifetime/whicht00/" + output_str + "/lifetime_wire" + ngroupwires + "_" + plane + "_" + side + "_" + recom_model + "_overlay.pdf";
  TString outfile_dir = gSystem->DirName(outfile_str);
    if (gSystem->AccessPathName(outfile_dir)) {
      std::cout << "Directory does not exist, creating: " << outfile_dir << std::endl;
      gSystem->mkdir(outfile_dir, kTRUE);
    }
    c -> SaveAs(outfile_str);
  c -> Close();
  outfile.close();
}

void draw_lifetime_fit_sce_efield_corr_fits(TString input_file_name_, TString sample_str_, TString output_str_){
  input_file_name = input_file_name_;
  sample_str = sample_str_;
  output_str = output_str_;

  setTDRStyle();

  TString sides[] = {"E", "W"};
  TString planes[] = {"plane0", "plane1", "plane2"};
  for(int i = 0; i < 3; i++){
    TString this_plane = planes[i];

    if(i != 2) continue; // FIXME: running only collection plane
    for(int j = 0; j < 2; j++){
      TString this_side = sides[j];
      //Fit_lifetime_overlay(this_side, "10", this_plane, "mb", 0., 1.3, 1070., 1200.);
      //Fit_lifetime_overlay(this_side, "10", this_plane, "emb", 0., 1.3, 1070., 1200.);

      //Fit_lifetime_overlay(this_side, "10", this_plane, "mb", 0., 1.3, 960., 1300.);
      Fit_lifetime_overlay(this_side, "10", this_plane, "emb", 0., 1.3, 1050., 1300.);
    }
  }

  
}
