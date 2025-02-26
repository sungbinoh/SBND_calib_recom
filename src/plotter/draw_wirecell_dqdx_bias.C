#include "canvas_margin.h"

void draw_wirecell_dqdx_bias(){
  setTDRStyle();
  gStyle->SetLineWidth(3);
  // Open the file
  std::ifstream file("./data/wirecell_dqdx_perf/muon_sp_bias.csv");
  if (!file.is_open()) {
    std::cerr << "Error: Cannot open the file!" << std::endl;
    return;
  }
  
  // Variables to store the data
  std::vector<double> angle_bin_start, angle_bin_stop, angle_bin_center, u_med, v_med, w_med;
  std::string line;
  
  // Skip the first line (header)
  std::getline(file, line);

  // Read data
  double a_start, a_stop, u, v, w;
  string elline;
  while(getline(file,elline)){
    std::istringstream is( elline );
    TString tstring_elline = elline;
    std::string token;
    cout << "[draw_wirecell_dqdx_bias] tstring_elline : " << tstring_elline << endl;
    if(tstring_elline.Contains("angle_bin_start")) continue;
    
    std::getline(is, token, ','); a_start = std::stod(token);
    std::getline(is, token, ','); a_stop = std::stod(token);
    std::getline(is, token, ','); u = std::stod(token);
    std::getline(is, token, ','); v = std::stod(token);
    std::getline(is, token, ','); w = std::stod(token);
    /*
    is >> a_start;
    is >> a_stop;
    is >> u;
    is >> v;
    is >> w;
    */
    //cout << "a_start: " << a_start << ", a_stop: " << a_stop << ", u: " << u << ", v: " << v << ", w: " << w << endl;

    angle_bin_start.push_back(a_start);
    angle_bin_stop.push_back(a_stop);
    angle_bin_center.push_back((a_start + a_stop) / 2.);
    u_med.push_back(u * 100.);
    v_med.push_back(v * 100.);
    w_med.push_back(w * 100.);
  }
  
  file.close();
  
  // Create a TGraph for visualization
  TCanvas *c1 = new TCanvas("c1", "Muon Bias", 1600, 1200);

  TH1D * template_h = new TH1D("", "", 1., 0, 90);
  template_h -> SetStats(0);
  template_h -> GetYaxis() -> SetRangeUser(-2., 2);
  template_h -> GetXaxis() -> SetTitle("Track angle [#circ]");
  template_h -> GetXaxis() -> SetTitleSize(0.037);
  template_h -> GetXaxis() -> SetTitleOffset(1.4);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("dQ/dx fractional bias [%]");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  //template_h -> SetBinContent(1, -1000.);
  template_h -> SetLineWidth(2);
  template_h -> SetLineStyle(7);
  template_h ->	SetLineColor(kRed);
  template_h -> Draw();

  /*
  TGraph *gU = new TGraph(angle_bin_start.size(), angle_bin_start.data(), u_med.data());
  gU->SetTitle("U Bias vs Angle;Angle Bin Start (deg);U Med");
  gU->SetMarkerStyle(20);
  gU->SetMarkerColor(kRed);
  gU->Draw("psame");
  
  TGraph *gV = new TGraph(angle_bin_start.size(), angle_bin_start.data(), v_med.data());
  gV->SetMarkerStyle(21);
  gV->SetMarkerColor(kBlue);
  gV->Draw("psame");
  */
  
  TGraph *gW = new TGraph(angle_bin_start.size(), angle_bin_center.data(), w_med.data());
  //gW->SetMarkerStyle(22);
  gW->SetMarkerSize(2.0);
  gW->SetMarkerColor(kGreen);
  gW->Draw("psame");

  
  TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]*x + [2]*x*x", 0, 90);
  fitFunc->SetParameters(0, 0.01, -0.01);
  gW->Fit(fitFunc, "RBOSQN");
  fitFunc->SetLineColor(kBlue);
  fitFunc->SetLineWidth(2);
  fitFunc->Draw("same");
  
  // Add legend
  auto legend = new TLegend(0.2, 0.7, 0.7, 0.9);
  /*
  legend->AddEntry(gU, "U Med", "p");
  legend->AddEntry(gV, "V Med", "p");
  */
  legend->AddEntry(gW, "Collection plane", "p");
  legend->AddEntry(fitFunc, "Fit with p_{0} + p_{1}x + p_{2}x^{2}", "l");
  legend->AddEntry(fitFunc, Form("p_{0} = %.2f #pm %.2e", fitFunc -> GetParameter(0), fitFunc -> GetParError(0)), "");
  legend->AddEntry(fitFunc, Form("p_{1} = %.2f #pm %.2e", fitFunc -> GetParameter(1), fitFunc -> GetParError(1)), "");
  legend->AddEntry(fitFunc, Form("p_{2} = %.2e #pm %.2e", fitFunc -> GetParameter(2), fitFunc -> GetParError(2)), "");
  legend->Draw("same");
  
  c1->SaveAs("./output/plots/muon_sp_bias.pdf");

}
