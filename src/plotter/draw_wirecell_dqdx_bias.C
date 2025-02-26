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
  std::vector<double> angle_bin_start, angle_bin_stop, u_med, v_med, w_med;
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
    u_med.push_back(u);
    v_med.push_back(v);
    w_med.push_back(w);
  }
  
  file.close();
  
  // Create a TGraph for visualization
  TCanvas *c1 = new TCanvas("c1", "Muon Bias", 800, 600);
  
  TGraph *gU = new TGraph(angle_bin_start.size(), angle_bin_start.data(), u_med.data());
  gU->SetTitle("U Bias vs Angle;Angle Bin Start (deg);U Med");
  gU->SetMarkerStyle(20);
  gU->SetMarkerColor(kRed);
  gU->Draw("AP");
  
  TGraph *gV = new TGraph(angle_bin_start.size(), angle_bin_start.data(), v_med.data());
  gV->SetMarkerStyle(21);
  gV->SetMarkerColor(kBlue);
  gV->Draw("P SAME");
  
  TGraph *gW = new TGraph(angle_bin_start.size(), angle_bin_start.data(), w_med.data());
  gW->SetMarkerStyle(22);
  gW->SetMarkerColor(kGreen);
  gW->Draw("P SAME");
  
  // Add legend
  auto legend = new TLegend(0.1, 0.7, 0.3, 0.9);
  legend->AddEntry(gU, "U Med", "p");
  legend->AddEntry(gV, "V Med", "p");
  legend->AddEntry(gW, "W Med", "p");
  legend->Draw();
  
  c1->SaveAs("muon_sp_bias.png");

}
