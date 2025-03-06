#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "mylib.h"
#include "Math/Vector3D.h"
#include "BetheBloch.h"

bool isdata = false;
TSpline3 * muon_sp_range_to_KE = Get_sp_range_KE(mass_muon);

double zprime_60deg(double y, double z, int pm = 1){
  double cos_60deg = 0.5;
  double sig_60deg = sqrt(3.) * 0.5;
  double zprime = z * cos_60deg - (pm + 0.) * sig_60deg * y;
  return zprime;
}

double Get_EndMediandQdx(const TTreeReaderArray<float>& rr, const TTreeReaderArray<float>& dqdx){

  vector<double> endp_dqdx;
  for (unsigned i = 0; i < dqdx.GetSize(); i++) {
    if(rr[i] > 0 && rr[i] < 5.){
      endp_dqdx.push_back(dqdx[i]);
    }
  }

  double med_dqdx = -1.;
  if (endp_dqdx.size()) {
    unsigned middle = endp_dqdx.size() / 2;
    std::nth_element(endp_dqdx.begin(), endp_dqdx.begin() + middle, endp_dqdx.end());
    med_dqdx = endp_dqdx[middle];

    // for even case
    if (endp_dqdx.size() % 2 == 0) {
      unsigned other_middle = middle - 1;
      std::nth_element(endp_dqdx.begin(), endp_dqdx.begin() + other_middle, endp_dqdx.end());
      med_dqdx = (med_dqdx + endp_dqdx[other_middle]) / 2.;
    }
  }

  return med_dqdx;
}

void Fill_track_plots(TString suffix, double dist_start, double dist_end, const TTreeReaderArray<float>& rr, const TTreeReaderArray<float>& dqdx){

  FillHist("dist_start_" + suffix, dist_start, 1., 1000., 0., 1000.);
  FillHist("dist_end_" + suffix, dist_end, 1., 1000., 0., 1000.);
  for (unsigned i = 0; i < dqdx.GetSize(); i++) {
    FillHist("rr_vs_dqdx_" + suffix, rr[i], dqdx[i], 1., 300., 0., 300., 3000., 0., 3000.);
  }
}

void Fill_hit_position_plots(TString suffix, const TTreeReaderArray<float>& sp_x, const TTreeReaderArray<float>& sp_y, const TTreeReaderArray<float>& sp_z){

  for(unsigned i = 0; i < sp_x.GetSize(); i++) {
    FillHist("x_vs_y_" + suffix, sp_x[i], sp_y[i], 1., 500., -250., 250., 500., -250., 250.);

    FillHist("z_vs_y_" + suffix, sp_z[i], sp_y[i], 1., 600., -50., 550., 500., -250., 250.);
    if(sp_x[i] > 0) FillHist("z_vs_y_west_" + suffix, sp_z[i], sp_y[i], 1., 600., -50., 550., 500., -250., 250.);
    else FillHist("z_vs_y_east_" + suffix, sp_z[i], sp_y[i], 1., 600., -50., 550., 500., -250., 250.);

    FillHist("x_vs_z_" + suffix, sp_x[i], sp_z[i], 1., 500., -250., 250., 600., -50., 550.);

    if(fabs(sp_x[i]) < 2.){
      FillHist("z_vs_y_smallx_" + suffix, sp_z[i], sp_y[i], 1., 600., -50., 550., 500., -250., 250.);
    }
  }
}

void Fill_hit_dqdx_angular_dependency_plots(TString suffix, const TTreeReaderArray<float>& sp_x, const TTreeReaderArray<float>& sp_y, const TTreeReaderArray<float>& sp_z, const TTreeReaderArray<float>& dqdx, double trk_vec_X, double trk_vec_Y, double trk_vec_Z, double cos_plus_zprimex, double cos_minus_zprimex){

  double cos_zenith = fabs(trk_vec_Y / sqrt( pow(trk_vec_X, 2.) + pow(trk_vec_Y, 2.) + pow(trk_vec_Z, 2.)));
  double zenith = acos(cos_zenith) * 180. / TMath::Pi();

  // == azimuthal angle with top down assumption
  double cos_azimuth = trk_vec_X / sqrt(pow(trk_vec_X, 2.) + pow(trk_vec_Z, 2.)); // == azimuth angle from x-axis
  double azimuth = acos(cos_azimuth) * 180. / TMath::Pi();
  if(trk_vec_Z < 0){
    azimuth = -1. * azimuth;
  }
  
  TString coming_from_str = "";
  if(trk_vec_X * trk_vec_Y > 0.) coming_from_str = "west";
  else coming_from_str = "east";

  bool not_filled_TPCcore = true;
  bool not_filled_NE = true;
  bool not_filled_NW = true;
  for (unsigned i = 0; i < dqdx.GetSize(); i++) {
    if(suffix.Contains("plane0")){
      if(sp_x[i] < 0. && fabs(cos_plus_zprimex) > 0.75) continue;
      if(sp_x[i] > 0. && fabs(cos_minus_zprimex) > 0.75) continue;
    }
    if(suffix.Contains("plane1")){
      if(sp_x[i] < 0. && fabs(cos_minus_zprimex) > 0.75) continue;
      if(sp_x[i] > 0. && fabs(cos_plus_zprimex) > 0.75) continue;
    }
    
    if(fabs(sp_x[i]) < 50 && fabs(sp_y[i]) < 50 && fabs(sp_z[i] - 250.) < 50.){ // consider only core part of TPC
      if(not_filled_TPCcore){
	FillHist("cos_zenith_TPCcore_" + suffix, cos_zenith, 1., 100., 0., 1.);
	FillHist("zenith_TPCcore_" + suffix, zenith, 1., 100., 0., 100.);
	FillHist("azimuth_TPCcore_" + suffix, azimuth, 1., 400., -200., 200.);
	not_filled_TPCcore = false;
      }
      FillHist("cos_zenith_vs_dqdx_TPCcore_" + suffix, cos_zenith, dqdx[i], 1., 100., 0., 1., 3000., 0., 3000.);
      FillHist("cos_zenith_vs_dqdx_from_" + coming_from_str + "_TPCcore_" + suffix, cos_zenith, dqdx[i], 1., 100., 0., 1., 3000., 0., 3000.);
      
      FillHist("zenith_vs_dqdx_TPCcore_" + suffix, zenith, dqdx[i], 1., 100., 0., 100., 3000., 0., 3000.);
      FillHist("zenith_vs_dqdx_from_" + coming_from_str + "_TPCcore_" + suffix, zenith, dqdx[i], 1., 100., 0., 100., 3000., 0., 3000.);
    }
    
    // == NW
    if(sp_z[i] > 250. && sp_x[i] > 0){
      if(not_filled_NW){
        FillHist("cos_zenith_NW_" + suffix, cos_zenith, 1., 100., 0., 1.);
        FillHist("zenith_NW_" + suffix, zenith, 1., 100., 0., 100.);
        FillHist("azimuth_NW_" + suffix, azimuth, 1., 400., -200., 200.);
        not_filled_NW = false;
      }
      FillHist("cos_zenith_vs_dqdx_NW_" + suffix, cos_zenith, dqdx[i], 1., 100., 0., 1., 3000., 0., 3000.);
      FillHist("cos_zenith_vs_dqdx_from_" + coming_from_str + "_NW_" + suffix, cos_zenith, dqdx[i], 1., 100., 0., 1., 3000., 0., 3000.);

      FillHist("zenith_vs_dqdx_NW_" + suffix, zenith, dqdx[i], 1., 100., 0., 100., 3000., 0., 3000.);
      FillHist("zenith_vs_dqdx_from_" + coming_from_str + "_NW_" + suffix, zenith, dqdx[i], 1., 100., 0., 100., 3000., 0., 3000.);
    }

    // == NE
    if(sp_z[i] > 250. && sp_x[i] < 0){
      if(not_filled_NE){
        FillHist("cos_zenith_NE_" + suffix, cos_zenith, 1., 100., 0., 1.);
        FillHist("zenith_NE_" + suffix, zenith, 1., 100., 0., 100.);
        FillHist("azimuth_NE_" + suffix, azimuth, 1., 400., -200., 200.);
        not_filled_NE = false;
      }
      FillHist("cos_zenith_vs_dqdx_NE_" + suffix, cos_zenith, dqdx[i], 1., 100., 0., 1., 3000., 0., 3000.);
      FillHist("cos_zenith_vs_dqdx_from_" + coming_from_str + "_NE_" + suffix, cos_zenith, dqdx[i], 1., 100., 0., 1., 3000., 0., 3000.);

      FillHist("zenith_vs_dqdx_NE_" + suffix, zenith, dqdx[i], 1., 100., 0., 100., 3000., 0., 3000.);
      FillHist("zenith_vs_dqdx_from_" + coming_from_str + "_NE_" + suffix, zenith, dqdx[i], 1., 100., 0., 100., 3000., 0., 3000.);
    }
  }
}

void Fill_hit_plots(TString suffix, const TTreeReaderArray<float>& rr, const TTreeReaderArray<float>& dqdx, const TTreeReaderArray<float>& sp_x, double cos_xy, double cos_yz, double cos_zx, double cos_plus_zprimex, double cos_minus_zprimex){

  for (unsigned i = 0; i < dqdx.GetSize(); i++) {
    FillHist("cos_xy_vs_dqdx_" + suffix, cos_xy, dqdx[i], 1., 220., -1.1, 1.1, 3000., 0., 3000.);
    FillHist("cos_yz_vs_dqdx_" + suffix, cos_yz, dqdx[i], 1., 220., -1.1, 1.1, 3000., 0., 3000.);
    FillHist("cos_zx_vs_dqdx_" + suffix, cos_zx, dqdx[i], 1., 220., -1.1, 1.1, 3000., 0., 3000.);
    FillHist("cos_plus_zprimex_vs_dqdx_" + suffix, cos_plus_zprimex, dqdx[i], 1., 220., -1.1, 1.1, 3000., 0., 3000.);
    FillHist("cos_minus_zprimex_vs_dqdx_" + suffix, cos_minus_zprimex, dqdx[i], 1., 220., -1.1, 1.1, 3000., 0., 3000.);
    if(sp_x[i] < 0.){
      FillHist("cos_plus_zprimex_vs_dqdx_" + suffix + "_East", cos_plus_zprimex, dqdx[i], 1., 220., -1.1, 1.1, 3000., 0., 3000.);
      FillHist("cos_minus_zprimex_vs_dqdx_" + suffix + "_East", cos_minus_zprimex, dqdx[i], 1., 220., -1.1, 1.1, 3000., 0., 3000.);
    }
    else{
      FillHist("cos_plus_zprimex_vs_dqdx_" + suffix + "_West", cos_plus_zprimex, dqdx[i], 1., 220., -1.1, 1.1, 3000., 0., 3000.);
      FillHist("cos_minus_zprimex_vs_dqdx_" + suffix + "_West", cos_minus_zprimex, dqdx[i], 1., 220., -1.1, 1.1, 3000., 0., 3000.);
    }
  }
}

void Fill_corrected_dqdx_plots(TString suffix, const TTreeReaderArray<float>& rr, const TTreeReaderArray<float>& dqdx, const TTreeReaderArray<float>& sp_x, const TTreeReaderArray<float>& sp_z, const TTreeReaderArray<float>& pitch, double last_x, double cos_plus_zprimex, double cos_minus_zprimex, bool do_ind_ang_cut = false){

  //for (unsigned i = 0; i < dqdx.GetSize(); i++) {
  if(dqdx.GetSize() < 1) return;
  for (unsigned i = 1; i < dqdx.GetSize() - 1; i++) {
    if(do_ind_ang_cut){
      if(suffix.Contains("plane0")){
        if(sp_x[i] < 0. && fabs(cos_plus_zprimex) > 0.75) continue;
        if(sp_x[i] > 0. && fabs(cos_minus_zprimex) > 0.75) continue;
      }
      if(suffix.Contains("plane1")){
        if(sp_x[i] < 0. && fabs(cos_minus_zprimex) > 0.75) continue;
        if(sp_x[i] > 0. && fabs(cos_plus_zprimex) > 0.75) continue;
      }
    }

    double this_lifetime_corr = Lifetime_Correction(sp_x[i], 10.0);
    if(isdata) this_lifetime_corr = 1.;
    double corrected_dqdx = dqdx[i] * this_lifetime_corr;
    FillHist("rr_vs_corr_dqdx_" + suffix, rr[i], corrected_dqdx, 1., 300., 0., 300., 3000., 0., 3000.);
    FillHist("rr_vs_pitch_" + suffix, rr[i], pitch[i], 1., 300., 0., 300., 200., 0., 2.);
    FillHist("pitch_" + suffix, pitch[i], 1., 200., 0., 2.);
    FillHist("pitch_x_vs_corr_dqdx_" + suffix, pitch[i], corrected_dqdx, 1., 200., 0., 2., 3000., 0., 3000.);
    if(i == 0) FillHist("last_x_" + suffix, last_x, 1., 500., -250., 250.);
    
    double this_KE= muon_sp_range_to_KE -> Eval(rr[i]); // == from rr

    double gamma = (this_KE/mass_muon)+1.0;
    double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
    double this_xi = Landau_xi(this_KE, pitch[i], mass_muon);
    double this_Wmax = Get_Wmax(this_KE, mass_muon);
    double this_kappa = this_xi / this_Wmax;
    double this_dEdx_BB = meandEdx(this_KE, mass_muon);
    double par[5] = {this_kappa, beta * beta, this_xi, this_dEdx_BB, pitch[i]};
    TF1 * this_dEdx_PDF = dEdx_PDF(par);
    double this_dEdx_MPV = this_dEdx_PDF -> GetMaximumX();

    if(rr[i] < 0.) continue;

    FillHist("dEdx_MPV_vs_corr_dqdx_" + suffix, this_dEdx_MPV, corrected_dqdx, 1., 3000., 0., 30., 3000., 0., 3000.);
    //cout << "[Fill_corrected_dqdx_plots] " << i << ", rr : " << rr[i] << ", KE : " << muon_sp_range_to_KE -> Eval(rr[i]) << ", this_kappa : " << this_kappa << ", this_dEdx_MPV : " << this_dEdx_MPV << endl;

    // == Divide into NE, NW, SE and SW
    if(sp_x[i] < 0.){
      if(sp_z[i] > 250.) FillHist("dEdx_MPV_vs_corr_dqdx_" + suffix + "_NE", this_dEdx_MPV, corrected_dqdx, 1., 3000., 0., 30., 3000., 0., 3000.);
      else FillHist("dEdx_MPV_vs_corr_dqdx_" + suffix + "_SE", this_dEdx_MPV, corrected_dqdx, 1., 3000., 0., 30., 3000., 0., 3000.);
    }
    else{
      if(sp_z[i] > 250.) FillHist("dEdx_MPV_vs_corr_dqdx_" + suffix + "_NW", this_dEdx_MPV, corrected_dqdx, 1., 3000., 0., 30., 3000., 0., 3000.);
      else FillHist("dEdx_MPV_vs_corr_dqdx_" + suffix + "_SW", this_dEdx_MPV, corrected_dqdx, 1., 3000., 0., 30., 3000., 0., 3000.);
    }
    
  }
}
/*
double Vav_MPV(float rr){

  //mass_muon
  double this_xi = dEdx.Get_Landau_xi(KE, width, mass);
  double this_Wmax = dEdx.Get_Wmax(KE, mass);
  double this_kappa = this_xi / this_Wmax;
  double this_dEdx_BB = dEdx.dEdx_Bethe_Bloch(KE, mass);
}
*/
void run_dqdx_uniform_study(int run_num = 0) {

  TString run_str = "";
  if(run_num != 0){
    isdata = true;
    run_str = TString::Format("%d", run_num);
  }
  
  /////////////////////////////////
  // == Define histograms
  /////////////////////////////////
  // == Histograms for overal events
  TH1F *hist_selected = new TH1F("selected", "selected", 3., -0.5, 2.5);

  /////////////////////////////////
  // == Call Trees
  /////////////////////////////////
  // Open the file containing the tree
  TChain *fChain = new TChain("caloskim/TrackCaloSkim");
  TString input_file_dir = getenv("DATA_PATH");
  //TString fileListPath = input_file_dir + "/sample_list/list_MCP2023B_corsika_1Dsim_1Dreco.txt";
  TString fileListPath = input_file_dir + "/sample_list/list_run_" + run_str + "_local.txt";
  if(!isdata) fileListPath = input_file_dir + "/sample_list/list_2023B_GENIE_CV_local.txt";
  cout << "Opening : " << fileListPath << endl;
  AddFilesToChain(fileListPath, fChain);

  TTreeReader myReader(fChain);

  // == Variables
  TTreeReaderValue<int> run(myReader, "trk.meta.run");
  TTreeReaderValue<int> evt(myReader, "trk.meta.evt");
  TTreeReaderValue<unsigned long> ts(myReader, "trk.meta.time");
  TTreeReaderValue<int> trkid(myReader, "trk.id");
  TTreeReaderValue<float> trklen(myReader, "trk.length");
  
  TTreeReaderValue<int> selected(myReader, "trk.selected");
  TTreeReaderArray<float> dqdx0(myReader, "trk.hits0.dqdx"); // hits on plane 0 (Induction)
  TTreeReaderArray<float> dqdx1(myReader, "trk.hits1.dqdx"); // hits on plane 1 (Induction)
  TTreeReaderArray<float> dqdx2(myReader, "trk.hits2.dqdx"); // hits on plane 2 (Collection)
  TTreeReaderArray<float> rr0(myReader, "trk.hits0.rr");
  TTreeReaderArray<float> rr1(myReader, "trk.hits1.rr");
  TTreeReaderArray<float> rr2(myReader, "trk.hits2.rr");
  TTreeReaderArray<float> pitch0(myReader, "trk.hits0.pitch");
  TTreeReaderArray<float> pitch1(myReader, "trk.hits1.pitch");
  TTreeReaderArray<float> pitch2(myReader, "trk.hits2.pitch");
  TTreeReaderArray<float> time0(myReader, "trk.hits0.h.time"); 
  TTreeReaderArray<float> time1(myReader, "trk.hits1.h.time");
  TTreeReaderArray<float> time2(myReader, "trk.hits2.h.time");
  TTreeReaderArray<float> sp_x0(myReader, "trk.hits0.h.sp.x");
  TTreeReaderArray<float> sp_y0(myReader, "trk.hits0.h.sp.y");
  TTreeReaderArray<float> sp_z0(myReader, "trk.hits0.h.sp.z");
  TTreeReaderArray<float> sp_x1(myReader, "trk.hits1.h.sp.x");
  TTreeReaderArray<float> sp_y1(myReader, "trk.hits1.h.sp.y");
  TTreeReaderArray<float> sp_z1(myReader, "trk.hits1.h.sp.z");
  TTreeReaderArray<float> sp_x2(myReader, "trk.hits2.h.sp.x");
  TTreeReaderArray<float> sp_y2(myReader, "trk.hits2.h.sp.y");
  TTreeReaderArray<float> sp_z2(myReader, "trk.hits2.h.sp.z");
  TTreeReaderArray<float> dir_x(myReader, "trk.dir.x"); 
  TTreeReaderArray<float> dir_y(myReader, "trk.dir.y");
  TTreeReaderArray<float> dir_z(myReader, "trk.dir.z");

  TTreeReaderArray<float> true_start_x(myReader, "trk.truth.p.start.x");
  TTreeReaderArray<float> true_start_y(myReader, "trk.truth.p.start.y");
  TTreeReaderArray<float> true_start_z(myReader, "trk.truth.p.start.z");
  TTreeReaderArray<float> true_end_x(myReader, "trk.truth.p.end.x");
  TTreeReaderArray<float> true_end_y(myReader, "trk.truth.p.end.y");
  TTreeReaderArray<float> true_end_z(myReader, "trk.truth.p.end.z");
  TTreeReaderArray<int> true_end_process(myReader, "trk.truth.p.end_process");
  TTreeReaderArray<float> true_hit_time(myReader, "trk.truth.p.truehits2.time");
  TTreeReaderArray<float> true_hit_tdrift(myReader, "trk.truth.p.truehits2.tdrift");


  /////////////////////////////////
  // == Loop for tracks
  /////////////////////////////////
  int N_entries = myReader.GetEntries();
  cout << "N_entries : " << N_entries << endl;
  int current_entry = 0;

  int N_run = 50000;
  double ADC_med_cut = 1600.; // == https://sbn-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=23472&filename=SBND%20Calib%20Workshop%202021.pdf&version=1
  double track_length_cut = 15.;
  // Loop over all entries of the TTree
  while (myReader.Next()) {
    if(current_entry > N_run) break;
   
    if(current_entry%100 == 0){
      cout << current_entry << " / " << N_entries << endl;
    }
    current_entry++;

    if(*run == 14480){
      if(*evt >= 748) continue;
    }
    if(*run == 14608){
      if(*evt >= 9695) continue;
    }
    
    hist_selected -> Fill(*selected);

    double end_meddqdx = Get_EndMediandQdx(rr2, dqdx2);
    FillHist("end_meddqdx", end_meddqdx, 1., 5000., 0., 5000.);
    FillHist(Form("end_meddqdx_selected%d", *selected), end_meddqdx, 1., 5000., 0., 5000.);
    // == Tracks selected as through-going
    if (*selected == 2 || *selected == 1) {
      
      unsigned N_reco_hits = rr2.GetSize();
      //cout << "N_reco_hits : " << N_reco_hits << endl;
      if(N_reco_hits < 2) continue;
      
      TVector3 this_reco_start(sp_x2[N_reco_hits - 1], sp_y2[N_reco_hits - 1], sp_z2[N_reco_hits - 1]);
      TVector3 this_reco_end(sp_x2[0], sp_y2[0], sp_z2[0]);
      TVector3 this_true_start(true_start_x[0], true_start_y[0], true_start_z[0]);
      TVector3 this_true_end(true_end_x[0], true_end_y[0], true_end_z[0]);

      double dist_start = (this_reco_start - this_true_start).Mag();
      double dist_end = (this_reco_end - this_true_end).Mag();

      double this_reco_trk_len = rr2[N_reco_hits - 1];

      // == Pick first and last hits to check if the track is passing the cathode
      double first_x = -999.;
      double first_y = -999.;
      double first_z = -999.;
      double last_x = -999.;
      double last_y = -999.;
      double last_z = -999.;

      first_x = sp_x2[N_reco_hits - 1];
      first_y = sp_y2[N_reco_hits - 1];
      first_z = sp_z2[N_reco_hits - 1];
      for (unsigned i = 0; i < dqdx2.GetSize(); i++) {
        if(rr2[i] > 0){
	  last_x = sp_x2[i];
	  last_y = sp_y2[i];
          last_z = sp_z2[i];
	  break;
	}	    
      }
      bool passing_cathode = false;
      ROOT::Math::XYZVector track_vec(last_x - first_x, last_y - first_y, last_z - first_z);
      double cos_xy = track_vec.X() / (sqrt(pow(track_vec.X(), 2.) + pow(track_vec.Y(), 2.)));
      double cos_yz = track_vec.Y() / (sqrt(pow(track_vec.Y(), 2.) + pow(track_vec.Z(), 2.)));
      double cos_zx = track_vec.X() / (sqrt(pow(track_vec.X(), 2.) + pow(track_vec.Z(), 2.)));
      double zprime_plus = zprime_60deg(track_vec.Y(), track_vec.Z(), 1);
      double zprime_minus = zprime_60deg(track_vec.Y(), track_vec.Z(), -1);
      double cos_plus_zprimex = track_vec.X() / (sqrt(pow(track_vec.X(), 2.) + pow(zprime_plus, 2.)));
      double cos_minus_zprimex = track_vec.X() / (sqrt(pow(track_vec.X(), 2.) + pow(zprime_minus, 2.)));
      
      if(first_x * last_x < 0.) passing_cathode= true;

      if(!passing_cathode) continue;
      FillHist("reco_trk_len_passing_cathode", this_reco_trk_len, 1., 1000., 0., 1000.);
      
      if(this_reco_trk_len < 400.) continue;
      FillHist("reco_trk_len_passing_cathode_trklen_400", this_reco_trk_len, 1., 1000., 0., 1000.);
      Fill_hit_position_plots("plane0_passing_cathode_trklen_400", sp_x0, sp_y0, sp_z0);
      Fill_hit_position_plots("plane1_passing_cathode_trklen_400", sp_x1, sp_y1, sp_z1);
      Fill_hit_position_plots("plane2_passing_cathode_trklen_400", sp_x2, sp_y2, sp_z2);

      Fill_hit_dqdx_angular_dependency_plots("plane0_passing_cathode_trklen_400", sp_x0, sp_y0, sp_z0, dqdx0, track_vec.X(), track_vec.Y(), track_vec.Z(), cos_plus_zprimex, cos_minus_zprimex);
      Fill_hit_dqdx_angular_dependency_plots("plane1_passing_cathode_trklen_400", sp_x1, sp_y1, sp_z1, dqdx1, track_vec.X(), track_vec.Y(), track_vec.Z(), cos_plus_zprimex, cos_minus_zprimex);
      if(fabs(cos_zx) < 0.75) Fill_hit_dqdx_angular_dependency_plots("plane2_passing_cathode_trklen_400", sp_x2, sp_y2, sp_z2, dqdx2, track_vec.X(), track_vec.Y(), track_vec.Z(), cos_plus_zprimex, cos_minus_zprimex);
    }
  }

  TString output_rootfile_dir = getenv("OUTPUTROOT_PATH");
  TString output_file_name = output_rootfile_dir + "/output_dqdx_uniform_" + run_str + ".root";
  if(!isdata) output_file_name = output_rootfile_dir + "/output_dqdx_uniform_2023B_GENIE_CV.root";
  out_rootfile = new TFile(output_file_name, "RECREATE");
  out_rootfile -> cd();
  
  hist_selected -> Write();
  WriteHist();

  out_rootfile -> Close();
}
