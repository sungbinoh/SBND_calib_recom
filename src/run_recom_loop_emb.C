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
BetheBloch *muon_BB = new BetheBloch(13);

double dqdx_scale_correction_angle(double theta){

  // theta shoud be deg
  double p0 = -1.55222;
  double p1 = 0.0500101;
  double p2 = -0.000211176;

  double this_bias = (p0 + p1 * theta + p2 * theta * theta) / 100.; // -- % to number
  double this_correction = 1. / (1. - this_bias);

  return this_correction;
}


double zprime_60deg(double y, double z, int pm = 1){
  // == This function is for angle cut for induction planes.
  // ==== Rotation on yz plane to provide z' for cos(theta_z'x) measurement
  double cos_60deg = 0.5;
  double sig_60deg = sqrt(3.) * 0.5;
  double zprime = z * cos_60deg - (pm + 0.) * sig_60deg * y;
  return zprime;
}

double Get_EndMediandQdx(const TTreeReaderArray<float>& rr, const TTreeReaderArray<float>& dqdx){
  // == This function is for a study with stronger cut on EndMediandQddx for stopping track selection
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

TString Get_theta_trk_x_bin_str(double theta_trk_x){
  // == This is for getting bin string for theta_trk_x, following ICARUS's binning (https://arxiv.org/pdf/2407.12969)
  if(theta_trk_x < 30.) return "0to30";
  else if(theta_trk_x < 80.){
    int tenth = theta_trk_x / 10;
    TString this_out = Form("%dto%d", tenth * 10, (tenth+1) * 10);
    return this_out;
  }
  else if(theta_trk_x < 85.) return "80to85";
  else if(theta_trk_x < 90.) return "85to90";
  else return "Error";

  return "Error";
}

void Fill_track_plots(TString suffix, double dist_start, double dist_end, const TTreeReaderArray<float>& rr, const TTreeReaderArray<float>& dqdx){
  // == Fill plots for varaibles related with track
  FillHist("dist_start_" + suffix, dist_start, 1., 1000., 0., 1000.);
  FillHist("dist_end_" + suffix, dist_end, 1., 1000., 0., 1000.);
  for (unsigned i = 0; i < dqdx.GetSize(); i++) {
    FillHist("rr_vs_dqdx_" + suffix, rr[i], dqdx[i], 1., 300., 0., 300., 5000., 0., 5000.);
  }
}

void Fill_hit_plots(TString suffix, const TTreeReaderArray<float>& rr, const TTreeReaderArray<float>& dqdx, const TTreeReaderArray<float>& sp_x, double cos_xy, double cos_yz, double cos_zx, double cos_plus_zprimex, double cos_minus_zprimex){
  // == Fill plots for variables related anguler cut for track selection
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

void Fill_corrected_dqdx_plots(TString suffix, const TTreeReaderArray<float>& rr, const TTreeReaderArray<float>& dqdx, const TTreeReaderArray<float>& sp_x, const TTreeReaderArray<float>& sp_z, const TTreeReaderArray<float>& pitch, double last_x, double cos_plus_zprimex, double cos_minus_zprimex, double theta_trk_x, TString theta_trk_x_str, bool do_ind_ang_cut = false){
  // == Fill plots for recombination fits
  if(dqdx.GetSize() < 1) return;
  for (unsigned i = 1; i < dqdx.GetSize() - 1; i++) { // == Not using first and last hits of a track
    if(rr[i] < 0.) continue;
    if(do_ind_ang_cut){ // == angle cut for induction planes
      if(suffix.Contains("plane0")){
        if(sp_x[i] < 0. && fabs(cos_plus_zprimex) > 0.75) continue;
        if(sp_x[i] > 0. && fabs(cos_minus_zprimex) > 0.75) continue;
      }
      if(suffix.Contains("plane1")){
        if(sp_x[i] < 0. && fabs(cos_minus_zprimex) > 0.75) continue;
        if(sp_x[i] > 0. && fabs(cos_plus_zprimex) > 0.75) continue;
      }
    }

    double this_lifetime_corr = Lifetime_Correction(sp_x[i], 100.0);
    if(isdata) this_lifetime_corr = 1.; // == FIXME, for data, do not apply lifetime correction. Should be updated in future to use different lifetime values for MC and data
    double corrected_dqdx = dqdx[i] * this_lifetime_corr;
    double this_dqdx_bias_corr = dqdx_scale_correction_angle(theta_trk_x);
    //cout << "corrected_dqdx: " << corrected_dqdx << endl;
    corrected_dqdx *= this_dqdx_bias_corr;
    //cout << "bias corrected_dqdx: " << corrected_dqdx << endl;

    FillHist("rr_vs_corr_dqdx_" + suffix, rr[i], corrected_dqdx, 1., 300., 0., 300., 5000., 0., 5000.);
    FillHist("rr_vs_pitch_" + suffix, rr[i], pitch[i], 1., 300., 0., 300., 200., 0., 2.);
    FillHist("pitch_" + suffix, pitch[i], 1., 200., 0., 2.);
    FillHist("pitch_x_vs_corr_dqdx_" + suffix, pitch[i], corrected_dqdx, 1., 200., 0., 2., 3000., 0., 3000.);
    if(i == 0) FillHist("last_x_" + suffix, last_x, 1., 500., -250., 250.);

    if(pitch[i] > 1.) continue;
    if(rr[i] < 0.) continue;
    
    double this_KE= muon_BB -> KEFromRangeSpline(rr[i]); // == from rr

    double gamma = (this_KE/mass_muon)+1.0;
    double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
    double this_xi = muon_BB -> Landau_xi(this_KE, pitch[i]);
    double this_Wmax = muon_BB -> Get_Wmax(this_KE);
    double this_kappa = this_xi / this_Wmax;
    double this_dEdx_BB = muon_BB -> meandEdx(this_KE);
    double par[5] = {this_kappa, beta * beta, this_xi, this_dEdx_BB, pitch[i]};
    TF1 * this_dEdx_PDF = muon_BB -> dEdx_PDF(this_KE, pitch[i]);
    double this_dEdx_MPV = this_dEdx_PDF -> GetMaximumX();
    delete this_dEdx_PDF;
    
    FillHist("dEdx_MPV_vs_corr_dqdx_" + suffix, this_dEdx_MPV, corrected_dqdx, 1., 3000., 0., 30., 3000., 0., 3000.);
    FillHist("dEdx_MPV_vs_corr_dqdx_" + suffix + "_phi" + theta_trk_x_str, this_dEdx_MPV, corrected_dqdx, 1., 3000., 0., 30., 3000., 0., 3000.);
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

void run_recom_loop_emb(int run_number = 0) {

  TString run_number_str = "";
  run_number_str = TString::Format("%d", run_number);
  if(run_number > 0){
    isdata = true;
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
  TString sample_list_dir = getenv("SAMPLE_PATH");
  TString sample_list_label = getenv("FILELIST_LABEL");

  TString fileListPath = sample_list_dir + "/list" + sample_list_label + run_number_str + ".txt";
  if(!isdata) fileListPath = sample_list_dir + "/list_2024B_MC_calib_ntuple.txt";
  //if(!isdata) fileListPath = sample_list_dir + "/list_2025A_Sprint25Dev_MC_bnbcosmics_calib_ntuple.txt";
  cout << "Opening : " << fileListPath << endl;
  // Check if the file exists
  std::ifstream file(fileListPath.Data());  // Convert TString to const char*
  if (!file) {
    cout << "File does not exist: " << fileListPath << endl;
    cout << "Exiting [run_recom_loop_emb]" << endl;
    return;
  }
      
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

  // Loop over all entries of the TTree
  int _run_to = 300000;
  //_run_to = 1000;
  while (myReader.Next()) {
    if(current_entry > _run_to) break;
   
    if(current_entry%100 == 0){
      cout << current_entry << " / " << N_entries << endl;
    }
    current_entry++;

    FillHist("selected", *selected, 1., 40, -0.5, 39.5);

    double end_meddqdx = Get_EndMediandQdx(rr2, dqdx2);
    FillHist("end_meddqdx", end_meddqdx, 1., 5000., 0., 5000.);
    FillHist(Form("end_meddqdx_selected%d", *selected), end_meddqdx, 1., 5000., 0., 5000.);

    // == Tracks selected as stopping
    if (*selected == 0) {
      
      unsigned N_reco_hits = rr2.GetSize();
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
      ROOT::Math::XYZVector x_unit(1.0, 0., 0.);
      ROOT::Math::XYZVector trk_cross_x = track_vec.Cross(x_unit);
      double sin_theta_trk_x = sqrt(trk_cross_x.Mag2()) / sqrt(track_vec.Mag2());
      double theta_trk_x = TMath::ASin(sin_theta_trk_x) * 180. / TMath::Pi(); // == [Deg.]
      TString theta_trk_x_str = Get_theta_trk_x_bin_str(theta_trk_x);
      //cout << "theta_trk_x : " << theta_trk_x << ", theta_trk_x_str : " << theta_trk_x_str << endl;
      
      if(first_x * last_x < 0.) passing_cathode= true;
      //cout << "passing_cathode : " << passing_cathode << ", this_reco_trk_len : " <<this_reco_trk_len << ", first_x : " << first_x << ", last_x : " << last_x << endl;
      
      // == Nocut
      /*
      FillHist("reco_trk_len_nocut", this_reco_trk_len, 1., 1000., 0., 1000.);
      FillHist("reco_trk_len_vs_passing_cathode_nocut", this_reco_trk_len, passing_cathode, 1., 1000., 0., 1000., 2., -0.5, 1.5);

      FillHist("dist_start_vs_reco_trk_len_nocut", dist_start, this_reco_trk_len, 1., 1000., 0., 1000., 1000., 0., 1000.);
      FillHist("dist_end_vs_reco_trk_len_nocut", dist_end, this_reco_trk_len, 1., 1000., 0., 1000., 1000., 0., 1000.);
      Fill_track_plots("nocut", dist_start, dist_end, rr, dqdx);

      // == Track length 5 cm cut
      if(this_reco_trk_len < 5.) continue;
      Fill_track_plots("trklen_05cm", dist_start, dist_end, rr, dqdx);

      // == Track length 10 cm cut
      if(this_reco_trk_len < 10.) continue;
      Fill_track_plots("trklen10cm", dist_start, dist_end, rr, dqdx);

      // == Track length 15 cm cut
      if(this_reco_trk_len < 15.) continue;
      Fill_track_plots("trklen_15cm", dist_start, dist_end, rr, dqdx);

      // == Track length 20 cm cut
      if(this_reco_trk_len < 20.) continue;
      Fill_track_plots("trklen_20cm", dist_start, dist_end, rr, dqdx);

      // == Track length 30 cm cut
      if(this_reco_trk_len < 30.) continue;
      Fill_track_plots("trklen_30cm", dist_start, dist_end, rr, dqdx);
      */

      // == Track length 60 cm cut
      if(this_reco_trk_len < 60. || !passing_cathode) continue;

      // == These are for checking impact of event selections
      /*
      Fill_track_plots("plane0_trklen_60cm_passing_cathode", dist_start, dist_end, rr0, dqdx0);
      Fill_track_plots("plane1_trklen_60cm_passing_cathode", dist_start, dist_end, rr1, dqdx1);
      Fill_track_plots("plane2_trklen_60cm_passing_cathode", dist_start, dist_end, rr2, dqdx2);
      Fill_hit_plots("plane0_trklen_60cm_passing_cathode", rr0, dqdx0, sp_x0, cos_xy, cos_yz, cos_zx, cos_plus_zprimex, cos_minus_zprimex);
      Fill_hit_plots("plane1_trklen_60cm_passing_cathode", rr1, dqdx1, sp_x1, cos_xy, cos_yz, cos_zx, cos_plus_zprimex, cos_minus_zprimex);
      Fill_hit_plots("plane2_trklen_60cm_passing_cathode", rr2, dqdx2, sp_x2, cos_xy, cos_yz, cos_zx, cos_plus_zprimex, cos_minus_zprimex);

      Fill_corrected_dqdx_plots("plane2_trklen_60cm_passing_cathode", rr2, dqdx2, sp_x2, sp_z2, pitch2, last_x, cos_plus_zprimex, cos_minus_zprimex, theta_trk_x_str, true);
      Fill_corrected_dqdx_plots("plane0_trklen_60cm_passing_cathode", rr0, dqdx0, sp_x0, sp_z0, pitch0, last_x, cos_plus_zprimex, cos_minus_zprimex, theta_trk_x_str, false);
      Fill_corrected_dqdx_plots("plane1_trklen_60cm_passing_cathode", rr1, dqdx1, sp_x1, sp_z1, pitch1, last_x, cos_plus_zprimex, cos_minus_zprimex, theta_trk_x_str, false);
      */

      // == These are the main plots
      Fill_corrected_dqdx_plots("plane0_trklen_60cm_passing_cathode_coszx", rr0, dqdx0, sp_x0, sp_z0, pitch0, last_x, cos_plus_zprimex, cos_minus_zprimex, theta_trk_x, theta_trk_x_str, true);
      Fill_corrected_dqdx_plots("plane1_trklen_60cm_passing_cathode_coszx", rr1, dqdx1, sp_x1, sp_z1, pitch1, last_x, cos_plus_zprimex, cos_minus_zprimex, theta_trk_x, theta_trk_x_str, true);
      if(fabs(cos_zx) < 0.75){
	Fill_track_plots("trklen_60cm_passing_cathode_coszx", dist_start, dist_end, rr2, dqdx2);
	FillHist("theta_trk_x_str_trklen_60cm_passing_cathode_coszx", theta_trk_x, 1., 100., 0., 100.);
        Fill_corrected_dqdx_plots("plane2_trklen_60cm_passing_cathode_coszx", rr2, dqdx2, sp_x2, sp_z2, pitch2, last_x, cos_plus_zprimex, cos_minus_zprimex, theta_trk_x, theta_trk_x_str, true);
      }
      
      // == Bellow are for plot approvals
      /*
      // == end dqdx median cut shift from 1600 to 1800
      if(end_meddqdx > 1800.){
	Fill_corrected_dqdx_plots("plane0_trklen_60cm_passing_cathode_coszx_meddqdx", rr0, dqdx0, sp_x0, sp_z0, pitch0, last_x, cos_plus_zprimex, cos_minus_zprimex, theta_trk_x_str, true);
        Fill_corrected_dqdx_plots("plane1_trklen_60cm_passing_cathode_coszx_meddqdx", rr1, dqdx1, sp_x1, sp_z1, pitch1, last_x, cos_plus_zprimex, cos_minus_zprimex, theta_trk_x_str, true);
        if(fabs(cos_zx) < 0.75) Fill_corrected_dqdx_plots("plane2_trklen_60cm_passing_cathode_coszx_meddqdx", rr2, dqdx2, sp_x2, sp_z2, pitch2, last_x, cos_plus_zprimex, cos_minus_zprimex, theta_trk_x_str, true);
      }

      // == End point away from cathode
      if(fabs(last_x) > 15. && fabs(last_x) < 185.){
	Fill_corrected_dqdx_plots("plane0_trklen_60cm_passing_cathode_coszx_cafv", rr0, dqdx0, sp_x0, sp_z0, pitch0, last_x, cos_plus_zprimex, cos_minus_zprimex, theta_trk_x_str, true);
	Fill_corrected_dqdx_plots("plane1_trklen_60cm_passing_cathode_coszx_cafv", rr1, dqdx1, sp_x1, sp_z1, pitch1, last_x, cos_plus_zprimex, cos_minus_zprimex, theta_trk_x_str, true);
	if(fabs(cos_zx) < 0.75) Fill_corrected_dqdx_plots("plane2_trklen_60cm_passing_cathode_coszx_cafv", rr2, dqdx2, sp_x2, sp_z2, pitch2, last_x, cos_plus_zprimex, cos_minus_zprimex, theta_trk_x_str, true);

	TString evt_id_str = Form("evt%d_id%d", *evt, *trkid);
	//cout << "evt_id_str : " << evt_id_str << ", length : " << this_reco_trk_len << " cm" << ", trklen : " << *trklen << ", ts : " << *ts << endl;
      }
      */   
    }
  }


  TString output_rootfile_dir = getenv("OUTPUTROOT_PATH");
  TString output_file_name = output_rootfile_dir + "/output_recom_loop_emb_run_" + run_number_str + ".root";
  //if(!isdata) output_file_name = output_rootfile_dir + "/output_recom_loop_emb_mc_2025a_spring.root";
  if(!isdata) output_file_name = output_rootfile_dir + "/output_recom_loop_emb_mc_2024b.root";
  out_rootfile = new TFile(output_file_name, "RECREATE");
  out_rootfile -> cd();
  
  hist_selected -> Write();
  WriteHist();

  out_rootfile -> Close();
}
