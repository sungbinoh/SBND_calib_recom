#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "mylib.h"

void run_lifetime_loop(int run_num = 0) {

  bool isdata = false;
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

  // == Histograms for stopping muons
  TH1F *hist_true_end_process = new TH1F("true_end_process", "true_end_process", 65., -0.5, 64.5);
  TH1F *hist_true_end_process_track_len_cut = new TH1F("true_end_process_track_len_cut", "true_end_process_track_len_cut", 65., -0.5, 64.5);

  // == Histograms for cathod-anode passing muons
  TH2F *hist_time_vs_dqdx = new TH2F("time_dqdx","time_dqdx", 4000., 0., 4000., 3000., 0., 3000.);
  TH2F *hist_sp_x_vs_dqdx = new TH2F("sp_x_dqdx","sp_x_dqdx", 500., -250., 250., 3000., 0., 3000.);
  TH2F *hist_sp_x_vs_corr_dqdx = new TH2F("sp_x_corr_dqdx","sp_x_corr_dqdx", 500., -250., 250., 3000., 0., 3000.);
  TH2F *hist_tdrift_vs_dqdx = new TH2F("tdrift_vs_dqdx", "tdrift_vs_dqdx", 150., 0., 1.5, 3000., 0., 3000.);
  TH2F *hist_tdrift_vs_corr_dqdx = new TH2F("tdrift_vs_corr_dqdx", "tdrift_vs_corr_dqdx", 150., 0., 1.5, 3000., 0., 3000.);

  /////////////////////////////////
  // == Call Trees
  /////////////////////////////////
  // Open the file containing the tree
  TChain *fChain = new TChain("caloskim/TrackCaloSkim");
  TString input_file_dir = getenv("DATA_PATH");
  TString fileListPath = input_file_dir + "/sample_list/list_run_" + run_str + "_local.txt";
  AddFilesToChain(fileListPath, fChain);

  TTreeReader myReader(fChain);

  // == Variables
  TTreeReaderValue<int> run(myReader, "trk.meta.run");
  TTreeReaderValue<int> evt(myReader, "trk.meta.evt");

  TTreeReaderValue<int> selected(myReader, "trk.selected");
  TTreeReaderArray<float> dqdx(myReader, "trk.hits2.dqdx"); // hits on plane 2 (Collection)
  TTreeReaderArray<float> rr(myReader, "trk.hits2.rr");
  TTreeReaderArray<float> time(myReader, "trk.hits2.h.time"); 
  TTreeReaderArray<float> sp_x(myReader, "trk.hits2.h.sp.x");
  TTreeReaderArray<float> sp_y(myReader, "trk.hits2.h.sp.y");
  TTreeReaderArray<float> sp_z(myReader, "trk.hits2.h.sp.z");
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

  int N_run = 100;
  double ADC_med_cut = 1600.; // == https://sbn-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=23472&filename=SBND%20Calib%20Workshop%202021.pdf&version=1
  double track_length_cut = 60.;
  // Loop over all entries of the TTree
  while (myReader.Next()) {
    //if(current_entry > N_run) break;
   
    if(current_entry%100 == 0){
      cout << current_entry << " / " << N_entries << endl;
    }
    current_entry++;

    if(*run == 14480){
      if(*evt >= 749) continue;
    }

    hist_selected -> Fill(*selected);

    // == Tracks selected as Anode+Cathode crossing
    if (*selected == 1) {

      unsigned N_reco_hits = rr.GetSize();
      double this_reco_trk_len = rr[N_reco_hits - 1];

      double first_x = -999.;
      double first_y = -999.;
      double first_z = -999.;
      double last_x = -999.;
      double last_y = -999.;
      double last_z = -999.;

      first_x = sp_x[N_reco_hits - 1];
      first_y = sp_y[N_reco_hits - 1];
      first_z = sp_z[N_reco_hits - 1];

      for (unsigned i = 0; i < dqdx.GetSize(); i++) {
        if(rr[i] > 0){
          last_x = sp_x[i];
          last_y = sp_y[i];
          last_z = sp_z[i];
          break;
	}
      }
      bool passing_cathode = false;
      ROOT::Math::XYZVector track_vec(last_x - first_x, last_y - first_y, last_z - first_z);
      double cos_xy = track_vec.X() / (sqrt(pow(track_vec.X(), 2.) + pow(track_vec.Y(), 2.)));
      double cos_yz = track_vec.Y() / (sqrt(pow(track_vec.Y(), 2.) + pow(track_vec.Z(), 2.)));
      double cos_zx = track_vec.X() / (sqrt(pow(track_vec.X(), 2.) + pow(track_vec.Z(), 2.)));
      if(first_x * last_x < 0.) passing_cathode= true;

      //cout << "passing_cathode : " << passing_cathode << ", this_reco_trk_len : " <<this_reco_trk_len << ", first_x : " << first_x << ", last_x : " << last_x << endl;
      if(rr[rr.GetSize() - 1] < track_length_cut) continue; // remove tracks shorter than 60 cm
      //cout << "[(*selected == 1)] dqdx.GetSize() : " << dqdx.GetSize() << ", true_hit_time.GetSize() : " << true_hit_time.GetSize() << endl;
      
      // == Collect maximum hit pick time of the track
      double t_min = 999.;
      for (unsigned i = 0; i < dqdx.GetSize(); i++) {
	double this_t = time[i] * 5.0e-4;
	if(this_t < t_min && rr[i] > 0) t_min = this_t;
	//cout << i << ", this_t :  " << this_t << ", sp_x[i] : " << sp_x[i] << ", rr[i] : " << rr[i] << endl;
      }

      // == No angular cut
      for (unsigned i = 0; i < dqdx.GetSize(); i++) {
	if(rr[i] < 0.) continue;
	
	FillHist("x_vs_y_acx", first_x, first_y, 1., 5000., -250., 250., 5000., -250., 250.);
	FillHist("x_vs_y_acx", last_x, last_y, 1., 5000., -250., 250., 5000., -250., 250.);
	
	FillHist("x_vs_z_acx", first_x, first_z, 1., 5000., -250., 250., 5000., -50., 550.);
	FillHist("x_vs_z_acx", last_x, last_z, 1., 5000., -250., 250., 5000., -50., 550.);
	
	FillHist("x_vs_y_trk_len", sp_x[i], sp_y[i], 1., 500., -250., 250., 500., -250., 250.);
        FillHist("x_vs_z_trk_len", sp_x[i], sp_z[i], 1., 500., -250., 250., 600., -50., 550.);

	//cout << i << ", time : " << time[i] << ", sp_x : " << sp_x[i] << endl;
	double t_drift_time = time[i] * 5.0e-4 - t_min;
        double t_drift_space = (200. - fabs(sp_x[i])) / v_drift;
	FillHist("sp_x_vs_dqdx_trk_len", sp_x[i], dqdx[i], 1., 500., -250., 250., 3000., 0., 3000.);
	FillHist("t_t0_vs_dqdx_trk_len", time[i], dqdx[i], 1., 3000., 0., 3000., 3000., 0., 3000.);
	FillHist("t_drift_time_vs_dqdx_trk_len", t_drift_time, dqdx[i], 1., 150., 0., 1.5, 3000., 0., 3000.);
        FillHist("t_drift_space_vs_dqdx_trk_len", t_drift_space, dqdx[i], 1., 150., 0., 1.5, 3000., 0., 3000.);

	double this_lifetime_corr = Lifetime_Correction(sp_x[i], 12.3);
	double corrected_dqdx = dqdx[i] * this_lifetime_corr;
	//cout << i << ", x : " << sp_x[i] << ", this_lifetime_corr : " << this_lifetime_corr << ", dqdx[i] : " << dqdx[i] << ", corrected_dqdx : " << corrected_dqdx << endl;
	FillHist("sp_x_vs_corr_dqdx_trk_len", sp_x[i], corrected_dqdx, 1., 500., -250., 250., 3000., 0., 3000.);
	FillHist("t_drift_time_vs_dqdx_corr_trk_len", t_drift_time, corrected_dqdx, 1., 150., 0., 1.5, 3000., 0., 3000.);
	FillHist("t_drift_space_vs_dqdx_corr_trk_len", t_drift_space, corrected_dqdx, 1.,150., 0., 1.5, 3000., 0., 3000.);

	FillHist("cos_xy_vs_dqdx_trk_len", cos_xy, dqdx[i], 1., 220., -1.1, 1.1, 3000., 0., 3000.);
	FillHist("cos_yz_vs_dqdx_trk_len", cos_yz, dqdx[i], 1., 220., -1.1, 1.1, 3000., 0., 3000.);
	FillHist("cos_zx_vs_dqdx_trk_len", cos_zx, dqdx[i], 1., 220., -1.1, 1.1, 3000., 0., 3000.);
      }

      // == With angular and passing cathode cut
      if(fabs(cos_zx) > 0.75 || !passing_cathode) continue;

      for(unsigned i = 0; i < dqdx.GetSize(); i++) {
	if(rr[i] < 0.) continue;

	FillHist("x_vs_y_angle_passing_cathode", sp_x[i], sp_y[i], 1., 500., -250., 250., 500., -250., 250.);
	FillHist("x_vs_z_angle_passing_cathode", sp_x[i], sp_z[i], 1., 500., -250., 250., 600., -50., 550.);

	double this_lifetime_corr = Lifetime_Correction(sp_x[i], 12.3);
        double corrected_dqdx = dqdx[i] * this_lifetime_corr;

	double t_drift_time = time[i] * 5.0e-4 - t_min;
	double t_drift_space = (200. - fabs(sp_x[i])) / v_drift;
	FillHist("t_drift_time_vs_dqdx_angle_passing_cathode", t_drift_time, dqdx[i], 1., 150., 0., 1.5, 3000., 0., 3000.);
	FillHist("t_drift_space_vs_dqdx_angle_passing_cathode", t_drift_space, dqdx[i], 1., 150., 0., 1.5, 3000., 0., 3000.);
	FillHist("t_drift_time_vs_dqdx_corr_angle_passing_cathode", t_drift_time, corrected_dqdx, 1., 150., 0., 1.5, 3000., 0., 3000.);
        FillHist("t_drift_space_vs_dqdx_corr_angle_passing_cathode", t_drift_space, corrected_dqdx, 1., 150., 0., 1.5, 3000., 0., 3000.);
      }
    }
  }

  TString output_rootfile_dir = getenv("OUTPUTROOT_PATH");
  out_rootfile = new TFile(output_rootfile_dir + "/output_lifetime_" + run_str + ".root", "RECREATE");
  out_rootfile -> cd();
  
  hist_selected -> Write();

  WriteHist();
  out_rootfile -> Close();
}
