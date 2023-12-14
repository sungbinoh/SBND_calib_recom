#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "mylib.h"

void run_lifetime_loop() {

  /////////////////////////////////
  // == Constants
  /////////////////////////////////
  double v_drift = 156.267; // -- [cm/ms]


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
  TH2F *hist_tdrift_vs_dqdx = new TH2F("tdrift_vs_dqdx", "tdrift_vs_dqdx", 150., 0., 1.5, 3000., 0., 3000.);


  /////////////////////////////////
  // == Call Trees
  /////////////////////////////////
  // Open the file containing the tree
  TChain *fChain = new TChain("caloskim/TrackCaloSkim");
  TString input_file_dir = getenv("DATA_PATH");
  TString fileListPath = input_file_dir + "/sample_list/list_MCP2023B_corsika_1Dsim_1Dreco.txt";
  AddFilesToChain(fileListPath, fChain);

  TTreeReader myReader(fChain);

  // == Variables
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

  int N_run = 3000;
  double ADC_med_cut = 1600.; // == https://sbn-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=23472&filename=SBND%20Calib%20Workshop%202021.pdf&version=1
  double track_length_cut = 15.;
  // Loop over all entries of the TTree
  while (myReader.Next()) {
    if(current_entry > N_run) break;
   
    if(current_entry%100 == 0){
      cout << current_entry << " / " << N_entries << endl;
    }
    current_entry++;

    hist_selected -> Fill(*selected);

    // == Tracks selected as Anode+Cathode crossing
    if (*selected == 1) {
      //cout << "[(*selected == 1)] dqdx.GetSize() : " << dqdx.GetSize() << ", true_hit_time.GetSize() : " << true_hit_time.GetSize() << endl;
      for (unsigned i = 0; i < dqdx.GetSize(); i++) {
	//cout << i << ", time : " << time[i] << ", sp_x : " << sp_x[i] << endl;
	double this_drift_time = (200. - fabs(sp_x[i])) / v_drift;
	hist_sp_x_vs_dqdx -> Fill(sp_x[i], dqdx[i]);
	hist_time_vs_dqdx -> Fill(time[i], dqdx[i]);
	hist_tdrift_vs_dqdx -> Fill(this_drift_time, dqdx[i]);
      }
    }
  }

  TString output_rootfile_dir = getenv("OUTPUTROOT_PATH");
  TFile *out_rootfile = new TFile(output_rootfile_dir + "/output_lifetime.root", "RECREATE");
  out_rootfile -> cd();
  
  hist_selected -> Write();
  hist_true_end_process -> Write();
  hist_true_end_process_track_len_cut -> Write();
  hist_sp_x_vs_dqdx -> Write();
  hist_time_vs_dqdx -> Write();
  hist_tdrift_vs_dqdx -> Write();

  out_rootfile -> Close();
}
