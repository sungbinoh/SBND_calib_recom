#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "mylib.h"


void Fill_track_plots(TString suffix, double dist_start, double dist_end){

  FillHist("dist_start_" + suffix, dist_start, 1., 1000., 0., 1000.);
  FillHist("dist_end_" + suffix, dist_end, 1., 1000., 0., 1000.);
}

void run_recom_loop() {

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

  int N_run = 1000;
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

    // == Tracks selected as stopping
    if (*selected == 0) {

      unsigned N_reco_hits = rr.GetSize();
      TVector3 this_reco_start(sp_x[N_reco_hits - 1], sp_y[N_reco_hits - 1], sp_z[N_reco_hits - 1]);
      TVector3 this_reco_end(sp_x[0], sp_y[0], sp_z[0]);
      TVector3 this_true_start(true_start_x[0], true_start_y[0], true_start_z[0]);
      TVector3 this_true_end(true_end_x[0], true_end_y[0], true_end_z[0]);

      double dist_start = (this_reco_start - this_true_start).Mag();
      double dist_end = (this_reco_end - this_true_end).Mag();

      double this_reco_trk_len = rr[N_reco_hits - 1];

      // == Nocut
      FillHist("dist_start_vs_reco_trk_len_nocut", dist_start, this_reco_trk_len, 1., 1000., 0., 1000., 1000., 0., 1000.);
      FillHist("dist_end_vs_reco_trk_len_nocut", dist_end, this_reco_trk_len, 1., 1000., 0., 1000., 1000., 0., 1000.);
      Fill_track_plots("nocut", dist_start, dist_end);


      // == Track length 5 cm cut
      if(this_reco_trk_len < 5.) continue;
      Fill_track_plots("trklen_05cm", dist_start, dist_end);

      // == Track length 10 cm cut
      if(this_reco_trk_len < 10.) continue;
      Fill_track_plots("trklen10cm", dist_start, dist_end);

      // == Track length 15 cm cut
      if(this_reco_trk_len < 15.) continue;
      Fill_track_plots("trklen_15cm", dist_start, dist_end);

      // == Track length 20 cm cut
      if(this_reco_trk_len < 20.) continue;
      Fill_track_plots("trklen_20cm", dist_start, dist_end);

      // == Track length 30 cm cut
      if(this_reco_trk_len < 30.) continue;
      Fill_track_plots("trklen_30cm", dist_start, dist_end);
      
    }
  }

  TString output_rootfile_dir = getenv("OUTPUTROOT_PATH");
  out_rootfile = new TFile(output_rootfile_dir + "/output_recom.root", "RECREATE");
  out_rootfile -> cd();
  
  hist_selected -> Write();
  hist_true_end_process -> Write();
  hist_true_end_process_track_len_cut -> Write();
  hist_sp_x_vs_dqdx -> Write();
  hist_time_vs_dqdx -> Write();
  hist_tdrift_vs_dqdx -> Write();
  hist_tdrift_vs_corr_dqdx -> Write();

  WriteHist();

  out_rootfile -> Close();
}
