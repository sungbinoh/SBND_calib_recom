#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "mylib.h"

void run_calib_Ntuple() {

  // == Histograms for overal events
  TH1F *hist_selected = new TH1F("selected", "selected", 3., -0.5, 2.5);

  // Open the file containing the tree.
  TChain *fChain = new TChain("caloskim/TrackCaloSkim");
  TString input_file_dir = getenv("DATA_PATH");
  TString fileListPath = input_file_dir + "/sample_list/list_run_14608_local.txt";
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
  

  int N_entries = myReader.GetEntries();
  cout << "N_entries : " << N_entries << endl;
  int current_entry = 0;

  int N_run = 100;
  double ADC_med_cut = 1600.; // == https://sbn-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=23472&filename=SBND%20Calib%20Workshop%202021.pdf&version=1
  double track_length_cut = 15.;
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
    if(*run == 14608){
      if(*evt >= 9695) continue;
    }

    hist_selected -> Fill(*selected);
    //cout << "selected : " << *selected << endl;
    // == Tracks selected as Anode+Cathode crossing
    if (*selected == 2) {
      unsigned N_reco_hits = rr.GetSize();
      if(N_reco_hits < 1) continue;
      double this_reco_trk_len = rr[N_reco_hits - 1];
      if(rr[rr.GetSize() - 1] < track_length_cut) continue;

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

      //cout << "first_x : " << first_x << ", last_x : " << last_x << endl;

      FillHist("x_vs_y_throu", first_x, first_y, 1., 500., -250., 250., 5000., -250., 250.);
      FillHist("x_vs_y_throu", last_x, last_y, 1., 500., -250., 250., 5000., -250., 250.);
      
      FillHist("x_vs_z_throu", first_x, first_z, 1., 500., -250., 250., 600., -50., 550.);
      FillHist("x_vs_z_throu", last_x, last_z, 1., 500., -250., 250., 600., -50., 550.);

    }
  }

  TString output_rootfile_dir = getenv("OUTPUTROOT_PATH");
  out_rootfile = new TFile(output_rootfile_dir + "/output_test_14608.root", "RECREATE");
  out_rootfile -> cd();
  
  hist_selected -> Write();

  WriteHist();
  out_rootfile -> Close();
}
