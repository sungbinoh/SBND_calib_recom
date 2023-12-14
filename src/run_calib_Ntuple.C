#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

void run_calib_Ntuple() {

  // == Histograms for overal events
  TH1F *hist_selected = new TH1F("selected", "selected", 3., -0.5, 2.5);

  // == Histograms for stopping muons
  TH2F *hist_dqdx_vs_rr = new TH2F("dqdx_rr","dqdx_rr", 300., 0., 300., 3000., 0., 3000.);
  TH2F *hist_dqdx_vs_rr_ADC_med_cut = new TH2F("dqdx_rr_ADC_med_cut","dqdx_rr_ADC_med_cut", 300., 0., 300., 3000., 0., 3000.);
  TH2F *hist_dqdx_vs_rr_track_len_cut = new TH2F("dqdx_rr_track_len_cut","dqdx_rr_track_len_cut", 300., 0., 300., 3000., 0., 3000.);
  TH1F *hist_true_end_process = new TH1F("true_end_process", "true_end_process", 65., -0.5, 64.5);
  TH1F *hist_true_end_process_track_len_cut = new TH1F("true_end_process_track_len_cut", "true_end_process_track_len_cut", 65., -0.5, 64.5);

  // == Histograms for cathod-anode passing muons
  TH2F *hist_time_vs_dqdx = new TH2F("time_dqdx","time_dqdx", 300., 0., 300., 3000., 0., 3000.);

  // Open the file containing the tree.
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
  

  int N_entries = myReader.GetEntries();
  cout << "N_entries : " << N_entries << endl;
  int current_entry = 0;

  int N_run = 100;
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
      // == With no cut
      hist_true_end_process -> Fill(true_end_process[0]);
      for (unsigned i = 0; i < dqdx.GetSize(); i++) {
        hist_dqdx_vs_rr -> Fill(rr[i], dqdx[i]);
      }

      // == With ADC Median Cut
      vector<double> vec_ADC_for_med;
      for(unsigned i = 0; i < rr.GetSize(); i++) {
        double this_rr = rr[i];
        if(this_rr > 5.) break;
        vec_ADC_for_med.push_back(dqdx[i]);
      }
      double ADC_med = TMath::Median(vec_ADC_for_med.size(), &vec_ADC_for_med[0]);
      //if(ADC_med < ADC_med_cut) continue;
      for (unsigned i = 0; i < dqdx.GetSize(); i++) {
        hist_dqdx_vs_rr_ADC_med_cut -> Fill(rr[i], dqdx[i]);
      }

      // == With Track Length cut
      if(rr[rr.GetSize() - 1] < track_length_cut) continue; // remove tracks shorter than 15 cm
      hist_true_end_process_track_len_cut -> Fill(true_end_process[0]);
      for (unsigned i = 0; i < dqdx.GetSize(); i++) {
	hist_dqdx_vs_rr_track_len_cut -> Fill(rr[i], dqdx[i]);
      }
    }

    // == Tracks selected as Anode+Cathode crossing
    if (*selected == 1) {
      //continue;
      cout << "[(*selected == 1)] dqdx.GetSize() : " << dqdx.GetSize() << ", true_hit_time.GetSize() : " << true_hit_time.GetSize() << endl;
      for (unsigned i = 0; i < true_hit_time.GetSize(); i++) {
        cout << i << ", true_hit_time : " << true_hit_time[i] << ", true_hit_tdrift : " << true_hit_tdrift[i] << endl;
      }
      for (unsigned i = 0; i < dqdx.GetSize(); i++) {
	cout << i << ", time : " << time[i] << ", sp_x : " << sp_x[i] << endl;
	hist_time_vs_dqdx -> Fill(time[i], dqdx[i]);
      }
    }
  }

  TString output_rootfile_dir = getenv("OUTPUTROOT_PATH");
  TFile *out_rootfile = new TFile(output_rootfile_dir + "/output.root", "RECREATE");
  out_rootfile -> cd();
  
  hist_selected -> Write();

  hist_true_end_process -> Write();
  hist_dqdx_vs_rr -> Write();
  hist_dqdx_vs_rr_ADC_med_cut -> Write();
  hist_true_end_process_track_len_cut -> Write();
  hist_dqdx_vs_rr_track_len_cut -> Write();

  hist_time_vs_dqdx -> Write();

  out_rootfile -> Close();
}
