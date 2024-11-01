#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "mylib.h"
#include "TrackCaloSkimmerObj.h"

double Modbox_dedx(double dqdx, double this_c_cal){

  double alpha = 0.93;
  double beta = 0.212;

  double dedx = ( exp( (beta * dqdx)/(29449.153 * this_c_cal) ) - alpha) / (1.4388489 * beta);

  return dedx;
}

void run_waveform() {

  //gSystem->Load("./lib/TrackCaloSkimmerObj_h.so");
  
  // == Histograms for overal events
  TH1F *hist_selected = new TH1F("selected", "selected", 3., -0.5, 2.5);

  // Open the file containing the tree.
  TChain *fChain = new TChain("caloskim/TrackCaloSkim");
  TString input_file_dir = getenv("DATA_PATH");
  TString fileListPath = input_file_dir + "/sample_list/list_run_14860_waveform.txt";
  AddFilesToChain(fileListPath, fChain);

  TTreeReader myReader(fChain);

  // == Variables
  TTreeReaderValue<int> run(myReader, "trk.meta.run");
  TTreeReaderValue<int> evt(myReader, "trk.meta.evt");
  TTreeReaderValue<unsigned long> ts(myReader, "trk.meta.time");
  TTreeReaderValue<int> trkid(myReader, "trk.id");
  TTreeReaderValue<float> trklen(myReader, "trk.length");
  
  TTreeReaderValue<int> selected(myReader, "trk.selected");
  TTreeReaderArray<float> dqdx0(myReader, "trk.hits0.dqdx");
  TTreeReaderArray<float> dqdx1(myReader, "trk.hits1.dqdx");
  TTreeReaderArray<float> dqdx2(myReader, "trk.hits2.dqdx"); // hits on plane 2 (Collection)
  TTreeReaderArray<float> rr0(myReader, "trk.hits0.rr");
  TTreeReaderArray<float> rr1(myReader, "trk.hits1.rr");
  TTreeReaderArray<float> rr2(myReader, "trk.hits2.rr");
  TTreeReaderArray<float> time(myReader, "trk.hits2.h.time"); 
  TTreeReaderArray<float> sp_x(myReader, "trk.hits2.h.sp.x");
  TTreeReaderArray<float> sp_y(myReader, "trk.hits2.h.sp.y");
  TTreeReaderArray<float> sp_z(myReader, "trk.hits2.h.sp.z");
  TTreeReaderArray<float> dir_x(myReader, "trk.dir.x"); 
  TTreeReaderArray<float> dir_y(myReader, "trk.dir.y");
  TTreeReaderArray<float> dir_z(myReader, "trk.dir.z");

  TTreeReaderArray<sbn::WireInfo> wires0(myReader, "trk.wires0");
  TTreeReaderArray<sbn::WireInfo> wires1(myReader, "trk.wires1");
  TTreeReaderArray<sbn::WireInfo> wires2(myReader, "trk.wires2");
   
  int N_entries = myReader.GetEntries();
  cout << "N_entries : " << N_entries << endl;
  int current_entry = 0;

  int N_run = 1000;
  double ADC_med_cut = 1600.; // == https://sbn-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=23472&filename=SBND%20Calib%20Workshop%202021.pdf&version=1
  double track_length_cut = 15.;
  // Loop over all entries of the TTree
  while (myReader.Next()) {
    //if(current_entry > N_run) break;

    if(current_entry%100 == 0){
      //cout << current_entry << " / " << N_entries << endl;
    }
    current_entry++;

    TString evt_id_str = Form("evt%d_id%d", *evt, *trkid);
    //cout << "evt_id_str : " << evt_id_str << endl;
    //if(*evt == 7160 && *trkid==0){
    
    if(*evt == 8473 && *trkid==2){
    //if(*evt == 10130 && *trkid==0){

      double this_trk_len0 = rr0[rr0.GetSize() - 1];
      double this_trk_len1 = rr1[rr1.GetSize() - 1];
      double this_trk_len2 = rr2[rr2.GetSize() - 1];
      cout << evt_id_str << " " << Form("trk length : %f, %f, %f, %f, timestamp : %lu", this_trk_len0, this_trk_len1, this_trk_len2, *trklen, *ts) << endl;
      
      if(wires0.GetSize() > 0){
	for(int i = 0; i < wires0.GetSize(); i++){
	  TString this_id = evt_id_str +  Form("_plane%d_wire%d", wires0[i].plane, wires0[i].wire);
	  int N_adcs = wires0[i].adcs.size();
	  maphist_TH1D[this_id] = new TH1D(this_id, this_id, N_adcs, 0., N_adcs + 0.);
	  for(int j = 0; j < N_adcs; j++){
	    maphist_TH1D[this_id] -> SetBinContent(j + 1, wires0[i].adcs.at(j));
	  }
	}
      }

      if(wires1.GetSize() > 0){
	for(int i = 0; i < wires1.GetSize(); i++){
	  TString this_id = evt_id_str + Form("_plane%d_wire%d", wires1[i].plane, wires1[i].wire);
	  int N_adcs = wires1[i].adcs.size();
	  maphist_TH1D[this_id] = new TH1D(this_id, this_id, N_adcs, 0., N_adcs + 0.);
	  for(int j = 0; j < N_adcs; j++){
	    maphist_TH1D[this_id] -> SetBinContent(j + 1, wires1[i].adcs.at(j));
	  }
	}
      }

      if(wires2.GetSize() > 0){
	for(int i = 0; i < wires2.GetSize(); i++){
	  TString this_id = evt_id_str + Form("_plane%d_wire%d", wires2[i].plane, wires2[i].wire);
	  int N_adcs = wires2[i].adcs.size();
	  maphist_TH1D[this_id] = new TH1D(this_id, this_id, N_adcs, 0., N_adcs + 0.);
	  for(int j = 0; j < N_adcs; j++){
	    maphist_TH1D[this_id] -> SetBinContent(j + 1, wires2[i].adcs.at(j));
	  }
	}
      }

      // == Save dE/dx vs rr
      vector<double> rr0_vec, rr1_vec, rr2_vec;
      vector<double> dqdx0_vec, dqdx1_vec, dqdx2_vec;
      vector<double> dedx0_vec, dedx1_vec, dedx2_vec;
      for(int i = 0; i < rr0.GetSize(); i++){
	rr0_vec.push_back(rr0[i]);
	dqdx0_vec.push_back(dqdx0[i]);
	double this_dedx = Modbox_dedx(dqdx0[i], 0.02038);
	dedx0_vec.push_back(this_dedx);
      }
      for(int i = 0; i < rr1.GetSize(); i++){
        rr1_vec.push_back(rr1[i]);
	dqdx1_vec.push_back(dqdx1[i]);
	double this_dedx = Modbox_dedx(dqdx1[i], 0.02101);
	dedx1_vec.push_back(this_dedx);
      }
      for(int i = 0; i < rr2.GetSize(); i++){
        rr2_vec.push_back(rr2[i]);
	dqdx2_vec.push_back(dqdx2[i]);
	double this_dedx = Modbox_dedx(dqdx2[i], 0.01928);
	dedx2_vec.push_back(this_dedx);
      }

      map_gr[evt_id_str + "_plane0_dedx"] = new TGraph(rr0_vec.size(), &rr0_vec[0], &dedx0_vec[0]);
      map_gr[evt_id_str + "_plane1_dedx"] = new TGraph(rr1_vec.size(), &rr1_vec[0], &dedx1_vec[0]);
      map_gr[evt_id_str + "_plane2_dedx"] = new TGraph(rr2_vec.size(), &rr2_vec[0], &dedx2_vec[0]);
      map_gr[evt_id_str + "_plane0_dedx"] -> SetName(evt_id_str + "_plane0_dedx");
      map_gr[evt_id_str + "_plane1_dedx"] -> SetName(evt_id_str + "_plane1_dedx");
      map_gr[evt_id_str + "_plane2_dedx"] -> SetName(evt_id_str + "_plane2_dedx");
      
      map_gr[evt_id_str + "_plane0_dqdx"] = new TGraph(rr0_vec.size(), &rr0_vec[0], &dqdx0_vec[0]);
      map_gr[evt_id_str + "_plane1_dqdx"] = new TGraph(rr1_vec.size(), &rr1_vec[0], &dqdx1_vec[0]);
      map_gr[evt_id_str + "_plane2_dqdx"] = new TGraph(rr2_vec.size(), &rr2_vec[0], &dqdx2_vec[0]);
      map_gr[evt_id_str + "_plane0_dqdx"] -> SetName(evt_id_str + "_plane0_dqdx");
      map_gr[evt_id_str + "_plane1_dqdx"] -> SetName(evt_id_str + "_plane1_dqdx");
      map_gr[evt_id_str + "_plane2_dqdx"] -> SetName(evt_id_str + "_plane2_dqdx");
    }
    
    hist_selected -> Fill(*selected);
    //cout << "selected : " << *selected << endl;
    // == Tracks selected as Anode+Cathode crossing

    /*
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
    }
    */
  }

  cout << "Writing output file" << endl;
  TString output_rootfile_dir = getenv("OUTPUTROOT_PATH");
  out_rootfile = new TFile(output_rootfile_dir + "/output_waveform_14860.root", "RECREATE");
  out_rootfile -> cd();
  
  hist_selected -> Write();

  WriteHist();

  for(std::map< TString, TGraph* >::iterator mapit = map_gr.begin(); mapit!=map_gr.end(); mapit++){
    mapit->second->Write();
  }
  out_rootfile -> Close();
}
