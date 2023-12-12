#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

void AddFilesToChain(const std::string& fileListPath, TChain* chain) {
  std::ifstream fileStream(fileListPath);

  if (!fileStream.is_open()) {
    std::cerr << "Error opening file: " << fileListPath << std::endl;
    return;
  }

  std::string fileName;
  while (std::getline(fileStream, fileName)) {
    // Add each file to the TChain
    chain->Add(fileName.c_str());
  }

  fileStream.close();
}

void run_calib_Ntuple() {
  // Create a histogram for the values we read.
  TH2F *hist_dqdx_vs_rr = new TH2F("dqdx_rr","dqdx_rr", 300., 0., 300., 3000., 0., 3000.);

  // Open the file containing the tree.
  TChain *fChain = new TChain("caloskim/TrackCaloSkim");
  std::string fileListPath = "./list.txt";
  AddFilesToChain(fileListPath, fChain);

  TTreeReader myReader(fChain);

  // Variables we are going to read
  TTreeReaderValue<int> selected(myReader, "trk.selected");
  TTreeReaderArray<float> dqdx(myReader, "trk.hits2.dqdx"); // hits on plane 2 (Collection)
  TTreeReaderArray<float> rr(myReader, "trk.hits2.rr"); // hits on plane 2 (Collection)
  TTreeReaderArray<float> dir_x(myReader, "trk.dir.x"); 
  TTreeReaderArray<float> dir_y(myReader, "trk.dir.y");
  TTreeReaderArray<float> dir_z(myReader, "trk.dir.z");

  int N_entries = myReader.GetEntries();
  cout << "N_entries : " << N_entries << endl;
  int current_entry = 0;

  int N_run = 100;
  double ADC_med_cut = 1600.; // == https://sbn-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=23472&filename=SBND%20Calib%20Workshop%202021.pdf&version=1
  // Loop over all entries of the TTree
  while (myReader.Next()) {
    //if(current_entry > N_run) break;
   
    if(current_entry%100 == 0){
      cout << current_entry << " / " << N_entries << endl;
    }
    current_entry++;

    if (*selected == 0) { // Tracks selected as stopping

      // == Remove fake tracks
      if(rr[rr.GetSize() - 1] < 5.) continue; // remove tracks shorter than 5 cm

      vector<double> vec_ADC_for_med;
      for(unsigned i = 0; i < rr.GetSize(); i++) {
	double this_rr = rr[i];
	if(this_rr > 5.) break;
	vec_ADC_for_med.push_back(dqdx[i]);
      }
      double ADC_med = TMath::Median(vec_ADC_for_med.size(), &vec_ADC_for_med[0]);
      //cout << "ADC_med : " << ADC_med << endl;
      if(ADC_med < ADC_med_cut) continue;

      // == Fill histograms
      for (unsigned i = 0; i < dqdx.GetSize(); i++) {
	hist_dqdx_vs_rr->Fill(rr[i], dqdx[i]);
      }
    }
  }

  TFile *out_rootfile = new TFile("./output.root", "RECREATE");
  out_rootfile -> cd();
  hist_dqdx_vs_rr -> Write();
  
  out_rootfile -> Close();
}
