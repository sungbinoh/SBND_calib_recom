#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "mylib.h"
#include "Math/Vector3D.h"
#include "BetheBloch.h"

TSpline3 * muon_sp_range_to_KE = Get_sp_range_KE(mass_muon);

void Fill_track_plots(TString suffix, double dist_start, double dist_end, const TTreeReaderArray<float>& rr, const TTreeReaderArray<float>& dqdx){

  FillHist("dist_start_" + suffix, dist_start, 1., 1000., 0., 1000.);
  FillHist("dist_end_" + suffix, dist_end, 1., 1000., 0., 1000.);
  for (unsigned i = 0; i < dqdx.GetSize(); i++) {
    FillHist("rr_vs_dqdx_" + suffix, rr[i], dqdx[i], 1., 300., 0., 300., 3000., 0., 3000.);
  }
}

void Fill_hit_plots(TString suffix, const TTreeReaderArray<float>& rr, const TTreeReaderArray<float>& dqdx, double cos_xy, double cos_yz, double cos_zx){

  for (unsigned i = 0; i < dqdx.GetSize(); i++) {
    FillHist("cos_xy_vs_dqdx_" + suffix, cos_xy, dqdx[i], 1., 220., -1.1, 1.1, 3000., 0., 3000.);
    FillHist("cos_yz_vs_dqdx_" + suffix, cos_yz, dqdx[i], 1., 220., -1.1, 1.1, 3000., 0., 3000.);
    FillHist("cos_zx_vs_dqdx_" + suffix, cos_zx, dqdx[i], 1., 220., -1.1, 1.1, 3000., 0., 3000.);
  }
}

void Fill_corrected_dqdx_plots(TString suffix, const TTreeReaderArray<float>& rr, const TTreeReaderArray<float>& dqdx, const TTreeReaderArray<float>& sp_x, const TTreeReaderArray<float>& pitch){

  for (unsigned i = 0; i < dqdx.GetSize(); i++) {
    double this_lifetime_corr = Lifetime_Correction(sp_x[i], 10.0);
    double corrected_dqdx = dqdx[i] * this_lifetime_corr;
    FillHist("rr_vs_corr_dqdx_" + suffix, rr[i], corrected_dqdx, 1., 300., 0., 300., 3000., 0., 3000.);
    FillHist("rr_vs_pitch_" + suffix, rr[i], pitch[i], 1., 300., 0., 300., 200., 0., 2.);

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
void run_recom_loop() {

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
  TString fileListPath = input_file_dir + "/sample_list/list_run_14480_local.txt";
  AddFilesToChain(fileListPath, fChain);

  TTreeReader myReader(fChain);

  // == Variables
  TTreeReaderValue<int> run(myReader, "trk.meta.run");
  TTreeReaderValue<int> evt(myReader, "trk.meta.evt");

  TTreeReaderValue<int> selected(myReader, "trk.selected");
  TTreeReaderArray<float> dqdx(myReader, "trk.hits2.dqdx"); // hits on plane 2 (Collection)
  TTreeReaderArray<float> rr(myReader, "trk.hits2.rr");
  TTreeReaderArray<float> pitch(myReader, "trk.hits2.pitch");
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

  int N_run = 300;
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
      if(*evt >= 748) continue;
    }
    if(*run == 14608){
      if(*evt >= 9695) continue;
    }
    
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

      // == Pick first and last hits to check if the track is passing the cathode
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
      Fill_track_plots("trklen_60cm_passing_cathode", dist_start, dist_end, rr, dqdx);
      Fill_hit_plots("trklen_60cm_passing_cathode", rr, dqdx, cos_xy, cos_yz, cos_zx);

      if(fabs(cos_zx) > 0.75) continue;
      Fill_track_plots("trklen_60cm_passing_cathode_coszx", dist_start, dist_end, rr, dqdx);
      Fill_corrected_dqdx_plots("trklen_60cm_passing_cathode_coszx", rr, dqdx, sp_x, pitch);

    }
  }

  TString output_rootfile_dir = getenv("OUTPUTROOT_PATH");
  out_rootfile = new TFile(output_rootfile_dir + "/output_recom_run14480.root", "RECREATE");
  out_rootfile -> cd();
  
  hist_selected -> Write();
  WriteHist();

  out_rootfile -> Close();
}
