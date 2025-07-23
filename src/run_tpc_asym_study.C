#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "mylib.h"
#include "Math/Vector3D.h"
#include "BetheBloch.h"
#include "SCECorr.h"

BetheBloch *muon_BB = new BetheBloch(13);

SCECorr *sce_corr_mc = new SCECorr(false);

bool isdata = false;

// == binnings for traj plots
const double x_min = -200, x_max = 200., y_min = -200, y_max = 200, z_min = 0, z_max = 500;
const double pixel_size = 0.1;
const int x_bins = (x_max - x_min) / pixel_size;
const int y_bins = (y_max - y_min) / pixel_size;
const int z_bins = (z_max - z_min) / pixel_size;

/*
TString get_y_str(double y){

  TString out = "";
  if(y < -200.) out = "lm200";
  else if(y < -150.) out = "m200tom150";
  else if(y < -100.) out = "m150tom100";
  else if(y < -50.) out = "m100tom50";
  else if(y < 0.) out = "m50to0";
  else if(y < 50.) out = "0top50";
  else if(y < 100.) out = "p50top100";
  else if(y < 150.) out = "p100top150";
  else if(y < 200.) out = "p150top200";
  else out = "gp200";

  return out;
}
*/
TString get_y_str(double y){

  TString out = "";
  if(y < -200.) out = "lm200";
  else if(y > 200.) out = "gp200";
  else if(y > -180. && y < -120.) out = "m180tom120";
  else if(y > -80. && y < -20.) out = "m80tom20";
  else if(y > 20. && y < 80.) out = "p20top80";
  else if(y > 120. && y < 200.) out = "p120top200";
  else out = "support";

  return out;
}


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

bool evt_sel(const TTreeReaderArray<float> &sp_x, const TTreeReaderArray<float> &sp_y, const TTreeReaderArray<float> &sp_z, const TTreeReaderArray<float> &rr, const TTreeReaderArray<float> &dqdx){

  bool out = false;

  unsigned N_reco_hits = rr.GetSize();
  if(N_reco_hits < 3) return out;

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
  double zprime_plus = zprime_60deg(track_vec.Y(), track_vec.Z(), 1);
  double zprime_minus = zprime_60deg(track_vec.Y(), track_vec.Z(), -1);
  double cos_plus_zprimex = track_vec.X() / (sqrt(pow(track_vec.X(), 2.) + pow(zprime_plus, 2.)));
  double cos_minus_zprimex = track_vec.X() / (sqrt(pow(track_vec.X(), 2.) + pow(zprime_minus, 2.)));

  if(first_x * last_x < 0.) passing_cathode = true;
  if(!passing_cathode) return out;

  // == aca
  bool is_ac = fabs(first_x) > 195. || fabs(last_x) > 195.;
  if(!is_ac) return out;
  
  out = true;
  return out;
}

void make_traj_plot(int current_entry, int plane, const TTreeReaderArray<float> &sp_x, const TTreeReaderArray<float> &sp_y, const TTreeReaderArray<float> &sp_z,
		    const TTreeReaderArray<float> &dir_x, const TTreeReaderArray<float> &dir_y, const TTreeReaderArray<float> &dir_z,
		    const TTreeReaderArray<float> &rr, const TTreeReaderArray<float> &dqdx, const TTreeReaderArray<float> &q_integ,
		    const TTreeReaderArray<float> &time
		    ){

  TString plane_str = Form("%d", plane);
  unsigned N_reco_hits = sp_x.GetSize();
    
  int idx_largest_x = -1;
  int idx_smallest_x = -1;
  double largest_x = -999.;
  double smallest_x = 999.;

  int idx_near_cathod_east = -1;
  int idx_near_cathod_west = -1;
  double min_abs_x_east = 999.;
  double min_abs_x_west = 999.;
  
  // == select only ACA tracks
  for (unsigned i = 0; i < sp_x.GetSize(); i++) {
    if(rr[i] > 0){
      double this_sp_x = sp_x[i];
      double this_sp_y = sp_y[i];
      double this_sp_z = sp_z[i];
      double this_dir_x = dir_x[i];
      double this_dir_y	= dir_y[i];
      double this_dir_z	= dir_z[i];
      double this_rr = rr[i];
      double this_dqdx = dqdx[i];
      double this_q_integ = q_integ[i];

      // == traj plots, post-SCE Corr.
      XYZVector sp_sce_uncorr(this_sp_x, this_sp_y, this_sp_z);
      XYZVector sp_sce_corr = sce_corr_mc -> WireToTrajectoryPosition(sp_sce_uncorr);
      double pitch_sce_uncorr = sce_corr_mc -> meas_pitch(this_sp_x, this_sp_y, this_sp_z, this_dir_x, this_dir_y, this_dir_z, plane, false);
      double pitch_sce_corr = sce_corr_mc -> meas_pitch(this_sp_x, this_sp_y, this_sp_z, this_dir_x, this_dir_y, this_dir_z, plane, true);
      double dqdx_sce_corr = this_dqdx * pitch_sce_uncorr / pitch_sce_corr;

      XYZVector EfieldOffset = sce_corr_mc -> GetEfieldOffsets(sp_sce_corr);
      double sce_efield = sce_corr_mc -> GetEfield(sp_sce_corr);

      /*
      if(pitch_sce_uncorr < 1.){
	cout << Form("%d: (x, y, z) = (%.2f, %.2f, %.2f), (rr, t) = (%.2f, %.2f)", i, sp_x[i], sp_y[i], sp_z[i], rr[i], time[i]) << endl;
      }
      */

      if(sp_x[i] > largest_x){
	largest_x = sp_x[i];
	idx_largest_x = i;
      }
      if(sp_x[i] < smallest_x){
	smallest_x = sp_x[i];
	idx_smallest_x = i;
      }

      if(sp_x[i] < 0.){
	if(fabs(sp_x[i]) < min_abs_x_east){
	  min_abs_x_east = fabs(sp_x[i]);
	  idx_near_cathod_east = i;
	}
      }
      if(sp_x[i] > 0.){
        if(fabs(sp_x[i]) < min_abs_x_west){
          min_abs_x_west = fabs(sp_x[i]);
          idx_near_cathod_west = i;
        }
      }
      
      if(fabs(sp_x[i]) < 5.){
	if(sp_x[i] < 0.) FillHist("aca_hit_zy_east_near_cathode_5cm_plane" + plane_str, sp_z[i], sp_y[i], 1., z_bins / 10, z_min, z_max, y_bins / 10, y_min, y_max);
	if(sp_x[i] > 0.) FillHist("aca_hit_zy_west_near_cathode_5cm_plane" + plane_str, sp_z[i], sp_y[i], 1., z_bins / 10, z_min, z_max, y_bins / 10, y_min, y_max);

	if(fabs(sp_x[i]) < 3.){
	  if(sp_x[i] < 0.) FillHist("aca_hit_zy_east_near_cathode_3cm_plane" + plane_str, sp_z[i], sp_y[i], 1., z_bins / 10, z_min, z_max, y_bins / 10, y_min, y_max);
	  if(sp_x[i] > 0.) FillHist("aca_hit_zy_west_near_cathode_3cm_plane" + plane_str, sp_z[i], sp_y[i], 1., z_bins / 10, z_min, z_max, y_bins / 10, y_min, y_max);

	  if(fabs(sp_x[i]) < 1.){
	    if(sp_x[i] < 0.) FillHist("aca_hit_zy_east_near_cathode_1cm_plane" + plane_str, sp_z[i], sp_y[i], 1., z_bins / 10, z_min, z_max, y_bins / 10, y_min, y_max);
	    if(sp_x[i] > 0.) FillHist("aca_hit_zy_west_near_cathode_1cm_plane" + plane_str, sp_z[i], sp_y[i], 1., z_bins / 10, z_min, z_max, y_bins / 10, y_min, y_max);
	  }
	}
      }
    }
  }

  if(idx_largest_x != -1 && idx_near_cathod_west != -1 && sp_x[idx_largest_x] > 195.){
    TString this_y_str = get_y_str(sp_y[idx_near_cathod_west]);
    FillHist("aca_hit_zx_west_near_anode_plane" + plane_str, sp_z[idx_near_cathod_west], sp_x[idx_largest_x], 1., z_bins / 40, z_min, z_max, 100., 190., 210.);
    FillHist("aca_hit_zx_west_near_anode_plane" + plane_str + "_y" + this_y_str, sp_z[idx_near_cathod_west], sp_x[idx_largest_x], 1., z_bins / 40, z_min, z_max, 100., 190., 210.);
    FillHist("aca_hit_zx_west_near_anode_plane" + plane_str + "_y" + this_y_str, sp_z[idx_near_cathod_west], sp_x[idx_largest_x], 1., 100, -100., 600., 500., -250., 250.);

  }
  if(idx_smallest_x != -1 && idx_near_cathod_east != -1 && sp_x[idx_smallest_x] < -195.){
    TString this_y_str = get_y_str(sp_y[idx_near_cathod_east]);
    FillHist("aca_hit_zx_east_near_anode_plane" + plane_str, sp_z[idx_near_cathod_east], sp_x[idx_smallest_x], 1., z_bins / 40, z_min, z_max, 100., -210., -190.);
    FillHist("aca_hit_zx_east_near_anode_plane" + plane_str + "_y" + this_y_str, sp_z[idx_near_cathod_east], sp_x[idx_smallest_x], 1., z_bins / 40, z_min, z_max, 100., -210., -190.);
    FillHist("aca_hit_zx_east_near_anode_plane" + plane_str + "_y" + this_y_str, sp_z[idx_near_cathod_east], sp_x[idx_smallest_x], 1., 100, -100., 600., 500., -250., 250.);
  }
}

void run_tpc_asym_study(TString list_file, TString out_suffix, bool IsData = false) {

  sce_corr_mc -> ReadHistograms();
  isdata = IsData;
  
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

  TString fileListPath = sample_list_dir + "/" + list_file;
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
  TTreeReaderValue<int> whicht0(myReader, "trk.whicht0");
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
  TTreeReaderArray<float> qinteg0(myReader, "trk.hits0.h.integral");
  TTreeReaderArray<float> qinteg1(myReader, "trk.hits1.h.integral");
  TTreeReaderArray<float> qinteg2(myReader, "trk.hits2.h.integral");
  TTreeReaderValue<float> dirx(myReader, "trk.dir.x"); 
  TTreeReaderValue<float> diry(myReader, "trk.dir.y");
  TTreeReaderValue<float> dirz(myReader, "trk.dir.z");
  TTreeReaderArray<float> dirx0(myReader, "trk.hits0.dir.x");
  TTreeReaderArray<float> diry0(myReader, "trk.hits0.dir.y");
  TTreeReaderArray<float> dirz0(myReader, "trk.hits0.dir.z");
  TTreeReaderArray<float> dirx1(myReader, "trk.hits1.dir.x");
  TTreeReaderArray<float> diry1(myReader, "trk.hits1.dir.y");
  TTreeReaderArray<float> dirz1(myReader, "trk.hits1.dir.z");
  TTreeReaderArray<float> dirx2(myReader, "trk.hits2.dir.x");
  TTreeReaderArray<float> diry2(myReader, "trk.hits2.dir.y");
  TTreeReaderArray<float> dirz2(myReader, "trk.hits2.dir.z");
  
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

  double ADC_med_cut = 1600.; // == https://sbn-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=23472&filename=SBND%20Calib%20Workshop%202021.pdf&version=1
  double track_length_cut = 15.;

  // Loop over all entries of the TTree
  //int N_run = 2000;
  //int N_run = 100000;
   
  while (myReader.Next()) {
    //if(current_entry > N_run) break;
   
    if(current_entry%100 == 0){
      cout << current_entry << " / " << N_entries << endl;
    }
    current_entry++;

    hist_selected -> Fill(*selected);

    double end_meddqdx = Get_EndMediandQdx(rr2, dqdx2);
    FillHist("end_meddqdx", end_meddqdx, 1., 5000., 0., 5000.);
    FillHist(Form("end_meddqdx_selected%d", *selected), end_meddqdx, 1., 5000., 0., 5000.);

    // == select TPC t0s
    if (*whicht0 == 0){
      // == Tracks selected as cathode passing through going
      if (*selected == 1) {
	
	//if(current_entry < N_run){
	  // == 1st ind plane
	  if(evt_sel(sp_x0, sp_y0, sp_z0, rr0, dqdx0)){
	    make_traj_plot(current_entry, 0, sp_x0, sp_y0, sp_z0, dirx0, diry0, dirz0, rr0, dqdx0, qinteg0, time0);
	  }
	  // == 2nd ind plane
	  if(evt_sel(sp_x1, sp_y1, sp_z1, rr1, dqdx1)){
	    make_traj_plot(current_entry, 1, sp_x1, sp_y1, sp_z1, dirx1, diry1, dirz1, rr1, dqdx1, qinteg1, time1);
	  }
	  // == collection plane
	  if(evt_sel(sp_x2, sp_y2, sp_z2, rr2, dqdx2)){
	    make_traj_plot(current_entry, 2, sp_x2, sp_y2, sp_z2, dirx2, diry2, dirz2, rr2, dqdx2, qinteg2, time2);
	  }
	  //}
      }
    }
  }
  TString output_rootfile_dir = getenv("OUTPUTROOT_PATH");
  TString output_file_name = output_rootfile_dir + "/output_tpc_asym_" + out_suffix + ".root";
  out_rootfile = new TFile(output_file_name, "RECREATE");
  out_rootfile -> cd();

  hist_selected -> Write();
  WriteHist();
  WriteCanvas();
  WriteGrs();
  out_rootfile -> Close();
}
