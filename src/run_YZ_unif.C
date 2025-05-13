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

const double y_min = -200, y_max = 200, z_min = 0, z_max = 500;
const double pixel_size = 2;
const int y_bins = (y_max - y_min) / pixel_size;
const int z_bins = (z_max - z_min) / pixel_size;

const double cos_min = -1., cos_max = 1.;
const double cos_size = 0.01;
const int cos_bins = (cos_max - cos_min) / cos_size;

bool InVeto_region_eastTPC_C(double y, double z) {
  // format considering {y_min, y_max, z_min, z_max}
  vector<vector<double>> vetoed_regions = {
    {-200, -185, 0, 15},
    {-200, -195, 240, 260},
    {-200, -195, 445, 455},
    {-200, -185, 495, 500},
    {190, 200, 460, 500},
    {185, 200, 230, 260},
    {180, 200, 0, 15},
    {-55, 45, 0, 10}
  };

  for (size_t i = 0; i < vetoed_regions.size(); ++i) {
    vector<double> region = vetoed_regions[i];

    double y_min = region[0];
    double y_max = region[1];
    double z_min = region[2];
    double z_max = region[3];

    if (y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
      return true;
    }
  }
  return false;
}

bool InVeto_region_westTPC_C(double y, double z) {
  // format considering {y_min, y_max, z_min, z_max}
  vector<vector<double>> vetoed_regions = {
    {-200, -185, 0, 15},
    {-200, -195, 230, 260},
    {-200, -180, 485, 500},
    {190, 200, 490, 500},
    {190, 200, 235, 275},
    {190, 200, 0, 15},
    //{-5, 5, 0, 15}
    {-50, 35, 0, 10}
  };

  for (size_t i = 0; i < vetoed_regions.size(); ++i) {
    vector<double> region = vetoed_regions[i];

    double y_min = region[0];
    double y_max = region[1];
    double z_min = region[2];
    double z_max = region[3];

    if (y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
      return true;
    }
  }
  return false;
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

  if(!Is_Edge(first_x, first_y, first_z) || !Is_Edge(last_x, last_y, last_z)) return out; // == through going track
  if(first_x * last_x < 0.) passing_cathode = true;
  if(!passing_cathode) return out;

  /*
  double end_meddqdx = Get_EndMediandQdx(rr, dqdx);
  if(end_meddqdx > 1500.) continue;
  */

  out = true;
  return out;
}


void get_cos_vals(const TTreeReaderArray<float> &sp_x, const TTreeReaderArray<float> &sp_y, const TTreeReaderArray<float> &sp_z, const TTreeReaderArray<float> &rr, Double_t *cos_vals){

  unsigned N_reco_hits = sp_x.GetSize();
  double first_x = -999.;
  double first_y = -999.;
  double first_z = -999.;
  double last_x = -999.;
  double last_y = -999.;
  double last_z = -999.;
  first_x = sp_x[N_reco_hits - 1];
  first_y = sp_y[N_reco_hits - 1];
  first_z = sp_z[N_reco_hits - 1];
  for (unsigned i = 0; i < sp_x.GetSize(); i++) {
    if(rr[i] > 0){
      last_x = sp_x[i];
      last_y = sp_y[i];
      last_z = sp_z[i];
      break;
    }
  }
  
  ROOT::Math::XYZVector track_vec(last_x - first_x, last_y - first_y, last_z - first_z);
  double cos_xy = track_vec.X() / (sqrt(pow(track_vec.X(), 2.) + pow(track_vec.Y(), 2.)));
  double cos_yz = track_vec.Y() / (sqrt(pow(track_vec.Y(), 2.) + pow(track_vec.Z(), 2.)));
  double cos_zx = track_vec.X() / (sqrt(pow(track_vec.X(), 2.) + pow(track_vec.Z(), 2.)));
  double zprime_plus = zprime_60deg(track_vec.Y(), track_vec.Z(), 1);
  double zprime_minus = zprime_60deg(track_vec.Y(), track_vec.Z(), -1);
  double cos_plus_zprimex = track_vec.X() / (sqrt(pow(track_vec.X(), 2.) + pow(zprime_plus, 2.)));
  double cos_minus_zprimex = track_vec.X() / (sqrt(pow(track_vec.X(), 2.) + pow(zprime_minus, 2.)));

  // == return four cos vals: cosyz, coszx, coszx+, coszx-
  cos_vals[0] = cos_yz;
  cos_vals[1] =	cos_zx;
  cos_vals[2] =	cos_plus_zprimex;
  cos_vals[3] =	cos_minus_zprimex;
}

void fill_meddqdx_cos(const std::vector<std::vector<std::vector<double>>> & dqdx, TString hist_name){

  FillHist(hist_name, 0., 0., 1., cos_bins, cos_min, cos_max, cos_bins, cos_min, cos_max);
  for (int xi = 0; xi < cos_bins; ++xi) {
    for (int yi = 0; yi < cos_bins; ++yi) {
      float local_median_dqdx=TMath::Median(dqdx[xi][yi].size(),&dqdx[xi][yi][0]); 
      maphist_TH2D[hist_name] -> SetBinContent(xi+1, yi+1, local_median_dqdx); 
    }
  }
}


void fill_meddqdx_yz(const std::vector<std::vector<std::vector<double>>> & dqdx, TString hist_name){

  cout << "[fill_meddqdx_yz] " << hist_name << endl;
  FillHist(hist_name + "_meddqdx", 1., 1., 1., z_bins, z_min, z_max, y_bins, y_min, y_max);
  vector<double> local_dqdx_medians;
  for (int yi = 0; yi < y_bins; ++yi) {
    for (int zi = 0; zi < z_bins; ++zi) {
      float local_median_dqdx = TMath::Median(dqdx[yi][zi].size(),&dqdx[yi][zi][0]);
      maphist_TH2D[hist_name + "_meddqdx"] -> SetBinContent(zi+1, yi+1, local_median_dqdx);

      local_dqdx_medians.push_back(local_median_dqdx);
    }
  }

  double global_median_dqdx = TMath::Median(local_dqdx_medians.size(), &local_dqdx_medians[0]);

  FillHist(hist_name + "_corr", 1., 1., 1., z_bins, z_min, z_max, y_bins, y_min, y_max);
  for (int yi = 0; yi < y_bins; ++yi) {
    for (int zi = 0; zi < z_bins; ++zi) {
      double this_median = maphist_TH2D[hist_name + "_meddqdx"] -> GetBinContent(zi+1, yi+1);
      double this_corr = global_median_dqdx / this_median;
      if(this_median < 1e-5) this_corr = 0.;
      maphist_TH2D[hist_name + "_corr"] -> SetBinContent(zi+1, yi+1, this_corr);
    }
  }
}

void run_YZ_unif(TString list_file, TString out_suffix, bool IsData = false) {

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
  int N_run = 200000;
   
  std::vector<std::vector<std::vector<double>>> cos_east_0(cos_bins, std::vector<std::vector<double>>(cos_bins));
  std::vector<std::vector<std::vector<double>>> cos_west_0(cos_bins, std::vector<std::vector<double>>(cos_bins));
  std::vector<std::vector<std::vector<double>>> cos_east_1(cos_bins, std::vector<std::vector<double>>(cos_bins));
  std::vector<std::vector<std::vector<double>>> cos_west_1(cos_bins, std::vector<std::vector<double>>(cos_bins));
  std::vector<std::vector<std::vector<double>>> cos_east_2(cos_bins, std::vector<std::vector<double>>(cos_bins));
  std::vector<std::vector<std::vector<double>>> cos_west_2(cos_bins, std::vector<std::vector<double>>(cos_bins));

  std::vector<std::vector<std::vector<double>>> dqdx_east_0(y_bins, std::vector<std::vector<double>>(z_bins)); // == x < 0
  std::vector<std::vector<std::vector<double>>> dqdx_west_0(y_bins, std::vector<std::vector<double>>(z_bins)); // == x > 0
  std::vector<std::vector<std::vector<double>>> dqdx_east_1(y_bins, std::vector<std::vector<double>>(z_bins)); // == x < 0
  std::vector<std::vector<std::vector<double>>> dqdx_west_1(y_bins, std::vector<std::vector<double>>(z_bins)); // == x > 0
  std::vector<std::vector<std::vector<double>>> dqdx_east_2(y_bins, std::vector<std::vector<double>>(z_bins)); // == x < 0
  std::vector<std::vector<std::vector<double>>> dqdx_west_2(y_bins, std::vector<std::vector<double>>(z_bins)); // == x > 0

  std::vector<std::vector<std::vector<double>>> dqdx_sce_east_0(y_bins, std::vector<std::vector<double>>(z_bins)); // == x < 0
  std::vector<std::vector<std::vector<double>>> dqdx_sce_west_0(y_bins, std::vector<std::vector<double>>(z_bins)); // == x > 0
  std::vector<std::vector<std::vector<double>>> dqdx_sce_east_1(y_bins, std::vector<std::vector<double>>(z_bins)); // == x < 0
  std::vector<std::vector<std::vector<double>>> dqdx_sce_west_1(y_bins, std::vector<std::vector<double>>(z_bins)); // == x > 0
  std::vector<std::vector<std::vector<double>>> dqdx_sce_east_2(y_bins, std::vector<std::vector<double>>(z_bins)); // == x < 0
  std::vector<std::vector<std::vector<double>>> dqdx_sce_west_2(y_bins, std::vector<std::vector<double>>(z_bins)); // == x > 0

  while (myReader.Next()) {
    if(current_entry > N_run) break;
   
    if(current_entry%100 == 0){
      cout << current_entry << " / " << N_entries << endl;
    }
    current_entry++;

    hist_selected -> Fill(*selected);

    double end_meddqdx = Get_EndMediandQdx(rr2, dqdx2);
    FillHist("end_meddqdx", end_meddqdx, 1., 5000., 0., 5000.);
    FillHist(Form("end_meddqdx_selected%d", *selected), end_meddqdx, 1., 5000., 0., 5000.);
    // == Tracks selected as cathode passing through going
    if (*selected == 2) {

      // == 1st ind plane
      if(evt_sel(sp_x0, sp_y0, sp_z0, rr0, dqdx0)){

	// -- cos vals: cosyz, coszx, coszx+, coszx-
	double cos_vals_0[4];
	get_cos_vals(sp_x0, sp_y0, sp_z0, rr0, cos_vals_0);
	
	for (size_t i = 0; i < sp_x0.GetSize(); ++i) {
	  if(rr0[i] < 0.) continue;

	  int y_index = (sp_y0[i] - y_min) / pixel_size;
	  int z_index = (sp_z0[i] - z_min) / pixel_size;

	  XYZVector sp_sce_uncorr(sp_x0[i], sp_y0[i], sp_z0[i]);
	  XYZVector sp_sce_corr = sce_corr_mc -> WireToTrajectoryPosition(sp_sce_uncorr);
	  int y_sce_index = (sp_sce_corr.Y() - y_min) / pixel_size;
	  int z_sce_index = (sp_sce_corr.Z() - z_min) / pixel_size;
	  double pitch_sce_corr = sce_corr_mc -> meas_pitch(sp_x0[i], sp_y0[i], sp_z0[i], dirx0[i], diry0[i], dirz0[i], 0, true);
	  double dqdx_sce_corr = dqdx0[i] * pitch0[i] / pitch_sce_corr;
	  
	  if(sp_x0[i] < 0){
	    int cos_x_index = (cos_vals_0[2] - cos_min) / cos_size;
	    int cos_y_index = (cos_vals_0[0] - cos_min) / cos_size;
	    cos_east_0[cos_x_index][cos_y_index].push_back(dqdx0[i]);
	    FillHist("plane_0_cos_zx_plus_vs_dqdx_east", cos_vals_0[2], dqdx0[i], 1., 100., -1., 1., 3000., 0., 3000.);

	    if(fabs(cos_vals_0[2]) < 0.75){// && !InVeto_region_eastTPC_C(sp_y0[i], sp_z0[i])){
	      if (y_index >= 0 && y_index < y_bins && z_index >= 0 && z_index < z_bins) {
		dqdx_east_0[y_index][z_index].push_back(dqdx0[i]); 
	      }
	      if (y_sce_index >= 0 && y_sce_index < y_bins && z_sce_index >= 0 && z_sce_index < z_bins) {
                dqdx_sce_east_0[y_sce_index][z_sce_index].push_back(dqdx_sce_corr);
              }
	    }
	  }
	  else{
	    int cos_x_index = (cos_vals_0[3] - cos_min) / cos_size;
            int cos_y_index = (cos_vals_0[0] - cos_min) / cos_size;
            cos_west_0[cos_x_index][cos_y_index].push_back(dqdx0[i]);
	    FillHist("plane_0_cos_zx_minus_vs_dqdx_west", cos_vals_0[3], dqdx0[i], 1., 100., -1., 1., 3000., 0., 3000.);

	    if(fabs(cos_vals_0[3]) < 0.75){// && !InVeto_region_westTPC_C(sp_y0[i], sp_z0[i])){
              if (y_index >= 0 && y_index < y_bins && z_index >= 0 && z_index < z_bins) {
                dqdx_west_0[y_index][z_index].push_back(dqdx0[i]);
              }
	      if (y_sce_index >= 0 && y_sce_index < y_bins && z_sce_index >= 0 && z_sce_index < z_bins) {
                dqdx_sce_west_0[y_sce_index][z_sce_index].push_back(dqdx_sce_corr);
              }
	    }
	  }
	}
      }

      // == 2nd ind plane
      if(evt_sel(sp_x1, sp_y1, sp_z1, rr1, dqdx1)){
	// -- cos vals: cosyz, coszx, coszx+, coszx-  
        double cos_vals_1[4];
        get_cos_vals(sp_x1, sp_y1, sp_z1, rr1, cos_vals_1);

        for (size_t i = 0; i < sp_x1.GetSize(); ++i) {
	  if(rr1[i] < 0.) continue;

	  int y_index = (sp_y1[i] - y_min) / pixel_size;
          int z_index = (sp_z1[i] - z_min) / pixel_size;

          XYZVector sp_sce_uncorr(sp_x1[i], sp_y1[i], sp_z1[i]);
          XYZVector sp_sce_corr = sce_corr_mc -> WireToTrajectoryPosition(sp_sce_uncorr);
          int y_sce_index = (sp_sce_corr.Y() - y_min) / pixel_size;
          int z_sce_index = (sp_sce_corr.Z() - z_min) / pixel_size;
          double pitch_sce_corr = sce_corr_mc -> meas_pitch(sp_x1[i], sp_y1[i], sp_z1[i], dirx1[i], diry1[i], dirz1[i], 1, true);
          double dqdx_sce_corr = dqdx1[i] * pitch1[i] / pitch_sce_corr;

	  if(sp_x1[i] < 0){
	    int cos_x_index = (cos_vals_1[3] - cos_min) / cos_size;
            int cos_y_index = (cos_vals_1[0] - cos_min) / cos_size;
            cos_east_1[cos_x_index][cos_y_index].push_back(dqdx1[i]);
	    FillHist("plane_1_cos_zx_minus_vs_dqdx_east", cos_vals_1[3], dqdx1[i], 1., 100., -1., 1., 3000., 0., 3000.);

	    if(fabs(cos_vals_1[3]) < 0.75){// && !InVeto_region_eastTPC_C(sp_y1[i], sp_z1[i])){
              if (y_index >= 0 && y_index < y_bins && z_index >= 0 && z_index < z_bins) {
                dqdx_east_1[y_index][z_index].push_back(dqdx1[i]);
              }
              if (y_sce_index >= 0 && y_sce_index < y_bins && z_sce_index >= 0 && z_sce_index < z_bins) {
		dqdx_sce_east_1[y_sce_index][z_sce_index].push_back(dqdx_sce_corr);
              }
            }
	  }
          else{
	    int cos_x_index = (cos_vals_1[2] - cos_min) / cos_size;
            int cos_y_index = (cos_vals_1[0] - cos_min) / cos_size;
            cos_west_1[cos_x_index][cos_y_index].push_back(dqdx1[i]);
	    FillHist("plane_1_cos_zx_plus_vs_dqdx_west", cos_vals_1[2], dqdx1[i], 1., 100., -1., 1., 3000., 0., 3000.);

	    if(fabs(cos_vals_1[2]) < 0.75){// && !InVeto_region_westTPC_C(sp_y1[i], sp_z1[i])){
              if (y_index >= 0 && y_index < y_bins && z_index >= 0 && z_index < z_bins) {
                dqdx_west_1[y_index][z_index].push_back(dqdx1[i]);
              }
              if (y_sce_index >= 0 && y_sce_index < y_bins && z_sce_index >= 0 && z_sce_index < z_bins) {
                dqdx_sce_west_1[y_sce_index][z_sce_index].push_back(dqdx_sce_corr);
              }
            }
	  }
        }
      }
      
      // == collection plane
      if(evt_sel(sp_x2, sp_y2, sp_z2, rr2, dqdx2)){
	// -- cos vals: cosyz, coszx, coszx+, coszx-                                           
        double cos_vals_2[4];
        get_cos_vals(sp_x2, sp_y2, sp_z2, rr2, cos_vals_2);

        for (size_t i = 0; i < sp_x2.GetSize(); ++i) {
	  if(rr2[i] < 0.) continue;

	  int y_index = (sp_y2[i] - y_min) / pixel_size;
          int z_index = (sp_z2[i] - z_min) / pixel_size;

          XYZVector sp_sce_uncorr(sp_x2[i], sp_y2[i], sp_z2[i]);
          XYZVector sp_sce_corr = sce_corr_mc -> WireToTrajectoryPosition(sp_sce_uncorr);
          int y_sce_index = (sp_sce_corr.Y() - y_min) / pixel_size;
          int z_sce_index = (sp_sce_corr.Z() - z_min) / pixel_size;
          double pitch_sce_corr = sce_corr_mc -> meas_pitch(sp_x2[i], sp_y2[i], sp_z2[i], dirx2[i], diry2[i], dirz2[i], 1, true);
          double dqdx_sce_corr = dqdx2[i] * pitch2[i] / pitch_sce_corr;

	  if(sp_x2[i] < 0){
	    int cos_x_index = (cos_vals_2[1] - cos_min) / cos_size;
            int cos_y_index = (cos_vals_2[0] - cos_min) / cos_size;
            cos_east_2[cos_x_index][cos_y_index].push_back(dqdx2[i]);
	    FillHist("plane_2_cos_zx_vs_dqdx_east", cos_vals_2[1], dqdx2[i], 1., 100., -1., 1., 3000., 0., 3000.);

	    if(fabs(cos_vals_2[1]) < 0.75){// && !InVeto_region_eastTPC_C(sp_y2[i], sp_z2[i])){
              if (y_index >= 0 && y_index < y_bins && z_index >= 0 && z_index < z_bins) {
                dqdx_east_2[y_index][z_index].push_back(dqdx2[i]);
              }
              if (y_sce_index >= 0 && y_sce_index < y_bins && z_sce_index >= 0 && z_sce_index < z_bins) {
                dqdx_sce_east_2[y_sce_index][z_sce_index].push_back(dqdx_sce_corr);
              }
            }
	  }
	  else{
	    int cos_x_index = (cos_vals_2[1] - cos_min) / cos_size;
            int cos_y_index = (cos_vals_2[0] - cos_min) / cos_size;
            cos_west_2[cos_x_index][cos_y_index].push_back(dqdx2[i]);
	    FillHist("plane_2_cos_zx_vs_dqdx_west", cos_vals_2[1], dqdx2[i], 1., 100., -1., 1., 3000., 0., 3000.);

	    if(fabs(cos_vals_2[1]) < 0.75){// && !InVeto_region_westTPC_C(sp_y2[i], sp_z2[i])){
              if (y_index >= 0 && y_index < y_bins && z_index >= 0 && z_index < z_bins) {
                dqdx_west_2[y_index][z_index].push_back(dqdx2[i]);
              }
              if (y_sce_index >= 0 && y_sce_index < y_bins && z_sce_index >= 0 && z_sce_index < z_bins) {
                dqdx_sce_west_2[y_sce_index][z_sce_index].push_back(dqdx_sce_corr);
              }
            }
	  }
	}
      }
    }
  }

  // == median dqdx for angles
  fill_meddqdx_cos(cos_east_0, "cos_meddqdx_plane0_east_zx_plus");
  fill_meddqdx_cos(cos_west_0, "cos_meddqdx_plane0_west_zx_minus");
  fill_meddqdx_cos(cos_east_1, "cos_meddqdx_plane1_east_zx_minus");
  fill_meddqdx_cos(cos_west_1, "cos_meddqdx_plane1_west_zx_plus");
  fill_meddqdx_cos(cos_east_2, "cos_meddqdx_plane2_east_zx");
  fill_meddqdx_cos(cos_west_2, "cos_meddqdx_plane2_west_zx");

  fill_meddqdx_yz(dqdx_east_0, "yz_plane0_east");
  fill_meddqdx_yz(dqdx_east_1, "yz_plane1_east");
  fill_meddqdx_yz(dqdx_east_2, "yz_plane2_east");
  fill_meddqdx_yz(dqdx_west_0, "yz_plane0_west");
  fill_meddqdx_yz(dqdx_west_1, "yz_plane1_west");
  fill_meddqdx_yz(dqdx_west_2, "yz_plane2_west");

  fill_meddqdx_yz(dqdx_sce_east_0, "yz_plane0_east_sce");
  fill_meddqdx_yz(dqdx_sce_east_1, "yz_plane1_east_sce");
  fill_meddqdx_yz(dqdx_sce_east_2, "yz_plane2_east_sce");
  fill_meddqdx_yz(dqdx_sce_west_0, "yz_plane0_west_sce");
  fill_meddqdx_yz(dqdx_sce_west_1, "yz_plane1_west_sce");
  fill_meddqdx_yz(dqdx_sce_west_2, "yz_plane2_west_sce");

  TString output_rootfile_dir = getenv("OUTPUTROOT_PATH");
  TString output_file_name = output_rootfile_dir + "/output_YZ_unif_" + out_suffix + ".root";
  out_rootfile = new TFile(output_file_name, "RECREATE");
  out_rootfile -> cd();

  hist_selected -> Write();
  WriteHist();

  out_rootfile -> Close();
}
