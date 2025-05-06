#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "mylib.h"
#include "Math/Vector3D.h"
#include "BetheBloch.h"
#include "Point3D.h"
#include "SCECorr.h"

BetheBloch *muon_BB = new BetheBloch(13);

SCECorr *sce_corr_mc = new SCECorr(false);

bool isdata = false;

const double y_min = -200, y_max = 200, z_min = 0, z_max = 500;
const double pixel_size = 5;
const int y_bins = (y_max - y_min) / pixel_size;
const int z_bins = (z_max - z_min) / pixel_size;

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
  TTreeReaderValue<float> dirx(myReader, "trk.dir.x"); 
  TTreeReaderValue<float> diry(myReader, "trk.dir.y");
  TTreeReaderValue<float> dirz(myReader, "trk.dir.z");
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
  int N_run = 10;
  std::vector<std::vector<std::vector<double>>> dqdx_east(y_bins, std::vector<std::vector<double>>(z_bins)); // == x < 0
  std::vector<std::vector<std::vector<double>>> dqdx_west(y_bins, std::vector<std::vector<double>>(z_bins)); // == x > 0

  while (myReader.Next()) {
    if(current_entry > N_run) break;
   
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

    double end_meddqdx = Get_EndMediandQdx(rr2, dqdx2);
    FillHist("end_meddqdx", end_meddqdx, 1., 5000., 0., 5000.);
    FillHist(Form("end_meddqdx_selected%d", *selected), end_meddqdx, 1., 5000., 0., 5000.);
    // == Tracks selected as cathode passing through going
    if (*selected == 2) {

      //cout << "(*selected == 2)" << endl;
      
      unsigned N_reco_hits = rr2.GetSize();

      if(N_reco_hits < 3) continue;
            
      TVector3 this_reco_start(sp_x2[N_reco_hits - 1], sp_y2[N_reco_hits - 1], sp_z2[N_reco_hits - 1]);
      TVector3 this_reco_end(sp_x2[0], sp_y2[0], sp_z2[0]);
      TVector3 this_true_start(true_start_x[0], true_start_y[0], true_start_z[0]);
      TVector3 this_true_end(true_end_x[0], true_end_y[0], true_end_z[0]);

      double dist_start = (this_reco_start - this_true_start).Mag();
      double dist_end = (this_reco_end - this_true_end).Mag();

      double this_reco_trk_len = rr2[N_reco_hits - 1];

      // == Pick first and last hits to check if the track is passing the cathode
      double first_x = -999.;
      double first_y = -999.;
      double first_z = -999.;
      double last_x = -999.;
      double last_y = -999.;
      double last_z = -999.;

      //cout << "get first and last hits" << endl;
      first_x = sp_x2[N_reco_hits - 1];
      first_y = sp_y2[N_reco_hits - 1];
      first_z = sp_z2[N_reco_hits - 1];
      for (unsigned i = 0; i < dqdx2.GetSize(); i++) {
        if(rr2[i] > 0){
	  last_x = sp_x2[i];
	  last_y = sp_y2[i];
          last_z = sp_z2[i];
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

      if(!Is_Edge(first_x, first_y, first_z) || !Is_Edge(last_x, last_y, last_z)) continue; // == through going track
      if(first_x * last_x < 0.) passing_cathode = true;
      if(!passing_cathode) continue;

      double thetaxz = acos(*dirz / sqrt(pow(*dirz,2)+pow(*dirx,2)))*180/TMath::Pi();
      if(*dirx<0) thetaxz = -thetaxz;
      double thetayz = acos(*dirz / sqrt(pow(*dirz,2)+pow(*diry,2)))*180/TMath::Pi();
      if(*diry<0) thetayz = -thetayz;

      if(abs(thetaxz)<115&&abs(thetaxz)>65) continue;//Angle
      if(abs(thetayz)<110&&abs(thetayz)>70) continue;//Angle
      
      if(end_meddqdx > 1500.) continue;

      FillHist("end_meddqdx_yz_sel", end_meddqdx, 1., 5000., 0., 5000.);
      FillHist("this_reco_trk_len", this_reco_trk_len, 1., 1000., 0., 1000.);

      for (size_t i = 0; i < sp_y2.GetSize(); ++i) {

	if(sp_x2[i] < 0 && InVeto_region_eastTPC_C(sp_y2[i], sp_z2[i])) continue;
	if(sp_x2[i] > 0 && InVeto_region_westTPC_C(sp_y2[i], sp_z2[i])) continue;
	
        int y_index = (sp_y2[i] - y_min) / pixel_size;
        int z_index = (sp_z2[i] - z_min) / pixel_size;
        if (y_index >= 0 && y_index < y_bins && z_index >= 0 && z_index < z_bins) {
	  if(sp_x2[i] < 0) dqdx_east[y_index][z_index].push_back(dqdx2[i]);
	  else dqdx_west[y_index][z_index].push_back(dqdx2[i]);

	  // == debugging pitch
	  //cout << Form("dir (x, y, z) = (%.2f, %.2f, %.2f)", dirx2[i], diry2[i], dirz2[i]) << endl;
	  double pitch_repro_sce_off = sce_corr_mc -> meas_pitch(sp_x2[i], sp_y2[i], sp_z2[i], dirx2[i], diry2[i], dirz2[i], 2, false);
	  double pitch_repro_sce_on = sce_corr_mc -> meas_pitch(sp_x2[i], sp_y2[i], sp_z2[i], dirx2[i], diry2[i], dirz2[i], 2, true);

	  cout << Form("pitch: %.5f, pitch_repro (SCE off): %.5f, pitch_repro (SCE on): %.5f", pitch2[i], pitch_repro_sce_off, pitch_repro_sce_on) << endl;

	}
      }
    }
  }


  vector<double> local_dqdx_medians_east, local_dqdx_medians_west;
  TH2D* h_yz_median_east = new TH2D("", "", z_bins, z_min, z_max, y_bins, y_min, y_max);
  TH2D* h_yz_median_west = new TH2D("", "", z_bins, z_min, z_max, y_bins, y_min, y_max);
  for (int yi = 0; yi < y_bins; ++yi) {
    for (int zi = 0; zi < z_bins; ++zi) {
      float local_median_dqdx_east=TMath::Median(dqdx_east[yi][zi].size(),&dqdx_east[yi][zi][0]);
      float local_median_dqdx_west=TMath::Median(dqdx_west[yi][zi].size(),&dqdx_west[yi][zi][0]);
      h_yz_median_east -> SetBinContent(zi+1, yi+1, local_median_dqdx_east);
      h_yz_median_west -> SetBinContent(zi+1, yi+1, local_median_dqdx_west);
      local_dqdx_medians_east.push_back(local_median_dqdx_east);
      local_dqdx_medians_west.push_back(local_median_dqdx_west);
    }
  }

  double global_median_dqdx_east = TMath::Median(local_dqdx_medians_east.size(), &local_dqdx_medians_east[0]);
  double global_median_dqdx_west = TMath::Median(local_dqdx_medians_west.size(), &local_dqdx_medians_west[0]);
  
  TH2D* h_yz_corr_east = (TH2D*)h_yz_median_east -> Clone();
  TH2D* h_yz_corr_west = (TH2D*)h_yz_median_west -> Clone();
  for (int yi = 0; yi < y_bins; ++yi) {
    for (int zi = 0; zi < z_bins; ++zi) {
      double this_east_median = h_yz_corr_east -> GetBinContent(zi+1, yi+1);
      double this_east_corr = global_median_dqdx_east / this_east_median;
      if(this_east_median < 1e-5) this_east_corr = 0.;
      h_yz_corr_east -> SetBinContent(zi+1, yi+1, this_east_corr);

      double this_west_median = h_yz_corr_west -> GetBinContent(zi+1, yi+1);
      double this_west_corr = global_median_dqdx_west /	this_west_median;
      if(this_west_median < 1e-5) this_west_corr = 0.;
      h_yz_corr_west ->	SetBinContent(zi+1, yi+1, this_west_corr);
    }
  }
  
  TString output_rootfile_dir = getenv("OUTPUTROOT_PATH");
  TString output_file_name = output_rootfile_dir + "/output_YZ_unif_" + out_suffix + ".root";
  out_rootfile = new TFile(output_file_name, "RECREATE");
  out_rootfile -> cd();

  hist_selected -> Write();
  WriteHist();

  h_yz_median_east -> SetName("yz_median_east");
  h_yz_median_west -> SetName("yz_median_west");
  h_yz_median_east -> Write();
  h_yz_median_west -> Write();

  h_yz_corr_east -> SetName("yz_corr_east");
  h_yz_corr_west -> SetName("yz_corr_west");
  h_yz_corr_east -> Write();
  h_yz_corr_west -> Write();
  
  out_rootfile -> Close();
}
