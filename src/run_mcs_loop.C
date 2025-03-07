#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "mylib.h"
#include "Math/Vector3D.h"
#include "BetheBloch.h"

bool isdata = false;
TSpline3 * muon_sp_range_to_KE = Get_sp_range_KE(mass_muon);

double muon_ek2p(double ek){

  double out = pow(pow(ek + mass_muon, 2.) - mass_muon * mass_muon, 0.5);
  return out;
}

// == mimic functions in https://github.com/SBNSoftware/sbncode/blob/develop/sbncode/LArRecoProducer/LArReco/TrajectoryMCSFitter.cxx
bool isInVolume(TVector3 point){

  bool out = false;
  if(fabs(point.X()) < 190. && fabs(point.Y()) < 190 && point.Z() > 10. && point.Z() < 490.) out = true;
  return out;
}


// == some default parameters based on https://github.com/LArSoft/larreco/blob/develop/larreco/TrackFinder/mcsfitproducer.fcl
const double segLen_ = 14.;
const int minNSegs_ = 3;
const int minHitsPerSegment_ = 2;
const double lar_radl_inv = 1./14.0;

// == order of hits is different from traj. Calibntuple's hit index starts from the last hit.
void breakTrajInSegments(const TTreeReaderArray<float>& sp_x, const TTreeReaderArray<float>& sp_y, const TTreeReaderArray<float>& sp_z, const TTreeReaderArray<float>& rr,
			 vector<size_t>& breakpoints, vector<float>& segradlengths, vector<float>& cumseglens, vector<bool>& breakpointsgood) {

  //cout << "[breakTrajInSegments] start" << endl;
  const unsigned N_reco_hits = rr.GetSize();
  if(N_reco_hits < 2) return;

  const double trajlen = rr[N_reco_hits - 1];
  const int nseg = std::max(minNSegs_,int(trajlen/segLen_));
  //const double thisSegLen = trajlen/double(nseg); // -- it does not make sense, we should keep segLen_. last segment would have length shorter then 14 cm and we can add flag if it is complete segLen_ long or not to each segment
  
  cumseglens.push_back(0.);
  int lastValid = 0;
  for(int i = N_reco_hits - 1; i >= 0; i--){
    if(rr[i] < 0.){
      //TVector3 this_pos(sp_x[i], sp_y[i], sp_z[i]);
      //if(isInVolume(this_pos)){
      lastValid = i + 1;
      break;
      //}
    }
  }

  int firstValid = 0;
  for(int i = N_reco_hits - 1; i >= 0; i--){
    if(rr[i] > 0.){
      TVector3 this_pos(sp_x[i], sp_y[i], sp_z[i]);
      if(isInVolume(this_pos)){
	firstValid = i;
	break;
      }
    }
  }

  //cout << "[breakTrajInSegments] N_reco_hits: " << N_reco_hits << ", firstValid: " << firstValid << ", lastValid: " << lastValid << endl;
  
  if(firstValid < 1) return;
  
  breakpoints.push_back(firstValid);
  TVector3 pos0(sp_x[firstValid], sp_y[firstValid], sp_z[firstValid]);
  bool pos0good = isInVolume(pos0); // -- FV cut only. No additional volume requiremet yet && !isInVolume(excludeVolumes, pos0);
  breakpointsgood.push_back(pos0good);
  double thislen = 0.;
  int npoints = 0;
  for(int i = firstValid; i >= lastValid ; i--){
    //cout << i << ", thislen : " << thislen << endl;
    TVector3 pos1(sp_x[i], sp_y[i], sp_z[i]);
    thislen += sqrt((pos1 - pos0).Mag2());
    pos0 = pos1;
    npoints++;
    if(thislen > segLen_){ // == use segLen_ instead of thisSegLen
      bool pos1good = isInVolume(pos1);
      breakpoints.push_back(i);
      if (npoints>=minHitsPerSegment_) segradlengths.push_back(thislen*lar_radl_inv);
      else segradlengths.push_back(-999.);
      cumseglens.push_back(cumseglens.back()+thislen);
      thislen = 0.;
      npoints = 0;
    }
    //cout << i << ", thislen : " << thislen << endl;
  }

  if (thislen>0.) {
    TVector3 endpointpos(sp_x[lastValid], sp_y[lastValid], sp_z[lastValid]);
    bool endpointposgood = isInVolume(endpointpos);
    breakpoints.push_back(N_reco_hits);
    breakpointsgood.push_back(endpointposgood);
    segradlengths.push_back(thislen*lar_radl_inv);
    cumseglens.push_back(cumseglens.back()+thislen);
  }
  return;
}

TVector3 avg_point_calculator(vector<TVector3> pos_vec){

  double avg_x = 0., avg_y = 0., avg_z = 0.;

  int n_pos = pos_vec.size();
  for(int i = 0; i < n_pos; i++){
    double this_x = pos_vec.at(i).X();
    double this_y = pos_vec.at(i).Y();
    double this_z = pos_vec.at(i).Z();

    avg_x += this_x;
    avg_y += this_y;
    avg_z += this_z;
  }

  avg_x = avg_x / (n_pos + 0.);
  avg_y = avg_y / (n_pos + 0.);
  avg_z = avg_z / (n_pos + 0.);

  TVector3 out(avg_x, avg_y, avg_z);
  return out;
}

void linearRegression(const TTreeReaderArray<float>& sp_x, const TTreeReaderArray<float>& sp_y, const TTreeReaderArray<float>& sp_z, const TTreeReaderArray<float>& rr, TVector3 track_vec,
					   const size_t firstPoint, const size_t lastPoint, TVector3& pcdir) {

  //cout << "[linearRegression] start" << endl;
  int npoints = 0;
  vector<TVector3> pos_vec;
  for(int i = firstPoint; i > lastPoint; i--){
    TVector3 this_pos(sp_x[i], sp_y[i], sp_z[i]);
    pos_vec.push_back(this_pos);
    npoints++;
  }

  TVector3 avgpos = avg_point_calculator(pos_vec);
  const double norm = 1./double(npoints);

  TMatrixDSym m(3);

  for(int i = firstPoint; i > lastPoint; i--){
    TVector3 p(sp_x[i], sp_y[i], sp_z[i]);
    const double xxw0 = p.X()-avgpos.X();
    const double yyw0 = p.Y()-avgpos.Y();
    const double zzw0 = p.Z()-avgpos.Z();
    m(0, 0) += xxw0*xxw0*norm;
    m(0, 1) += xxw0*yyw0*norm;
    m(0, 2) += xxw0*zzw0*norm;
    m(1, 0) += yyw0*xxw0*norm;
    m(1, 1) += yyw0*yyw0*norm;
    m(1, 2) += yyw0*zzw0*norm;
    m(2, 0) += zzw0*xxw0*norm;
    m(2, 1) += zzw0*yyw0*norm;
    m(2, 2) += zzw0*zzw0*norm;
  }

  const TMatrixDSymEigen me(m);
  const auto& eigenval = me.GetEigenValues();
  const auto& eigenvec = me.GetEigenVectors();
  //
  int maxevalidx = 0;
  double maxeval = eigenval(0);
  for (int i=1; i<3; ++i) {
    if (eigenval(i)>maxeval) {
      maxevalidx = i;
      maxeval = eigenval(i);
    }
  } 
  //
  pcdir = TVector3(eigenvec(0, maxevalidx),eigenvec(1, maxevalidx),eigenvec(2, maxevalidx));
  if (track_vec.Dot(pcdir)<0.) pcdir*=-1.;
}

TVector3 RotateToZaxis(TVector3 reference_vec, TVector3 input_vec){
  double theta = reference_vec.Theta();
  TVector3 ortho_line_on_xy(1., -1. * reference_vec.X() / reference_vec.Y(), 0.);
  TVector3 cross = reference_vec.Cross(ortho_line_on_xy);
  if(cross.Z() < 0.) theta = -1. * theta;
  TVector3 out(input_vec);
  out.Rotate(-1 * theta, ortho_line_on_xy);
  return out;
}

double dqdx_scale_correction_angle(double theta){

  // theta shoud be deg
  double p0 = -1.55222;
  double p1 = 0.0500101;
  double p2 = -0.000211176;

  double this_bias = (p0 + p1 * theta + p2 * theta * theta) / 100.; // -- % to number
  double this_correction = 1. / (1. - this_bias);

  return this_correction;
}


double zprime_60deg(double y, double z, int pm = 1){
  // == This function is for angle cut for induction planes.
  // ==== Rotation on yz plane to provide z' for cos(theta_z'x) measurement
  double cos_60deg = 0.5;
  double sig_60deg = sqrt(3.) * 0.5;
  double zprime = z * cos_60deg - (pm + 0.) * sig_60deg * y;
  return zprime;
}

double Get_EndMediandQdx(const TTreeReaderArray<float>& rr, const TTreeReaderArray<float>& dqdx){
  // == This function is for a study with stronger cut on EndMediandQddx for stopping track selection
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

TString Get_theta_trk_x_bin_str(double theta_trk_x){
  // == This is for getting bin string for theta_trk_x, following ICARUS's binning (https://arxiv.org/pdf/2407.12969)
  if(theta_trk_x < 30.) return "0to30";
  else if(theta_trk_x < 80.){
    int tenth = theta_trk_x / 10;
    TString this_out = Form("%dto%d", tenth * 10, (tenth+1) * 10);
    return this_out;
  }
  else if(theta_trk_x < 85.) return "80to85";
  else if(theta_trk_x < 90.) return "85to90";
  else return "Error";

  return "Error";
}

void Fill_corrected_dqdx_plots(TString suffix, const TTreeReaderArray<float>& rr, const TTreeReaderArray<float>& dqdx, const TTreeReaderArray<float>& sp_x, const TTreeReaderArray<float>& sp_z, const TTreeReaderArray<float>& pitch, double last_x, double cos_plus_zprimex, double cos_minus_zprimex, double theta_trk_x, TString theta_trk_x_str, bool do_ind_ang_cut = false){
  // == Fill plots for recombination fits
  if(dqdx.GetSize() < 1) return;
  for (unsigned i = 1; i < dqdx.GetSize() - 1; i++) { // == Not using first and last hits of a track
    if(rr[i] < 0.) continue;
    if(do_ind_ang_cut){ // == angle cut for induction planes
      if(suffix.Contains("plane0")){
        if(sp_x[i] < 0. && fabs(cos_plus_zprimex) > 0.75) continue;
        if(sp_x[i] > 0. && fabs(cos_minus_zprimex) > 0.75) continue;
      }
      if(suffix.Contains("plane1")){
        if(sp_x[i] < 0. && fabs(cos_minus_zprimex) > 0.75) continue;
        if(sp_x[i] > 0. && fabs(cos_plus_zprimex) > 0.75) continue;
      }
    }

    double this_lifetime_corr = Lifetime_Correction(sp_x[i], 100.0);
    if(isdata) this_lifetime_corr = 1.; // == FIXME, for data, do not apply lifetime correction. Should be updated in future to use different lifetime values for MC and data
    double corrected_dqdx = dqdx[i] * this_lifetime_corr;
    double this_dqdx_bias_corr = dqdx_scale_correction_angle(theta_trk_x);
    //cout << "corrected_dqdx: " << corrected_dqdx << endl;
    corrected_dqdx *= this_dqdx_bias_corr;
    //cout << "bias corrected_dqdx: " << corrected_dqdx << endl;

    FillHist("rr_vs_corr_dqdx_" + suffix, rr[i], corrected_dqdx, 1., 300., 0., 300., 3000., 0., 3000.);
    FillHist("rr_vs_pitch_" + suffix, rr[i], pitch[i], 1., 300., 0., 300., 200., 0., 2.);
    FillHist("pitch_" + suffix, pitch[i], 1., 200., 0., 2.);
    FillHist("pitch_x_vs_corr_dqdx_" + suffix, pitch[i], corrected_dqdx, 1., 200., 0., 2., 3000., 0., 3000.);
    if(i == 0) FillHist("last_x_" + suffix, last_x, 1., 500., -250., 250.);

    if(pitch[i] > 1.) continue;
    
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
    FillHist("dEdx_MPV_vs_corr_dqdx_" + suffix + "_phi" + theta_trk_x_str, this_dEdx_MPV, corrected_dqdx, 1., 3000., 0., 30., 3000., 0., 3000.);
    //cout << "[Fill_corrected_dqdx_plots] " << i << ", rr : " << rr[i] << ", KE : " << muon_sp_range_to_KE -> Eval(rr[i]) << ", this_kappa : " << this_kappa << ", this_dEdx_MPV : " << this_dEdx_MPV << endl;

    // == Divide into NE, NW, SE and SW
    if(sp_x[i] < 0.){
      if(sp_z[i] > 250.) FillHist("dEdx_MPV_vs_corr_dqdx_" + suffix + "_NE", this_dEdx_MPV, corrected_dqdx, 1., 3000., 0., 30., 3000., 0., 3000.);
      else FillHist("dEdx_MPV_vs_corr_dqdx_" + suffix + "_SE", this_dEdx_MPV, corrected_dqdx, 1., 3000., 0., 30., 3000., 0., 3000.);
    }
    else{
      if(sp_z[i] > 250.) FillHist("dEdx_MPV_vs_corr_dqdx_" + suffix + "_NW", this_dEdx_MPV, corrected_dqdx, 1., 3000., 0., 30., 3000., 0., 3000.);
      else FillHist("dEdx_MPV_vs_corr_dqdx_" + suffix + "_SW", this_dEdx_MPV, corrected_dqdx, 1., 3000., 0., 30., 3000., 0., 3000.);
    }
  }
}

void run_mcs_loop(int run_number = 0) {

  TString run_number_str = "";
  run_number_str = TString::Format("%d", run_number);
  if(run_number > 0){
    isdata = true;
  }
  
  /////////////////////////////////
  // == Call Trees
  /////////////////////////////////
  // Open the file containing the tree
  TChain *fChain = new TChain("caloskim/TrackCaloSkim");
  TString input_file_dir = getenv("DATA_PATH");
  TString sample_list_dir = getenv("SAMPLE_PATH");
  TString sample_list_label = getenv("FILELIST_LABEL");

  TString fileListPath = sample_list_dir + "/list" + sample_list_label + run_number_str + ".txt";
  if(!isdata) fileListPath = sample_list_dir + "/calib_ntuple_moon_v10_04_1.list";
  //if(!isdata) fileListPath = sample_list_dir + "/calib_ntuple_moon_v10_04_1_single_file.list";
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

  // Loop over all entries of the TTree
  int _run_to = 370000;
  //_run_to = 1000;
  while (myReader.Next()) {
    if(current_entry > _run_to) break;
   
    if(current_entry%100 == 0){
      cout << current_entry << " / " << N_entries << endl;
    }
    current_entry++;

    FillHist("selected", *selected, 1., 40, -0.5, 39.5);

    double end_meddqdx = Get_EndMediandQdx(rr2, dqdx2);
    FillHist("end_meddqdx", end_meddqdx, 1., 5000., 0., 5000.);
    FillHist(Form("end_meddqdx_selected%d", *selected), end_meddqdx, 1., 5000., 0., 5000.);

    // == Tracks selected as stopping
    if (*selected == 0) {
      
      unsigned N_reco_hits = rr2.GetSize();
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
      TVector3 track_vec(last_x - first_x, last_y - first_y, last_z - first_z);
      double cos_xy = track_vec.X() / (sqrt(pow(track_vec.X(), 2.) + pow(track_vec.Y(), 2.)));
      double cos_yz = track_vec.Y() / (sqrt(pow(track_vec.Y(), 2.) + pow(track_vec.Z(), 2.)));
      double cos_zx = track_vec.X() / (sqrt(pow(track_vec.X(), 2.) + pow(track_vec.Z(), 2.)));
      double zprime_plus = zprime_60deg(track_vec.Y(), track_vec.Z(), 1);
      double zprime_minus = zprime_60deg(track_vec.Y(), track_vec.Z(), -1);
      double cos_plus_zprimex = track_vec.X() / (sqrt(pow(track_vec.X(), 2.) + pow(zprime_plus, 2.)));
      double cos_minus_zprimex = track_vec.X() / (sqrt(pow(track_vec.X(), 2.) + pow(zprime_minus, 2.)));
      TVector3 x_unit(1.0, 0., 0.);
      TVector3 trk_cross_x = track_vec.Cross(x_unit);
      double sin_theta_trk_x = sqrt(trk_cross_x.Mag2()) / sqrt(track_vec.Mag2());
      double theta_trk_x = TMath::ASin(sin_theta_trk_x) * 180. / TMath::Pi(); // == [Deg.]
      TString theta_trk_x_str = Get_theta_trk_x_bin_str(theta_trk_x);
      
      if(first_x * last_x < 0.) passing_cathode= true;
      if(this_reco_trk_len < 100. || !passing_cathode) continue;

      vector<size_t> breakpoints;
      vector<float> segradlengths;
      vector<float> cumseglens;
      vector<bool> breakpointsgood;
      breakTrajInSegments(sp_x2, sp_y2, sp_z2, rr2, breakpoints, segradlengths, cumseglens, breakpointsgood);
      //cout << "segradlengths.size(): " << segradlengths.size() << endl;
      if(segradlengths.size()<3) continue;
      vector<float> theta_yz_vec;
      vector<float> theta_xz_vec;
      vector<float> dtheta;
      vector<float> pmuon_vec;
      TVector3 pcdir0;
      TVector3 pcdir1;
      for (unsigned int p = 0; p<segradlengths.size() - 1; p++) {
	linearRegression(sp_x2, sp_y2, sp_z2, rr2, track_vec, breakpoints[p], breakpoints[p+1], pcdir1);


	if (p>0) {
	  /*
	  if (segradlengths[p]<-100. || segradlengths[p-1]<-100.) {
	    dtheta.push_back(-999.);
	  } 
	  else if (!breakpointsgood[p] || !breakpointsgood[p-1]) {
	    dtheta.push_back(-999.);
	  }
	  else {
	  */ 
	    TVector3 rotated_pcdir0 = RotateToZaxis(pcdir0, pcdir0);
	    TVector3 rotated_pcdir1 = RotateToZaxis(pcdir0, pcdir1);
	    double theta_yz = 1000. * TMath::ATan(rotated_pcdir1.Y() / rotated_pcdir1.Z());
	    double theta_xz = 1000. * TMath::ATan(rotated_pcdir1.X() / rotated_pcdir1.Z());
	    double theta_3D = 1000. * rotated_pcdir1.Theta();

	    const double cosval = pcdir0.X()*pcdir1.X()+pcdir0.Y()*pcdir1.Y()+pcdir0.Z()*pcdir1.Z();
	    double dt = 1000.*acos(cosval);
	    theta_yz_vec.push_back(theta_yz);
	    theta_xz_vec.push_back(theta_xz);
	    dtheta.push_back(dt);
	    double this_rr = rr2[breakpoints[p-1]]; // -- momentuem before the track enter the segment
	    //double this_rr = rr2[breakpoints[p]]; // -- momentum after the track travelled the segment
	    double this_KE_CSDA = muon_sp_range_to_KE -> Eval(this_rr);
	    double this_P_CSDA = muon_ek2p(this_KE_CSDA);
	    pmuon_vec.push_back(this_P_CSDA);
	}
	pcdir0 = pcdir1;
      }
      //    FillHist("dEdx_MPV_vs_corr_dqdx_" + suffix, this_dEdx_MPV, corrected_dqdx, 1., 3000., 0., 30., 3000., 0., 3000.);

      for(unsigned int i = 0; i < pmuon_vec.size(); i++){
	double theta_yz = theta_yz_vec.at(i);
	double theta_xz = theta_xz_vec.at(i);
	double theta_3d = dtheta.at(i);
	double P_muon = pmuon_vec.at(i); 
	FillHist("muon_p_vs_theta_xz", P_muon, theta_xz, 1., 3000., 0., 3000., 1000., -500., 500.);
	FillHist("muon_p_vs_theta_yz", P_muon, theta_yz, 1., 3000., 0., 3000., 1000., -500., 500.);
	FillHist("muon_p_vs_theta_2d", P_muon, theta_xz, 1., 3000., 0., 3000., 1000., -500., 500.);
	FillHist("muon_p_vs_theta_2d", P_muon, theta_yz, 1., 3000., 0., 3000., 1000., -500., 500.);
	FillHist("muon_p_vs_theta_3d", P_muon, theta_3d, 1., 3000., 0., 3000., 1000., 0., 1000.);
      }

    }
  }


  TString output_rootfile_dir = getenv("OUTPUTROOT_PATH");
  TString output_file_name = output_rootfile_dir + "/output_mcs_loop_run_" + run_number_str + ".root";
  if(!isdata) output_file_name = output_rootfile_dir + "/output_mcs_loop_mc.root";
  out_rootfile = new TFile(output_file_name, "RECREATE");
  out_rootfile -> cd();
  WriteHist();

  out_rootfile -> Close();
}
