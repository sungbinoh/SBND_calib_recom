#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "mylib.h"
#include "SCECorr.h"
#include "recom.h"

SCECorr *sce_corr_mc = new SCECorr(false);
recom *recom_fns = new recom();
bool isdata = false;

//int nGroupedWires = 10;
int NBinsX = 100;
int NBinsT = 100;
int NBinsdQdx = 300;
double minX = -200.;
double halfX = 0.;
double maxX = 200.;
double mindQdx = 0.;
double maxdQdx = 3000.;
double minT = 0.;
double maxT = 1.3;

double zprime_60deg(double y, double z, int pm = 1){
  double cos_60deg = 0.5;
  double sig_60deg = sqrt(3.) * 0.5;
  double zprime = z * cos_60deg - (pm + 0.) * sig_60deg * y;
  return zprime;
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
  cos_vals[1] = cos_zx;
  cos_vals[2] = cos_plus_zprimex;
  cos_vals[3] = cos_minus_zprimex;
}

void fill_lifetime_hists(int nGroupedWires, int plane, const TTreeReaderArray<float> &sp_x, const TTreeReaderArray<float> &sp_y, const TTreeReaderArray<float> &sp_z,
			 const TTreeReaderArray<float> &dirx, const TTreeReaderArray<float> &diry, const TTreeReaderArray<float> &dirz,
			 const TTreeReaderArray<uint16_t> &wire, const TTreeReaderArray<float> &dqdx, const TTreeReaderArray<float> &time, float trk_t0){
  // Refering to https://github.com/annab101/ElectronLifetimeSBND/blob/main/Main/dQdx_hist.cpp
  auto it_start_x = std::find_if(sp_x.begin(), sp_x.end(), [](float f){return !std::isnan(f);});
  int start_index = std::distance(sp_x.begin(), it_start_x);

  if(start_index == sp_x.GetSize()){return;}
  auto it_end_x = std::find_if(it_start_x, sp_x.end(), [](float f){return std::isnan(f);}) - 1;
  int end_index = std::distance(sp_x.begin(), it_end_x);
  if(end_index < 0){return;}

  int minWire = wire[start_index];
  int maxWire = wire[end_index];

  if(minWire > maxWire){
    std::swap(minWire, maxWire);
  }

  double dQdx_sum=0;
  double x_sum=0;
  double t_sum=0;
  int count=0;

  double dQdx_sce_sum=0;
  double x_sce_sum=0;

  vector<double> dedx_assumes = {1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3};
  vector<TString> dedx_assume_str = {"1p7", "1p8", "1p9", "2p0", "2p1", "2p2", "2p3"};
  vector<double> mb_dQdx_sum_vec_dedx_assumes(dedx_assumes.size(), 0.0);
  vector<double> mb_x_sce_sum_vec_dedx_assumes(dedx_assumes.size(), 0.0);
  vector<double> emb_dQdx_sum_vec_dedx_assumes(dedx_assumes.size(), 0.0);
  vector<double> emb_x_sce_sum_vec_dedx_assumes(dedx_assumes.size(), 0.0);

  for(int i = start_index; i <= end_index; i++){

    dQdx_sum += dqdx[i];
    x_sum += sp_x[i];
    t_sum += time[i];

    XYZVector sp_sce_uncorr(sp_x[i], sp_y[i], sp_z[i]);
    XYZVector sp_sce_corr = sce_corr_mc -> WireToTrajectoryPosition(sp_sce_uncorr);
    double pitch_sce_uncorr = sce_corr_mc -> meas_pitch(sp_x[i], sp_y[i], sp_z[i], dirx[i], diry[i], dirz[i], plane, false);
    double pitch_sce_corr = sce_corr_mc -> meas_pitch(sp_x[i], sp_y[i], sp_z[i], dirx[i], diry[i], dirz[i], plane, true);
    double dqdx_sce_corr = dqdx[i] * pitch_sce_uncorr / pitch_sce_corr;

    dQdx_sce_sum += dqdx_sce_corr;
    x_sce_sum += sp_sce_corr.X();

    double sce_efield = sce_corr_mc -> GetEfield(sp_sce_corr);
    vector<double> mb_recom_factor_nom_efield;
    vector<double> mb_recom_factor_sce_efield;
    vector<double> emb_recom_factor_nom_efield;
    vector<double> emb_recom_factor_sce_efield;
    
    for(unsigned int j = 0; j < dedx_assumes.size(); j++){
      double this_dedx = dedx_assumes.at(j);
      double mb_recom_fac_nom_field = recom_fns -> dedx2recomfactor_modbox(this_dedx, 0.5);
      double mb_recom_fac_sce_field = recom_fns -> dedx2recomfactor_modbox(this_dedx, sce_efield);
      double mb_recom_fac_corr_dqdx = dqdx_sce_corr * mb_recom_fac_nom_field / mb_recom_fac_sce_field;
      mb_dQdx_sum_vec_dedx_assumes.at(j) = mb_dQdx_sum_vec_dedx_assumes.at(j) + mb_recom_fac_corr_dqdx;

      double emb_recom_fac_nom_field = recom_fns -> dedx2recomfactor_emb(this_dedx, 0.5);
      double emb_recom_fac_sce_field = recom_fns -> dedx2recomfactor_emb(this_dedx, sce_efield);
      double emb_recom_fac_corr_dqdx = dqdx_sce_corr * emb_recom_fac_nom_field / emb_recom_fac_sce_field;
      emb_dQdx_sum_vec_dedx_assumes.at(j) = emb_dQdx_sum_vec_dedx_assumes.at(j) + emb_recom_fac_corr_dqdx;
    }

    count += 1;

    //Fill hist if swap TPC
    if( i < end_index && (sp_x[i] * sp_x[i+1]) < 0.){
      FillHist(Form("h_dQdx_xDrift_%dwires_plane%d", nGroupedWires, plane), x_sum/count, dQdx_sum/count, 1., NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);
      FillHist(Form("h_dQdx_xDrift_%dwires_plane%d_sce_corr", nGroupedWires, plane), x_sce_sum/count, dQdx_sce_sum/count, 1., NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);

      for(unsigned int j = 0; j < dedx_assumes.size(); j++){
	FillHist(Form("h_dQdx_xDrift_%dwires_plane%d_sce_corr_e_distor_mb_dedx", nGroupedWires, plane) + dedx_assume_str.at(j), x_sce_sum/count, mb_dQdx_sum_vec_dedx_assumes.at(j) / count, 1., NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);
	FillHist(Form("h_dQdx_xDrift_%dwires_plane%d_sce_corr_e_distor_emb_dedx", nGroupedWires, plane) + dedx_assume_str.at(j), x_sce_sum/count, emb_dQdx_sum_vec_dedx_assumes.at(j) / count, 1., NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);
      }

      if(-200. < x_sum/count && x_sum/count < 0.){
	FillHist(Form("h_dQdx_tDriftE_%dwires_plane%d", nGroupedWires, plane), t_sum/(count*2000) - 0.2 - trk_t0/1000000, dQdx_sum/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
	FillHist(Form("h_dQdx_tDriftE_%dwires_plane%d_sce_corr", nGroupedWires, plane), t_sum/(count*2000) - 0.2 - trk_t0/1000000, dQdx_sce_sum/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);

	for(unsigned int j = 0; j < dedx_assumes.size(); j++){
	  FillHist(Form("h_dQdx_tDriftE_%dwires_plane%d_sce_corr_e_distor_mb_dedx", nGroupedWires, plane) + dedx_assume_str.at(j), t_sum/(count*2000) - 0.2 - trk_t0/1000000, mb_dQdx_sum_vec_dedx_assumes.at(j)/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
	  FillHist(Form("h_dQdx_tDriftE_%dwires_plane%d_sce_corr_e_distor_emb_dedx", nGroupedWires, plane) + dedx_assume_str.at(j), t_sum/(count*2000) - 0.2 - trk_t0/1000000, emb_dQdx_sum_vec_dedx_assumes.at(j)/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
	}
      }
      if(0. < x_sum/count && x_sum/count < 200.){
	FillHist(Form("h_dQdx_tDriftW_%dwires_plane%d", nGroupedWires, plane), t_sum/(count*2000) - 0.2 - trk_t0/1000000, dQdx_sum/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
	FillHist(Form("h_dQdx_tDriftW_%dwires_plane%d_sce_corr", nGroupedWires, plane), t_sum/(count*2000) - 0.2 - trk_t0/1000000, dQdx_sce_sum/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);

	for(unsigned int j = 0; j < dedx_assumes.size(); j++){
          FillHist(Form("h_dQdx_tDriftW_%dwires_plane%d_sce_corr_e_distor_mb_dedx", nGroupedWires, plane) + dedx_assume_str.at(j), t_sum/(count*2000) - 0.2 - trk_t0/1000000, mb_dQdx_sum_vec_dedx_assumes.at(j)/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
          FillHist(Form("h_dQdx_tDriftW_%dwires_plane%d_sce_corr_e_distor_emb_dedx", nGroupedWires, plane) + dedx_assume_str.at(j), t_sum/(count*2000) - 0.2 - trk_t0/1000000, emb_dQdx_sum_vec_dedx_assumes.at(j)/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
        }
      }				
      
      dQdx_sum = 0.;
      x_sum = 0.;
      t_sum = 0.;
      count = 0;

      dQdx_sce_sum = 0.;
      x_sce_sum = 0.;

      for(unsigned int j = 0; j < dedx_assumes.size(); j++){
	mb_dQdx_sum_vec_dedx_assumes.at(j) = 0.;
	emb_dQdx_sum_vec_dedx_assumes.at(j) = 0.;
      }
    }
    else if(i == end_index){

      FillHist(Form("h_dQdx_xDrift_%dwires_plane%d", nGroupedWires, plane), x_sum/count, dQdx_sum/count, 1., NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);
      FillHist(Form("h_dQdx_xDrift_%dwires_plane%d_sce_corr", nGroupedWires, plane), x_sce_sum/count, dQdx_sce_sum/count, 1., NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);

      for(unsigned int j = 0; j < dedx_assumes.size(); j++){
        FillHist(Form("h_dQdx_xDrift_%dwires_plane%d_sce_corr_e_distor_mb_dedx", nGroupedWires, plane) + dedx_assume_str.at(j), x_sce_sum/count, mb_dQdx_sum_vec_dedx_assumes.at(j) / count, 1., NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);
        FillHist(Form("h_dQdx_xDrift_%dwires_plane%d_sce_corr_e_distor_emb_dedx", nGroupedWires, plane) + dedx_assume_str.at(j), x_sce_sum/count, emb_dQdx_sum_vec_dedx_assumes.at(j) / count, 1., NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);
      }

      if(-200. < x_sum/count && x_sum/count < 0.){
	FillHist(Form("h_dQdx_tDriftE_%dwires_plane%d", nGroupedWires, plane), t_sum/(count*2000) - 0.2 - trk_t0/1000000, dQdx_sum/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
	FillHist(Form("h_dQdx_tDriftE_%dwires_plane%d_sce_corr", nGroupedWires, plane), t_sum/(count*2000) - 0.2 - trk_t0/1000000, dQdx_sce_sum/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);

	for(unsigned int j = 0; j < dedx_assumes.size(); j++){
          FillHist(Form("h_dQdx_tDriftE_%dwires_plane%d_sce_corr_e_distor_mb_dedx", nGroupedWires, plane) + dedx_assume_str.at(j), t_sum/(count*2000) - 0.2 - trk_t0/1000000, mb_dQdx_sum_vec_dedx_assumes.at(j)/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
          FillHist(Form("h_dQdx_tDriftE_%dwires_plane%d_sce_corr_e_distor_emb_dedx", nGroupedWires, plane) + dedx_assume_str.at(j), t_sum/(count*2000) - 0.2 - trk_t0/1000000, emb_dQdx_sum_vec_dedx_assumes.at(j)/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
        }
      }
      
      if(0. < x_sum/count && x_sum/count < 200.){
	FillHist(Form("h_dQdx_tDriftW_%dwires_plane%d", nGroupedWires, plane), t_sum/(count*2000) - 0.2 - trk_t0/1000000, dQdx_sum/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
        FillHist(Form("h_dQdx_tDriftW_%dwires_plane%d_sce_corr", nGroupedWires, plane), t_sum/(count*2000) - 0.2 - trk_t0/1000000, dQdx_sce_sum/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);

	for(unsigned int j = 0; j < dedx_assumes.size(); j++){
          FillHist(Form("h_dQdx_tDriftW_%dwires_plane%d_sce_corr_e_distor_mb_dedx", nGroupedWires, plane) + dedx_assume_str.at(j), t_sum/(count*2000) - 0.2 - trk_t0/1000000, mb_dQdx_sum_vec_dedx_assumes.at(j)/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
          FillHist(Form("h_dQdx_tDriftW_%dwires_plane%d_sce_corr_e_distor_emb_dedx", nGroupedWires, plane) + dedx_assume_str.at(j), t_sum/(count*2000) - 0.2 - trk_t0/1000000, emb_dQdx_sum_vec_dedx_assumes.at(j)/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
        }
      }
      
      dQdx_sum = 0.;
      x_sum = 0.;
      t_sum = 0.;
      count = 0;

      dQdx_sce_sum = 0.;
      x_sce_sum = 0.;

      for(unsigned int j = 0; j < dedx_assumes.size(); j++){
        mb_dQdx_sum_vec_dedx_assumes.at(j) = 0.;
        emb_dQdx_sum_vec_dedx_assumes.at(j) = 0.;
      }
    }
    else{
      
      if((wire[i] - minWire)/nGroupedWires != (wire[i+1] - minWire)/nGroupedWires){
	
	FillHist(Form("h_dQdx_xDrift_%dwires_plane%d", nGroupedWires, plane), x_sum/count, dQdx_sum/count, 1., NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);
	FillHist(Form("h_dQdx_xDrift_%dwires_plane%d_sce_corr", nGroupedWires, plane), x_sce_sum/count, dQdx_sce_sum/count, 1., NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);

	for(unsigned int j = 0; j < dedx_assumes.size(); j++){
	  FillHist(Form("h_dQdx_xDrift_%dwires_plane%d_sce_corr_e_distor_mb_dedx", nGroupedWires, plane) + dedx_assume_str.at(j), x_sce_sum/count, mb_dQdx_sum_vec_dedx_assumes.at(j) / count, 1., NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);
	  FillHist(Form("h_dQdx_xDrift_%dwires_plane%d_sce_corr_e_distor_emb_dedx", nGroupedWires, plane) + dedx_assume_str.at(j), x_sce_sum/count, emb_dQdx_sum_vec_dedx_assumes.at(j) / count, 1., NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);
	}

	if(-200. < x_sum/count && x_sum/count < 0.){
	  FillHist(Form("h_dQdx_tDriftE_%dwires_plane%d", nGroupedWires, plane), t_sum/(count*2000) - 0.2 - trk_t0/1000000, dQdx_sum/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
	  FillHist(Form("h_dQdx_tDriftE_%dwires_plane%d_sce_corr", nGroupedWires, plane), t_sum/(count*2000) - 0.2 - trk_t0/1000000, dQdx_sce_sum/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);

	  for(unsigned int j = 0; j < dedx_assumes.size(); j++){
	    FillHist(Form("h_dQdx_tDriftE_%dwires_plane%d_sce_corr_e_distor_mb_dedx", nGroupedWires, plane) + dedx_assume_str.at(j), t_sum/(count*2000) - 0.2 - trk_t0/1000000, mb_dQdx_sum_vec_dedx_assumes.at(j)/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
	    FillHist(Form("h_dQdx_tDriftE_%dwires_plane%d_sce_corr_e_distor_emb_dedx", nGroupedWires, plane) + dedx_assume_str.at(j), t_sum/(count*2000) - 0.2 - trk_t0/1000000, emb_dQdx_sum_vec_dedx_assumes.at(j)/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
	  }
	}
	
	if(0. < x_sum/count && x_sum/count < 200.){
	  FillHist(Form("h_dQdx_tDriftW_%dwires_plane%d", nGroupedWires, plane), t_sum/(count*2000) - 0.2 - trk_t0/1000000, dQdx_sum/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
          FillHist(Form("h_dQdx_tDriftW_%dwires_plane%d_sce_corr", nGroupedWires, plane), t_sum/(count*2000) - 0.2 - trk_t0/1000000, dQdx_sce_sum/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);

	  for(unsigned int j = 0; j < dedx_assumes.size(); j++){
	    FillHist(Form("h_dQdx_tDriftW_%dwires_plane%d_sce_corr_e_distor_mb_dedx", nGroupedWires, plane) + dedx_assume_str.at(j), t_sum/(count*2000) - 0.2 - trk_t0/1000000, mb_dQdx_sum_vec_dedx_assumes.at(j)/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
	    FillHist(Form("h_dQdx_tDriftW_%dwires_plane%d_sce_corr_e_distor_emb_dedx", nGroupedWires, plane) + dedx_assume_str.at(j), t_sum/(count*2000) - 0.2 - trk_t0/1000000, emb_dQdx_sum_vec_dedx_assumes.at(j)/count, 1., NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
	  }
	}
	
	dQdx_sum = 0.;
	x_sum = 0.;
	t_sum = 0.;
	count = 0;

	dQdx_sce_sum = 0.;
	x_sce_sum = 0.;

	for(unsigned int j = 0; j < dedx_assumes.size(); j++){
	  mb_dQdx_sum_vec_dedx_assumes.at(j) = 0.;
	  emb_dQdx_sum_vec_dedx_assumes.at(j) = 0.;
	}
      }
    }
  }
}

void run_lifetime_loop(TString list_file, TString out_suffix, bool IsData = false) {

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
  // Check if the file exits
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
  TTreeReaderValue<Float_t> trk_t0(myReader, "trk.t0");
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
  TTreeReaderArray<uint16_t> wire0(myReader, "trk.hits0.h.wire");
  TTreeReaderArray<uint16_t> wire1(myReader, "trk.hits1.h.wire");
  TTreeReaderArray<uint16_t> wire2(myReader, "trk.hits2.h.wire");

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

  int N_run = 200000;
  double track_length_cut = 60.;
  // Loop over all entries of the TTree
  while (myReader.Next()) {
    if(current_entry > N_run) break;

    if(current_entry%100 == 0){
      cout << current_entry << " / " << N_entries << endl;
    }
    current_entry++;

    hist_selected -> Fill(*selected);

    // == Tracks selected as Anode+Cathode crossing
    if (*selected == 1) {
      // == 1st ind plane
      if(evt_sel(sp_x0, sp_y0, sp_z0, rr0, dqdx0)){
	// -- cos vals: cosyz, coszx, coszx+, coszx-
	double cos_vals_0[4];
        get_cos_vals(sp_x0, sp_y0, sp_z0, rr0, cos_vals_0);
	fill_lifetime_hists(10, 0, sp_x0, sp_y0, sp_z0, dirx0, diry0, dirz0, wire0, dqdx0, time0, *trk_t0);
      }

      if(evt_sel(sp_x1, sp_y1, sp_z1, rr1, dqdx1)){
        // -- cos vals: cosyz, coszx, coszx+, coszx-
        double cos_vals_1[4];
        get_cos_vals(sp_x1, sp_y1, sp_z1, rr1, cos_vals_1);
        fill_lifetime_hists(10, 1, sp_x1, sp_y1, sp_z1, dirx1, diry1, dirz1, wire1, dqdx1, time1, *trk_t0);
      }

      if(evt_sel(sp_x2, sp_y2, sp_z2, rr2, dqdx2)){
        // -- cos vals: cosyz, coszx, coszx+, coszx-
        double cos_vals_2[4];
        get_cos_vals(sp_x2, sp_y2, sp_z2, rr2, cos_vals_2);
        fill_lifetime_hists(10, 2, sp_x2, sp_y2, sp_z2, dirx2, diry2, dirz2, wire2, dqdx2, time2, *trk_t0);
      }
      
      /*
      // == No angular cut
      for (unsigned i = 0; i < dqdx.GetSize(); i++) {
	if(rr[i] < 0.) continue;
	
	FillHist("x_vs_y_acx", first_x, first_y, 1., 5000., -250., 250., 5000., -250., 250.);
	FillHist("x_vs_y_acx", last_x, last_y, 1., 5000., -250., 250., 5000., -250., 250.);
	
	FillHist("x_vs_z_acx", first_x, first_z, 1., 5000., -250., 250., 5000., -50., 550.);
	FillHist("x_vs_z_acx", last_x, last_z, 1., 5000., -250., 250., 5000., -50., 550.);
	
	FillHist("x_vs_y_trk_len", sp_x[i], sp_y[i], 1., 500., -250., 250., 500., -250., 250.);
        FillHist("x_vs_z_trk_len", sp_x[i], sp_z[i], 1., 500., -250., 250., 600., -50., 550.);
	
	//cout << i << ", time : " << time[i] << ", sp_x : " << sp_x[i] << endl;
	double t_drift_time = time[i] * 5.0e-4 - t_min;
        double t_drift_space = (200. - fabs(sp_x[i])) / v_drift;
	FillHist("sp_x_vs_dqdx_trk_len", sp_x[i], dqdx[i], 1., 500., -250., 250., 3000., 0., 3000.);
	FillHist("t_t0_vs_dqdx_trk_len", time[i], dqdx[i], 1., 3000., 0., 3000., 3000., 0., 3000.);
	FillHist("t_drift_time_vs_dqdx_trk_len", t_drift_time, dqdx[i], 1., 150., 0., 1.5, 3000., 0., 3000.);
        FillHist("t_drift_space_vs_dqdx_trk_len", t_drift_space, dqdx[i], 1., 150., 0., 1.5, 3000., 0., 3000.);

	double this_lifetime_corr = Lifetime_Correction(sp_x[i], 12.3);
	double corrected_dqdx = dqdx[i] * this_lifetime_corr;
	//cout << i << ", x : " << sp_x[i] << ", this_lifetime_corr : " << this_lifetime_corr << ", dqdx[i] : " << dqdx[i] << ", corrected_dqdx : " << corrected_dqdx << endl;
	FillHist("sp_x_vs_corr_dqdx_trk_len", sp_x[i], corrected_dqdx, 1., 500., -250., 250., 3000., 0., 3000.);
	FillHist("t_drift_time_vs_dqdx_corr_trk_len", t_drift_time, corrected_dqdx, 1., 150., 0., 1.5, 3000., 0., 3000.);
	FillHist("t_drift_space_vs_dqdx_corr_trk_len", t_drift_space, corrected_dqdx, 1.,150., 0., 1.5, 3000., 0., 3000.);

	FillHist("cos_xy_vs_dqdx_trk_len", cos_xy, dqdx[i], 1., 220., -1.1, 1.1, 3000., 0., 3000.);
	FillHist("cos_yz_vs_dqdx_trk_len", cos_yz, dqdx[i], 1., 220., -1.1, 1.1, 3000., 0., 3000.);
	FillHist("cos_zx_vs_dqdx_trk_len", cos_zx, dqdx[i], 1., 220., -1.1, 1.1, 3000., 0., 3000.);
      }

      // == With angular and passing cathode cut
      if(fabs(cos_zx) > 0.75 || !passing_cathode) continue;

      for(unsigned i = 0; i < dqdx.GetSize(); i++) {
	if(rr[i] < 0.) continue;

	FillHist("x_vs_y_angle_passing_cathode", sp_x[i], sp_y[i], 1., 500., -250., 250., 500., -250., 250.);
	FillHist("x_vs_z_angle_passing_cathode", sp_x[i], sp_z[i], 1., 500., -250., 250., 600., -50., 550.);

	double this_lifetime_corr = Lifetime_Correction(sp_x[i], 12.3);
        double corrected_dqdx = dqdx[i] * this_lifetime_corr;

	double t_drift_time = time[i] * 5.0e-4 - t_min;
	double t_drift_space = (200. - fabs(sp_x[i])) / v_drift;
	FillHist("t_drift_time_vs_dqdx_angle_passing_cathode", t_drift_time, dqdx[i], 1., 150., 0., 1.5, 3000., 0., 3000.);
	FillHist("t_drift_space_vs_dqdx_angle_passing_cathode", t_drift_space, dqdx[i], 1., 150., 0., 1.5, 3000., 0., 3000.);
	FillHist("t_drift_time_vs_dqdx_corr_angle_passing_cathode", t_drift_time, corrected_dqdx, 1., 150., 0., 1.5, 3000., 0., 3000.);
        FillHist("t_drift_space_vs_dqdx_corr_angle_passing_cathode", t_drift_space, corrected_dqdx, 1., 150., 0., 1.5, 3000., 0., 3000.);
      }
      */
    }
  }

  TString output_rootfile_dir = getenv("OUTPUTROOT_PATH");
  TString output_file_name = output_rootfile_dir + "/output_lifetime_" + out_suffix + ".root";
  out_rootfile = new TFile(output_file_name, "RECREATE");
  out_rootfile -> cd();
  
  hist_selected -> Write();

  WriteHist();
  out_rootfile -> Close();
}
