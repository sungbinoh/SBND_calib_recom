#ifndef SCECORR_H
#define SCECORR_H

#include "Math/Vector3D.h"
#include <iostream>
#include <vector>

using namespace std;
using ROOT::Math::XYZVector;

class SCECorr {

 public:
  SCECorr();
  SCECorr(bool isdata);
  ~SCECorr();

  void ReadHistograms();
  
  void SetIsData(bool isdata);
  bool GetIsData(){ return _isdata; };

  double meas_pitch(double spx,  double spy, double spz, double dirx, double diry, double dirz, int plane, bool apply_sce = true);
  XYZVector GetCalPosOffsets(const XYZVector& in);
  XYZVector WireToTrajectoryPosition(const XYZVector& in);

 private:
  bool _isdata;
  double _ind_pl_pitch = 0.3;
  double _col_pl_pitch = 0.3;
  double _theta_z_east_1st = TMath::Pi() / 3.;
  double _theta_z_east_2nd = -1. * TMath::Pi() / 3.;
  double _theta_z_west_1st = -1. * TMath::Pi() / 3.;
  double _theta_z_west_2nd = TMath::Pi() / 3.;
  double _theta_z_col = 0.;

  // == SCE Corr Maps
  // E TPC
  TH3F* hTrueFwdX_E;
  TH3F* hTrueFwdY_E;
  TH3F* hTrueFwdZ_E;
  TH3F* hTrueBkwdX_E;
  TH3F* hTrueBkwdY_E;
  TH3F* hTrueBkwdZ_E;
  TH3F* hTrueEFieldX_E;
  TH3F* hTrueEFieldY_E;
  TH3F* hTrueEFieldZ_E;
  // W TPC
  TH3F* hTrueFwdX_W;
  TH3F* hTrueFwdY_W;
  TH3F* hTrueFwdZ_W;
  TH3F* hTrueBkwdX_W;
  TH3F* hTrueBkwdY_W;
  TH3F* hTrueBkwdZ_W;
  TH3F* hTrueEFieldX_W;
  TH3F* hTrueEFieldY_W;
  TH3F* hTrueEFieldZ_W;
};


SCECorr::SCECorr() : _isdata(false) {}

SCECorr::SCECorr(bool isdata) : SCECorr() { SetIsData(isdata); }


void SCECorr::ReadHistograms(){

  cout << "[SCECorr::ReadHistograms] Reading SCE histograms" << endl;
  TString datapath = getenv("SCEMAP_PATH");
  TDirectory* origDir = gDirectory;

  // == Refering
  // ==== https://github.com/SBNSoftware/sbndcode/blob/develop/sbndcode/SpaceCharge/SpaceChargeSBND.cxx
  // ==== https://github.com/LArSoft/larreco/blob/develop/larreco/Calorimetry/GnocchiCalorimetry_module.cc
  // == SCE map
  TString sce_map_file_path = datapath + "/SCEoffsets_SBND_E500_dualmap_voxelTH3.root";
  TFile *infile = TFile::Open(sce_map_file_path);
  // E TPC
  hTrueFwdX_E = (TH3F*) infile->Get("TrueFwd_Displacement_X_E");
  hTrueFwdY_E = (TH3F*) infile->Get("TrueFwd_Displacement_Y_E");
  hTrueFwdZ_E = (TH3F*) infile->Get("TrueFwd_Displacement_Z_E");
  hTrueBkwdX_E = (TH3F*) infile->Get("TrueBkwd_Displacement_X_E");
  hTrueBkwdY_E = (TH3F*) infile->Get("TrueBkwd_Displacement_Y_E");
  hTrueBkwdZ_E = (TH3F*) infile->Get("TrueBkwd_Displacement_Z_E");
  hTrueEFieldX_E = (TH3F*) infile->Get("True_ElecField_X_E");
  hTrueEFieldY_E = (TH3F*) infile->Get("True_ElecField_Y_E");
  hTrueEFieldZ_E = (TH3F*) infile->Get("True_ElecField_Z_E");
  // W TPC
  hTrueFwdX_W = (TH3F*) infile->Get("TrueFwd_Displacement_X_W");
  hTrueFwdY_W = (TH3F*) infile->Get("TrueFwd_Displacement_Y_W");
  hTrueFwdZ_W = (TH3F*) infile->Get("TrueFwd_Displacement_Z_W");
  hTrueBkwdX_W = (TH3F*) infile->Get("TrueBkwd_Displacement_X_W");
  hTrueBkwdY_W = (TH3F*) infile->Get("TrueBkwd_Displacement_Y_W");
  hTrueBkwdZ_W = (TH3F*) infile->Get("TrueBkwd_Displacement_Z_W");
  hTrueEFieldX_W = (TH3F*) infile->Get("True_ElecField_X_W");
  hTrueEFieldY_W = (TH3F*) infile->Get("True_ElecField_Y_W");
  hTrueEFieldZ_W = (TH3F*) infile->Get("True_ElecField_Z_W");

  hTrueFwdX_E->SetDirectory(0);
  hTrueFwdY_E->SetDirectory(0);
  hTrueFwdZ_E->SetDirectory(0);
  hTrueBkwdX_E->SetDirectory(0);
  hTrueBkwdY_E->SetDirectory(0);
  hTrueBkwdZ_E->SetDirectory(0);
  hTrueEFieldX_E->SetDirectory(0);
  hTrueEFieldY_E->SetDirectory(0);
  hTrueEFieldZ_E->SetDirectory(0);
  hTrueFwdX_W->SetDirectory(0);
  hTrueFwdY_W->SetDirectory(0);
  hTrueFwdZ_W->SetDirectory(0);
  hTrueBkwdX_W->SetDirectory(0);
  hTrueBkwdY_W->SetDirectory(0);
  hTrueBkwdZ_W->SetDirectory(0);
  hTrueEFieldX_W->SetDirectory(0);
  hTrueEFieldY_W->SetDirectory(0);
  hTrueEFieldZ_W->SetDirectory(0);
}

void SCECorr::SetIsData(bool isdata) {
  _isdata = isdata;
}

XYZVector SCECorr::GetCalPosOffsets(const XYZVector& in){

  std::vector<double> theCalPosOffsets;
  double xx = in.X(), yy = in.Y(), zz = in.Z();

  if(xx < -199.999){ xx = -199.999;}
  else if(xx > 199.999){ xx = 199.999;}
  if(yy < -199.999){ yy = -199.999;}
  else if(yy > 199.999){ yy = 199.999;}
  if(zz < 0.001){ zz = 0.001;}
  else if(zz > 499.999){ zz = 499.999;}

  //correct for charge drifted across cathode
  if ((xx > -2.5) && (xx < 0.)) { xx = -2.5; }
  if ((xx < 2.5) && (xx > 0.)) { xx = 2.5; }
  double offset_x = 0., offset_y = 0., offset_z = 0.;
  if(xx < 0){
    offset_x = hTrueBkwdX_E -> Interpolate(xx,yy,zz);
    offset_y = hTrueBkwdY_E -> Interpolate(xx,yy,zz);
    offset_z = hTrueBkwdZ_E -> Interpolate(xx,yy,zz);
  }
  else{
    offset_x = hTrueBkwdX_W -> Interpolate(xx,yy,zz);
    offset_y = hTrueBkwdY_W -> Interpolate(xx,yy,zz);
    offset_z = hTrueBkwdZ_W -> Interpolate(xx,yy,zz);
  }
  theCalPosOffsets = {offset_x, offset_y, offset_z};
  
  return { theCalPosOffsets[0], theCalPosOffsets[1], theCalPosOffsets[2] };
}

XYZVector SCECorr::WireToTrajectoryPosition(const XYZVector& in){
  XYZVector out = in;
  XYZVector offset = GetCalPosOffsets(in);

  // == FieldDistortionCorrectionXSign: 1 in SBND pandoraCalo
  out.SetX(in.X() + offset.X());
  out.SetY(in.Y() + offset.Y());
  out.SetZ(in.Z() + offset.Z());

  return out;
}

double SCECorr::meas_pitch(double spx, double spy, double spz, double dirx, double diry, double dirz, int plane, bool apply_sce){
  // == Calibration ntuple is made of tracks without SCE correction
  // == With this function, pitch after applying correction is measured
  // == For a double-checking, measured pitch after setting apply_sce to false is compared with the pitch branch that showed exactly the same results
  double this_angleToVert = 0.;
  double this_pitch = _ind_pl_pitch;
  if(plane == 0){
    if(spx < 0) this_angleToVert = _theta_z_east_1st;
    else this_angleToVert = _theta_z_west_1st;
  }
  else if(plane == 1){
    if(spx < 0) this_angleToVert = _theta_z_east_2nd;
    else this_angleToVert = _theta_z_west_2nd;
  }
  else if(plane == 2){
    this_pitch = _col_pl_pitch;
  }

  double cosgamma = fabs(sin(this_angleToVert) * diry + cos(this_angleToVert) * dirz);
  double pitch = 0.;
  if(cosgamma) pitch = this_pitch / cosgamma; 
  
  if(apply_sce){
    XYZVector dir(dirx, diry, dirz);
    XYZVector loc_w(spx, spy, spz);
    XYZVector locw_pdx_traj = WireToTrajectoryPosition(loc_w + pitch * dir);
    XYZVector loc = WireToTrajectoryPosition(loc_w);
    pitch = (locw_pdx_traj - loc).R();
  }
  
  return pitch;
}

#endif
