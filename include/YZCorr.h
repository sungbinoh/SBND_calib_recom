#ifndef YZCORR_H
#define YZCORR_H

#include "Math/Vector3D.h"
#include <iostream>
#include <vector>

using namespace std;
using ROOT::Math::XYZVector;

class YZCorr {

 public:
  YZCorr();
  YZCorr(TString yz_unif_file_str);
  ~YZCorr();

  void ReadHistograms();
  
  void SetFileStr(TString yz_unif_file_str);
  bool GetIsData(){ return _isdata; };

  double meas_pitch(double spx,  double spy, double spz, double dirx, double diry, double dirz, int plane, bool apply_sce = true);
  XYZVector GetCalPosOffsets(const XYZVector& in);
  XYZVector WireToTrajectoryPosition(const XYZVector& in);
  XYZVector GetEfieldOffsets(const XYZVector& in);
  double GetEfield(const XYZVector& in);

 private:
  TString _yz_unif_file_str;

  // == YZ Corr Maps
  // E TPC
  TH2D* hyzcorr_plane0_E;
  TH2D* hyzcorr_plane1_E;
  TH2D* hyzcorr_plane2_E;
  
  // W TPC
  TH2D* hyzcorr_plane0_W;
  TH2D* hyzcorr_plane1_W;
  TH2D* hyzcorr_plane2_W;
};


YZCorr::YZCorr() : __yz_unif_file_str;("") {}

YZCorr::YZCorr(TString yz_unif_file_str) : YZCorr() { SetFileStr(yz_unif_file_str); }

void YZCorr::ReadHistograms(){

  cout << "[YZCorr::ReadHistograms] Reading SCE histograms" << endl;
  TString datapath = getenv("SBND_YZCORR_PATH");
  TString data_v = getenv("SBND_YZCORR_VERSION");
  datapath = datapath + data_v + "/";
  TDirectory* origDir = gDirectory;

  // == Refering
  TString yzcorr_file_path = datapath + "/" + _yz_unif_file_str;
  TFile *infile = TFile::Open(yzcorr_file_path);
  // E TPC
  hyzcorr_plane0_E = (TH2D*) infile->Get("yz_plane0_east_sce_corr");
  hyzcorr_plane1_E = (TH2D*) infile->Get("yz_plane1_east_sce_corr");
  hyzcorr_plane2_E = (TH2D*) infile->Get("yz_plane2_east_sce_corr");
  
  // W TPC
  hyzcorr_plane0_W = (TH2D*) infile->Get("yz_plane0_west_sce_corr");
  hyzcorr_plane1_W = (TH2D*) infile->Get("yz_plane1_west_sce_corr");
  hyzcorr_plane2_W = (TH2D*) infile->Get("yz_plane2_west_sce_corr");
}

void YZCorr::SetFileStr(TString yz_unif_file_str) {
  _yz_unif_file_str = yz_unif_file_str;
}

double YZCorr::GetYZCorr(const XYZVector& in, int plane){
  // == sce spactial distortion corrected position should be input
  double xx = in.X(), yy = in.Y(), zz = in.Z();
  if(xx < -199.999){ xx = -199.999;}
  else if(xx > 199.999){ xx = 199.999;}
  if(yy < -199.999){ yy = -199.999;}
  else if(yy > 199.999){ yy = 199.999;}
  if(zz < 0.001){ zz = 0.001;}
  else if(zz > 499.999){ zz = 499.999;}

  double out_corr = 1.;
  if(plane == 0){
    if(xx < 0) out_corr = hyzcorr_plane0_E->GetBinContent(hyzcorr_plane0_E->FindBin(zz, yy));
    else out_corr = hyzcorr_plane0_W->GetBinContent(hyzcorr_plane0_W->FindBin(zz, yy));
  }
  else if(plane == 1){
    if(xx < 0) out_corr = hyzcorr_plane1_E->GetBinContent(hyzcorr_plane1_E->FindBin(zz, yy));
    else out_corr = hyzcorr_plane1_W->GetBinContent(hyzcorr_plane1_W->FindBin(zz, yy));
  }
  else if(plane == 2){
    if(xx < 0) out_corr = hyzcorr_plane2_E->GetBinContent(hyzcorr_plane2_E->FindBin(zz, yy));
    else out_corr = hyzcorr_plane2_W->GetBinContent(hyzcorr_plane2_W->FindBin(zz, yy));
  }

  if(out_corr < 0.8 || out_corr > 1.2) out_corr = 1.;
  return out_corr;
}

XYZVector YZCorr::GetCalPosOffsets(const XYZVector& in){
  // == Get SCE corrections to space points
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
  //cout << "[YZCorr::GetCalPosOffsets] xx: " << xx << ", yy: " << yy << ", zz: " << zz << endl;
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

XYZVector YZCorr::WireToTrajectoryPosition(const XYZVector& in){
  // == Apply SCE correction to space points
  XYZVector out = in;
  XYZVector offset = GetCalPosOffsets(in);

  // == FieldDistortionCorrectionXSign: 1 in SBND pandoraCalo
  out.SetX(in.X() + offset.X());
  out.SetY(in.Y() + offset.Y());
  out.SetZ(in.Z() + offset.Z());

  return out;
}

double YZCorr::meas_pitch(double spx, double spy, double spz, double dirx, double diry, double dirz, int plane, bool apply_sce){
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

XYZVector YZCorr::GetEfieldOffsets(const XYZVector& in){
  // == This is for getting offset in E-field for a space point. Therefore, the input space point should be SCE corrected one.
  std::vector<double> theEfieldOffsets;
  double xx=in.X(), yy=in.Y(), zz=in.Z();
  double offset_x=0., offset_y=0., offset_z=0.;

  if(xx<-199.999){xx=-199.999;}
  else if(xx>199.999){xx=199.999;}
  if(yy<-199.999){yy=-199.999;}
  else if(yy>199.999){yy=199.999;}
  if(zz<0.001){zz=0.001;}
  else if(zz>499.999){zz=499.999;}
  if(xx < 0){
    offset_x = hTrueEFieldX_E->Interpolate(xx, yy, zz);
    offset_y = hTrueEFieldY_E->Interpolate(xx, yy, zz);
    offset_z = hTrueEFieldZ_E->Interpolate(xx, yy, zz);
  }
  else{
    offset_x = hTrueEFieldX_W->Interpolate(xx, yy, zz);
    offset_y = hTrueEFieldY_W->Interpolate(xx, yy, zz);
    offset_z = hTrueEFieldZ_W->Interpolate(xx, yy, zz);
  }
  theEfieldOffsets = {offset_x, offset_y, offset_z};

  return { theEfieldOffsets[0], theEfieldOffsets[1], theEfieldOffsets[2] };
}

double YZCorr::GetEfield(const XYZVector& in){
  XYZVector EfieldOffsets = GetEfieldOffsets(in);
  XYZVector xunit{1., 0., 0.};
  EfieldOffsets = EfieldOffsets + xunit;
  EfieldOffsets = EfieldOffsets * 0.5;
  double EField = EfieldOffsets.R();

  return EField;
}


#endif
