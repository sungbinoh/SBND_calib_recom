#ifndef SCECORR_H
#define SCECORR_H

using namespace std;

class SCECorr {

 public:
  SCECorr();
  SCECorr(bool isdata);
  ~SCECorr();

  void SetIsData(bool isdata);
  bool GetIsData(){ return _isdata; };

  double meas_pitch(double spx, double dirx, double diry, double dirz, int plane);
  
 private:
  bool _isdata;
  double _ind_pl_pitch = 0.3;
  double _col_pl_pitch = 0.3;
  double _theta_z_east_1st = TMath::Pi() / 3.;
  double _theta_z_east_2nd = -1. * TMath::Pi() / 3.;
  double _theta_z_west_1st = -1. * TMath::Pi() / 3.;
  double _theta_z_west_2nd = TMath::Pi() / 3.;
  double _theta_z_col = 0.;
  
};


SCECorr::SCECorr() : _isdata(false) {}

SCECorr::SCECorr(bool isdata) : SCECorr() { SetIsData(isdata); }

void SCECorr::SetIsData(bool isdata) {
  _isdata = isdata;
}

double SCECorr::meas_pitch(double spx, double dirx, double diry, double dirz, int plane){

  double this_angleToVert = 0.;
  double this_pitch = _ind_pl_pitch;
  if(plane == 0){
    if(spx < 0) this_angleToVert = _theta_z_east_1st;
    else this_angleToVert = _theta_z_west_1st;
  }
  else if(plane == 1){
    if(spx < 0) this_angleToVert = _theta_z_east_2nd;
    else this_angleToVert =	_theta_z_west_2nd;
  }
  else if(plane == 2){
    this_pitch = _col_pl_pitch;
  }

  double cosgamma = fabs(sin(this_angleToVert) * diry + cos(this_angleToVert) * dirz);
  double pitch = 0.;
  if(cosgamma) pitch = this_pitch / cosgamma; 

  return pitch;
}

#endif
