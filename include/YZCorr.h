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
  double GetYZCorr(const XYZVector& in, int plane);

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


YZCorr::YZCorr() : _yz_unif_file_str("") {}

YZCorr::YZCorr(TString yz_unif_file_str) : YZCorr() { SetFileStr(yz_unif_file_str); }

YZCorr::~YZCorr() {}

void YZCorr::ReadHistograms(){

  cout << "[YZCorr::ReadHistograms] Reading YZ unif. corr. histograms" << endl;
  TString datapath = getenv("SBND_YZCORR_PATH");
  TString data_v = getenv("SBND_YZCORR_VERSION");
  datapath = datapath + data_v + "/";
  TDirectory* origDir = gDirectory;

  // == Refering
  TString yzcorr_file_path = datapath + "/" + _yz_unif_file_str;
  //TString yzcorr_file_path = datapath + "/";
  TFile *infile = TFile::Open(yzcorr_file_path);
  // E TPC
  hyzcorr_plane0_E = (TH2D*) infile->Get("CzyHist_0_0");
  hyzcorr_plane1_E = (TH2D*) infile->Get("CzyHist_1_0");
  hyzcorr_plane2_E = (TH2D*) infile->Get("CzyHist_2_0");
  
  // W TPC
  hyzcorr_plane0_W = (TH2D*) infile->Get("CzyHist_0_1");
  hyzcorr_plane1_W = (TH2D*) infile->Get("CzyHist_1_1");
  hyzcorr_plane2_W = (TH2D*) infile->Get("CzyHist_2_1");
}

void YZCorr::SetFileStr(TString yz_unif_file_str) {
  _yz_unif_file_str = yz_unif_file_str;
}

double YZCorr::GetYZCorr(const XYZVector& in, int plane){
  // == get yz corr
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

  if(out_corr < 1e-3) out_corr = 1.;
  return out_corr;
}

#endif
