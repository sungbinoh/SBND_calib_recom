#ifndef recom_h
#define recom_h

#include "mylib.h"

using namespace std;

class recom {

 public:
  recom();
  ~recom();

  // == Mod box functions 
  double dedx2dqdx_modbox(double dedx, double Ef = 0.5);
  double dqdx2dedx_modbox(double dqdx, double Ef = 0.5);
  double dedx2recomfactor_modbox(double dedx, double Ef = 0.5);

  // == EMB functions 
  double B_phi(double Ef = 0.5, double phi = 90.); // phi in deg
  double dedx2dqdx_emb(double dedx, double Ef = 0.5, double phi = 90.);
  double dedx2recomfactor_emb(double dedx, double Ef = 0.5, double phi = 90.);

 private:
  // == LAr properties
  double Rho = 1.39;
  double Wion = 23.6e-6;

  // == Mod. box par
  double alpha = 0.93;
  double betap = 0.212;

  // == EMB par
  double alpha_emb = 0.904;
  double beta_90 = 0.204;
  double R_emb = 1.25;

};

recom::recom(){}

// == Mod box functions
double recom::dedx2dqdx_modbox(double dedx, double Ef = 0.5){
  return ((Rho*Ef)/(betap*Wion)) * log(alpha + (betap*dedx) / (Rho*Ef));
}

double recom::dqdx2dedx_modbox(double dqdx, double Ef = 0.5){
  return (exp(dqdx*(betap/(Rho*Ef)*Wion))-alpha)/(betap/(Rho*Ef));
}

double recom::dedx2recomfactor_modbox(double dedx, double Ef = 0.5){
  double this_dqdx = dedx2dqdx_modbox(dedx, Ef);
  return (this_dqdx * Wion) / dedx;
}

// == EMB functions
double recom::B_phi(double Ef = 0.5, double phi = 90.){ // phi in deg
  double phi_rad = phi * TMath::Pi() / 180.;
  double out = beta_90 / (sqrt(pow(sin(phi_rad), 2.) + pow(cos(phi_rad) / R_emb, 2.)));
  return out;
}

double recom::dedx2dqdx_emb(double dedx, double Ef = 0.5, double phi = 90.){
  double this_B_phi = B_phi(Ef, phi);
  return ((Rho*Ef)/(this_B_phi*Wion)) * log(alpha_emb + (this_B_phi*dedx) / (Rho*Ef));
}

double recom::dedx2recomfactor_emb(double dedx, double Ef = 0.5, double phi = 90.){
  double this_dqdx = dedx2dqdx_emb(dedx, Ef, phi);
  return (this_dqdx * Wion) / dedx;
}

recom::~recom() {}

#endif
