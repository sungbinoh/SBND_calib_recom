#ifndef BETHEBLOCH_H
#define BETHEBLOCH_H

#include <map>
#include "Math/VavilovAccurate.h"
#include "TF1.h"
#include "TSpline.h"
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

ROOT::Math::VavilovAccurate vav;

// == Bethe-Bloch parameters, https://indico.fnal.gov/event/14933/contributions/28526/attachments/17961/22583/Final_SIST_Paper.pdf
const double rho = 1.39; // [g/cm3], density of LAr
const double K = 0.307075; // [MeV cm2 / mol]
const double Z = 18.; // atomic number of Ar
const double A = 39.948; // [g / mol], atomic mass of Ar
const double I = 188.0e-6; // [MeV], mean excitation energy
const double me = 0.511; // [Mev], mass of electron
// == Parameters for the density correction
const double density_C = 5.2146;
const double density_y0 = 0.2;
const double density_y1 = 3.0;
const double density_a = 0.19559;
const double density_k = 3.0;

double densityEffect(double beta, double gamma, double mass){

  double density_y = TMath::Log10(beta * gamma);
  double ln10 = TMath::Log(10);
  double this_delta = 0.;
  if(density_y > density_y1){
    this_delta = 2.0 * ln10 * density_y - density_C;
  }
  else if (density_y < density_y0){
    this_delta = 0.;
  }
  else{
    this_delta = 2.0 * ln10 * density_y - density_C + density_a * pow(density_y1 - density_y, density_k);
  }

  return this_delta;
}

double betaGamma(double KE, double mass){

  double gamma, beta;
  gamma = (KE + mass) / mass;
  beta = sqrt( 1 - 1/pow(gamma,2));

  return beta*gamma;

}

double Landau_xi(double KE, double pitch, double mass){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double xi = rho * pitch * 0.5 * K * (Z / A) * pow(1. / beta, 2);
  return xi;
}

double Get_Wmax(double KE, double mass){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * me * pow(beta * gamma, 2)) / (1.0 + 2.0 * me * (gamma / mass) + pow((me / mass),2));

  return Wmax;
}

double meandEdx(double KE, double mass){

  double gamma = (KE + mass) / mass;
  double beta = sqrt( 1 - 1/pow(gamma,2));
  double wmax = Get_Wmax(KE, mass);
  double dEdX = (rho*K*Z*pow(1.,2))/(A*pow(beta,2))*(0.5*log(2*me*pow(gamma,2)*pow(beta,2)*wmax/pow(I,2)) - pow(beta,2) - densityEffect( beta, gamma, mass )/2 ); // charge = 1

  return dEdX;
}

double MPVdEdx(double KE, double pitch, double mass){

  double gamma = (KE + mass) / mass;
  double beta = sqrt( 1 - 1/pow(gamma,2));

  double xi = Landau_xi(KE, pitch, mass);

  double eloss_mpv = xi*(log( 2*me*pow(gamma,2)*pow(beta,2) / I ) + log( xi / I ) + 0.2 - pow(beta,2) - densityEffect( beta, gamma, mass ) )/pitch;

  return eloss_mpv;
}

double IntegratedEdx(double mass, double KE0, double KE1, int n = 10000){

  if (KE0>KE1) swap(KE0, KE1);

  double step = (KE1-KE0)/n;

  double area = 0;

  for (int i = 0; i<n; ++i){
    double dEdx = meandEdx(KE0 + (i+0.5)*step, mass);
    if (dEdx)
      area += 1/dEdx*step;
  }
  return area;
}

double RangeFromKE(double KE, double mass){

  return IntegratedEdx(mass, 0, KE);
}

TSpline3 * Get_sp_range_KE(double mass, int np = 1000, double minke = .01, double maxke = 2e5){

  double *KE = new double[np];
  double *Range = new double[np];

  for (int i = 0; i<np; ++i){
    double ke = pow(10, log10(minke)+i*log10(maxke/minke)/np);
    KE[i] = ke;
    Range[i] = RangeFromKE(ke, mass);
  }

  TSpline3 *out = new TSpline3("sp_range_KE", Range, KE, np, "b2e2", 0, 0);
  delete[] KE;
  delete[] Range;

  return out;
}

double dEdx_PDF_setting(double *x, double *par){

  // == par[5] = {kappa, beta^2, xi, <dE/dx>BB, width}
  double a = par[2] / par[4];
  double b = (0.422784 + par[1] + log(par[0])) * par[2] / par[4] + par[3];
  double y = (x[0] - b) / a;

  double this_vav = 0.;

  if(par[0] < 0.01){ // == Landau
    this_vav = TMath::Landau(y);
    this_vav =  this_vav / a;
  }
  else if(par[0] > 10.){ // == Gaussian
    double mu = vav.Mean(par[0], par[1]);
    double sigma = sqrt(vav.Variance(par[0], par[1]));
    this_vav =  TMath::Gaus(y, mu, sigma);
  }
  else{ // == Vavilov
    this_vav =  vav.Pdf(y, par[0], par[1]);
    this_vav =  this_vav / a;
  }

  return this_vav;
}

TF1 *dEdx_PDF(double *par){

  TF1 *out = new TF1("", dEdx_PDF_setting, 0., 1000., 5);
  out -> SetParameter(0, par[0]);
  out -> SetParameter(1, par[1]);
  out -> SetParameter(2, par[2]);
  out -> SetParameter(3, par[3]);
  out -> SetParameter(4, par[4]);
  out -> SetParameters(par[0], par[1], par[2], par[3], par[4]);

  return out;
}

#endif
