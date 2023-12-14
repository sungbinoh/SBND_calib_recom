#ifndef VAVGAUSFIT_H
#define VAVGAUSFIT_H

#include "TH1.h"
#include "TF1.h"
#include "Math/VavilovAccurate.h"

ROOT::Math::VavilovAccurate vav;

// == Bethe-Bloch parameters, https://indico.fnal.gov/event/14933/contributions/28526/attachments/17961/22583/Final_SIST_Paper.pdf
const double rho = 1.39; // [g/cm3], densityof Lar
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

double densityEffect(double beta, double gamma){
  // == Estimate the density correction
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

double Get_Wmax(double KE, double mass){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * me * pow(beta * gamma, 2)) / (1.0 + 2.0 * me * (gamma / mass) + pow((me / mass),2));

  return Wmax;
}

double Landau_xi(double KE, double pitch, double mass){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double xi = rho * pitch * 0.5 * K * (Z / A) * pow(1. / beta, 2);
  return xi;
}

double meandEdx(double KE, double mass){

  double gamma = (KE + mass) / mass;
  double beta = sqrt( 1 - 1/pow(gamma,2));
  double wmax = Get_Wmax(KE, mass);
  double dEdX = (rho*K*Z*pow(1.,2))/(A*pow(beta,2))*(0.5*log(2*me*pow(gamma,2)*pow(beta,2)*wmax/pow(I,2)) - pow(beta,2) - densityEffect( beta, gamma )/2 );

  return dEdX;
}

double ResLength_to_KE_BB(double ResLength, double mass){
  // == KE to ResLength using Bethe-Bloch (BB)
  double KE_BB = 0.1; // [MeV], starting with non-zero energy
  double this_length = 0.;
  double step = 0.01; // [cm]
  bool first = true;
  while(this_length < ResLength){
    double this_dEdx = meandEdx(KE_BB, mass);
    KE_BB += this_dEdx * step;
    this_length += step;
  }
  return KE_BB;
}

double MPVdEdxfun(Double_t *x, Double_t *par){

  Double_t this_pitch = par[0];
  Double_t this_mass = par[1];
  Double_t this_range = x[0];
  Double_t this_KE = ResLength_to_KE_BB(this_range, this_mass);

  double this_gamma = (this_KE + this_mass) / this_mass;
  double this_beta = sqrt( 1 - 1 / pow(this_gamma,2));

  double this_xi = Landau_xi(this_KE, this_pitch, this_mass);

  double eloss_mpv = this_xi*(log( 2*me*pow(this_gamma,2)*pow(this_beta,2) / I ) + log( this_xi / I ) + 0.2 - pow(this_beta,2) - densityEffect( this_beta, this_gamma ) ) / this_pitch;

  return eloss_mpv;
}

Double_t vavgaufun(Double_t *x, Double_t *par) {
  
  Double_t this_pitch = par[0];
  Double_t this_KE = par[1];
  Double_t this_mass = par[2];
  Double_t this_Gaus_sigma = par[3];
  Double_t this_norm_factor = par[4];
  
  Double_t this_Wmax = Get_Wmax(this_KE, this_mass);
  Double_t this_Landau_xi = Landau_xi(this_KE, this_pitch, this_mass);

  Double_t this_kappa = this_Landau_xi / this_Wmax;

  Double_t this_gamma = (this_KE/this_mass)+1.0;
  Double_t this_beta_sqaure = 1.0 - 1.0/(this_gamma * this_gamma);

  Double_t this_a = this_Landau_xi / this_pitch;
  Double_t this_mean_dedx = meandEdx(this_KE, this_mass);
  Double_t this_b = (0.422784 + this_beta_sqaure + log(this_kappa)) * this_a + this_mean_dedx;
  //Double_t this_y = (x[0] - this_a) / this_b;

  Double_t np = 500.0;
  Double_t sc = 5.0;// convolution extends to +-sc Gaussian sigmas
  Double_t xx;
  Double_t yy;
  Double_t f_vav;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;

  xlow = x[0] - sc * this_Gaus_sigma;
  xupp = x[0] + sc * this_Gaus_sigma;
  step = (xupp-xlow)/np;

  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    yy = (xx - this_b) / this_a;
    if(this_kappa > 10.){
      double this_mu = vav.Mean(this_kappa, this_beta_sqaure);
      double this_sigma = sqrt(vav.Variance(this_kappa, this_beta_sqaure));
      f_vav = TMath::Gaus(yy, this_mu, this_sigma);
    }
    else if(this_kappa < 0.01){
      f_vav = TMath::Landau(yy);
      f_vav = f_vav / this_a;
    }
    else f_vav = vav.Pdf(yy, this_kappa, this_beta_sqaure) / this_a;
    sum += f_vav * TMath::Gaus(x[0],xx,this_Gaus_sigma);

    xx = xupp - (i-.5) * step;
    yy = (xx - this_b) / this_a;
    if(this_kappa > 10.){
      double this_mu = vav.Mean(this_kappa, this_beta_sqaure);
      double this_sigma = sqrt(vav.Variance(this_kappa, this_beta_sqaure));
      f_vav = TMath::Gaus(yy, this_mu, this_sigma);
    }
    else if(this_kappa < 0.01){
      f_vav = TMath::Landau(yy);
      f_vav = f_vav / this_a;
    }
    else f_vav = vav.Pdf(yy, this_kappa, this_beta_sqaure) / this_a;
    sum += f_vav * TMath::Gaus(x[0],xx,this_Gaus_sigma);
  }

  Double_t invsq2pi = 0.398942280401;
  return (this_norm_factor * step * sum * invsq2pi / this_Gaus_sigma);
}

TF1 *vavgaufit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF, Int_t *Status, TString FunName){

  Int_t i;
  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;

  TF1 *ffit = new TF1(FunName,vavgaufun,fitrange[0],fitrange[1],5);
  ffit->SetParameters(startvalues);
  cout << "[vavgaufit] startvalues[2] : " << startvalues[2] << endl;
  cout << "[vavgaufit] startvalues[4] : " << startvalues[4] << endl;
  ffit->SetParNames("pitch","KE","mass","GSigma","NormFactor");

  for (i=0; i<5; i++) {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }
  //ffit->FixParameter(0, 0.65);
  ffit->FixParameter(2, startvalues[2]);


  TFitResultPtr fitres = his->Fit(FunName,"RSN"); // fit within specified range, use ParLimits, do not plot 
  ffit->GetParameters(fitparams);    // obtain fit parameters
  for (i=0; i<5; i++) {
    fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
  }

  ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
  NDF[0] = ffit->GetNDF();           // obtain pdf
  Status[0] = fitres->CovMatrixStatus();

  return (ffit);              // return fit function
}

#endif
