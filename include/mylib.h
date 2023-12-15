#ifndef mylib_h
#define mylib_h

ofstream outfile;

using namespace std;
const double alpha = 1 - 0.6827;
// == Debugging Mode
bool debug = true;
bool blind_SR = false;
TString Simulator;
// == Call all needed maps
map<TString, TH1D*> maphist;
map<TString, TH1F*> mapTH1F;
map<TString, TH2D*> maphist2D;
map<TString, TGraph*> map_gr;
map<TString, TGraphErrors*> map_err_gr;
map<TString, TGraphAsymmErrors*> map_asym_gr;
map<TString, TFile*> mapfile;
map<TString, TCanvas*> mapcanvas;
map<TString, TPad*> mappad;
map<TString, THStack*> maphstack;
map<TString, TLegend*> maplegend;
map<TString, double> map_overflow;
map<TString, TLine*> mapline;
map<TString, TKey*> maphistcheck;
map<TString, TList*> maplist;
map<TString, std::vector<double> > map_bin_vector;
map<TString, std::vector<TString> > map_sample_names;
map<TString, std::vector<double> > map_syst_array;
map<TString, std::vector<double> > map_syst_table;

double v_drift = 156.267;
double mass_muon = 105.658; // [MeV]
double mass_pion = 139.57; // [MeV]
double mass_proton = 938.272; // [MeV]
const int N_pi_type = 10;
TString pi_type_str[N_pi_type] = {"Fake Data",
				  "#pi^{+}_{Inel.}",
				  "#pi^{+}_{Elas.}",
				  "#mu",
				  "Misid. : cosmic",
				  "Misid. : p",
				  "Misid. : #pi^{+}",
				  "Misid. : #mu",
				  "Misid. : e/#gamma",
				  "Misid. : other"};
const int N_p_type = 9;
TString p_type_str[N_p_type] = {"Fake Data",
				"PInel",
				"PElas",
				"misID:cosmic",
				"misID:p",
				"misID:pi",
				"misID:mu",
				"misID:e/#gamma",
				"misID:other"};


const int N_MC = 6;
TString MC_category[N_MC] = {"NC", "NumuCC", "External_NC", "External_NueCC", "External_NumuCC", "NueCC"};

const int N_smear_bit = 8;
TString smear_flags[N_smear_bit] = {"NONE", "P", "Theta", "P_Theta", "Phi", "P_Phi", "Phi_Theta", "All"};

// == TChain
void AddFilesToChain(TString fileListPath, TChain* chain) {
  ifstream in(fileListPath);

  string fileName;
  while(getline(in,fileName)){
    chain->Add(fileName.c_str());
  }
}

// == Lifetime Correction
double Lifetime_Correction(double x, double tau){
  
  double out = 1.;
  if(fabs(x) > 200.) return out;

  double this_tdrift = (200. - fabs(x)) / v_drift;
  out = 1. / exp(-1. * this_tdrift / tau);

  return out;
}

// == HL parameters
double HL_kappa_a = 0.022;
double HL_kappa_c = 9.078;
double HL_sigma_res = 0.001908;
double HL_epsilon = 0.038;
double MCS_Get_HL_Sigma(double segment_size, double P, double mass){
  double kappa = ( (HL_kappa_a / (P*P)) + HL_kappa_c );
  double one_over_pbeta = pow(P*P + mass * mass, 0.5) / (P*P);
  double root_term = pow(segment_size / 14., 0.5);

  double sigma_HL = kappa * one_over_pbeta * root_term * (1 + HL_epsilon * log(segment_size / 14.));
  double out = pow(sigma_HL * sigma_HL + HL_sigma_res * HL_sigma_res, 0.5);
  return out;
}

vector<double> vx, vy, vexl, vexh, veyl, veyh;

TH1D * GetHist(TString hname){

  TH1D *h = NULL;
  std::map<TString, TH1D*>::iterator mapit = maphist.find(hname);
  if(mapit != maphist.end()) return mapit-> second;

  return h;

}

void Rebin_with_overflow(TString histname, int N_bin, double binx[]){
  if(debug) cout << "[Rebin_with_overflow]" << endl;
  maphist[histname + "rebin"] = (TH1D*)maphist[histname] -> Rebin(N_bin - 1, histname + "rebin", binx);
  double last_bin = 0.;
  last_bin =  maphist[histname + "rebin"] -> GetBinContent(N_bin - 1) + maphist[histname + "rebin"] -> GetBinContent(N_bin);
  maphist[histname + "rebin"] -> SetBinContent(N_bin - 1, last_bin);
}

TH1D* Add_Under_and_Overflow(TH1D* hist){

  TH1D* out = (TH1D*)hist -> Clone();
  double this_underflow_content = hist -> GetBinContent(0);
  double this_underflow_err = hist -> GetBinError(0);
  double this_first_bin_content = hist -> GetBinContent(1);
  double this_first_bin_err = hist -> GetBinError(1);
  double new_first_bin_content = this_underflow_content + this_first_bin_content;
  double new_first_bin_err = sqrt(pow(this_underflow_err, 2.) + pow(this_first_bin_err, 2.));
  out -> SetBinContent(1, new_first_bin_content);
  out -> SetBinError(1, new_first_bin_err);

  int N_bins = hist -> GetNbinsX();
  double this_overflow_content = hist -> GetBinContent(N_bins + 1);
  double this_overflow_err = hist -> GetBinError(N_bins + 1);
  double this_last_bin_content =  hist -> GetBinContent(N_bins);
  double this_last_bin_err =  hist -> GetBinError(N_bins);
  double new_last_bin_content = this_overflow_content + this_last_bin_content;
  double new_last_bin_err = sqrt(pow(this_overflow_err, 2.) + pow(this_last_bin_err, 2.));
  out -> SetBinContent(N_bins, new_last_bin_content);
  out -> SetBinError(N_bins, new_last_bin_err);

  return out;
  
}

void change_to_pseudo_data(TString current_histname){
  
  int N_bin = maphist[current_histname] -> GetNbinsX();
  for(int i = 1; i <= N_bin; i++){
    int current_bin_int =  maphist[current_histname] -> GetBinContent(i);
    double current_error = pow(current_bin_int, 0.5);
    maphist[current_histname] -> SetBinContent(i, current_bin_int);
    maphist[current_histname] -> SetBinError(i, current_error);
  }
}

void Proper_error_data(TString nameofhistogram, int N_bin, double binx[]){
  
  map_asym_gr[nameofhistogram + "correct_error"] = new TGraphAsymmErrors(GetHist(nameofhistogram + "rebin"));
  
  for(int i = 0; i < N_bin; i++){
    int N = GetHist(nameofhistogram + "rebin") -> GetBinContent(i + 1);
    double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
    double U =  (N==0) ? ( ROOT::Math::gamma_quantile_c(alpha,N+1,1) ) : ( ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) );
    if( N!=0 ){
      map_asym_gr[nameofhistogram + "correct_error"] -> SetPointEYlow(i, (N-L) ); // / (binx[i + 1] - binx[i]) );
      map_asym_gr[nameofhistogram + "correct_error"] -> SetPointEXlow(i, 0);
      map_asym_gr[nameofhistogram + "correct_error"] -> SetPointEYhigh(i, (U-N) );// / (binx[i + 1] - binx[i]) );
      map_asym_gr[nameofhistogram + "correct_error"] -> SetPointEXhigh(i, 0);
      
      double current_x = -1., current_y = -1.;
      map_asym_gr[nameofhistogram + "correct_error"] -> GetPoint(i, current_x, current_y);
      double modified_y = -1.;
      modified_y = current_y; // ( binx[i + 1] - binx[i] );
      map_asym_gr[nameofhistogram + "correct_error"] -> SetPoint(i, current_x, modified_y);

      
      //if(debug) cout << "i : " << i << ", current_x : " << current_x << ", current_y : " << current_y << ", modified_y : " << modified_y << endl;

      veyl.push_back( (N-L)); // / (binx[i + 1] - binx[i]) );
      veyh.push_back( (U-N)); // / (binx[i + 1] - binx[i]));
    }
    else{
      double zerodata_err_low = 0.1 ;// / (binx[i + 1] - binx[i]);
      double zerodata_err_high = 1.8 ;// / (binx[i + 1] - binx[i]);

      veyl.push_back(zerodata_err_low);
      veyh.push_back(zerodata_err_high);
      
      double current_x = GetHist(nameofhistogram + "rebin") -> GetBinCenter(i + 1);
      double current_ex = GetHist(nameofhistogram + "rebin") -> GetBinWidth(i + 1) /2.;

      map_asym_gr[nameofhistogram + "correct_error"] -> SetPointEYlow(i, zerodata_err_low);
      map_asym_gr[nameofhistogram + "correct_error"] -> SetPointEXlow(i, 0.);
      map_asym_gr[nameofhistogram + "correct_error"] -> SetPointEYhigh(i, zerodata_err_high);
      map_asym_gr[nameofhistogram + "correct_error"] -> SetPointEXhigh(i, 0.);

      vx.push_back(current_x);
      vexl.push_back(current_ex);
      vexh.push_back(current_ex);
    }
  }//end for i on N_bin
    
}

TGraphAsymmErrors* hist_to_graph(TH1D* hist, bool YErrorZero=false){

  TH1::SetDefaultSumw2(true);

  int Nbins = hist->GetXaxis()->GetNbins();
  double x[Nbins], y[Nbins], xlow[Nbins], xup[Nbins], ylow[Nbins], yup[Nbins];
  TAxis *xaxis = hist->GetXaxis();
  for(Int_t i=0; i<Nbins; i++){
    x[i] = xaxis->GetBinCenter(i+1);
    y[i] = hist->GetBinContent(i+1);
    xlow[i] = xaxis->GetBinCenter(i+1)-xaxis->GetBinLowEdge(i+1);
    xup[i] = xaxis->GetBinUpEdge(i+1)-xaxis->GetBinCenter(i+1);
    ylow[i] = hist->GetBinError(i+1);
    yup[i] = hist->GetBinError(i+1);
    if(YErrorZero){
      ylow[i] = 0;
      yup[i] = 0;
    }
  }
  TGraphAsymmErrors *out = new TGraphAsymmErrors(Nbins, x, y, xlow, xup, ylow, yup);
  out->SetTitle("");
  return out;

}

TGraphAsymmErrors* hist_to_graph(TH1D* hist, int n_skip_x_left){

  TH1::SetDefaultSumw2(true);

  int Nbins = hist->GetXaxis()->GetNbins()-n_skip_x_left;
  double x[Nbins], y[Nbins], xlow[Nbins], xup[Nbins], ylow[Nbins], yup[Nbins];
  TAxis *xaxis = hist->GetXaxis();
  for(Int_t i=1; i<=Nbins; i++){
    x[i-1] = xaxis->GetBinCenter(i+n_skip_x_left);
    y[i-1] = hist->GetBinContent(i+n_skip_x_left);
    xlow[i-1] = xaxis->GetBinCenter(i+n_skip_x_left)-xaxis->GetBinLowEdge(i+n_skip_x_left);
    xup[i-1] = xaxis->GetBinUpEdge(i+n_skip_x_left)-xaxis->GetBinCenter(i+n_skip_x_left);
    ylow[i-1] = hist->GetBinError(i+n_skip_x_left);
    yup[i-1] = hist->GetBinError(i+n_skip_x_left);
  }
  TGraphAsymmErrors *out = new TGraphAsymmErrors(Nbins, x, y, xlow, xup, ylow, yup);
  out->SetTitle("");
  return out;

}

TGraphAsymmErrors* hist_to_graph(TH1D* hist, int n_skip_x_left, int n_x_shift, int i_x_shift){

  TH1::SetDefaultSumw2(true);

  int Nbins = hist->GetXaxis()->GetNbins()-n_skip_x_left;
  double x[Nbins], y[Nbins], xlow[Nbins], xup[Nbins], ylow[Nbins], yup[Nbins];
  TAxis *xaxis = hist->GetXaxis();
  for(Int_t i=1; i<=Nbins; i++){
    x[i-1] = xaxis->GetBinCenter(i+n_skip_x_left);
    y[i-1] = hist->GetBinContent(i+n_skip_x_left);
    xlow[i-1] = xaxis->GetBinCenter(i+n_skip_x_left)-xaxis->GetBinLowEdge(i+n_skip_x_left);
    xup[i-1] = xaxis->GetBinUpEdge(i+n_skip_x_left)-xaxis->GetBinCenter(i+n_skip_x_left);
    ylow[i-1] = hist->GetBinError(i+n_skip_x_left);
    yup[i-1] = hist->GetBinError(i+n_skip_x_left);

    double dx = (xaxis->GetBinUpEdge(i+n_skip_x_left)-xaxis->GetBinLowEdge(i+n_skip_x_left))/double(n_x_shift+1);
    x[i-1] = xaxis->GetBinLowEdge(i+n_skip_x_left) + double(i_x_shift+1) * dx;
    xlow[i-1] = double(i_x_shift+1) * dx;
    xup[i-1] = xaxis->GetBinUpEdge(i+n_skip_x_left)-x[i-1];
  }
  TGraphAsymmErrors *out = new TGraphAsymmErrors(Nbins, x, y, xlow, xup, ylow, yup);
  out->SetTitle("");
  return out;

}


TGraphAsymmErrors* GraphSubtract(TGraphAsymmErrors *a, TGraphAsymmErrors *b, bool Rel){

  //==== do a-b

  int NX = a->GetN();
  TGraphAsymmErrors* gr_out = (TGraphAsymmErrors*)a->Clone();

  for(int i=0; i<NX; i++){

    double a_x, a_y, b_x, b_y;

    a->GetPoint(i, a_x, a_y);
    b->GetPoint(i, b_x, b_y);

    if(Rel==true){
      gr_out->SetPoint( i, a_x, (a_y-b_y)/b_y );
    }
    else{
      gr_out->SetPoint( i, a_x, a_y-b_y );
    }

  }

  return gr_out;

}

void RemoveLargeError(TGraphAsymmErrors *a){

  int NX = a->GetN();

  for(int i=0; i<NX; i++){

    double x, y, yerr_low, yerr_high;

    a->GetPoint(i, x, y);
    yerr_low  = a->GetErrorYlow(i);
    yerr_high = a->GetErrorYhigh(i);

    if(y+yerr_high==1.){
      a->SetPointEYhigh( i, yerr_low );
    }

  }

}

void ScaleGraph(TGraphAsymmErrors *a, double c){

  int NX = a->GetN();

  for(int i=0; i<NX; i++){

    double x, y, yerr_low, yerr_high;

    a->GetPoint(i, x, y);
    yerr_low  = a->GetErrorYlow(i);
    yerr_high = a->GetErrorYhigh(i);

    a->SetPoint(i, x, c*y);
    a->SetPointEYlow(i, c*yerr_low);
    a->SetPointEYhigh(i, c*yerr_high);

  }

}

double GetMaximum(TH1D* hist){

  TAxis *xaxis = hist->GetXaxis();

  double maxval(-1.);
  for(int i=1; i<=xaxis->GetNbins(); i++){
    if( hist->GetBinContent(i) + hist->GetBinError(i) > maxval ){
      maxval = hist->GetBinContent(i) + hist->GetBinError(i);
    }
  }

  return maxval;

}

double GetMaximum(TGraphAsymmErrors *a){

  int NX = a->GetN();

  double maxval(-9999.);
  for(int i=0; i<NX; i++){

    double x, y, yerr_low, yerr_high;

    a->GetPoint(i, x, y);
    yerr_low  = a->GetErrorYlow(i);
    yerr_high = a->GetErrorYhigh(i);

    if( y+yerr_high > maxval ){
      maxval = y+yerr_high;
    }

  }

  return maxval;

}

double GetYieldSystematics(TH1D *hist){

  int n_syst = hist->GetXaxis()->GetNbins();
  int n_source = (n_syst-1)/2;

  //==== Bin index
  //==== i=1 : central
  //==== i=2 : _MuonEn_up
  //==== i=3 : _MuonEn_down
  //==== -> n_syst = 3
  //==== -> n_source = (n_syst-1)/2 = (3-1)/2 = 1

  double y_central = hist->GetBinContent(1);

  double sum_syst = 0.;
  for(int i=1; i<=n_source; i++){
    double yield_up = hist->GetBinContent(i*2);
    double yield_down = hist->GetBinContent(i*2+1);

    double syst_up = fabs(yield_up-y_central);
    double syst_down = fabs(yield_down-y_central);

    sum_syst += 0.5*(syst_up*syst_up+syst_down*syst_down);

  }

  return sqrt(sum_syst);

}

TDirectory *MakeTemporaryDirectory(){

  gROOT->cd();
  TDirectory* tempDir = 0;
  int counter = 0;
  while (not tempDir) {
    // First, let's find a directory name that doesn't exist yet:
    std::stringstream dirname;
    dirname << "HNCommonLeptonFakes_%i" << counter;
    if (gROOT->GetDirectory((dirname.str()).c_str())) {
      ++counter;
      continue;
    }
    // Let's try to make this directory:
    tempDir = gROOT->mkdir((dirname.str()).c_str());

  }

  return tempDir;

}

void AddPhantomZero(double a, TString align, int digit_int, int digit_frac){

  if(align=="r"){
    int number_maxdigit = 0;
    for(int i=10; i>=0; i--){
      if(a/pow(10,i)>=1.){
        number_maxdigit = i;
        break;
      }
    }
    for(int i=0; i<digit_int-(number_maxdigit+1); i++) cout << "\\phantom{0}";
    cout << std::fixed<<std::setprecision(digit_frac) << a;
  }

  else if(align=="l"){
    int target_total_digit = digit_int+digit_frac;
    int number_maxdigit = 0;
    for(int i=10; i>=0; i--){
      if(a/pow(10,i)>=1.){
        number_maxdigit = i;
        break;
      }
    }
    cout << std::fixed<<std::setprecision(digit_frac) << a;
    for(int i=0; i<target_total_digit-(number_maxdigit+1)-digit_frac; i++) cout << "\\phantom{0}";
  }

}

TString Get_KE_Range_Str(double KE, int KE_step){
  
  int N_KE_step = KE / KE_step;
  TString out = Form("%dto%d", KE_step * N_KE_step, KE_step * (N_KE_step + 1));
  return out;
}

TH1D* GetHist1D(TString histname){

  TH1D *h = NULL;
  std::map<TString, TH1D*>::iterator mapit = maphist.find(histname);
  if(mapit != maphist.end()) return mapit->second;

  return h;

}

void FillHist(TString histname, double value, double weight, int n_bin, double x_min, double x_max){

  TH1D *this_hist = GetHist1D(histname);
  if( !this_hist ){
    this_hist = new TH1D(histname, "", n_bin, x_min, x_max);
    this_hist->SetDirectory(NULL);
    maphist[histname] = this_hist;
  }

  this_hist->Fill(value, weight);

}

#endif
