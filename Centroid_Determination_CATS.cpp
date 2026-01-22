#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TH2F.h>
#include <TGraph.h>

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <array>
#include <iostream>
using namespace std;


void citire_calibrare(std::ifstream &file, double slope[28], double intercept[28],double pedestal[28])
{

    int strip;
    for (int i = 0; i < 28; ++i){
        file >>strip>> slope[i] >> intercept[i]>>pedestal[i];
        cout<<strip<<" "<<slope[i]<<" "<<intercept[i]<<" "<<pedestal[i]<<'\n';
    }
    cout<<'\n';
}

// double calculate_centroid(double strip_energies[28]){
//     double cluster_energy_sum = 0;
//     double weighted_cluster_energy_sum = 0;
//     int cluster_strip_count = 0;

//     std::vector<double> centroid;
//     std::vector<double> energy;

//     for(int strip_index = 0; strip_index < 28; strip_index++){

//         if(strip_energies[strip_index] == 0){
//             if(cluster_strip_count > 0){
//                 if(cluster_strip_count > 1){
//                     centroid.push_back(weighted_cluster_energy_sum / cluster_energy_sum);
//                     energy.push_back(cluster_energy_sum);
//                 }
//                 cluster_energy_sum = 0;
//                 weighted_cluster_energy_sum = 0;
//                 cluster_strip_count = 0;
//             }
//             continue;
//         }

//         weighted_cluster_energy_sum += strip_index * strip_energies[strip_index];
//         cluster_energy_sum += strip_energies[strip_index];
//         cluster_strip_count++;
//     }

//     if(cluster_strip_count > 1){
//         centroid.push_back(weighted_cluster_energy_sum / cluster_energy_sum);
//         energy.push_back(cluster_energy_sum);
//     }

    
//     if(centroid.empty()){ return -1000;}

//     int index_max_energy = std::max_element(energy.begin(), energy.end()) - energy.begin();
//     double centroid_max = centroid[index_max_energy];

//     return ((centroid_max - 13.5) * 2.54);
// }



std::vector<int> build_fstrip(const double strip_energies[28]){
    std::vector<int> fired;
    fired.reserve(28);
    for(int strip_index = 0; strip_index < 28; strip_index++){
        if(strip_energies[strip_index] > 0.0) fired.push_back(strip_index);
    }
    std::sort(fired.begin(), fired.end(),
              [&](int left_index, int right_index){
                  return strip_energies[left_index] > strip_energies[right_index];
              });
    return fired;
}

bool check_neighbours(std::vector<int>& fstrip, const double strip_energies[28]){
    if((int)fstrip.size() < 3) return false;

    int center_strip = fstrip[0];

    auto is_neighbour = [&](int strip_index){
        return (std::abs(strip_index - center_strip) == 1) && (strip_energies[strip_index] > 0.0);
    };

    if(!is_neighbour(fstrip[1])) return false;
    if(is_neighbour(fstrip[2])) return true;

    if((int)fstrip.size() > 3 && is_neighbour(fstrip[3])){
        int tmp = fstrip[2]; fstrip[2] = fstrip[3]; fstrip[3] = tmp;
        return true;
    }

    if((int)fstrip.size() > 4 && is_neighbour(fstrip[4])){
        int tmp = fstrip[2]; fstrip[2] = fstrip[4]; fstrip[4] = tmp;
        return true;
    }

    return false;
}

double weighted_average(const double strip_energies[28], const std::vector<int>& fstrip, int nstrips){
    if(nstrips <= 0) return -1000.0;
    if((int)fstrip.size() < nstrips) return -1000.0;

    double numerator = 0.0;
    double denominator = 0.0;

    for(int rank_index = 0; rank_index < nstrips; rank_index++){
        int strip_index = fstrip[rank_index];
        double energy = strip_energies[strip_index];
        numerator += (double)strip_index * energy;
        denominator += energy;
    }

    if(denominator <= 0.0) return -1000.0;
    return numerator / denominator;
}

double sechip(const double strip_energies[28], const std::vector<int>& fstrip){
    if((int)fstrip.size() < 3) return -1000.0;

    int s0 = fstrip[0];
    int S2 = fstrip[1];
    int s2 = fstrip[2];

    double q0 = strip_energies[s0];
    double q1 = strip_energies[S2];
    double q2 = strip_energies[s2];

    if(q0 <= 0.0 || q1 <= 0.0 || q2 <= 0.0) return -1000.0;

    double v0 = std::sqrt(q0 / q2);
    double v1 = std::sqrt(q0 / q1);
    double v2 = 0.5 * (v0 + v1);
    if(v2 < 1.0) return -1000.0;

    double v3 = std::log(v2 + std::sqrt(v2 * v2 - 1.0));
    if(v3 == 0.0) return -1000.0;

    double v4 = (v0 - v1) / (2.0 * std::sinh(v3));
    if(std::abs(v4) >= 1.0) return -1000.0;

    double v5 = 0.5 * std::log((1.0 + v4) / (1.0 - v4));

    double xs = (double)s0 - (double)(s0 - S2) * v5 / v3;
    if(xs <= 0.0 || xs >= 28.0) return -1000.0;

    return xs;
}

double fit_gauss_array(const double strip_energies[28]){
    static TH1D* h1 = new TH1D("hQCats_tmp","hQCats_tmp",29,-0.5,28.5);
    static TF1* f1 = new TF1("f1_tmp","gaus(0)",0,28);

    h1->Reset();

    int nmax = -1;
    double qmax = 0.0;

    for(int strip_index = 0; strip_index < 28; strip_index++){
        double q = strip_energies[strip_index];
        if(q <= 0.0) continue;

        h1->Fill((double)strip_index, q);

        if(q > qmax){
            qmax = q;
            nmax = strip_index;
        }
    }

    if(nmax < 0) return -1000.0;

    f1->SetParameter(0,1000);
    f1->SetParameter(1,14);
    f1->SetParameter(2,1.0);

    f1->SetParLimits(1,0,28);
    f1->SetParLimits(2,0,3);

    double fit_min = (double)nmax - 4.0;
    double fit_max = (double)nmax + 4.0;
    if(fit_min < 0.0) fit_min = 0.0;
    if(fit_max > 28.0) fit_max = 28.0;

    h1->Fit(f1,"QNLB","",fit_min,fit_max);

    double xtmp = f1->GetParameter(1);
    if(xtmp > 0.0 && xtmp < 28.0) return xtmp;
    return -1000.0;
}

double calculate_centroid(double strip_energies[28]){
    std::vector<int> fstrip = build_fstrip(strip_energies);
    int m = (int)fstrip.size();
    if(m < 1) return -1000.0;

    double x = -1000.0;

    if(m >= 3){
        if(check_neighbours(fstrip, strip_energies)){
            x = sechip(strip_energies, fstrip);
        } else {
            x = fit_gauss_array(strip_energies);
        }
    } else if(m == 2){
        x = weighted_average(strip_energies, fstrip, 2);
    } else {
        if(strip_energies[fstrip[0]] > 500.0) x = (double)fstrip[0];
        else x = -1000.0;
    }

    if(x == -1000.0) return -1000.0;
    return ((x - 13.5) * 2.54);
}

void Centroid_Determination_CATS(){

ifstream fisier_coeficienti_CATS2X("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS2X.txt");
ifstream fisier_coeficienti_CATS2Y("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS2Y.txt");
TFile *inputRootFile = TFile::Open("/media/olivia/Partition1/CATS/r0422_000a_good_order_nocal.root", "READ");

double slope_CATS2X[28];
double slope_CATS2Y[28];
double pedestal_CATS2X[28];

double intercept_CATS2X[28];
double intercept_CATS2Y[28];
double pedestal_CATS2Y[28];

citire_calibrare(fisier_coeficienti_CATS2X, slope_CATS2X, intercept_CATS2X, pedestal_CATS2X);
citire_calibrare(fisier_coeficienti_CATS2Y, slope_CATS2Y, intercept_CATS2Y, pedestal_CATS2Y);

TH2F *hCATS2XY = new TH2F("hCATS2XY", "CATSXV vs CATS YV", 1000,-40,40, 1000,-40,40);

Int_t   CATS2XVM;
Float_t CATS2XV[28];
UShort_t CATS2XVN[28];

Int_t   CATS2YVM;
Float_t CATS2YV[28];
UShort_t CATS2YVN[28];
//event multiplicity, energy depositied, strip number

Float_t Id_6;
Float_t Id_11;

TTree *catsTree = (TTree*)inputRootFile->Get("AD");

catsTree->SetBranchAddress("Id_6", &Id_6);
catsTree->SetBranchAddress("Id_11", &Id_11);

catsTree->SetBranchAddress("CATS2XVM", &CATS2XVM);
catsTree->SetBranchAddress("CATS2XV",  &CATS2XV);
catsTree->SetBranchAddress("CATS2XVN", CATS2XVN);

catsTree->SetBranchAddress("CATS2YVM", &CATS2YVM);
catsTree->SetBranchAddress("CATS2YV",  &CATS2YV);
catsTree->SetBranchAddress("CATS2YVN", &CATS2YVN);

catsTree->SetBranchAddress("Id_6", &Id_6);
catsTree->SetBranchAddress("Id_11", &Id_11);

Long64_t entries = catsTree->GetEntries();

//normalizam raportat la primul strip
double slope_1_X     = slope_CATS2X[1];
double intercept_1_X = intercept_CATS2X[1];

for(int i = 0; i < 28; i++){
    slope_CATS2X[i]     = slope_CATS2X[i] / slope_1_X;
    intercept_CATS2X[i] = (intercept_CATS2X[i] - intercept_1_X) / slope_1_X;
}

double slope_1_Y     = slope_CATS2Y[1];
double intercept_1_Y = intercept_CATS2Y[1];

for(int i = 0; i < 28; i++){
    slope_CATS2Y[i]     = slope_CATS2Y[i] / slope_1_Y;
    intercept_CATS2Y[i] = (intercept_CATS2Y[i] - intercept_1_Y) / slope_1_Y;
}

//std::ofstream outFile("CATS2_centroid.txt");

for (Long64_t entry = 0; entry < entries; ++entry) {

  catsTree->GetEntry(entry);
    //trebuie sa facem x1 vs y1, calculam centroizii si dupa plotam

  //if(!(Id_6>5 && Id_6<6.5 && Id_11>269 && Id_11<273))continue;

  double energy_CATS2X_byStrip[28];
  double energy_CATS2Y_byStrip[28];
  double centroid_Y =0, centroid_X=0;
  for (int strip_index = 0; strip_index < 28; ++strip_index){ energy_CATS2X_byStrip[strip_index] = 0.0; energy_CATS2Y_byStrip[strip_index] = 0.0;}
 
  for (int j = 0; j < CATS2XVM; ++j) {

    //calibram si salvam energiile per strip
    if (CATS2XVN[j] ==0 || CATS2XVN[j] > 27) { continue; }
    if(CATS2XV[j]< pedestal_CATS2X[CATS2XVN[j]]+100)continue; //threshold
    double calibrated_energy = (CATS2XV[j] - pedestal_CATS2X[CATS2XVN[j]]-100) * slope_CATS2X[CATS2XVN[j]] + intercept_CATS2X[CATS2XVN[j]];
    energy_CATS2X_byStrip[CATS2XVN[j]]+=calibrated_energy;
  }

  for (int j = 0; j < CATS2YVM; ++j) {

    if (CATS2YVN[j] ==0 || CATS2YVN[j] > 27) { continue; }
    if(CATS2YV[j]< pedestal_CATS2Y[CATS2YVN[j]]+100)continue;
    double calibrated_energy = (CATS2YV[j] - pedestal_CATS2Y[CATS2YVN[j]]) * slope_CATS2Y[CATS2YVN[j]] + intercept_CATS2Y[CATS2YVN[j]];
    energy_CATS2Y_byStrip[CATS2YVN[j]]+=calibrated_energy;
  }

  //construim histograma cu centroizi
  centroid_X = calculate_centroid(energy_CATS2X_byStrip);
  centroid_Y = calculate_centroid(energy_CATS2Y_byStrip);
  //luam doar evenimentele nenule   
  if(centroid_X==-1000 || centroid_Y==-1000) continue;
  hCATS2XY->Fill(centroid_X, centroid_Y);
  outFile << centroid_X << " " << centroid_Y << "\n";
 }

outFile.close();

TCanvas *cXY = new TCanvas("cX","CATS2XV vs XVN",800,600);
hCATS2XY->Draw("COLZ");
cXY->Update();

// TFile *f = new TFile("CATS2_Centroid.root","RECREATE");
// hCATS2XY->Write();
// cXY->Write();
// f->Close();

}
