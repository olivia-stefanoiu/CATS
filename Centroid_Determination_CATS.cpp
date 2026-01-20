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

double calculate_centroid(double strip_energies[28]){
    double cluster_energy_sum = 0;
    double weighted_cluster_energy_sum = 0;
    int cluster_strip_count = 0;

    std::vector<double> centroid;
    std::vector<double> energy;

    for(int strip_index = 0; strip_index < 28; strip_index++){

        if(strip_energies[strip_index] == 0){
            if(cluster_strip_count > 0){
                if(cluster_strip_count > 1){
                    centroid.push_back(weighted_cluster_energy_sum / cluster_energy_sum);
                    energy.push_back(cluster_energy_sum);
                }
                cluster_energy_sum = 0;
                weighted_cluster_energy_sum = 0;
                cluster_strip_count = 0;
            }
            continue;
        }

        weighted_cluster_energy_sum += strip_index * strip_energies[strip_index];
        cluster_energy_sum += strip_energies[strip_index];
        cluster_strip_count++;
    }

    if(cluster_strip_count > 1){
        centroid.push_back(weighted_cluster_energy_sum / cluster_energy_sum);
        energy.push_back(cluster_energy_sum);
    }

    if(centroid.empty()) return -1000;

    int index_max_energy = std::max_element(energy.begin(), energy.end()) - energy.begin();
    double centroid_max = centroid[index_max_energy];

    return ((centroid_max - 13.5) * 2.54);
}



void Centroid_Determination_CATS(){

ifstream fisier_coeficienti_CATS1X("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS1X.txt");
ifstream fisier_coeficienti_CATS1Y("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS1Y.txt");
TFile *inputRootFile = TFile::Open("/media/olivia/Partition1/CATS/r0421_000a_mycal_nth_gord.root", "READ");

double slope_CATS1X[28];
double slope_CATS1Y[28];
double pedestal_CATS1X[28];

double intercept_CATS1X[28];
double intercept_CATS1Y[28];
double pedestal_CATS1Y[28];

citire_calibrare(fisier_coeficienti_CATS1X, slope_CATS1X, intercept_CATS1X, pedestal_CATS1X);
citire_calibrare(fisier_coeficienti_CATS1Y, slope_CATS1Y, intercept_CATS1Y, pedestal_CATS1Y);

TH2F *hCATS1XY = new TH2F("hCATS1XY", "CATSXV vs CATS YV", 1000,-25,25, 1000,-25,25);

Int_t   CATS1XVM;
Float_t CATS1XV[28];
UShort_t CATS1XVN[28];

Int_t   CATS1YVM;
Float_t CATS1YV[28];
UShort_t CATS1YVN[28];
//event multiplicity, energy depositied, strip number

Float_t Id_6;
Float_t Id_11;

TTree *catsTree = (TTree*)inputRootFile->Get("AD");

catsTree->SetBranchAddress("Id_6", &Id_6);
catsTree->SetBranchAddress("Id_11", &Id_11);

catsTree->SetBranchAddress("CATS1XVM", &CATS1XVM);
catsTree->SetBranchAddress("CATS1XV",  &CATS1XV);
catsTree->SetBranchAddress("CATS1XVN", CATS1XVN);

catsTree->SetBranchAddress("CATS1YVM", &CATS1YVM);
catsTree->SetBranchAddress("CATS1YV",  &CATS1YV);
catsTree->SetBranchAddress("CATS1YVN", &CATS1YVN);

catsTree->SetBranchAddress("Id_6", &Id_6);
catsTree->SetBranchAddress("Id_11", &Id_11);

Long64_t entries = catsTree->GetEntries();

//normalizam raportat la primul strip
double slope_1_X     = slope_CATS1X[1];
double intercept_1_X = intercept_CATS1X[1];

for(int i = 0; i < 28; i++){
    slope_CATS1X[i]     = slope_CATS1X[i] / slope_1_X;
    intercept_CATS1X[i] = (intercept_CATS1X[i] - intercept_1_X) / slope_1_X;
}

double slope_1_Y     = slope_CATS1Y[1];
double intercept_1_Y = intercept_CATS1Y[1];

for(int i = 0; i < 28; i++){
    slope_CATS1Y[i]     = slope_CATS1Y[i] / slope_1_Y;
    intercept_CATS1Y[i] = (intercept_CATS1Y[i] - intercept_1_Y) / slope_1_Y;
}

std::ofstream outFile("CATS1_centroid.txt");
for (Long64_t entry = 0; entry < entries; ++entry) {

  catsTree->GetEntry(entry);
    //trebuie sa facem x1 vs y1, calculam centroizii si dupa plotam

  if(!(Id_6>5 && Id_6<6.5 && Id_11>269 && Id_11<273))continue;

  double energy_CATS1X_byStrip[28];
  double energy_CATS1Y_byStrip[28];
  double centroid_Y =0, centroid_X=0;
  for (int strip_index = 0; strip_index < 28; ++strip_index){ energy_CATS1X_byStrip[strip_index] = 0.0; energy_CATS1Y_byStrip[strip_index] = 0.0;}
 
  for (int j = 0; j < CATS1XVM; ++j) {

    //calibram si salvam energiile per strip
    if (CATS1XVN[j] ==0 || CATS1XVN[j] > 27) { continue; }
    if(CATS1XV[j]< pedestal_CATS1X[CATS1XVN[j]]+100)continue; //threshold
    double calibrated_energy = (CATS1XV[j] - pedestal_CATS1X[CATS1XVN[j]]-100) * slope_CATS1X[CATS1XVN[j]] + intercept_CATS1X[CATS1XVN[j]];
    energy_CATS1X_byStrip[CATS1XVN[j]]+=calibrated_energy;
  }

  for (int j = 0; j < CATS1YVM; ++j) {

    if (CATS1YVN[j] ==0 || CATS1YVN[j] > 27) { continue; }
    if(CATS1YV[j]< pedestal_CATS1Y[CATS1YVN[j]]+100)continue;
    double calibrated_energy = (CATS1YV[j] - pedestal_CATS1Y[CATS1YVN[j]]) * slope_CATS1Y[CATS1YVN[j]] + intercept_CATS1Y[CATS1YVN[j]];
    energy_CATS1Y_byStrip[CATS1YVN[j]]+=calibrated_energy;
  }

  //construim histograma cu centroizi
  centroid_X = calculate_centroid(energy_CATS1X_byStrip);
  centroid_Y = calculate_centroid(energy_CATS1Y_byStrip);
  //luam doar evenimentele nenule   
  if(centroid_X==-1000 || centroid_Y==-1000) continue;
  hCATS1XY->Fill(centroid_X, centroid_Y);
  outFile << centroid_X << " " << centroid_Y << "\n";
 }

outFile.close();

TCanvas *cXY = new TCanvas("cX","CATS1XV vs XVN",800,600);
hCATS1XY->Draw("COLZ");
cXY->Update();

TFile *f = new TFile("CATS1_Centroid.root","RECREATE");
hCATS1XY->Write();
cXY->Write();
f->Close();

}
