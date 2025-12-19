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

   double cluster_energy_sum = 0; //adunam clusterele adiacente
   double weighted_cluster_energy_sum = 0;
   std::vector<double> centroid;

   for(int strip_index=0;strip_index<28;i++){
      
       if( cluster_energy_sum == 0 && strip_energies[i] ==0) continue;
       if( cluster_energy_sum != 0 && strip_energies[i] ==0){ // asta inseamna ca am ajuns la finalul unui cluster si salvez valorile calculate
            centroid.push_back(weighted_cluster_energy_sum/cluster_energy_sum);
            cluster_energy_sum = 0; 
            weighted_cluster_energy_sum = 0;
            continue;
       }

       weighted_cluster_energy_sum += strip_index * strip_energies;
       cluster_energy_sum += strip_energies;
   }
   //c++ este un jeg de limbaj si vectorii sai nu au functie de max =)
   double centroid_max = *std::max_element(centroid.begin(), centroid.end());
   return ((centroid_max-14.5)*2.54) //convertim in mm 
}


void Centroid_Determination_CATS(){

ifstream fisier_coeficienti_CATS1X("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_Cats2X.txt");
ifstream fisier_coeficienti_CATS1Y("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_Cats2Y.txt");
TFile *inputRootFile = TFile::Open("/media/olivia/Partition1/CATS/r0421_000a.root", "READ");

double slope_CATS1X[28];
double slope_CATS1Y[28];
double pedestal_CATS1X[28];

double intercept_CATS1X[28];
double intercept_CATS1Y[28];
double pedestal_CATS1Y[28];

citire_calibrare(fisier_coeficienti_CATS1X, slope_CATS1X, intercept_CATS1X, pedestal_CATS1X);
citire_calibrare(fisier_coeficienti_CATS1Y, slope_CATS1Y, intercept_CATS1Y, pedestal_CATS1Y);

TH1F *hCATS1X = new TH1F("hCATS1X", "CATS1X calibrated;Channel;Energy", 28,0,28, 1000, 0, 7000);
TH2F *hCATS1Y = new TH2F("hCATS1Y", "CATS1Y calibrated;Channel;Energy", 28, 0, 28, 1000, 0.0, 7000);

Int_t   CATS1XVM;
Float_t CATS1XV[28];
UShort_t CATS1XVN[28];

Int_t   CATS1YVM;
Float_t CATS1YV[28];
UShort_t CATS1YVN[28];
//event multiplicity, energy depositied, strip number

TTree *catsTree = (TTree*)inputRootFile->Get("AD");
catsTree->SetBranchAddress("CATS1XVM", &CATS1XVM);
catsTree->SetBranchAddress("CATS1XV",  &CATS1XV);
catsTree->SetBranchAddress("CATS1XVN", CATS1XVN);

catsTree->SetBranchAddress("CATS1YVM", &CATS1YVM);
catsTree->SetBranchAddress("CATS1YV",  &CATS1YV);
catsTree->SetBranchAddress("CATS1YVN", &CATS1YVN);

Long64_t entries = catsTree->GetEntries();

for (Long64_t entry = 0; entry < entries; ++entry) {

  catsTree->GetEntry(entry);
    //trebuie sa facem x1 vs y1, calculam centroizii si dupa plotam

  double energy_CATS1X_byStrip[28];
  double energy_CATS1Y_byStrip[28];
  double centroid_Y =0, centroid_X=0;
  for (int strip_index = 0; strip_index < 28; ++strip_index){ energy_CATS1X_byStrip[strip_index] = 0.0; energy_CATS1Y_byStrip[strip_index] = 0.0;}

  for (int j = 0; j < CATS1XVM; ++j) {

    //calibram si salvam energiile per strip
    if (CATS1XVN[j] == 0 || CATS1XVN[j] > 27) { continue; }
    if(CATS1XV[j]< pedestal_CATS1X[CATS1XVN[j]]+100)continue; //threshold
    double calibrated_energy = (CATS1XV[j] - pedestal_CATS1X[CATS1XVN[j]]) * slope_CATS1X[CATS1XVN[j]] + intercept_CATS1X[CATS1XVN[j]];
    energy_CATS1X_byStrip[CATS1XVN[j]]+=calibrated_energy;
  }

  for (int j = 0; j < CATS1YVM; ++j) {

    if (CATS1YVN[j] == 0 || CATS1YVN[j] > 27) { continue; }
    if(CATS1YV[j]< pedestal_CATS1Y[CATS1YVN[j]]+100)continue;
    double calibrated_energy = (CATS1YV[j] - pedestal_CATS1Y[CATS1YVN[j]]) * slope_CATS1Y[CATS1YVN[j]] + intercept_CATS1Y[CATS1YVN[j]];
    energy_CATS1Y_byStrip[CATS1XVN[j]]+=calibrated_energy;
  }

  centroid_X = calculate_centroid(energy_CATS1X_byStrip);
  centroid_Y = calculate_centroid(energy_CATS1Y_byStrip);

}

TCanvas *cX = new TCanvas("cX","CATS1XV vs XVN",800,600);
hCATS1X->Draw("COLZ");
cX->Update();


}
