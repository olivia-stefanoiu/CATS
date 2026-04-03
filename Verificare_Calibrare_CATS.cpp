#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TGraph.h>

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <array>
#include <iostream>
using namespace std;

void citire_fisier(std::ifstream &file, double slope[28], double intercept[28],double pedestal[28])
{
    int strip;
    for (int i = 0; i < 28; ++i){
        file >>strip>> slope[i] >> intercept[i]>>pedestal[i];
        cout<<strip<<" "<<slope[i]<<" "<<intercept[i]<<" "<<pedestal[i]<<'\n';
    }
    cout<<'\n';
}


void Verificare_Calibrare_CATS(){
//TODO si pentru catS1

double slope_CATS1X[28];
double slope_CATS1Y[28];
double pedestal_CATS1X[28];

double intercept_CATS1X[28];
double intercept_CATS1Y[28];
double pedestal_CATS1Y[28];

ifstream fisier_coeficienti_CATS1X("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_proliant_CATS1X.txt");
ifstream fisier_coeficienti_CATS1Y("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_proliant_CATS1Y.txt");
TFile *inputRootFile = TFile::Open("/media/olivia/Partition1/CATS/Remerged/r0193_000a.root", "READ");
TTree *catsTree = (TTree*)inputRootFile->Get("AD");

citire_fisier(fisier_coeficienti_CATS1X, slope_CATS1X, intercept_CATS1X, pedestal_CATS1X);
citire_fisier(fisier_coeficienti_CATS1Y, slope_CATS1Y, intercept_CATS1Y, pedestal_CATS1Y);

TH2F *hCATS1X = new TH2F("hCATS1X", "CATS1X calibrated;Channel;Energy", 29, 0, 29, 2000, 0, 10000);
TH2F *hCATS1Y = new TH2F("hCATS1Y", "CATS1Y calibrated;Channel;Energy", 29, 0, 29, 2000, 0, 10000);

Int_t   CATS1XVM;
Float_t CATS1XV[28];
UShort_t CATS1XVN[28];

Int_t   CATS1YVM;
Float_t CATS1YV[28];
UShort_t CATS1YVN[28];
//event multiplicity, energy depositied, strip number
catsTree->SetBranchAddress("CATS1XVM", &CATS1XVM);
catsTree->SetBranchAddress("CATS1XV",  &CATS1XV);
catsTree->SetBranchAddress("CATS1XVN", CATS1XVN);

catsTree->SetBranchAddress("CATS1YVM", &CATS1YVM);
catsTree->SetBranchAddress("CATS1YV",  &CATS1YV);
catsTree->SetBranchAddress("CATS1YVN", &CATS1YVN);

Long64_t entries = catsTree->GetEntries();
//TODO make this into a function
for (Long64_t ie = 0; ie < entries; ++ie) {

        catsTree->GetEntry(ie);
        //cout<<"MULT "<<CATS1XVM<<'\n';
        for (int j = 0; j <CATS1XVM; ++j){
           //  cout<<"STRIP"<<CATS1XVN[j]<<'\n';
            if(CATS1XVN[j] == 0 || CATS1XVN[j]>27) {continue;} //acest strip este defect
           
           // cout<<"slope "<<slope_CATS1X[CATS1XVN[j]]<<'\n';
           // cout<<"en"<<CATS1XV[j]<<'\n';
            double calibrated_energy = (CATS1XV[j]-pedestal_CATS1X[CATS1XVN[j]])*slope_CATS1X[CATS1XVN[j]]+intercept_CATS1X[CATS1XVN[j]];
            hCATS1X->Fill(CATS1XVN[j], calibrated_energy);

        }

    }   

for (Long64_t ie = 0; ie < entries; ++ie) {

        catsTree->GetEntry(ie);

        for (int j = 0; j <CATS1YVM; ++j){

            if(CATS1YVN[j]==0 || CATS1YVN[j]>27) {continue;} //acest strip este defect

            double calibrated_energy = (CATS1YV[j]-pedestal_CATS1Y[CATS1YVN[j]])*slope_CATS1Y[CATS1YVN[j]]+intercept_CATS1Y[CATS1YVN[j]];
            hCATS1Y->Fill(CATS1YVN[j],calibrated_energy);

        }

    }  

TCanvas *cX = new TCanvas("cX", "CATS1X calibrated", 1000, 700);
hCATS1X->Draw();
cX->Update();

TCanvas *cY = new TCanvas("cY", "CATS1Y calibrated", 1000, 700);
hCATS1Y->Draw();
cY->Update();

}