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
//TODO si pentru cats2

double slope_CATS2X[28];
double slope_CATS2Y[28];
double pedestal_CATS2X[28];

double intercept_CATS2X[28];
double intercept_CATS2Y[28];
double pedestal_CATS2Y[28];

ifstream fisier_coeficienti_CATS2X("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS2X.txt");
ifstream fisier_coeficienti_CATS2Y("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS2Y.txt");
TFile *inputRootFile = TFile::Open("/media/olivia/Partition1/CATS/r0193_000a.root", "READ");
TTree *catsTree = (TTree*)inputRootFile->Get("AD");

citire_fisier(fisier_coeficienti_CATS2X, slope_CATS2X, intercept_CATS2X, pedestal_CATS2X);
citire_fisier(fisier_coeficienti_CATS2Y, slope_CATS2Y, intercept_CATS2Y, pedestal_CATS2Y);

TH2F *hCATS2X = new TH2F("hCATS2X", "CATS2X calibrated;Channel;Energy", 28, 0.5, 28.5, 120, 0.0, 120);
TH2F *hCATS2Y = new TH2F("hCATS2Y", "CATS2Y calibrated;Channel;Energy", 28, 0.5, 28.5, 120, 0.0, 120);

Int_t   CATS2XVM;
Float_t CATS2XV[28];
UShort_t CATS2XVN[28];

Int_t   CATS2YVM;
Float_t CATS2YV[28];
UShort_t CATS2YVN[28];
//event multiplicity, energy depositied, strip number
catsTree->SetBranchAddress("CATS2XVM", &CATS2XVM);
catsTree->SetBranchAddress("CATS2XV",  &CATS2XV);
catsTree->SetBranchAddress("CATS2XVN", CATS2XVN);

catsTree->SetBranchAddress("CATS2YVM", &CATS2YVM);
catsTree->SetBranchAddress("CATS2YV",  &CATS2YV);
catsTree->SetBranchAddress("CATS2YVN", &CATS2YVN);

Long64_t entries = catsTree->GetEntries();
//TODO make this into a function
for (Long64_t ie = 0; ie < entries; ++ie) {

        catsTree->GetEntry(ie);
        //cout<<"MULT "<<CATS1XVM<<'\n';
        for (int j = 0; j <CATS2XVM; ++j){
           //  cout<<"STRIP"<<CATS1XVN[j]<<'\n';
            if(CATS2XVN[j] == 0 || CATS2XVN[j]>27) {continue;} //acest strip este defect
           
           // cout<<"slope "<<slope_CATS1X[CATS1XVN[j]]<<'\n';
           // cout<<"en"<<CATS1XV[j]<<'\n';
            double calibrated_energy = (CATS2XV[j]-pedestal_CATS2X[CATS2XVN[j]])*slope_CATS2X[CATS2XVN[j]]+intercept_CATS2X[CATS2XVN[j]];
            hCATS2X->Fill(CATS2XVN[j], calibrated_energy);

        }

    }   

for (Long64_t ie = 0; ie < entries; ++ie) {

        catsTree->GetEntry(ie);

        for (int j = 0; j <CATS2YVM; ++j){

            if(CATS2YVN[j]==0 || CATS2YVN[j]>27) {continue;} //acest strip este defect

            double calibrated_energy = (CATS2YV[j]-pedestal_CATS2Y[CATS2YVN[j]])*slope_CATS2Y[CATS2YVN[j]]+intercept_CATS2Y[CATS2YVN[j]];
            hCATS2Y->Fill(CATS2YVN[j],calibrated_energy);

        }

    }  

TCanvas *cX = new TCanvas("cX", "CATS2X calibrated", 1000, 700);
hCATS2X->Draw();
cX->Update();

TCanvas *cY = new TCanvas("cY", "CATS2Y calibrated", 1000, 700);
hCATS2Y->Draw();
cY->Update();

}