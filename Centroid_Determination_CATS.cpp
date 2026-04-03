#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TGraph.h>
#include <TChain.h>
#include <TH1F.h>
#include <TSystem.h>

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <array>
#include <iostream>
#include <math.h>

using namespace std;


void citire_calibrare(std::ifstream &file, double slope[28], double intercept[28], double pedestal[28])
{
    int strip;
    for (int i = 0; i < 28; ++i){
        file >> strip >> slope[i] >> intercept[i] >> pedestal[i];
        cout << strip << " " << slope[i] << " " << intercept[i] << " " << pedestal[i] << '\n';
    }
    cout << '\n';
}


void NormalizeToFirstStrip(double* slopeArray, double* interceptArray, int arrayLength = 28, int referenceIndex = 1)
{
    double referenceSlope = slopeArray[referenceIndex];
    double referenceIntercept = interceptArray[referenceIndex];

    for(int index = 0; index < arrayLength; index++)
    {
        slopeArray[index] = slopeArray[index] / referenceSlope;
        interceptArray[index] = (interceptArray[index] - referenceIntercept) / referenceSlope;
    }
}


void FillEnergiesByStrip(const UShort_t* stripNumberArray,
                         const Float_t* rawValueArray,
                         int multiplicity,
                         const double* pedestalArray,
                         const double* slopeArray,
                         const double* interceptArray,
                         double* energyByStripArray,
                         double thresholdOffset = 0.0)
{
    for(int strip_index = 0; strip_index < 28; ++strip_index) energyByStripArray[strip_index] = 0.0;

    for(int hit_index = 0; hit_index < multiplicity; ++hit_index)
    {
        int stripIndex = (int)stripNumberArray[hit_index];

        if(stripIndex == 0 || stripIndex > 27) continue;
        if((double)rawValueArray[hit_index] < pedestalArray[stripIndex] + thresholdOffset) continue;

        double calibratedEnergy = ((double)rawValueArray[hit_index] - pedestalArray[stripIndex]) * slopeArray[stripIndex] + interceptArray[stripIndex];
        energyByStripArray[stripIndex] += calibratedEnergy;
    }
}


std::vector<int> build_ordered_strip_energies(const double strip_energies[28])
{
    std::vector<int> fired;
    fired.reserve(28);
    for(int strip_index = 1; strip_index < 27; strip_index++){
        if(strip_energies[strip_index] > 0.0) fired.push_back(strip_index);
    }
    std::sort(fired.begin(), fired.end(),
              [&](int left_index, int right_index){
                  return strip_energies[left_index] > strip_energies[right_index];
              });
    return fired;
}


double weighted_average(const double strip_energies[28], const std::vector<int>& fstrip)
{
    if(fstrip.size() < 1) return -1000.0;

    int max_strip_index = fstrip[0];
    if(max_strip_index < 1 || max_strip_index > 26) return -1000.0;

    int min_strip_index = max_strip_index - 1;
    int max_window_index = max_strip_index + 1;

    if(min_strip_index < 1) min_strip_index = 1;
    if(max_window_index > 26) max_window_index = 26;

    double weighted_sum = 0.0;
    double energy_sum = 0.0;

    for(int strip_index = min_strip_index; strip_index <= max_window_index; ++strip_index)
    {
        double energy = strip_energies[strip_index];
        if(energy <= 0.0) continue;

        weighted_sum += energy * (double)strip_index;
        energy_sum += energy;
    }

    if(energy_sum <= 0.0) return -1000.0;
    return weighted_sum / energy_sum;
}


double calculate_centroid(double strip_energies[28], const std::vector<int>& fstrip)
{
    int m = (int)fstrip.size();
    if(m <= 0) return -1000.0;

    double x = weighted_average(strip_energies, fstrip);
    if(x == -1000.0) return -1000.0;

    double random_value = ((double)rand() / (double)RAND_MAX) * 1.0 - 0.51;
    return ((x - 13.5) * 2.54) + random_value;
}


void Centroid_Determination_CATS()
{
    ifstream fisier_coeficienti_CATS1X("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS1X (copy).txt");
    ifstream fisier_coeficienti_CATS1Y("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS1Y (copy).txt");

    TFile *inputRootFile = TFile::Open("/media/olivia/Partition1/CATS/r0421_000a_mycal_nth_gord.root", "READ");

    double slope_CATS1X[28];
    double slope_CATS1Y[28];
    double intercept_CATS1X[28];
    double intercept_CATS1Y[28];
    double pedestal_CATS1X[28];
    double pedestal_CATS1Y[28];

    citire_calibrare(fisier_coeficienti_CATS1X, slope_CATS1X, intercept_CATS1X, pedestal_CATS1X);
    citire_calibrare(fisier_coeficienti_CATS1Y, slope_CATS1Y, intercept_CATS1Y, pedestal_CATS1Y);

    NormalizeToFirstStrip(slope_CATS1X, intercept_CATS1X);
    NormalizeToFirstStrip(slope_CATS1Y, intercept_CATS1Y);

    TH2F *hCATS1XY = new TH2F("hCATS1XY", "CATS1 centroid;X (mm);Y (mm)", 1000, -40, 40, 1000, -40, 40);

    Int_t    CATS1XVM;
    Float_t  CATS1XV[28];
    UShort_t CATS1XVN[28];

    Int_t    CATS1YVM;
    Float_t  CATS1YV[28];
    UShort_t CATS1YVN[28];

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
    catsTree->SetBranchAddress("CATS1YVN", CATS1YVN);

    Long64_t entries = catsTree->GetEntries();

    double energy_CATS1X_byStrip[28];
    double energy_CATS1Y_byStrip[28];

    for (Long64_t entry = 0; entry < entries; ++entry)
    {
        catsTree->GetEntry(entry);

        FillEnergiesByStrip(CATS1XVN, CATS1XV, CATS1XVM, pedestal_CATS1X, slope_CATS1X, intercept_CATS1X,
                            energy_CATS1X_byStrip, 100.0);

        FillEnergiesByStrip(CATS1YVN, CATS1YV, CATS1YVM, pedestal_CATS1Y, slope_CATS1Y, intercept_CATS1Y,
                            energy_CATS1Y_byStrip, 0.0);

        std::vector<int> fstrip_CATS1X = build_ordered_strip_energies(energy_CATS1X_byStrip);
        std::vector<int> fstrip_CATS1Y = build_ordered_strip_energies(energy_CATS1Y_byStrip);

        double centroid_X = calculate_centroid(energy_CATS1X_byStrip, fstrip_CATS1X);
        double centroid_Y = calculate_centroid(energy_CATS1Y_byStrip, fstrip_CATS1Y);

        if(centroid_X == -1000.0 || centroid_Y == -1000.0) continue;

        hCATS1XY->Fill(centroid_X, centroid_Y);
    }

    TCanvas *cXY = new TCanvas("cXY","CATS1 centroid",800,600);
    hCATS1XY->Draw("COLZ");
    cXY->Update();

    gSystem->ProcessEvents();
}