#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
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


void citire_calibrare(std::ifstream &file, double slope[28], double intercept[28],double pedestal[28])
{

    int strip;
    for (int i = 0; i < 28; ++i){
        file >>strip>> slope[i] >> intercept[i]>>pedestal[i];
        cout<<strip<<" "<<slope[i]<<" "<<intercept[i]<<" "<<pedestal[i]<<'\n';
    }
    cout<<'\n';
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
                         double* energyByStripArray)
{
    for(int strip_index = 0; strip_index < 28; ++strip_index) energyByStripArray[strip_index] = 0.0;

    for(int hit_index = 0; hit_index < multiplicity; ++hit_index)
    {
        int stripIndex = (int)stripNumberArray[hit_index];

        if(stripIndex == 0 || stripIndex > 27) continue;
        if((double)rawValueArray[hit_index] < pedestalArray[stripIndex] + 100.0) continue;

        double calibratedEnergy = ((double)rawValueArray[hit_index] - pedestalArray[stripIndex] - 100.0) * slopeArray[stripIndex] + interceptArray[stripIndex];
        energyByStripArray[stripIndex] += calibratedEnergy;
    }
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


void rotate_and_shift_centroid = [](double input_x, double input_y,
                           double x_shift, double slope_shift, double intercept_shift, double center_position,
                           double& output_x, double& output_y)
{
    double rotation_angle = atan(slope_shift);
    double cos_angle = cos(-rotation_angle);
    double sin_angle = sin(-rotation_angle);

    double x_new = input_x + x_shift;
    double y_new = input_y - (slope_shift * center_position + intercept_shift);

    output_x = x_new * cos_angle - y_new * sin_angle;
    output_y = x_new * sin_angle + y_new * cos_angle;
};


void Target_reconstruction(){

ifstream fisier_coeficienti_CATS1X("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS1X.txt");
ifstream fisier_coeficienti_CATS1Y("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS1Y.txt");
ifstream fisier_coeficienti_CATS2X("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS2X.txt");
ifstream fisier_coeficienti_CATS2Y("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS2Y.txt");
TFile *inputRootFile = TFile::Open("/media/olivia/Partition1/CATS/r0421_000a_mycal_nth_gord.root", "READ");

double slope_CATS1X[28];
double slope_CATS1Y[28];

double intercept_CATS1X[28];
double intercept_CATS1Y[28];

double pedestal_CATS1X[28];
double pedestal_CATS1Y[28];

double slope_CATS2X[28];
double slope_CATS2Y[28];

double intercept_CATS2X[28];
double intercept_CATS2Y[28];

double pedestal_CATS2X[28];
double pedestal_CATS2Y[28];

double energy_CATS1X_byStrip[28];
double energy_CATS1Y_byStrip[28];
double energy_CATS2X_byStrip[28];
double energy_CATS2Y_byStrip[28];

double centroid_CATS1X;
double centroid_CATS1Y;
double centroid_CATS2X;
double centroid_CATS2Y;

double slope_shift_CATS1 = 0.02250442496464144;
double intercept_shift_CATS1 = 2.7132336197449174;
double x_shift_CATS1 = 7.046736916766848;
double center_position_CATS1 = -6.5;

double slope_shift_CATS2 = 0.005910027638249295;
double intercept_shift_CATS2 = -0.7013276374635101;
double x_shift_CATS2 = 0.8083158756358418;
double center_position_CATS2 = -9.5;

double X_f, Y_f;

citire_calibrare(fisier_coeficienti_CATS1X, slope_CATS1X, intercept_CATS1X, pedestal_CATS1X);
citire_calibrare(fisier_coeficienti_CATS1Y, slope_CATS1Y, intercept_CATS1Y, pedestal_CATS1Y);
citire_calibrare(fisier_coeficienti_CATS2X, slope_CATS2X, intercept_CATS2X, pedestal_CATS2X);
citire_calibrare(fisier_coeficienti_CATS2Y, slope_CATS2Y, intercept_CATS2Y, pedestal_CATS2Y);

NormalizeToFirstStrip(slope_CATS1X, intercept_CATS1X);
NormalizeToFirstStrip(slope_CATS1Y, intercept_CATS1Y);
NormalizeToFirstStrip(slope_CATS2X, intercept_CATS2X);
NormalizeToFirstStrip(slope_CATS2Y, intercept_CATS2Y);

Int_t   CATS1XVM;
Float_t CATS1XV[28];
UShort_t CATS1XVN[28];

Int_t   CATS1YVM;
Float_t CATS1YV[28];
UShort_t CATS1YVN[28];

Int_t    CATS2XVM;
Float_t  CATS2XV[28];
UShort_t CATS2XVN[28];

Int_t    CATS2YVM;
Float_t  CATS2YV[28];
UShort_t CATS2YVN[28];

ULong64_t DATATRIG_CATS1TS;
ULong64_t DATATRIG_CATS2TS;

Float_t Id_6;
Float_t Id_11;

double deltaZ_CATS21 = 632; //lucram in mm
double deltaZ_targetFromCATS2 = 743;

TTree *catsTree = (TTree*)inputRootFile->Get("AD");
catsTree->SetBranchAddress("CATS1XVM", &CATS1XVM);
catsTree->SetBranchAddress("CATS1XV",  &CATS1XV);
catsTree->SetBranchAddress("CATS1XVN", CATS1XVN);

catsTree->SetBranchAddress("CATS1YVM", &CATS1YVM);
catsTree->SetBranchAddress("CATS1YV",  &CATS1YV);
catsTree->SetBranchAddress("CATS1YVN", &CATS1YVN);

catsTree->SetBranchAddress("CATS2XVM", &CATS2XVM);
catsTree->SetBranchAddress("CATS2XV",  &CATS2XV);
catsTree->SetBranchAddress("CATS2XVN", CATS2XVN);

catsTree->SetBranchAddress("CATS2YVM", &CATS2YVM);
catsTree->SetBranchAddress("CATS2YV",  &CATS2YV);
catsTree->SetBranchAddress("CATS2YVN", &CATS2YVN);


Long64_t entries = catsTree->GetEntries();
TH2F* hTargetXY = new TH2F("hTargetXY", "Target (X_f,Y_f);X_{f};Y_{f}", 1000, -25, 25, 1000, -25, 25);

for (Long64_t entry = 0; entry < entries; ++entry) {

    catsTree->GetEntry(entry);

    if(!(Id_6>5 && Id_6<6.5 && Id_11>269 && Id_11<273))continue; //verificam ca am selectat doar anumiti nuclei 

    //if ( DATATRIG_CATS2TS-DATATRIG_CATS1TS>2000) continue; not needed all events are in the interval

    FillEnergiesByStrip(CATS1XVN, CATS1XV, CATS1XVM, pedestal_CATS1X, slope_CATS1X, intercept_CATS1X, energy_CATS1X_byStrip);
    FillEnergiesByStrip(CATS1YVN, CATS1YV, CATS1YVM, pedestal_CATS1Y, slope_CATS1Y, intercept_CATS1Y, energy_CATS1Y_byStrip);
    FillEnergiesByStrip(CATS2XVN, CATS2XV, CATS2XVM, pedestal_CATS2X, slope_CATS2X, intercept_CATS2X, energy_CATS2X_byStrip);
    FillEnergiesByStrip(CATS2YVN, CATS2YV, CATS2YVM, pedestal_CATS2Y, slope_CATS2Y, intercept_CATS2Y, energy_CATS2Y_byStrip);

    centroid_CATS1X = calculate_centroid(energy_CATS1X_byStrip);
    centroid_CATS1Y = calculate_centroid(energy_CATS1Y_byStrip);
    centroid_CATS2X = calculate_centroid(energy_CATS2X_byStrip);
    centroid_CATS2Y = calculate_centroid(energy_CATS2Y_byStrip);

    if(centroid_CATS1X == -1000 || centroid_CATS1Y == -1000 || centroid_CATS2X == -1000 || centroid_CATS2Y == -1000) continue;

    double corrected_CATS1X, corrected_CATS1Y;
    double corrected_CATS2X, corrected_CATS2Y;

    rotate_and_shift(centroid_CATS1X, centroid_CATS1Y,
                    x_shift_CATS1, slope_shift_CATS1, intercept_shift_CATS1, center_position_CATS1,
                    corrected_CATS1X, corrected_CATS1Y);

    rotate_and_shift(centroid_CATS2X, centroid_CATS2Y,
                    x_shift_CATS2, slope_shift_CATS2, intercept_shift_CATS2, center_position_CATS2,
                    corrected_CATS2X, corrected_CATS2Y);

    centroid_CATS1X = corrected_CATS1X;
    centroid_CATS1Y = corrected_CATS1Y;
    centroid_CATS2X = corrected_CATS2X;
    centroid_CATS2Y = corrected_CATS2Y;

    X_f = centroid_CATS2X + ( (centroid_CATS2X - centroid_CATS1X) / deltaZ_CATS21 ) * deltaZ_targetFromCATS2;
    Y_f = centroid_CATS2Y + ( (centroid_CATS2Y - centroid_CATS1Y) / deltaZ_CATS21 ) * deltaZ_targetFromCATS2;

   // hTargetXY->Fill(X_f, Y_f);
    hTargetXY->Fill(centroid_CATS1X,centroid_CATS1Y);

}

TCanvas* canvasTargetXY = new TCanvas("canvasTargetXY", "Target XY", 900, 800);
hTargetXY->Draw("COLZ");
canvasTargetXY->Update();
gPad->WaitPrimitive();

}