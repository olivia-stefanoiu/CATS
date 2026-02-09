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
#include <math.h>


void citire_calibrare(std::ifstream &file, double slope[28], double intercept[28],double pedestal[28])
{
    //Calibram per strip pana ajung toate la acelasi nivel. Pedestalul este foarte important pentru eliminarea zgomotolui
    int strip;
    for (int i = 0; i < 28; ++i){
        file >>strip>> slope[i] >> intercept[i]>>pedestal[i];
        cout<<strip<<" "<<slope[i]<<" "<<intercept[i]<<" "<<pedestal[i]<<'\n';
    }
    cout<<'\n';
}

//Normalizam la un strip la alegere valorile celorlalte
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

//selectam doar stripurile care au valori peste pedestal si calibram 
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
        if((double)rawValueArray[hit_index] < pedestalArray[stripIndex]+200) continue;
        //TODO Formula completa contine si inercept dar noi l-am luat zero (raw-ped)^2
        double calibratedEnergy = ((double)rawValueArray[hit_index] - pedestalArray[stripIndex]) * slopeArray[stripIndex] + interceptArray[stripIndex];
        energyByStripArray[stripIndex] += calibratedEnergy;
    }
}
// for the reconstruction algorithms we need to get the ordered values of the strip energies and select only the nstrip highest ones
std::vector<int> build_ordered_strip_energies(const double strip_energies[28]){
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

double sechip(const double strip_energies[28], const std::vector<int>& fstrip){
    if((int)fstrip.size() < 3) return -1000.0;

    int s0 = fstrip[0];
    int s1 = fstrip[1];
    int s2 = fstrip[2];

    double q0 = strip_energies[s0];
    double q1 = strip_energies[s1];
    double q2 = strip_energies[s2];

    if(q0 <= 0.0 || q1 <= 0.0 || q2 <= 0.0) return -1000.0;

    double v0 = std::sqrt(q0 / q2);
    double v1 = std::sqrt(q0 / q1);
    double v2 = 0.5 * (v0 + v1);
    if(v2 < 1.0) return -1000.0;

    double v3 = std::log(v2 + std::sqrt(v2 * v2 - 1.0));
    if(v3 == 0.0) return -1000.0;

    double v4 = (v0 - v1) / (2.0 * std::sinh(v3));
    if(std::fabs(v4) >= 1.0) return -1000.0;

    double v5 = 0.5 * std::log((1.0 + v4) / (1.0 - v4));

    double xs = (double)s0 - (double)(s0 - s1) * v5 / v3;
    if(xs <= 0.0 || xs >= 28.0) return -1000.0;

    return xs;
}


//Weighted avg has the best results for nstrip = 3 respectively 5. Adding more strips does not increase performance 
// double weighted_average(const double strip_energies[28], const std::vector<int>& ordered_strip_energies, int nstrips=3){
//     //verific deja in fct calculate_centroid
//     if((int)ordered_strip_energies.size()<=1) return -1000.0;
//     if((int)ordered_strip_energies.size() < nstrips) nstrips = (int)ordered_strip_energies.size();

//     double numerator = 0.0;
//     double denominator = 0.0;

   
    
//     for(int rank_index = 0; rank_index < nstrips; rank_index++){
//         int strip_index = ordered_strip_energies[rank_index];
//         double energy = strip_energies[strip_index];
//         numerator += (double)strip_index * energy;
//         denominator += energy;
//     }

//     if(denominator <= 0.0) return -1000.0;
//     return numerator / denominator;
// }

// bool check_neighbours(std::vector<int>& fstrip, const double strip_energies[28]){
//     if((int)fstrip.size() < 3) return false;

//     int center_strip = fstrip[0];

//     auto is_neighbour = [&](int strip_index){
//         return (std::abs(strip_index - center_strip) == 1) && (strip_energies[strip_index] > 0.0);
//     };

//     if(!is_neighbour(fstrip[1])) return false;
//     if(is_neighbour(fstrip[2])) return true;

//     if((int)fstrip.size() > 3 && is_neighbour(fstrip[3])){
//         int tmp = fstrip[2]; fstrip[2] = fstrip[3]; fstrip[3] = tmp;
//         return true;
//     }

//     if((int)fstrip.size() > 4 && is_neighbour(fstrip[4])){
//         int tmp = fstrip[2]; fstrip[2] = fstrip[4]; fstrip[4] = tmp;
//         return true;
//     }

//     return false;
// }


// double calculate_centroid(double strip_energies[28]){
//     std::vector<int> ordered_strip_energies = build_ordered_strip_energies(strip_energies);
//     int multiplicity = (int)ordered_strip_energies.size();
  
//     double x = -1000.0;
//    // x = weighted_average(strip_energies,ordered_strip_energies, 3);
//     if (!check_neighbours(ordered_strip_energies,strip_energies))return -1000;
//     x = sechip(strip_energies,ordered_strip_energies);
//     if(x == -1000.0) return -1000.0;
//     //TODO verifica daca de aici e problema cu gridul 
//     return ((x - 13.5) * 2.54);
// }

void rotate_and_shift_centroid(double input_x, double input_y,
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
}

bool check_neighbours(std::vector<int>& fstrip)
{
    if((int)fstrip.size() < 3) return false;

    int tmp = 0;
    int s0 = fstrip[0];

    if(std::abs(s0 - fstrip[1]) == 1 &&
       std::abs(s0 - fstrip[2]) == 1)
        return true;

    bool rval = false;

    if((int)fstrip.size() > 3)
    {
        if(std::abs(s0 - fstrip[1]) == 1 &&
           std::abs(s0 - fstrip[3]) == 1)
        {
            tmp = fstrip[2];
            fstrip[2] = fstrip[3];
            fstrip[3] = tmp;
            rval = true;
        }
    }

    if(!rval && (int)fstrip.size() > 4)
    {
        if(std::abs(s0 - fstrip[1]) == 1 &&
           std::abs(s0 - fstrip[4]) == 1)
        {
            tmp = fstrip[2];
            fstrip[2] = fstrip[4];
            fstrip[4] = tmp;
            rval = true;
        }
    }

    return rval;
}

double weighted_average(const double strip_energies[28],
                        const std::vector<int>& ordered_strip_energies,
                        int nstrips,
                        int number_of_detectors = 28,
                        double invalid_value = -1000.0)
{
    if(nstrips <= 0) return invalid_value;
    if((int)ordered_strip_energies.size() < nstrips) return invalid_value;

    double v0 = 0.0;
    double v1 = 0.0;

    for(int k = 0; k < nstrips; k++)
    {
        int strip_number = ordered_strip_energies[k];
        if(strip_number < 0 || strip_number >= number_of_detectors) continue;

        double Q = strip_energies[strip_number];
        double N = (double)strip_number;

        v0 += Q * N;
        v1 += Q;
    }

    if(v1 == 0.0) return invalid_value;

    double xwa = v0 / v1;

    if(xwa > 0.0 && xwa < (double)number_of_detectors) return xwa;
    return invalid_value;
}

double calculate_centroid(double strip_energies[28],
                          const std::vector<int>& fstrip_input,
                          int number_of_detectors = 28,
                          double invalid_value = -1000.0,
                          double single_hit_threshold = 500.0,
                          double ref_position_mm = 0.0)
{
    int m = (int)fstrip_input.size();
    if(m <= 0) return invalid_value;

    double x_temp_strip = invalid_value;

    if(m >= 3)
    {
        std::vector<int> fstrip = fstrip_input;

        if(check_neighbours(fstrip))
        {
            x_temp_strip = sechip(strip_energies, fstrip);
        }
        else
        {
            x_temp_strip = invalid_value;
        }
    }
//     else if(m == 2)
//     {
//         x_temp_strip = weighted_average(strip_energies, fstrip_input, 2, number_of_detectors, invalid_value);
//     }
//     else
//     {
//         int strip0 = fstrip_input[0];
//         if(strip0 >= 0 && strip0 < number_of_detectors)
//         {
//             if(strip_energies[strip0] > single_hit_threshold)
//                 x_temp_strip = (double)strip0;
//         }
//     }

    if(!(x_temp_strip > 0.0 && x_temp_strip < (double)number_of_detectors))
        return invalid_value;

    double centered_strip = x_temp_strip - ((double)number_of_detectors / 2.0 - 0.5);
    return centered_strip * 2.54 + ref_position_mm;
}


void Target_reconstruction(){

ifstream fisier_coeficienti_CATS1X("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS1X.txt");
ifstream fisier_coeficienti_CATS1Y("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS1Y.txt");
ifstream fisier_coeficienti_CATS2X("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS2X.txt");
ifstream fisier_coeficienti_CATS2Y("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS2Y.txt");
const char* file_address = "/media/olivia/Partition1/CATS/r0948_00*a.root";

//TFile *inputRootFile = TFile::Open("/media/olivia/Partition1/CATS/r1163_000a.root", "READ");

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

//TODO no more hardcoding, read from file
//TODO recalculeaza ptr Cats1 coeficientii ca e stramb
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

//TTree *catsTree = (TTree*)inputRootFile->Get("AD");
TChain* catsTree = new TChain("AD");
catsTree->Add(file_address);


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

catsTree->SetBranchAddress("Id_6", &Id_6);
catsTree->SetBranchAddress("Id_11", &Id_11);

// TCanvas* cId = new TCanvas("cId","Id_6 vs Id_11",900,800);
// catsTree->Draw("Id_11:Id_6>>hId11Id6(1000,0,10,1000,225,235)","","COLZ");
// cId->Update();
// std::cout << "Press Enter..." << std::flush;
// std::cin.get();

// TCanvas* cSame = new TCanvas("cSame","CATS1XVN and CATS2XVN",900,800);
// catsTree->Draw("CATS1XVN","Id_6 > 6.4 && Id_6 < 8 && Id_11>226.5 && Id_11<229.5 && CATS1XV>1000");
// catsTree->Draw("CATS2XVN","Id_6 > 6.4 && Id_6 < 8 && Id_11>226.5 && Id_11<229.5 && CATS2XV>1000","SAME");
// cSame->Update();
// std::cout << "Press Enter..." << std::flush;
// std::cin.get();


Long64_t entries = catsTree->GetEntries();
TH2F* hTargetXY = new TH2F("hTargetXY", "Target (X_f,Y_f);X_{f};Y_{f}", 1000, -40,40, 1000, -40, 40);

for (Long64_t entry = 0; entry < entries; ++entry) {

    catsTree->GetEntry(entry);

   //if(!(Id_6>6.4 && Id_6<6.7 && Id_11>225.3 && Id_11<225.9))continue; //verificam ca am selectat doar anumiti nuclei 
   //if(!(Id_6>4 && Id_6<22 && Id_11>220 && Id_11<235))continue;
    //if ( DATATRIG_CATS2TS-DATATRIG_CATS1TS>2000) continue; not needed all events are in the interval

    //948 SULF
   if(!(Id_6>6.4 && Id_6<8 && Id_11>226.5 && Id_11<229.5))continue;
   //Beamul se shifteaza la SI
   //if(!(Id_6>0&& Id_6<10 && Id_11>220 && Id_11<230))continue;

    // FillEnergiesByStrip(CATS1XVN, CATS1XV, CATS1XVM, pedestal_CATS1X, slope_CATS1X, intercept_CATS1X, energy_CATS1X_byStrip);
    // FillEnergiesByStrip(CATS1YVN, CATS1YV, CATS1YVM, pedestal_CATS1Y, slope_CATS1Y, intercept_CATS1Y, energy_CATS1Y_byStrip);
    FillEnergiesByStrip(CATS2XVN, CATS2XV, CATS2XVM, pedestal_CATS2X, slope_CATS2X, intercept_CATS2X, energy_CATS2X_byStrip);
    FillEnergiesByStrip(CATS2YVN, CATS2YV, CATS2YVM, pedestal_CATS2Y, slope_CATS2Y, intercept_CATS2Y, energy_CATS2Y_byStrip);

//     std::vector<int> fstrip_CATS1X = build_ordered_strip_energies(energy_CATS1X_byStrip);
// std::vector<int> fstrip_CATS1Y = build_ordered_strip_energies(energy_CATS1Y_byStrip);
std::vector<int> fstrip_CATS2X = build_ordered_strip_energies(energy_CATS2X_byStrip);
std::vector<int> fstrip_CATS2Y = build_ordered_strip_energies(energy_CATS2Y_byStrip);


//     centroid_CATS1X = calculate_centroid(energy_CATS1X_byStrip);
//     centroid_CATS1Y = calculate_centroid(energy_CATS1Y_byStrip);
//    centroid_CATS2X = calculate_centroid(energy_CATS2X_byStrip);
//    centroid_CATS2Y = calculate_centroid(energy_CATS2Y_byStrip);

// centroid_CATS1X = calculate_centroid(energy_CATS1X_byStrip, fstrip_CATS1X);
// centroid_CATS1Y = calculate_centroid(energy_CATS1Y_byStrip, fstrip_CATS1Y);
centroid_CATS2X = calculate_centroid(energy_CATS2X_byStrip, fstrip_CATS2X);
centroid_CATS2Y = calculate_centroid(energy_CATS2Y_byStrip, fstrip_CATS2Y);


   // if(centroid_CATS1X == -1000 || centroid_CATS1Y == -1000 || centroid_CATS2X == -1000 || centroid_CATS2Y == -1000) continue;
   if(centroid_CATS2X == -1000 || centroid_CATS2Y == -1000 ) continue;

    double corrected_CATS1X, corrected_CATS1Y;
    double corrected_CATS2X, corrected_CATS2Y;

    // rotate_and_shift_centroid(centroid_CATS1X, centroid_CATS1Y,
    //                       x_shift_CATS1, slope_shift_CATS1, intercept_shift_CATS1, center_position_CATS1,
    //                       corrected_CATS1X, corrected_CATS1Y);

    rotate_and_shift_centroid(centroid_CATS2X, centroid_CATS2Y,
                          x_shift_CATS2, slope_shift_CATS2, intercept_shift_CATS2, center_position_CATS2,
                          corrected_CATS2X, corrected_CATS2Y);
    
//    if(centroid_CATS1X>5 && centroid_CATS1X<6 )
//     hTargetXY->Fill(centroid_CATS2X,centroid_CATS2Y);

    // centroid_CATS1X = corrected_CATS1X;
    // centroid_CATS1Y = corrected_CATS1Y;
    centroid_CATS2X = corrected_CATS2X;
    centroid_CATS2Y = corrected_CATS2Y;

    // X_f = centroid_CATS2X + ( (centroid_CATS2X - centroid_CATS1X) / deltaZ_CATS21 ) * deltaZ_targetFromCATS2;
    // Y_f = centroid_CATS2Y + ( (centroid_CATS2Y - centroid_CATS1Y) / deltaZ_CATS21 ) * deltaZ_targetFromCATS2;

   hTargetXY->Fill(centroid_CATS2X,centroid_CATS2Y);
   //hTargetXY->Fill(X_f, Y_f);
  //std:cout<<centroid_CATS2X<<" "<<centroid_CATS2Y;
    
}

TCanvas* canvasTargetXY = new TCanvas("canvasTargetXY", "Target XY", 900, 800);
hTargetXY->Draw("COLZ");
canvasTargetXY->Update();
gPad->WaitPrimitive();

std::cout<<"FINALIZAT";

}