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

    for(int index = 2; index < arrayLength-1; index++)
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
{//TODO scapa de jegul de hardcodare pentru stripurile care nu merg
    for(int strip_index = 1; strip_index < 27; ++strip_index) energyByStripArray[strip_index] = 0.0;

    for(int hit_index = 0; hit_index < multiplicity; ++hit_index)
    {
        int stripIndex = (int)stripNumberArray[hit_index];

        if(stripIndex == 0 || stripIndex > 26) continue;
        if((double)rawValueArray[hit_index] < pedestalArray[stripIndex]) continue;
        //TODO Formula completa contine si inercept dar noi l-am luat zero (raw-ped)^2
        //TODO de ce se adauga acel random intre 0.5 si de ce isi iau ei intervalele constante?
        double calibratedEnergy = ((double)rawValueArray[hit_index] - pedestalArray[stripIndex]) * slopeArray[stripIndex] + interceptArray[stripIndex];
        double random_value = ((double)rand() / (double)RAND_MAX) - 0.5;
        calibratedEnergy+=random_value;

        energyByStripArray[stripIndex] += calibratedEnergy;
    }
}


//TODO nu prea conteaza daca dau skip la 0 si 27 ca oricum sunt zero din calibrare 
// for the reconstruction algorithms we need to get the ordered values of the strip energies and select only the nstrip highest ones
std::vector<int> build_ordered_strip_energies(const double strip_energies[28]){
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

void shift_centroid(double input_x, double input_y,
                    double x_shift, double y_shift,
                    double& output_x, double& output_y)
{
    output_x = input_x + x_shift;
    output_y = input_y + y_shift;
}



double weighted_average(const double strip_energies[28], const std::vector<int>& fstrip)
{
    //TODO scapa de hardcodare ptr n strip
    if(fstrip.size() < 3) return -1000.0;

    int max_strip_index = fstrip[0];
    if(max_strip_index < 1 || max_strip_index > 26) return -1000.0;

    int min_strip_index = max_strip_index - 2;
    int max_window_index = max_strip_index + 2;

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


double calculate_centroid(double strip_energies[28], vector<int>fstrip){
    int m = (int)fstrip.size();
    if(m < 1) return -1000.0;

    double x = -1000.0;

    x = weighted_average(strip_energies, fstrip);
  
    if(x == -1000.0) return -1000.0;
    return ((x - 13.5) * 2.54);
}

void Target_reconstruction(){

ifstream fisier_coeficienti_CATS1X("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS1X.txt");
ifstream fisier_coeficienti_CATS1Y("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS1Y.txt");
ifstream fisier_coeficienti_CATS2X("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS2X.txt");
ifstream fisier_coeficienti_CATS2Y("/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS2Y.txt");

//const char* file_address = "/media/olivia/Partition1/CATS/r1159_00*a.root";
TString runTag = "r1159";
TString fileAddressPattern = Form("/media/olivia/Partition1/CATS/%s_00*a.root", runTag.Data());
const char* file_address = fileAddressPattern.Data();

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
// double slope_shift_CATS1 = 0.02250442496464144;
// double intercept_shift_CATS1 = 2.7132336197449174;
// double x_shift_CATS1 = 7.046736916766848;
// double center_position_CATS1 = -6.5;

// double slope_shift_CATS2 = 0.005910027638249295;
// double intercept_shift_CATS2 = -0.7013276374635101;
// double x_shift_CATS2 = 0.8083158756358418;
// double center_position_CATS2 = -9.5;

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
catsTree->SetBranchAddress("CATS2YVN", CATS2YVN);

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

TH2F* hTargetXY = new TH2F("hTargetXY", "Target (X_f,Y_f);X_{f};Y_{f}", 1000, -40, 40, 1000, -40, 40);

TH2F* hC1_centroid = new TH2F("hC1_centroid", "CATS1 centroid;X (mm);Y (mm)", 1000, -40, 40, 1000, -40, 40);
TH2F* hC1_shifted  = new TH2F("hC1_shifted",  "CATS1 shifted;X (mm);Y (mm)", 1000, -40, 40, 1000, -40, 40);
TH2F* hC1_rotated  = new TH2F("hC1_rotated",  "CATS1 rotated;X (mm);Y (mm)", 1000, -40, 40, 1000, -40, 40);

TH2F* hC2_centroid = new TH2F("hC2_centroid", "CATS2 centroid;X (mm);Y (mm)", 1000, -40, 40, 1000, -40, 40);
TH2F* hC2_shifted  = new TH2F("hC2_shifted",  "CATS2 shifted;X (mm);Y (mm)", 1000, -40, 40, 1000, -40, 40);

TH1F* hTargetXY2 = new TH1F("h2", "CATS1 X centroid;X (mm);Counts", 1000, -40, 40);
TH1F* hTargetXY3 = new TH1F("h3", "CATS1 Y centroid;Y (mm);Counts", 1000, -40, 40);
TH1F* hTargetXY4 = new TH1F("h4", "CATS2 X centroid;X (mm);Counts", 1000, -40, 40);
TH1F* hTargetXY5 = new TH1F("h5", "CATS2 Y centroid;Y (mm);Counts", 1000, -40, 40);

std::cout << entries << '\n';
for (Long64_t entry = 0; entry < entries; ++entry) {

    catsTree->GetEntry(entry);

    //990
   // if(!(Id_6>4.68 && Id_6<5.45 && Id_11>222.2 && Id_11<225))continue; //verificam ca am selectat doar anumiti nuclei 
    //948
  //if(!(Id_6>6.55 && Id_6<7.65 && Id_11>226.8 && Id_11<228.96))continue; //verificam ca am selectat doar anumiti nuclei 
   //1159
   if(!(Id_6>6.14 && Id_6<7.05 && Id_11>224.49 && Id_11<226.71))continue;

    FillEnergiesByStrip(CATS1XVN, CATS1XV, CATS1XVM, pedestal_CATS1X, slope_CATS1X, intercept_CATS1X, energy_CATS1X_byStrip);
    FillEnergiesByStrip(CATS1YVN, CATS1YV, CATS1YVM, pedestal_CATS1Y, slope_CATS1Y, intercept_CATS1Y, energy_CATS1Y_byStrip);
    FillEnergiesByStrip(CATS2XVN, CATS2XV, CATS2XVM, pedestal_CATS2X, slope_CATS2X, intercept_CATS2X, energy_CATS2X_byStrip);
    FillEnergiesByStrip(CATS2YVN, CATS2YV, CATS2YVM, pedestal_CATS2Y, slope_CATS2Y, intercept_CATS2Y, energy_CATS2Y_byStrip);

    std::vector<int> fstrip_CATS1X = build_ordered_strip_energies(energy_CATS1X_byStrip);
    std::vector<int> fstrip_CATS1Y = build_ordered_strip_energies(energy_CATS1Y_byStrip);
    std::vector<int> fstrip_CATS2X = build_ordered_strip_energies(energy_CATS2X_byStrip);
    std::vector<int> fstrip_CATS2Y = build_ordered_strip_energies(energy_CATS2Y_byStrip);


    //CENTROID DETERMINATION 
    double corrected_CATS1X, corrected_CATS1Y;
   

    //scoase din elog si din codul lor
    double rotation_angle_cats1_degrees = -1.2;
    const double Pi = 3.14159265358979323846;
    double x_shift_CATS1 =6.5;
    double y_shift_CATS1 =-1.705;
    double x_shift_CATS2 = 0.188;
    double y_shift_CATS2 =1.509;

    //din teza QUentin 6584.8-6037.7
    //TODO pare ciudat ca da negativ dar asa e formula 
    double deltaZ_CATS21 = -547.1; //lucram in mm
    //din fp ref + teza quentin facut diferenta 7600-6584.8
    double deltaZ_targetFromCATS2 = 1015.2;
    double target_angle = 36.7; //tinta nu era dreapta, trebuie corectate distantele
    //din ganil sheet


    centroid_CATS1X = calculate_centroid(energy_CATS1X_byStrip, fstrip_CATS1X);
    centroid_CATS1Y = calculate_centroid(energy_CATS1Y_byStrip, fstrip_CATS1Y);
    centroid_CATS2X = calculate_centroid(energy_CATS2X_byStrip, fstrip_CATS2X);
    centroid_CATS2Y = calculate_centroid(energy_CATS2Y_byStrip, fstrip_CATS2Y);

    if(centroid_CATS1X == -1000 || centroid_CATS1Y == -1000 || centroid_CATS2X == -1000 || centroid_CATS2Y == -1000) continue;

    double shifted_CATS1X,shifted_CATS1Y,shifted_CATS2X,shifted_CATS2Y;

    shift_centroid(centroid_CATS1X, centroid_CATS1Y, x_shift_CATS1, y_shift_CATS1, shifted_CATS1X, shifted_CATS1Y);
    shift_centroid(centroid_CATS2X, centroid_CATS2Y, x_shift_CATS2, y_shift_CATS2, shifted_CATS2X, shifted_CATS2Y);

    corrected_CATS1X = shifted_CATS1X * cos(rotation_angle_cats1_degrees * Pi / 180.0)
                    - shifted_CATS1Y * sin(rotation_angle_cats1_degrees * Pi / 180.0);

    X_f = shifted_CATS2X + ((corrected_CATS1X - shifted_CATS2X) / deltaZ_CATS21) * deltaZ_targetFromCATS2;

    corrected_CATS1Y = shifted_CATS1Y * cos(rotation_angle_cats1_degrees * Pi / 180.0)
                    + shifted_CATS1X * sin(rotation_angle_cats1_degrees * Pi / 180.0);

    Y_f = shifted_CATS2Y + ((corrected_CATS1Y - shifted_CATS2Y) / deltaZ_CATS21) * deltaZ_targetFromCATS2;

    double Z_f = X_f * tan(target_angle * Pi / 180.0);

    double X_f_corr = shifted_CATS2X
                    + ((corrected_CATS1X - shifted_CATS2X) / deltaZ_CATS21)
                    * (deltaZ_targetFromCATS2 + Z_f);

    double Y_f_corr = shifted_CATS2Y
                    + ((corrected_CATS1Y - shifted_CATS2Y) / deltaZ_CATS21)
                    * (deltaZ_targetFromCATS2 + Z_f);

    hTargetXY->Fill(X_f_corr, Y_f_corr);

    hC1_centroid->Fill(centroid_CATS1X, centroid_CATS1Y);
    hC1_shifted->Fill(shifted_CATS1X, shifted_CATS1Y);
    hC1_rotated->Fill(corrected_CATS1X, corrected_CATS1Y);

    hC2_centroid->Fill(centroid_CATS2X, centroid_CATS2Y);
    hC2_shifted->Fill(shifted_CATS2X, shifted_CATS2Y);

    hTargetXY2->Fill(centroid_CATS1X);
    hTargetXY3->Fill(centroid_CATS1Y);
    hTargetXY4->Fill(centroid_CATS2X);
    hTargetXY5->Fill(centroid_CATS2Y);

    //std:cout<<centroid_CATS2X<<" "<<centroid_CATS2Y;
}

TCanvas* canvasTargetXY = new TCanvas("canvasTargetXY", "Target (c2,c2)", 900, 800);
hTargetXY->Draw("COLZ");
canvasTargetXY->Update();

TCanvas* canvasC1X = new TCanvas("canvasC1X", "C1X", 900, 800);
hTargetXY2->Draw();
canvasC1X->Update();

TCanvas* canvasC1Y = new TCanvas("canvasC1Y", "C1Y", 900, 800);
hTargetXY3->Draw();
canvasC1Y->Update();

TCanvas* canvasC2X = new TCanvas("canvasC2X", "C2X", 900, 800);
hTargetXY4->Draw();
canvasC2X->Update();

TCanvas* canvasC2Y = new TCanvas("canvasC2Y", "C2Y", 900, 800);
hTargetXY5->Draw();
canvasC2Y->Update();

gSystem->ProcessEvents();

TFile* outFile = new TFile(Form("Target_Reconstruction_CATS_%s.root", runTag.Data()), "RECREATE");
outFile->cd();

hTargetXY->SetName("h_XfYf");
hTargetXY->SetTitle("XfYf;X_{f} (mm);Y_{f} (mm)");
hTargetXY->Write();

hC1_centroid->Write();
hC1_shifted->Write();
hC1_rotated->Write();
hC2_centroid->Write();
hC2_shifted->Write();

hTargetXY2->SetName("h_C1X");
hTargetXY2->SetTitle("C1X;CATS1 X (mm);Counts");
hTargetXY2->Write();

hTargetXY3->SetName("h_C1Y");
hTargetXY3->SetTitle("C1Y;CATS1 Y (mm);Counts");
hTargetXY3->Write();

hTargetXY4->SetName("h_C2X");
hTargetXY4->SetTitle("C2X;CATS2 X (mm);Counts");
hTargetXY4->Write();

hTargetXY5->SetName("h_C2Y");
hTargetXY5->SetTitle("C2Y;CATS2 Y (mm);Counts");
hTargetXY5->Write();

outFile->Close();

std::cin.get();

}