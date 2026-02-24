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


//TString tags[] = {"r0927","r0944","r0948","r1159"};
TString tags[] = {"r0964","r0989","r1010","r1041","r1048","r1064","r1095","r1099","r1112"};
int number_of_tags = sizeof(tags) / sizeof(tags[0]);

for (int tag_index = 0; tag_index < number_of_tags; tag_index++) {

  TString runTag = tags[tag_index];
  std::cout << runTag << " " << std::endl;

  TString fileAddressPattern = Form("/media/olivia/Partition1/CATS/Remerged/%s_00*a.root", runTag.Data());

  TChain* catsTree = new TChain("AD");
  catsTree->Add(fileAddressPattern.Data());

  
    TCanvas* cId = new TCanvas(Form("cId_%s", runTag.Data()), "Id_6 vs Id_11", 900, 800);
    catsTree->Draw(Form("Id_11:Id_6>>hId11Id6_%s(1000,2,10,1000,220,230)", runTag.Data()), "", "COLZ");
    cId->Update();

    TH2* h2 = (TH2*)gDirectory->Get(Form("hId11Id6_%s", runTag.Data()));
    if (!h2) { std::cerr << "missing histogram for " << runTag << "\n"; continue; }

    Int_t globalBin = h2->GetMaximumBin();
    Int_t binx = 0, biny = 0, binz = 0;
    h2->GetBinXYZ(globalBin, binx, biny, binz);

    Double_t xCenter = h2->GetXaxis()->GetBinCenter(binx);
    Double_t yCenter = h2->GetYaxis()->GetBinCenter(biny);

    TCanvas* cSame = new TCanvas(Form("cSame_%s", runTag.Data()), "Window around max bin", 900, 800);
    catsTree->Draw("Id_11:Id_6", Form("Id_6 > %g && Id_6 < %g && Id_11 > %g && Id_11 < %g && Max$(CATS1XV) > %g", xCenter-0.75, xCenter+0.75, yCenter-1.75, yCenter+1.75, 300.0), "COLZ");
    cSame->Update();

    std::cout<< xCenter<< " "<<yCenter<<'\n';
  

  delete catsTree;
}

std::cout << "Press Enter..." << std::flush;
std::cin.get();
// Long64_t entries = catsTree->GetEntries();

// TH2F* hTargetXY = new TH2F("hTargetXY", "Target (X_f,Y_f);X_{f};Y_{f}", 1000, -40, 40, 1000, -40, 40);

// TH2F* hTargetXY_simple = new TH2F("hTargetXY_simple", "Target simple (X_f,Y_f);X_{f};Y_{f}", 1000, -40, 40, 1000, -40, 40);

// TH1F* hXf = new TH1F("hXf", "X_{f} (corr);X_{f} (mm);Counts", 1000, -40, 40);
// TH1F* hYf = new TH1F("hYf", "Y_{f} (corr);Y_{f} (mm);Counts", 1000, -40, 40);

// TH1F* hXf_simple = new TH1F("hXf_simple", "X_{f} simple;X_{f} (mm);Counts", 1000, -40, 40);
// TH1F* hYf_simple = new TH1F("hYf_simple", "Y_{f} simple;Y_{f} (mm);Counts", 1000, -40, 40);

// TH2F* hC1_centroid = new TH2F("hC1_centroid", "CATS1 centroid;X (mm);Y (mm)", 1000, -40, 40, 1000, -40, 40);
// TH2F* hC1_shifted  = new TH2F("hC1_shifted",  "CATS1 shifted;X (mm);Y (mm)", 1000, -40, 40, 1000, -40, 40);
// TH2F* hC1_rotated  = new TH2F("hC1_rotated",  "CATS1 rotated;X (mm);Y (mm)", 1000, -40, 40, 1000, -40, 40);

// TH2F* hC2_centroid = new TH2F("hC2_centroid", "CATS2 centroid;X (mm);Y (mm)", 1000, -40, 40, 1000, -40, 40);
// TH2F* hC2_shifted  = new TH2F("hC2_shifted",  "CATS2 shifted;X (mm);Y (mm)", 1000, -40, 40, 1000, -40, 40);

// TH1F* hTargetXY2 = new TH1F("h2", "CATS1 X centroid;X (mm);Counts", 1000, -40, 40);
// TH1F* hTargetXY3 = new TH1F("h3", "CATS1 Y centroid;Y (mm);Counts", 1000, -40, 40);
// TH1F* hTargetXY4 = new TH1F("h4", "CATS2 X centroid;X (mm);Counts", 1000, -40, 40);
// TH1F* hTargetXY5 = new TH1F("h5", "CATS2 Y centroid;Y (mm);Counts", 1000, -40, 40);

// std::cout << entries << '\n';
// for (Long64_t entry = 0; entry < entries/10; ++entry) {

//     catsTree->GetEntry(entry);

//     //990
//   //  if(!(Id_6>4.68 && Id_6<5.45 && Id_11>222.2 && Id_11<225))continue; //verificam ca am selectat doar anumiti nuclei 
//     //948
//   //  if(!(Id_6>6.55 && Id_6<7.65 && Id_11>226.8 && Id_11<228.96))continue; //verificam ca am selectat doar anumiti nuclei 
//    //1159
//    //if(!(Id_6>6.14 && Id_6<7.05 && Id_11>224.49 && Id_11<226.71))continue;
//    //1010
//     if(!(Id_6>3.83 && Id_6<5.06 && Id_11>221.57 && Id_11<224.71))continue;

//     FillEnergiesByStrip(CATS1XVN, CATS1XV, CATS1XVM, pedestal_CATS1X, slope_CATS1X, intercept_CATS1X, energy_CATS1X_byStrip);
//     FillEnergiesByStrip(CATS1YVN, CATS1YV, CATS1YVM, pedestal_CATS1Y, slope_CATS1Y, intercept_CATS1Y, energy_CATS1Y_byStrip);
//     FillEnergiesByStrip(CATS2XVN, CATS2XV, CATS2XVM, pedestal_CATS2X, slope_CATS2X, intercept_CATS2X, energy_CATS2X_byStrip);
//     FillEnergiesByStrip(CATS2YVN, CATS2YV, CATS2YVM, pedestal_CATS2Y, slope_CATS2Y, intercept_CATS2Y, energy_CATS2Y_byStrip);

//     std::vector<int> fstrip_CATS1X = build_ordered_strip_energies(energy_CATS1X_byStrip);
//     std::vector<int> fstrip_CATS1Y = build_ordered_strip_energies(energy_CATS1Y_byStrip);
//     std::vector<int> fstrip_CATS2X = build_ordered_strip_energies(energy_CATS2X_byStrip);
//     std::vector<int> fstrip_CATS2Y = build_ordered_strip_energies(energy_CATS2Y_byStrip);

//     //CENTROID DETERMINATION 
//     double corrected_CATS1X, corrected_CATS1Y;
   
//     //scoase din elog si din codul lor
//     double rotation_angle_cats1_degrees = -1.2;
//     const double Pi = 3.14159265358979323846;
//     double x_shift_CATS1 =7.03;
//     double y_shift_CATS1 =-2.28;
//     double x_shift_CATS2 = 0.805+2;  ///0.188; +2 sa incerc sa scap de 2 struct
//     double y_shift_CATS2 =0.627;

//     //din teza QUentin 6584.8-6037.7
//     //TODO pare ciudat ca da negativ dar asa e formula 
//     double deltaZ_CATS21 = -547.1; //lucram in mm
//     //din fp ref + teza quentin facut diferenta 7600-6584.8
//     double deltaZ_targetFromCATS2 = 1015.2;
//     double target_angle = 36.7; //tinta nu era dreapta, trebuie corectate distantele
//     //din ganil sheet

//     centroid_CATS1X = calculate_centroid(energy_CATS1X_byStrip, fstrip_CATS1X);
//     centroid_CATS1Y = calculate_centroid(energy_CATS1Y_byStrip, fstrip_CATS1Y);
//     centroid_CATS2X = calculate_centroid(energy_CATS2X_byStrip, fstrip_CATS2X);
//     centroid_CATS2Y = calculate_centroid(energy_CATS2Y_byStrip, fstrip_CATS2Y);

//     if(centroid_CATS1X == -1000 || centroid_CATS1Y == -1000 || centroid_CATS2X == -1000 || centroid_CATS2Y == -1000) continue;

//     // if(centroid_CATS1X>5 && centroid_CATS1X<6 )
//     //   hTargetXY->Fill(centroid_CATS2X,centroid_CATS2Y);


//     double shifted_CATS1X,shifted_CATS1Y,shifted_CATS2X,shifted_CATS2Y;

//     shift_centroid(centroid_CATS1X, centroid_CATS1Y, x_shift_CATS1, y_shift_CATS1, shifted_CATS1X, shifted_CATS1Y);
//     shift_centroid(centroid_CATS2X, centroid_CATS2Y, x_shift_CATS2, y_shift_CATS2, shifted_CATS2X, shifted_CATS2Y);

//     corrected_CATS1X = shifted_CATS1X * cos(rotation_angle_cats1_degrees * Pi / 180.0)
//                     - shifted_CATS1Y * sin(rotation_angle_cats1_degrees * Pi / 180.0);


//     X_f = shifted_CATS2X + ((corrected_CATS1X - shifted_CATS2X) / deltaZ_CATS21) * deltaZ_targetFromCATS2;

//     corrected_CATS1Y = shifted_CATS1Y * cos(rotation_angle_cats1_degrees * Pi / 180.0)
//                     + shifted_CATS1X * sin(rotation_angle_cats1_degrees * Pi / 180.0);

//    Y_f = shifted_CATS2Y + ((corrected_CATS1Y - shifted_CATS2Y) / deltaZ_CATS21) * deltaZ_targetFromCATS2;

//       //  if (corrected_CATS1X>2 && corrected_CATS1X<10)
//       //   hTargetXY->Fill(corrected_CATS1X,corrected_CATS1Y);
//       //   hXf->Fill(corrected_CATS1X);
//       //   hYf->Fill(corrected_CATS1Y);
   
//     double Z_f = X_f * tan(target_angle * Pi / 180.0);

//     double X_f_corr = shifted_CATS2X
//                     + ((corrected_CATS1X - shifted_CATS2X) / deltaZ_CATS21)
//                     * (deltaZ_targetFromCATS2 + Z_f);

//     double Y_f_corr = shifted_CATS2Y
//                     + ((corrected_CATS1Y - shifted_CATS2Y) / deltaZ_CATS21)
//                     * (deltaZ_targetFromCATS2 + Z_f);

//     hTargetXY->Fill(X_f_corr, Y_f_corr);
//     hXf->Fill(X_f_corr);
//     hYf->Fill(Y_f_corr);

//     hTargetXY_simple->Fill(X_f, Y_f);
//     hXf_simple->Fill(X_f);
//     hYf_simple->Fill(Y_f);

//     hC1_centroid->Fill(centroid_CATS1X, centroid_CATS1Y);
//     hC1_shifted->Fill(shifted_CATS1X, shifted_CATS1Y);
//     hC1_rotated->Fill(corrected_CATS1X, corrected_CATS1Y);

//     hC2_centroid->Fill(centroid_CATS2X, centroid_CATS2Y);
//     hC2_shifted->Fill(shifted_CATS2X, shifted_CATS2Y);

//     hTargetXY2->Fill(centroid_CATS1X);
//     hTargetXY3->Fill(centroid_CATS1Y);
//     hTargetXY4->Fill(centroid_CATS2X);
//     hTargetXY5->Fill(centroid_CATS2Y);

//     //std:cout<<centroid_CATS2X<<" "<<centroid_CATS2Y;
// }

// TCanvas* canvasTargetXY = new TCanvas("canvasTargetXY", "Target (c2,c2)", 900, 800);
// hTargetXY->Draw("COLZ");
// canvasTargetXY->Update();

// TCanvas* canvasXf = new TCanvas("canvasXf", "Xf (corr)", 900, 800);
// hXf->Draw("HIST");
// canvasXf->Update();

// TCanvas* canvasYf = new TCanvas("canvasYf", "Yf (corr)", 900, 800);
// hYf->Draw("HIST");
// canvasYf->Update();

// // TCanvas* canvasTargetXY_simple = new TCanvas("canvasTargetXY_simple", "Target simple (X_f,Y_f)", 900, 800);
// // hTargetXY_simple->Draw("COLZ");
// // canvasTargetXY_simple->Update();

// // TCanvas* canvasXf_simple = new TCanvas("canvasXf_simple", "Xf (simple)", 900, 800);
// // hXf_simple->Draw("HIST");
// // canvasXf_simple->Update();

// // TCanvas* canvasYf_simple = new TCanvas("canvasYf_simple", "Yf (simple)", 900, 800);
// // hYf_simple->Draw("HIST");
// // canvasYf_simple->Update();

// TCanvas* canvasC1X = new TCanvas("canvasC1X", "C1X", 900, 800);
// hTargetXY2->Draw();
// canvasC1X->Update();

// TCanvas* canvasC1Y = new TCanvas("canvasC1Y", "C1Y", 900, 800);
// hTargetXY3->Draw();
// canvasC1Y->Update();

// TCanvas* canvasC2X = new TCanvas("canvasC2X", "C2X", 900, 800);
// hTargetXY4->Draw();
// canvasC2X->Update();

// TCanvas* canvasC2Y = new TCanvas("canvasC2Y", "C2Y", 900, 800);
// hTargetXY5->Draw();
// canvasC2Y->Update();

// TCanvas* canvasC1_centroid = new TCanvas("canvasC1_centroid", "CATS1 centroid (2D)", 900, 800);
// hC1_centroid->Draw("COLZ");
// canvasC1_centroid->Update();

// TCanvas* canvasC1_shifted = new TCanvas("canvasC1_shifted", "CATS1 shifted (2D)", 900, 800);
// hC1_shifted->Draw("COLZ");
// canvasC1_shifted->Update();

// TCanvas* canvasC1_rotated = new TCanvas("canvasC1_rotated", "CATS1 rotated (2D)", 900, 800);
// hC1_rotated->Draw("COLZ");
// canvasC1_rotated->Update();

// TCanvas* canvasC2_centroid = new TCanvas("canvasC2_centroid", "CATS2 centroid (2D)", 900, 800);
// hC2_centroid->Draw("COLZ");
// canvasC2_centroid->Update();

// TCanvas* canvasC2_shifted = new TCanvas("canvasC2_shifted", "CATS2 shifted (2D)", 900, 800);
// hC2_shifted->Draw("COLZ");
// canvasC2_shifted->Update();

// gSystem->ProcessEvents();

// TFile* outFile = new TFile(Form("Target_Reconstruction_CATS_remerged_%s.root", runTag.Data()), "RECREATE");
// outFile->cd();

// hTargetXY->SetName("h_XfYf");
// hTargetXY->SetTitle("XfYf;X_{f} (mm);Y_{f} (mm)");
// hTargetXY->Write();

// hTargetXY_simple->SetName("h_XfYf_simple");
// hTargetXY_simple->SetTitle("XfYf simple;X_{f} (mm);Y_{f} (mm)");
// hTargetXY_simple->Write();

// hXf->SetName("h_Xf");
// hXf->SetTitle("Xf;X_{f} (mm);Counts");
// hXf->Write();

// hYf->SetName("h_Yf");
// hYf->SetTitle("Yf;Y_{f} (mm);Counts");
// hYf->Write();

// hXf_simple->SetName("h_Xf_simple");
// hXf_simple->SetTitle("Xf simple;X_{f} (mm);Counts");
// hXf_simple->Write();

// hYf_simple->SetName("h_Yf_simple");
// hYf_simple->SetTitle("Yf simple;Y_{f} (mm);Counts");
// hYf_simple->Write();

// hC1_centroid->Write();
// hC1_shifted->Write();
// hC1_rotated->Write();
// hC2_centroid->Write();
// hC2_shifted->Write();

// hTargetXY2->SetName("h_C1X");
// hTargetXY2->SetTitle("C1X;CATS1 X (mm);Counts");
// hTargetXY2->Write();

// hTargetXY3->SetName("h_C1Y");
// hTargetXY3->SetTitle("C1Y;CATS1 Y (mm);Counts");
// hTargetXY3->Write();

// hTargetXY4->SetName("h_C2X");
// hTargetXY4->SetTitle("C2X;CATS2 X (mm);Counts");
// hTargetXY4->Write();

// hTargetXY5->SetName("h_C2Y");
// hTargetXY5->SetTitle("C2Y;CATS2 Y (mm);Counts");
// hTargetXY5->Write();

// outFile->Close();

// std::cout<<"FINALIZAT";
// std::cin.get();












}
