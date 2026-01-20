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

std::string input_file = "/media/olivia/Partition1/CATS/r0193_000a.root";
std::string output_file_name = "coeficienti_calibrare_Cats2_Y.txt";
//std::ofstream output_lin_reg("coeficienti_regresie_Cats2Y.txt"); Nu merge declaratia in afara main


constexpr int number_of_wires = 28;
constexpr int number_of_peaks = 6;
double pulser_values[5]={10,25,50,75,100};
double pedestal[28];


//Salvam histogramele si fiturile
std::vector<TH1F*> stripHistograms;
std::vector<TCanvas*> stripCanvases;
std::array<TF1*, number_of_peaks> gaussianFits;


//Fiind pointer putem salva direct in el valorile schimbate
void Gaussian_Fit(Double_t *PeakPositions, TH1F *stripHistogram, int strip_index){
    for (int peakIndex = 0; peakIndex < 6; ++peakIndex)
    {
        double mean = PeakPositions[peakIndex];

        TF1 *fit = new TF1(
            Form("fit_strip_%d_peak_%d", strip_index, peakIndex),
            "gaus",
            mean - 200.0,
            mean + 200.0
        );

        stripHistogram->Fit(fit, "RQ+");
        PeakPositions[peakIndex] = fit->GetParameter(1);
    }
}

void Linear_Regression(double *pulser_values, double *PeakPositions, int stripIndex, std::ofstream &output_lin_reg){

   double shifted_Peaks[5];
   for(int i=1;i<=5;i++)shifted_Peaks[i-1]=PeakPositions[i]-PeakPositions[0];
   TGraph *graph = new TGraph(5, &shifted_Peaks[0], &pulser_values[0]);

    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.2);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);

    TF1 *line = new TF1(
        Form("lin_strip_%d", stripIndex),
        "pol1",
       shifted_Peaks[0],
       shifted_Peaks[4]
    );

    graph->Fit(line, "Q");

    graph->Draw("AP");      // A = axes, P = points
    line->Draw("same");

    // std::cin.get();

    double intercept = line->GetParameter(0);
    double slope     = line->GetParameter(1);
    //output_lin_reg<<stripIndex<<" "<<slope<<" "<<intercept<<" "<<PeakPositions[0]<<'\n';

}

void CATS_calibration()
{
    //I/O 
    TFile *inputRootFile = TFile::Open(input_file.c_str(), "READ");
    if (!inputRootFile || inputRootFile->IsZombie()) { printf("Cannot open file\n"); return; }
    std::ofstream output_lin_reg("coeficienti_regresie_Cats2Y.txt");

    TTree *catsTree = (TTree*)inputRootFile->Get("AD");
    if (!catsTree) { printf("Cannot find tree AD\n"); return; }

    //std::ofstream output_file(output_file_name);

    const int histogramBins = 4000;
    const double histogramMin = 0.0;
    const double histogramMax = 20000.0; // y este doar pana la 6000 

    const int maxPeaks = 6;
    const double searchSigma = 2.0;
    const double searchThreshold = 0.3; //peakurile noastre au dimensiuni similare
    const double fitHalfWindow = 100.0;

    for (int stripIndex = 1; stripIndex < 27; ++stripIndex) // zero este defect
    {
        //Generam canvasurile si histogramele
        TCanvas *stripCanvas = new TCanvas(Form("c_CATS2Y_strip_%d", stripIndex), Form("CATS2Y strip %d", stripIndex), 900, 650);
        stripCanvases.push_back(stripCanvas);

        catsTree->Draw(
            Form("CATS2YV>>h_CATS2Y_strip_%d(%d,%f,%f)", stripIndex, histogramBins, histogramMin, histogramMax),
            Form("CATS2YVN==%d", stripIndex),
            "goff"
        );

        TH1F *stripHistogram = (TH1F*)gDirectory->Get(Form("h_CATS2Y_strip_%d", stripIndex));
        if (!stripHistogram) continue;
        stripHistograms.push_back(stripHistogram); //salvam histograma pentru export


        //Gasim peakurile si le aranjam crescator
        TSpectrum *peakFinder = new TSpectrum(maxPeaks);
        Int_t numberOfPeaksFound = peakFinder->Search(stripHistogram, searchSigma, "", searchThreshold);
       // cout<<numberOfPeaksFound<<endl;
        Double_t *PeakPositions = peakFinder->GetPositionX();
        std::sort(PeakPositions, PeakPositions + numberOfPeaksFound); // nu este nevoie sa le salvam, salvam direct valorile din fit 

      
        Gaussian_Fit(PeakPositions,stripHistogram,stripIndex);
       // Linear_Regression(pulser_values,PeakPositions,stripIndex,output_lin_reg);
       // Verify_Calibration();

        // output_file << stripIndex << " ";
        // for (int i = 0; i < 6; i++){

        //     std::cout << PeakPositions[i] << " ";
        //     output_file << PeakPositions[i] << " ";
        // }
        // output_file << "\n";
       
     }
        
   
}
