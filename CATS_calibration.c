#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <TDirectory.h>
#include <TF1.h>

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <cstdio>

std::string input_file = "/media/olivia/Partition1/CATS/r0193_000a.root";
std::string output_file_name = "coeficienti_calibrare_Cats1_Y.txt";
std::ofstream output_lin_reg("coeficienti_regresie_Cats1Y.txt");


constexpr int number_of_wires = 28;
constexpr int number_of_peaks = 6;
double pulser_values[6]={10,25,50,75,100};
double pedestal[28];


//Salvam histogramele si fiturile
std::vector<TH1F*> stripHistograms(number_of_wires);
std::vector<TCanvas*> stripCanvases(number_of_wires);
std::array<TF1*, number_of_peaks> gaussianFits;



//Fiind pointer putem salva direct in el valorile schimbate
void Gaussian_Fit(Double_t *PeakPositions, TH1F *stripHistogram, int strip_index){
    for (int peakIndex = 0; peakIndex < 6; ++peakIndex)
    {
        double mean = PeakPositions[peakIndex];

        TF1 *fit = new TF1(
            Form("fit_strip_%d_peak_%d", strip_index, strip_index),
            "gaus",
            mean - 200.0,
            mean + 200.0
        );

        stripHistogram->Fit(fit, "RQ+");
        double intercept = line->GetParameter(0);
        double slope     = line->GetParameter(1);

        PeakPositions[peakIndex] = fit->GetParameter(1);
    }
}

void Linear_Regression(double *pulser_values, double *PeakPositions, int stripIndex){

TGraph *graph = new TGraph(5, &pulser_values[1], &PeakPositions[1]);
TF1 *line = new TF1(Form("lin_strip_%d", stripIndex), "pol1", pulser_values[1], pulser_values[5]);

graph->Fit(line, "Q");
graph->Draw("AP");
line->Draw("same");

cin.get()

delete graph;
delete line;

}


void CATS_calibration()
{

    //I/O 
    TFile *inputRootFile = TFile::Open(input_file.c_str(), "READ");
    if (!inputRootFile || inputRootFile->IsZombie()) { printf("Cannot open file\n"); return; }

    TTree *catsTree = (TTree*)inputRootFile->Get("AD");
    if (!catsTree) { printf("Cannot find tree AD\n"); return; }

    std::ofstream output_file(output_file_name);

    const int histogramBins = 2000;
    const double histogramMin = 0.0;
    const double histogramMax = 12000.0; // y este doar pana la 6000 

    const int maxPeaks = 6;
    const double searchSigma = 2.0;
    const double searchThreshold = 0.3; //peakurile noastre au dimensiuni similare
    const double fitHalfWindow = 100.0;


    for (int stripIndex = 0; stripIndex < number_of_wires; ++stripIndex)
    {

        //Generam canvasurile si histogramele
        TCanvas *stripCanvas = new TCanvas(Form("c_CATS1Y_strip_%d", stripIndex), Form("CATS1Y strip %d", stripIndex), 900, 650);
        stripCanvases.push_back(stripCanvas);

        catsTree->Draw(
            Form("CATS1YV>>h_CATS1Y_strip_%d(%d,%f,%f)", stripIndex, histogramBins, histogramMin, histogramMax),
            Form("CATS1YVN==%d", stripIndex),
            "goff"
        );

        TH1F *stripHistogram = (TH1F*)gDirectory->Get(Form("h_CATS1Y_strip_%d", stripIndex));
        if (!stripHistogram) continue;
        stripHistograms.push_back(stripHistogram); //salvam histograma pentru export


        //Gasim peakurile si le aranjam crescator
        TSpectrum *peakFinder = new TSpectrum(maxPeaks);
        Int_t numberOfPeaksFound = peakFinder->Search(stripHistogram, searchSigma, "", searchThreshold);
        cout<<numberOfPeaksFound<<endl;
        Double_t *PeakPositions = peakFinder->GetPositionX();
        std::sort(PeakPositions, PeakPositions + numberOfPeaksFound); // nu este nevoie sa le salvam, salvam direct valorile din fit 

        Gaussian_Fit(PeakPositions,stripHistogram,stripIndex);

  double pulser_values[6]={10,25,50,75,100};
    output_file << stripIndex << " ";
    for (int i = 0; i < 6; i++)
    {
        std::cout << PeakPositions[i] << " ";
        output_file << PeakPositions[i] << " ";
    }
    output_file << "\n";
       
    }
        
    // std::ofstream outputTextFile(output_file.c_str());
    // for (int stripIndex = 0; stripIndex < (int)number_of_wires; ++stripIndex)
    // {
    //     for (int peakIndex = 0; peakIndex < 6; ++peakIndex)
    //     {
    //         if ((int)stripPeakPositions[stripIndex].size() > peakIndex) pulser_values[stripIndex][peakIndex] = stripPeakPositions[stripIndex][peakIndex];
    //         else pulser_values[stripIndex][peakIndex] = 0.0;
    //     }

    //     outputTextFile << "strip " << stripIndex << " pedestal " << pedestal[stripIndex];
    //     for (int peakIndex = 0; peakIndex < 6; ++peakIndex) outputTextFile << " " << pulser_values[stripIndex][peakIndex];
    //     outputTextFile << "\n";
    // }
    // outputTextFile.close();
}
