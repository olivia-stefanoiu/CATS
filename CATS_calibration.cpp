#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TSpectrum.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TGraph.h>
#include <TString.h>

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <TSystem.h>


void Read_Reference_Values(const char *file_name, double reference_values[28])
{
    std::ifstream input_stream(file_name);
    if (!input_stream.is_open()) {
        std::cout << "Cannot open reference file: " << file_name << '\n';
        return;
    }

    int strip_index;
    double slope_value;
    double intercept_value;
    double last_column_value;

    while (input_stream >> strip_index >> slope_value >> intercept_value >> last_column_value) {
        if (strip_index >= 0 && strip_index < 28) {
            reference_values[strip_index] = last_column_value;
        }
    }

    input_stream.close();
}

void Gaussian_Fit(Double_t *PeakPositions, TH1F *stripHistogram, int strip_index, int numberOfPeaksFound)
{
    for (int peak_index = 0; peak_index < numberOfPeaksFound; ++peak_index)
    {
        double mean_value = PeakPositions[peak_index];

        TF1 *fit = new TF1(
            Form("fit_strip_%d_peak_%d", strip_index, peak_index),
            "gaus",
            mean_value - 200.0,
            mean_value + 200.0
        );

        stripHistogram->Fit(fit, "RQ+");
        PeakPositions[peak_index] = fit->GetParameter(1);
    }
}
void Linear_Regression(double pulser_values[6], double *PeakPositions, int stripIndex, double reference_values[28])
{
    double shifted_Peaks[6];

    for (int peak_index = 0; peak_index <= 5; ++peak_index) {
        shifted_Peaks[peak_index] = PeakPositions[peak_index] - reference_values[stripIndex];
    }

    TCanvas *canvas = new TCanvas(
        Form("c_lin_strip_%d", stripIndex),
        Form("Linear regression strip %d", stripIndex),
        800,
        600
    );

    TGraph *graph = new TGraph(6, shifted_Peaks, pulser_values);
    graph->SetTitle(Form("Strip %d;Shifted peak position;Pulser value", stripIndex));
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.2);

    TF1 *line = new TF1(
        Form("lin_strip_%d", stripIndex),
        "pol1",
        shifted_Peaks[0],
        shifted_Peaks[5]
    );

    graph->Draw("AP");
    graph->Fit(line, "Q");
    line->Draw("same");
    canvas->Modified();
    canvas->Update();
    gSystem->ProcessEvents();

    double intercept = line->GetParameter(0);
    double slope = line->GetParameter(1);

    std::cout << stripIndex << " " << slope << " " << intercept << " " << reference_values[stripIndex] << '\n';
    std::cout << "Press Enter for next strip..." << std::endl;
    std::cin.get();
}
void CATS_calibration()
{
    const char *input_file = "/media/olivia/Partition1/CATS/r0193_000a.root";
    const char *reference_file = "/home/olivia/Desktop/scripts/CATS/coeficienti_regresie_CATS2Y.txt";

    double pulser_values[6] = {0,10,25,50,75,100};
    double reference_values[28] = {0};

    Read_Reference_Values(reference_file, reference_values);

    TFile *inputRootFile = TFile::Open(input_file, "READ");
    if (!inputRootFile || inputRootFile->IsZombie()) {
        printf("Cannot open file\n");
        return;
    }

    TTree *catsTree = (TTree*)inputRootFile->Get("AD");
    if (!catsTree) {
        printf("Cannot find tree AD\n");
        inputRootFile->Close();
        return;
    }

    const int histogramBins = 4000;
    const double histogramMin = 0.0;
    const double histogramMax = 20000.0;

    const int maxPeaks = 6;
    const double searchSigma = 2.0;
    const double searchThreshold = 0.3;

    for (int stripIndex = 1; stripIndex <= 27; ++stripIndex)
    {
        gDirectory->Delete(Form("h_CATS2Y_strip_%d;*", stripIndex));

        catsTree->Draw(
            Form("CATS2YV>>h_CATS2Y_strip_%d(%d,%f,%f)", stripIndex, histogramBins, histogramMin, histogramMax),
            Form("CATS2YVN==%d", stripIndex),
            "goff"
        );

        TH1F *stripHistogram = (TH1F*)gDirectory->Get(Form("h_CATS2Y_strip_%d", stripIndex));
        if (!stripHistogram) {
            std::cout << "Strip " << stripIndex << " histogram not found\n";
            continue;
        }

        TSpectrum *peakFinder = new TSpectrum(maxPeaks);
        Int_t numberOfPeaksFound = peakFinder->Search(stripHistogram, searchSigma, "", searchThreshold);

        if (numberOfPeaksFound < 6) {
            std::cout << "Strip " << stripIndex << " has only " << numberOfPeaksFound << " peaks found\n";
            delete peakFinder;
            continue;
        }

        Double_t *PeakPositions = peakFinder->GetPositionX();
        std::sort(PeakPositions, PeakPositions + numberOfPeaksFound);

        Gaussian_Fit(PeakPositions, stripHistogram, stripIndex, numberOfPeaksFound);
        std::sort(PeakPositions, PeakPositions + numberOfPeaksFound);

        Linear_Regression(pulser_values, PeakPositions, stripIndex, reference_values);

        delete peakFinder;
    }

    gSystem->ProcessEvents();
std::cin.get();
    inputRootFile->Close();
}