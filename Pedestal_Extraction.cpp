#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include <iostream>
using namespace std;

#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TDirectory.h>
#include <TString.h>
#include <iostream>

int Pedestal_Extraction() {

    TChain chain("AD");
    chain.Add("/media/olivia/Partition1/CATS/r0948_*.root");

    TTree* catsTree = &chain;

    Float_t CATS2YV[28];
    UShort_t CATS2YVN[28];

    catsTree->SetBranchAddress("CATS2YV", CATS2YV);
    catsTree->SetBranchAddress("CATS2YVN", CATS2YVN);

    TCanvas* canvas_by_strip[28];
    TH1F* histogram_by_strip[28];

    TString pdf_name = "CATS2YV_vs_CATS2YVN.pdf";

    for (int strip_nr = 1; strip_nr <= 27; strip_nr++) {

        canvas_by_strip[strip_nr] = new TCanvas(Form("canvas_CATS2YVN_%d", strip_nr),
                                                Form("canvas_CATS2YVN_%d", strip_nr),
                                                900, 700);

        catsTree->Draw(Form("CATS2YV>>hist_CATS2YVN_%d(1000,0,2000)", strip_nr),
                       Form("CATS2YVN==%d", strip_nr));

        histogram_by_strip[strip_nr] = (TH1F*)gDirectory->Get(Form("hist_CATS2YVN_%d", strip_nr));
        if (!histogram_by_strip[strip_nr] || histogram_by_strip[strip_nr]->GetEntries() < 50) {
            std::cout << strip_nr << " nan nan" << std::endl;

            // if (strip_nr == 1) canvas_by_strip[strip_nr]->Print(pdf_name + "(");
            // else if (strip_nr == 27) canvas_by_strip[strip_nr]->Print(pdf_name + ")");
            // else canvas_by_strip[strip_nr]->Print(pdf_name);

            continue;
        }

        int maximum_bin = histogram_by_strip[strip_nr]->GetMaximumBin();
        double peak_x = histogram_by_strip[strip_nr]->GetBinCenter(maximum_bin);

        double fit_min = peak_x - 600;
        double fit_max = peak_x + 600;
        if (fit_min < 0.0) fit_min = 0.0;
        if (fit_max > 2000.0) fit_max = 2000.0;

        TF1* gaussian_fit = new TF1(Form("gaus_fit_%d", strip_nr), "gaus", fit_min, fit_max);
        histogram_by_strip[strip_nr]->Fit(gaussian_fit, "QR");

        double mean_value = gaussian_fit->GetParameter(1);
        double sigma_value = gaussian_fit->GetParameter(2);

        std::cout << strip_nr << " " << mean_value << " " << sigma_value << std::endl;

        canvas_by_strip[strip_nr]->SetLogy(1);
        canvas_by_strip[strip_nr]->Update();

        if (strip_nr == 1) canvas_by_strip[strip_nr]->Print(pdf_name + "(");
        else if (strip_nr == 27) canvas_by_strip[strip_nr]->Print(pdf_name + ")");
        else canvas_by_strip[strip_nr]->Print(pdf_name);
    }

    return 0;
}
//TODO e posibil ca ptr cats 2 sa nu fie factorul 3.4 ci mai mult i guess