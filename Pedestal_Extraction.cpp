#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include <iostream>
using namespace std;

int Pedestal_Extraction() {

    TChain chain("AD");
    chain.Add("/media/olivia/Partition1/CATS/r1163_*.root");

    TTree* catsTree = &chain;

    Float_t CATS1YV[28];
    UShort_t CATS1YVN[28];

    catsTree->SetBranchAddress("CATS1YV", CATS1YV);
    catsTree->SetBranchAddress("CATS1YVN", CATS1YVN);

    TCanvas* canvas_by_strip[28];
    TH1F* histogram_by_strip[28];

    TString pdf_name = "CATS1YV_vs_CATS1YVN";

 for (int strip_nr = 1; strip_nr <=27 ; strip_nr++) {
    // cout<<"STRIP NR "<<strip_nr<<endl;
    canvas_by_strip[strip_nr] = new TCanvas(Form("canvas_CATS1YVN_%d", strip_nr), Form("canvas_CATS1YVN_%d", strip_nr), 900, 700);
    catsTree->Draw(Form("CATS1YV>>hist_CATS1YVN_%d(1000,0,2000)", strip_nr), Form("CATS1YVN==%d", strip_nr));
    canvas_by_strip[strip_nr]->SetLogy(1);
    canvas_by_strip[strip_nr]->Update();

     if (strip_nr == 0) canvas_by_strip[strip_nr]->Print(pdf_name + "(");
        else if (strip_nr == 27) canvas_by_strip[strip_nr]->Print(pdf_name + ")");
        else canvas_by_strip[strip_nr]->Print(pdf_name);
}


    return 0;
}
