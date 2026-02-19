#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TROOT.h>

static TH1* GetH1(TFile* filePointer, const char* histogramName)
{
    TH1* histogramPointer = (TH1*)filePointer->Get(histogramName);
    TH1* clonedHistogramPointer = (TH1*)histogramPointer->Clone(Form("%s_%s", histogramName, filePointer->GetName()));
    clonedHistogramPointer->SetDirectory(nullptr);
    return clonedHistogramPointer;
}

static TH2* GetH2(TFile* filePointer, const char* histogramName)
{
    TH2* histogramPointer = (TH2*)filePointer->Get(histogramName);
    TH2* clonedHistogramPointer = (TH2*)histogramPointer->Clone(Form("%s_%s", histogramName, filePointer->GetName()));
    clonedHistogramPointer->SetDirectory(nullptr);
    return clonedHistogramPointer;
}

void Histogram_target_recon_display()
{
    TString filePath1 = "/home/olivia/Desktop/scripts/CATS/Target_Reconstruction_CATS_r0990.root";
    TString filePath2 = "/home/olivia/Desktop/scripts/CATS/Target_Reconstruction_CATS_r0948.root";
    TString filePath3 = "/home/olivia/Desktop/scripts/CATS/Target_Reconstruction_CATS_r1159.root";

    TString label1 = "run 990";
    TString label2 = "run 948";
    TString label3 = "run 1159";

    TFile* filePointer1 = TFile::Open(filePath1, "READ");
    TFile* filePointer2 = TFile::Open(filePath2, "READ");
    TFile* filePointer3 = TFile::Open(filePath3, "READ");

    TH2* hXfYf_1 = GetH2(filePointer1, "h_XfYf");
    TH2* hXfYf_2 = GetH2(filePointer2, "h_XfYf");
    TH2* hXfYf_3 = GetH2(filePointer3, "h_XfYf");

    TH1* hC1X_1 = GetH1(filePointer1, "h_C1X");
    TH1* hC1X_2 = GetH1(filePointer2, "h_C1X");
    TH1* hC1X_3 = GetH1(filePointer3, "h_C1X");

    TH1* hC1Y_1 = GetH1(filePointer1, "h_C1Y");
    TH1* hC1Y_2 = GetH1(filePointer2, "h_C1Y");
    TH1* hC1Y_3 = GetH1(filePointer3, "h_C1Y");

    TH2* hC1_rotated_1 = GetH2(filePointer1, "hC1_rotated");
    TH2* hC1_rotated_2 = GetH2(filePointer2, "hC1_rotated");
    TH2* hC1_rotated_3 = GetH2(filePointer3, "hC1_rotated");

    TH2* hC2_shifted_1 = GetH2(filePointer1, "hC2_shifted");
    TH2* hC2_shifted_2 = GetH2(filePointer2, "hC2_shifted");
    TH2* hC2_shifted_3 = GetH2(filePointer3, "hC2_shifted");

    TCanvas* canvasXf = new TCanvas("canvasXf", "XfYf (3 runs)", 1800, 600);
    canvasXf->Divide(3, 1);
    canvasXf->cd(1); hXfYf_1->SetTitle(Form("XfYf - %s;X_{f} (mm);Y_{f} (mm)", label1.Data())); hXfYf_1->Draw("COLZ");
    canvasXf->cd(2); hXfYf_2->SetTitle(Form("XfYf - %s;X_{f} (mm);Y_{f} (mm)", label2.Data())); hXfYf_2->Draw("COLZ");
    canvasXf->cd(3); hXfYf_3->SetTitle(Form("XfYf - %s;X_{f} (mm);Y_{f} (mm)", label3.Data())); hXfYf_3->Draw("COLZ");
    canvasXf->Update();

    TCanvas* canvasC1_1D = new TCanvas("canvasC1_1D", "CATS1 centroids (overlay)", 1400, 600);
    canvasC1_1D->Divide(2, 1);

    canvasC1_1D->cd(1);
    hC1X_1->SetLineColor(kRed); hC1X_2->SetLineColor(kBlue); hC1X_3->SetLineColor(kGreen+2);
    hC1X_1->SetLineWidth(2); hC1X_2->SetLineWidth(2); hC1X_3->SetLineWidth(2);
    hC1X_1->SetTitle("CATS1 X centroid;X (mm);Counts"); hC1X_1->Draw("HIST"); hC1X_2->Draw("HIST SAME"); hC1X_3->Draw("HIST SAME");
    TLegend* legendC1X = new TLegend(0.65, 0.70, 0.88, 0.88);
    legendC1X->AddEntry(hC1X_1, label1, "l"); legendC1X->AddEntry(hC1X_2, label2, "l"); legendC1X->AddEntry(hC1X_3, label3, "l");
    legendC1X->Draw();

    canvasC1_1D->cd(2);
    hC1Y_1->SetLineColor(kRed); hC1Y_2->SetLineColor(kBlue); hC1Y_3->SetLineColor(kGreen+2);
    hC1Y_1->SetLineWidth(2); hC1Y_2->SetLineWidth(2); hC1Y_3->SetLineWidth(2);
    hC1Y_1->SetTitle("CATS1 Y centroid;Y (mm);Counts"); hC1Y_1->Draw("HIST"); hC1Y_2->Draw("HIST SAME"); hC1Y_3->Draw("HIST SAME");
    TLegend* legendC1Y = new TLegend(0.65, 0.70, 0.88, 0.88);
    legendC1Y->AddEntry(hC1Y_1, label1, "l"); legendC1Y->AddEntry(hC1Y_2, label2, "l"); legendC1Y->AddEntry(hC1Y_3, label3, "l");
    legendC1Y->Draw();

    canvasC1_1D->Update();

    TCanvas* canvasC1_2D = new TCanvas("canvasC1_2D", "CATS1 2D (rotated) (3 runs)", 1800, 600);
    canvasC1_2D->Divide(3, 1);
    canvasC1_2D->cd(1); hC1_rotated_1->SetTitle(Form("CATS1 rotated - %s;X (mm);Y (mm)", label1.Data())); hC1_rotated_1->Draw("COLZ");
    canvasC1_2D->cd(2); hC1_rotated_2->SetTitle(Form("CATS1 rotated - %s;X (mm);Y (mm)", label2.Data())); hC1_rotated_2->Draw("COLZ");
    canvasC1_2D->cd(3); hC1_rotated_3->SetTitle(Form("CATS1 rotated - %s;X (mm);Y (mm)", label3.Data())); hC1_rotated_3->Draw("COLZ");
    canvasC1_2D->Update();

    TCanvas* canvasC2_2D = new TCanvas("canvasC2_2D", "CATS2 2D (shifted) (3 runs)", 1800, 600);
    canvasC2_2D->Divide(3, 1);
    canvasC2_2D->cd(1); hC2_shifted_1->SetTitle(Form("CATS2 shifted - %s;X (mm);Y (mm)", label1.Data())); hC2_shifted_1->Draw("COLZ");
    canvasC2_2D->cd(2); hC2_shifted_2->SetTitle(Form("CATS2 shifted - %s;X (mm);Y (mm)", label2.Data())); hC2_shifted_2->Draw("COLZ");
    canvasC2_2D->cd(3); hC2_shifted_3->SetTitle(Form("CATS2 shifted - %s;X (mm);Y (mm)", label3.Data())); hC2_shifted_3->Draw("COLZ");
    canvasC2_2D->Update();

    gSystem->ProcessEvents();
    std::cin.get();

    filePointer1->Close();
    filePointer2->Close();
    filePointer3->Close();
}
