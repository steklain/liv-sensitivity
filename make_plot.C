#include "TH1.h"
#include "TH2.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH2S.h"


void make_plot(){

    gROOT->SetBatch(true);

    TCanvas *canvas = new TCanvas("canvas", "Chi2", 800, 600);

    auto grnp  = new TGraph("BSMchi2_NP.dat");
    //auto grstat  = new TGraph("BSMchi2_statisticsonly.dat");
    //auto grsys  = new TGraph("BSMchi2_systematics.dat");

    //gr->GetXaxis()->SetRangeUser(0.25,3.0);
    grnp->GetYaxis()->SetRangeUser(0,10);

    grnp->SetTitle("#chi^{2} for Statistics only and systematics;|a_{#mu e}| (10^{-23}GeV);#chi^{2}");
    grnp->GetXaxis()->CenterTitle();
    grnp->GetYaxis()->CenterTitle();
    grnp->GetYaxis()->SetNdivisions(505);

    //grstat->SetLineColor(kRed);
    //grsys->SetLineColor(kBlue);
    grnp->SetLineColor(kBlue);
    //grnp->SetLineStyle(7);

    //grstat->Draw();
    //grsys->Draw("SAME");
    grnp->Draw();

    auto legend = new TLegend(0.7,0.7,0.9,0.9);
    //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    //legend->AddEntry(grstat,"Statistics only","l");
    //legend->AddEntry(grsys,"Systematics","l");
    legend->AddEntry(grnp,"Projection","l");
   
    legend->Draw();

    canvas->SaveAs("chi2_NP.png");
}