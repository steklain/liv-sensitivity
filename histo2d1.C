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


void histo2d1(){

    TGraph2D *T13DCP = new TGraph2D("BSMchi2_NP.dat");
    TH2 *HistoT13DCP1 = T13DCP->GetHistogram();
    TH2 *HistoT13DCP2 = T13DCP->GetHistogram();
    TH2 *HistoT13DCP3 = T13DCP->GetHistogram();

    gStyle->SetPadLeftMargin(0.23);
    gStyle->SetPadBottomMargin(0.23);
    gStyle->SetPadRightMargin(0.03);
    gStyle->SetPadTopMargin(0.03);
    gStyle->SetLineWidth(3);
    gStyle->SetTextFont(132);

    HistoT13DCP1->SetLineWidth(4);
    //HistoT13DCP1->SetMaximum(7.0);
    HistoT13DCP1->SetTitle("");
    HistoT13DCP1->GetXaxis()->CenterTitle();
    HistoT13DCP1->GetYaxis()->CenterTitle();
    HistoT13DCP1->GetXaxis()->SetTitle("log sin^{2} (2#theta_{13})");
    HistoT13DCP1->GetYaxis()->SetTitle("#delta_{CP} (degrees)");
    HistoT13DCP1->GetXaxis()->SetTitleOffset(1.5);
    HistoT13DCP1->GetYaxis()->SetTitleOffset(1.7);
    HistoT13DCP1->GetXaxis()->SetTitleSize(0.05);
    HistoT13DCP1->GetYaxis()->SetTitleSize(0.05);
    HistoT13DCP1->GetXaxis()->SetLabelSize(0.05);
    HistoT13DCP1->GetYaxis()->SetLabelSize(0.05);
    HistoT13DCP1->GetXaxis()->SetTickLength(0.04);
    HistoT13DCP1->GetYaxis()->SetTickLength(0.04);
    //HistoT13DCP1->SetNdivisions(1005, "X"); 
    //HistoT13DCP1->SetNdivisions(1005, "Y"); 
    //HistoT13DCP1->SetNdivisions(1005, "Z"); 
    //HistoT13DCP1->GetXaxis()->SetDecimals();
    //HistoT13DCP1->GetYaxis()->SetDecimals();
    //HistoT13DCP1->GetXaxis()->SetNoExponent();
    //HistoT13DCP1->GetYaxis()->SetNoExponent();

    //HistoT13DCP1->SetAxisRange(-4, -2.2);
    //HistoT13DCP1->SetAxisRange(150, 350, "y");
    
    
    Double_t contours[1];
     //contours[0]=2.3;
     //contours[1]=4.61;
     contours[0]=5.99;
    HistoT13DCP1->SetContour(1, contours);

    Int_t palette[1];
     palette[0] = 2;
     //palette[1] = 3;
     //palette[2] = 4;
    gStyle->SetPalette(1,palette);


 TCanvas *plot3 = new TCanvas("plot3","plot3",600,500);
    plot3->SetTickx(1);
    plot3->SetTicky(1);

    HistoT13DCP1->Draw("CONT1");
    plot3->Update();

   auto legend = new TLegend(0.75,0.85,0.95,0.95);
   legend->SetHeader("#chi^{2}=2.3","C"); // option "C" allows to center the header
  
   legend->Draw();

    //gr->SetMarkerColor(kBlue);
    //gr->SetMarkerStyle(kFullCircle);

//    gr1->SetMarkerColor(kRed);
//    gr1->SetMarkerStyle(kFullCircle);

   // gr->Draw("P");
    
  //  gr1->Draw("P");

    
//    HistoT13DCP1->SetLineColorAlpha(kWhite, 0.00);
//    HistoT13DCP1->Draw("SAME CONT2 LIST");
//    plot3->SaveAs("T13DCP1.png");


    HistoT13DCP1->SetLineWidth(1);
    //HistoT13DCP1->GetXaxis()->SetRangeUser(-4,-2);

    //TCanvas *plot4 = new TCanvas("plot4","plot4") ;  
    //HistoT13DCP1->Draw("surf1");


}