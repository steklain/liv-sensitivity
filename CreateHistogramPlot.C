#include <iostream>
#include <TCanvas.h>
#include <TGraphErrors.h>

void CreateHistogramPlot() {
    // Binned data with values and corresponding errors
    Double_t binValues[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    Double_t binErrors[] = {0.2, 0.3, 0.1, 0.4, 0.25};
    Int_t numBins = sizeof(binValues) / sizeof(binValues[0]);

    // Create a TGraphErrors object
    TGraphErrors *graph = new TGraphErrors(numBins);

    // Set the values and errors for each bin
    for (Int_t i = 0; i < numBins; ++i) {
        graph->SetPoint(i, i + 0.5, binValues[i]);  // Shift X-position by 0.5 for centering bins
        graph->SetPointError(i, 0.5, binErrors[i]);  // Set symmetric error in X and Y directions
    }

    // Create a canvas
    TCanvas *canvas = new TCanvas("canvas", "Histogram with Error Bars", 800, 600);

    // Set graph options
    graph->SetTitle("Histogram with Error Bars");
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.5);

    // Draw the graph
    graph->Draw("AP");

    // Save the plot as an image file
    canvas->SaveAs("histogram_with_error_bars.png");

    // Clean up
    delete graph;
    delete canvas;
}

int main() {
    CreateHistogramPlot();
    return 0;
}
