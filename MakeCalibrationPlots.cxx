// ----------------------------------------------------------------------------
// 'MakeCalibrationPlots.C'
// Derek Anderson
// 09.13.2023
//
// Macro to plot histograms from 'TrainAndEvaluateBHCalCalibration.C'
// ----------------------------------------------------------------------------

// c utilties
#include <vector>
#include <string>
#include <iostream>
// root classes
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"

using namespace std;

// global constants
static const UInt_t NTxt(3);
static const UInt_t NVtx(4);
static const UInt_t NHist(4);
static const UInt_t NRange(2);
static const UInt_t NEneBins(10);
static const UInt_t NTmvaVar(28);
static const UInt_t NTmvaSpec(1);



void MakeCalibrationPlots() {

  const uint32_t fColEneBin[NEneBins] = {923, 799, 809, 899, 909, 879, 889, 859, 869, 839};
  const uint32_t fMarEneBin[NEneBins] = {25,  27,  32,  26,  29,  30,  25,  27,  32,  26};
  const string   sHCalEne[NEneBins] = {
    "hHCalEne_ene2",
    "hHCalEne_ene3",
    "hHCalEne_ene4",
    "hHCalEne_ene5",
    "hHCalEne_ene6",
    "hHCalEne_ene8",
    "hHCalEne_ene10",
    "hHCalEne_ene12",
    "hHCalEne_ene16",
    "hHCalEne_ene20"
  };
  const string sHCalDiff[NEneBins] = {
    "hHCalDiff_ene2",
    "hHCalDiff_ene3",
    "hHCalDiff_ene4",
    "hHCalDiff_ene5",
    "hHCalDiff_ene6",
    "hHCalDiff_ene8",
    "hHCalDiff_ene10",
    "hHCalDiff_ene12",
    "hHCalDiff_ene16",
    "hHCalDiff_ene20"
  };
  const string sEneTitleX("E_{lead}^{BHCal} [GeV]");
  const string sDiffTitleX("#DeltaE / E_{par}");
  const string sTitleY("arbitrary units");

  // generic resolution parameters
  const double enePar[NEneBins]    = {2.,  3.,  4.,  5.,  6.,  8,   10.,  12.,  16.,  20.};
  const double eneParMin[NEneBins] = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 9.5,  11.5, 13.5, 18.5};
  const double eneParMax[NEneBins] = {2.5, 3.5, 4.5, 5.5, 6.5, 9.5, 11.5, 13.5, 18.5, 21.5};

  // style parameters
  const UInt_t fFil(0);
  const UInt_t fLin(1);
  const UInt_t fTxt(42);
  const UInt_t fAln(12);
  const UInt_t fCenter(1);
  const float  fOffX(1.2);
  const float  fOffY(1.3);
  const string sTitle("");

  // text parameters
  const string sHeader("#bf{Tower Clusters} (reduced energy param.s)");
  const string sTxt[NTxt] = {
    "#bf{ePIC} simulation [23.05.0]",
    "single #pi^{-}",
    "#bf{Imaging Configuration}"
  };
  const string sLabel[NEneBins] = {
    "E_{par} = 2 GeV",
    "E_{par} = 3 GeV",
    "E_{par} = 4 GeV",
    "E_{par} = 5 GeV",
    "E_{par} = 6 GeV",
    "E_{par} = 8 GeV",
    "E_{par} = 10 GeV",
    "E_{par} = 12 GeV",
    "E_{par} = 16 GeV",
    "E_{par} = 20 GeV"
  };

  // set histogram styles
  for (UInt_t iEneBin = 0; iEneBin < NEneBins; iEneBin++) {
    hHCalEneBin[iEneBin]  -> SetMarkerColor(fColEneBin[iEneBin]);
    hHCalEneBin[iEneBin]  -> SetMarkerStyle(fMarEneBin[iEneBin]);
    hHCalEneBin[iEneBin]  -> SetLineColor(fColEneBin[iEneBin]);
    hHCalEneBin[iEneBin]  -> SetLineStyle(fLin);
    hHCalEneBin[iEneBin]  -> SetFillColor(fColEneBin[iEneBin]);
    hHCalEneBin[iEneBin]  -> SetFillStyle(fFil);
    hHCalEneBin[iEneBin]  -> SetTitle(sTitle.data());
    hHCalEneBin[iEneBin]  -> SetTitleFont(fTxt);
    hHCalEneBin[iEneBin]  -> GetXaxis() -> SetTitle(sEneTitleX.data());
    hHCalEneBin[iEneBin]  -> GetXaxis() -> SetTitleFont(fTxt);
    hHCalEneBin[iEneBin]  -> GetXaxis() -> SetTitleOffset(fOffX);
    hHCalEneBin[iEneBin]  -> GetXaxis() -> CenterTitle(fCenter);
    hHCalEneBin[iEneBin]  -> GetYaxis() -> SetTitle(sTitleY.data());
    hHCalEneBin[iEneBin]  -> GetYaxis() -> SetTitleFont(fTxt);
    hHCalEneBin[iEneBin]  -> GetYaxis() -> SetTitleOffset(fOffY);
    hHCalEneBin[iEneBin]  -> GetYaxis() -> CenterTitle(fCenter);
    hHCalDiffBin[iEneBin] -> SetMarkerColor(fColEneBin[iEneBin]);
    hHCalDiffBin[iEneBin] -> SetMarkerStyle(fMarEneBin[iEneBin]);
    hHCalDiffBin[iEneBin] -> SetLineColor(fColEneBin[iEneBin]);
    hHCalDiffBin[iEneBin] -> SetLineStyle(fLin);
    hHCalDiffBin[iEneBin] -> SetFillColor(fColEneBin[iEneBin]);
    hHCalDiffBin[iEneBin] -> SetFillStyle(fFil);
    hHCalDiffBin[iEneBin] -> SetTitle(sTitle.data());
    hHCalDiffBin[iEneBin] -> SetTitleFont(fTxt);
    hHCalDiffBin[iEneBin] -> GetXaxis() -> SetTitle(sDiffTitleX.data());
    hHCalDiffBin[iEneBin] -> GetXaxis() -> SetTitleFont(fTxt);
    hHCalDiffBin[iEneBin] -> GetXaxis() -> SetTitleOffset(fOffX);
    hHCalDiffBin[iEneBin] -> GetXaxis() -> CenterTitle(fCenter);
    hHCalDiffBin[iEneBin] -> GetYaxis() -> SetTitle(sTitleY.data());
    hHCalDiffBin[iEneBin] -> GetYaxis() -> SetTitleFont(fTxt);
    hHCalDiffBin[iEneBin] -> GetYaxis() -> SetTitleOffset(fOffY);
    hHCalDiffBin[iEneBin] -> GetYaxis() -> CenterTitle(fCenter);
  }
  cout << "    Set styles of resolution histograms." << endl;

  // make legend
  const UInt_t  fColLeg      = 0;
  const UInt_t  fFilLeg      = 0;
  const UInt_t  fLinLeg      = 0;
  const float hObjLeg      = NEneBins * 0.05;
  const float yObjLeg      = 0.1 + hObjLeg;
  const float fLegXY[NVtx] = {0.1, 0.1, 0.3, yObjLeg};

  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], sHeader.data());
  leg -> SetFillColor(fColLeg);
  leg -> SetFillStyle(fFilLeg);
  leg -> SetLineColor(fColLeg);
  leg -> SetLineStyle(fLinLeg);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  for (UInt_t iEneBin = 0; iEneBin < NEneBins; iEneBin++) {
    leg -> AddEntry(hHCalEneBin[iEneBin], sLabel[iEneBin].data(), "pf");
  }
  cout << "    Made legend." << endl;

  // make text
  const UInt_t fColTxt      = 0;
  const UInt_t fFilTxt      = 0;
  const UInt_t fLinTxt      = 0;
  const float  hObjTxt      = NTxt * 0.05;
  const float  yObjTxt      = 0.1 + hObjTxt;
  const float  fTxtXY[NVtx] = {0.3, 0.1, 0.5, yObjTxt};

  TPaveText *txt = new TPaveText(fTxtXY[0], fTxtXY[1], fTxtXY[2], fTxtXY[3], "NDC NB");
  txt -> SetFillColor(fColTxt);
  txt -> SetFillStyle(fFilTxt);
  txt -> SetLineColor(fColTxt);
  txt -> SetLineStyle(fLinTxt);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAln);
  for (UInt_t iTxt = 0; iTxt < NTxt; iTxt++) {
    txt -> AddText(sTxt[iTxt].data());
  }
  cout << "    Made text." << endl;

  // plot fit distributions
  const UInt_t  width(750);
  const UInt_t  height(750);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY(1);
  const UInt_t  fFrame(0);
  const float fMarginL(0.15);
  const float fMarginR(0.02);
  const float fMarginT(0.02);
  const float fMarginB(0.15);

  TCanvas *cResoEne = new TCanvas("cResoEne", "", width, height);
  cResoEne       -> SetGrid(fGrid, fGrid);
  cResoEne       -> SetTicks(fTick, fTick);
  cResoEne       -> SetBorderMode(fMode);
  cResoEne       -> SetBorderSize(fBord);
  cResoEne       -> SetFrameBorderMode(fFrame);
  cResoEne       -> SetLeftMargin(fMarginL);
  cResoEne       -> SetRightMargin(fMarginR);
  cResoEne       -> SetTopMargin(fMarginT);
  cResoEne       -> SetBottomMargin(fMarginB);
  cResoEne       -> SetLogx(fLogX);
  cResoEne       -> SetLogy(fLogY);
  cResoEne       -> cd();
  hHCalEneBin[0] -> Draw();
  for (UInt_t iEneBin = 1; iEneBin < NEneBins; iEneBin++) {
    hHCalEneBin[iEneBin] -> Draw("same");
  }
  leg      -> Draw();
  txt      -> Draw();
  fOutput  -> cd();
  cResoEne -> Write();
  cResoEne -> Close();

  TCanvas *cResoDiff = new TCanvas("cResoDiff", "", width, height);
  cResoDiff       -> SetGrid(fGrid, fGrid);
  cResoDiff       -> SetTicks(fTick, fTick);
  cResoDiff       -> SetBorderMode(fMode);
  cResoDiff       -> SetBorderSize(fBord);
  cResoDiff       -> SetFrameBorderMode(fFrame);
  cResoDiff       -> SetLeftMargin(fMarginL);
  cResoDiff       -> SetRightMargin(fMarginR);
  cResoDiff       -> SetTopMargin(fMarginT);
  cResoDiff       -> SetBottomMargin(fMarginB);
  cResoDiff       -> SetLogx(fLogX);
  cResoDiff       -> SetLogy(fLogY);
  cResoDiff       -> cd();
  cResoDiff       -> cd();
  hHCalDiffBin[0] -> Draw();
  for (UInt_t iEneBin = 1; iEneBin < NEneBins; iEneBin++) {
    hHCalDiffBin[iEneBin] -> Draw("same");
  }
  leg       -> Draw();
  leg       -> Draw();
  fOutput   -> cd();
  cResoDiff -> Write();
  cResoDiff -> Close();
  cout << "    Made resolution plots." << endl;

}

// end ------------------------------------------------------------------------
