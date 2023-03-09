// ----------------------------------------------------------------------------
// 'MakeResolutionComparison.C'
// Derek Anderson
// 04.18.2021
//
// Use this quickly plot the calculated
// resolutions from 'DoHCalCalibration.C'
// and 'TMVARegressionApplication.C'.
// ----------------------------------------------------------------------------

#include <iostream>
#include "TH2.h"
#include "TPad.h"
#include "TFile.h"
#include "TLine.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGraphErrors.h"

using namespace std;

// global constants
static const UInt_t NReso(3);
static const UInt_t NPlot(2);
static const UInt_t NPad(2);
static const UInt_t NVtx(4);
static const UInt_t NTxt(2);



void MakeResolutionComparison() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning resolution comparison plot-maker..." << endl;

  // output and denominator parameters
  const TString sOutput("resoComparison.testBeamVsSciGlassVsImage.d8m3y2023.root");
  const TString sReso[NReso]      = {"calibration_output/forTestBeamReso.e1t20th35145n5KeaPim.d8m3y2023.root",
                                     "calibration_output/forSciGlassReso.application.e2t20th35145n5KeaPim.d8m3y2023.root",
                                     "calibration_output/forImagingReso.application.e2t20th35145n5KeaPim.d8m3y2023.root"};
  const TString sHistReso[NReso]  = {"Resolution/grResoEne",
                                     "grResoEne_LD",
                                     "grResoEne_LD"};
  const TString sNameReso[NReso]  = {"grTestBeamReso",
                                     "grSciGlassReso",
                                     "grImagingReso"};

  // plot parameters
  const UInt_t  nFrameX(51);
  const UInt_t  nFrameY(102);
  const TString sOptReso[NReso]    = {"LP", "LP", "LP"};
  const Float_t xyFrameRange[NVtx] = {-1., -1., 50., 50.};
  const Float_t xyPlotRange[NVtx]  = {0.,  0.,  35., 1.2};

  // style parameters
  const TString sTitle("");
  const TString sTitleX("E_{par} [GeV]");
  const TString sTitleY("Resolution (#sigma_{E} / <E_{reco}>)");
  const UInt_t  fColRes[NReso] = {923, 634, 602};
  const UInt_t  fMarRes[NReso] = {20,  22,  23};

  // text parameters
  const TString sHeader("");
  const TString sTxt[NTxt] = {"ePIC simulation [23.01.0]", "single #pi^{-}"};
  const TString sLabelReso[NReso] = {"Only BHCal",
                                     "Full detector (SciGlass)",
                                     "Full detector (Imaging)"};

  // open output file
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  if (!fOutput) {
    cerr << "PANIC: couldn't open output file!\n" << endl;
    return;
  }

  // open resolution file
  TFile *fReso[NReso];
  for (UInt_t iReso = 0; iReso < NReso; iReso++) {
    fReso[iReso] = new TFile(sReso[iReso].Data(), "read");
    if (!fReso[iReso]) {
      cerr << "PANIC: couldn't open resolution file #" << iReso << "!" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab input graphs
  TGraphErrors *grReso[NReso];
  for (UInt_t iReso = 0; iReso < NReso; iReso++) {
    grReso[iReso] = (TGraphErrors*) fReso[iReso] -> Get(sHistReso[iReso]);
    if (!grReso[iReso]) {
      cerr << "PANIC: couldn't grab resolution graph #" << iReso << "!" << endl;
      return;
    }
    grReso[iReso] -> SetName(sNameReso[iReso].Data());
  }
  cout << "    Grabbed graphs." << endl;

  // set styles
  const UInt_t  fFil(0);
  const UInt_t  fLin(1);
  const UInt_t  fWid(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fLab(0.04);
  const Float_t fTit(0.04);
  const Float_t fOffX(1.1);
  const Float_t fOffY(1.3);
  for (UInt_t iReso = 0; iReso < NReso; iReso++) {
    grReso[iReso] -> SetMarkerColor(fColRes[iReso]);
    grReso[iReso] -> SetMarkerStyle(fMarRes[iReso]);
    grReso[iReso] -> SetFillColor(fColRes[iReso]);
    grReso[iReso] -> SetFillStyle(fFil);
    grReso[iReso] -> SetLineColor(fColRes[iReso]);
    grReso[iReso] -> SetLineStyle(fLin);
    grReso[iReso] -> SetLineWidth(fWid);
    grReso[iReso] -> SetTitle(sTitle.Data());
    grReso[iReso] -> GetXaxis() -> SetRangeUser(xyPlotRange[0], xyPlotRange[2]);
    grReso[iReso] -> GetXaxis() -> SetTitle(sTitleX.Data());
    grReso[iReso] -> GetXaxis() -> SetTitleFont(fTxt);
    grReso[iReso] -> GetXaxis() -> SetTitleSize(fTit);
    grReso[iReso] -> GetXaxis() -> SetTitleOffset(fOffX);
    grReso[iReso] -> GetXaxis() -> SetLabelFont(fTxt);
    grReso[iReso] -> GetXaxis() -> SetLabelSize(fLab);
    grReso[iReso] -> GetXaxis() -> CenterTitle(fCnt);
    grReso[iReso] -> GetYaxis() -> SetRangeUser(xyPlotRange[1], xyPlotRange[3]);
    grReso[iReso] -> GetYaxis() -> SetTitle(sTitleY.Data());
    grReso[iReso] -> GetYaxis() -> SetTitleFont(fTxt);
    grReso[iReso] -> GetYaxis() -> SetTitleSize(fTit);
    grReso[iReso] -> GetYaxis() -> SetTitleOffset(fOffY);
    grReso[iReso] -> GetYaxis() -> SetLabelFont(fTxt);
    grReso[iReso] -> GetYaxis() -> SetLabelSize(fLab);
    grReso[iReso] -> GetYaxis() -> CenterTitle(fCnt);
  }

  // make frame histogram
  TH2D *hFrame = new TH2D("hFrame", "", nFrameX, xyFrameRange[0], xyFrameRange[2], nFrameY, xyFrameRange[1], xyFrameRange[3]);
  hFrame -> SetTitle(sTitle.Data());
  hFrame -> SetTitleFont(fTxt);
  hFrame -> GetXaxis() -> SetRangeUser(xyPlotRange[0], xyPlotRange[2]);
  hFrame -> GetXaxis() -> SetTitle(sTitleX.Data());
  hFrame -> GetXaxis() -> SetTitleFont(fTxt);
  hFrame -> GetXaxis() -> SetTitleSize(fTit);
  hFrame -> GetXaxis() -> SetTitleOffset(fOffX);
  hFrame -> GetXaxis() -> SetLabelFont(fTxt);
  hFrame -> GetXaxis() -> SetLabelSize(fLab);
  hFrame -> GetXaxis() -> CenterTitle(fCnt);
  hFrame -> GetYaxis() -> SetRangeUser(xyPlotRange[1], xyPlotRange[3]);
  hFrame -> GetYaxis() -> SetTitle(sTitleY.Data());
  hFrame -> GetYaxis() -> SetTitleFont(fTxt);
  hFrame -> GetYaxis() -> SetTitleSize(fTit);
  hFrame -> GetYaxis() -> SetTitleOffset(fOffY);
  hFrame -> GetYaxis() -> SetLabelFont(fTxt);
  hFrame -> GetYaxis() -> SetLabelSize(fLab);
  hFrame -> GetYaxis() -> CenterTitle(fCnt);
  cout << "    Set styles." << endl;

  // make legend
  const UInt_t  fColLeg      = 0;
  const UInt_t  fFilLeg      = 0;
  const UInt_t  fLinLeg      = 0;
  const Float_t hObjLeg      = NReso * 0.05;
  const Float_t yObjLeg      = 0.1 + hObjLeg;
  const Float_t fLegXY[NVtx] = {0.1, 0.1, 0.3, yObjLeg};

  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], sHeader.Data());
  leg -> SetFillColor(fColLeg);
  leg -> SetFillStyle(fFilLeg);
  leg -> SetLineColor(fColLeg);
  leg -> SetLineStyle(fLinLeg);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  for (UInt_t iReso = 0; iReso < NReso; iReso++) {
    leg -> AddEntry(grReso[iReso], sLabelReso[iReso], "p");
  }
  cout << "    Made legend." << endl;

  // make text
  const UInt_t  fColTxt      = 0;
  const UInt_t  fFilTxt      = 0;
  const UInt_t  fLinTxt      = 0;
  const Float_t hObjTxt      = NTxt * 0.05;
  const Float_t yObjTxt      = 0.1 + hObjTxt;
  const Float_t fTxtXY[NVtx] = {0.3, 0.1, 0.5, yObjTxt};

  TPaveText *txt = new TPaveText(fTxtXY[0], fTxtXY[1], fTxtXY[2], fTxtXY[3], "NDC NB");
  txt -> SetFillColor(fColTxt);
  txt -> SetFillStyle(fFilTxt);
  txt -> SetLineColor(fColTxt);
  txt -> SetLineStyle(fLinTxt);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAln);
  for (UInt_t iTxt = 0; iTxt < NTxt; iTxt++) {
    txt -> AddText(sTxt[iTxt].Data());
  }
  cout << "    Made text." << endl;

  // make plot
  const UInt_t  width(750);
  const UInt_t  height(750);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY(0);
  const UInt_t  fFrame(0);
  const Float_t fMarginL(0.15);
  const Float_t fMarginR(0.02);
  const Float_t fMarginT(0.02);
  const Float_t fMarginB(0.15);

  TCanvas *cPlot = new TCanvas("cPlot", "", width, height);
  cPlot  -> SetGrid(fGrid, fGrid);
  cPlot  -> SetTicks(fTick, fTick);
  cPlot  -> SetBorderMode(fMode);
  cPlot  -> SetBorderSize(fBord);
  cPlot  -> SetFrameBorderMode(fFrame);
  cPlot  -> SetLeftMargin(fMarginL);
  cPlot  -> SetRightMargin(fMarginR);
  cPlot  -> SetTopMargin(fMarginT);
  cPlot  -> SetBottomMargin(fMarginB);
  cPlot  -> SetLogx(fLogX);
  cPlot  -> SetLogy(fLogY);
  cPlot  -> cd();
  hFrame -> Draw();
  for(UInt_t iReso = 0; iReso < NReso; iReso++) {
    grReso[iReso] -> Draw(sOptReso[iReso].Data());
  }
  leg     -> Draw();
  txt     -> Draw();
  fOutput -> cd();
  cPlot   -> Write();
  cPlot   -> Close();
  cout << "    Made plot." << endl;

  // save histograms
  fOutput -> cd();
  hFrame  -> Write();
  for (UInt_t iReso = 0; iReso < NReso; iReso++) {
    grReso[iReso] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  for (UInt_t iReso = 0; iReso < NReso; iReso++) {
    fReso[iReso] -> cd();
    fReso[iReso] -> Close();
  }
  cout << "  Finished plot!\n" << endl;

}

// end ------------------------------------------------------------------------
