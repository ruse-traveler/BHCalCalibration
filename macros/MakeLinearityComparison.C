// ----------------------------------------------------------------------------
// 'MakeLinearityComparison.C'
// Derek Anderson
// 03.12.2023
//
// Use this quickly plot the calculated
// linearities from 'DoHCalCalibration.C'.
// ----------------------------------------------------------------------------

#include <iostream>
#include "TH2.h"
#include "TPad.h"
#include "TFile.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGraphErrors.h"

using namespace std;

// global constants
static const UInt_t NLine(2);
static const UInt_t NPlot(2);
static const UInt_t NPad(2);
static const UInt_t NVtx(4);
static const UInt_t NTxt(2);



void MakeLinearityComparison() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning linearity comparison plot-maker..." << endl;

  // output and denominator parameters
  const TString sOutput("lineComparison.sciGlassVsImaging.d12m3y2023.root");
  const TString sLine[NLine]      = {"forSciGlassReso.application_forLinearityAnd2dPlots.e2t20th35145n5KeaPim.d12m3y2023.root",
                                     "forImagingReso.application_forLinearityAnd2dPlots.e2t20th35145n5KeaPim.d12m3y2023.root"};
  const TString sHistLine[NLine]  = {"grLineEne_LD",
                                     "grLineEne_LD"};
  const TString sNameLine[NLine]  = {"grSciGlassLine",
                                     "grImagingLine"};

  // plot parameters
  const UInt_t  nFrameX(51);
  const UInt_t  nFrameY(102);
  const TString sOptLine[NLine]    = {"LP", "LP"};
  const Float_t xyFrameRange[NVtx] = {-1., -1., 50., 50.};
  const Float_t xyPlotRange[NVtx]  = {0.,  0.,  35., 25.};

  // style parameters
  const TString sTitle("");
  const TString sTitleX("E_{par} [GeV]");
  const TString sTitleY("Linearity");
  const UInt_t  fColRes[NLine] = {634, 602};
  const UInt_t  fMarRes[NLine] = {22,  23};

  // text parameters
  const TString sHeader("");
  const TString sTxt[NTxt] = {"ePIC simulation [23.01.0]", "single #pi^{-}"};
  const TString sLabelLine[NLine] = {"Full detector (SciGlass)",
                                     "Full detector (Imaging)"};

  // open output file
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  if (!fOutput) {
    cerr << "PANIC: couldn't open output file!\n" << endl;
    return;
  }

  // open linearity file
  TFile *fLine[NLine];
  for (UInt_t iLine = 0; iLine < NLine; iLine++) {
    fLine[iLine] = new TFile(sLine[iLine].Data(), "read");
    if (!fLine[iLine]) {
      cerr << "PANIC: couldn't open linearity file #" << iLine << "!" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab input graphs
  TGraphErrors *grLine[NLine];
  for (UInt_t iLine = 0; iLine < NLine; iLine++) {
    grLine[iLine] = (TGraphErrors*) fLine[iLine] -> Get(sHistLine[iLine]);
    if (!grLine[iLine]) {
      cerr << "PANIC: couldn't grab linearity graph #" << iLine << "!" << endl;
      return;
    }
    grLine[iLine] -> SetName(sNameLine[iLine].Data());
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
  for (UInt_t iLine = 0; iLine < NLine; iLine++) {
    grLine[iLine] -> SetMarkerColor(fColRes[iLine]);
    grLine[iLine] -> SetMarkerStyle(fMarRes[iLine]);
    grLine[iLine] -> SetFillColor(fColRes[iLine]);
    grLine[iLine] -> SetFillStyle(fFil);
    grLine[iLine] -> SetLineColor(fColRes[iLine]);
    grLine[iLine] -> SetLineStyle(fLin);
    grLine[iLine] -> SetLineWidth(fWid);
    grLine[iLine] -> SetTitle(sTitle.Data());
    grLine[iLine] -> GetXaxis() -> SetRangeUser(xyPlotRange[0], xyPlotRange[2]);
    grLine[iLine] -> GetXaxis() -> SetTitle(sTitleX.Data());
    grLine[iLine] -> GetXaxis() -> SetTitleFont(fTxt);
    grLine[iLine] -> GetXaxis() -> SetTitleSize(fTit);
    grLine[iLine] -> GetXaxis() -> SetTitleOffset(fOffX);
    grLine[iLine] -> GetXaxis() -> SetLabelFont(fTxt);
    grLine[iLine] -> GetXaxis() -> SetLabelSize(fLab);
    grLine[iLine] -> GetXaxis() -> CenterTitle(fCnt);
    grLine[iLine] -> GetYaxis() -> SetRangeUser(xyPlotRange[1], xyPlotRange[3]);
    grLine[iLine] -> GetYaxis() -> SetTitle(sTitleY.Data());
    grLine[iLine] -> GetYaxis() -> SetTitleFont(fTxt);
    grLine[iLine] -> GetYaxis() -> SetTitleSize(fTit);
    grLine[iLine] -> GetYaxis() -> SetTitleOffset(fOffY);
    grLine[iLine] -> GetYaxis() -> SetLabelFont(fTxt);
    grLine[iLine] -> GetYaxis() -> SetLabelSize(fLab);
    grLine[iLine] -> GetYaxis() -> CenterTitle(fCnt);
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
  const Float_t hObjLeg      = NLine * 0.05;
  const Float_t yObjLeg      = 0.1 + hObjLeg;
  const Float_t fLegXY[NVtx] = {0.1, 0.1, 0.3, yObjLeg};

  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], sHeader.Data());
  leg -> SetFillColor(fColLeg);
  leg -> SetFillStyle(fFilLeg);
  leg -> SetLineColor(fColLeg);
  leg -> SetLineStyle(fLinLeg);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  for (UInt_t iLine = 0; iLine < NLine; iLine++) {
    leg -> AddEntry(grLine[iLine], sLabelLine[iLine], "p");
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
  for(UInt_t iLine = 0; iLine < NLine; iLine++) {
    grLine[iLine] -> Draw(sOptLine[iLine].Data());
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
  for (UInt_t iLine = 0; iLine < NLine; iLine++) {
    grLine[iLine] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  for (UInt_t iLine = 0; iLine < NLine; iLine++) {
    fLine[iLine] -> cd();
    fLine[iLine] -> Close();
  }
  cout << "  Finished plot!\n" << endl;

}

// end ------------------------------------------------------------------------
