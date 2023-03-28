// ----------------------------------------------------------------------------
// 'MakeEnergyComparisonPlot.C'
// Derek Anderson
// 03.12.2023
//
// Use this quickly plot the calibrated energies
// from 'TMVARegressionApplication.C' (or
// otherwise) from both configurations on the
// same canvas.
// ----------------------------------------------------------------------------

#include <iostream>
#include "TH1.h"
#include "TPad.h"
#include "TFile.h"
#include "TLine.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;

// global constants
static const UInt_t NHist(4);
static const UInt_t NPlot(2);
static const UInt_t NPad(2);
static const UInt_t NVtx(4);
static const UInt_t NTxt(2);



void MakeEnergyComparisonPlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning ernergy comparsion plot macro..." << endl;

  // file parameters
  const TString sOutput("calibEnergyComparison.sciGlassVsImaging.d12m3y2023.root");
  const TString sInSciGlass("forSciGlassReso.application_forLinearityAnd2dPlots.e2t20th35145n5KeaPim.d12m3y2023.root");
  const TString sInImaging("forImagingReso.application_forLinearityAnd2dPlots.e2t20th35145n5KeaPim.d12m3y2023.root");

  // denominator parameters
  const TString sLabelSciGlass("Full detector (SciGlass)");
  const TString sHistSciGlass[NHist]  = {"hHCalEne_ene2_LD",
                                         "hHCalEne_ene5_LD",
                                         "hHCalEne_ene10_LD",
                                         "hHCalEne_ene20_LD"};
  const TString sNameSciGlass[NHist]  = {"hSciGlassEne_ePar2",
                                         "hSciGlassEne_ePar5",
                                         "hSciGlassEne_ePar10",
                                         "hSciGlassEne_ePar2"};

  // numerator parameters
  const TString sLabelImaging("Full detector (Imaging)");
  const TString sHistImaging[NHist]  = {"hHCalEne_ene2_LD",
                                        "hHCalEne_ene5_LD",
                                        "hHCalEne_ene10_LD",
                                        "hHCalEne_ene20_LD"};
  const TString sNameImaging[NHist]  = {"hImagingEne_ePar2",
                                        "hImagingEne_ePar5",
                                        "hImagingEne_ePar10",
                                        "hImagingEne_ePar20"};

  // text parameters
  const TString sTitle("");
  const TString sTitleX("E_{par}^{reco} [GeV]");
  const TString sTitleY("arbitrary units");
  const TString sTitleR("Imaging / SciGlass");
  const TString sText[NTxt]    = {"ePIC simulation [23.01.0]", "single #pi^{-}"};
  const TString sLabels[NHist] = {"E_{par} = 2 GeV", "E_{par} = 5 GeV", "E_{par} = 10 GeV", "E_{par} = 20 GeV"};

  // plot parameters
  const TString sNameRatio[NHist]   = {"hRatio_ePar2",  "hRatio_ePar5",  "hRatio_ePar10",  "hRatio_ePar20"};
  const TString sNameLegend[NHist]  = {"hLegend_ePar2", "hLegend_ePar5", "hLegend_ePar10", "hLegend_ePar20"};
  const TString sOptSciGlass[NHist] = {"",              "same",          "same",           "same"};
  const TString sOptImaging[NHist]  = {"same",          "same",          "same",           "same"};
  const TString sOptRatio[NHist]    = {"",              "same",          "same",           "same"};
  const Float_t xPlotRange[NPlot]   = {-1., 40.};
  const UInt_t  fColSci[NHist]      = {809, 909, 889, 869};
  const UInt_t  fColIma[NHist]      = {808, 908, 888, 868};
  const UInt_t  fMarSci[NHist]      = {22,  23,  20,  21};
  const UInt_t  fMarIma[NHist]      = {26,  32,  24,  25};
  const UInt_t  fColLegSci(923);
  const UInt_t  fColLegIma(921);
  const UInt_t  fMarLegSci(20);
  const UInt_t  fMarLegIma(24);

  // norm/rebin parameters
  const Bool_t doIntNorm(false);
  const Bool_t doRebinSciGlass[NHist] = {false, false, false, false};
  const Bool_t doRebinImaging[NHist]  = {false, false, false, false};
  const UInt_t nRebinSciGlass[NHist]  = {2, 2, 2, 2};
  const UInt_t nRebinImaging[NHist]   = {2, 2, 2, 2};

  // open files
  TFile *fOutput   = new TFile(sOutput.Data(),     "recreate");
  TFile *fSciGlass = new TFile(sInSciGlass.Data(), "read");
  TFile *fImaging  = new TFile(sInImaging.Data(),  "read");
  if (!fOutput || !fSciGlass || !fImaging) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fOutput = " << fOutput << ", fSciGlass = " << fSciGlass << ", fImaging = " << fImaging << "\n"
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hSciGlass[NHist];
  TH1D *hImaging[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hSciGlass[iHist] = (TH1D*) fSciGlass -> Get(sHistSciGlass[iHist].Data());
    hImaging[iHist]  = (TH1D*) fImaging  -> Get(sHistImaging[iHist].Data());
    if (!hSciGlass[iHist] || !hImaging[iHist]) {
      cerr << "PANIC: couldn't grab SciGlass or Imaging histogram # " << iHist << "!\n"
           << "       hSciGlass = " << hSciGlass[iHist] << ", hImaging = " << hImaging[iHist] << "\n"
           << endl;
      return;
    }
    hSciGlass[iHist] -> SetName(sNameSciGlass[iHist].Data());
    hImaging[iHist]  -> SetName(sNameImaging[iHist].Data());
  }
  cout << "    Grabbed histograms." << endl;

  // rebin histograms
  Bool_t doRebin(false);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    doRebin = (doRebinSciGlass[iHist] || doRebinImaging[iHist]);
    if (doRebin) break;
  }

  if (doRebin) {
    for (UInt_t iHist = 0; iHist < NHist; iHist++) {
      if (doRebinSciGlass[iHist]) hSciGlass[iHist] -> Rebin(nRebinSciGlass[iHist]);
      if (doRebinImaging[iHist])  hImaging[iHist]  -> Rebin(nRebinImaging[iHist]);
    }
    cout << "    Rebinned histograms." << endl;
  }

  // normalize by integrals )if needed)
  if (doIntNorm) {
    for (UInt_t iHist = 0; iHist < NHist; iHist++) {
      const Double_t intSciGlass = hSciGlass[iHist] -> Integral();
      const Double_t intImaging = hImaging[iHist] -> Integral();
      hSciGlass[iHist] -> Scale(1. / intSciGlass);
      hImaging[iHist]  -> Scale(1. / intImaging);
    }
    cout << "    Normalized histograms by integral." << endl;
  }

  // calculate ratios
  TH1D *hRatio[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hRatio[iHist] = (TH1D*) hSciGlass[iHist] -> Clone();
    hRatio[iHist] -> Reset("ICE");
    hRatio[iHist] -> Divide(hImaging[iHist], hSciGlass[iHist], 1., 1.);
    hRatio[iHist] -> SetName(sNameRatio[iHist]);
  }
  cout << "    Calculated ratios." << endl;

  // set styles
  const UInt_t  fFil(0);
  const UInt_t  fLin(1);
  const UInt_t  fWid(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fLab[NPad]  = {0.074, 0.04};
  const Float_t fTit[NPad]  = {0.074, 0.04};
  const Float_t fOffX[NPad] = {1.1, 1.};
  const Float_t fOffY[NPad] = {0.7, 1.3};
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hSciGlass[iHist] -> SetMarkerColor(fColSci[iHist]);
    hSciGlass[iHist] -> SetMarkerStyle(fMarSci[iHist]);
    hSciGlass[iHist] -> SetFillColor(fColSci[iHist]);
    hSciGlass[iHist] -> SetFillStyle(fFil);
    hSciGlass[iHist] -> SetLineColor(fColSci[iHist]);
    hSciGlass[iHist] -> SetLineStyle(fLin);
    hSciGlass[iHist] -> SetLineWidth(fWid);
    hSciGlass[iHist] -> SetTitle(sTitle.Data());
    hSciGlass[iHist] -> SetTitleFont(fTxt);
    hSciGlass[iHist] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hSciGlass[iHist] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hSciGlass[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    hSciGlass[iHist] -> GetXaxis() -> SetTitleSize(fTit[1]);
    hSciGlass[iHist] -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hSciGlass[iHist] -> GetXaxis() -> SetLabelFont(fTxt);
    hSciGlass[iHist] -> GetXaxis() -> SetLabelSize(fLab[1]);
    hSciGlass[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    hSciGlass[iHist] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hSciGlass[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    hSciGlass[iHist] -> GetYaxis() -> SetTitleSize(fTit[1]);
    hSciGlass[iHist] -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hSciGlass[iHist] -> GetYaxis() -> SetLabelFont(fTxt);
    hSciGlass[iHist] -> GetYaxis() -> SetLabelSize(fLab[1]);
    hSciGlass[iHist] -> GetYaxis() -> CenterTitle(fCnt);
    hImaging[iHist]  -> SetMarkerColor(fColIma[iHist]);
    hImaging[iHist]  -> SetMarkerStyle(fMarIma[iHist]);
    hImaging[iHist]  -> SetFillColor(fColIma[iHist]);
    hImaging[iHist]  -> SetFillStyle(fFil);
    hImaging[iHist]  -> SetLineColor(fColIma[iHist]);
    hImaging[iHist]  -> SetLineStyle(fLin);
    hImaging[iHist]  -> SetLineWidth(fWid);
    hImaging[iHist]  -> SetTitle(sTitle.Data());
    hImaging[iHist]  -> SetTitleFont(fTxt);
    hImaging[iHist]  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hImaging[iHist]  -> GetXaxis() -> SetTitle(sTitleX.Data());
    hImaging[iHist]  -> GetXaxis() -> SetTitleFont(fTxt);
    hImaging[iHist]  -> GetXaxis() -> SetTitleSize(fTit[1]);
    hImaging[iHist]  -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hImaging[iHist]  -> GetXaxis() -> SetLabelFont(fTxt);
    hImaging[iHist]  -> GetXaxis() -> SetLabelSize(fLab[1]);
    hImaging[iHist]  -> GetXaxis() -> CenterTitle(fCnt);
    hImaging[iHist]  -> GetYaxis() -> SetTitle(sTitleY.Data());
    hImaging[iHist]  -> GetYaxis() -> SetTitleFont(fTxt);
    hImaging[iHist]  -> GetYaxis() -> SetTitleSize(fTit[1]);
    hImaging[iHist]  -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hImaging[iHist]  -> GetYaxis() -> SetLabelFont(fTxt);
    hImaging[iHist]  -> GetYaxis() -> SetLabelSize(fLab[1]);
    hImaging[iHist]  -> GetYaxis() -> CenterTitle(fCnt);
    hRatio[iHist]    -> SetMarkerColor(fColIma[iHist]);
    hRatio[iHist]    -> SetMarkerStyle(fMarIma[iHist]);
    hRatio[iHist]    -> SetFillColor(fColIma[iHist]);
    hRatio[iHist]    -> SetFillStyle(fFil);
    hRatio[iHist]    -> SetLineColor(fColIma[iHist]);
    hRatio[iHist]    -> SetLineStyle(fLin);
    hRatio[iHist]    -> SetLineWidth(fWid);
    hRatio[iHist]    -> SetTitle(sTitle.Data());
    hRatio[iHist]    -> SetTitleFont(fTxt);
    hRatio[iHist]    -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hRatio[iHist]    -> GetXaxis() -> SetTitle(sTitleX.Data());
    hRatio[iHist]    -> GetXaxis() -> SetTitleFont(fTxt);
    hRatio[iHist]    -> GetXaxis() -> SetTitleSize(fTit[0]);
    hRatio[iHist]    -> GetXaxis() -> SetTitleOffset(fOffX[0]);
    hRatio[iHist]    -> GetXaxis() -> SetLabelFont(fTxt);
    hRatio[iHist]    -> GetXaxis() -> SetLabelSize(fLab[0]);
    hRatio[iHist]    -> GetXaxis() -> CenterTitle(fCnt);
    hRatio[iHist]    -> GetYaxis() -> SetTitle(sTitleR.Data());
    hRatio[iHist]    -> GetYaxis() -> SetTitleFont(fTxt);
    hRatio[iHist]    -> GetYaxis() -> SetTitleSize(fTit[0]);
    hRatio[iHist]    -> GetYaxis() -> SetTitleOffset(fOffY[0]);
    hRatio[iHist]    -> GetYaxis() -> SetLabelFont(fTxt);
    hRatio[iHist]    -> GetYaxis() -> SetLabelSize(fLab[0]);
    hRatio[iHist]    -> GetYaxis() -> CenterTitle(fCnt);
  }

  // create histograms for legend
  const UInt_t fFilHist(0);
  const UInt_t fLinHist(1);

  TH1D *hLegSci = (TH1D*) hSciGlass[0] -> Clone();
  TH1D *hLegIma = (TH1D*) hImaging[0]  -> Clone();
  hLegSci -> SetName("hLegendSciGlass");
  hLegIma -> SetName("hLegendImaging");
  hLegSci -> SetMarkerColor(fColLegSci);
  hLegSci -> SetMarkerStyle(fMarLegSci);
  hLegSci -> SetFillColor(fColLegSci);
  hLegSci -> SetFillStyle(fFilHist);
  hLegSci -> SetLineColor(fColLegSci);
  hLegSci -> SetLineStyle(fLinHist);
  hLegIma -> SetMarkerColor(fColLegIma);
  hLegIma -> SetMarkerStyle(fMarLegIma);
  hLegIma -> SetFillColor(fColLegIma);
  hLegIma -> SetFillStyle(fFilHist);
  hLegIma -> SetLineColor(fColLegIma);
  hLegIma -> SetLineStyle(fLinHist);

  TH1D *hLegHist[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hLegHist[iHist] = (TH1D*) hSciGlass[iHist] -> Clone();
    hLegHist[iHist] -> SetName(sNameLegend[iHist].Data());
    hLegHist[iHist] -> SetFillStyle(fFilHist);
    hLegHist[iHist] -> SetLineStyle(fLinHist);
  }
  cout << "    Set styles." << endl;

  // make legend
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const UInt_t  nColumn(2);
  const UInt_t  nObjLe((NHist + 2) / nColumn);
  const Float_t hObjLe(nObjLe * 0.05);
  const Float_t yObjLe(0.1 + hObjLe);
  const Float_t fLegXY[NVtx] = {0.1, 0.1, 0.7, yObjLe};

  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3]);
  leg -> SetNColumns(nColumn);
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  leg -> AddEntry(hLegSci, sLabelSciGlass.Data(), "pf");
  leg -> AddEntry(hLegIma, sLabelImaging.Data(),  "pf");
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    leg -> AddEntry(hLegHist[iHist], sLabels[iHist].Data(), "pf");
  }
  cout << "    Made legend." << endl;

  // make text
  const UInt_t  fColTx(0);
  const UInt_t  fFilTx(0);
  const UInt_t  fLinTx(0);
  const UInt_t  nObjTx(NTxt);
  const Float_t hObjTx(nObjTx * 0.05);
  const Float_t yObjTx(0.1 + hObjTx);
  const Float_t fTxtXY[NVtx] = {0.7, 0.1, 0.9, yObjTx};

  TPaveText *txt = new TPaveText(fTxtXY[0], fTxtXY[1], fTxtXY[2], fTxtXY[3], "NDC NB");
  txt -> SetFillColor(fColTx);
  txt -> SetFillStyle(fFilTx);
  txt -> SetLineColor(fColTx);
  txt -> SetLineStyle(fLinTx);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAln);
  for (UInt_t iTxt = 0; iTxt < NTxt; iTxt++) {
    txt -> AddText(sText[iTxt].Data());
  }
  cout << "    Made text." << endl;

  // make line
  const UInt_t  fColLi(923);
  const UInt_t  fLinLi(9);
  const UInt_t  fWidLi(1);
  const Float_t fLinXY[NVtx] = {xPlotRange[0], 1., xPlotRange[1], 1.};

  TLine *line = new TLine(fLinXY[0], fLinXY[1], fLinXY[2], fLinXY[3]);
  line -> SetLineColor(fColLi);
  line -> SetLineStyle(fLinLi);
  line -> SetLineWidth(fWidLi);
  cout << "    Made line." << endl;

  // make plot
  const UInt_t  width(750);
  const UInt_t  height(950);
  const UInt_t  heightNR(750);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY1(0);
  const UInt_t  fLogY2(1);
  const UInt_t  fFrame(0);
  const Float_t fMarginL(0.15);
  const Float_t fMarginR(0.02);
  const Float_t fMarginT1(0.005);
  const Float_t fMarginT2(0.02);
  const Float_t fMarginB1(0.25);
  const Float_t fMarginB2(0.005);
  const Float_t fMarginBNR(0.15);
  const Float_t fPadXY1[NVtx] = {0., 0., 1., 0.35};
  const Float_t fPadXY2[NVtx] = {0., 0.35, 1., 1.};

  TCanvas *cPlot = new TCanvas("cPlot", "", width, height);
  TPad    *pPad1 = new TPad("pPad1", "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
  TPad    *pPad2 = new TPad("pPad2", "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  cPlot     -> SetGrid(fGrid, fGrid);
  cPlot     -> SetTicks(fTick, fTick);
  cPlot     -> SetBorderMode(fMode);
  cPlot     -> SetBorderSize(fBord);
  pPad1     -> SetGrid(fGrid, fGrid);
  pPad1     -> SetTicks(fTick, fTick);
  pPad1     -> SetLogx(fLogX);
  pPad1     -> SetLogy(fLogY1);
  pPad1     -> SetBorderMode(fMode);
  pPad1     -> SetBorderSize(fBord);
  pPad1     -> SetFrameBorderMode(fFrame);
  pPad1     -> SetLeftMargin(fMarginL);
  pPad1     -> SetRightMargin(fMarginR);
  pPad1     -> SetTopMargin(fMarginT1);
  pPad1     -> SetBottomMargin(fMarginB1);
  pPad2     -> SetGrid(fGrid, fGrid);
  pPad2     -> SetTicks(fTick, fTick);
  pPad2     -> SetLogx(fLogX);
  pPad2     -> SetLogy(fLogY2);
  pPad2     -> SetBorderMode(fMode);
  pPad2     -> SetBorderSize(fBord);
  pPad2     -> SetFrameBorderMode(fFrame);
  pPad2     -> SetLeftMargin(fMarginL);
  pPad2     -> SetRightMargin(fMarginR);
  pPad2     -> SetTopMargin(fMarginT2);
  pPad2     -> SetBottomMargin(fMarginB2);
  cPlot     -> cd();
  pPad1     -> Draw();
  pPad2     -> Draw();
  pPad1     -> cd();
  hRatio[0] -> Draw(sOptRatio[0].Data());
  for (UInt_t iHist = 1; iHist < NHist; iHist++) {
    hRatio[iHist] -> Draw(sOptRatio[iHist].Data());
  }
  line         -> Draw();
  pPad2        -> cd();
  hSciGlass[0] -> Draw(sOptSciGlass[0].Data());
  hImaging[0]  -> Draw(sOptImaging[0].Data());
  for(UInt_t iHist = 1; iHist < NHist; iHist++) {
    hSciGlass[iHist] -> Draw(sOptSciGlass[iHist].Data());
    hImaging[iHist]  -> Draw(sOptImaging[iHist].Data());
  }
  leg     -> Draw();
  txt     -> Draw();
  fOutput -> cd();
  cPlot   -> Write();
  cPlot   -> Close();

  TCanvas *cPlotNR = new TCanvas("cPlot_NoRatio", "", width, heightNR);
  cPlotNR      -> SetGrid(fGrid, fGrid);
  cPlotNR      -> SetTicks(fTick, fTick);
  cPlotNR      -> SetLogx(fLogX);
  cPlotNR      -> SetLogy(fLogY2);
  cPlotNR      -> SetBorderMode(fMode);
  cPlotNR      -> SetBorderSize(fBord);
  cPlotNR      -> SetFrameBorderMode(fFrame);
  cPlotNR      -> SetLeftMargin(fMarginL);
  cPlotNR      -> SetRightMargin(fMarginR);
  cPlotNR      -> SetTopMargin(fMarginT1);
  cPlotNR      -> SetBottomMargin(fMarginBNR);
  hSciGlass[0] -> Draw(sOptSciGlass[0].Data());
  hImaging[0]  -> Draw(sOptImaging[0].Data());
  for(UInt_t iHist = 1; iHist < NHist; iHist++) {
    hSciGlass[iHist] -> Draw(sOptSciGlass[iHist].Data());
    hImaging[iHist]  -> Draw(sOptImaging[iHist].Data());
  }
  leg     -> Draw();
  txt     -> Draw();
  fOutput -> cd();
  cPlotNR -> Write();
  cPlotNR -> Close();
  cout << "    Made plot." << endl;

  // save histograms
  fOutput -> cd();
  hLegSci -> Write();
  hLegIma -> Write();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hLegHist[iHist]  -> Write();
    hSciGlass[iHist] -> Write();
    hImaging[iHist]  -> Write();
    hRatio[iHist]    -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOutput   -> cd();
  fOutput   -> Close();
  fSciGlass -> cd();
  fSciGlass -> Close();
  fImaging  -> cd();
  fImaging  -> Close();
  cout << "  Finished plot!\n" << endl;

}

// end ------------------------------------------------------------------------
