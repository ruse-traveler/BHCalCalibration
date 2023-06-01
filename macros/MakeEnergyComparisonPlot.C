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
static const UInt_t NParBins(4);
static const UInt_t NInputs(3);
static const UInt_t NPlot(2);
static const UInt_t NPad(2);
static const UInt_t NVtx(4);
static const UInt_t NTxt(3);



void MakeEnergyComparisonPlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning energy comparsion plot macro..." << endl;

  // i/o parameters
  const TString sOutput("uncalibEnergyComparison.forTileVsTowerCheck_varyTileEneMin.e220th45pim.d1m6y2023.root");
  const TString sInput[NInputs] = {
    "../ecal_study/calibration_output/mar/forImagingReso.training_noNClustAndWithNHits_withGraphicUpdate.e2t20th35145n5KeaPim.d9m3y2023.root",
    "tileVsTowerCalibCheck.tileTraining_emin06ecen6.e220th35145n30Kpim.d1m6y2023.tmva.root",
    "tileVsTowerCalibCheck.tileTraining_emin3ecen30.e220th35145n30Kpim.d1m6y2023.tmva.root"
  };
  const TString sHistEne[NInputs][NParBins] = {
    {"Resolution/hHCalEne_ene2", "Resolution/hHCalEne_ene5", "Resolution/hHCalEne_ene10", "Resolution/hHCalEne_ene20"},
    {"Resolution/hHCalEne_ene2", "Resolution/hHCalEne_ene5", "Resolution/hHCalEne_ene10", "Resolution/hHCalEne_ene20"},
    {"Resolution/hHCalEne_ene2", "Resolution/hHCalEne_ene5", "Resolution/hHCalEne_ene10", "Resolution/hHCalEne_ene20"},
  };
  const TString sNameEne[NInputs][NParBins] = {
    {"hRawTowerClustEne_ePar2",            "hRawTowerClustEne_ePar5",            "hRawTowerClustEne_ePar10",            "hRawTowerClustEne_ePar20"},
    {"hRawTileClustEne_emin06ecen6_ePar2", "hRawTileClustEne_emin06ecen6_ePar5", "hRawTileClustEne_emin06ecen6_ePar10", "hRawTileClustEne_emin06ecen6_epar20"},
    {"hRawTileClustEne_emin3ecen30_ePar2", "hRawTileClustEne_emin3ecen30_ePar5", "hRawTileClustEne_emin3ecen30_ePar10", "hRawTileClustEne_emin3ecen30_epar20"}
  };

  // plot parameters
  const TString sOptEne[NInputs][NParBins] = {
    {"",     "same", "same", "same"},
    {"same", "same", "same", "same"},
    {"same", "same", "same", "same"}
  };
  const Float_t xPlotRange[NPlot] = {-1., 40.};

  // style parameters
  const TString sTitle("");
  const TString sTitleX("E_{clust}^{reco} [GeV]");
  const TString sTitleY("arbitrary units");
  const UInt_t  fColEne[NInputs][NParBins] = {
    {803, 893, 883, 863},
    {809, 899, 889, 869},
    {806, 896, 886, 866}
  };
  const UInt_t  fMarEne[NInputs][NParBins] = {
    {20, 20, 20, 20},
    {26, 26, 26, 26},
    {32, 32, 32, 32}
  };

  // legend parameters
  const UInt_t  fColIns[NInputs]  = {923, 922, 921};
  const UInt_t  fMarIns[NInputs]  = {20,  26,  32};
  const UInt_t  fColPar[NParBins] = {809, 909, 889, 869};
  const UInt_t  fMarPar[NParBins] = {20,  20,  20,  20};
  const TString sNameIns[NInputs] = {
    "hLegTowerClusters",
    "hLegTileClusters_emin06ecen6",
    "hLegTileClusters_emin3ecen30"
  };
  const TString sNamePar[NParBins] = {
    "hLegEnePar2",
    "hLegEnePar5",
    "hLegEnePar10",
    "hLegEnePar20"
  };
  const TString sInputLabels[NInputs] = {
    "Tower clusters",
    "Tile clusters: E_{min} = 0.6 MeV, E_{min}^{cent} = 6 MeV",
    "Tile clusters: E_{min} = 3 MeV, E_{min}^{cent} = 30 MeV"
  };
  const TString sParBinLabels[NParBins] = {
    "E_{par} = 2 GeV",
    "E_{par} = 5 GeV",
    "E_{par} = 10 GeV",
    "E_{par} = 20 GeV"
  };

  // text parameters
  const TString sText[NTxt] = {
    "#bf{ePIC} Simulation [23.05.0]",
    "single #pi^{-}, #theta #in (45^{#circ}, 145^{#circ})",
    "#bf{Imaging configuration}"
  };

  // norm/rebin parameters
  const Bool_t doIntNorm[NInputs][NParBins] = {
    {false, false, false, false},
    {false, false, false, false},
    {false, false, false, false}
  };
  const Bool_t doRebin[NInputs][NParBins] = {
    {false, false, false, false},
    {false, false, false, false},
    {false, false, false, false}
  };
  const UInt_t nRebin[NInputs][NParBins] = {
    {2, 2, 2, 2},
    {2, 2, 2, 2},
    {2, 2, 2, 2}
  };

  // open files
  TFile *fOutput = new TFile(sOutput.Data(),     "recreate");
  if (!fOutput) {
    cerr << "PANIC: couldn't open output file!\n" << endl;
    return;
  }

  TFile *fInput[NInputs];
  for (UInt_t iInput = 0; iInput < NInputs; iInput++) {
    fInput[iInput] = new TFile(sInput[iInput].Data(),  "read");
    if (!fInput[iInput]) {
      cerr << "PANIC: couldn't open input file #" << iInput << "!\n" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hEne[NInputs][NParBins];
  for (UInt_t iInput = 0; iInput < NInputs; iInput++) {
    for (UInt_t iParBin = 0; iParBin < NParBins; iParBin++) {
      hEne[iInput][iParBin] = (TH1D*) fInput[iInput] -> Get(sHistEne[iInput][iParBin].Data());
      if (!hEne[iInput][iParBin]) {
        cerr << "PANIC: couldn't grab input histogram # " << iParBin << " from input file #" << iInput << "!\n" << endl;
        return;
      }
      hEne[iInput][iParBin] -> SetName(sNameEne[iInput][iParBin].Data());
    }
  }
  cout << "    Grabbed histograms." << endl;

  // rebin histograms
  Bool_t didRebin(false);
  for (UInt_t iInput = 0; iInput < NInputs; iInput++) {
    for (UInt_t iParBin = 0; iParBin < NParBins; iParBin++) {
      if (doRebin[iInput][iParBin]) {
        hEne[iInput][iParBin] -> Rebin(nRebin[iInput][iParBin]);
        didRebin = true;
      }
    }
  }
  if (didRebin) cout << "    Rebinned histograms." << endl;

  // normalize by integrals (if needed)
  Bool_t didIntNorm(false);
  for (UInt_t iInput = 0; iInput < NInputs; iInput++) {
    for (UInt_t iParBin = 0; iParBin < NParBins; iParBin++) {
      const Double_t intEne = hEne[iInput][iParBin] -> Integral();
      hEne[iInput][iParBin] -> Scale(1. / intEne);
    }
  }
  if (didIntNorm)  cout << "    Normalized histograms by integral." << endl;

  // set styles
  const UInt_t  fFil(0);
  const UInt_t  fLin(1);
  const UInt_t  fWid(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fLab(0.04);
  const Float_t fTit(0.04);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.3);
  for (UInt_t iInput = 0; iInput < NInputs; iInput++) {
    for (UInt_t iParBin = 0; iParBin < NParBins; iParBin++) {
      hEne[iInput][iParBin] -> SetMarkerColor(fColEne[iInput][iParBin]);
      hEne[iInput][iParBin] -> SetMarkerStyle(fMarEne[iInput][iParBin]);
      hEne[iInput][iParBin] -> SetFillColor(fColEne[iInput][iParBin]);
      hEne[iInput][iParBin] -> SetFillStyle(fFil);
      hEne[iInput][iParBin] -> SetLineColor(fColEne[iInput][iParBin]);
      hEne[iInput][iParBin] -> SetLineStyle(fLin);
      hEne[iInput][iParBin] -> SetLineWidth(fWid);
      hEne[iInput][iParBin] -> SetTitle(sTitle.Data());
      hEne[iInput][iParBin] -> SetTitleFont(fTxt);
      hEne[iInput][iParBin] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
      hEne[iInput][iParBin] -> GetXaxis() -> SetTitle(sTitleX.Data());
      hEne[iInput][iParBin] -> GetXaxis() -> SetTitleFont(fTxt);
      hEne[iInput][iParBin] -> GetXaxis() -> SetTitleSize(fTit);
      hEne[iInput][iParBin] -> GetXaxis() -> SetTitleOffset(fOffX);
      hEne[iInput][iParBin] -> GetXaxis() -> SetLabelFont(fTxt);
      hEne[iInput][iParBin] -> GetXaxis() -> SetLabelSize(fLab);
      hEne[iInput][iParBin] -> GetXaxis() -> CenterTitle(fCnt);
      hEne[iInput][iParBin] -> GetYaxis() -> SetTitle(sTitleY.Data());
      hEne[iInput][iParBin] -> GetYaxis() -> SetTitleFont(fTxt);
      hEne[iInput][iParBin] -> GetYaxis() -> SetTitleSize(fTit);
      hEne[iInput][iParBin] -> GetYaxis() -> SetTitleOffset(fOffY);
      hEne[iInput][iParBin] -> GetYaxis() -> SetLabelFont(fTxt);
      hEne[iInput][iParBin] -> GetYaxis() -> SetLabelSize(fLab);
      hEne[iInput][iParBin] -> GetYaxis() -> CenterTitle(fCnt);
    }
  }

  // create histograms for legend
  const UInt_t fFilHist(0);
  const UInt_t fLinHist(1);

  TH1D *hLegIns[NInputs];
  TH1D *hLegPar[NParBins];
  for (UInt_t iInput = 0; iInput < NInputs; iInput++) {
    hLegIns[iInput] = (TH1D*) hEne[iInput][0] -> Clone();
    hLegIns[iInput] -> SetName(sNameIns[iInput].Data());
    hLegIns[iInput] -> SetMarkerColor(fColIns[iInput]);
    hLegIns[iInput] -> SetMarkerStyle(fMarIns[iInput]);
    hLegIns[iInput] -> SetFillColor(fColIns[iInput]);
    hLegIns[iInput] -> SetFillStyle(fFilHist);
    hLegIns[iInput] -> SetLineColor(fColIns[iInput]);
    hLegIns[iInput] -> SetLineStyle(fLinHist);
  }
  for (UInt_t iParBin = 0; iParBin < NParBins; iParBin++) {
    hLegPar[iParBin] = (TH1D*) hEne[0][iParBin] -> Clone();
    hLegPar[iParBin] -> SetName(sNamePar[iParBin].Data());
    hLegPar[iParBin] -> SetMarkerColor(fColPar[iParBin]);
    hLegPar[iParBin] -> SetMarkerStyle(fMarPar[iParBin]);
    hLegPar[iParBin] -> SetFillColor(fColPar[iParBin]);
    hLegPar[iParBin] -> SetFillStyle(fFilHist);
    hLegPar[iParBin] -> SetLineColor(fColPar[iParBin]);
    hLegPar[iParBin] -> SetLineStyle(fLinHist);
  }
  cout << "    Set styles." << endl;

  // determine the number of columns
  const Bool_t isNInMoreThanNPar = (NInputs >= NParBins);

  UInt_t nColumn(1);
  if (isNInMoreThanNPar) {
    nColumn = NInputs;
  } else {
    nColumn = NParBins;
  }

  // make legend
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const UInt_t  nObjLe(2);
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
  if (isNInMoreThanNPar) {
    for (UInt_t iInput = 0; iInput < NInputs; iInput++) {
      leg -> AddEntry(hLegIns[iInput], sInputLabels[iInput].Data(), "pf");
    }
    for (UInt_t iParBin = 0; iParBin < NParBins; iParBin++) {
      leg -> AddEntry(hLegPar[iParBin], sParBinLabels[iParBin].Data(), "pf");
    }
  } else {
    for (UInt_t iParBin = 0; iParBin < NParBins; iParBin++) {
      leg -> AddEntry(hLegPar[iParBin], sParBinLabels[iParBin].Data(), "pf");
    }
    for (UInt_t iInput = 0; iInput < NInputs; iInput++) {
      leg -> AddEntry(hLegIns[iInput], sInputLabels[iInput].Data(), "pf");
    }
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

  // make plot
  const UInt_t  width(750);
  const UInt_t  height(750);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY(1);
  const UInt_t  fFrame(0);
  const Float_t fMarginL(0.15);
  const Float_t fMarginR(0.02);
  const Float_t fMarginT(0.005);
  const Float_t fMarginB(0.15);

  TCanvas *cPlot = new TCanvas("cPlot", "", width, height);
  cPlot   -> SetGrid(fGrid, fGrid);
  cPlot   -> SetTicks(fTick, fTick);
  cPlot   -> SetLogx(fLogX);
  cPlot   -> SetLogy(fLogY);
  cPlot   -> SetBorderMode(fMode);
  cPlot   -> SetBorderSize(fBord);
  cPlot   -> SetFrameBorderMode(fFrame);
  cPlot   -> SetLeftMargin(fMarginL);
  cPlot   -> SetRightMargin(fMarginR);
  cPlot   -> SetTopMargin(fMarginT);
  cPlot   -> SetBottomMargin(fMarginB);
  for (UInt_t iInput = 0; iInput < NInputs; iInput++) {
    for (UInt_t iParBin = 0; iParBin < NParBins; iParBin++) {
      hEne[iInput][iParBin] -> Draw(sOptEne[iInput][iParBin].Data());
    }
  }
  leg     -> Draw();
  txt     -> Draw();
  fOutput -> cd();
  cPlot   -> Write();
  cPlot   -> Close();
  cout << "    Made plot." << endl;

  // save histograms
  fOutput -> cd();
  for (UInt_t iInput = 0; iInput < NInputs; iInput++) {
    for (UInt_t iParBin = 0; iParBin < NParBins; iParBin++) {
      hEne[iInput][iParBin] -> Write();
    }
  }
  for (UInt_t iInput = 0; iInput < NInputs; iInput++) {
    hLegIns[iInput] -> Write();
  }
  for (UInt_t iParBin = 0; iParBin < NParBins; iParBin++) {
    hLegPar[iParBin] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  for (UInt_t iInput = 0; iInput < NInputs; iInput++) {
    fInput[iInput] -> cd();
    fInput[iInput] -> Close();
  }
  cout << "  Finished plot!\n" << endl;

}

// end ------------------------------------------------------------------------
