// ----------------------------------------------------------------------------
// 'MakeResolutionAndLinearityComparison.C'
// Derek Anderson
// 10.19.2023
//
// Use this quickly plot the calculated
// resolutions and linearities from
// 'DoHCalCalibration.C' and
// 'TMVARegressionApplication.C'.
// ----------------------------------------------------------------------------

#include <iostream>
#include "TH2.h"
#include "TPad.h"
#include "TFile.h"
#include "TError.h"
#include "TGraph.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGraphErrors.h"

using namespace std;

// global constants
static const UInt_t NHist(4);
static const UInt_t NPlot(2);
static const UInt_t NTest(7);
static const UInt_t NPad(2);
static const UInt_t NVtx(4);
static const UInt_t NTxt(3);



void MakeResolutionAndLinearityComparison() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning resolution and linearity comparison plot-maker..." << endl;

  // output and denominator parameters
  const TString sOutput("resoComparison_hist.oldVsNewThresholds_noBECalLayers_ddsim.e220th45pim.d19m10y2023.root");
  const TString sInput[NHist] = {
    "tmva_output/forLowTresholdCheck.withDDSim_withTowers_noBECalLayers_emin3ecen30.epic23050image.e220th45n120Kpim.d26m9y2023.tmva.root",
    "tmva_output/forLowTresholdCheck.withDDSim_noBECalLayers_emin3ecen30.epic23050image.e220th45n120Kpim.d26m9y2023.tmva.root",
    "tmva_output/forLowTresholdCheck.withDDSim_noBECalLayers_emin06ecen6.epic23050image.e220th45n120Kpim.d26m9y2023.tmva.root",
    "tmva_output/forLowTresholdCheck.withDDSim_noBECalLayers.epic23080image.e220th45n250Kpim.d5m10y2023.tmva.root"
  };
  const TString sHistReso[NHist] = {
    "resolution/grResoCalibHist_LD",
    "resolution/grResoCalibHist_LD",
    "resolution/grResoCalibHist_LD",
    "resolution/grResoCalibHist_LD"
  };
  const TString sHistLine[NHist] = {
    "resolution/grLineCalibHist_LD",
    "resolution/grLineCalibHist_LD",
    "resolution/grLineCalibHist_LD",
    "resolution/grLineCalibHist_LD"
  };
  const TString sNameReso[NHist] = {
    "grTowerReso_emin3ecen30_hist",
    "grTileReso_emin3ecen30_hist",
    "grTileReso_emin06ecen6_hist",
    "grTileReso_emin05ecen30_hist"
  };
  const TString sNameLine[NHist] = {
    "grTowerLine_emin3ecen30_hist",
    "grTileLine_emin3ecen30_hist",
    "grTileLine_emin06ecen6_hist",
    "grTileLine_emin05ecen30_hist"
  };

  // plot parameters
  const UInt_t  nFrameX(51);
  const UInt_t  nFrameY(102);
  const TString sOptReso[NHist]    = {"LP", "LP", "LP", "LP"};
  const TString sOptLine[NHist]    = {"LP", "LP", "LP", "LP"};
  const Float_t xyFrameRange[NVtx] = {-1., -1., 50., 50.};
  const Float_t xyResoRange[NVtx]  = {0.,  0.,  35., 1.2};
  const Float_t xyLineRange[NVtx]  = {0.,  0.,  23., 23.};

  // style parameters
  const TString sTitle("");
  const TString sTitleX("E_{par} [GeV]");
  const TString sTitleYR("Resolution (#sigma_{E} / #mu(E_{reco}))");
  const TString sTitleYL("Linearity");
  const UInt_t  fColHist[NHist] = {923, 859, 879, 899};
  const UInt_t  fMarHist[NHist] = {20,  26,  32,  24};

  // text parameters
  const TString sHeader("#bf{Hist. Reso.}");
  const TString sTxt[NTxt] = {
    "#bf{ePIC} simulation [23.05.0 vs. 23.08.0]",
    "single #pi^{-}, #theta #in (45^{#circ}, 135^{#circ})",
    "#bf{Imaging configuration}"
  };
  const TString sLabelHist[NHist] = {
    "E_{min} = 3 MeV, E_{min}^{cent} = 30 MeV (tower clust.)",
    "E_{min} = 3 MeV, E_{min}^{cent} = 30 MeV (tile clust.)",
    "E_{min} = 0.6 MeV, E_{min}^{cent} = 6 MeV (tile clust.)",
    "E_{min} = 5 MeV, E_{min}^{cent} = 30 MeV (tile clust.)" 
  };

  // test beam points & parameters
  const Bool_t   addTestBeam(false);
  const UInt_t   fColTest(618);
  const UInt_t   fMarTest(29);
  const TString  sOptTest("LP");
  const TString  sLabelTest("sPHENIX test beam data");
  const TString  sTestRef("[IEEE Transactions on Nuc. Sci., Vol. 65, Iss. 12, pp. 2901-2919, Dec. 2018]");
  const Double_t xValsTestBeam[NTest] = {4.14959877108356, 6.14450880383323, 8.1692122326946,  12.15563223082159, 16.20408511280676, 24.14495469139409, 32.17897143943406};
  const Double_t yValsTestBeam[NTest] = {0.47719893154717, 0.34697739951106, 0.30316859721537, 0.26110700323024,  0.23476189744027,  0.20405296417384,  0.19063440434873};

  // open output file
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  if (!fOutput) {
    cerr << "PANIC: couldn't open output file!\n" << endl;
    return;
  }

  // open resolution file
  TFile *fInput[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fInput[iHist] = new TFile(sInput[iHist].Data(), "read");
    if (!fInput[iHist]) {
      cerr << "PANIC: couldn't open resolution file #" << iHist << "!" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab input graphs
  TGraphErrors *grReso[NHist];
  TGraphErrors *grLine[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    grReso[iHist] = (TGraphErrors*) fInput[iHist] -> Get(sHistReso[iHist]);
    grLine[iHist] = (TGraphErrors*) fInput[iHist] -> Get(sHistLine[iHist]);
    if (!grReso[iHist]) {
      cerr << "PANIC: couldn't grab resolution graph #" << iHist << "!" << endl;
      return;
    }
    if (!grLine[iHist]) {
      cerr << "PANIC: couldn't grab linearity graph #" << iHist << "!" << endl;
      return;
    }
    grReso[iHist] -> SetName(sNameReso[iHist].Data());
    grLine[iHist] -> SetName(sNameLine[iHist].Data());
  }
  cout << "    Grabbed graphs." << endl;

  // create test beam curve
  TGraph *grTest = new TGraph(NTest, xValsTestBeam, yValsTestBeam);
  grTest -> SetName("grFromTestBeamPaper");
  if (addTestBeam) cout << "    Made test beam graph." << endl;

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
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    grReso[iHist] -> SetMarkerColor(fColHist[iHist]);
    grReso[iHist] -> SetMarkerStyle(fMarHist[iHist]);
    grReso[iHist] -> SetFillColor(fColHist[iHist]);
    grReso[iHist] -> SetFillStyle(fFil);
    grReso[iHist] -> SetLineColor(fColHist[iHist]);
    grReso[iHist] -> SetLineStyle(fLin);
    grReso[iHist] -> SetLineWidth(fWid);
    grReso[iHist] -> SetTitle(sTitle.Data());
    grReso[iHist] -> GetXaxis() -> SetRangeUser(xyResoRange[0], xyResoRange[2]);
    grReso[iHist] -> GetXaxis() -> SetTitle(sTitleX.Data());
    grReso[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    grReso[iHist] -> GetXaxis() -> SetTitleSize(fTit);
    grReso[iHist] -> GetXaxis() -> SetTitleOffset(fOffX);
    grReso[iHist] -> GetXaxis() -> SetLabelFont(fTxt);
    grReso[iHist] -> GetXaxis() -> SetLabelSize(fLab);
    grReso[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    grReso[iHist] -> GetYaxis() -> SetRangeUser(xyResoRange[1], xyResoRange[3]);
    grReso[iHist] -> GetYaxis() -> SetTitle(sTitleYR.Data());
    grReso[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    grReso[iHist] -> GetYaxis() -> SetTitleSize(fTit);
    grReso[iHist] -> GetYaxis() -> SetTitleOffset(fOffY);
    grReso[iHist] -> GetYaxis() -> SetLabelFont(fTxt);
    grReso[iHist] -> GetYaxis() -> SetLabelSize(fLab);
    grReso[iHist] -> GetYaxis() -> CenterTitle(fCnt);
    grLine[iHist] -> SetMarkerColor(fColHist[iHist]);
    grLine[iHist] -> SetMarkerStyle(fMarHist[iHist]);
    grLine[iHist] -> SetFillColor(fColHist[iHist]);
    grLine[iHist] -> SetFillStyle(fFil);
    grLine[iHist] -> SetLineColor(fColHist[iHist]);
    grLine[iHist] -> SetLineStyle(fLin);
    grLine[iHist] -> SetLineWidth(fWid);
    grLine[iHist] -> SetTitle(sTitle.Data());
    grLine[iHist] -> GetXaxis() -> SetRangeUser(xyLineRange[0], xyLineRange[2]);
    grLine[iHist] -> GetXaxis() -> SetTitle(sTitleX.Data());
    grLine[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    grLine[iHist] -> GetXaxis() -> SetTitleSize(fTit);
    grLine[iHist] -> GetXaxis() -> SetTitleOffset(fOffX);
    grLine[iHist] -> GetXaxis() -> SetLabelFont(fTxt);
    grLine[iHist] -> GetXaxis() -> SetLabelSize(fLab);
    grLine[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    grLine[iHist] -> GetYaxis() -> SetRangeUser(xyLineRange[1], xyLineRange[3]);
    grLine[iHist] -> GetYaxis() -> SetTitle(sTitleYL.Data());
    grLine[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    grLine[iHist] -> GetYaxis() -> SetTitleSize(fTit);
    grLine[iHist] -> GetYaxis() -> SetTitleOffset(fOffY);
    grLine[iHist] -> GetYaxis() -> SetLabelFont(fTxt);
    grLine[iHist] -> GetYaxis() -> SetLabelSize(fLab);
    grLine[iHist] -> GetYaxis() -> CenterTitle(fCnt);
  }
  grTest -> SetMarkerColor(fColTest);
  grTest -> SetMarkerStyle(fMarTest);
  grTest -> SetFillColor(fColTest);
  grTest -> SetFillStyle(fFil);
  grTest -> SetLineColor(fColTest);
  grTest -> SetLineStyle(fLin);
  grTest -> SetLineWidth(fWid);
  grTest -> SetTitle(sTitle.Data());
  grTest -> GetXaxis() -> SetRangeUser(xyResoRange[0], xyResoRange[2]);
  grTest -> GetXaxis() -> SetTitle(sTitleX.Data());
  grTest -> GetXaxis() -> SetTitleFont(fTxt);
  grTest -> GetXaxis() -> SetTitleSize(fTit);
  grTest -> GetXaxis() -> SetTitleOffset(fOffX);
  grTest -> GetXaxis() -> SetLabelFont(fTxt);
  grTest -> GetXaxis() -> SetLabelSize(fLab);
  grTest -> GetXaxis() -> CenterTitle(fCnt);
  grTest -> GetYaxis() -> SetRangeUser(xyResoRange[1], xyResoRange[3]);
  grTest -> GetYaxis() -> SetTitle(sTitleYR.Data());
  grTest -> GetYaxis() -> SetTitleFont(fTxt);
  grTest -> GetYaxis() -> SetTitleSize(fTit);
  grTest -> GetYaxis() -> SetTitleOffset(fOffY);
  grTest -> GetYaxis() -> SetLabelFont(fTxt);
  grTest -> GetYaxis() -> SetLabelSize(fLab);
  grTest -> GetYaxis() -> CenterTitle(fCnt);

  // make frame histograms
  TH2D *hResoFrame = new TH2D("hResoFrame", "", nFrameX, xyFrameRange[0], xyFrameRange[2], nFrameY, xyFrameRange[1], xyFrameRange[3]);
  hResoFrame -> SetTitle(sTitle.Data());
  hResoFrame -> SetTitleFont(fTxt);
  hResoFrame -> GetXaxis() -> SetRangeUser(xyResoRange[0], xyResoRange[2]);
  hResoFrame -> GetXaxis() -> SetTitle(sTitleX.Data());
  hResoFrame -> GetXaxis() -> SetTitleFont(fTxt);
  hResoFrame -> GetXaxis() -> SetTitleSize(fTit);
  hResoFrame -> GetXaxis() -> SetTitleOffset(fOffX);
  hResoFrame -> GetXaxis() -> SetLabelFont(fTxt);
  hResoFrame -> GetXaxis() -> SetLabelSize(fLab);
  hResoFrame -> GetXaxis() -> CenterTitle(fCnt);
  hResoFrame -> GetYaxis() -> SetRangeUser(xyResoRange[1], xyResoRange[3]);
  hResoFrame -> GetYaxis() -> SetTitle(sTitleYR.Data());
  hResoFrame -> GetYaxis() -> SetTitleFont(fTxt);
  hResoFrame -> GetYaxis() -> SetTitleSize(fTit);
  hResoFrame -> GetYaxis() -> SetTitleOffset(fOffY);
  hResoFrame -> GetYaxis() -> SetLabelFont(fTxt);
  hResoFrame -> GetYaxis() -> SetLabelSize(fLab);
  hResoFrame -> GetYaxis() -> CenterTitle(fCnt);

  TH2D *hLineFrame = new TH2D("hLineFrame", "", nFrameX, xyFrameRange[0], xyFrameRange[2], nFrameY, xyFrameRange[1], xyFrameRange[3]);
  hLineFrame -> SetTitle(sTitle.Data());
  hLineFrame -> SetTitleFont(fTxt);
  hLineFrame -> GetXaxis() -> SetRangeUser(xyLineRange[0], xyLineRange[2]);
  hLineFrame -> GetXaxis() -> SetTitle(sTitleX.Data());
  hLineFrame -> GetXaxis() -> SetTitleFont(fTxt);
  hLineFrame -> GetXaxis() -> SetTitleSize(fTit);
  hLineFrame -> GetXaxis() -> SetTitleOffset(fOffX);
  hLineFrame -> GetXaxis() -> SetLabelFont(fTxt);
  hLineFrame -> GetXaxis() -> SetLabelSize(fLab);
  hLineFrame -> GetXaxis() -> CenterTitle(fCnt);
  hLineFrame -> GetYaxis() -> SetRangeUser(xyLineRange[1], xyLineRange[3]);
  hLineFrame -> GetYaxis() -> SetTitle(sTitleYL.Data());
  hLineFrame -> GetYaxis() -> SetTitleFont(fTxt);
  hLineFrame -> GetYaxis() -> SetTitleSize(fTit);
  hLineFrame -> GetYaxis() -> SetTitleOffset(fOffY);
  hLineFrame -> GetYaxis() -> SetLabelFont(fTxt);
  hLineFrame -> GetYaxis() -> SetLabelSize(fLab);
  hLineFrame -> GetYaxis() -> CenterTitle(fCnt);
  cout << "    Set styles." << endl;

  // make legend
  UInt_t nObjLeg(NHist);
  if (addTestBeam) nObjLeg += 2;

  const UInt_t  fColLeg      = 0;
  const UInt_t  fFilLeg      = 0;
  const UInt_t  fLinLeg      = 0;
  const Float_t hObjLeg      = nObjLeg * 0.05;
  const Float_t yObjLeg      = 0.1 + hObjLeg;
  const Float_t fLegXY[NVtx] = {0.1, 0.1, 0.3, yObjLeg};

  TLegend *legR = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], sHeader.Data());
  legR -> SetFillColor(fColLeg);
  legR -> SetFillStyle(fFilLeg);
  legR -> SetLineColor(fColLeg);
  legR -> SetLineStyle(fLinLeg);
  legR -> SetTextFont(fTxt);
  legR -> SetTextAlign(fAln);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    legR -> AddEntry(grReso[iHist], sLabelHist[iHist].Data(), "fp");
  }
  if (addTestBeam) {
    legR -> AddEntry(grTest,      sLabelTest.Data(), "fp");
    legR -> AddEntry((TObject*)0, sTestRef.Data(),   "");
  }

  TLegend *legL = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], sHeader.Data());
  legL -> SetFillColor(fColLeg);
  legL -> SetFillStyle(fFilLeg);
  legL -> SetLineColor(fColLeg);
  legL -> SetLineStyle(fLinLeg);
  legL -> SetTextFont(fTxt);
  legL -> SetTextAlign(fAln);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    legL -> AddEntry(grLine[iHist], sLabelHist[iHist], "pf");
  }
  cout << "    Made legends." << endl;

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

  TCanvas *cReso = new TCanvas("cReso", "", width, height);
  cReso      -> SetGrid(fGrid, fGrid);
  cReso      -> SetTicks(fTick, fTick);
  cReso      -> SetBorderMode(fMode);
  cReso      -> SetBorderSize(fBord);
  cReso      -> SetFrameBorderMode(fFrame);
  cReso      -> SetLeftMargin(fMarginL);
  cReso      -> SetRightMargin(fMarginR);
  cReso      -> SetTopMargin(fMarginT);
  cReso      -> SetBottomMargin(fMarginB);
  cReso      -> SetLogx(fLogX);
  cReso      -> SetLogy(fLogY);
  cReso      -> cd();
  hResoFrame -> Draw();
  for(UInt_t iHist = 0; iHist < NHist; iHist++) {
    grReso[iHist] -> Draw(sOptReso[iHist].Data());
  }
  if (addTestBeam) grTest -> Draw(sOptTest.Data());
  legR    -> Draw();
  txt     -> Draw();
  fOutput -> cd();
  cReso   -> Write();
  cReso   -> Close();

  TCanvas *cLine = new TCanvas("cLine", "", width, height);
  cLine      -> SetGrid(fGrid, fGrid);
  cLine      -> SetTicks(fTick, fTick);
  cLine      -> SetBorderMode(fMode);
  cLine      -> SetBorderSize(fBord);
  cLine      -> SetFrameBorderMode(fFrame);
  cLine      -> SetLeftMargin(fMarginL);
  cLine      -> SetRightMargin(fMarginR);
  cLine      -> SetTopMargin(fMarginT);
  cLine      -> SetBottomMargin(fMarginB);
  cLine      -> SetLogx(fLogX);
  cLine      -> SetLogy(fLogY);
  cLine      -> cd();
  hLineFrame -> Draw();
  for(UInt_t iHist = 0; iHist < NHist; iHist++) {
    grLine[iHist] -> Draw(sOptLine[iHist].Data());
  }
  legL    -> Draw();
  txt     -> Draw();
  fOutput -> cd();
  cLine   -> Write();
  cLine   -> Close();
  cout << "    Made plots." << endl;

  // save histograms
  fOutput    -> cd();
  hResoFrame -> Write();
  hLineFrame -> Write();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    grReso[iHist] -> Write();
  }
  if (addTestBeam) grTest -> Write();
  cout << "    Saved histograms." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fInput[iHist] -> cd();
    fInput[iHist] -> Close();
  }
  cout << "  Finished plot!\n" << endl;

}

// end ------------------------------------------------------------------------

