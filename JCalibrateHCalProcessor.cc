// 'JCalibrateHCalProcessor.cc'
// Derek Anderson
// 11.02.2022
//
// A simple JANA plugin to compare the
// reconstructed hit and cluster energy
// in the HCal to simulated particles.

// user includes
#include "JCalibrateHCalProcessor.h"
#include <services/rootfile/RootFile_service.h>

// The following just makes this a JANA plugin
extern "C" {
  void InitPlugin(JApplication *app) {
  InitJANAPlugin(app);
  app -> Add(new JCalibrateHCalProcessor);
  }
}




//-------------------------------------------
// InitWithGlobalRootLock
//-------------------------------------------
void JCalibrateHCalProcessor::InitWithGlobalRootLock(){

  // create directory in output file
  auto rootfile_svc = GetApplication() -> GetService<RootFile_service>();
  auto rootfile     = rootfile_svc     -> GetHistFile();
  rootfile -> mkdir("JCalibrateHCal") -> cd();

  // initialize histograms
  const unsigned long nNumBin(200);
  const unsigned long nChrgBin(6);
  const unsigned long nMassBin(1000);
  const unsigned long nEneBin(500);
  const unsigned long nMomBin(500);
  const unsigned long nPosBin(8000);
  const unsigned long nDiffBin(500);
  const unsigned long rNumBin[NRange]  = {0,      200};
  const double        rChrgBin[NRange] = {-3.,    3.};
  const double        rMassBin[NRange] = {0.,     5.};
  const double        rEneBin[NRange]  = {0.,     100.};
  const double        rMomBin[NRange]  = {-50.,   50.};
  const double        rPosBin[NRange]  = {-4000., 4000.};
  const double        rDiffBin[NRange] = {-50.,   50.};
  // particle histograms
  hParChrg           = new TH1D("hParChrg",           "", nChrgBin, rChrgBin[0], rChrgBin[1]);
  hParMass           = new TH1D("hParMass",           "", nMassBin, rMassBin[0], rMassBin[1]);
  hParEne            = new TH1D("hParEne",            "", nEneBin,  rEneBin[0],  rEneBin[1]);
  hParMom            = new TH1D("hParMom",            "", nEneBin,  rEneBin[0],  rEneBin[1]);
  hParMomX           = new TH1D("hParMomX",           "", nMomBin,  rMomBin[0],  rMomBin[1]);
  hParMomY           = new TH1D("hParMomY",           "", nMomBin,  rMomBin[0],  rMomBin[1]);
  hParMomZ           = new TH1D("hParMomZ",           "", nMomBin,  rMomBin[0],  rMomBin[1]);
  // reconstructed hit histograms
  hRecHitEne         = new TH1D("hRecHitEne",         "", nEneBin,  rEneBin[0],  rEneBin[1]);
  hRecHitPosZ        = new TH1D("hRecHitPosZ",        "", nPosBin,  rPosBin[0],  rPosBin[1]);
  hRecHitParDiff     = new TH1D("hRecHitParDiff",     "", nDiffBin, rDiffBin[0], rDiffBin[1]);
  hRecHitPosYvsX     = new TH2D("hRecHitPosYvsX",     "", nPosBin,  rPosBin[0],  rPosBin[1], nPosBin, rPosBin[0], rPosBin[1]);
  hRecHitVsParEne    = new TH2D("hRecHitVsParEne",    "", nEneBin,  rEneBin[0],  rEneBin[1], nEneBin, rEneBin[0], rEneBin[1]);
  // cluster histograms
  hClusterEne        = new TH1D("hClusterEne",        "", nEneBin,  rEneBin[0],  rEneBin[1]);
  hClusterPosZ       = new TH1D("hClusterPosZ",       "", nPosBin,  rPosBin[0],  rPosBin[1]);
  hClusterNumHit     = new TH1I("hClusterNumHit",     "", nNumBin,  rNumBin[0],  rNumBin[1]);
  hClusterParDiff    = new TH1D("hClusterParDiff",    "", nDiffBin, rDiffBin[0], rDiffBin[1]);
  hClusterPosYvsX    = new TH2D("hClusterPosYvsX",    "", nPosBin,  rPosBin[0],  rPosBin[1], nPosBin, rPosBin[0], rPosBin[1]);
  hClusterVsParEne   = new TH2D("hClusterVsParEne",   "", nEneBin,  rEneBin[0],  rEneBin[1], nEneBin, rEneBin[0], rEneBin[1]);
  // event-wise histograms
  hEvtNumPar         = new TH1I("hEvtNumPar",         "", nNumBin,  rNumBin[0],  rNumBin[1]);
  hEvtNumHit         = new TH1I("hEvtNumHit",         "", nNumBin,  rNumBin[0],  rNumBin[1]);
  hEvtNumClust       = new TH1I("hEvtNumClust",       "", nNumBin,  rNumBin[0],  rNumBin[1]);
  hEvtSumHitEne      = new TH1D("hEvtSumHitEne",      "", nEneBin,  rEneBin[0],  rEneBin[1]);
  hEvtSumClustEne    = new TH1D("hEvtSumClustEne",    "", nEneBin,  rEneBin[0],  rEneBin[1]);
  hEvtLeadClustEne   = new TH1D("hEvtLeadClustEne",   "", nEneBin,  rEneBin[0],  rEneBin[1]);
  hEvtSumHitDiff     = new TH1D("hEvtSumHitDiff",     "", nDiffBin, rDiffBin[0], rDiffBin[1]);
  hEvtSumClustDiff   = new TH1D("hEvtSumClustDiff",   "", nDiffBin, rDiffBin[0], rDiffBin[1]);
  hEvtLeadClustDiff  = new TH1D("hEvtLeadClustDiff",  "", nDiffBin, rDiffBin[0], rDiffBin[1]);
  hEvtSumHitVsPar    = new TH2D("hEvtSumHitVsPar",    "", nEneBin,  rEneBin[0],  rEneBin[1], nEneBin, rEneBin[0], rEneBin[1]);
  hEvtSumClustVsPar  = new TH2D("hEvtSumClustVsPar",  "", nEneBin,  rEneBin[0],  rEneBin[1], nEneBin, rEneBin[0], rEneBin[1]);
  hEvtLeadClustVsPar = new TH2D("hEvtLeadClustVsPar", "", nEneBin,  rEneBin[0],  rEneBin[1], nEneBin, rEneBin[0], rEneBin[1]);
  // errors
  hParChrg           -> Sumw2();
  hParMass           -> Sumw2();
  hParEne            -> Sumw2();
  hParMom            -> Sumw2();
  hParMomX           -> Sumw2();
  hParMomY           -> Sumw2();
  hParMomZ           -> Sumw2();
  hRecHitEne         -> Sumw2();
  hRecHitPosZ        -> Sumw2();
  hRecHitParDiff     -> Sumw2();
  hRecHitPosYvsX     -> Sumw2();
  hRecHitVsParEne    -> Sumw2();
  hClusterEne        -> Sumw2();
  hClusterPosZ       -> Sumw2();
  hClusterNumHit     -> Sumw2();
  hClusterParDiff    -> Sumw2();
  hClusterPosYvsX    -> Sumw2();
  hClusterVsParEne   -> Sumw2();
  hEvtNumPar         -> Sumw2();
  hEvtNumHit         -> Sumw2();
  hEvtNumClust       -> Sumw2();
  hEvtSumHitEne      -> Sumw2();
  hEvtSumClustEne    -> Sumw2();
  hEvtLeadClustEne   -> Sumw2();
  hEvtSumHitDiff     -> Sumw2();
  hEvtSumClustDiff   -> Sumw2();
  hEvtLeadClustDiff  -> Sumw2();
  hEvtSumHitVsPar    -> Sumw2();
  hEvtSumClustVsPar  -> Sumw2();
  hEvtLeadClustVsPar -> Sumw2();
  return;

}  // end 'InitWithGlobalRootLock()'




//-------------------------------------------
// ProcessSequential
//-------------------------------------------
void JCalibrateHCalProcessor::ProcessSequential(const std::shared_ptr<const JEvent>& event) {

  // MC particle properties
  float  cMcPar(0.);
  double mMcPar(0.);
  double eMcPar(0.);
  double pTotMcPar(0.);
  double pMcPar[NComp] = {0., 0., 0.};

  // particle loop
  UInt_t nPar(0);
  for (auto par : genParticles()) {

    // grab particle properties
    const auto cPar  = par -> getCharge();
    const auto mPar  = par -> getMass();
    const auto ePar  = par -> getEnergy();
    const auto pParX = par -> getMomentum().x;
    const auto pParY = par -> getMomentum().y;
    const auto pParZ = par -> getMomentum().z;
    const auto pPar  = std::sqrt((pParX * pParX) + (pParY * pParY) + (pParZ * pParZ));

    // select MC particle
    const bool isRightCharge   = (cPar == CPar);
    const bool isRightMass     = ((mPar >= MParMin) && (mPar <= MParMax));
    const bool isRightMomentum = ((pPar >= PParMin) && (pPar <= PParMax));
    const bool isMcParticle    = (isRightCharge && isRightMass && isRightMomentum);
    if (isMcParticle) {
      cMcPar    = cPar;
      mMcPar    = mPar;
      eMcPar    = ePar;
      pMcPar[0] = pParX;
      pMcPar[1] = pParY;
      pMcPar[2] = pParZ;
      pTotMcPar = pPar;
    }
    ++nPar;
  }  // end particle loop

  // fill particle histograms
  hParChrg -> Fill(cMcPar);
  hParMass -> Fill(mMcPar);
  hParEne  -> Fill(eMcPar);
  hParMom  -> Fill(pTotMcPar);
  hParMomX -> Fill(pMcPar[0]);
  hParMomY -> Fill(pMcPar[1]);
  hParMomZ -> Fill(pMcPar[2]);

  // hit and cluster sums
  double eHitSum(0.);
  double eClustSum(0.);

  // for highest energy cluster
  double eLeadClust(0.);
  double diffLeadClust(0.);

  // reconstructed hit loop
  UInt_t nHit(0);
  for (auto hit : bhcalRecHits()) {

    // grab hit properties
    const auto rHitX   = hit -> getPosition().x;
    const auto rHitY   = hit -> getPosition().y;
    const auto rHitZ   = hit -> getPosition().z;
    const auto eHit    = hit -> getEnergy();
    const auto diffHit = (eHit - eMcPar) / eHit;

    // fill hit histograms and increment sums/counters
    hRecHitEne      -> Fill(eHit);
    hRecHitPosZ     -> Fill(rHitZ);
    hRecHitParDiff  -> Fill(diffHit);
    hRecHitPosYvsX  -> Fill(rHitX, rHitY);
    hRecHitVsParEne -> Fill(eMcPar, eHit);
    eHitSum += eHit;
    ++nHit;
  }  // end hit loop

  // cluster loop
  UInt_t nClust(0);
  for (auto cluster : bhcalClusters()) {

    // grab hit properties
    const auto rClustX   = cluster -> getPosition().x;
    const auto rClustY   = cluster -> getPosition().y;
    const auto rClustZ   = cluster -> getPosition().z;
    const auto eClust    = cluster -> getEnergy();
    const auto nHitClust = cluster -> hits_size();
    const auto diffClust = (eClust - eMcPar) / eClust;

    // fill cluster histograms and increment sums/counters
    hClusterEne      -> Fill(eClust);
    hClusterPosZ     -> Fill(rClustZ);
    hClusterNumHit   -> Fill(nHitClust);
    hClusterParDiff  -> Fill(diffClust);
    hClusterPosYvsX  -> Fill(rClustX, rClustY);
    hClusterVsParEne -> Fill(eMcPar, eClust);
    eClustSum += eClust;
    ++nClust;

    // select leading cluster
    const bool isBiggerEne = (eClust > eLeadClust);
    if (isBiggerEne) {
      eLeadClust    = eClust;
      diffLeadClust = diffClust;
    }
  }  // end cluster loop

  // do event-wise calculations
  const auto diffHitSum   = (eHitSum - eMcPar) / eHitSum;
  const auto diffClustSum = (eClustSum - eMcPar) / eClustSum;

  // fill event-wise histograms
  hEvtNumPar         -> Fill(nPar);
  hEvtNumHit         -> Fill(nHit);
  hEvtNumClust       -> Fill(nClust);
  hEvtSumHitEne      -> Fill(eHitSum);
  hEvtSumClustEne    -> Fill(eClustSum);
  hEvtLeadClustEne   -> Fill(eLeadClust);
  hEvtSumHitDiff     -> Fill(diffHitSum);
  hEvtSumClustDiff   -> Fill(diffClustSum);
  hEvtLeadClustDiff  -> Fill(diffLeadClust);
  hEvtSumHitVsPar    -> Fill(eMcPar, eHitSum);
  hEvtSumClustVsPar  -> Fill(eMcPar, eClustSum);
  hEvtLeadClustVsPar -> Fill(eMcPar, eLeadClust);
  return;

}  // end 'ProcessSequential(std::shared_ptr<JEvent>&)'



//-------------------------------------------
// FinishWithGlobalRootLock
//-------------------------------------------
void JCalibrateHCalProcessor::FinishWithGlobalRootLock() {

  // axis titles
  const TString sCount("counts");
  const TString sCharge("charge");
  const TString sMass("m_{par} [GeV/c^{2}]");
  const TString sEnePar("E_{par} [GeV]");
  const TString sEneHit("e_{hit} [GeV]");
  const TString sEneClust("e_{clust} [GeV]");
  const TString sEneHitDiff("#Deltae_{hit} / e_{hit} = (e_{hit} - E_{par}) / e_{hit} [GeV]");
  const TString sEneClustDiff("#Deltae_{clust} / e_{clust} = (e_{clust} - E_{par}) / e_{clust} [GeV]");
  const TString sEneHitSum("E_{sum, hit} = #Sigmae_{hit} [GeV]");
  const TString sEneClustSum("E_{sum, clust} = #Sigmae_{clust} [GeV]");
  const TString sEneClustLead("E_{lead, clust} [GeV]");
  const TString sEneHitSumDiff("#DeltaE_{sum, hit} / E_{sum, hit} = (E_{sum, hit} - E_{par}) / E_{sum, hit} [GeV]");
  const TString sEneClustSumDiff("#DeltaE_{sum, clust} / E_{sum, clust} = (E_{sum, clust} - E_{par}) / E_{sum, clust} [GeV]");
  const TString sEneClustLeadDiff("#DeltaE_{lead, clust} / E_{lead, clust} = (E_{lead, clust} - E_{par}) / E_{lead, clust} [GeV]");
  const TString sMomPar("p_{par} [GeV/c]");
  const TString sMomParX("p_{x, par} [GeV/c]");
  const TString sMomParY("p_{y, par} [GeV/c]");
  const TString sMomParZ("p_{z, par} [GeV/c]");
  const TString sPosHitX("x_{hit} [mm]");
  const TString sPosHitY("y_{hit} [mm]");
  const TString sPosHitZ("z_{hit} [mm]");
  const TString sPosClustX("x_{clust} [mm]");
  const TString sPosClustY("y_{clust} [mm]");
  const TString sPosClustZ("z_{clust} [mm]");
  const TString sNumHitClust("N_{hit} per cluster");
  const TString sNumParEvt("N_{par} per event");
  const TString sNumHitEvt("N_{hit} per event");
  const TString sNumClustEvt("N_{clust} per event");

  // set particle axis titles
  hParChrg           -> GetXaxis() -> SetTitle(sCharge.Data());
  hParChrg           -> GetYaxis() -> SetTitle(sCount.Data());
  hParMass           -> GetXaxis() -> SetTitle(sMass.Data());
  hParMass           -> GetYaxis() -> SetTitle(sCount.Data());
  hParEne            -> GetXaxis() -> SetTitle(sEnePar.Data());
  hParEne            -> GetYaxis() -> SetTitle(sCount.Data());
  hParMom            -> GetXaxis() -> SetTitle(sMomPar.Data());
  hParMom            -> GetYaxis() -> SetTitle(sCount.Data());
  hParMomX           -> GetXaxis() -> SetTitle(sMomParX.Data());
  hParMomX           -> GetYaxis() -> SetTitle(sCount.Data());
  hParMomY           -> GetXaxis() -> SetTitle(sMomParY.Data());
  hParMomY           -> GetYaxis() -> SetTitle(sCount.Data());
  hParMomZ           -> GetXaxis() -> SetTitle(sMomParZ.Data());
  hParMomZ           -> GetYaxis() -> SetTitle(sCount.Data());
  // set reconstructed hit axis titles
  hRecHitEne         -> GetXaxis() -> SetTitle(sEneHit.Data());
  hRecHitEne         -> GetYaxis() -> SetTitle(sCount.Data());
  hRecHitPosZ        -> GetXaxis() -> SetTitle(sPosHitZ.Data());
  hRecHitPosZ        -> GetYaxis() -> SetTitle(sCount.Data());
  hRecHitParDiff     -> GetXaxis() -> SetTitle(sEneHitDiff.Data());
  hRecHitParDiff     -> GetYaxis() -> SetTitle(sCount.Data());
  hRecHitPosYvsX     -> GetXaxis() -> SetTitle(sPosHitX.Data());
  hRecHitPosYvsX     -> GetYaxis() -> SetTitle(sPosHitY.Data());
  hRecHitPosYvsX     -> GetZaxis() -> SetTitle(sCount.Data());
  hRecHitVsParEne    -> GetXaxis() -> SetTitle(sEnePar.Data());
  hRecHitVsParEne    -> GetYaxis() -> SetTitle(sEneHit.Data());
  hRecHitVsParEne    -> GetZaxis() -> SetTitle(sCount.Data());
  // set cluster axis titles
  hClusterEne        -> GetXaxis() -> SetTitle(sEneClust.Data());
  hClusterEne        -> GetYaxis() -> SetTitle(sCount.Data());
  hClusterPosZ       -> GetXaxis() -> SetTitle(sPosClustZ.Data());
  hClusterPosZ       -> GetYaxis() -> SetTitle(sCount.Data());
  hClusterNumHit     -> GetXaxis() -> SetTitle(sNumHitClust.Data());
  hClusterNumHit     -> GetYaxis() -> SetTitle(sCount.Data());
  hClusterParDiff    -> GetXaxis() -> SetTitle(sEneClustDiff.Data());
  hClusterParDiff    -> GetYaxis() -> SetTitle(sCount.Data());
  hClusterPosYvsX    -> GetXaxis() -> SetTitle(sPosClustX.Data());
  hClusterPosYvsX    -> GetYaxis() -> SetTitle(sPosClustY.Data());
  hClusterPosYvsX    -> GetZaxis() -> SetTitle(sCount.Data());
  hClusterVsParEne   -> GetXaxis() -> SetTitle(sEnePar.Data());
  hClusterVsParEne   -> GetYaxis() -> SetTitle(sEneClust.Data());
  hClusterVsParEne   -> GetZaxis() -> SetTitle(sCount.Data());
  // set event-wise axis titles
  hEvtNumPar         -> GetXaxis() -> SetTitle(sNumParEvt.Data());
  hEvtNumPar         -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtNumHit         -> GetXaxis() -> SetTitle(sNumHitEvt.Data());
  hEvtNumHit         -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtNumClust       -> GetXaxis() -> SetTitle(sNumClustEvt.Data());
  hEvtNumClust       -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtSumHitEne      -> GetXaxis() -> SetTitle(sEneHitSum.Data());
  hEvtSumHitEne      -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtSumClustEne    -> GetXaxis() -> SetTitle(sEneClustSum.Data());
  hEvtSumClustEne    -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtLeadClustEne   -> GetXaxis() -> SetTitle(sEneClustLead.Data());
  hEvtLeadClustEne   -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtSumHitDiff     -> GetXaxis() -> SetTitle(sEneHitSumDiff.Data());
  hEvtSumHitDiff     -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtSumClustDiff   -> GetXaxis() -> SetTitle(sEneClustSumDiff.Data());
  hEvtSumClustDiff   -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtLeadClustDiff  -> GetXaxis() -> SetTitle(sEneClustLeadDiff.Data());
  hEvtLeadClustDiff  -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtSumHitVsPar    -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtSumHitVsPar    -> GetYaxis() -> SetTitle(sEneHitSum.Data());
  hEvtSumHitVsPar    -> GetZaxis() -> SetTitle(sCount.Data());
  hEvtSumClustVsPar  -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtSumClustVsPar  -> GetYaxis() -> SetTitle(sEneClustSum.Data());
  hEvtSumClustVsPar  -> GetZaxis() -> SetTitle(sCount.Data());
  hEvtLeadClustVsPar -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtLeadClustVsPar -> GetYaxis() -> SetTitle(sEneClustLead.Data());
  hEvtLeadClustVsPar -> GetZaxis() -> SetTitle(sCount.Data());
  return;

}  // end 'FinishWithGlobalRootLock()'

// end ------------------------------------------------------------------------