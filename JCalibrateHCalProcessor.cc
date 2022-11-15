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
  const unsigned long nEneBin(200);
  const unsigned long nMomBin(200);
  const unsigned long nPosTrBin(800);
  const unsigned long nPosLoBin(30);
  const unsigned long nDiffBin(200);
  const unsigned long rNumBin[NRange]   = {0,      200};
  const double        rChrgBin[NRange]  = {-3.,    3.};
  const double        rMassBin[NRange]  = {0.,     5.};
  const double        rEneBin[NRange]   = {0.,     100.};
  const double        rMomBin[NRange]   = {-50.,   50.};
  const double        rPosTrBin[NRange] = {-4000., 4000.};
  const double        rPosLoBin[NRange] = {-3000., 3000.};
  const double        rDiffBin[NRange]  = {-50.,   50.};
  // particle histograms
  hParChrg               = new TH1D("hParChrg",               "", nChrgBin,  rChrgBin[0],  rChrgBin[1]);
  hParMass               = new TH1D("hParMass",               "", nMassBin,  rMassBin[0],  rMassBin[1]);
  hParEne                = new TH1D("hParEne",                "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hParMom                = new TH1D("hParMom",                "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hParMomX               = new TH1D("hParMomX",               "", nMomBin,   rMomBin[0],   rMomBin[1]);
  hParMomY               = new TH1D("hParMomY",               "", nMomBin,   rMomBin[0],   rMomBin[1]);
  hParMomZ               = new TH1D("hParMomZ",               "", nMomBin,   rMomBin[0],   rMomBin[1]);
  // reconstructed hit histograms
  hRecHitEne             = new TH1D("hRecHitEne",             "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hRecHitPosZ            = new TH1D("hRecHitPosZ",            "", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hRecHitParDiff         = new TH1D("hRecHitParDiff",         "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hRecHitPosYvsX         = new TH2D("hRecHitPosYvsX",         "", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hRecHitVsParEne        = new TH2D("hRecHitVsParEne",        "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // reconstructed cluster histograms
  hClusterEne            = new TH1D("hClusterEne",            "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hClusterPosZ           = new TH1D("hClusterPosZ",           "", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hClusterNumHit         = new TH1I("hClusterNumHit",         "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hClusterParDiff        = new TH1D("hClusterParDiff",        "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hClusterPosYvsX        = new TH2D("hClusterPosYvsX",        "", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hClusterVsParEne       = new TH2D("hClusterVsParEne",       "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // reconstructed cluster debug histograms
  hDebugClusterSum5      = new TH1D("hDebugClusterSum5",      "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hDebugClusterSum10     = new TH1D("hDebugClusterSum10",     "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hDebugClusterSum100    = new TH1D("hDebugClusterSum100",    "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hDebugClusterSum1000   = new TH1D("hDebugClusterSum1000",   "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hDebugClusterDiff5     = new TH1D("hDebugClusterDiff5",     "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hDebugClusterDiff10    = new TH1D("hDebugClusterDiff10",    "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hDebugClusterDiff100   = new TH1D("hDebugClusterDiff100",   "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hDebugClusterDiff1000  = new TH1D("hDebugClusterDiff1000",  "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  // truth cluster histograms
  hTruClustEne           = new TH1D("hTruClustEne",           "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hTruClustPosZ          = new TH1D("hTruClustPosZ",          "", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hTruClustNumHit        = new TH1I("hTruClustNumHit",        "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hTruClustParDiff       = new TH1D("hTruClustParDiff",       "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hTruClustPosYvsX       = new TH2D("hTruClustPosYvsX",       "", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hTruClustVsParEne      = new TH2D("hTruClustVsParEne",      "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // truth cluster debug histograms
  hDebugTruClustSum5     = new TH1D("hDebugTruClustSum5",     "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hDebugTruClustSum10    = new TH1D("hDebugTruClustSum10",    "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hDebugTruClustSum100   = new TH1D("hDebugTruClustSum100",   "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hDebugTruClustSum1000  = new TH1D("hDebugTruClustSum1000",  "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hDebugTruClustDiff5    = new TH1D("hDebugTruClustDiff5",    "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hDebugTruClustDiff10   = new TH1D("hDebugTruClustDiff10",   "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hDebugTruClustDiff100  = new TH1D("hDebugTruClustDiff100",  "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hDebugTruClustDiff1000 = new TH1D("hDebugTruClustDiff1000", "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  // event-wise histograms
  hEvtNumPar             = new TH1I("hEvtNumPar",             "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtNumHit             = new TH1I("hEvtNumHit",             "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtNumClust           = new TH1I("hEvtNumClust",           "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtNumTruClust        = new TH1I("hEvtNumTruClust",        "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtSumHitEne          = new TH1D("hEvtSumHitEne",          "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtSumClustEne        = new TH1D("hEvtSumClustEne",        "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtSumTruClustEne     = new TH1D("hEvtSumTruClustEne",     "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtLeadClustEne       = new TH1D("hEvtLeadClustEne",       "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtLeadTruClustEne    = new TH1D("hEvtLeadTruClustEne",    "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtSumHitDiff         = new TH1D("hEvtSumHitDiff",         "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtSumClustDiff       = new TH1D("hEvtSumClustDiff",       "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtSumTruClustDiff    = new TH1D("hEvtSumTruClustDiff",    "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtLeadClustDiff      = new TH1D("hEvtLeadClustDiff",      "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtLeadTruClustDiff   = new TH1D("hEvtLeadTruClustDiff",   "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtNumClustVsHit      = new TH2I("hEvtNumClustVsHit",      "", nNumBin,   rNumBin[0],   rNumBin[1],   nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtNumTruClustVsClust = new TH2I("hEvtNumTruClustVsClust", "", nNumBin,   rNumBin[0],   rNumBin[1],   nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtSumHitVsPar        = new TH2D("hEvtSumHitVsPar",        "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtSumClustVsPar      = new TH2D("hEvtSumClustVsPar",      "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtSumTruClustVsPar   = new TH2D("hEvtSumTruClustVsPar",   "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtLeadClustVsPar     = new TH2D("hEvtLeadClustVsPar",     "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtLeadTruClustVsPar  = new TH2D("hEvtLeadTruClustVsPar",  "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // errors
  hParChrg               -> Sumw2();
  hParMass               -> Sumw2();
  hParEne                -> Sumw2();
  hParMom                -> Sumw2();
  hParMomX               -> Sumw2();
  hParMomY               -> Sumw2();
  hParMomZ               -> Sumw2();
  hRecHitEne             -> Sumw2();
  hRecHitPosZ            -> Sumw2();
  hRecHitParDiff         -> Sumw2();
  hRecHitPosYvsX         -> Sumw2();
  hRecHitVsParEne        -> Sumw2();
  hClusterEne            -> Sumw2();
  hClusterPosZ           -> Sumw2();
  hClusterNumHit         -> Sumw2();
  hClusterParDiff        -> Sumw2();
  hClusterPosYvsX        -> Sumw2();
  hClusterVsParEne       -> Sumw2();
  hDebugClusterSum5      -> Sumw2();
  hDebugClusterSum10     -> Sumw2();
  hDebugClusterSum100    -> Sumw2();
  hDebugClusterSum1000   -> Sumw2();
  hDebugClusterDiff5     -> Sumw2();
  hDebugClusterDiff10    -> Sumw2();
  hDebugClusterDiff100   -> Sumw2();
  hDebugClusterDiff1000  -> Sumw2();
  hTruClustEne           -> Sumw2();
  hTruClustPosZ          -> Sumw2();
  hTruClustNumHit        -> Sumw2();
  hTruClustParDiff       -> Sumw2();
  hTruClustPosYvsX       -> Sumw2();
  hTruClustVsParEne      -> Sumw2();
  hDebugTruClustSum5     -> Sumw2();
  hDebugTruClustSum10    -> Sumw2();
  hDebugTruClustSum100   -> Sumw2();
  hDebugTruClustSum1000  -> Sumw2();
  hDebugTruClustDiff5    -> Sumw2();
  hDebugTruClustDiff10   -> Sumw2();
  hDebugTruClustDiff100  -> Sumw2();
  hDebugTruClustDiff1000 -> Sumw2();
  hEvtNumPar             -> Sumw2();
  hEvtNumHit             -> Sumw2();
  hEvtNumClust           -> Sumw2();
  hEvtNumTruClust        -> Sumw2();
  hEvtSumHitEne          -> Sumw2();
  hEvtSumClustEne        -> Sumw2();
  hEvtSumTruClustEne     -> Sumw2();
  hEvtLeadClustEne       -> Sumw2();
  hEvtLeadTruClustEne    -> Sumw2();
  hEvtSumHitDiff         -> Sumw2();
  hEvtSumClustDiff       -> Sumw2();
  hEvtSumTruClustDiff    -> Sumw2();
  hEvtLeadClustDiff      -> Sumw2();
  hEvtLeadTruClustDiff   -> Sumw2();
  hEvtNumClustVsHit      -> Sumw2();
  hEvtNumTruClustVsClust -> Sumw2();
  hEvtSumHitVsPar        -> Sumw2();
  hEvtSumClustVsPar      -> Sumw2();
  hEvtLeadClustVsPar     -> Sumw2();
  hEvtLeadTruClustVsPar  -> Sumw2();
  return;

}  // end 'InitWithGlobalRootLock()'




//-------------------------------------------
// ProcessSequential
//-------------------------------------------
void JCalibrateHCalProcessor::ProcessSequential(const std::shared_ptr<const JEvent>& event) {

  // hit and cluster sums
  double eHitSum(0.);
  double eClustSum(0.);
  double eTruClustSum(0.);

  // sum hit energy
  for (auto hit : bhcalRecHits()) {
    eHitSum += hit -> getEnergy();
  }  // end 1st hit loop

  // if hit sum is 0, skip event
  const bool isHitSumNonzero = (eHitSum > 0.);
  if (!isHitSumNonzero) {
    return;
  }

  // MC particle properties
  float  cMcPar(0.);
  double mMcPar(0.);
  double eMcPar(0.);
  double pTotMcPar(0.);
  double pMcPar[NComp] = {0., 0., 0.};

  // particle loop
  unsigned long nPar(0);
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

  // reconstructed hit loop
  unsigned long nHit(0);
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
    ++nHit;
  }  // end 2nd hit loop

  // for highest energy clusters
  int    iLeadClust(-1);
  int    iLeadTruClust(-1);
  double eLeadClust(0.);
  double eLeadTruClust(0.);
  double diffLeadClust(0.);
  double diffLeadTruClust(0.);

  // reconstructed cluster loop
  unsigned long iClust(0);
  unsigned long nClust(0);
  for (auto cluster : bhcalClusters()) {

    // grab cluster properties
    const auto rClustX   = cluster -> getPosition().x;
    const auto rClustY   = cluster -> getPosition().y;
    const auto rClustZ   = cluster -> getPosition().z;
    const auto eClust    = cluster -> getEnergy();
    const auto nHitClust = cluster -> hits_size();
    const auto diffClust = (eClust - eMcPar) / eClust;

    // fill cluster histograms and increment counters
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
      iLeadClust    = iClust;
      eLeadClust    = eClust;
      diffLeadClust = diffClust;
    }
    ++iClust;
  }  // end reconstructed cluster loop

  // for debugging reco clusters
  double eDebugSumClust5(0.);
  double eDebugSumClust10(0.);
  double eDebugSumClust100(0.);
  double eDebugSumClust1000(0.);

  // debug reconstructed cluster loop
  unsigned long iDebugClust(0);
  for (auto debugCluster : bhcalClusters()) {

    // select leading cluster
    const bool isLeadCluster = (iDebugClust == iLeadClust);
    if (isLeadCluster) {

      // grab lead cluster properties
      const auto rLeadClustX = debugCluster -> getPosition().x;
      const auto rLeadClustY = debugCluster -> getPosition().y;
      const auto eLeadClust  = debugCluster -> getEnergy();

      unsigned long iOtherClust(0);
      double        toAddClust5(0.);
      double        toAddClust10(0.);
      double        toAddClust100(0.);
      double        toAddClust1000(0.);
      for (auto otherCluster : bhcalClusters()) {

        // ignore same cluster
        const bool isSameClust = (iOtherClust == iDebugClust);
        if (isSameClust) continue;

        // grab other cluster properties
        const auto rOtherClustX = otherCluster -> getPosition().x;
        const auto rOtherClustY = otherCluster -> getPosition().y;
        const auto eOtherClust  = otherCluster -> getEnergy();
        const auto drLeadOtherX = rOtherClustX - rLeadClustX;
        const auto drLeadOtherY = rOtherClustY - rLeadClustY;
        const auto drLeadOther  = std::sqrt((drLeadOtherX * drLeadOtherX) + (drLeadOtherY * drLeadOtherY));

        // increment relevant sums and counters
        const bool isIn5mm    = (drLeadOther < 5.);
        const bool isIn10mm   = (drLeadOther < 10.);
        const bool isIn100mm  = (drLeadOther < 100.);
        const bool isIn1000mm = (drLeadOther < 1000.);
        if (isIn5mm)    toAddClust5    += eOtherClust;
        if (isIn10mm)   toAddClust10   += eOtherClust;
        if (isIn100mm)  toAddClust100  += eOtherClust;
        if (isIn1000mm) toAddClust1000 += eOtherClust;
        ++iOtherClust;
      }  // end other reco. cluster loo

      // add sums to lead energy
      eDebugSumClust5     = eLeadClust;
      eDebugSumClust10    = eLeadClust;
      eDebugSumClust100   = eLeadClust;
      eDebugSumClust1000  = eLeadClust;
      eDebugSumClust5    += toAddClust5;
      eDebugSumClust10   += toAddClust10;
      eDebugSumClust100  += toAddClust100;
      eDebugSumClust1000 += toAddClust1000;

    }
    ++iDebugClust;
  }  // end debug reco. cluster loop

  // true cluster loop
  unsigned long iTruClust(0);
  unsigned long nTruClust(0);
  for (auto truthCluster : bhcalTruthClusters()) {

    // grab cluster properties
    const auto rTruClustX   = truthCluster -> getPosition().x;
    const auto rTruClustY   = truthCluster -> getPosition().y;
    const auto rTruClustZ   = truthCluster -> getPosition().z;
    const auto eTruClust    = truthCluster -> getEnergy();
    const auto nHitTruClust = truthCluster -> hits_size();
    const auto diffTruClust = (eTruClust - eMcPar) / eTruClust;

    // fill cluster histograms and increment counters
    hTruClustEne      -> Fill(eTruClust);
    hTruClustPosZ     -> Fill(rTruClustZ);
    hTruClustNumHit   -> Fill(nHitTruClust);
    hTruClustParDiff  -> Fill(diffTruClust);
    hTruClustPosYvsX  -> Fill(rTruClustX, rTruClustY);
    hTruClustVsParEne -> Fill(eMcPar, eTruClust);
    eTruClustSum += eTruClust;
    ++nTruClust;

    // select leading cluster
    const bool isBiggerEne = (eTruClust > eLeadTruClust);
    if (isBiggerEne) {
      iLeadTruClust    = iTruClust;
      eLeadTruClust    = eTruClust;
      diffLeadTruClust = diffTruClust;
    }
    ++iTruClust;
  }  // end true cluster loop

  // for debugging truth clusters
  double eDebugSumTruClust5(0.);
  double eDebugSumTruClust10(0.);
  double eDebugSumTruClust100(0.);
  double eDebugSumTruClust1000(0.);

  // debug reconstructed cluster loop
  unsigned long iDebugTruClust(0);
  for (auto debugTruthCluster : bhcalTruthClusters()) {

    // select leading cluster
    const bool isLeadTruCluster = (iDebugTruClust == iLeadTruClust);
    if (isLeadTruCluster) {

      // grab lead cluster properties
      const auto rLeadTruClustX = debugTruthCluster -> getPosition().x;
      const auto rLeadTruClustY = debugTruthCluster -> getPosition().y;
      const auto eLeadTruClust  = debugTruthCluster -> getEnergy();

      unsigned long iOtherTruClust(0);
      double        toAddTruClust5(0.);
      double        toAddTruClust10(0.);
      double        toAddTruClust100(0.);
      double        toAddTruClust1000(0.);
      for (auto otherTruthCluster : bhcalTruthClusters()) {

        // ignore same cluster
        const bool isSameTruClust = (iOtherTruClust == iDebugTruClust);
        if (isSameTruClust) continue;

        // grab other cluster properties
        const auto rOtherTruClustX = otherTruthCluster -> getPosition().x;
        const auto rOtherTruClustY = otherTruthCluster -> getPosition().y;
        const auto eOtherTruClust  = otherTruthCluster -> getEnergy();
        const auto drLeadOtherX    = rOtherTruClustX - rLeadTruClustX;
        const auto drLeadOtherY    = rOtherTruClustY - rLeadTruClustY;
        const auto drLeadOther     = std::sqrt((drLeadOtherX * drLeadOtherX) + (drLeadOtherY * drLeadOtherY));

        // increment relevant sums and counters
        const bool isIn5mm    = (drLeadOther < 5);
        const bool isIn10mm   = (drLeadOther < 10);
        const bool isIn100mm  = (drLeadOther < 100);
        const bool isIn1000mm = (drLeadOther < 1000);
        if (isIn5mm)    toAddTruClust5    += eOtherTruClust;
        if (isIn10mm)   toAddTruClust10   += eOtherTruClust;
        if (isIn100mm)  toAddTruClust100  += eOtherTruClust;
        if (isIn1000mm) toAddTruClust1000 += eOtherTruClust;
        ++iOtherTruClust;
      }  // end other true cluster loo

      // add sums to lead energy
      eDebugSumTruClust5     = eLeadTruClust;
      eDebugSumTruClust10    = eLeadTruClust;
      eDebugSumTruClust100   = eLeadTruClust;
      eDebugSumTruClust1000  = eLeadTruClust;
      eDebugSumTruClust5    += toAddTruClust5;
      eDebugSumTruClust10   += toAddTruClust10;
      eDebugSumTruClust100  += toAddTruClust100;
      eDebugSumTruClust1000 += toAddTruClust1000;

    }
    ++iDebugTruClust;
  }  // end debug true cluster loop

  // do event-wise calculations
  const auto diffHitSum      = (eHitSum - eMcPar) / eHitSum;
  const auto diffClustSum    = (eClustSum - eMcPar) / eClustSum;
  const auto diffTruClustSum = (eTruClustSum - eMcPar) / eTruClustSum;

  // fill event-wise histograms
  hEvtNumPar             -> Fill(nPar);
  hEvtNumHit             -> Fill(nHit);
  hEvtNumClust           -> Fill(nClust);
  hEvtNumTruClust        -> Fill(nTruClust);
  hEvtSumHitEne          -> Fill(eHitSum);
  hEvtSumClustEne        -> Fill(eClustSum);
  hEvtSumTruClustEne     -> Fill(eTruClustSum);
  hEvtLeadClustEne       -> Fill(eLeadClust);
  hEvtLeadTruClustEne    -> Fill(eLeadTruClust);
  hEvtSumHitDiff         -> Fill(diffHitSum);
  hEvtSumClustDiff       -> Fill(diffClustSum);
  hEvtSumTruClustDiff    -> Fill(diffTruClustSum);
  hEvtLeadClustDiff      -> Fill(diffLeadClust);
  hEvtLeadTruClustDiff   -> Fill(diffLeadTruClust);
  hEvtNumClustVsHit      -> Fill(nHit, nClust);
  hEvtNumTruClustVsClust -> Fill(nClust, nTruClust);
  hEvtSumHitVsPar        -> Fill(eMcPar, eHitSum);
  hEvtSumClustVsPar      -> Fill(eMcPar, eClustSum);
  hEvtSumTruClustVsPar   -> Fill(eMcPar, eTruClustSum);
  hEvtLeadClustVsPar     -> Fill(eMcPar, eLeadClust);
  hEvtLeadTruClustVsPar  -> Fill(eMcPar, eLeadTruClust);

  // do debugging calculations
  const auto diffDebugClustSum5       = (eDebugSumClust5 - eMcPar) / eDebugSumClust5;
  const auto diffDebugClustSum10      = (eDebugSumClust10 - eMcPar) / eDebugSumClust5;
  const auto diffDebugClustSum100     = (eDebugSumClust100 - eMcPar) / eDebugSumClust5;
  const auto diffDebugClustSum1000    = (eDebugSumClust1000 - eMcPar) / eDebugSumClust5;
  const auto diffDebugTruClustSum5    = (eDebugSumTruClust5 - eMcPar) / eDebugSumTruClust5;
  const auto diffDebugTruClustSum10   = (eDebugSumTruClust10 - eMcPar) / eDebugSumTruClust5;
  const auto diffDebugTruClustSum100  = (eDebugSumTruClust100 - eMcPar) / eDebugSumTruClust5;
  const auto diffDebugTruClustSum1000 = (eDebugSumTruClust1000 - eMcPar) / eDebugSumTruClust5;

  // fill reco. cluster debug histograms
  hDebugClusterSum5       -> Fill(eDebugSumClust5);
  hDebugClusterSum10      -> Fill(eDebugSumClust10);
  hDebugClusterSum100     -> Fill(eDebugSumClust100);
  hDebugClusterSum1000    -> Fill(eDebugSumClust1000);
  hDebugTruClustSum5      -> Fill(eDebugSumClust5);
  hDebugTruClustSum10     -> Fill(eDebugSumClust10);
  hDebugTruClustSum100    -> Fill(eDebugSumClust100);
  hDebugTruClustSum1000   -> Fill(eDebugSumClust1000);
  hDebugClusterDiff5      -> Fill(diffDebugClustSum5);
  hDebugClusterDiff10     -> Fill(diffDebugClustSum10);
  hDebugClusterDiff100    -> Fill(diffDebugClustSum100);
  hDebugClusterDiff1000   -> Fill(diffDebugClustSum1000);
  hDebugTruClustDiff5     -> Fill(diffDebugClustSum5);
  hDebugTruClustDiff10    -> Fill(diffDebugClustSum10);
  hDebugTruClustDiff100   -> Fill(diffDebugClustSum100);
  hDebugTruClustDiff1000  -> Fill(diffDebugClustSum1000);
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
  const TString sEneTruClust("e^{truth}_{clust} [GeV]");
  const TString sEneHitDiff("#Deltae_{hit} / e_{hit} = (e_{hit} - E_{par}) / e_{hit} [GeV]");
  const TString sEneClustDiff("#Deltae_{clust} / e_{clust} = (e_{clust} - E_{par}) / e_{clust} [GeV]");
  const TString sEneTruClustDiff("#Deltae^{truth}_{clust} / e^{truth}_{clust} / (e^{truth}_{clust} - E_{par}) / e^{truth}_{clust} [GeV]");
  const TString sEneHitSum("E^{sum}_{hit} = #Sigmae_{hit} [GeV]");
  const TString sEneClustSum("E^{sum}_{clust} = #Sigmae_{clust} [GeV]");
  const TString sEneTruClustSum("E^{sum/truth}_{clust} = #Sigmae^{truth}_{clust} [GeV]");
  const TString sEneClustLead("E^{lead}_{clust} [GeV]");
  const TString sEneTruClustLead("E^{lead/truth}_{clust} [GeV]");
  const TString sEneHitSumDiff("#DeltaE^{sum}_{hit} / E^{sum}_{hit} = (E^{sum}_{hit} - E_{par}) / E^{sum}_{hit} [GeV]");
  const TString sEneClustSumDiff("#DeltaE^{sum}_{clust} / E^{sum}_{clust} = (E^{sum}_{clust} - E_{par}) / E^{sum}_{clust} [GeV]");
  const TString sEneTruClustSumDiff("#DeltaE^{sum/truth}_{clust} / E^{sum/truth}_{clust} = (E^{sum/truth}_{clust} - E_{par}) / E^{sum/truth}_{clust} [GeV]");
  const TString sEneClustLeadDiff("#DeltaE^{lead}_{clust} / E^{lead}_{clust} = (E^{lead}_{clust} - E_{par}) / E^{lead}_{clust} [GeV]");
  const TString sEneTruClustLeadDiff("#DeltaE^{lead/truth}_{clust} / E^{lead/truth}_{clust} = (E^{lead/truth} _{clust} - E_{par}) / E^{lead/truth}_{clust} [GeV]");
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
  const TString sPosTruClustX("x_{truth clust} [mm]");
  const TString sPosTruClustY("y_{truth clust} [mm]");
  const TString sPosTruClustZ("z_{truth clust} [mm]");
  const TString sNumHitClust("N_{hit} per cluster");
  const TString sNumHitTruClust("N_{hit} per truth cluster");
  const TString sNumParEvt("N_{par} per event");
  const TString sNumHitEvt("N_{hit} per event");
  const TString sNumClustEvt("N_{clust} per event");
  const TString sNumTruClustEvt("N_{truth clust} per event");

  // set particle axis titles
  hParChrg              -> GetXaxis() -> SetTitle(sCharge.Data());
  hParChrg              -> GetYaxis() -> SetTitle(sCount.Data());
  hParMass              -> GetXaxis() -> SetTitle(sMass.Data());
  hParMass              -> GetYaxis() -> SetTitle(sCount.Data());
  hParEne               -> GetXaxis() -> SetTitle(sEnePar.Data());
  hParEne               -> GetYaxis() -> SetTitle(sCount.Data());
  hParMom               -> GetXaxis() -> SetTitle(sMomPar.Data());
  hParMom               -> GetYaxis() -> SetTitle(sCount.Data());
  hParMomX              -> GetXaxis() -> SetTitle(sMomParX.Data());
  hParMomX              -> GetYaxis() -> SetTitle(sCount.Data());
  hParMomY              -> GetXaxis() -> SetTitle(sMomParY.Data());
  hParMomY              -> GetYaxis() -> SetTitle(sCount.Data());
  hParMomZ              -> GetXaxis() -> SetTitle(sMomParZ.Data());
  hParMomZ              -> GetYaxis() -> SetTitle(sCount.Data());
  // set reconstructed hit axis titles
  hRecHitEne            -> GetXaxis() -> SetTitle(sEneHit.Data());
  hRecHitEne            -> GetYaxis() -> SetTitle(sCount.Data());
  hRecHitPosZ           -> GetXaxis() -> SetTitle(sPosHitZ.Data());
  hRecHitPosZ           -> GetYaxis() -> SetTitle(sCount.Data());
  hRecHitParDiff        -> GetXaxis() -> SetTitle(sEneHitDiff.Data());
  hRecHitParDiff        -> GetYaxis() -> SetTitle(sCount.Data());
  hRecHitPosYvsX        -> GetXaxis() -> SetTitle(sPosHitX.Data());
  hRecHitPosYvsX        -> GetYaxis() -> SetTitle(sPosHitY.Data());
  hRecHitPosYvsX        -> GetZaxis() -> SetTitle(sCount.Data());
  hRecHitVsParEne       -> GetXaxis() -> SetTitle(sEnePar.Data());
  hRecHitVsParEne       -> GetYaxis() -> SetTitle(sEneHit.Data());
  hRecHitVsParEne       -> GetZaxis() -> SetTitle(sCount.Data());
  // set reconstructed cluster axis titles
  hClusterEne           -> GetXaxis() -> SetTitle(sEneClust.Data());
  hClusterEne           -> GetYaxis() -> SetTitle(sCount.Data());
  hClusterPosZ          -> GetXaxis() -> SetTitle(sPosClustZ.Data());
  hClusterPosZ          -> GetYaxis() -> SetTitle(sCount.Data());
  hClusterNumHit        -> GetXaxis() -> SetTitle(sNumHitClust.Data());
  hClusterNumHit        -> GetYaxis() -> SetTitle(sCount.Data());
  hClusterParDiff       -> GetXaxis() -> SetTitle(sEneClustDiff.Data());
  hClusterParDiff       -> GetYaxis() -> SetTitle(sCount.Data());
  hClusterPosYvsX       -> GetXaxis() -> SetTitle(sPosClustX.Data());
  hClusterPosYvsX       -> GetYaxis() -> SetTitle(sPosClustY.Data());
  hClusterPosYvsX       -> GetZaxis() -> SetTitle(sCount.Data());
  hClusterVsParEne      -> GetXaxis() -> SetTitle(sEnePar.Data());
  hClusterVsParEne      -> GetYaxis() -> SetTitle(sEneClust.Data());
  hClusterVsParEne      -> GetZaxis() -> SetTitle(sCount.Data());
  // set truth cluster axis titles
  hTruClustEne          -> GetXaxis() -> SetTitle(sEneTruClust.Data());
  hTruClustEne          -> GetYaxis() -> SetTitle(sCount.Data());
  hTruClustPosZ         -> GetXaxis() -> SetTitle(sPosTruClustZ.Data());
  hTruClustPosZ         -> GetYaxis() -> SetTitle(sCount.Data());
  hTruClustNumHit       -> GetXaxis() -> SetTitle(sNumHitTruClust.Data());
  hTruClustNumHit       -> GetYaxis() -> SetTitle(sCount.Data());
  hTruClustParDiff      -> GetXaxis() -> SetTitle(sEneTruClustDiff.Data());
  hTruClustParDiff      -> GetYaxis() -> SetTitle(sCount.Data());
  hTruClustPosYvsX      -> GetXaxis() -> SetTitle(sPosTruClustX.Data());
  hTruClustPosYvsX      -> GetYaxis() -> SetTitle(sPosTruClustY.Data());
  hTruClustPosYvsX      -> GetZaxis() -> SetTitle(sCount.Data());
  hTruClustVsParEne     -> GetXaxis() -> SetTitle(sEnePar.Data());
  hTruClustVsParEne     -> GetYaxis() -> SetTitle(sEneTruClust.Data());
  hTruClustVsParEne     -> GetZaxis() -> SetTitle(sCount.Data());
  // set event-wise axis titles
  hEvtNumPar            -> GetXaxis() -> SetTitle(sNumParEvt.Data());
  hEvtNumPar            -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtNumHit            -> GetXaxis() -> SetTitle(sNumHitEvt.Data());
  hEvtNumHit            -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtNumClust          -> GetXaxis() -> SetTitle(sNumClustEvt.Data());
  hEvtNumClust          -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtNumTruClust       -> GetXaxis() -> SetTitle(sNumTruClustEvt.Data());
  hEvtNumTruClust       -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtSumHitEne         -> GetXaxis() -> SetTitle(sEneHitSum.Data());
  hEvtSumHitEne         -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtSumClustEne       -> GetXaxis() -> SetTitle(sEneClustSum.Data());
  hEvtSumClustEne       -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtSumTruClustEne    -> GetXaxis() -> SetTitle(sEneTruClustSum.Data());
  hEvtSumTruClustEne    -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtLeadClustEne      -> GetXaxis() -> SetTitle(sEneClustLead.Data());
  hEvtLeadClustEne      -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtLeadTruClustEne   -> GetXaxis() -> SetTitle(sEneTruClustLead.Data());
  hEvtLeadTruClustEne   -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtSumHitDiff        -> GetXaxis() -> SetTitle(sEneHitSumDiff.Data());
  hEvtSumHitDiff        -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtSumClustDiff      -> GetXaxis() -> SetTitle(sEneClustSumDiff.Data());
  hEvtSumClustDiff      -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtSumTruClustDiff   -> GetXaxis() -> SetTitle(sEneTruClustSumDiff.Data());
  hEvtSumTruClustDiff   -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtLeadClustDiff     -> GetXaxis() -> SetTitle(sEneClustLeadDiff.Data());
  hEvtLeadClustDiff     -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtLeadTruClustDiff  -> GetXaxis() -> SetTitle(sEneTruClustLeadDiff.Data());
  hEvtLeadTruClustDiff  -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtSumHitVsPar       -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtSumHitVsPar       -> GetYaxis() -> SetTitle(sEneHitSum.Data());
  hEvtSumHitVsPar       -> GetZaxis() -> SetTitle(sCount.Data());
  hEvtSumClustVsPar     -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtSumClustVsPar     -> GetYaxis() -> SetTitle(sEneClustSum.Data());
  hEvtSumClustVsPar     -> GetZaxis() -> SetTitle(sCount.Data());
  hEvtSumTruClustVsPar  -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtSumTruClustVsPar  -> GetYaxis() -> SetTitle(sEneTruClustSum.Data());
  hEvtSumTruClustVsPar  -> GetZaxis() -> SetTitle(sCount.Data());
  hEvtLeadClustVsPar    -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtLeadClustVsPar    -> GetYaxis() -> SetTitle(sEneClustLead.Data());
  hEvtLeadClustVsPar    -> GetZaxis() -> SetTitle(sCount.Data());
  hEvtLeadTruClustVsPar -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtLeadTruClustVsPar -> GetYaxis() -> SetTitle(sEneTruClustLead.Data());
  hEvtLeadTruClustVsPar -> GetZaxis() -> SetTitle(sCount.Data());
  return;

}  // end 'FinishWithGlobalRootLock()'

// end ------------------------------------------------------------------------
