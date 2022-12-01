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
  hHCalRecHitEne             = new TH1D("hHCalRecHitEne",             "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalRecHitPosZ            = new TH1D("hHCalRecHitPosZ",            "", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hHCalRecHitParDiff         = new TH1D("hHCalRecHitParDiff",         "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalRecHitPosYvsX         = new TH2D("hHCalRecHitPosYvsX",         "", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hHCalRecHitVsParEne        = new TH2D("hHCalRecHitVsParEne",        "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // reconstructed cluster histograms
  hHCalClustEne            = new TH1D("hHCalClustEne",            "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalClustPosZ           = new TH1D("hHCalClustPosZ",           "", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hHCalClustNumHit         = new TH1I("hHCalClustNumHit",         "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hHCalClustParDiff        = new TH1D("hHCalClustParDiff",        "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalClustPosYvsX        = new TH2D("hHCalClustPosYvsX",        "", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hHCalClustVsParEne       = new TH2D("hHCalClustVsParEne",       "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // reconstructed cluster debug histograms
  hHCalDebugClustSum5      = new TH1D("hHCalDebugClustSum5",      "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugClustSum10     = new TH1D("hHCalDebugClustSum10",     "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugClustSum100    = new TH1D("hHCalDebugClustSum100",    "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugClustSum1000   = new TH1D("hHCalDebugClustSum1000",   "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugClustDiff5     = new TH1D("hHCalDebugClustDiff5",     "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalDebugClustDiff10    = new TH1D("hHCalDebugClustDiff10",    "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalDebugClustDiff100   = new TH1D("hHCalDebugClustDiff100",   "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalDebugClustDiff1000  = new TH1D("hHCalDebugClustDiff1000",  "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  // truth cluster histograms
  hHCalTruClustEne           = new TH1D("hHCalTruClustEne",           "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalTruClustPosZ          = new TH1D("hHCalTruClustPosZ",          "", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hHCalTruClustNumHit        = new TH1I("hHCalTruClustNumHit",        "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hHCalTruClustParDiff       = new TH1D("hHCalTruClustParDiff",       "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalTruClustPosYvsX       = new TH2D("hHCalTruClustPosYvsX",       "", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hHCalTruClustVsParEne      = new TH2D("hHCalTruClustVsParEne",      "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // truth cluster debug histograms
  hHCalDebugTruClustSum5     = new TH1D("hHCalDebugTruClustSum5",     "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugTruClustSum10    = new TH1D("hHCalDebugTruClustSum10",    "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugTruClustSum100   = new TH1D("hHCalDebugTruClustSum100",   "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugTruClustSum1000  = new TH1D("hHCalDebugTruClustSum1000",  "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugTruClustDiff5    = new TH1D("hHCalDebugTruClustDiff5",    "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalDebugTruClustDiff10   = new TH1D("hHCalDebugTruClustDiff10",   "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalDebugTruClustDiff100  = new TH1D("hHCalDebugTruClustDiff100",  "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalDebugTruClustDiff1000 = new TH1D("hHCalDebugTruClustDiff1000", "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  // event-wise histograms
  hEvtHCalNumPar             = new TH1I("hEvtHCalNumPar",             "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtHCalNumHit             = new TH1I("hEvtHCalNumHit",             "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtHCalNumClust           = new TH1I("hEvtHCalNumClust",           "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtHCalNumTruClust        = new TH1I("hEvtHCalNumTruClust",        "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtHCalSumHitEne          = new TH1D("hEvtHCalSumHitEne",          "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalSumClustEne        = new TH1D("hEvtHCalSumClustEne",        "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalSumTruClustEne     = new TH1D("hEvtHCalSumTruClustEne",     "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalLeadClustEne       = new TH1D("hEvtHCalLeadClustEne",       "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalLeadTruClustEne    = new TH1D("hEvtHCalLeadTruClustEne",    "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalSumHitDiff         = new TH1D("hEvtHCalSumHitDiff",         "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtHCalSumClustDiff       = new TH1D("hEvtHCalSumClustDiff",       "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtHCalSumTruClustDiff    = new TH1D("hEvtHCalSumTruClustDiff",    "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtHCalLeadClustDiff      = new TH1D("hEvtHCalLeadClustDiff",      "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtHCalLeadTruClustDiff   = new TH1D("hEvtHCalLeadTruClustDiff",   "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtHCalNumClustVsHit      = new TH2I("hEvtHCalNumClustVsHit",      "", nNumBin,   rNumBin[0],   rNumBin[1],   nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtHCalNumTruClustVsClust = new TH2I("hEvtHCalNumTruClustVsClust", "", nNumBin,   rNumBin[0],   rNumBin[1],   nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtHCalSumHitVsPar        = new TH2D("hEvtHCalSumHitVsPar",        "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalSumClustVsPar      = new TH2D("hEvtHCalSumClustVsPar",      "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalSumTruClustVsPar   = new TH2D("hEvtHCalSumTruClustVsPar",   "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalLeadClustVsPar     = new TH2D("hEvtHCalLeadClustVsPar",     "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalLeadTruClustVsPar  = new TH2D("hEvtHCalLeadTruClustVsPar",  "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // errors
  hParChrg               -> Sumw2();
  hParMass               -> Sumw2();
  hParEne                -> Sumw2();
  hParMom                -> Sumw2();
  hParMomX               -> Sumw2();
  hParMomY               -> Sumw2();
  hParMomZ               -> Sumw2();
  hHCalRecHitEne             -> Sumw2();
  hHCalRecHitPosZ            -> Sumw2();
  hHCalRecHitParDiff         -> Sumw2();
  hHCalRecHitPosYvsX         -> Sumw2();
  hHCalRecHitVsParEne        -> Sumw2();
  hHCalClustEne            -> Sumw2();
  hHCalClustPosZ           -> Sumw2();
  hHCalClustNumHit         -> Sumw2();
  hHCalClustParDiff        -> Sumw2();
  hHCalClustPosYvsX        -> Sumw2();
  hHCalClustVsParEne       -> Sumw2();
  hHCalDebugClustSum5      -> Sumw2();
  hHCalDebugClustSum10     -> Sumw2();
  hHCalDebugClustSum100    -> Sumw2();
  hHCalDebugClustSum1000   -> Sumw2();
  hHCalDebugClustDiff5     -> Sumw2();
  hHCalDebugClustDiff10    -> Sumw2();
  hHCalDebugClustDiff100   -> Sumw2();
  hHCalDebugClustDiff1000  -> Sumw2();
  hHCalTruClustEne           -> Sumw2();
  hHCalTruClustPosZ          -> Sumw2();
  hHCalTruClustNumHit        -> Sumw2();
  hHCalTruClustParDiff       -> Sumw2();
  hHCalTruClustPosYvsX       -> Sumw2();
  hHCalTruClustVsParEne      -> Sumw2();
  hHCalDebugTruClustSum5     -> Sumw2();
  hHCalDebugTruClustSum10    -> Sumw2();
  hHCalDebugTruClustSum100   -> Sumw2();
  hHCalDebugTruClustSum1000  -> Sumw2();
  hHCalDebugTruClustDiff5    -> Sumw2();
  hHCalDebugTruClustDiff10   -> Sumw2();
  hHCalDebugTruClustDiff100  -> Sumw2();
  hHCalDebugTruClustDiff1000 -> Sumw2();
  hEvtHCalNumPar             -> Sumw2();
  hEvtHCalNumHit             -> Sumw2();
  hEvtHCalNumClust           -> Sumw2();
  hEvtHCalNumTruClust        -> Sumw2();
  hEvtHCalSumHitEne          -> Sumw2();
  hEvtHCalSumClustEne        -> Sumw2();
  hEvtHCalSumTruClustEne     -> Sumw2();
  hEvtHCalLeadClustEne       -> Sumw2();
  hEvtHCalLeadTruClustEne    -> Sumw2();
  hEvtHCalSumHitDiff         -> Sumw2();
  hEvtHCalSumClustDiff       -> Sumw2();
  hEvtHCalSumTruClustDiff    -> Sumw2();
  hEvtHCalLeadClustDiff      -> Sumw2();
  hEvtHCalLeadTruClustDiff   -> Sumw2();
  hEvtHCalNumClustVsHit      -> Sumw2();
  hEvtHCalNumTruClustVsClust -> Sumw2();
  hEvtHCalSumHitVsPar        -> Sumw2();
  hEvtHCalSumClustVsPar      -> Sumw2();
  hEvtHCalLeadClustVsPar     -> Sumw2();
  hEvtHCalLeadTruClustVsPar  -> Sumw2();
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

  // sum hcal hit energy
  for (auto hit : bhcalRecHits()) {
    eHitSum += hit -> getEnergy();
  }  // end 1st hcal hit loop

  // if hit sum is 0, skip event
  const bool isHCalHitSumNonzero = (eHitSum > 0.);
  if (!isHCalHitSumNonzero) {
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

  // reco. hcal hit loop
  unsigned long nHit(0);
  for (auto hit : bhcalRecHits()) {

    // grab hit properties
    const auto rHitX   = hit -> getPosition().x;
    const auto rHitY   = hit -> getPosition().y;
    const auto rHitZ   = hit -> getPosition().z;
    const auto eHit    = hit -> getEnergy();
    const auto diffHit = (eHit - eMcPar) / eHit;

    // fill hit histograms and increment sums/counters
    hHCalRecHitEne      -> Fill(eHit);
    hHCalRecHitPosZ     -> Fill(rHitZ);
    hHCalRecHitParDiff  -> Fill(diffHit);
    hHCalRecHitPosYvsX  -> Fill(rHitX, rHitY);
    hHCalRecHitVsParEne -> Fill(eMcPar, eHit);
    ++nHit;
  }  // end 2nd hcal hit loop

  // for highest energy clusters
  int    iLeadClust(-1);
  int    iLeadTruClust(-1);
  double eLeadClust(0.);
  double eLeadTruClust(0.);
  double diffLeadClust(0.);
  double diffLeadTruClust(0.);

  // reco. hcal cluster loop
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
    hHCalClustEne      -> Fill(eClust);
    hHCalClustPosZ     -> Fill(rClustZ);
    hHCalClustNumHit   -> Fill(nHitClust);
    hHCalClustParDiff  -> Fill(diffClust);
    hHCalClustPosYvsX  -> Fill(rClustX, rClustY);
    hHCalClustVsParEne -> Fill(eMcPar, eClust);
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
  }  // end reco. hcal cluster loop

  // for debugging reco. hcal clusters
  double eDebugSumClust5(0.);
  double eDebugSumClust10(0.);
  double eDebugSumClust100(0.);
  double eDebugSumClust1000(0.);

  // debug reco. hcal cluster loop
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
  }  // end debug reco. hcal cluster loop

  // true hcal cluster loop
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
    hHCalTruClustEne      -> Fill(eTruClust);
    hHCalTruClustPosZ     -> Fill(rTruClustZ);
    hHCalTruClustNumHit   -> Fill(nHitTruClust);
    hHCalTruClustParDiff  -> Fill(diffTruClust);
    hHCalTruClustPosYvsX  -> Fill(rTruClustX, rTruClustY);
    hHCalTruClustVsParEne -> Fill(eMcPar, eTruClust);
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
  }  // end true hcal cluster loop

  // for debugging truth hcal clusters
  double eDebugSumTruClust5(0.);
  double eDebugSumTruClust10(0.);
  double eDebugSumTruClust100(0.);
  double eDebugSumTruClust1000(0.);

  // debug truth hcal cluster loop
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
      }  // end other true hcal cluster loo

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

  // fill event-wise hcal histograms
  hEvtHCalNumPar             -> Fill(nPar);
  hEvtHCalNumHit             -> Fill(nHit);
  hEvtHCalNumClust           -> Fill(nClust);
  hEvtHCalNumTruClust        -> Fill(nTruClust);
  hEvtHCalSumHitEne          -> Fill(eHitSum);
  hEvtHCalSumClustEne        -> Fill(eClustSum);
  hEvtHCalSumTruClustEne     -> Fill(eTruClustSum);
  hEvtHCalLeadClustEne       -> Fill(eLeadClust);
  hEvtHCalLeadTruClustEne    -> Fill(eLeadTruClust);
  hEvtHCalSumHitDiff         -> Fill(diffHitSum);
  hEvtHCalSumClustDiff       -> Fill(diffClustSum);
  hEvtHCalSumTruClustDiff    -> Fill(diffTruClustSum);
  hEvtHCalLeadClustDiff      -> Fill(diffLeadClust);
  hEvtHCalLeadTruClustDiff   -> Fill(diffLeadTruClust);
  hEvtHCalNumClustVsHit      -> Fill(nHit, nClust);
  hEvtHCalNumTruClustVsClust -> Fill(nClust, nTruClust);
  hEvtHCalSumHitVsPar        -> Fill(eMcPar, eHitSum);
  hEvtHCalSumClustVsPar      -> Fill(eMcPar, eClustSum);
  hEvtHCalSumTruClustVsPar   -> Fill(eMcPar, eTruClustSum);
  hEvtHCalLeadClustVsPar     -> Fill(eMcPar, eLeadClust);
  hEvtHCalLeadTruClustVsPar  -> Fill(eMcPar, eLeadTruClust);

  // do debugging calculations
  const auto diffDebugClustSum5       = (eDebugSumClust5 - eMcPar) / eDebugSumClust5;
  const auto diffDebugClustSum10      = (eDebugSumClust10 - eMcPar) / eDebugSumClust5;
  const auto diffDebugClustSum100     = (eDebugSumClust100 - eMcPar) / eDebugSumClust5;
  const auto diffDebugClustSum1000    = (eDebugSumClust1000 - eMcPar) / eDebugSumClust5;
  const auto diffDebugTruClustSum5    = (eDebugSumTruClust5 - eMcPar) / eDebugSumTruClust5;
  const auto diffDebugTruClustSum10   = (eDebugSumTruClust10 - eMcPar) / eDebugSumTruClust5;
  const auto diffDebugTruClustSum100  = (eDebugSumTruClust100 - eMcPar) / eDebugSumTruClust5;
  const auto diffDebugTruClustSum1000 = (eDebugSumTruClust1000 - eMcPar) / eDebugSumTruClust5;

  // fill reco. cluster hcal debug histograms
  hHCalDebugClustSum5        -> Fill(eDebugSumClust5);
  hHCalDebugClustSum10       -> Fill(eDebugSumClust10);
  hHCalDebugClustSum100      -> Fill(eDebugSumClust100);
  hHCalDebugClustSum1000     -> Fill(eDebugSumClust1000);
  hHCalDebugTruClustSum5     -> Fill(eDebugSumClust5);
  hHCalDebugTruClustSum10    -> Fill(eDebugSumClust10);
  hHCalDebugTruClustSum100   -> Fill(eDebugSumClust100);
  hHCalDebugTruClustSum1000  -> Fill(eDebugSumClust1000);
  hHCalDebugClustDiff5       -> Fill(diffDebugClustSum5);
  hHCalDebugClustDiff10      -> Fill(diffDebugClustSum10);
  hHCalDebugClustDiff100     -> Fill(diffDebugClustSum100);
  hHCalDebugClustDiff1000    -> Fill(diffDebugClustSum1000);
  hHCalDebugTruClustDiff5    -> Fill(diffDebugClustSum5);
  hHCalDebugTruClustDiff10   -> Fill(diffDebugClustSum10);
  hHCalDebugTruClustDiff100  -> Fill(diffDebugClustSum100);
  hHCalDebugTruClustDiff1000 -> Fill(diffDebugClustSum1000);
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
  hParChrg                  -> GetXaxis() -> SetTitle(sCharge.Data());
  hParChrg                  -> GetYaxis() -> SetTitle(sCount.Data());
  hParMass                  -> GetXaxis() -> SetTitle(sMass.Data());
  hParMass                  -> GetYaxis() -> SetTitle(sCount.Data());
  hParEne                   -> GetXaxis() -> SetTitle(sEnePar.Data());
  hParEne                   -> GetYaxis() -> SetTitle(sCount.Data());
  hParMom                   -> GetXaxis() -> SetTitle(sMomPar.Data());
  hParMom                   -> GetYaxis() -> SetTitle(sCount.Data());
  hParMomX                  -> GetXaxis() -> SetTitle(sMomParX.Data());
  hParMomX                  -> GetYaxis() -> SetTitle(sCount.Data());
  hParMomY                  -> GetXaxis() -> SetTitle(sMomParY.Data());
  hParMomY                  -> GetYaxis() -> SetTitle(sCount.Data());
  hParMomZ                  -> GetXaxis() -> SetTitle(sMomParZ.Data());
  hParMomZ                  -> GetYaxis() -> SetTitle(sCount.Data());
  // set reco. hit hcal axis titles
  hHCalRecHitEne            -> GetXaxis() -> SetTitle(sEneHit.Data());
  hHCalRecHitEne            -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalRecHitPosZ           -> GetXaxis() -> SetTitle(sPosHitZ.Data());
  hHCalRecHitPosZ           -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalRecHitParDiff        -> GetXaxis() -> SetTitle(sEneHitDiff.Data());
  hHCalRecHitParDiff        -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalRecHitPosYvsX        -> GetXaxis() -> SetTitle(sPosHitX.Data());
  hHCalRecHitPosYvsX        -> GetYaxis() -> SetTitle(sPosHitY.Data());
  hHCalRecHitPosYvsX        -> GetZaxis() -> SetTitle(sCount.Data());
  hHCalRecHitVsParEne       -> GetXaxis() -> SetTitle(sEnePar.Data());
  hHCalRecHitVsParEne       -> GetYaxis() -> SetTitle(sEneHit.Data());
  hHCalRecHitVsParEne       -> GetZaxis() -> SetTitle(sCount.Data());
  // set reco. cluster hcal axis titles
  hHCalClustEne             -> GetXaxis() -> SetTitle(sEneClust.Data());
  hHCalClustEne             -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustPosZ            -> GetXaxis() -> SetTitle(sPosClustZ.Data());
  hHCalClustPosZ            -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustNumHit          -> GetXaxis() -> SetTitle(sNumHitClust.Data());
  hHCalClustNumHit          -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustParDiff         -> GetXaxis() -> SetTitle(sEneClustDiff.Data());
  hHCalClustParDiff         -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustPosYvsX         -> GetXaxis() -> SetTitle(sPosClustX.Data());
  hHCalClustPosYvsX         -> GetYaxis() -> SetTitle(sPosClustY.Data());
  hHCalClustPosYvsX         -> GetZaxis() -> SetTitle(sCount.Data());
  hHCalClustVsParEne        -> GetXaxis() -> SetTitle(sEnePar.Data());
  hHCalClustVsParEne        -> GetYaxis() -> SetTitle(sEneClust.Data());
  hHCalClustVsParEne        -> GetZaxis() -> SetTitle(sCount.Data());
  // set truth cluster hcal axis titles
  hHCalTruClustEne          -> GetXaxis() -> SetTitle(sEneTruClust.Data());
  hHCalTruClustEne          -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalTruClustPosZ         -> GetXaxis() -> SetTitle(sPosTruClustZ.Data());
  hHCalTruClustPosZ         -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalTruClustNumHit       -> GetXaxis() -> SetTitle(sNumHitTruClust.Data());
  hHCalTruClustNumHit       -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalTruClustParDiff      -> GetXaxis() -> SetTitle(sEneTruClustDiff.Data());
  hHCalTruClustParDiff      -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalTruClustPosYvsX      -> GetXaxis() -> SetTitle(sPosTruClustX.Data());
  hHCalTruClustPosYvsX      -> GetYaxis() -> SetTitle(sPosTruClustY.Data());
  hHCalTruClustPosYvsX      -> GetZaxis() -> SetTitle(sCount.Data());
  hHCalTruClustVsParEne     -> GetXaxis() -> SetTitle(sEnePar.Data());
  hHCalTruClustVsParEne     -> GetYaxis() -> SetTitle(sEneTruClust.Data());
  hHCalTruClustVsParEne     -> GetZaxis() -> SetTitle(sCount.Data());
  // set event-wise hcal axis titles
  hEvtHCalNumPar            -> GetXaxis() -> SetTitle(sNumParEvt.Data());
  hEvtHCalNumPar            -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalNumHit            -> GetXaxis() -> SetTitle(sNumHitEvt.Data());
  hEvtHCalNumHit            -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalNumClust          -> GetXaxis() -> SetTitle(sNumClustEvt.Data());
  hEvtHCalNumClust          -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalNumTruClust       -> GetXaxis() -> SetTitle(sNumTruClustEvt.Data());
  hEvtHCalNumTruClust       -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalSumHitEne         -> GetXaxis() -> SetTitle(sEneHitSum.Data());
  hEvtHCalSumHitEne         -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalSumClustEne       -> GetXaxis() -> SetTitle(sEneClustSum.Data());
  hEvtHCalSumClustEne       -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalSumTruClustEne    -> GetXaxis() -> SetTitle(sEneTruClustSum.Data());
  hEvtHCalSumTruClustEne    -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalLeadClustEne      -> GetXaxis() -> SetTitle(sEneClustLead.Data());
  hEvtHCalLeadClustEne      -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalLeadTruClustEne   -> GetXaxis() -> SetTitle(sEneTruClustLead.Data());
  hEvtHCalLeadTruClustEne   -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalSumHitDiff        -> GetXaxis() -> SetTitle(sEneHitSumDiff.Data());
  hEvtHCalSumHitDiff        -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalSumClustDiff      -> GetXaxis() -> SetTitle(sEneClustSumDiff.Data());
  hEvtHCalSumClustDiff      -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalSumTruClustDiff   -> GetXaxis() -> SetTitle(sEneTruClustSumDiff.Data());
  hEvtHCalSumTruClustDiff   -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalLeadClustDiff     -> GetXaxis() -> SetTitle(sEneClustLeadDiff.Data());
  hEvtHCalLeadClustDiff     -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalLeadTruClustDiff  -> GetXaxis() -> SetTitle(sEneTruClustLeadDiff.Data());
  hEvtHCalLeadTruClustDiff  -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtHCalSumHitVsPar       -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtHCalSumHitVsPar       -> GetYaxis() -> SetTitle(sEneHitSum.Data());
  hEvtHCalSumHitVsPar       -> GetZaxis() -> SetTitle(sCount.Data());
  hEvtHCalSumClustVsPar     -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtHCalSumClustVsPar     -> GetYaxis() -> SetTitle(sEneClustSum.Data());
  hEvtHCalSumClustVsPar     -> GetZaxis() -> SetTitle(sCount.Data());
  hEvtHCalSumTruClustVsPar  -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtHCalSumTruClustVsPar  -> GetYaxis() -> SetTitle(sEneTruClustSum.Data());
  hEvtHCalSumTruClustVsPar  -> GetZaxis() -> SetTitle(sCount.Data());
  hEvtHCalLeadClustVsPar    -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtHCalLeadClustVsPar    -> GetYaxis() -> SetTitle(sEneClustLead.Data());
  hEvtHCalLeadClustVsPar    -> GetZaxis() -> SetTitle(sCount.Data());
  hEvtHCalLeadTruClustVsPar -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtHCalLeadTruClustVsPar -> GetYaxis() -> SetTitle(sEneTruClustLead.Data());
  hEvtHCalLeadTruClustVsPar -> GetZaxis() -> SetTitle(sCount.Data());
  return;

}  // end 'FinishWithGlobalRootLock()'

// end ------------------------------------------------------------------------
