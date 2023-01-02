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
#include <boost/math/special_functions/sign.hpp>


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
  const unsigned long nPhiBin(60);
  const unsigned long nEtaBin(40);
  const unsigned long nEneBin(200);
  const unsigned long nMomBin(200);
  const unsigned long nPosTrBin(800);
  const unsigned long nPosLoBin(30);
  const unsigned long nDiffBin(200);
  const unsigned long rNumBin[NRange]   = {0,      200};
  const double        rChrgBin[NRange]  = {-3.,    3.};
  const double        rMassBin[NRange]  = {0.,     5.};
  const double        rPhiBin[NRange]   = {-3.15,  3.15};
  const double        rEtaBin[NRange]   = {-2.,    2.};
  const double        rEneBin[NRange]   = {0.,     100.};
  const double        rMomBin[NRange]   = {-50.,   50.};
  const double        rPosTrBin[NRange] = {-4000., 4000.};
  const double        rPosLoBin[NRange] = {-3000., 3000.};
  const double        rDiffBin[NRange]  = {-50.,   50.};
  // particle histograms
  hParChrg                   = new TH1D("hParChrg",                   "Gen. Particles",         nChrgBin,  rChrgBin[0],  rChrgBin[1]);
  hParMass                   = new TH1D("hParMass",                   "Gen. Particles",         nMassBin,  rMassBin[0],  rMassBin[1]);
  hParPhi                    = new TH1D("hParPhi",                    "Gen. Particles",         nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  hParEta                    = new TH1D("hParEta",                    "Gen. Particles",         nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hParEne                    = new TH1D("hParEne",                    "Gen. Particles",         nEneBin,   rEneBin[0],   rEneBin[1]);
  hParMom                    = new TH1D("hParMom",                    "Gen. Particles",         nEneBin,   rEneBin[0],   rEneBin[1]);
  hParMomX                   = new TH1D("hParMomX",                   "Gen. Particles",         nMomBin,   rMomBin[0],   rMomBin[1]);
  hParMomY                   = new TH1D("hParMomY",                   "Gen. Particles",         nMomBin,   rMomBin[0],   rMomBin[1]);
  hParMomZ                   = new TH1D("hParMomZ",                   "Gen. Particles",         nMomBin,   rMomBin[0],   rMomBin[1]);
  hParEtaVsPhi               = new TH2D("hParEtaVsPhi",               "Gen. Particles",         nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,  rEtaBin[0],   rEtaBin[1]);
  // reco. hcal hit histograms
  hHCalRecHitPhi             = new TH1D("hHCalRecHitPhi",             "Barrel HCal",            nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  hHCalRecHitEta             = new TH1D("hHCalRecHitEta",             "Barrel HCal",            nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalRecHitEne             = new TH1D("hHCalRecHitEne",             "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalRecHitPosZ            = new TH1D("hHCalRecHitPosZ",            "Barrel HCal",            nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hHCalRecHitParDiff         = new TH1D("hHCalRecHitParDiff",         "Barrel HCal",            nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalRecHitPosYvsX         = new TH2D("hHCalRecHitPosYvsX",         "Barrel HCal",            nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hHCalRecHitEtaVsPhi        = new TH2D("hHCalRecHitEtaVsPhi",        "Barrel HCal",            nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalRecHitVsParEne        = new TH2D("hHCalRecHitVsParEne",        "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // reco. hcal cluster histograms
  hHCalClustPhi              = new TH1D("hHCalClustPhi",              "Barrel HCal",            nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  hHCalClustEta              = new TH1D("hHCalClustEta",              "Barrel HCal",            nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalClustEne              = new TH1D("hHCalClustEne",              "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalClustPosZ             = new TH1D("hHCalClustPosZ",             "Barrel HCal",            nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hHCalClustNumHit           = new TH1I("hHCalClustNumHit",           "Barrel HCal",            nNumBin,   rNumBin[0],   rNumBin[1]);
  hHCalClustParDiff          = new TH1D("hHCalClustParDiff",          "Barrel HCal",            nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalClustPosYvsX          = new TH2D("hHCalClustPosYvsX",          "Barrel HCal",            nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hHCalClustEtaVsPhi         = new TH2D("hHCalClustEtaVsPhi",         "Barrel HCal",            nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalClustVsParEne         = new TH2D("hHCalClustVsParEne",         "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // reco. hcal cluster debug histograms
  hHCalDebugClustSum5        = new TH1D("hHCalDebugClustSum5",        "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugClustSum10       = new TH1D("hHCalDebugClustSum10",       "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugClustSum100      = new TH1D("hHCalDebugClustSum100",      "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugClustSum1000     = new TH1D("hHCalDebugClustSum1000",     "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugClustDiff5       = new TH1D("hHCalDebugClustDiff5",       "Barrel HCal",            nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalDebugClustDiff10      = new TH1D("hHCalDebugClustDiff10",      "Barrel HCal",            nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalDebugClustDiff100     = new TH1D("hHCalDebugClustDiff100",     "Barrel HCal",            nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalDebugClustDiff1000    = new TH1D("hHCalDebugClustDiff1000",    "Barrel HCal",            nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  // truth hcal cluster histograms
  hHCalTruClustPhi           = new TH1D("hHCalTruClustPhi",           "Barrel HCal",            nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  hHCalTruClustEta           = new TH1D("hHCalTruClustEta",           "Barrel HCal",            nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalTruClustEne           = new TH1D("hHCalTruClustEne",           "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalTruClustPosZ          = new TH1D("hHCalTruClustPosZ",          "Barrel HCal",            nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hHCalTruClustNumHit        = new TH1I("hHCalTruClustNumHit",        "Barrel HCal",            nNumBin,   rNumBin[0],   rNumBin[1]);
  hHCalTruClustParDiff       = new TH1D("hHCalTruClustParDiff",       "Barrel HCal",            nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalTruClustPosYvsX       = new TH2D("hHCalTruClustPosYvsX",       "Barrel HCal",            nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hHCalTruClustEtaVsPhi      = new TH2D("hHCalTruClustEtaVsPhi",      "Barrel HCal",            nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalTruClustVsParEne      = new TH2D("hHCalTruClustVsParEne",      "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // truth hcal cluster debug histograms
  hHCalDebugTruClustSum5     = new TH1D("hHCalDebugTruClustSum5",     "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugTruClustSum10    = new TH1D("hHCalDebugTruClustSum10",    "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugTruClustSum100   = new TH1D("hHCalDebugTruClustSum100",   "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugTruClustSum1000  = new TH1D("hHCalDebugTruClustSum1000",  "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugTruClustDiff5    = new TH1D("hHCalDebugTruClustDiff5",    "Barrel HCal",            nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalDebugTruClustDiff10   = new TH1D("hHCalDebugTruClustDiff10",   "Barrel HCal",            nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalDebugTruClustDiff100  = new TH1D("hHCalDebugTruClustDiff100",  "Barrel HCal",            nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalDebugTruClustDiff1000 = new TH1D("hHCalDebugTruClustDiff1000", "Barrel HCal",            nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  // event-wise hcal histograms
  hEvtHCalNumPar             = new TH1I("hEvtHCalNumPar",             "Barrel HCal",            nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtHCalNumHit             = new TH1I("hEvtHCalNumHit",             "Barrel HCal",            nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtHCalNumClust           = new TH1I("hEvtHCalNumClust",           "Barrel HCal",            nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtHCalNumTruClust        = new TH1I("hEvtHCalNumTruClust",        "Barrel HCal",            nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtHCalSumHitEne          = new TH1D("hEvtHCalSumHitEne",          "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalSumClustEne        = new TH1D("hEvtHCalSumClustEne",        "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalSumTruClustEne     = new TH1D("hEvtHCalSumTruClustEne",     "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalLeadClustEne       = new TH1D("hEvtHCalLeadClustEne",       "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalLeadTruClustEne    = new TH1D("hEvtHCalLeadTruClustEne",    "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalSumHitDiff         = new TH1D("hEvtHCalSumHitDiff",         "Barrel HCal",            nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtHCalSumClustDiff       = new TH1D("hEvtHCalSumClustDiff",       "Barrel HCal",            nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtHCalSumTruClustDiff    = new TH1D("hEvtHCalSumTruClustDiff",    "Barrel HCal",            nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtHCalLeadClustDiff      = new TH1D("hEvtHCalLeadClustDiff",      "Barrel HCal",            nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtHCalLeadTruClustDiff   = new TH1D("hEvtHCalLeadTruClustDiff",   "Barrel HCal",            nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtHCalNumClustVsHit      = new TH2I("hEvtHCalNumClustVsHit",      "Barrel HCal",            nNumBin,   rNumBin[0],   rNumBin[1],   nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtHCalNumTruClustVsClust = new TH2I("hEvtHCalNumTruClustVsClust", "Barrel HCal",            nNumBin,   rNumBin[0],   rNumBin[1],   nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtHCalSumHitVsPar        = new TH2D("hEvtHCalSumHitVsPar",        "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalSumClustVsPar      = new TH2D("hEvtHCalSumClustVsPar",      "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalSumTruClustVsPar   = new TH2D("hEvtHCalSumTruClustVsPar",   "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalLeadClustVsPar     = new TH2D("hEvtHCalLeadClustVsPar",     "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalLeadTruClustVsPar  = new TH2D("hEvtHCalLeadTruClustVsPar",  "Barrel HCal",            nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // errors
  hParChrg                   -> Sumw2();
  hParMass                   -> Sumw2();
  hParPhi                    -> Sumw2();
  hParEta                    -> Sumw2();
  hParEne                    -> Sumw2();
  hParMom                    -> Sumw2();
  hParMomX                   -> Sumw2();
  hParMomY                   -> Sumw2();
  hParMomZ                   -> Sumw2();
  hParEtaVsPhi               -> Sumw2();
  hHCalRecHitPhi             -> Sumw2();
  hHCalRecHitEta             -> Sumw2();
  hHCalRecHitEne             -> Sumw2();
  hHCalRecHitPosZ            -> Sumw2();
  hHCalRecHitParDiff         -> Sumw2();
  hHCalRecHitPosYvsX         -> Sumw2();
  hHCalRecHitEtaVsPhi        -> Sumw2();
  hHCalRecHitVsParEne        -> Sumw2();
  hHCalClustPhi              -> Sumw2();
  hHCalClustEta              -> Sumw2();
  hHCalClustEne              -> Sumw2();
  hHCalClustPosZ             -> Sumw2();
  hHCalClustNumHit           -> Sumw2();
  hHCalClustParDiff          -> Sumw2();
  hHCalClustPosYvsX          -> Sumw2();
  hHCalClustEtaVsPhi         -> Sumw2();
  hHCalClustVsParEne         -> Sumw2();
  hHCalDebugClustSum5        -> Sumw2();
  hHCalDebugClustSum10       -> Sumw2();
  hHCalDebugClustSum100      -> Sumw2();
  hHCalDebugClustSum1000     -> Sumw2();
  hHCalDebugClustDiff5       -> Sumw2();
  hHCalDebugClustDiff10      -> Sumw2();
  hHCalDebugClustDiff100     -> Sumw2();
  hHCalDebugClustDiff1000    -> Sumw2();
  hHCalTruClustPhi           -> Sumw2();
  hHCalTruClustEta           -> Sumw2();
  hHCalTruClustEne           -> Sumw2();
  hHCalTruClustPosZ          -> Sumw2();
  hHCalTruClustNumHit        -> Sumw2();
  hHCalTruClustParDiff       -> Sumw2();
  hHCalTruClustPosYvsX       -> Sumw2();
  hHCalTruClustEtaVsPhi      -> Sumw2();
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
  double eHCalHitSum(0.);
  double eECalHitSum(0.);
  double eHCalClustSum(0.);
  double eECalClustSum(0.);
  double eTruHCalClustSum(0.);
  double eTruECalClustSum(0.);

  // sum hcal hit energy
  for (auto bhCalHit : bhcalRecHits()) {
    eHCalHitSum += bhCalHit -> getEnergy();
  }  // end 1st hcal hit loop

  // if both hit sums are 0, skip event
  const bool isHCalHitSumNonzero = (eHCalHitSum > 0.);
  if (!isHCalHitSumNonzero) {
    return;
  }

  // MC particle properties
  float  cMcPar(0.);
  double mMcPar(0.);
  double fMcPar(0.);
  double hMcPar(0.);
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
    const auto fPar  = std::atan(pParY / pParX);
    const auto hPar  = std::atanh(pParZ / pPar);

    // select MC particle
    const bool isRightCharge   = (cPar == CPar);
    const bool isRightMass     = ((mPar >= MParMin) && (mPar <= MParMax));
    const bool isRightMomentum = ((pPar >= PParMin) && (pPar <= PParMax));
    const bool isMcParticle    = (isRightCharge && isRightMass && isRightMomentum);
    if (isMcParticle) {
      cMcPar    = cPar;
      mMcPar    = mPar;
      fMcPar    = fPar;
      hMcPar    = hPar;
      eMcPar    = ePar;
      pMcPar[0] = pParX;
      pMcPar[1] = pParY;
      pMcPar[2] = pParZ;
      pTotMcPar = pPar;
    }
    ++nPar;
  }  // end particle loop

  // fill particle histograms
  hParChrg     -> Fill(cMcPar);
  hParMass     -> Fill(mMcPar);
  hParPhi      -> Fill(fMcPar);
  hParEta      -> Fill(hMcPar);
  hParEne      -> Fill(eMcPar);
  hParMom      -> Fill(pTotMcPar);
  hParMomX     -> Fill(pMcPar[0]);
  hParMomY     -> Fill(pMcPar[1]);
  hParMomZ     -> Fill(pMcPar[2]);
  hParEtaVsPhi -> Fill(fMcPar, hMcPar);

  // reco. hcal hit loop
  unsigned long nHCalHit(0);
  for (auto bhCalHit : bhcalRecHits()) {

    // grab hit properties
    const auto rHCalHitX   = bhCalHit -> getPosition().x;
    const auto rHCalHitY   = bhCalHit -> getPosition().y;
    const auto rHCalHitZ   = bhCalHit -> getPosition().z;
    const auto eHCalHit    = bhCalHit -> getEnergy();
    const auto rHCalHitS   = std::sqrt((rHCalHitX * rHCalHitX) + (rHCalHitY * rHCalHitY));
    const auto rHCalHitR   = std::sqrt((rHCalHitS * rHCalHitS) + (rHCalHitZ * rHCalHitZ));
    const auto fHCalHit    = boost::math::sign(rHCalHitY) * acos(rHCalHitX / rHCalHitS);
    const auto tHCalHit    = std::acos(rHCalHitZ / rHCalHitR);
    const auto hHCalHit    = (-1.) * std::log(std::atan(tHCalHit / 2.));
    const auto diffHCalHit = (eHCalHit - eMcPar) / eHCalHit;

    // fill hit histograms and increment sums/counters
    hHCalRecHitPhi      -> Fill(fHCalHit);
    hHCalRecHitEta      -> Fill(hHCalHit);
    hHCalRecHitEne      -> Fill(eHCalHit);
    hHCalRecHitPosZ     -> Fill(rHCalHitZ);
    hHCalRecHitParDiff  -> Fill(diffHCalHit);
    hHCalRecHitPosYvsX  -> Fill(rHCalHitX, rHCalHitY);
    hHCalRecHitEtaVsPhi -> Fill(fHCalHit, hHCalHit);
    hHCalRecHitVsParEne -> Fill(eMcPar, eHCalHit);
    ++nHCalHit;
  }  // end 2nd hcal hit loop

  // for highest energy clusters
  int    iLeadHCalClust(-1);
  int    iLeadTruHCalClust(-1);
  double eLeadHCalClust(0.);
  double eLeadTruHCalClust(0.);
  double diffLeadHCalClust(0.);
  double diffLeadTruHCalClust(0.);

  // reco. hcal cluster loop
  unsigned long iHCalClust(0);
  unsigned long nHCalClust(0);
  for (auto bhCalClust : bhcalClusters()) {

    // grab cluster properties
    const auto rHCalClustX   = bhCalClust -> getPosition().x;
    const auto rHCalClustY   = bhCalClust -> getPosition().y;
    const auto rHCalClustZ   = bhCalClust -> getPosition().z;
    const auto eHCalClust    = bhCalClust -> getEnergy();
    const auto nHitHCalClust = bhCalClust -> getNhits();
    const auto fHCalClust    = bhCalClust -> getIntrinsicPhi();
    const auto tHCalClust    = bhCalClust -> getIntrinsicTheta();
    const auto hHCalClust    = (1.) * std::log(std::atan(tHCalcClust / 2.));
    const auto diffHCalClust = (eHCalClust - eMcPar) / eHCalClust;

    // fill cluster histograms and increment counters
    hHCalClustPhi      -> Fill(fHCalClust);
    hHCalClustEta      -> Fill(hHCalClust);
    hHCalClustEne      -> Fill(eHCalClust);
    hHCalClustPosZ     -> Fill(rHCalClustZ);
    hHCalClustNumHit   -> Fill(nHitHCalClust);
    hHCalClustParDiff  -> Fill(diffHCalClust);
    hHCalClustPosYvsX  -> Fill(rHCalClustX, rHCalClustY);
    hHCalClustEtaVsPhi -> Fill(fHCalClust, hHCalClust);
    hHCalClustVsParEne -> Fill(eMcPar, eHCalClust);
    eHCalClustSum += eHCalClust;
    ++nHCalClust;

    // select leading cluster
    const bool isBiggerEne = (eHCalClust > eLeadHCalClust);
    if (isBiggerEne) {
      iLeadHCalClust    = iHCalClust;
      eLeadHCalClust    = eHCalClust;
      diffLeadHCalClust = diffHCalClust;
    }
    ++iHCalClust;
  }  // end reco. hcal cluster loop

  // for debugging reco. clusters
  double eDebugSumHCalClust5(eLeadHCalClust);
  double eDebugSumHCalClust10(eLeadHCalClust);
  double eDebugSumHCalClust100(eLeadHCalClust);
  double eDebugSumHCalClust1000(eLeadHCalClust);

  // debug reco. hcal cluster loop
  unsigned long iDebugHCalClust(0);
  for (auto debugHCalClust : bhcalClusters()) {

    // select leading cluster
    const bool isLeadHCalClust = (iDebugHCalClust == iLeadHCalClust);
    if (isLeadHCalClust) {

      // grab lead cluster properties
      const auto rLeadHCalClustX = debugHCalClust -> getPosition().x;
      const auto rLeadHCalClustY = debugHCalClust -> getPosition().y;
      const auto eLeadHCalClust  = debugHCalClust -> getEnergy();

      unsigned long iOtherHCalClust(0);
      double        toAddHCalClust5(0.);
      double        toAddHCalClust10(0.);
      double        toAddHCalClust100(0.);
      double        toAddHCalClust1000(0.);
      for (auto otherHCalClust : bhcalClusters()) {

        // ignore same cluster
        const bool isSameHCalClust = (iOtherHCalClust == iDebugHCalClust);
        if (isSameHCalClust) continue;

        // grab other cluster properties
        const auto rOtherHCalClustX = otherHCalClust -> getPosition().x;
        const auto rOtherHCalClustY = otherHCalClust -> getPosition().y;
        const auto eOtherHCalClust  = otherHCalClust -> getEnergy();
        const auto drLeadOtherX     = rOtherHCalClustX - rLeadHCalClustX;
        const auto drLeadOtherY     = rOtherHCalClustY - rLeadHCalClustY;
        const auto drLeadOther      = std::sqrt((drLeadOtherX * drLeadOtherX) + (drLeadOtherY * drLeadOtherY));

        // increment relevant sums and counters
        const bool isIn5mm    = (drLeadOther < 5.);
        const bool isIn10mm   = (drLeadOther < 10.);
        const bool isIn100mm  = (drLeadOther < 100.);
        const bool isIn1000mm = (drLeadOther < 1000.);
        if (isIn5mm)    toAddHCalClust5    += eOtherHCalClust;
        if (isIn10mm)   toAddHCalClust10   += eOtherHCalClust;
        if (isIn100mm)  toAddHCalClust100  += eOtherHCalClust;
        if (isIn1000mm) toAddHCalClust1000 += eOtherHCalClust;
        ++iOtherHCalClust;
      }  // end other reco. cluster loo

      // add sums to lead energy
      eDebugSumHCalClust5    += toAddHCalClust5;
      eDebugSumHCalClust10   += toAddHCalClust10;
      eDebugSumHCalClust100  += toAddHCalClust100;
      eDebugSumHCalClust1000 += toAddHCalClust1000;
    }
    ++iDebugHCalClust;
  }  // end debug reco. hcal cluster loop

  // true hcal cluster loop
  unsigned long iTruHCalClust(0);
  unsigned long nTruHCalClust(0);
  for (auto truthHCalClust : bhcalTruthClusters()) {

    // grab cluster properties
    const auto rTruHCalClustX   = truthHCalClust -> getPosition().x;
    const auto rTruHCalClustY   = truthHCalClust -> getPosition().y;
    const auto rTruHCalClustZ   = truthHCalClust -> getPosition().z;
    const auto eTruHCalClust    = truthHCalClust -> getEnergy();
    const auto nHitTruHCalClust = truthHCalClust -> getNhits();
    const auto fTruHCalClust    = truthHCalClust -> getIntrinsicPhi();
    const auto tTruHCalClust    = truthHCalClust -> getIntrinsicTheta();
    const auto hTruHCalClust    = (-1.) * std::log(std::atan(tTruHCalClust / 2.));
    const auto diffTruHCalClust = (eTruHCalClust - eMcPar) / eTruHCalClust;

    // fill cluster histograms and increment counters
    hHCalTruClustPhi      -> Fill(fTruHCalClust);
    hHCalTruClustEta      -> Fill(hTruHCalClust);
    hHCalTruClustEne      -> Fill(eTruHCalClust);
    hHCalTruClustPosZ     -> Fill(rTruHCalClustZ);
    hHCalTruClustNumHit   -> Fill(nHitTruHCalClust);
    hHCalTruClustParDiff  -> Fill(diffTruHCalClust);
    hHCalTruClustPosYvsX  -> Fill(rTruHCalClustX, rTruHCalClustY);
    hHCalTruClustEtaVsPhi -> Fill(fTruHCalClust, hTruHCalClust);
    hHCalTruClustVsParEne -> Fill(eMcPar, eTruHCalClust);
    eTruHCalClustSum += eTruHCalClust;
    ++nTruHCalClust;

    // select leading cluster
    const bool isBiggerEne = (eTruHCalClust > eLeadTruHCalClust);
    if (isBiggerEne) {
      iLeadTruHCalClust    = iTruHCalClust;
      eLeadTruHCalClust    = eTruHCalClust;
      diffLeadTruHCalClust = diffTruHCalClust;
    }
    ++iTruHCalClust;
  }  // end true hcal cluster loop

  // for debugging truth clusters
  double eDebugSumTruHCalClust5(eLeadTruHCalClust);
  double eDebugSumTruHCalClust10(eLeadTruHCalClust);
  double eDebugSumTruHCalClust100(eLeadTruHCalClust);
  double eDebugSumTruHCalClust1000(eLeadTruHCalClust);

  // debug truth hcal cluster loop
  unsigned long iDebugTruHCalClust(0);
  for (auto debugTruthHCalClust : bhcalTruthClusters()) {

    // select leading cluster
    const bool isLeadTruHCalClust = (iDebugTruHCalClust == iLeadTruHCalClust);
    if (isLeadTruHCalClust) {

      // grab lead cluster properties
      const auto rLeadTruHCalClustX = debugTruthHCalClust -> getPosition().x;
      const auto rLeadTruHCalClustY = debugTruthHCalClust -> getPosition().y;
      const auto eLeadTruHCalClust  = debugTruthHCalClust -> getEnergy();

      unsigned long iOtherTruHCalClust(0);
      double        toAddTruHCalClust5(0.);
      double        toAddTruHCalClust10(0.);
      double        toAddTruHCalClust100(0.);
      double        toAddTruHCalClust1000(0.);
      for (auto otherTruthHCalClust : bhcalTruthClusters()) {

        // ignore same cluster
        const bool isSameTruHCalClust = (iOtherTruHCalClust == iDebugTruHCalClust);
        if (isSameTruHCalClust) continue;

        // grab other cluster properties
        const auto rOtherTruHCalClustX = otherTruthHCalClust -> getPosition().x;
        const auto rOtherTruHCalClustY = otherTruthHCalClust -> getPosition().y;
        const auto eOtherTruHCalClust  = otherTruthHCalClust -> getEnergy();
        const auto drLeadOtherX        = rOtherTruHCalClustX - rLeadTruHCalClustX;
        const auto drLeadOtherY        = rOtherTruHCalClustY - rLeadTruHCalClustY;
        const auto drLeadOther         = std::sqrt((drLeadOtherX * drLeadOtherX) + (drLeadOtherY * drLeadOtherY));

        // increment relevant sums and counters
        const bool isIn5mm    = (drLeadOther < 5);
        const bool isIn10mm   = (drLeadOther < 10);
        const bool isIn100mm  = (drLeadOther < 100);
        const bool isIn1000mm = (drLeadOther < 1000);
        if (isIn5mm)    toAddTruHCalClust5    += eOtherTruHCalClust;
        if (isIn10mm)   toAddTruHCalClust10   += eOtherTruHCalClust;
        if (isIn100mm)  toAddTruHCalClust100  += eOtherTruHCalClust;
        if (isIn1000mm) toAddTruHCalClust1000 += eOtherTruHCalClust;
        ++iOtherTruHCalClust;
      }  // end other true hcal cluster loo

      // add sums to lead energy
      eDebugSumTruHCalClust5    += toAddTruHCalClust5;
      eDebugSumTruHCalClust10   += toAddTruHCalClust10;
      eDebugSumTruHCalClust100  += toAddTruHCalClust100;
      eDebugSumTruHCalClust1000 += toAddTruHCalClust1000;
    }
    ++iDebugTruHCalClust;
  }  // end debug truth hcal cluster loop

  // do event-wise calculations
  const auto diffHCalHitSum      = (eHCalHitSum - eMcPar) / eHCalHitSum;
  const auto diffHCalClustSum    = (eHCalClustSum - eMcPar) / eHCalClustSum;
  const auto diffTruHCalClustSum = (eTruHCalClustSum - eMcPar) / eTruHCalClustSum;

  // fill event-wise hcal histograms
  hEvtHCalNumPar             -> Fill(nPar);
  hEvtHCalNumHit             -> Fill(nHCalHit);
  hEvtHCalNumClust           -> Fill(nHCalClust);
  hEvtHCalNumTruClust        -> Fill(nTruHCalClust);
  hEvtHCalSumHitEne          -> Fill(eHCalHitSum);
  hEvtHCalSumClustEne        -> Fill(eHCalClustSum);
  hEvtHCalSumTruClustEne     -> Fill(eTruHCalClustSum);
  hEvtHCalLeadClustEne       -> Fill(eLeadHCalClust);
  hEvtHCalLeadTruClustEne    -> Fill(eLeadTruHCalClust);
  hEvtHCalSumHitDiff         -> Fill(diffHCalHitSum);
  hEvtHCalSumClustDiff       -> Fill(diffHCalClustSum);
  hEvtHCalSumTruClustDiff    -> Fill(diffTruHCalClustSum);
  hEvtHCalLeadClustDiff      -> Fill(diffLeadHCalClust);
  hEvtHCalLeadTruClustDiff   -> Fill(diffLeadTruHCalClust);
  hEvtHCalNumClustVsHit      -> Fill(nHCalHit, nHCalClust);
  hEvtHCalNumTruClustVsClust -> Fill(nHCalClust, nTruHCalClust);
  hEvtHCalSumHitVsPar        -> Fill(eMcPar, eHCalHitSum);
  hEvtHCalSumClustVsPar      -> Fill(eMcPar, eHCalClustSum);
  hEvtHCalSumTruClustVsPar   -> Fill(eMcPar, eTruHCalClustSum);
  hEvtHCalLeadClustVsPar     -> Fill(eMcPar, eLeadHCalClust);
  hEvtHCalLeadTruClustVsPar  -> Fill(eMcPar, eLeadTruHCalClust);

  // do debugging calculations
  const auto diffDebugHCalClustSum5       = (eDebugSumHCalClust5 - eMcPar) / eDebugSumHCalClust5;
  const auto diffDebugHCalClustSum10      = (eDebugSumHCalClust10 - eMcPar) / eDebugSumHCalClust5;
  const auto diffDebugHCalClustSum100     = (eDebugSumHCalClust100 - eMcPar) / eDebugSumHCalClust5;
  const auto diffDebugHCalClustSum1000    = (eDebugSumHCalClust1000 - eMcPar) / eDebugSumHCalClust5;
  const auto diffDebugTruHCalClustSum5    = (eDebugSumTruHCalClust5 - eMcPar) / eDebugSumTruHCalClust5;
  const auto diffDebugTruHCalClustSum10   = (eDebugSumTruHCalClust10 - eMcPar) / eDebugSumTruHCalClust5;
  const auto diffDebugTruHCalClustSum100  = (eDebugSumTruHCalClust100 - eMcPar) / eDebugSumTruHCalClust5;
  const auto diffDebugTruHCalClustSum1000 = (eDebugSumTruHCalClust1000 - eMcPar) / eDebugSumTruHCalClust5;

  // fill cluster hcal debug histograms
  hHCalDebugClustSum5        -> Fill(eDebugSumHCalClust5);
  hHCalDebugClustSum10       -> Fill(eDebugSumHCalClust10);
  hHCalDebugClustSum100      -> Fill(eDebugSumHCalClust100);
  hHCalDebugClustSum1000     -> Fill(eDebugSumHCalClust1000);
  hHCalDebugTruClustSum5     -> Fill(eDebugSumHCalClust5);
  hHCalDebugTruClustSum10    -> Fill(eDebugSumHCalClust10);
  hHCalDebugTruClustSum100   -> Fill(eDebugSumHCalClust100);
  hHCalDebugTruClustSum1000  -> Fill(eDebugSumHCalClust1000);
  hHCalDebugClustDiff5       -> Fill(diffDebugHCalClustSum5);
  hHCalDebugClustDiff10      -> Fill(diffDebugHCalClustSum10);
  hHCalDebugClustDiff100     -> Fill(diffDebugHCalClustSum100);
  hHCalDebugClustDiff1000    -> Fill(diffDebugHCalClustSum1000);
  hHCalDebugTruClustDiff5    -> Fill(diffDebugTruHCalClustSum5);
  hHCalDebugTruClustDiff10   -> Fill(diffDebugTruHCalClustSum10);
  hHCalDebugTruClustDiff100  -> Fill(diffDebugTruHCalClustSum100);
  hHCalDebugTruClustDiff1000 -> Fill(diffDebugTruHCalClustSum1000);
  return;

}  // end 'ProcessSequential(std::shared_ptr<JEvent>&)'



//-------------------------------------------
// FinishWithGlobalRootLock
//-------------------------------------------
void JCalibrateHCalProcessor::FinishWithGlobalRootLock() {

  // generic axis titles
  const TString sCount("counts");

  // particle axis titles
  const TString sMass("m_{par} [GeV/c^{2}]");
  const TString sCharge("charge");
  const TString sPhiPar("#varphi_{par}");
  const TString sEtaPar("#eta_{Par}");
  const TString sEnePar("E_{par} [GeV]");
  const TString sMomPar("p_{par} [GeV/c]");
  const TString sMomParX("p_{x, par} [GeV/c]");
  const TString sMomParY("p_{y, par} [GeV/c]");
  const TString sMomParZ("p_{z, par} [GeV/c]");
  const TString sNumParEvt("N_{par} per event");

  // hit axis titles
  const TString sPosHitX("x_{hit} [mm]");
  const TString sPosHitY("y_{hit} [mm]");
  const TString sPosHitZ("z_{hit} [mm]");
  const TString sPhiHit("#varphi_{hit}");
  const TString sEtaHit("#eta_{hit}");
  const TString sEneHit("e_{hit} [GeV]");
  const TString sEneHitSum("E^{sum}_{hit} = #Sigmae_{hit} [GeV]");
  const TString sEneHitDiff("#Deltae_{hit} / e_{hit} = (e_{hit} - E_{par}) / e_{hit} [GeV]");
  const TString sEneHitSumDiff("#DeltaE^{sum}_{hit} / E^{sum}_{hit} = (E^{sum}_{hit} - E_{par}) / E^{sum}_{hit} [GeV]");
  const TString sNumHitEvt("N_{hit} per event");

  // reco. cluster axis titles
  const TString sPosClustX("x_{clust} [mm]");
  const TString sPosClustY("y_{clust} [mm]");
  const TString sPosClustZ("z_{clust} [mm]");
  const TString sEneClust("e_{clust} [GeV]");
  const TString sPhiClust("#varphi_{clust}");
  const TString sEtaClust("#eta_{clust}");
  const TString sEneClustSum("E^{sum}_{clust} = #Sigmae_{clust} [GeV]");
  const TString sEneClustDiff("#Deltae_{clust} / e_{clust} = (e_{clust} - E_{par}) / e_{clust} [GeV]");
  const TString sEneClustLead("E^{lead}_{clust} [GeV]");
  const TString sEneClustSumDiff("#DeltaE^{sum}_{clust} / E^{sum}_{clust} = (E^{sum}_{clust} - E_{par}) / E^{sum}_{clust} [GeV]");
  const TString sEneClustLeadDiff("#DeltaE^{lead}_{clust} / E^{lead}_{clust} = (E^{lead}_{clust} - E_{par}) / E^{lead}_{clust} [GeV]");
  const TString sNumHitClust("N_{hit} per cluster");
  const TString sNumClustEvt("N_{clust} per event");

  // truth cluster axis titles
  const TString sPosTruClustX("x_{truth clust} [mm]");
  const TString sPosTruClustY("y_{truth clust} [mm]");
  const TString sPosTruClustZ("z_{truth clust} [mm]");
  const TString sPhiTruClust("#varphi^{truth}_{clust}");
  const TString sEtaTruClust("#eta^{truth}_{clust}");
  const TString sEneTruClust("e^{truth}_{clust} [GeV]");
  const TString sEneTruClustDiff("#Deltae^{truth}_{clust} / e^{truth}_{clust} / (e^{truth}_{clust} - E_{par}) / e^{truth}_{clust} [GeV]");
  const TString sEneTruClustSum("E^{sum/truth}_{clust} = #Sigmae^{truth}_{clust} [GeV]");
  const TString sEneTruClustLead("E^{lead/truth}_{clust} [GeV]");
  const TString sEneTruClustSumDiff("#DeltaE^{sum/truth}_{clust} / E^{sum/truth}_{clust} = (E^{sum/truth}_{clust} - E_{par}) / E^{sum/truth}_{clust} [GeV]");
  const TString sEneTruClustLeadDiff("#DeltaE^{lead/truth}_{clust} / E^{lead/truth}_{clust} = (E^{lead/truth} _{clust} - E_{par}) / E^{lead/truth}_{clust} [GeV]");
  const TString sNumHitTruClust("N_{hit} per truth cluster");
  const TString sNumTruClustEvt("N_{truth clust} per event");

  // set particle axis titles
  hParChrg                  -> GetXaxis() -> SetTitle(sCharge.Data());
  hParChrg                  -> GetYaxis() -> SetTitle(sCount.Data());
  hParMass                  -> GetXaxis() -> SetTitle(sMass.Data());
  hParMass                  -> GetYaxis() -> SetTitle(sCount.Data());
  hParPhi                   -> GetXaxis() -> SetTitle(sPhiPar.Data());
  hParPhi                   -> GetYaxis() -> SetTitle(sCount.Data());
  hParEta                   -> GetXaxis() -> SetTitle(sEtaPar.Data());
  hParEta                   -> GetYaxis() -> SetTitle(sCount.Data());
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
  hParEtaVsPhi              -> GetXaxis() -> SetTitle(sPhiPar.Data());
  hParEtaVsPhi              -> GetYaxis() -> SetTitle(sEtaPar.Data());
  hParEtaVsPhi              -> GetZaxis() -> SetTitle(sCount.Data());
  // set reco. hit hcal axis titles
  hHCalRecHitPhi            -> GetXaxis() -> SetTitle(sPhiHit.Data());
  hHCalRecHitPhi            -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalRecHitEta            -> GetXaxis() -> SetTitle(sEtaHit.Data());
  hHCalRecHitEta            -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalRecHitEne            -> GetXaxis() -> SetTitle(sEneHit.Data());
  hHCalRecHitEne            -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalRecHitPosZ           -> GetXaxis() -> SetTitle(sPosHitZ.Data());
  hHCalRecHitPosZ           -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalRecHitParDiff        -> GetXaxis() -> SetTitle(sEneHitDiff.Data());
  hHCalRecHitParDiff        -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalRecHitPosYvsX        -> GetXaxis() -> SetTitle(sPosHitX.Data());
  hHCalRecHitPosYvsX        -> GetYaxis() -> SetTitle(sPosHitY.Data());
  hHCalRecHitPosYvsX        -> GetZaxis() -> SetTitle(sCount.Data());
  hHCalRecHitEtaVsPhi       -> GetXaxis() -> SetTitle(sPhiHit.Data());
  hHCalRecHitEtaVsPhi       -> GetYaxis() -> SetTitle(sEtaHit.Data());
  hHCalRecHitEtaVsPhi       -> GetZaxis() -> SetTitle(sCount.Data());
  hHCalRecHitVsParEne       -> GetXaxis() -> SetTitle(sEnePar.Data());
  hHCalRecHitVsParEne       -> GetYaxis() -> SetTitle(sEneHit.Data());
  hHCalRecHitVsParEne       -> GetZaxis() -> SetTitle(sCount.Data());
  // set reco. cluster hcal axis titles
  hHCalClustPhi             -> GetXaxis() -> SetTitle(sPhiClust.Data());
  hHCalClustPhi             -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustEta             -> GetXaxis() -> SetTitle(sEtaClust.Data());
  hHCalClustEta             -> GetYaxis() -> SetTitle(sCount.Data());
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
  hHCalClustEtaVsPhi        -> GetXaxis() -> SetTitle(sPhiClust.Data());
  hHCalClustEtaVsPhi        -> GetYaxis() -> SetTitle(sEtaClust.Data());
  hHCalClustEtaVsPhi        -> GetZaxis() -> SetTitle(sCount.Data());
  hHCalClustVsParEne        -> GetXaxis() -> SetTitle(sEnePar.Data());
  hHCalClustVsParEne        -> GetYaxis() -> SetTitle(sEneClust.Data());
  hHCalClustVsParEne        -> GetZaxis() -> SetTitle(sCount.Data());
  // set truth cluster hcal axis titles
  hHCalTruClustPhi          -> GetXaxis() -> SetTitle(sPhiTruClust.Data());
  hHCalTruClustPhi          -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalTruClustEta          -> GetXaxis() -> SetTitle(sEtaTruClust.Data());
  hHCalTruClustEta          -> GetYaxis() -> SetTitle(sCount.Data());
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
  hHCalTruClustEtaVsPhi     -> GetXaxis() -> SetTitle(sPhiTruClust.Data());
  hHCalTruClustEtaVsPhi     -> GetYaxis() -> SetTitle(sEtaTruClust.Data());
  hHCalTruClustEtaVsPhi     -> GetZaxis() -> SetTitle(sCount.Data());
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

  // debug sum axis titles
  const TString sDebugClustSum5("E^{sum}_{clust}(5 mm) = E^{lead}_{clust} + #Sigma_{#Deltar < 5 mm} e_{clust}");
  const TString sDebugClustSum10("E^{sum}_{clust}(10 mm) = E^{lead}_{clust} + #Sigma_{#Deltar < 10 mm} e_{clust}");
  const TString sDebugClustSum100("E^{sum}_{clust}(1 cm) = E^{lead}_{clust} + #Sigma_{#Deltar < 1 cm} e_{clust}");
  const TString sDebugClustSum1000("E^{sum}_{clust}(1 m) = E^{lead}_{clust} + #Sigma_{#Deltar < 1 m} e_{clust}");
  const TString sDebugTruClustSum5("E^{sum/truth}_{clust}(5 mm) = E^{lead}_{clust} + #Sigma_{#Deltar < 5 mm} e_{clust}");
  const TString sDebugTruClustSum10("E^{sum/truth}_{clust}(10 mm) = E^{lead/truth}_{clust} + #Sigma_{#Deltar < 10 mm} e^{truth}_{clust}");
  const TString sDebugTruClustSum100("E^{sum/truth}_{clust}(1 cm) = E^{lead/truth}_{clust} + #Sigma_{#Deltar < 1 cm} e^{truth}_{clust}");
  const TString sDebugTruClustSum1000("E^{sum/truth}_{clust}(1 m) = E^{lead/truth}_{clust} + #Sigma_{#Deltar < 1 m} e^{truth}_{clust}");

  // debug difference axis titles
  const TString sDebugClustDiff5("#DeltaE^{sum}_{clust}(5 mm) = [E^{sum}_{clust}(5 mm) - E_{par}] / E^{sum}_{clust}(5 mm)");
  const TString sDebugClustDiff10("#DeltaE#^{sum}_{clust}(10 mm) = [E^{sum}_{clust}(10 mm) - E_{par}] / E^{sum}_{clust}(10 mm)");
  const TString sDebugClustDiff100("#DeltaE^{sum}_{clust}(1 cm) = [E^{sum}_{clust}(1 cm) - E_{par}] / E^{sum}_{clust}(1 cm)");
  const TString sDebugClustDiff1000("#DeltaE^{sum}_{clust}(1 m) = [E^{sum}_{clust}(1 m) - E_{par}] / E^{sum}_{clust}(1 m)");
  const TString sDebugTruClustDiff5("#DeltaE^{sum/truth}_{clust}(5 mm) = [E^{sum/truth}_{clust}(5 mm) - E_{par}] / E^{sum/truth}_{clust}(5 mm)");
  const TString sDebugTruClustDiff10("#DeltaE^{sum/truth}_{clust}(10 mm) = [E^{sum/truth}_{clust}(10 mm) - E_{par}] / E^{sum/truth}_{clust}(10 mm)");
  const TString sDebugTruClustDiff100("#DeltaE^{sum/truth}_{clust}(1 cm) = [E^{sum/truth}_{clust}(1 cm) - E_{par}] / E^{sum/truth}_{clust}(1 cm)");
  const TString sDebugTruClustDiff1000("#DeltaE^{sum/truth}_{clust}(1 m) = [E^{sum/truth}_{clust}(1 m) - E_{par}] / E^{sum/truth}_{clust}(1 m)");

  // set hcal debug axis titles
  hHCalDebugClustSum5        -> GetXaxis() -> SetTitle(sDebugClustSum5.Data());
  hHCalDebugClustSum5        -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalDebugClustSum10       -> GetXaxis() -> SetTitle(sDebugClustSum10.Data());
  hHCalDebugClustSum10       -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalDebugClustSum100      -> GetXaxis() -> SetTitle(sDebugClustSum100.Data());
  hHCalDebugClustSum100      -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalDebugClustSum1000     -> GetXaxis() -> SetTitle(sDebugClustSum1000.Data());
  hHCalDebugClustSum1000     -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalDebugTruClustSum5     -> GetXaxis() -> SetTitle(sDebugTruClustSum5.Data());
  hHCalDebugTruClustSum5     -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalDebugTruClustSum10    -> GetXaxis() -> SetTitle(sDebugTruClustSum10.Data());
  hHCalDebugTruClustSum10    -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalDebugTruClustSum100   -> GetXaxis() -> SetTitle(sDebugTruClustSum100.Data());
  hHCalDebugTruClustSum100   -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalDebugTruClustSum1000  -> GetXaxis() -> SetTitle(sDebugTruClustSum1000.Data());
  hHCalDebugTruClustSum1000  -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalDebugClustDiff5       -> GetXaxis() -> SetTitle(sDebugClustDiff5.Data());
  hHCalDebugClustDiff5       -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalDebugClustDiff10      -> GetXaxis() -> SetTitle(sDebugClustDiff10.Data());
  hHCalDebugClustDiff10      -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalDebugClustDiff100     -> GetXaxis() -> SetTitle(sDebugClustDiff100.Data());
  hHCalDebugClustDiff100     -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalDebugClustDiff1000    -> GetXaxis() -> SetTitle(sDebugClustDiff1000.Data());
  hHCalDebugClustDiff1000    -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalDebugTruClustDiff5    -> GetXaxis() -> SetTitle(sDebugTruClustDiff5.Data());
  hHCalDebugTruClustDiff5    -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalDebugTruClustDiff10   -> GetXaxis() -> SetTitle(sDebugTruClustDiff10.Data());
  hHCalDebugTruClustDiff10   -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalDebugTruClustDiff100  -> GetXaxis() -> SetTitle(sDebugTruClustDiff100.Data());
  hHCalDebugTruClustDiff100  -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalDebugTruClustDiff1000 -> GetXaxis() -> SetTitle(sDebugTruClustDiff1000.Data());
  hHCalDebugTruClustDiff1000 -> GetYaxis() -> SetTitle(sCount.Data());
  return;

}  // end 'FinishWithGlobalRootLock()'

// end ------------------------------------------------------------------------
