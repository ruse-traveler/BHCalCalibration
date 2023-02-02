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
  hParChrg                   = new TH1D("hParChrg",     "Gen. Particles", nChrgBin,  rChrgBin[0], rChrgBin[1]);
  hParMass                   = new TH1D("hParMass",     "Gen. Particles", nMassBin,  rMassBin[0], rMassBin[1]);
  hParPhi                    = new TH1D("hParPhi",      "Gen. Particles", nPhiBin,   rPhiBin[0],  rPhiBin[1]);
  hParEta                    = new TH1D("hParEta",      "Gen. Particles", nEtaBin,   rEtaBin[0],  rEtaBin[1]);
  hParEne                    = new TH1D("hParEne",      "Gen. Particles", nEneBin,   rEneBin[0],  rEneBin[1]);
  hParMom                    = new TH1D("hParMom",      "Gen. Particles", nEneBin,   rEneBin[0],  rEneBin[1]);
  hParMomX                   = new TH1D("hParMomX",     "Gen. Particles", nMomBin,   rMomBin[0],  rMomBin[1]);
  hParMomY                   = new TH1D("hParMomY",     "Gen. Particles", nMomBin,   rMomBin[0],  rMomBin[1]);
  hParMomZ                   = new TH1D("hParMomZ",     "Gen. Particles", nMomBin,   rMomBin[0],  rMomBin[1]);
  hParEtaVsPhi               = new TH2D("hParEtaVsPhi", "Gen. Particles", nPhiBin,   rPhiBin[0],  rPhiBin[1],   nEtaBin,  rEtaBin[0],   rEtaBin[1]);
  // reco. hcal hit histograms
  hHCalRecHitPhi             = new TH1D("hHCalRecHitPhi",      "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  hHCalRecHitEta             = new TH1D("hHCalRecHitEta",      "Barrel HCal", nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalRecHitEne             = new TH1D("hHCalRecHitEne",      "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalRecHitPosZ            = new TH1D("hHCalRecHitPosZ",     "Barrel HCal", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hHCalRecHitParDiff         = new TH1D("hHCalRecHitParDiff",  "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalRecHitPosYvsX         = new TH2D("hHCalRecHitPosYvsX",  "Barrel HCal", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hHCalRecHitEtaVsPhi        = new TH2D("hHCalRecHitEtaVsPhi", "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalRecHitVsParEne        = new TH2D("hHCalRecHitVsParEne", "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // hcal cluster hit histograms
  hHCalClustHitPhi           = new TH1D("hHCalClustHitPhi",      "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  hHCalClustHitEta           = new TH1D("hHCalClustHitEta",      "Barrel HCal", nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalClustHitEne           = new TH1D("hHCalClustHitEne",      "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalClustHitPosZ          = new TH1D("hHCalClustHitPosZ",     "Barrel HCal", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hHCalClustHitParDiff       = new TH1D("hHCalClustHitParDiff",  "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalClustHitPosYvsX       = new TH2D("hHCalClustHitPosYvsX",  "Barrel HCal", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hHCalClustHitEtaVsPhi      = new TH2D("hHCalClustHitEtaVsPhi", "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalClustHitVsParEne      = new TH2D("hHCalClustHitVsParEne", "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // reco. hcal cluster histograms
  hHCalClustPhi              = new TH1D("hHCalClustPhi",      "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  hHCalClustEta              = new TH1D("hHCalClustEta",      "Barrel HCal", nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalClustEne              = new TH1D("hHCalClustEne",      "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalClustPosZ             = new TH1D("hHCalClustPosZ",     "Barrel HCal", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hHCalClustNumHit           = new TH1I("hHCalClustNumHit",   "Barrel HCal", nNumBin,   rNumBin[0],   rNumBin[1]);
  hHCalClustParDiff          = new TH1D("hHCalClustParDiff",  "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalClustPosYvsX          = new TH2D("hHCalClustPosYvsX",  "Barrel HCal", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hHCalClustEtaVsPhi         = new TH2D("hHCalClustEtaVsPhi", "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalClustVsParEne         = new TH2D("hHCalClustVsParEne", "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // hcal cluster hit histograms
  hHCalTruClustHitPhi        = new TH1D("hHCalTruClustHitPhi",      "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  hHCalTruClustHitEta        = new TH1D("hHCalTruClustHitEta",      "Barrel HCal", nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalTruClustHitEne        = new TH1D("hHCalTruClustHitEne",      "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalTruClustHitPosZ       = new TH1D("hHCalTruClustHitPosZ",     "Barrel HCal", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hHCalTruClustHitParDiff    = new TH1D("hHCalTruClustHitParDiff",  "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalTruClustHitPosYvsX    = new TH2D("hHCalTruClustHitPosYvsX",  "Barrel HCal", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hHCalTruClustHitEtaVsPhi   = new TH2D("hHCalTruClustHitEtaVsPhi", "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalTruClustHitVsParEne   = new TH2D("hHCalTruClustHitVsParEne", "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // truth hcal cluster histograms
  hHCalTruClustPhi           = new TH1D("hHCalTruClustPhi",      "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  hHCalTruClustEta           = new TH1D("hHCalTruClustEta",      "Barrel HCal", nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalTruClustEne           = new TH1D("hHCalTruClustEne",      "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalTruClustPosZ          = new TH1D("hHCalTruClustPosZ",     "Barrel HCal", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hHCalTruClustNumHit        = new TH1I("hHCalTruClustNumHit",   "Barrel HCal", nNumBin,   rNumBin[0],   rNumBin[1]);
  hHCalTruClustParDiff       = new TH1D("hHCalTruClustParDiff",  "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalTruClustPosYvsX       = new TH2D("hHCalTruClustPosYvsX",  "Barrel HCal", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hHCalTruClustEtaVsPhi      = new TH2D("hHCalTruClustEtaVsPhi", "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  hHCalTruClustVsParEne      = new TH2D("hHCalTruClustVsParEne", "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // event-wise hcal histograms
  hEvtHCalNumPar             = new TH1I("hEvtHCalNumPar",             "Barrel HCal", nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtHCalNumHit             = new TH1I("hEvtHCalNumHit",             "Barrel HCal", nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtHCalNumClust           = new TH1I("hEvtHCalNumClust",           "Barrel HCal", nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtHCalNumTruClust        = new TH1I("hEvtHCalNumTruClust",        "Barrel HCal", nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtHCalSumHitEne          = new TH1D("hEvtHCalSumHitEne",          "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalSumClustEne        = new TH1D("hEvtHCalSumClustEne",        "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalSumTruClustEne     = new TH1D("hEvtHCalSumTruClustEne",     "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalLeadClustEne       = new TH1D("hEvtHCalLeadClustEne",       "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalLeadTruClustEne    = new TH1D("hEvtHCalLeadTruClustEne",    "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalSumHitDiff         = new TH1D("hEvtHCalSumHitDiff",         "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtHCalSumClustDiff       = new TH1D("hEvtHCalSumClustDiff",       "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtHCalSumTruClustDiff    = new TH1D("hEvtHCalSumTruClustDiff",    "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtHCalLeadClustDiff      = new TH1D("hEvtHCalLeadClustDiff",      "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtHCalLeadTruClustDiff   = new TH1D("hEvtHCalLeadTruClustDiff",   "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtHCalNumClustVsHit      = new TH2I("hEvtHCalNumClustVsHit",      "Barrel HCal", nNumBin,   rNumBin[0],   rNumBin[1],   nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtHCalNumTruClustVsClust = new TH2I("hEvtHCalNumTruClustVsClust", "Barrel HCal", nNumBin,   rNumBin[0],   rNumBin[1],   nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtHCalSumHitVsPar        = new TH2D("hEvtHCalSumHitVsPar",        "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalSumClustVsPar      = new TH2D("hEvtHCalSumClustVsPar",      "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalSumTruClustVsPar   = new TH2D("hEvtHCalSumTruClustVsPar",   "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalLeadClustVsPar     = new TH2D("hEvtHCalLeadClustVsPar",     "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtHCalLeadTruClustVsPar  = new TH2D("hEvtHCalLeadTruClustVsPar",  "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
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
  hHCalClustHitPhi           -> Sumw2();
  hHCalClustHitEta           -> Sumw2();
  hHCalClustHitEne           -> Sumw2();
  hHCalClustHitPosZ          -> Sumw2();
  hHCalClustHitParDiff       -> Sumw2();
  hHCalClustHitPosYvsX       -> Sumw2();
  hHCalClustHitEtaVsPhi      -> Sumw2();
  hHCalClustHitVsParEne      -> Sumw2();
  hHCalClustPhi              -> Sumw2();
  hHCalClustEta              -> Sumw2();
  hHCalClustEne              -> Sumw2();
  hHCalClustPosZ             -> Sumw2();
  hHCalClustNumHit           -> Sumw2();
  hHCalClustParDiff          -> Sumw2();
  hHCalClustPosYvsX          -> Sumw2();
  hHCalClustEtaVsPhi         -> Sumw2();
  hHCalClustVsParEne         -> Sumw2();
  hHCalTruClustHitPhi        -> Sumw2();
  hHCalTruClustHitEta        -> Sumw2();
  hHCalTruClustHitEne        -> Sumw2();
  hHCalTruClustHitPosZ       -> Sumw2();
  hHCalTruClustHitParDiff    -> Sumw2();
  hHCalTruClustHitPosYvsX    -> Sumw2();
  hHCalTruClustHitEtaVsPhi   -> Sumw2();
  hHCalTruClustHitVsParEne   -> Sumw2();
  hHCalTruClustPhi           -> Sumw2();
  hHCalTruClustEta           -> Sumw2();
  hHCalTruClustEne           -> Sumw2();
  hHCalTruClustPosZ          -> Sumw2();
  hHCalTruClustNumHit        -> Sumw2();
  hHCalTruClustParDiff       -> Sumw2();
  hHCalTruClustPosYvsX       -> Sumw2();
  hHCalTruClustEtaVsPhi      -> Sumw2();
  hHCalTruClustVsParEne      -> Sumw2();
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
  double eHCalClustSum(0.);
  double eTruHCalClustSum(0.);

  // sum hcal hit energy
  for (auto bhCalHit : bhcalRecHits()) {
    eHCalHitSum += bhCalHit -> getEnergy();
  }  // end 1st hcal hit loop

  // if hit sum is 0, skip event
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
    const bool isRightMomentum = ((ePar >= EParMin) && (ePar <= EParMax));
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

  // get protoclusters
  auto bhCalProtoClusters = event -> Get<edm4eic::ProtoCluster>("HcalBarrelIslandProtoClusters");

  // reco. hcal cluster loop
  unsigned long iHCalClust(0);
  unsigned long nHCalProto(0);
  unsigned long nHCalClust(0);
  for (auto bhCalClust : bhcalClusters()) {
    
    // loop over protoclusters
    unsigned long iHCalProto(0);
    unsigned long nProtoHits(0);
    for (auto bhCalProto : bhCalProtoClusters) {

      // check if proto index is same as cluster index
      // FIXME: this might not be the correct way to match reco and proto clusters
      const bool isSameIndex = (iHCalProto == iHCalClust);
      if (!isSameIndex) continue;

      // loop over hits
      nProtoHits = bhCalProto -> hits_size();
      for (uint32_t iProtoHit = 0; iProtoHit < nProtoHits; iProtoHit++) {

        // get hit
        const auto bhCalProtoHit = bhCalProto -> getHits(iProtoHit);

        // grab hit properties
        const auto rHCalProtoHitX   = bhCalProtoHit.getPosition().x;
        const auto rHCalProtoHitY   = bhCalProtoHit.getPosition().y;
        const auto rHCalProtoHitZ   = bhCalProtoHit.getPosition().z;
        const auto eHCalProtoHit    = bhCalProtoHit.getEnergy();
        const auto rHCalProtoHitS   = std::sqrt((rHCalProtoHitX * rHCalProtoHitX) + (rHCalProtoHitY * rHCalProtoHitY));
        const auto rHCalProtoHitR   = std::sqrt((rHCalProtoHitS * rHCalProtoHitS) + (rHCalProtoHitZ * rHCalProtoHitZ));
        const auto fHCalProtoHit    = boost::math::sign(rHCalProtoHitY) * acos(rHCalProtoHitX / rHCalProtoHitS);
        const auto tHCalProtoHit    = std::acos(rHCalProtoHitZ / rHCalProtoHitR);
        const auto hHCalProtoHit    = (-1.) * std::log(std::atan(tHCalProtoHit / 2.));
        const auto diffHCalProtoHit = (eHCalProtoHit - eMcPar) / eHCalProtoHit;

        // fill hit histograms and increment sums/counters
        hHCalClustHitPhi      -> Fill(fHCalProtoHit);
        hHCalClustHitEta      -> Fill(hHCalProtoHit);
        hHCalClustHitEne      -> Fill(eHCalProtoHit);
        hHCalClustHitPosZ     -> Fill(rHCalProtoHitZ);
        hHCalClustHitParDiff  -> Fill(diffHCalProtoHit);
        hHCalClustHitPosYvsX  -> Fill(rHCalProtoHitX, rHCalProtoHitY);
        hHCalClustHitEtaVsPhi -> Fill(fHCalProtoHit, hHCalProtoHit);
        hHCalClustHitVsParEne -> Fill(eMcPar, eHCalProtoHit);
      }
      ++nHCalProto;
    }  // end protocluster loop

    // grab cluster properties
    const auto rHCalClustX   = bhCalClust -> getPosition().x;
    const auto rHCalClustY   = bhCalClust -> getPosition().y;
    const auto rHCalClustZ   = bhCalClust -> getPosition().z;
    const auto eHCalClust    = bhCalClust -> getEnergy();
    const auto nHitHCalClust = bhCalClust -> getNhits();
    const auto fHCalClust    = bhCalClust -> getIntrinsicPhi();
    const auto tHCalClust    = bhCalClust -> getIntrinsicTheta();
    const auto hHCalClust    = (1.) * std::log(std::atan(tHCalClust / 2.));
    const auto diffHCalClust = (eHCalClust - eMcPar) / eHCalClust;

    // fill cluster histograms and increment counters
    hHCalClustPhi      -> Fill(fHCalClust);
    hHCalClustEta      -> Fill(hHCalClust);
    hHCalClustEne      -> Fill(eHCalClust);
    hHCalClustPosZ     -> Fill(rHCalClustZ);
    hHCalClustNumHit   -> Fill(nProtoHits);
    hHCalClustParDiff  -> Fill(diffHCalClust);
    hHCalClustPosYvsX  -> Fill(rHCalClustX, rHCalClustY);
    hHCalClustEtaVsPhi -> Fill(fHCalClust, hHCalClust);
    hHCalClustVsParEne -> Fill(eMcPar, eHCalClust);
    eHCalClustSum += eHCalClust;
    ++nHCalClust;
    ++iHCalClust;

    // select leading cluster
    const bool isBiggerEne = (eHCalClust > eLeadHCalClust);
    if (isBiggerEne) {
      iLeadHCalClust    = iHCalClust;
      eLeadHCalClust    = eHCalClust;
      diffLeadHCalClust = diffHCalClust;
    }
  }  // end reco. hcal cluster loop

  // get truth protoclusters
  auto bhCalTruProtoClusters = event -> Get<edm4eic::ProtoCluster>("HcalBarrelTruthProtoClusters");

  // true hcal cluster loop
  unsigned long iTruHCalClust(0);
  unsigned long nTruHCalProto(0);
  unsigned long nTruHCalClust(0);
  for (auto truthHCalClust : bhcalTruthClusters()) {
    
    // loop over protoclusters
    unsigned long iTruHCalProto(0);
    unsigned long nTruProtoHits(0);
    for (auto bhCalTruProto : bhCalTruProtoClusters) {

      // check if truth proto index is same as truth cluster index
      // FIXME: this might not be the correct way to match reco and proto clusters
      const bool isSameIndex = (iTruHCalProto == iTruHCalClust);
      if (!isSameIndex) continue;

      // loop over hits
      nTruProtoHits = bhCalTruProto -> hits_size();
      for (uint32_t iTruProtoHit = 0; iTruProtoHit < nTruProtoHits; iTruProtoHit++) {

        // get hit
        const auto bhCalTruProtoHit = bhCalTruProto -> getHits(iTruProtoHit);

        // grab hit properties
        const auto rTruHCalProtoHitX   = bhCalTruProtoHit.getPosition().x;
        const auto rTruHCalProtoHitY   = bhCalTruProtoHit.getPosition().y;
        const auto rTruHCalProtoHitZ   = bhCalTruProtoHit.getPosition().z;
        const auto eTruHCalProtoHit    = bhCalTruProtoHit.getEnergy();
        const auto rTruHCalProtoHitS   = std::sqrt((rTruHCalProtoHitX * rTruHCalProtoHitX) + (rTruHCalProtoHitY * rTruHCalProtoHitY));
        const auto rTruHCalProtoHitR   = std::sqrt((rTruHCalProtoHitS * rTruHCalProtoHitS) + (rTruHCalProtoHitZ * rTruHCalProtoHitZ));
        const auto fTruHCalProtoHit    = boost::math::sign(rTruHCalProtoHitY) * acos(rTruHCalProtoHitX / rTruHCalProtoHitS);
        const auto tTruHCalProtoHit    = std::acos(rTruHCalProtoHitZ / rTruHCalProtoHitR);
        const auto hTruHCalProtoHit    = (-1.) * std::log(std::atan(tTruHCalProtoHit / 2.));
        const auto diffTruHCalProtoHit = (eTruHCalProtoHit - eMcPar) / eTruHCalProtoHit;

        // fill hit histograms and increment sums/counters
        hHCalTruClustHitPhi      -> Fill(fTruHCalProtoHit);
        hHCalTruClustHitEta      -> Fill(hTruHCalProtoHit);
        hHCalTruClustHitEne      -> Fill(eTruHCalProtoHit);
        hHCalTruClustHitPosZ     -> Fill(rTruHCalProtoHitZ);
        hHCalTruClustHitParDiff  -> Fill(diffTruHCalProtoHit);
        hHCalTruClustHitPosYvsX  -> Fill(rTruHCalProtoHitX, rTruHCalProtoHitY);
        hHCalTruClustHitEtaVsPhi -> Fill(fTruHCalProtoHit, hTruHCalProtoHit);
        hHCalTruClustHitVsParEne -> Fill(eMcPar, eTruHCalProtoHit);
      }
      ++nTruHCalProto;
    }  // end protocluster loop

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
  // set cluster hit hcal axis titles
  hHCalClustHitPhi          -> GetXaxis() -> SetTitle(sPhiHit.Data());
  hHCalClustHitPhi          -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustHitEta          -> GetXaxis() -> SetTitle(sEtaHit.Data());
  hHCalClustHitEta          -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustHitEne          -> GetXaxis() -> SetTitle(sEneHit.Data());
  hHCalClustHitEne          -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustHitPosZ         -> GetXaxis() -> SetTitle(sPosHitZ.Data());
  hHCalClustHitPosZ         -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustHitParDiff      -> GetXaxis() -> SetTitle(sEneHitDiff.Data());
  hHCalClustHitParDiff      -> GetYaxis() -> SetTitle(sCount.Data());
  hHCalClustHitPosYvsX      -> GetXaxis() -> SetTitle(sPosHitX.Data());
  hHCalClustHitPosYvsX      -> GetYaxis() -> SetTitle(sPosHitY.Data());
  hHCalClustHitPosYvsX      -> GetZaxis() -> SetTitle(sCount.Data());
  hHCalClustHitEtaVsPhi     -> GetXaxis() -> SetTitle(sPhiHit.Data());
  hHCalClustHitEtaVsPhi     -> GetYaxis() -> SetTitle(sEtaHit.Data());
  hHCalClustHitEtaVsPhi     -> GetZaxis() -> SetTitle(sCount.Data());
  hHCalClustHitVsParEne     -> GetXaxis() -> SetTitle(sEnePar.Data());
  hHCalClustHitVsParEne     -> GetYaxis() -> SetTitle(sEneHit.Data());
  hHCalClustHitVsParEne     -> GetZaxis() -> SetTitle(sCount.Data());
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
  return;

}  // end 'FinishWithGlobalRootLock()'

// end ------------------------------------------------------------------------
