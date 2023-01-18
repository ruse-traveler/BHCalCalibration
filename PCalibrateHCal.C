// 'PCalibrateHCal.C'
// Derek Anderson
// 01.18.2023
//
// A macro to read in PODIO
// collections relating to
// the Barrel HCal from
// EICrecon and produce
// several histograms.

R__LOAD_LIBRARY(podioRootIO)
R__LOAD_LIBRARY(podioDict)
R__LOAD_LIBRARY(edm4eic)

// c includes
#include <cmath>
#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>
// root includes
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TError.h>
#include <TString.h>
#include <TDirectory.h>
// podio includes
#include <podio/EventStore.h>
#include <podio/ROOTReader.h>
#include <podio/CollectionBase.h>
// misc includes
#include <boost/math/special_functions/sign.hpp>

using namespace std;

// global constants
static const size_t NRange       = 2;
static const size_t NComp        = 3;
static const size_t NType        = 7;
static const float  MParMinDef   = 0.135;
static const float  MParMaxDef   = 0.145;
static const float  EParMinDef   = 4.9;
static const float  EParMaxDef   = 5.1;
static const bool   IsInBatchDef = false;
static const string SOutputDef   = "test_out.root";
static const string SInputDef    = "../forPodioReaderTest_fromEicRecon.e5th70n500pip.d18m1y2023.podio.root";



void PCalibrateHCal(const string sOutput = SOutputDef, const string sInput = SInputDef, const float mParMin = MParMinDef, const float mParMax = MParMaxDef, const float eParMin = EParMinDef, const float eParMax = EParMaxDef, const bool isInBatchMode = IsInBatchDef) {

  // announce start
  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning PodioReader-based calibration script..." << endl;

  // open input and store events
  podio::ROOTReader reader;
  podio::EventStore store;
  reader.openFile(sInput);
  store.setReader(&reader);
  cout << "    Grabbed input events." << endl;

  // create output file
  TFile *output = new TFile(sOutput.data(), "recreate");
  if (!output) {
    cerr << "PANIC: couldn't open output file!\n" << endl;
    assert(output);
  } else {
    cout << "    Opened output file." << endl;
  }

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
  TH1D *hParChrg                   = new TH1D("hParChrg",     "Gen. Particles", nChrgBin,  rChrgBin[0], rChrgBin[1]);
  TH1D *hParMass                   = new TH1D("hParMass",     "Gen. Particles", nMassBin,  rMassBin[0], rMassBin[1]);
  TH1D *hParPhi                    = new TH1D("hParPhi",      "Gen. Particles", nPhiBin,   rPhiBin[0],  rPhiBin[1]);
  TH1D *hParEta                    = new TH1D("hParEta",      "Gen. Particles", nEtaBin,   rEtaBin[0],  rEtaBin[1]);
  TH1D *hParEne                    = new TH1D("hParEne",      "Gen. Particles", nEneBin,   rEneBin[0],  rEneBin[1]);
  TH1D *hParMom                    = new TH1D("hParMom",      "Gen. Particles", nEneBin,   rEneBin[0],  rEneBin[1]);
  TH1D *hParMomX                   = new TH1D("hParMomX",     "Gen. Particles", nMomBin,   rMomBin[0],  rMomBin[1]);
  TH1D *hParMomY                   = new TH1D("hParMomY",     "Gen. Particles", nMomBin,   rMomBin[0],  rMomBin[1]);
  TH1D *hParMomZ                   = new TH1D("hParMomZ",     "Gen. Particles", nMomBin,   rMomBin[0],  rMomBin[1]);
  TH2D *hParEtaVsPhi               = new TH2D("hParEtaVsPhi", "Gen. Particles", nPhiBin,   rPhiBin[0],  rPhiBin[1],   nEtaBin,  rEtaBin[0],   rEtaBin[1]);
  // reco. hcal hit histograms
  TH1D *hHCalRecHitPhi             = new TH1D("hHCalRecHitPhi",      "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  TH1D *hHCalRecHitEta             = new TH1D("hHCalRecHitEta",      "Barrel HCal", nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  TH1D *hHCalRecHitEne             = new TH1D("hHCalRecHitEne",      "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  TH1D *hHCalRecHitPosZ            = new TH1D("hHCalRecHitPosZ",     "Barrel HCal", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  TH1D *hHCalRecHitParDiff         = new TH1D("hHCalRecHitParDiff",  "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  TH2D *hHCalRecHitPosYvsX         = new TH2D("hHCalRecHitPosYvsX",  "Barrel HCal", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  TH2D *hHCalRecHitEtaVsPhi        = new TH2D("hHCalRecHitEtaVsPhi", "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  TH2D *hHCalRecHitVsParEne        = new TH2D("hHCalRecHitVsParEne", "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // hcal cluster hit histograms
  TH1D *hHCalClustHitPhi           = new TH1D("hHCalClustHitPhi",      "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  TH1D *hHCalClustHitEta           = new TH1D("hHCalClustHitEta",      "Barrel HCal", nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  TH1D *hHCalClustHitEne           = new TH1D("hHCalClustHitEne",      "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  TH1D *hHCalClustHitPosZ          = new TH1D("hHCalClustHitPosZ",     "Barrel HCal", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  TH1D *hHCalClustHitParDiff       = new TH1D("hHCalClustHitParDiff",  "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  TH2D *hHCalClustHitPosYvsX       = new TH2D("hHCalClustHitPosYvsX",  "Barrel HCal", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  TH2D *hHCalClustHitEtaVsPhi      = new TH2D("hHCalClustHitEtaVsPhi", "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  TH2D *hHCalClustHitVsParEne      = new TH2D("hHCalClustHitVsParEne", "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // reco. hcal cluster histograms
  TH1D *hHCalClustPhi              = new TH1D("hHCalClustPhi",      "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  TH1D *hHCalClustEta              = new TH1D("hHCalClustEta",      "Barrel HCal", nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  TH1D *hHCalClustEne              = new TH1D("hHCalClustEne",      "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  TH1D *hHCalClustPosZ             = new TH1D("hHCalClustPosZ",     "Barrel HCal", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  TH1I *hHCalClustNumHit           = new TH1I("hHCalClustNumHit",   "Barrel HCal", nNumBin,   rNumBin[0],   rNumBin[1]);
  TH1D *hHCalClustParDiff          = new TH1D("hHCalClustParDiff",  "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  TH2D *hHCalClustPosYvsX          = new TH2D("hHCalClustPosYvsX",  "Barrel HCal", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  TH2D *hHCalClustEtaVsPhi         = new TH2D("hHCalClustEtaVsPhi", "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  TH2D *hHCalClustVsParEne         = new TH2D("hHCalClustVsParEne", "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // hcal cluster hit histograms
  TH1D *hHCalTruClustHitPhi        = new TH1D("hHCalTruClustHitPhi",      "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  TH1D *hHCalTruClustHitEta        = new TH1D("hHCalTruClustHitEta",      "Barrel HCal", nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  TH1D *hHCalTruClustHitEne        = new TH1D("hHCalTruClustHitEne",      "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  TH1D *hHCalTruClustHitPosZ       = new TH1D("hHCalTruClustHitPosZ",     "Barrel HCal", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  TH1D *hHCalTruClustHitParDiff    = new TH1D("hHCalTruClustHitParDiff",  "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  TH2D *hHCalTruClustHitPosYvsX    = new TH2D("hHCalTruClustHitPosYvsX",  "Barrel HCal", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  TH2D *hHCalTruClustHitEtaVsPhi   = new TH2D("hHCalTruClustHitEtaVsPhi", "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  TH2D *hHCalTruClustHitVsParEne   = new TH2D("hHCalTruClustHitVsParEne", "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // truth hcal cluster histograms
  TH1D *hHCalTruClustPhi           = new TH1D("hHCalTruClustPhi",      "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1]);
  TH1D *hHCalTruClustEta           = new TH1D("hHCalTruClustEta",      "Barrel HCal", nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  TH1D *hHCalTruClustEne           = new TH1D("hHCalTruClustEne",      "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  TH1D *hHCalTruClustPosZ          = new TH1D("hHCalTruClustPosZ",     "Barrel HCal", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  TH1I *hHCalTruClustNumHit        = new TH1I("hHCalTruClustNumHit",   "Barrel HCal", nNumBin,   rNumBin[0],   rNumBin[1]);
  TH1D *hHCalTruClustParDiff       = new TH1D("hHCalTruClustParDiff",  "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  TH2D *hHCalTruClustPosYvsX       = new TH2D("hHCalTruClustPosYvsX",  "Barrel HCal", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  TH2D *hHCalTruClustEtaVsPhi      = new TH2D("hHCalTruClustEtaVsPhi", "Barrel HCal", nPhiBin,   rPhiBin[0],   rPhiBin[1],   nEtaBin,   rEtaBin[0],   rEtaBin[1]);
  TH2D *hHCalTruClustVsParEne      = new TH2D("hHCalTruClustVsParEne", "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // event-wise hcal histograms
  TH1I *hEvtHCalNumPar             = new TH1I("hEvtHCalNumPar",             "Barrel HCal", nNumBin,   rNumBin[0],   rNumBin[1]);
  TH1I *hEvtHCalNumHit             = new TH1I("hEvtHCalNumHit",             "Barrel HCal", nNumBin,   rNumBin[0],   rNumBin[1]);
  TH1I *hEvtHCalNumClust           = new TH1I("hEvtHCalNumClust",           "Barrel HCal", nNumBin,   rNumBin[0],   rNumBin[1]);
  TH1I *hEvtHCalNumTruClust        = new TH1I("hEvtHCalNumTruClust",        "Barrel HCal", nNumBin,   rNumBin[0],   rNumBin[1]);
  TH1D *hEvtHCalSumHitEne          = new TH1D("hEvtHCalSumHitEne",          "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  TH1D *hEvtHCalSumClustEne        = new TH1D("hEvtHCalSumClustEne",        "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  TH1D *hEvtHCalSumTruClustEne     = new TH1D("hEvtHCalSumTruClustEne",     "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  TH1D *hEvtHCalLeadClustEne       = new TH1D("hEvtHCalLeadClustEne",       "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  TH1D *hEvtHCalLeadTruClustEne    = new TH1D("hEvtHCalLeadTruClustEne",    "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1]);
  TH1D *hEvtHCalSumHitDiff         = new TH1D("hEvtHCalSumHitDiff",         "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  TH1D *hEvtHCalSumClustDiff       = new TH1D("hEvtHCalSumClustDiff",       "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  TH1D *hEvtHCalSumTruClustDiff    = new TH1D("hEvtHCalSumTruClustDiff",    "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  TH1D *hEvtHCalLeadClustDiff      = new TH1D("hEvtHCalLeadClustDiff",      "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  TH1D *hEvtHCalLeadTruClustDiff   = new TH1D("hEvtHCalLeadTruClustDiff",   "Barrel HCal", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  TH2I *hEvtHCalNumClustVsHit      = new TH2I("hEvtHCalNumClustVsHit",      "Barrel HCal", nNumBin,   rNumBin[0],   rNumBin[1],   nNumBin,   rNumBin[0],   rNumBin[1]);
  TH2I *hEvtHCalNumTruClustVsClust = new TH2I("hEvtHCalNumTruClustVsClust", "Barrel HCal", nNumBin,   rNumBin[0],   rNumBin[1],   nNumBin,   rNumBin[0],   rNumBin[1]);
  TH2D *hEvtHCalSumHitVsPar        = new TH2D("hEvtHCalSumHitVsPar",        "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  TH2D *hEvtHCalSumClustVsPar      = new TH2D("hEvtHCalSumClustVsPar",      "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  TH2D *hEvtHCalSumTruClustVsPar   = new TH2D("hEvtHCalSumTruClustVsPar",   "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  TH2D *hEvtHCalLeadClustVsPar     = new TH2D("hEvtHCalLeadClustVsPar",     "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  TH2D *hEvtHCalLeadTruClustVsPar  = new TH2D("hEvtHCalLeadTruClustVsPar",  "Barrel HCal", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
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

  // event loop
  const auto nEvts = reader.getEntries();
  cout << "    Beginning event loop: " << nEvts << " events to process." << endl;

  for (unsigned int iEvt = 0; iEvt < nEvts; iEvt++) {

    // announce progress
    const auto iProg = iEvt + 1;
    if (isInBatchMode) {
      cout << "      Processing event " << iProg << "/" << nEvts << "..." << endl;
    } else {
      cout << "      Processing event " << iProg << "/" << nEvts << "...\r" << flush;
      if (iProg == nEvts) cout << endl;
    }

    // grab relevant collections
    const auto &genParticles       = store.get<edm4eic::ReconstructedParticleCollection>("GeneratedParticles");
    const auto &bhcalRecHits       = store.get<edm4eic::CalorimeterHitCollection>("HcalBarrelRecHits");
    const auto &bhcalClusters      = store.get<edm4eic::ClusterCollection>("HcalBarrelClusters");
    const auto &bhcalTruthClusters = store.get<edm4eic::ClusterCollection>("HcalBarrelTruthClusters");

    // for hit and cluster sums
    float eHCalHitSum(0.);
    float eHCalClustSum(0.);
    float eTruHCalClustSum(0.);

    // sum hcal hit energy
    for (auto bhCalHit : bhcalRecHits) {
      eHCalHitSum += bhCalHit.getEnergy();
    }  // end 1st hcal hit loop

    // if hit sum is 0, skip event
    const bool isHCalHitSumNonzero = (eHCalHitSum > 0.);
    if (!isHCalHitSumNonzero) {
      store.clear();
      reader.endOfEvent();
    }

    // clear store and prepare for next event
    store.clear();
    reader.endOfEvent();

  }  // end event loop
  cout << "    Finished event loop!" << endl;

  // create output directories
  const string directNames[NType] = {"GenParticles", "RecoHits", "RecoClustHits", "RecoClusters", "TruthClustHits", "TruthClusters", "EventInfo"};

  TDirectory *outDir[NType];
  for (unsigned int iDir = 0; iDir < NType; iDir++) {
    output       -> cd();
    outDir[iDir] = (TDirectory*) output -> mkdir(directNames[iDir].data());
  }
  cout << "    Made output directories." << endl;

  // save histograms
  outDir[0]                  -> cd();
  hParChrg                   -> Write();
  hParMass                   -> Write();
  hParPhi                    -> Write();
  hParEta                    -> Write();
  hParEne                    -> Write();
  hParMom                    -> Write();
  hParMomX                   -> Write();
  hParMomY                   -> Write();
  hParMomZ                   -> Write();
  hParEtaVsPhi               -> Write();
  outDir[1]                  -> cd();
  hHCalRecHitPhi             -> Write();
  hHCalRecHitEta             -> Write();
  hHCalRecHitEne             -> Write();
  hHCalRecHitPosZ            -> Write();
  hHCalRecHitParDiff         -> Write();
  hHCalRecHitPosYvsX         -> Write();
  hHCalRecHitEtaVsPhi        -> Write();
  hHCalRecHitVsParEne        -> Write();
  outDir[2]                 -> cd();
  hHCalClustHitPhi           -> Write();
  hHCalClustHitEta           -> Write();
  hHCalClustHitEne           -> Write();
  hHCalClustHitPosZ          -> Write();
  hHCalClustHitParDiff       -> Write();
  hHCalClustHitPosYvsX       -> Write();
  hHCalClustHitEtaVsPhi      -> Write();
  hHCalClustHitVsParEne      -> Write();
  outDir[3]                  -> cd();
  hHCalClustPhi              -> Write();
  hHCalClustEta              -> Write();
  hHCalClustEne              -> Write();
  hHCalClustPosZ             -> Write();
  hHCalClustNumHit           -> Write();
  hHCalClustParDiff          -> Write();
  hHCalClustPosYvsX          -> Write();
  hHCalClustEtaVsPhi         -> Write();
  hHCalClustVsParEne         -> Write();
  outDir[4]                  -> cd();
  hHCalTruClustHitPhi        -> Write();
  hHCalTruClustHitEta        -> Write();
  hHCalTruClustHitEne        -> Write();
  hHCalTruClustHitPosZ       -> Write();
  hHCalTruClustHitParDiff    -> Write();
  hHCalTruClustHitPosYvsX    -> Write();
  hHCalTruClustHitEtaVsPhi   -> Write();
  hHCalTruClustHitVsParEne   -> Write();
  outDir[5]                  -> cd();
  hHCalTruClustPhi           -> Write();
  hHCalTruClustEta           -> Write();
  hHCalTruClustEne           -> Write();
  hHCalTruClustPosZ          -> Write();
  hHCalTruClustNumHit        -> Write();
  hHCalTruClustParDiff       -> Write();
  hHCalTruClustPosYvsX       -> Write();
  hHCalTruClustEtaVsPhi      -> Write();
  hHCalTruClustVsParEne      -> Write();
  outDir[6]                  -> cd();
  hEvtHCalNumPar             -> Write();
  hEvtHCalNumHit             -> Write();
  hEvtHCalNumClust           -> Write();
  hEvtHCalNumTruClust        -> Write();
  hEvtHCalSumHitEne          -> Write();
  hEvtHCalSumClustEne        -> Write();
  hEvtHCalSumTruClustEne     -> Write();
  hEvtHCalLeadClustEne       -> Write();
  hEvtHCalLeadTruClustEne    -> Write();
  hEvtHCalSumHitDiff         -> Write();
  hEvtHCalSumClustDiff       -> Write();
  hEvtHCalSumTruClustDiff    -> Write();
  hEvtHCalLeadClustDiff      -> Write();
  hEvtHCalLeadTruClustDiff   -> Write();
  hEvtHCalNumClustVsHit      -> Write();
  hEvtHCalNumTruClustVsClust -> Write();
  hEvtHCalSumHitVsPar        -> Write();
  hEvtHCalSumClustVsPar      -> Write();
  hEvtHCalLeadClustVsPar     -> Write();
  hEvtHCalLeadTruClustVsPar  -> Write();
  cout << "    Saved histograms." << endl;

  // close files and exit
  output -> cd();
  output -> Close();
  reader.closeFile();
  cout << "  Finished calibration script!\n" << endl;

}

// end ------------------------------------------------------------------------
