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
  hParChrg                   = new TH1D("hParChrg",                   "", nChrgBin,  rChrgBin[0],  rChrgBin[1]);
  hParMass                   = new TH1D("hParMass",                   "", nMassBin,  rMassBin[0],  rMassBin[1]);
  hParEne                    = new TH1D("hParEne",                    "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hParMom                    = new TH1D("hParMom",                    "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hParMomX                   = new TH1D("hParMomX",                   "", nMomBin,   rMomBin[0],   rMomBin[1]);
  hParMomY                   = new TH1D("hParMomY",                   "", nMomBin,   rMomBin[0],   rMomBin[1]);
  hParMomZ                   = new TH1D("hParMomZ",                   "", nMomBin,   rMomBin[0],   rMomBin[1]);
  // reco. hcal hit histograms
  hHCalRecHitEne             = new TH1D("hHCalRecHitEne",             "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalRecHitPosZ            = new TH1D("hHCalRecHitPosZ",            "", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hHCalRecHitParDiff         = new TH1D("hHCalRecHitParDiff",         "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalRecHitPosYvsX         = new TH2D("hHCalRecHitPosYvsX",         "", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hHCalRecHitVsParEne        = new TH2D("hHCalRecHitVsParEne",        "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // reco. ecal hit histograms
  hECalRecHitEne             = new TH1D("hECalRecHitEne",             "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hECalRecHitPosZ            = new TH1D("hECalRecHitPosZ",            "", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hECalRecHitParDiff         = new TH1D("hECalRecHitParDiff",         "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hECalRecHitPosYvsX         = new TH2D("hECalRecHitPosYvsX",         "", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hECalRecHitVsParEne        = new TH2D("hECalRecHitVsParEne",        "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // reco. hcal cluster histograms
  hHCalClustEne              = new TH1D("hHCalClustEne",              "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalClustPosZ             = new TH1D("hHCalClustPosZ",             "", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hHCalClustNumHit           = new TH1I("hHCalClustNumHit",           "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hHCalClustParDiff          = new TH1D("hHCalClustParDiff",          "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalClustPosYvsX          = new TH2D("hHCalClustPosYvsX",          "", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hHCalClustVsParEne         = new TH2D("hHCalClustVsParEne",         "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // reco. ecal cluster histograms
  hECalClustEne              = new TH1D("hECalClustEne",              "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hECalClustPosZ             = new TH1D("hECalClustPosZ",             "", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hECalClustNumHit           = new TH1I("hECalClustNumHit",           "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hECalClustParDiff          = new TH1D("hECalClustParDiff",          "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hECalClustPosYvsX          = new TH2D("hECalClustPosYvsX",          "", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hECalClustVsParEne         = new TH2D("hECalClustVsParEne",         "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // reco. hcal cluster debug histograms
  hHCalDebugClustSum5        = new TH1D("hHCalDebugClustSum5",        "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugClustSum10       = new TH1D("hHCalDebugClustSum10",       "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugClustSum100      = new TH1D("hHCalDebugClustSum100",      "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugClustSum1000     = new TH1D("hHCalDebugClustSum1000",     "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugClustDiff5       = new TH1D("hHCalDebugClustDiff5",       "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalDebugClustDiff10      = new TH1D("hHCalDebugClustDiff10",      "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalDebugClustDiff100     = new TH1D("hHCalDebugClustDiff100",     "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalDebugClustDiff1000    = new TH1D("hHCalDebugClustDiff1000",    "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  // reco. ecal cluster debug histograms
  hECalDebugClustSum5        = new TH1D("hECalDebugClustSum5",        "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hECalDebugClustSum10       = new TH1D("hECalDebugClustSum10",       "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hECalDebugClustSum100      = new TH1D("hECalDebugClustSum100",      "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hECalDebugClustSum1000     = new TH1D("hECalDebugClustSum1000",     "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hECalDebugClustDiff5       = new TH1D("hECalDebugClustDiff5",       "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hECalDebugClustDiff10      = new TH1D("hECalDebugClustDiff10",      "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hECalDebugClustDiff100     = new TH1D("hECalDebugClustDiff100",     "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hECalDebugClustDiff1000    = new TH1D("hECalDebugClustDiff1000",    "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  // truth hcal cluster histograms
  hHCalTruClustEne           = new TH1D("hHCalTruClustEne",           "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalTruClustPosZ          = new TH1D("hHCalTruClustPosZ",          "", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hHCalTruClustNumHit        = new TH1I("hHCalTruClustNumHit",        "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hHCalTruClustParDiff       = new TH1D("hHCalTruClustParDiff",       "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalTruClustPosYvsX       = new TH2D("hHCalTruClustPosYvsX",       "", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hHCalTruClustVsParEne      = new TH2D("hHCalTruClustVsParEne",      "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // truth ecal cluster histograms
  hECalTruClustEne           = new TH1D("hECalTruClustEne",           "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hECalTruClustPosZ          = new TH1D("hECalTruClustPosZ",          "", nPosLoBin, rPosLoBin[0], rPosLoBin[1]);
  hECalTruClustNumHit        = new TH1I("hECalTruClustNumHit",        "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hECalTruClustParDiff       = new TH1D("hECalTruClustParDiff",       "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hECalTruClustPosYvsX       = new TH2D("hECalTruClustPosYvsX",       "", nPosTrBin, rPosTrBin[0], rPosTrBin[1], nPosTrBin, rPosTrBin[0], rPosTrBin[1]);
  hECalTruClustVsParEne      = new TH2D("hECalTruClustVsParEne",      "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // truth hcal cluster debug histograms
  hHCalDebugTruClustSum5     = new TH1D("hHCalDebugTruClustSum5",     "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugTruClustSum10    = new TH1D("hHCalDebugTruClustSum10",    "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugTruClustSum100   = new TH1D("hHCalDebugTruClustSum100",   "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugTruClustSum1000  = new TH1D("hHCalDebugTruClustSum1000",  "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hHCalDebugTruClustDiff5    = new TH1D("hHCalDebugTruClustDiff5",    "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalDebugTruClustDiff10   = new TH1D("hHCalDebugTruClustDiff10",   "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalDebugTruClustDiff100  = new TH1D("hHCalDebugTruClustDiff100",  "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hHCalDebugTruClustDiff1000 = new TH1D("hHCalDebugTruClustDiff1000", "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  // truth ecal cluster debug histograms
  hECalDebugTruClustSum5     = new TH1D("hECalDebugTruClustSum5",     "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hECalDebugTruClustSum10    = new TH1D("hECalDebugTruClustSum10",    "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hECalDebugTruClustSum100   = new TH1D("hECalDebugTruClustSum100",   "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hECalDebugTruClustSum1000  = new TH1D("hECalDebugTruClustSum1000",  "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hECalDebugTruClustDiff5    = new TH1D("hECalDebugTruClustDiff5",    "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hECalDebugTruClustDiff10   = new TH1D("hECalDebugTruClustDiff10",   "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hECalDebugTruClustDiff100  = new TH1D("hECalDebugTruClustDiff100",  "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hECalDebugTruClustDiff1000 = new TH1D("hECalDebugTruClustDiff1000", "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  // event-wise hcal histograms
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
  // event-wise ecal histograms
  hEvtECalNumPar             = new TH1I("hEvtECalNumPar",             "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtECalNumHit             = new TH1I("hEvtECalNumHit",             "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtECalNumClust           = new TH1I("hEvtECalNumClust",           "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtECalNumTruClust        = new TH1I("hEvtECalNumTruClust",        "", nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtECalSumHitEne          = new TH1D("hEvtECalSumHitEne",          "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtECalSumClustEne        = new TH1D("hEvtECalSumClustEne",        "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtECalSumTruClustEne     = new TH1D("hEvtECalSumTruClustEne",     "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtECalLeadClustEne       = new TH1D("hEvtECalLeadClustEne",       "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtECalLeadTruClustEne    = new TH1D("hEvtECalLeadTruClustEne",    "", nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtECalSumHitDiff         = new TH1D("hEvtECalSumHitDiff",         "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtECalSumClustDiff       = new TH1D("hEvtECalSumClustDiff",       "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtECalSumTruClustDiff    = new TH1D("hEvtECalSumTruClustDiff",    "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtECalLeadClustDiff      = new TH1D("hEvtECalLeadClustDiff",      "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtECalLeadTruClustDiff   = new TH1D("hEvtECalLeadTruClustDiff",   "", nDiffBin,  rDiffBin[0],  rDiffBin[1]);
  hEvtECalNumClustVsHit      = new TH2I("hEvtECalNumClustVsHit",      "", nNumBin,   rNumBin[0],   rNumBin[1],   nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtECalNumTruClustVsClust = new TH2I("hEvtECalNumTruClustVsClust", "", nNumBin,   rNumBin[0],   rNumBin[1],   nNumBin,   rNumBin[0],   rNumBin[1]);
  hEvtECalSumHitVsPar        = new TH2D("hEvtECalSumHitVsPar",        "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtECalSumClustVsPar      = new TH2D("hEvtECalSumClustVsPar",      "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtECalSumTruClustVsPar   = new TH2D("hEvtECalSumTruClustVsPar",   "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtECalLeadClustVsPar     = new TH2D("hEvtECalLeadClustVsPar",     "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  hEvtECalLeadTruClustVsPar  = new TH2D("hEvtECalLeadTruClustVsPar",  "", nEneBin,   rEneBin[0],   rEneBin[1],   nEneBin,   rEneBin[0],   rEneBin[1]);
  // errors
  hParChrg                   -> Sumw2();
  hParMass                   -> Sumw2();
  hParEne                    -> Sumw2();
  hParMom                    -> Sumw2();
  hParMomX                   -> Sumw2();
  hParMomY                   -> Sumw2();
  hParMomZ                   -> Sumw2();
  hHCalRecHitEne             -> Sumw2();
  hHCalRecHitPosZ            -> Sumw2();
  hHCalRecHitParDiff         -> Sumw2();
  hHCalRecHitPosYvsX         -> Sumw2();
  hHCalRecHitVsParEne        -> Sumw2();
  hECalRecHitEne             -> Sumw2();
  hECalRecHitPosZ            -> Sumw2();
  hECalRecHitParDiff         -> Sumw2();
  hECalRecHitPosYvsX         -> Sumw2();
  hECalRecHitVsParEne        -> Sumw2();
  hHCalClustEne              -> Sumw2();
  hHCalClustPosZ             -> Sumw2();
  hHCalClustNumHit           -> Sumw2();
  hHCalClustParDiff          -> Sumw2();
  hHCalClustPosYvsX          -> Sumw2();
  hHCalClustVsParEne         -> Sumw2();
  hECalClustEne              -> Sumw2();
  hECalClustPosZ             -> Sumw2();
  hECalClustNumHit           -> Sumw2();
  hECalClustParDiff          -> Sumw2();
  hECalClustPosYvsX          -> Sumw2();
  hECalClustVsParEne         -> Sumw2();
  hHCalDebugClustSum5        -> Sumw2();
  hHCalDebugClustSum10       -> Sumw2();
  hHCalDebugClustSum100      -> Sumw2();
  hHCalDebugClustSum1000     -> Sumw2();
  hHCalDebugClustDiff5       -> Sumw2();
  hHCalDebugClustDiff10      -> Sumw2();
  hHCalDebugClustDiff100     -> Sumw2();
  hHCalDebugClustDiff1000    -> Sumw2();
  hECalDebugClustSum5        -> Sumw2();
  hECalDebugClustSum10       -> Sumw2();
  hECalDebugClustSum100      -> Sumw2();
  hECalDebugClustSum1000     -> Sumw2();
  hECalDebugClustDiff5       -> Sumw2();
  hECalDebugClustDiff10      -> Sumw2();
  hECalDebugClustDiff100     -> Sumw2();
  hECalDebugClustDiff1000    -> Sumw2();
  hHCalTruClustEne           -> Sumw2();
  hHCalTruClustPosZ          -> Sumw2();
  hHCalTruClustNumHit        -> Sumw2();
  hHCalTruClustParDiff       -> Sumw2();
  hHCalTruClustPosYvsX       -> Sumw2();
  hHCalTruClustVsParEne      -> Sumw2();
  hECalTruClustEne           -> Sumw2();
  hECalTruClustPosZ          -> Sumw2();
  hECalTruClustNumHit        -> Sumw2();
  hECalTruClustParDiff       -> Sumw2();
  hECalTruClustPosYvsX       -> Sumw2();
  hECalTruClustVsParEne      -> Sumw2();
  hHCalDebugTruClustSum5     -> Sumw2();
  hHCalDebugTruClustSum10    -> Sumw2();
  hHCalDebugTruClustSum100   -> Sumw2();
  hHCalDebugTruClustSum1000  -> Sumw2();
  hHCalDebugTruClustDiff5    -> Sumw2();
  hHCalDebugTruClustDiff10   -> Sumw2();
  hHCalDebugTruClustDiff100  -> Sumw2();
  hHCalDebugTruClustDiff1000 -> Sumw2();
  hECalDebugTruClustSum5     -> Sumw2();
  hECalDebugTruClustSum10    -> Sumw2();
  hECalDebugTruClustSum100   -> Sumw2();
  hECalDebugTruClustSum1000  -> Sumw2();
  hECalDebugTruClustDiff5    -> Sumw2();
  hECalDebugTruClustDiff10   -> Sumw2();
  hECalDebugTruClustDiff100  -> Sumw2();
  hECalDebugTruClustDiff1000 -> Sumw2();
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
  hEvtECalNumPar             -> Sumw2();
  hEvtECalNumHit             -> Sumw2();
  hEvtECalNumClust           -> Sumw2();
  hEvtECalNumTruClust        -> Sumw2();
  hEvtECalSumHitEne          -> Sumw2();
  hEvtECalSumClustEne        -> Sumw2();
  hEvtECalSumTruClustEne     -> Sumw2();
  hEvtECalLeadClustEne       -> Sumw2();
  hEvtECalLeadTruClustEne    -> Sumw2();
  hEvtECalSumHitDiff         -> Sumw2();
  hEvtECalSumClustDiff       -> Sumw2();
  hEvtECalSumTruClustDiff    -> Sumw2();
  hEvtECalLeadClustDiff      -> Sumw2();
  hEvtECalLeadTruClustDiff   -> Sumw2();
  hEvtECalNumClustVsHit      -> Sumw2();
  hEvtECalNumTruClustVsClust -> Sumw2();
  hEvtECalSumHitVsPar        -> Sumw2();
  hEvtECalSumClustVsPar      -> Sumw2();
  hEvtECalLeadClustVsPar     -> Sumw2();
  hEvtECalLeadTruClustVsPar  -> Sumw2();
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

  // ecal hit energy
  for (auto beCalHit : becalRecHits()) {
    eECalHitSum += beCalHit -> getEnergy();
  }  // end 1st ecal hit loop

  // if both hit sums are 0, skip event
  const bool isHCalHitSumNonzero = (eHCalHitSum > 0.);
  const bool isECalHitSumNonzero = (eECalHitSum > 0.);
  const bool areHitSumsNonzero   = (isHCalHitSumNonzero || isECalHitSumNonzero);
  if (!areHitSumsNonzero) {
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
  unsigned long nHCalHit(0);
  for (auto bhCalHit : bhcalRecHits()) {

    // grab hit properties
    const auto rHCalHitX   = bhCalHit -> getPosition().x;
    const auto rHCalHitY   = bhCalHit -> getPosition().y;
    const auto rHCalHitZ   = bhCalHit -> getPosition().z;
    const auto eHCalHit    = bhCalHit -> getEnergy();
    const auto diffHCalHit = (eHCalHit - eMcPar) / eHCalHit;

    // fill hit histograms and increment sums/counters
    hHCalRecHitEne      -> Fill(eHCalHit);
    hHCalRecHitPosZ     -> Fill(rHCalHitZ);
    hHCalRecHitParDiff  -> Fill(diffHCalHit);
    hHCalRecHitPosYvsX  -> Fill(rHCalHitX, rHCalHitY);
    hHCalRecHitVsParEne -> Fill(eMcPar, eHCalHit);
    ++nHCalHit;
  }  // end 2nd hcal hit loop

  // reco. ecal hit loop
  unsigned long nECalHit(0);
  for (auto beCalHit : becalRecHits()) {

    // grab hit properties
    const auto rECalHitX   = beCalHit -> getPosition().x;
    const auto rECalHitY   = beCalHit -> getPosition().y;
    const auto rECalHitZ   = beCalHit -> getPosition().z;
    const auto eECalHit    = beCalHit -> getEnergy();
    const auto diffECalHit = (eECalHit - eMcPar) / eECalHit;

    // fill hit histograms and increment sums/counters
    hECalRecHitEne      -> Fill(eECalHit);
    hECalRecHitPosZ     -> Fill(rECalHitZ);
    hECalRecHitParDiff  -> Fill(diffECalHit);
    hECalRecHitPosYvsX  -> Fill(rECalHitX, rECalHitY);
    hECalRecHitVsParEne -> Fill(eMcPar, eECalHit);
    ++nECalHit;
  }  // end 2nd hcal hit loop

  // for highest energy clusters
  int    iLeadHCalClust(-1);
  int    iLeadECalClust(-1);
  int    iLeadTruHCalClust(-1);
  int    iLeadTruECalClust(-1);
  double eLeadHCalClust(0.);
  double eLeadECalClust(0.);
  double eLeadTruHCalClust(0.);
  double eLeadTruECalClust(0.);
  double diffLeadHCalClust(0.);
  double diffLeadECalClust(0.);
  double diffLeadTruHCalClust(0.);
  double diffLeadTruECalClust(0.);

  // reco. hcal cluster loop
  unsigned long iHCalClust(0);
  unsigned long nHCalClust(0);
  for (auto bhCalClust : bhcalClusters()) {

    // grab cluster properties
    const auto rHCalClustX   = bhCalClust -> getPosition().x;
    const auto rHCalClustY   = bhCalClust -> getPosition().y;
    const auto rHCalClustZ   = bhCalClust -> getPosition().z;
    const auto eHCalClust    = bhCalClust -> getEnergy();
    const auto nHitHCalClust = bhCalClust -> hits_size();
    const auto diffHCalClust = (eHCalClust - eMcPar) / eHCalClust;

    // fill cluster histograms and increment counters
    hHCalClustEne      -> Fill(eHCalClust);
    hHCalClustPosZ     -> Fill(rHCalClustZ);
    hHCalClustNumHit   -> Fill(nHitHCalClust);
    hHCalClustParDiff  -> Fill(diffHCalClust);
    hHCalClustPosYvsX  -> Fill(rHCalClustX, rHCalClustY);
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

  // reco. ecal cluster loop
  unsigned long iECalClust(0);
  unsigned long nECalClust(0);
  for (auto beCalClust : becalClusters()) {

    // grab cluster properties
    const auto rECalClustX   = beCalClust -> getPosition().x;
    const auto rECalClustY   = beCalClust -> getPosition().y;
    const auto rECalClustZ   = beCalClust -> getPosition().z;
    const auto eECalClust    = beCalClust -> getEnergy();
    const auto nHitECalClust = beCalClust -> hits_size();
    const auto diffECalClust = (eECalClust - eMcPar) / eECalClust;

    // fill cluster histograms and increment counters
    hECalClustEne      -> Fill(eECalClust);
    hECalClustPosZ     -> Fill(rECalClustZ);
    hECalClustNumHit   -> Fill(nHitECalClust);
    hECalClustParDiff  -> Fill(diffECalClust);
    hECalClustPosYvsX  -> Fill(rECalClustX, rECalClustY);
    hECalClustVsParEne -> Fill(eMcPar, eECalClust);
    eECalClustSum += eECalClust;
    ++nECalClust;

    // select leading cluster
    const bool isBiggerEne = (eECalClust > eLeadECalClust);
    if (isBiggerEne) {
      iLeadECalClust    = iECalClust;
      eLeadECalClust    = eECalClust;
      diffLeadECalClust = diffECalClust;
    }
    ++iECalClust;
  }  // end reco. ecal cluster loop

  // for debugging reco. clusters
  double eDebugSumHCalClust5(eLeadHCalClust);
  double eDebugSumECalClust5(eLeadECalClust);
  double eDebugSumHCalClust10(eLeadHCalClust);
  double eDebugSumECalClust10(eLeadECalClust);
  double eDebugSumHCalClust100(eLeadHCalClust);
  double eDebugSumECalClust100(eLeadECalClust);
  double eDebugSumHCalClust1000(eLeadHCalClust);
  double eDebugSumECalClust1000(eLeadECalClust);

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

  // debug reco. ecal cluster loop
  unsigned long iDebugECalClust(0);
  for (auto debugECalClust : becalClusters()) {

    // select leading cluster
    const bool isLeadECalClust = (iDebugECalClust == iLeadECalClust);
    if (isLeadECalClust) {

      // grab lead cluster properties
      const auto rLeadECalClustX = debugECalClust -> getPosition().x;
      const auto rLeadECalClustY = debugECalClust -> getPosition().y;
      const auto eLeadECalClust  = debugECalClust -> getEnergy();

      unsigned long iOtherECalClust(0);
      double        toAddECalClust5(0.);
      double        toAddECalClust10(0.);
      double        toAddECalClust100(0.);
      double        toAddECalClust1000(0.);
      for (auto otherECalClust : becalClusters()) {

        // ignore same cluster
        const bool isSameECalClust = (iOtherECalClust == iDebugECalClust);
        if (isSameECalClust) continue;

        // grab other cluster properties
        const auto rOtherECalClustX = otherECalClust -> getPosition().x;
        const auto rOtherECalClustY = otherECalClust -> getPosition().y;
        const auto eOtherECalClust  = otherECalClust -> getEnergy();
        const auto drLeadOtherX     = rOtherECalClustX - rLeadECalClustX;
        const auto drLeadOtherY     = rOtherECalClustY - rLeadECalClustY;
        const auto drLeadOther      = std::sqrt((drLeadOtherX * drLeadOtherX) + (drLeadOtherY * drLeadOtherY));

        // increment relevant sums and counters
        const bool isIn5mm    = (drLeadOther < 5.);
        const bool isIn10mm   = (drLeadOther < 10.);
        const bool isIn100mm  = (drLeadOther < 100.);
        const bool isIn1000mm = (drLeadOther < 1000.);
        if (isIn5mm)    toAddECalClust5    += eOtherECalClust;
        if (isIn10mm)   toAddECalClust10   += eOtherECalClust;
        if (isIn100mm)  toAddECalClust100  += eOtherECalClust;
        if (isIn1000mm) toAddECalClust1000 += eOtherECalClust;
        ++iOtherECalClust;
      }  // end other reco. cluster loo

      // add sums to lead energy
      eDebugSumECalClust5    += toAddECalClust5;
      eDebugSumECalClust10   += toAddECalClust10;
      eDebugSumECalClust100  += toAddECalClust100;
      eDebugSumECalClust1000 += toAddECalClust1000;
    }
    ++iDebugECalClust;
  }  // end debug reco. ecal cluster loop

  // true hcal cluster loop
  unsigned long iTruHCalClust(0);
  unsigned long nTruHCalClust(0);
  for (auto truthHCalClust : bhcalTruthClusters()) {

    // grab cluster properties
    const auto rTruHCalClustX   = truthHCalClust -> getPosition().x;
    const auto rTruHCalClustY   = truthHCalClust -> getPosition().y;
    const auto rTruHCalClustZ   = truthHCalClust -> getPosition().z;
    const auto eTruHCalClust    = truthHCalClust -> getEnergy();
    const auto nHitTruHCalClust = truthHCalClust -> hits_size();
    const auto diffTruHCalClust = (eTruHCalClust - eMcPar) / eTruHCalClust;

    // fill cluster histograms and increment counters
    hHCalTruClustEne      -> Fill(eTruHCalClust);
    hHCalTruClustPosZ     -> Fill(rTruHCalClustZ);
    hHCalTruClustNumHit   -> Fill(nHitTruHCalClust);
    hHCalTruClustParDiff  -> Fill(diffTruHCalClust);
    hHCalTruClustPosYvsX  -> Fill(rTruHCalClustX, rTruHCalClustY);
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

  // true ecal cluster loop
  unsigned long iTruECalClust(0);
  unsigned long nTruECalClust(0);
  for (auto truthECalClust : becalTruthClusters()) {

    // grab cluster properties
    const auto rTruECalClustX   = truthECalClust -> getPosition().x;
    const auto rTruECalClustY   = truthECalClust -> getPosition().y;
    const auto rTruECalClustZ   = truthECalClust -> getPosition().z;
    const auto eTruECalClust    = truthECalClust -> getEnergy();
    const auto nHitTruECalClust = truthECalClust -> hits_size();
    const auto diffTruECalClust = (eTruECalClust - eMcPar) / eTruECalClust;

    // fill cluster histograms and increment counters
    hECalTruClustEne      -> Fill(eTruECalClust);
    hECalTruClustPosZ     -> Fill(rTruECalClustZ);
    hECalTruClustNumHit   -> Fill(nHitTruECalClust);
    hECalTruClustParDiff  -> Fill(diffTruECalClust);
    hECalTruClustPosYvsX  -> Fill(rTruECalClustX, rTruECalClustY);
    hECalTruClustVsParEne -> Fill(eMcPar, eTruECalClust);
    eTruECalClustSum += eTruECalClust;
    ++nTruECalClust;

    // select leading cluster
    const bool isBiggerEne = (eTruECalClust > eLeadTruECalClust);
    if (isBiggerEne) {
      iLeadTruECalClust    = iTruECalClust;
      eLeadTruECalClust    = eTruECalClust;
      diffLeadTruECalClust = diffTruECalClust;
    }
    ++iTruECalClust;
  }  // end true ecal cluster loop

  // for debugging truth clusters
  double eDebugSumTruHCalClust5(eLeadTruHCalClust);
  double eDebugSumTruECalClust5(eLeadTruECalClust);
  double eDebugSumTruHCalClust10(eLeadTruHCalClust);
  double eDebugSumTruECalClust10(eLeadTruECalClust);
  double eDebugSumTruHCalClust100(eLeadTruHCalClust);
  double eDebugSumTruECalClust100(eLeadTruECalClust);
  double eDebugSumTruHCalClust1000(eLeadTruHCalClust);
  double eDebugSumTruECalClust1000(eLeadTruECalClust);

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

  // debug truth ecal cluster loop
  unsigned long iDebugTruECalClust(0);
  for (auto debugTruthECalClust : becalTruthClusters()) {

    // select leading cluster
    const bool isLeadTruECalClust = (iDebugTruECalClust == iLeadTruECalClust);
    if (isLeadTruECalClust) {

      // grab lead cluster properties
      const auto rLeadTruECalClustX = debugTruthECalClust -> getPosition().x;
      const auto rLeadTruECalClustY = debugTruthECalClust -> getPosition().y;
      const auto eLeadTruECalClust  = debugTruthECalClust -> getEnergy();

      unsigned long iOtherTruECalClust(0);
      double        toAddTruECalClust5(0.);
      double        toAddTruECalClust10(0.);
      double        toAddTruECalClust100(0.);
      double        toAddTruECalClust1000(0.);
      for (auto otherTruthECalClust : becalTruthClusters()) {

        // ignore same cluster
        const bool isSameTruECalClust = (iOtherTruECalClust == iDebugTruECalClust);
        if (isSameTruECalClust) continue;

        // grab other cluster properties
        const auto rOtherTruECalClustX = otherTruthECalClust -> getPosition().x;
        const auto rOtherTruECalClustY = otherTruthECalClust -> getPosition().y;
        const auto eOtherTruECalClust  = otherTruthECalClust -> getEnergy();
        const auto drLeadOtherX        = rOtherTruECalClustX - rLeadTruECalClustX;
        const auto drLeadOtherY        = rOtherTruECalClustY - rLeadTruECalClustY;
        const auto drLeadOther         = std::sqrt((drLeadOtherX * drLeadOtherX) + (drLeadOtherY * drLeadOtherY));

        // increment relevant sums and counters
        const bool isIn5mm    = (drLeadOther < 5);
        const bool isIn10mm   = (drLeadOther < 10);
        const bool isIn100mm  = (drLeadOther < 100);
        const bool isIn1000mm = (drLeadOther < 1000);
        if (isIn5mm)    toAddTruECalClust5    += eOtherTruECalClust;
        if (isIn10mm)   toAddTruECalClust10   += eOtherTruECalClust;
        if (isIn100mm)  toAddTruECalClust100  += eOtherTruECalClust;
        if (isIn1000mm) toAddTruECalClust1000 += eOtherTruECalClust;
        ++iOtherTruECalClust;
      }  // end other true ecal cluster loo

      // add sums to lead energy
      eDebugSumTruECalClust5    += toAddTruECalClust5;
      eDebugSumTruECalClust10   += toAddTruECalClust10;
      eDebugSumTruECalClust100  += toAddTruECalClust100;
      eDebugSumTruECalClust1000 += toAddTruECalClust1000;
    }
    ++iDebugTruECalClust;
  }  // end debug truth ecal cluster loop

  // do event-wise calculations
  const auto diffHCalHitSum      = (eHCalHitSum - eMcPar) / eHCalHitSum;
  const auto diffECalHitSum      = (eECalHitSum - eMcPar) / eECalHitSum;
  const auto diffHCalClustSum    = (eHCalClustSum - eMcPar) / eHCalClustSum;
  const auto diffECalClustSum    = (eECalClustSum - eMcPar) / eECalClustSum;
  const auto diffTruHCalClustSum = (eTruHCalClustSum - eMcPar) / eTruHCalClustSum;
  const auto diffTruECalClustSum = (eTruECalClustSum - eMcPar) / eTruECalClustSum;

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

  // fill event-wise ecal histograms
  hEvtECalNumPar             -> Fill(nPar);
  hEvtECalNumHit             -> Fill(nECalHit);
  hEvtECalNumClust           -> Fill(nECalClust);
  hEvtECalNumTruClust        -> Fill(nTruECalClust);
  hEvtECalSumHitEne          -> Fill(eECalHitSum);
  hEvtECalSumClustEne        -> Fill(eECalClustSum);
  hEvtECalSumTruClustEne     -> Fill(eTruECalClustSum);
  hEvtECalLeadClustEne       -> Fill(eLeadECalClust);
  hEvtECalLeadTruClustEne    -> Fill(eLeadTruECalClust);
  hEvtECalSumHitDiff         -> Fill(diffECalHitSum);
  hEvtECalSumClustDiff       -> Fill(diffECalClustSum);
  hEvtECalSumTruClustDiff    -> Fill(diffTruECalClustSum);
  hEvtECalLeadClustDiff      -> Fill(diffLeadECalClust);
  hEvtECalLeadTruClustDiff   -> Fill(diffLeadTruECalClust);
  hEvtECalNumClustVsHit      -> Fill(nECalHit, nECalClust);
  hEvtECalNumTruClustVsClust -> Fill(nECalClust, nTruECalClust);
  hEvtECalSumHitVsPar        -> Fill(eMcPar, eECalHitSum);
  hEvtECalSumClustVsPar      -> Fill(eMcPar, eECalClustSum);
  hEvtECalSumTruClustVsPar   -> Fill(eMcPar, eTruECalClustSum);
  hEvtECalLeadClustVsPar     -> Fill(eMcPar, eLeadECalClust);
  hEvtECalLeadTruClustVsPar  -> Fill(eMcPar, eLeadTruECalClust);

  // do debugging calculations
  const auto diffDebugHCalClustSum5       = (eDebugSumHCalClust5 - eMcPar) / eDebugSumHCalClust5;
  const auto diffDebugECalClustSum5       = (eDebugSumECalClust5 - eMcPar) / eDebugSumECalClust5;
  const auto diffDebugHCalClustSum10      = (eDebugSumHCalClust10 - eMcPar) / eDebugSumHCalClust5;
  const auto diffDebugECalClustSum10      = (eDebugSumECalClust10 - eMcPar) / eDebugSumECalClust5;
  const auto diffDebugHCalClustSum100     = (eDebugSumHCalClust100 - eMcPar) / eDebugSumHCalClust5;
  const auto diffDebugECalClustSum100     = (eDebugSumECalClust100 - eMcPar) / eDebugSumECalClust5;
  const auto diffDebugHCalClustSum1000    = (eDebugSumHCalClust1000 - eMcPar) / eDebugSumHCalClust5;
  const auto diffDebugECalClustSum1000    = (eDebugSumECalClust1000 - eMcPar) / eDebugSumECalClust5;
  const auto diffDebugTruHCalClustSum5    = (eDebugSumTruHCalClust5 - eMcPar) / eDebugSumTruHCalClust5;
  const auto diffDebugTruECalClustSum5    = (eDebugSumTruECalClust5 - eMcPar) / eDebugSumTruECalClust5;
  const auto diffDebugTruHCalClustSum10   = (eDebugSumTruHCalClust10 - eMcPar) / eDebugSumTruHCalClust5;
  const auto diffDebugTruECalClustSum10   = (eDebugSumTruECalClust10 - eMcPar) / eDebugSumTruECalClust5;
  const auto diffDebugTruHCalClustSum100  = (eDebugSumTruHCalClust100 - eMcPar) / eDebugSumTruHCalClust5;
  const auto diffDebugTruECalClustSum100  = (eDebugSumTruECalClust100 - eMcPar) / eDebugSumTruECalClust5;
  const auto diffDebugTruHCalClustSum1000 = (eDebugSumTruHCalClust1000 - eMcPar) / eDebugSumTruHCalClust5;
  const auto diffDebugTruECalClustSum1000 = (eDebugSumTruECalClust1000 - eMcPar) / eDebugSumTruECalClust5;

  // fill reco. cluster hcal debug histograms
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

  // fill reco. cluster ecal debug histograms
  hECalDebugClustSum5        -> Fill(eDebugSumECalClust5);
  hECalDebugClustSum10       -> Fill(eDebugSumECalClust10);
  hECalDebugClustSum100      -> Fill(eDebugSumECalClust100);
  hECalDebugClustSum1000     -> Fill(eDebugSumECalClust1000);
  hECalDebugTruClustSum5     -> Fill(eDebugSumECalClust5);
  hECalDebugTruClustSum10    -> Fill(eDebugSumECalClust10);
  hECalDebugTruClustSum100   -> Fill(eDebugSumECalClust100);
  hECalDebugTruClustSum1000  -> Fill(eDebugSumECalClust1000);
  hECalDebugClustDiff5       -> Fill(diffDebugECalClustSum5);
  hECalDebugClustDiff10      -> Fill(diffDebugECalClustSum10);
  hECalDebugClustDiff100     -> Fill(diffDebugECalClustSum100);
  hECalDebugClustDiff1000    -> Fill(diffDebugECalClustSum1000);
  hECalDebugTruClustDiff5    -> Fill(diffDebugTruECalClustSum5);
  hECalDebugTruClustDiff10   -> Fill(diffDebugTruECalClustSum10);
  hECalDebugTruClustDiff100  -> Fill(diffDebugTruECalClustSum100);
  hECalDebugTruClustDiff1000 -> Fill(diffDebugTruECalClustSum1000);
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
  // set reco. hit ecal axis titles
  hECalRecHitEne            -> GetXaxis() -> SetTitle(sEneHit.Data());
  hECalRecHitEne            -> GetYaxis() -> SetTitle(sCount.Data());
  hECalRecHitPosZ           -> GetXaxis() -> SetTitle(sPosHitZ.Data());
  hECalRecHitPosZ           -> GetYaxis() -> SetTitle(sCount.Data());
  hECalRecHitParDiff        -> GetXaxis() -> SetTitle(sEneHitDiff.Data());
  hECalRecHitParDiff        -> GetYaxis() -> SetTitle(sCount.Data());
  hECalRecHitPosYvsX        -> GetXaxis() -> SetTitle(sPosHitX.Data());
  hECalRecHitPosYvsX        -> GetYaxis() -> SetTitle(sPosHitY.Data());
  hECalRecHitPosYvsX        -> GetZaxis() -> SetTitle(sCount.Data());
  hECalRecHitVsParEne       -> GetXaxis() -> SetTitle(sEnePar.Data());
  hECalRecHitVsParEne       -> GetYaxis() -> SetTitle(sEneHit.Data());
  hECalRecHitVsParEne       -> GetZaxis() -> SetTitle(sCount.Data());
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
  // set reco. cluster ecal axis titles
  hECalClustEne             -> GetXaxis() -> SetTitle(sEneClust.Data());
  hECalClustEne             -> GetYaxis() -> SetTitle(sCount.Data());
  hECalClustPosZ            -> GetXaxis() -> SetTitle(sPosClustZ.Data());
  hECalClustPosZ            -> GetYaxis() -> SetTitle(sCount.Data());
  hECalClustNumHit          -> GetXaxis() -> SetTitle(sNumHitClust.Data());
  hECalClustNumHit          -> GetYaxis() -> SetTitle(sCount.Data());
  hECalClustParDiff         -> GetXaxis() -> SetTitle(sEneClustDiff.Data());
  hECalClustParDiff         -> GetYaxis() -> SetTitle(sCount.Data());
  hECalClustPosYvsX         -> GetXaxis() -> SetTitle(sPosClustX.Data());
  hECalClustPosYvsX         -> GetYaxis() -> SetTitle(sPosClustY.Data());
  hECalClustPosYvsX         -> GetZaxis() -> SetTitle(sCount.Data());
  hECalClustVsParEne        -> GetXaxis() -> SetTitle(sEnePar.Data());
  hECalClustVsParEne        -> GetYaxis() -> SetTitle(sEneClust.Data());
  hECalClustVsParEne        -> GetZaxis() -> SetTitle(sCount.Data());
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
  // set truth cluster ecal axis titles
  hECalTruClustEne          -> GetXaxis() -> SetTitle(sEneTruClust.Data());
  hECalTruClustEne          -> GetYaxis() -> SetTitle(sCount.Data());
  hECalTruClustPosZ         -> GetXaxis() -> SetTitle(sPosTruClustZ.Data());
  hECalTruClustPosZ         -> GetYaxis() -> SetTitle(sCount.Data());
  hECalTruClustNumHit       -> GetXaxis() -> SetTitle(sNumHitTruClust.Data());
  hECalTruClustNumHit       -> GetYaxis() -> SetTitle(sCount.Data());
  hECalTruClustParDiff      -> GetXaxis() -> SetTitle(sEneTruClustDiff.Data());
  hECalTruClustParDiff      -> GetYaxis() -> SetTitle(sCount.Data());
  hECalTruClustPosYvsX      -> GetXaxis() -> SetTitle(sPosTruClustX.Data());
  hECalTruClustPosYvsX      -> GetYaxis() -> SetTitle(sPosTruClustY.Data());
  hECalTruClustPosYvsX      -> GetZaxis() -> SetTitle(sCount.Data());
  hECalTruClustVsParEne     -> GetXaxis() -> SetTitle(sEnePar.Data());
  hECalTruClustVsParEne     -> GetYaxis() -> SetTitle(sEneTruClust.Data());
  hECalTruClustVsParEne     -> GetZaxis() -> SetTitle(sCount.Data());
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
  // set event-wise ecal axis titles
  hEvtECalNumPar            -> GetXaxis() -> SetTitle(sNumParEvt.Data());
  hEvtECalNumPar            -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalNumHit            -> GetXaxis() -> SetTitle(sNumHitEvt.Data());
  hEvtECalNumHit            -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalNumClust          -> GetXaxis() -> SetTitle(sNumClustEvt.Data());
  hEvtECalNumClust          -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalNumTruClust       -> GetXaxis() -> SetTitle(sNumTruClustEvt.Data());
  hEvtECalNumTruClust       -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalSumHitEne         -> GetXaxis() -> SetTitle(sEneHitSum.Data());
  hEvtECalSumHitEne         -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalSumClustEne       -> GetXaxis() -> SetTitle(sEneClustSum.Data());
  hEvtECalSumClustEne       -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalSumTruClustEne    -> GetXaxis() -> SetTitle(sEneTruClustSum.Data());
  hEvtECalSumTruClustEne    -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalLeadClustEne      -> GetXaxis() -> SetTitle(sEneClustLead.Data());
  hEvtECalLeadClustEne      -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalLeadTruClustEne   -> GetXaxis() -> SetTitle(sEneTruClustLead.Data());
  hEvtECalLeadTruClustEne   -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalSumHitDiff        -> GetXaxis() -> SetTitle(sEneHitSumDiff.Data());
  hEvtECalSumHitDiff        -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalSumClustDiff      -> GetXaxis() -> SetTitle(sEneClustSumDiff.Data());
  hEvtECalSumClustDiff      -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalSumTruClustDiff   -> GetXaxis() -> SetTitle(sEneTruClustSumDiff.Data());
  hEvtECalSumTruClustDiff   -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalLeadClustDiff     -> GetXaxis() -> SetTitle(sEneClustLeadDiff.Data());
  hEvtECalLeadClustDiff     -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalLeadTruClustDiff  -> GetXaxis() -> SetTitle(sEneTruClustLeadDiff.Data());
  hEvtECalLeadTruClustDiff  -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtECalSumHitVsPar       -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtECalSumHitVsPar       -> GetYaxis() -> SetTitle(sEneHitSum.Data());
  hEvtECalSumHitVsPar       -> GetZaxis() -> SetTitle(sCount.Data());
  hEvtECalSumClustVsPar     -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtECalSumClustVsPar     -> GetYaxis() -> SetTitle(sEneClustSum.Data());
  hEvtECalSumClustVsPar     -> GetZaxis() -> SetTitle(sCount.Data());
  hEvtECalSumTruClustVsPar  -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtECalSumTruClustVsPar  -> GetYaxis() -> SetTitle(sEneTruClustSum.Data());
  hEvtECalSumTruClustVsPar  -> GetZaxis() -> SetTitle(sCount.Data());
  hEvtECalLeadClustVsPar    -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtECalLeadClustVsPar    -> GetYaxis() -> SetTitle(sEneClustLead.Data());
  hEvtECalLeadClustVsPar    -> GetZaxis() -> SetTitle(sCount.Data());
  hEvtECalLeadTruClustVsPar -> GetXaxis() -> SetTitle(sEnePar.Data());
  hEvtECalLeadTruClustVsPar -> GetYaxis() -> SetTitle(sEneTruClustLead.Data());
  hEvtECalLeadTruClustVsPar -> GetZaxis() -> SetTitle(sCount.Data());
  return;

}  // end 'FinishWithGlobalRootLock()'

// end ------------------------------------------------------------------------
