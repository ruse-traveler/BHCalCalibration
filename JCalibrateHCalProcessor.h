// 'JCalibrateHCalProcessor.h'
// Derek Anderson
// 11.02.2022
//
// A simple JANA plugin to compare the
// reconstructed hit and cluster energy
// in the HCal to simulated particles.

// C includes
#include <cmath>
// ROOT includes
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TString.h>
// JANA includes
#include <JANA/JEventProcessorSequentialRoot.h>
// EDM includes
#include <edm4eic/CalorimeterHit.h>
#include <edm4eic/ReconstructedParticle.h>
#include <edm4eic/Cluster.h>

// global constants
static const size_t NRange(2);
static const size_t NComp(3);
static const float  CPar(1.);
static const float  MParMin(0.135);
static const float  MParMax(0.145);
static const float  PParMin(1.9);
static const float  PParMax(2.1);



class JCalibrateHCalProcessor : public JEventProcessorSequentialRoot {

  private:

    // Data objects we will need from JANA
    PrefetchT<edm4eic::ReconstructedParticle> genParticles       = {this, "GeneratedParticles"};
    PrefetchT<edm4eic::CalorimeterHit>        bhcalRecHits       = {this, "HcalBarrelRecHits"};
    PrefetchT<edm4eic::Cluster>               bhcalClusters      = {this, "HcalBarrelClusters"};
    PrefetchT<edm4eic::Cluster>               bhcalTruthClusters = {this, "HcalBarrelTruthClusters"};

    // particle histograms
    TH1D *hParChrg               = nullptr;
    TH1D *hParMass               = nullptr;
    TH1D *hParEne                = nullptr;
    TH1D *hParMom                = nullptr;
    TH1D *hParMomX               = nullptr;
    TH1D *hParMomY               = nullptr;
    TH1D *hParMomZ               = nullptr;
    // reconstructed hit histograms
    TH1D *hRecHitEne             = nullptr;
    TH1D *hRecHitPosZ            = nullptr;
    TH1D *hRecHitParDiff         = nullptr;
    TH2D *hRecHitPosYvsX         = nullptr;
    TH2D *hRecHitVsParEne        = nullptr;
    // reconstructed cluster histograms
    TH1D *hClusterEne            = nullptr;
    TH1D *hClusterPosZ           = nullptr;
    TH1I *hClusterNumHit         = nullptr;
    TH1D *hClusterParDiff        = nullptr;
    TH2D *hClusterPosYvsX        = nullptr;
    TH2D *hClusterVsParEne       = nullptr;
    // reconstructed cluster debug histograms
    TH1D *hDebugClusterSum5      = nullptr;
    TH1D *hDebugClusterSum10     = nullptr;
    TH1D *hDebugClusterSum100    = nullptr;
    TH1D *hDebugClusterSum1000   = nullptr;
    TH1D *hDebugClusterDiff5     = nullptr;
    TH1D *hDebugClusterDiff10    = nullptr;
    TH1D *hDebugClusterDiff100   = nullptr;
    TH1D *hDebugClusterDiff1000  = nullptr;
    // truth cluster histograms
    TH1D *hTruClustEne           = nullptr;
    TH1D *hTruClustPosZ          = nullptr;
    TH1I *hTruClustNumHit        = nullptr;
    TH1D *hTruClustParDiff       = nullptr;
    TH2D *hTruClustPosYvsX       = nullptr;
    TH2D *hTruClustVsParEne      = nullptr;
    // truth cluster debug histograms
    TH1D *hDebugTruClustSum5     = nullptr;
    TH1D *hDebugTruClustSum10    = nullptr;
    TH1D *hDebugTruClustSum100   = nullptr;
    TH1D *hDebugTruClustSum1000  = nullptr;
    TH1D *hDebugTruClustDiff5    = nullptr;
    TH1D *hDebugTruClustDiff10   = nullptr;
    TH1D *hDebugTruClustDiff100  = nullptr;
    TH1D *hDebugTruClustDiff1000 = nullptr;
    // event-wise histograms
    TH1I *hEvtNumPar             = nullptr;
    TH1I *hEvtNumHit             = nullptr;
    TH1I *hEvtNumClust           = nullptr;
    TH1I *hEvtNumTruClust        = nullptr;
    TH1D *hEvtSumHitEne          = nullptr;
    TH1D *hEvtSumClustEne        = nullptr;
    TH1D *hEvtSumTruClustEne     = nullptr;
    TH1D *hEvtLeadClustEne       = nullptr;
    TH1D *hEvtLeadTruClustEne    = nullptr;
    TH1D *hEvtSumHitDiff         = nullptr;
    TH1D *hEvtSumClustDiff       = nullptr;
    TH1D *hEvtSumTruClustDiff    = nullptr;
    TH1D *hEvtLeadClustDiff      = nullptr;
    TH1D *hEvtLeadTruClustDiff   = nullptr;
    TH2I *hEvtNumClustVsHit      = nullptr;
    TH2I *hEvtNumTruClustVsClust = nullptr;
    TH2D *hEvtSumHitVsPar        = nullptr;
    TH2D *hEvtSumClustVsPar      = nullptr;
    TH2D *hEvtSumTruClustVsPar   = nullptr;
    TH2D *hEvtLeadClustVsPar     = nullptr;
    TH2D *hEvtLeadTruClustVsPar  = nullptr;

  public:

    // ctor
    JCalibrateHCalProcessor() { SetTypeName(NAME_OF_THIS); }

    // inherited methods
    void InitWithGlobalRootLock()                                      override;
    void ProcessSequential(const std::shared_ptr<const JEvent>& event) override;
    void FinishWithGlobalRootLock()                                    override;
};

// end ------------------------------------------------------------------------
