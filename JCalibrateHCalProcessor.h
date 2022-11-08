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
static const float  PParMin(2.);
static const float  PParMax(5.);



class JCalibrateHCalProcessor : public JEventProcessorSequentialRoot {

  private:

    // Data objects we will need from JANA
    PrefetchT<edm4eic::ReconstructedParticle> genParticles  = {this, "GeneratedParticles"};
    PrefetchT<edm4eic::CalorimeterHit>        bhcalRecHits  = {this, "HcalBarrelRecHits"};
    PrefetchT<edm4eic::Cluster>               bhcalClusters = {this, "HcalBarrelClusters"};

    // particle histograms
    TH1D *hParChrg           = nullptr;
    TH1D *hParMass           = nullptr;
    TH1D *hParEne            = nullptr;
    TH1D *hParMom            = nullptr;
    TH1D *hParMomX           = nullptr;
    TH1D *hParMomY           = nullptr;
    TH1D *hParMomZ           = nullptr;
    // reconstructed hit histograms
    TH1D *hRecHitEne         = nullptr;
    TH1D *hRecHitPosZ        = nullptr;
    TH1D *hRecHitParDiff     = nullptr;
    TH2D *hRecHitPosYvsX     = nullptr;
    TH2D *hRecHitVsParEne    = nullptr;
    // cluster histograms
    TH1D *hClusterEne        = nullptr;
    TH1D *hClusterPosZ       = nullptr;
    TH1I *hClusterNumHit     = nullptr;
    TH1D *hClusterParDiff    = nullptr;
    TH2D *hClusterPosYvsX    = nullptr;
    TH2D *hClusterVsParEne   = nullptr;
    // event-wise histograms
    TH1I *hEvtNumPar         = nullptr;
    TH1I *hEvtNumHit         = nullptr;
    TH1I *hEvtNumClust       = nullptr;
    TH1D *hEvtSumHitEne      = nullptr;
    TH1D *hEvtSumClustEne    = nullptr;
    TH1D *hEvtLeadClustEne   = nullptr;
    TH1D *hEvtSumHitDiff     = nullptr;
    TH1D *hEvtSumClustDiff   = nullptr;
    TH1D *hEvtLeadClustDiff  = nullptr;
    TH2D *hEvtSumHitVsPar    = nullptr;
    TH2D *hEvtSumClustVsPar  = nullptr;
    TH2D *hEvtLeadClustVsPar = nullptr;

  public:

    // ctor
    JCalibrateHCalProcessor() { SetTypeName(NAME_OF_THIS); }

    // inherited methods
    void InitWithGlobalRootLock()                                      override;
    void ProcessSequential(const std::shared_ptr<const JEvent>& event) override;
    void FinishWithGlobalRootLock()                                    override;
};

// end ------------------------------------------------------------------------
