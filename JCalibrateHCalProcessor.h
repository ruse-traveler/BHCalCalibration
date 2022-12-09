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
static const float  PParMin(9.9);
static const float  PParMax(10.1);



class JCalibrateHCalProcessor : public JEventProcessorSequentialRoot {

  private:

    // Data objects we will need from JANA
    PrefetchT<edm4eic::ReconstructedParticle> genParticles       = {this, "GeneratedParticles"};
    PrefetchT<edm4eic::CalorimeterHit>        bhcalRecHits       = {this, "HcalBarrelRecHits"};
    PrefetchT<edm4eic::CalorimeterHit>        becalRecHits       = {this, "EcalBarrelSciGlassRecHits"};
    PrefetchT<edm4eic::Cluster>               bhcalClusters      = {this, "HcalBarrelClusters"};
    PrefetchT<edm4eic::Cluster>               becalClusters      = {this, "EcalBarrelSciGlassClusters"};
    PrefetchT<edm4eic::Cluster>               bhcalTruthClusters = {this, "HcalBarrelTruthClusters"};
    PrefetchT<edm4eic::Cluster>               becalTruthClusters = {this, "EcalBarrelSciGlassTruthClusters"};

    // particle histograms
    TH1D *hParChrg                   = nullptr;
    TH1D *hParMass                   = nullptr;
    TH1D *hParEta                    = nullptr;
    TH1D *hParPhi                    = nullptr;
    TH1D *hParEne                    = nullptr;
    TH1D *hParMom                    = nullptr;
    TH1D *hParMomX                   = nullptr;
    TH1D *hParMomY                   = nullptr;
    TH1D *hParMomZ                   = nullptr;
    TH2D *hParEtaVsPhi               = nullptr;
    // hcal reconstructed hit histograms
    TH1D *hHCalRecHitEta             = nullptr;
    TH1D *hHCalRecHitPhi             = nullptr;
    TH1D *hHCalRecHitEne             = nullptr;
    TH1D *hHCalRecHitPosZ            = nullptr;
    TH1D *hHCalRecHitParDiff         = nullptr;
    TH2D *hHCalRecHitPosYvsX         = nullptr;
    TH2D *hHCalRecHitEtaVsPhi        = nullptr;
    TH2D *hHCalRecHitVsParEne        = nullptr;
    // ecal reconstructed hit histograms
    TH1D *hECalRecHitEta             = nullptr;
    TH1D *hECalRecHitPhi             = nullptr;
    TH1D *hECalRecHitEne             = nullptr;
    TH1D *hECalRecHitPosZ            = nullptr;
    TH1D *hECalRecHitParDiff         = nullptr;
    TH2D *hECalRecHitPosYvsX         = nullptr;
    TH2D *hECalRecHitEtaVsPhi        = nullptr;
    TH2D *hECalRecHitVsParEne        = nullptr;
    // hcal reconstructed cluster histograms
    TH1D *hHCalClustEta              = nullptr;
    TH1D *hHCalClustPhi              = nullptr;
    TH1D *hHCalClustEne              = nullptr;
    TH1D *hHCalClustPosZ             = nullptr;
    TH1I *hHCalClustNumHit           = nullptr;
    TH1D *hHCalClustParDiff          = nullptr;
    TH2D *hHCalClustPosYvsX          = nullptr;
    TH2D *hHCalClustEtaVsPhi         = nullptr;
    TH2D *hHCalClustVsParEne         = nullptr;
    // ecal reconstructed cluster histograms
    TH1D *hECalClustEta              = nullptr;
    TH1D *hECalClustPhi              = nullptr;
    TH1D *hECalClustEne              = nullptr;
    TH1D *hECalClustPosZ             = nullptr;
    TH1I *hECalClustNumHit           = nullptr;
    TH1D *hECalClustParDiff          = nullptr;
    TH2D *hECalClustPosYvsX          = nullptr;
    TH2D *hECalClustEtaVsPhi         = nullptr;
    TH2D *hECalClustVsParEne         = nullptr;
    // hcal reco. cluster debug histograms
    TH1D *hHCalDebugClustSum5        = nullptr;
    TH1D *hHCalDebugClustSum10       = nullptr;
    TH1D *hHCalDebugClustSum100      = nullptr;
    TH1D *hHCalDebugClustSum1000     = nullptr;
    TH1D *hHCalDebugClustDiff5       = nullptr;
    TH1D *hHCalDebugClustDiff10      = nullptr;
    TH1D *hHCalDebugClustDiff100     = nullptr;
    TH1D *hHCalDebugClustDiff1000    = nullptr;
    // ecal reco. cluster debug histograms
    TH1D *hECalDebugClustSum5        = nullptr;
    TH1D *hECalDebugClustSum10       = nullptr;
    TH1D *hECalDebugClustSum100      = nullptr;
    TH1D *hECalDebugClustSum1000     = nullptr;
    TH1D *hECalDebugClustDiff5       = nullptr;
    TH1D *hECalDebugClustDiff10      = nullptr;
    TH1D *hECalDebugClustDiff100     = nullptr;
    TH1D *hECalDebugClustDiff1000    = nullptr;
    // hcal truth cluster histograms
    TH1D *hHCalTruClustEta           = nullptr;
    TH1D *hHCalTruClustPhi           = nullptr;
    TH1D *hHCalTruClustEne           = nullptr;
    TH1D *hHCalTruClustPosZ          = nullptr;
    TH1I *hHCalTruClustNumHit        = nullptr;
    TH1D *hHCalTruClustParDiff       = nullptr;
    TH2D *hHCalTruClustPosYvsX       = nullptr;
    TH2D *hHCalTruClustEtaVsPhi      = nullptr;
    TH2D *hHCalTruClustVsParEne      = nullptr;
    // ecal truth cluster histograms
    TH1D *hECalTruClustEta           = nullptr;
    TH1D *hECalTruClustPhi           = nullptr;
    TH1D *hECalTruClustEne           = nullptr;
    TH1D *hECalTruClustPosZ          = nullptr;
    TH1I *hECalTruClustNumHit        = nullptr;
    TH1D *hECalTruClustParDiff       = nullptr;
    TH2D *hECalTruClustPosYvsX       = nullptr;
    TH2D *hECalTruClustEtaVsPhi      = nullptr;
    TH2D *hECalTruClustVsParEne      = nullptr;
    // hcal truth cluster debug histograms
    TH1D *hHCalDebugTruClustSum5     = nullptr;
    TH1D *hHCalDebugTruClustSum10    = nullptr;
    TH1D *hHCalDebugTruClustSum100   = nullptr;
    TH1D *hHCalDebugTruClustSum1000  = nullptr;
    TH1D *hHCalDebugTruClustDiff5    = nullptr;
    TH1D *hHCalDebugTruClustDiff10   = nullptr;
    TH1D *hHCalDebugTruClustDiff100  = nullptr;
    TH1D *hHCalDebugTruClustDiff1000 = nullptr;
    // ecal truth cluster debug histograms
    TH1D *hECalDebugTruClustSum5     = nullptr;
    TH1D *hECalDebugTruClustSum10    = nullptr;
    TH1D *hECalDebugTruClustSum100   = nullptr;
    TH1D *hECalDebugTruClustSum1000  = nullptr;
    TH1D *hECalDebugTruClustDiff5    = nullptr;
    TH1D *hECalDebugTruClustDiff10   = nullptr;
    TH1D *hECalDebugTruClustDiff100  = nullptr;
    TH1D *hECalDebugTruClustDiff1000 = nullptr;
    // hcal event-wise histograms
    TH1I *hEvtHCalNumPar             = nullptr;
    TH1I *hEvtHCalNumHit             = nullptr;
    TH1I *hEvtHCalNumClust           = nullptr;
    TH1I *hEvtHCalNumTruClust        = nullptr;
    TH1D *hEvtHCalSumHitEne          = nullptr;
    TH1D *hEvtHCalSumClustEne        = nullptr;
    TH1D *hEvtHCalSumTruClustEne     = nullptr;
    TH1D *hEvtHCalLeadClustEne       = nullptr;
    TH1D *hEvtHCalLeadTruClustEne    = nullptr;
    TH1D *hEvtHCalSumHitDiff         = nullptr;
    TH1D *hEvtHCalSumClustDiff       = nullptr;
    TH1D *hEvtHCalSumTruClustDiff    = nullptr;
    TH1D *hEvtHCalLeadClustDiff      = nullptr;
    TH1D *hEvtHCalLeadTruClustDiff   = nullptr;
    TH2I *hEvtHCalNumClustVsHit      = nullptr;
    TH2I *hEvtHCalNumTruClustVsClust = nullptr;
    TH2D *hEvtHCalSumHitVsPar        = nullptr;
    TH2D *hEvtHCalSumClustVsPar      = nullptr;
    TH2D *hEvtHCalSumTruClustVsPar   = nullptr;
    TH2D *hEvtHCalLeadClustVsPar     = nullptr;
    TH2D *hEvtHCalLeadTruClustVsPar  = nullptr;
    // ecal event-wise histograms
    TH1I *hEvtECalNumPar             = nullptr;
    TH1I *hEvtECalNumHit             = nullptr;
    TH1I *hEvtECalNumClust           = nullptr;
    TH1I *hEvtECalNumTruClust        = nullptr;
    TH1D *hEvtECalSumHitEne          = nullptr;
    TH1D *hEvtECalSumClustEne        = nullptr;
    TH1D *hEvtECalSumTruClustEne     = nullptr;
    TH1D *hEvtECalLeadClustEne       = nullptr;
    TH1D *hEvtECalLeadTruClustEne    = nullptr;
    TH1D *hEvtECalSumHitDiff         = nullptr;
    TH1D *hEvtECalSumClustDiff       = nullptr;
    TH1D *hEvtECalSumTruClustDiff    = nullptr;
    TH1D *hEvtECalLeadClustDiff      = nullptr;
    TH1D *hEvtECalLeadTruClustDiff   = nullptr;
    TH2I *hEvtECalNumClustVsHit      = nullptr;
    TH2I *hEvtECalNumTruClustVsClust = nullptr;
    TH2D *hEvtECalSumHitVsPar        = nullptr;
    TH2D *hEvtECalSumClustVsPar      = nullptr;
    TH2D *hEvtECalSumTruClustVsPar   = nullptr;
    TH2D *hEvtECalLeadClustVsPar     = nullptr;
    TH2D *hEvtECalLeadTruClustVsPar  = nullptr;

  public:

    // ctor
    JCalibrateHCalProcessor() { SetTypeName(NAME_OF_THIS); }

    // inherited methods
    void InitWithGlobalRootLock()                                      override;
    void ProcessSequential(const std::shared_ptr<const JEvent>& event) override;
    void FinishWithGlobalRootLock()                                    override;
};

// end ------------------------------------------------------------------------
