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
#include <JANA/JEvent.h>
// EDM includes
#include <edm4eic/CalorimeterHit.h>
#include <edm4eic/ReconstructedParticle.h>
#include <edm4eic/ProtoCluster.h>
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
    PrefetchT<edm4eic::Cluster>               bhcalClusters      = {this, "HcalBarrelClusters"};
    PrefetchT<edm4eic::Cluster>               bhcalTruthClusters = {this, "HcalBarrelTruthClusters"};

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
    // hcal cluster hit histograms
    TH1D *hHCalClustHitEta           = nullptr;
    TH1D *hHCalClustHitPhi           = nullptr;
    TH1D *hHCalClustHitEne           = nullptr;
    TH1D *hHCalClustHitPosZ          = nullptr;
    TH1D *hHCalClustHitParDiff       = nullptr;
    TH2D *hHCalClustHitPosYvsX       = nullptr;
    TH2D *hHCalClustHitEtaVsPhi      = nullptr;
    TH2D *hHCalClustHitVsParEne      = nullptr;
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
    // hcal truth cluster hit histograms
    TH1D *hHCalTruClustHitEta        = nullptr;
    TH1D *hHCalTruClustHitPhi        = nullptr;
    TH1D *hHCalTruClustHitEne        = nullptr;
    TH1D *hHCalTruClustHitPosZ       = nullptr;
    TH1D *hHCalTruClustHitParDiff    = nullptr;
    TH2D *hHCalTruClustHitPosYvsX    = nullptr;
    TH2D *hHCalTruClustHitEtaVsPhi   = nullptr;
    TH2D *hHCalTruClustHitVsParEne   = nullptr;
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

  public:

    // ctor
    JCalibrateHCalProcessor() { SetTypeName(NAME_OF_THIS); }

    // inherited methods
    void InitWithGlobalRootLock()                                      override;
    void ProcessSequential(const std::shared_ptr<const JEvent>& event) override;
    void FinishWithGlobalRootLock()                                    override;
};

// end ------------------------------------------------------------------------
