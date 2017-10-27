#ifndef Ntuplizer_H
#define Ntuplizer_H

// CMSSW
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <DataFormats/HepMCCandidate/interface/GenParticleFwd.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include <FWCore/Framework/interface/ESHandle.h>
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
//
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
// ROOT
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

// C++
#include<memory>
#include<vector>

using namespace std;

class Ntuplizer : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
public:
    explicit Ntuplizer(const edm::ParameterSet&);
    ~Ntuplizer();

    typedef math::XYZTLorentzVector LorentzVector ;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    void Init();
    void FillEvent(const edm::Event&, const edm::EventSetup&);
    void FillElectrons(const edm::Event&, const edm::EventSetup&);
    void FillVertices(const edm::Event&, const edm::EventSetup&);
    void FillTruth(const edm::Event&, const edm::EventSetup&);
    void FillMET (const edm::Event& iEvent, const edm::EventSetup& iSetup);

    //void setMomentum(TLorentzVector & myvector, const LorentzVector & mom);
    void setMomentum(TLorentzVector & myvector, const LorentzVector & mom) ;
    virtual void beginRun(edm::Run const &iEvent, edm::EventSetup const &) override;
    virtual void endRun(edm::Run const &iEvent, edm::EventSetup const &) override;

private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------
    //inputTag
    edm::InputTag EleTag_;
    //edm::InputTag MuonTag_;
    //edm::InputTag JetTag_;
    //edm::InputTag PhotonTag_;
    edm::InputTag VerticesTag_;
    edm::InputTag genPartInputTag_;
    // Trigger Stuff
    //edm::InputTag triggerEventTag_;
    //edm::InputTag MCTag_ ;
    bool isMC_;
    //int lepton_setup;

    //edm::InputTag MuRhoCorrection_;
    //edm::InputTag EleRhoCorrection_;
    //edm::InputTag SigmaRhoCorrection_;
    edm::InputTag PileupSrc_;

    edm::EDGetTokenT<std::vector<reco::Vertex>> vertices_;
    edm::EDGetTokenT<std::vector<reco::GsfElectron>> electrons_;
    edm::EDGetTokenT<reco::ConversionCollection> conversions_;
    edm::EDGetTokenT<reco::BeamSpot> beamspot_;
    edm::EDGetTokenT<reco::PFCandidateCollection> pfcandidates_;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genparticles_;
    //std::vector<edm::InputTag > HLT_Filters_;
    //edm::InputTag SCTag_;

    //tree
    TTree *_mytree;
    TLorentzVector myvector ;

    //global variables
    int _nEvent, _nRun, _nLumi;
    //pile-up
    int _PU_N;

    //vertices
    int _vtx_N;

    //trigger fired names
    char trig_fired_names[10000];

    // METr
    double _met_pf_et,_met_pf_px, _met_pf_py, _met_pf_phi, _met_pf_set, _met_pf_sig;

    //electrons
    int ele_N;
//    TClonesArray * m_electrons;

    std::vector<int> ele_echarge;
    std::vector<double>  ele_he;
    std::vector<double> ele_hebc;
    std::vector<double> ele_eta;
    std::vector<double> ele_phi;
    std::vector<double> ele_pt;
    std::vector<double> ele_e;

    //
    std::vector<double> ele_eseedpout;
    std::vector<double> ele_ep;
    std::vector<double> ele_eseedp;
    std::vector<double> ele_eelepout;
    std::vector<double>  ele_deltaetaseed;
    std::vector<double> ele_deltaetaele;
    std::vector<double> ele_deltaphiseed;
    std::vector<double> ele_deltaphiele;
    std::vector<double> ele_deltaetain;
    std::vector<double> ele_deltaphiin;
    //
    std::vector<double>  ele_sigmaietaieta;
    std::vector<double> ele_sigmaetaeta;
    std::vector<double> ele_sigmaiphiiphi;
    std::vector<double> ele_e15;
    std::vector<double> ele_e25max;
    std::vector<double> ele_e55;
    std::vector<double> ele_e1;
    std::vector<double> ele_r9;
    //
    std::vector<double>  ele_oldsigmaetaeta;
    std::vector<double> ele_oldsigmaietaieta;
    std::vector<double> ele_oldsigmaiphiiphi;
    std::vector<double>ele_oldsigmaietaiphi;
    std::vector<double>ele_oldr9;
    std::vector<double>ele_olde15;
    std::vector<double>ele_olde25max;
    std::vector<double>ele_olde55;
    std::vector<double> ele_oldhe;
    std::vector<double> ele_oldhebc;

    //, ele_e33; ele_e2overe9;
    std::vector<double>  ele_pin_mode;
    std::vector<double> ele_pout_mode;
    std::vector<double> ele_pTin_mode;
    std::vector<double> ele_pTout_mode;
    //
    std::vector<double>  ele_fbrem;
    std::vector<double> ele_SCfbrem;
    std::vector<double> ele_pfSCfbrem;
    std::vector<double> ele_trackfbrem;
    std::vector<int>  ele_nbrem;
    //
    std::vector<double>  ele_mva;
    //
    std::vector<int>  ele_isbarrel;
    std::vector<int>  ele_isendcap;
    std::vector<int>  ele_isEBetaGap;
    std::vector<int>  ele_isEBphiGap;
    std::vector<int>  ele_isEEdeeGap;
    std::vector<int>  ele_isEEringGap;
    std::vector<int>  ele_isecalDriven;
    std::vector<int>  ele_istrackerDriven;
    std::vector<int>  ele_eClass;
    //
    std::vector<double>  ele_gsfchi2;
    std::vector<double>  ele_kfchi2;
    std::vector<int>  ele_kfhits;
    std::vector<int>  ele_gsfhits;
    std::vector<int>  ele_valid_hits;
    std::vector<int>  ele_lost_hits;
    //
    std::vector<double>  ele_dxyB;
    std::vector<double> ele_dxy;
    std::vector<double> ele_dzB;
    std::vector<double> ele_dz;
    std::vector<double> ele_dszB;
    std::vector<double> ele_dsz;
    //
    //std::vector<double>  ele_tkSumPt_dr03; ele_ecalRecHitSumEt_dr03; ele_hcalDepth1TowerSumEt_dr03; ele_hcalDepth2TowerSumEt_dr03;
    //ele_tkSumPt_dr04; ele_ecalRecHitSumEt_dr04; ele_hcalDepth1TowerSumEt_dr04; ele_hcalDepth2TowerSumEt_dr04;
    //
    std::vector<double>  ele_conv_dcot;
    std::vector<double>  ele_conv_dist;
    std::vector<double>  ele_conv_radius;
    std::vector<int>  ele_expected_inner_hits;
    std::vector<int>  ele_vtxconv;
    //
    // std::vector<double>  ele_eidVeryLoose; ele_eidLoose; ele_eidMedium; ele_eidTight;
    //std::vector<double>  ele_eidHZZVeryLoose; ele_eidHZZLoose; ele_eidHZZMedium; ele_eidHZZTight; ele_eidHZZSuperTight;
    std::vector<double>  ele_eidMVATrig;
    std::vector<double> ele_eidMVANoTrig;
    //std::vector<double>  ele_HZZisoTk;ele_HZZisoTk5;  ele_HZZisoEcal; ele_HZZisoHcal;ele_HZZisoComb;
    //
    std::vector<double>  ele_pfChargedHadIso;
    std::vector<double> ele_pfNeutralHadIso;
    std::vector<double> ele_pfPhotonIso;
    std::vector<double> ele_pfChargedIso;
    std::vector<double> ele_pfSumPUIso;
    //ele_pfChargedHadPUIso;ele_pfCombRelIso;
    //
    std::vector<double>  ele_dzPV;
    std::vector<double> ele_d0;
    std::vector<double> ele_d0err;
    std::vector<double>  ele_IP;
    std::vector<double> ele_IPError;
    std::vector<double> ele_SIP;
    //
    std::vector<double>  ele_sclE;
    std::vector<double> ele_sclEt;
    std::vector<double> ele_sclEta;
    std::vector<double> ele_sclPhi;
    std::vector<double>ele_sclRawE;
    std::vector<int>  ele_sclNclus;

    double ele_sclsubE[50][20];
    double ele_sclsubEta[50][20], ele_sclsubPhi[50][20];
    int ele_sclsubisseed[50][20];

    //std::vector<std::vector<double> > ele_sclsubE;//[50][20],

    std::vector<double>  ele_ecalE;
    std::vector<double> ele_ecalErr;
    std::vector<double> ele_trackErr;
    std::vector<double> ele_combErr;
    std::vector<double> ele_PFcombErr;
    std::vector<double>  ele_ecalRegressionEnergy;
    std::vector<double> ele_ecalRegressionError;
    std::vector<double> ele_ecalTrackRegressionEnergy;
    std::vector<double> ele_ecalTrackRegressionError;
    std::vector<double> ele_ecalScale;
    std::vector<double> ele_ecalSmear;
    std::vector<double> ele_ecalRegressionScale;
    std::vector<double> ele_ecalRegressionSmear;
    std::vector<double> ele_ecalTrackRegressionScale;
    std::vector<double> ele_ecalTrackRegressionSmear;

    std::vector<double>  ele_HCALFullConeSum;
    //
    std::vector<double>  ele_sclphiwidth;
    std::vector<double>  ele_scletawidth;
    //
    std::vector<double>  ele_psE;



    std::vector<double>  ele_mvaphys14;
    std::vector<double>  ele_mvaphys14fix;

    //
    // PF Candidates around electrons variables
    std::vector<int> ele_pf_number;


    int  ele_pf_id[50][100];
    double ele_pf_E[50][100];
    double  ele_pf_eta[50][100];
    double  ele_pf_phi[50][100];
    double  ele_pf_pt[50][100];
    double  ele_pf_dz[50][100];
    double  ele_pf_dxy[50][100];
    double  ele_pf_vx[50][100];
    double  ele_pf_vy[50][100];
    double  ele_pf_vz[50][100];
    double  ele_pf_mva_nog[50][100];
    double  ele_pf_mva_epi[50][100];

    /* std::vector<int> ele_pf_id; */
    /* std::vector<float> ele_pf_eta; */
    /* std::vector<float> ele_pf_phi;  */
    /* std::vector<float> ele_pf_pt;  */
    /* std::vector<float> ele_pf_dz;  */
    /* std::vector<float> ele_pf_dxy; */
    /* std::vector<float> ele_pf_vx;  */
    /* std::vector<float> ele_pf_vy;  */
    /* std::vector<float> ele_pf_vz; */
    /* std::vector<float> ele_pf_mva_nog;  */
    /* std::vector<float> ele_pf_mva_epi; */
    /* std::vector<double>  ele_mvafbrem; ele_mvadetain; ele_mvadphiin; ele_mvasieie; ele_mvahoe; ele_mvaeop;  */
    /* 	ele_mvae1x5e5x5; ele_mvaeleopout; ele_mvakfchi2; ele_mvadist;  ele_mvadcot; ele_mvaeta; */
    /* 	ele_mvapt; */
    /*       int ele_mvakfhits; ele_mvamishits; ele_mvaecalseed	; */

    // MC
    //MC
    ////////////////////
    // reco::GenParticles
    //
    std::vector<float> gen_eta_;
    std::vector<float> gen_phi_;
    std::vector<float> gen_pt_;
    std::vector<float> gen_energy_;
    std::vector<int> gen_charge_;
    std::vector<int> gen_pdgid_;
    std::vector<int> gen_status_;
    std::vector<std::vector<int>> gen_daughters_;



};
#endif
