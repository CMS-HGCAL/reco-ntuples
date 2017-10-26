// -*- C++ -*-
//
// Package:    Analyzer/Ntuplizer
// Class:      Ntuplizer
//
/**\class Ntuplizer Ntuplizer.cc Analyzer/Ntuplizer/plugins/Ntuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Christophe ochando
//         Created:  Mon, 10 Mar 2014 14:51:20 GMT
//
//


// MY include
#include "Ntuplizer.h"

// C++ include files
#include <memory>

// CMSSW include
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include <DataFormats/Common/interface/MergeableCounter.h>
//#include <DataFormats/Common/interface/View.h>
//#include <DataFormats/Candidate/interface/Candidate.h>
//#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
//#include <DataFormats/MuonReco/interface/Muon.h>
//#include <DataFormats/MuonReco/interface/MuonFwd.h>
//


//
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
//CMSSW/ RecoEcal/ EgammaCoreTools/ interface/ EcalClusterLazyTools.h
#include "FWCore/Utilities/interface/isFinite.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "DataFormats/PatCandidates/interface/MET.h"
//#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "DataFormats/PatCandidates/interface/Jet.h"
//#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>

#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"

//"

//
//#include "DataFormats/Math/interface/LorentzVector.h"
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"


#include "DataFormats/TrackReco/interface/Track.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

// Trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

//
// class declaration
//

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//using namespace std;
using namespace reco;
using namespace edm;

// =============================================================================================
// constructor
Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig) :
// ==============================================================================================
EleTag_ (iConfig.getParameter<edm::InputTag> ("EleTag")),
VerticesTag_(iConfig.getParameter<edm::InputTag> ("VerticesTag")),
genPartInputTag_(iConfig.getParameter<edm::InputTag> ("GenParticles")),
isMC_ (iConfig.getParameter<bool>("isMC")),
PileupSrc_ ("addPileupInfo")
// std::vector<edm::InputTag> MVAidCollection_;
{
    vertices_ = consumes<std::vector<reco::Vertex>>(VerticesTag_);
    electrons_ = consumes<std::vector<reco::GsfElectron>>(EleTag_);
    conversions_ = consumes<reco::ConversionCollection>(edm::InputTag("allConversions"));
    beamspot_ = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
    pfcandidates_ = consumes<reco::PFCandidateCollection>(edm::InputTag("particleFlow"));
    genparticles_ = consumes<std::vector<reco::GenParticle>>(genPartInputTag_);
}

// =============================================================================================
// destructor
Ntuplizer::~Ntuplizer()
// =============================================================================================
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
//  delete m_electrons ;
//   if (fill_L1trigger) {
//     delete m_L1emIso;
//     delete m_L1emNonIso;
//   }
//   delete m_muons;
//   delete _m_jets_pf;
//   delete m_photons;

  if(isMC_ ) {

  } // if MC


}

// =============================================================================================
// ------------ method called once each job just before starting event loop  ------------
void Ntuplizer::beginJob()
//=============================================================================================
{
  // Book histograms
  usesResource(TFileService::kSharedResource);
  edm::Service<TFileService> fs ;
  _mytree  = fs->make <TTree>("simpleRootTree","simpleRootTree");

  //// Counters
  //_mytree->Branch("Nevt_Gen",&Nevt_Gen,"Nevt_Gen/I");
  //_mytree->Branch("Nevt_Skim",&Nevt_afterSkim,"Nevt_Skim/I");

  // Global
  _mytree->Branch("nEvent",&_nEvent,"nEvent/I");
  _mytree->Branch("nRun",&_nRun,"nRun/I");
  _mytree->Branch("nLumi",&_nLumi,"nLumi/I");

  // Pile UP
  _mytree->Branch("PU_N",&_PU_N,"PU_N/I");

  // Trigger
  _mytree->Branch("trig_fired_names",&trig_fired_names,"trig_fired_names[10000]/C");

  // Vertices
  _mytree->Branch("vtx_N",&_vtx_N,"vtx_N/I");

  // Electrons
  _mytree->Branch("ele_N",&ele_N,"ele_N/I");
//  m_electrons = new TClonesArray ("TLorentzVector");
//  _mytree->Branch ("electrons", "TClonesArray", &m_electrons, 256000,0);
  _mytree->Branch("ele_eta,&ele_eta");
  _mytree->Branch("ele_phi,&ele_eta");
  _mytree->Branch("ele_pt,&ele_eta");

  _mytree->Branch("ele_echarge",&ele_echarge); // ,"ele_echarge[50]/I");
  //
  _mytree->Branch("ele_he",&ele_he); //,"ele_he[50]/D");
  _mytree->Branch("ele_hebc",&ele_hebc); //,"ele_hebc[50]/D");
  _mytree->Branch("ele_oldhe",&ele_oldhe); //,"ele_oldhe[50]/D");
  _mytree->Branch("ele_oldhebc",&ele_oldhebc); //,"ele_oldhebc[50]/D");
  //
  _mytree->Branch("ele_eseedpout",&ele_eseedpout); //,"ele_eseedpout[50]/D");
  _mytree->Branch("ele_ep",&ele_ep); //,"ele_ep[50]/D");
  _mytree->Branch("ele_eseedp",&ele_eseedp); //,"ele_eseedp[50]/D");
  _mytree->Branch("ele_eelepout",&ele_eelepout); //,"ele_eelepout[50]/D");
  //
  _mytree->Branch("ele_pin_mode",&ele_pin_mode); //,"ele_pin_mode[50]/D");
  _mytree->Branch("ele_pout_mode",&ele_pout_mode); //,"ele_pout_mode[50]/D");
  _mytree->Branch("ele_pTin_mode",&ele_pTin_mode); //,"ele_pTin_mode[50]/D");
  _mytree->Branch("ele_pTout_mode",&ele_pTout_mode); // ,"ele_pTout_mode[50]/D");
  //
  _mytree->Branch("ele_deltaetaseed",&ele_deltaetaseed); //,"ele_deltaetaseed[50]/D");
  _mytree->Branch("ele_deltaphiseed",&ele_deltaphiseed); //,"ele_deltaphiseed[50]/D");
  _mytree->Branch("ele_deltaetaele",&ele_deltaetaele); //,"ele_deltaetaele[50]/D");
  _mytree->Branch("ele_deltaphiele",&ele_deltaphiele); //,"ele_deltaphiele[50]/D");
  _mytree->Branch("ele_deltaetain",&ele_deltaetain); //,"ele_deltaetain[50]/D");
  _mytree->Branch("ele_deltaphiin",&ele_deltaphiin); //,"ele_deltaphiin[50]/D");
  //
  _mytree->Branch("ele_sigmaietaieta",&ele_sigmaietaieta); //,"ele_sigmaietaieta[50]/D");
  _mytree->Branch("ele_sigmaetaeta",&ele_sigmaetaeta); //,"ele_sigmaetaeta[50]/D");
  _mytree->Branch("ele_sigmaiphiiphi",&ele_sigmaiphiiphi); //,"ele_sigmaiphiiphi[50]/D");
  _mytree->Branch("ele_e15",&ele_e15); //,"ele_e15[50]/D");
  _mytree->Branch("ele_e25max",&ele_e25max); //,"ele_e25max[50]/D");
  _mytree->Branch("ele_e55",&ele_e55); //,"ele_e55[50]/D");
  //_mytree->Branch("ele_e1",&ele_e1,"ele_e1[50]/D");
  _mytree->Branch("ele_r9",&ele_r9); //,"ele_r9[50]/D");
  //
  _mytree->Branch("ele_oldsigmaietaieta",&ele_oldsigmaietaieta); //,"ele_oldsigmaietaieta[50]/D");
  _mytree->Branch("ele_oldsigmaetaeta",&ele_oldsigmaetaeta); //,"ele_oldsigmaetaeta[50]/D");
  _mytree->Branch("ele_oldsigmaiphiiphi",&ele_oldsigmaiphiiphi); //,"ele_oldsigmaiphiiphi[50]/D");
  _mytree->Branch("ele_oldsigmaietaiphi",&ele_oldsigmaietaiphi); //,"ele_oldsigmaietaiphi[50]/D");
  _mytree->Branch("ele_olde15",&ele_olde15); //,"ele_olde15[50]/D");
  _mytree->Branch("ele_olde25max",&ele_olde25max); //,"ele_olde25max[50]/D");
  _mytree->Branch("ele_olde55",&ele_olde55); //,"ele_olde55[50]/D");
  //_mytree->Branch("ele_olde1",&ele_olde1,"ele_olde1[50]/D");
  _mytree->Branch("ele_oldr9",&ele_oldr9); //,"ele_oldr9[50]/D");
  //
  //_mytree->Branch("ele_e33",&ele_e33,"ele_e33[50]/D");  ---> No ECAL Reduced Collection
  //_mytree->Branch("ele_e2overe9",&ele_e2overe9,"ele_e2overe9[50]/D");  ---> No ECAL Reduced Collection
  //
  _mytree->Branch("ele_fbrem",&ele_fbrem); //,"ele_fbrem[50]/D");
  _mytree->Branch("ele_trackfbrem",&ele_trackfbrem); //,"ele_trackfbrem[50]/D");
  _mytree->Branch("ele_SCfbrem",&ele_SCfbrem); //,"ele_SCfbrem[50]/D");
  _mytree->Branch("ele_pfSCfbrem",&ele_pfSCfbrem); //,"ele_pfSCfbrem[50]/D");
  _mytree->Branch("ele_nbrem",&ele_nbrem); //,"ele_nbrem[50]/I");
  _mytree->Branch("ele_eClass",&ele_eClass); //,"ele_eClass[50]/I");
  //
  _mytree->Branch("ele_mva",&ele_mva); //,"ele_mva[50]/D");
  //
  _mytree->Branch("ele_isbarrel",&ele_isbarrel); //,"ele_isbarrel[50]/I");
  _mytree->Branch("ele_isendcap",&ele_isendcap); //,"ele_isendcap[50]/I");
  _mytree->Branch("ele_isEBetaGap",&ele_isEBetaGap); //,"ele_isEBetaGap[50]/I");
  _mytree->Branch("ele_isEBphiGap",&ele_isEBphiGap); //p,"ele_isEBphiGap[50]/I");
  _mytree->Branch("ele_isEEdeeGap",&ele_isEEdeeGap); //,"ele_isEEdeeGap[50]/I");
  _mytree->Branch("ele_isEEringGap",&ele_isEEringGap); //,"ele_isEEringGap[50]/I");
  _mytree->Branch("ele_isecalDriven",&ele_isecalDriven); //,"ele_isecalDriven[50]/I");
  _mytree->Branch("ele_istrackerDriven",&ele_istrackerDriven); //,"ele_istrackerDriven[50]/I");
  //
  _mytree->Branch("ele_valid_hits",&ele_valid_hits); // _valid_hits[50]/I");
  _mytree->Branch("ele_lost_hits",&ele_lost_hits); // _lost_hits[50]/I");
  _mytree->Branch("ele_gsfchi2",&ele_gsfchi2); // _gsfchi2[50]/D");
  //
  _mytree->Branch("ele_dxy",&ele_dxy); // _dxy[50]/D");
  _mytree->Branch("ele_dxyB",&ele_dxyB); // _dxyB[50]/D");
  _mytree->Branch("ele_dz",&ele_dz); // _dz[50]/D");
  _mytree->Branch("ele_dzB",&ele_dzB); // _dzB[50]/D");
  _mytree->Branch("ele_dsz",&ele_dsz); // _dsz[50]/D");
  _mytree->Branch("ele_dszB",&ele_dszB); // _dszB[50]/D");
  //
  _mytree->Branch("ele_dzPV", &ele_dzPV, "ele_dzPV[50]/D");
  _mytree->Branch("ele_d0", &ele_d0, "ele_d0[50]/D");
  _mytree->Branch("ele_d0err",&ele_d0err, "ele_d0err[50]/D");
  //
  _mytree->Branch("ele_IP",&ele_IP); // _IP[50]/D");
  _mytree->Branch("ele_IPError",&ele_IPError); // _IPError[50]/D");
  _mytree->Branch("ele_SIP",&ele_SIP); // _SIP[50]/D");
  //  _mytree->Branch("ele_tkSumPt_dr03",&ele_tkSumPt_dr03); // _tkSumPt_dr03[50]/D");
  //   _mytree->Branch("ele_ecalRecHitSumEt_dr03",&ele_ecalRecHitSumEt_dr03); // _ecalRecHitSumEt_dr03[50]/D");
  //   _mytree->Branch("ele_hcalDepth1TowerSumEt_dr03",&ele_hcalDepth1TowerSumEt_dr03); // _hcalDepth1TowerSumEt_dr03[50]/D");
  //   _mytree->Branch("ele_hcalDepth2TowerSumEt_dr03",&ele_hcalDepth2TowerSumEt_dr03); // _hcalDepth2TowerSumEt_dr03[50]/D");
  //   _mytree->Branch("ele_tkSumPt_dr04",&ele_tkSumPt_dr04); // _tkSumPt_dr04[50]/D");
  //   _mytree->Branch("ele_ecalRecHitSumEt_dr04",&ele_ecalRecHitSumEt_dr04); // _ecalRecHitSumEt_dr04[50]/D");
  //   _mytree->Branch("ele_hcalDepth1TowerSumEt_dr04",&ele_hcalDepth1TowerSumEt_dr04); // _hcalDepth1TowerSumEt_dr04[50]/D");
  //   _mytree->Branch("ele_hcalDepth2TowerSumEt_dr04",&ele_hcalDepth2TowerSumEt_dr04); // _hcalDepth2TowerSumEt_dr04[50]/D");
  //
  _mytree->Branch("ele_conv_dist",&ele_conv_dist); // _conv_dist[50]/D");
  _mytree->Branch("ele_conv_dcot",&ele_conv_dcot); // _conv_dcot[50]/D");
  _mytree->Branch("ele_conv_radius",&ele_conv_radius); // _conv_radius[50]/D");
  _mytree->Branch("ele_expected_inner_hits",&ele_expected_inner_hits); // _expected_inner_hits[50]/I");
  _mytree->Branch("ele_vtxconv",&ele_vtxconv); // _vtxconv[50]/I");
  //
  _mytree->Branch("ele_pfChargedHadIso",&ele_pfChargedHadIso); // _pfChargedHadIso[50]/D");
  _mytree->Branch("ele_pfNeutralHadIso",&ele_pfNeutralHadIso); // _pfNeutralHadIso[50]/D");
  _mytree->Branch("ele_pfPhotonIso",&ele_pfPhotonIso); // _pfPhotonIso[50]/D");

  _mytree->Branch("ele_pfChargedIso", &ele_pfChargedIso); // "ele_pfChargedIso[50]/D");
  _mytree->Branch("ele_pfSumPUIso", &ele_pfSumPUIso); //"ele_pfSumPUIso[50]/D");
  //  _mytree->Branch("ele_pfChargedHadPUIso",&ele_pfChargedHadPUIso); // _pfChargedHadPUIso[50]/D");
  //   _mytree->Branch("ele_pfCombRelIso",&ele_pfCombRelIso); // _pfCombRelIso[50]/D");
  //

  //
  _mytree->Branch("ele_sclRawE", &ele_sclRawE); //"ele_sclRawE[50]/D");
  _mytree->Branch("ele_sclE",    &ele_sclE); //    "ele_sclE[50]/D");
  _mytree->Branch("ele_sclEt",   &ele_sclEt); //   "ele_sclEt[50]/D");
  _mytree->Branch("ele_sclEta",  &ele_sclEta); //  "ele_sclEta[50]/D");
  _mytree->Branch("ele_sclPhi",  &ele_sclPhi); //  "ele_sclPhi[50]/D");
  _mytree->Branch("ele_sclNclus",  &ele_sclNclus); //  "ele_sclNclus[50]/I");
  _mytree->Branch("ele_sclphiwidth", &ele_sclphiwidth); //"ele_sclphiwidth[50]/D");
  _mytree->Branch("ele_scletawidth", &ele_scletawidth); // "ele_scletawidth[50]/D");
  //
  _mytree->Branch("ele_sclsubE", &ele_sclsubE, "ele_sclsubE[50][20]/D");
  _mytree->Branch("ele_sclsubEta", &ele_sclsubEta, "ele_sclsubEta[50][20]/D");
  _mytree->Branch("ele_sclsubPhi", &ele_sclsubPhi, "ele_sclsubPhi[50][20]/D");
  _mytree->Branch("ele_sclsubisseed", &ele_sclsubisseed, "ele_sclsubisseed[50][20]/I");
  //
  _mytree->Branch("ele_psE", &ele_psE); // "ele_psE[50]/D");
  //
  _mytree->Branch("ele_ecalE",   &ele_ecalE); //  "ele_ecalE[50]/D");
  _mytree->Branch("ele_ecalErr",   &ele_ecalErr); //  "ele_ecalErr[50]/D");
  _mytree->Branch("ele_trackErr",  &ele_trackErr); //  "ele_trackErr[50]/D");
  _mytree->Branch("ele_combErr",   &ele_combErr); //  "ele_combErr[50]/D");
  _mytree->Branch("ele_PFcombErr", &ele_PFcombErr); // "ele_PFcombErr[50]/D");

  // Electron ID
  _mytree->Branch("ele_mvaphys14",   &ele_mvaphys14); //  "ele_mvaphys14[50]/D");
  _mytree->Branch("ele_mvaphys14fix",   &ele_mvaphys14fix); //   "ele_mvaphys14fix[50]/D");


  //  _mytree->Branch("ele_ecalRegressionEnergy",   ele_ecalRegressionEnergy,   "ele_ecalRegressionEnergy[50]/D");
  //   _mytree->Branch("ele_ecalRegressionError", ele_ecalRegressionError, "ele_ecalRegressionError[50]/D");
  //   _mytree->Branch("ele_ecalTrackRegressionEnergy",&ele_ecalTrackRegressionEnergy); // _ecalTrackRegressionEnergy[50]/D");
  //   _mytree->Branch("ele_ecalTrackRegressionError",&ele_ecalTrackRegressionError); // _ecalTrackRegressionError[50]/D");
  //   _mytree->Branch("ele_ecalScale",ele_ecalScale); // _ecalScale[50]/D");
  //   _mytree->Branch("ele_ecalSmear",ele_ecalSmear); // _ecalSmear[50]/D");
  //   _mytree->Branch("ele_ecalRegressionScale",ele_ecalRegressionScale); // _ecalRegressionScale[50]/D");
  //   _mytree->Branch("ele_ecalRegressionSmear",ele_ecalRegressionSmear); // _ecalRegressionSmear[50]/D");
  //   _mytree->Branch("ele_ecalTrackRegressionScale",ele_ecalTrackRegressionScale); // _ecalTrackRegressionScale[50]/D");
  //   _mytree->Branch("ele_ecalTrackRegressionSmear",ele_ecalTrackRegressionSmear); // _ecalTrackRegressionSmear[50]/D");

  // Variables for the mva. Most of them are duplicated, but since they are corrected at analysis level, it could be dangerous
  //  _mytree->Branch("ele_mvafbrem", ele_mvafbrem); // _mvafbrem[50]/D");
  //   _mytree->Branch("ele_mvadetain", ele_mvadetain); // _mvadetain[50]/D");
  //   _mytree->Branch("ele_mvadphiin", ele_mvadphiin); // _mvadphiin[50]/D");
  //   _mytree->Branch("ele_mvasieie", ele_mvasieie); // _mvasiesie[50]/D");
  //   _mytree->Branch("ele_mvahoe", ele_mvahoe); // _mvahoe[50]/D");
  //   _mytree->Branch("ele_mvaeop", ele_mvaeop); // _mvaeop[50]/D");
  //   _mytree->Branch("ele_mvae1x5e5x5", ele_mvae1x5e5x5); // _mvae1x5e5x5[50]/D");
  //   _mytree->Branch("ele_mvaeleopout", ele_mvaeleopout); // _mvaeleopout[50]/D");
  _mytree->Branch("ele_kfchi2",&ele_kfchi2); // _kfchi2[50]/D");
  _mytree->Branch("ele_kfhits",&ele_kfhits); // _kfhits[50]/I");
  _mytree->Branch("ele_gsfhits",&ele_gsfhits); // _gsfhits[50]/I");

  //   _mytree->Branch("ele_mvamishits", ele_mvamishits); // _mvamisthits[50]/I");
  //   _mytree->Branch("ele_mvadist", ele_mvadist); // _mvadist[50]/D");
  //   _mytree->Branch("ele_mvadcot", ele_mvadcot); // _mvadcot[50]/D");
  //   _mytree->Branch("ele_mvaeta", ele_mvaeta); // _mvaeta[50]/D");
  //   _mytree->Branch("ele_mvapt", ele_mvapt); // _mvapt[50]/D");
  //   _mytree->Branch("ele_mvaecalseed", ele_mvaecalseed); // _mvaecalseed[50]/I");

  // PF Candidates around electrons...
  //_mytree->Branch("ele_sclsubE", &ele_sclsubE, "ele_sclsubE[50][20]/D");
  _mytree->Branch("ele_pf_number", &ele_pf_number);
  _mytree->Branch("ele_pf_id", &ele_pf_id,  "ele_pf_id[50][100]/I");
  _mytree->Branch("ele_pf_E", &ele_pf_E,  "ele_pf_E[50][100]/D");
  _mytree->Branch("ele_pf_phi", &ele_pf_phi,  "ele_pf_phi[50][100]/D");
  _mytree->Branch("ele_pf_eta", &ele_pf_eta,  "ele_pf_eta[50][100]/D");
  _mytree->Branch("ele_pf_pt", &ele_pf_pt,  "ele_pf_pt[50][100]/D");
  _mytree->Branch("ele_pf_dz", &ele_pf_dz,  "ele_pf_dz[50][100]/D");
  _mytree->Branch("ele_pf_dxy", &ele_pf_dxy,  "ele_pf_dxy[50][100]/D");
  _mytree->Branch("ele_pf_vx", &ele_pf_vx,  "ele_pf_vx[50][100]/D");
  _mytree->Branch("ele_pf_vy", &ele_pf_vy,  "ele_pf_vy[50][100]/D");
  _mytree->Branch("ele_pf_vz", &ele_pf_vz,  "ele_pf_vz[50][100]/D");
  _mytree->Branch("ele_pf_mva_nog", &ele_pf_mva_nog,  "ele_pf_mva_nog[50][100]/D");
  _mytree->Branch("ele_pf_mva_epi", &ele_pf_mva_epi,  "ele_pf_mva_epi[50][100]/D");

  // PFMET
  //_mytree->Branch("met_pf_et",&_met_pf_et); //,"met_pf_et/D");
  //_mytree->Branch("met_pf_px",&_met_pf_px); //,"met_pf_px/D");
  //_mytree->Branch("met_pf_py",&_met_pf_py); //,"met_pf_py/D");
  //_mytree->Branch("met_pf_phi",&_met_pf_phi); //,"met_pf_phi/D");
  //_mytree->Branch("met_pf_set",&_met_pf_set); //,"met_pf_set/D");
  //_mytree->Branch("met_pf_sig",&_met_pf_sig); //,"met_pf_sig/D");

  // Truth Leptons
  _mytree->Branch("gen_eta", &gen_eta_);
  _mytree->Branch("gen_phi", &gen_phi_);
  _mytree->Branch("gen_pt", &gen_pt_);
  _mytree->Branch("gen_energy", &gen_energy_);
  _mytree->Branch("gen_charge", &gen_charge_);
  _mytree->Branch("gen_pdgid", &gen_pdgid_);
  _mytree->Branch("gen_status", &gen_status_);
  _mytree->Branch("gen_daughters", &gen_daughters_);
  //cout << "truth leptons" << endl;
  //_m_MC_gen_V = new TClonesArray ("TLorentzVector");
  //_mytree->Branch ("MC_gen_V", "TClonesArray", &_m_MC_gen_V, 256000,0);
  //_mytree->Branch ("MC_gen_V_pdgid",&_MC_gen_V_pdgid, "MC_gen_V_pdgid[10]/D");
  ////
  //_m_MC_gen_Higgs = new TClonesArray ("TLorentzVector");
  //_mytree->Branch ("MC_gen_Higgs", "TClonesArray", &_m_MC_gen_Higgs, 256000,0);
  ////_mytree->Branch ("MC_gen_Higgs_pdgid",&_MC_gen_Higgs_pdgid, "MC_gen_Higgs_pdgid[10]/D");
  ////
  //_m_MC_gen_leptons         = new TClonesArray ("TLorentzVector");
  //_m_MC_gen_leptons_status1 = new TClonesArray ("TLorentzVector");
  //_m_MC_gen_leptons_status2 = new TClonesArray ("TLorentzVector");
  //_mytree->Branch ("MC_gen_leptons", "TClonesArray", &_m_MC_gen_leptons, 256000,0);
  //_mytree->Branch ("MC_gen_leptons_status1", "TClonesArray", &_m_MC_gen_leptons_status1, 256000,0);
  //_mytree->Branch ("MC_gen_leptons_status2", "TClonesArray", &_m_MC_gen_leptons_status2, 256000,0);
  //_mytree->Branch ("MC_gen_leptons_pdgid",&_MC_gen_leptons_pdgid, "MC_gen_leptons_pdgid[30]/D");
  //_mytree->Branch ("MC_gen_leptons_status1_pdgid",&_MC_gen_leptons_status1_pdgid, "MC_gen_leptons_status1_pdgid[30]/D");
  //_mytree->Branch ("MC_gen_leptons_status1_FromWZ",&_MC_gen_leptons_status1_FromWZ, "MC_gen_leptons_status1_FromWZ[30]/I");
  //_mytree->Branch ("MC_gen_leptons_status1_FromTaus",&_MC_gen_leptons_status1_FromTaus, "MC_gen_leptons_status1_FromTaus[30]/I");
  //_mytree->Branch ("MC_gen_leptons_status1_FromNonPrompt",&_MC_gen_leptons_status1_FromNonPrompt, "MC_gen_leptons_status1_FromNonPrompt[30]/I");
  //_mytree->Branch ("MC_gen_leptons_status1_FromBC",&_MC_gen_leptons_status1_FromBC, "MC_gen_leptons_status1x_FromBC[30]/I");
  //_mytree->Branch ("MC_gen_leptons_status2_pdgid",&_MC_gen_leptons_status2_pdgid, "MC_gen_leptons_status2_pdgid[30]/D");
  ////_mytree->Branch ("MC_pthat",&_MC_pthat,"MC_pthat/D");
  //_mytree->Branch ("MC_flavor",&_MC_flavor,"MC_flavor[2]/I");
  //_m_MC_gen_photons         = new TClonesArray ("TLorentzVector");
  //_mytree->Branch ("MC_gen_photons", "TClonesArray", &_m_MC_gen_photons, 256000,0);
  //_mytree->Branch ("MC_gen_photons_isFSR",&_MC_gen_photons_isFSR,"MC_gen_photons_isFSR[5000]/I");


}


//
// member functions
//
// =============================================================================================
// ------------ method called for each event  ------------
void Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// =============================================================================================
{
   using namespace edm;

   Init();

   //cout << "salut" << endl;

   FillEvent(iEvent, iSetup);

   FillVertices(iEvent, iSetup);

//   m_electrons -> Clear();
   FillElectrons(iEvent, iSetup);

//   FillMET (iEvent, iSetup);


   if(isMC_ ) {

     FillTruth(iEvent, iSetup);
   }

   _mytree->Fill();

// #ifdef THIS_IS_AN_EVENT_EXAMPLE
//    Handle<ExampleData> pIn;
//    iEvent.getByLabel("example",pIn);
// #endif

// #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//    ESHandle<SetupData> pSetup;
//    iSetup.get<SetupRecord>().get(pSetup);
// #endif
}

// =============================================================================================
void Ntuplizer::FillEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup)
//=============================================================================================
{

  _nEvent = iEvent.id().event();
  _nRun   = iEvent.id().run();
  _nLumi  = iEvent.luminosityBlock();

  // ----------------
  // Fired Triggers
  // ----------------
  //cout << "fired" << endl;
  /*
  Handle<edm::TriggerResults> triggerResultsHandle;
  iEvent.getByLabel (HLTTag_,triggerResultsHandle);
  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResultsHandle);

  //Get List of available Triggers
  //for (int in=0;in<(int)triggerNames.size();in++) {
  //cout << " Trigger Names " << in << " = " << triggerNames.triggerName(in) << endl;
  //} // for loop in triggernames

  // LOOP Over Trigger Results
  char trig_fired_names_local[10000];
  strcpy(trig_fired_names_local,"*");
  for (int iHLT = 0 ;
       iHLT<static_cast<int>(triggerResultsHandle->size());
       ++iHLT) {

    if (triggerResultsHandle->accept (iHLT)) {
      if ( strlen(trig_fired_names_local) <= 9950) {
	{
	  const char* c_str();
	  string hlt_string = triggerNames.triggerName(iHLT);
	  strcat(trig_fired_names_local,hlt_string.c_str());
	  strcat(trig_fired_names_local,"*");
	}
      }
    } // if HLT
  }
  strcpy(trig_fired_names,trig_fired_names_local);
*/
  //cout << "trig = " << trig_fired_names << endl;
  //cout << "" << endl;

} // end of FillEvent


// =============================================================================================
void Ntuplizer::FillVertices(const edm::Event& iEvent, const edm::EventSetup& iSetup)
//=============================================================================================
{

  Handle<vector<reco::Vertex> >  recoPrimaryVertexCollection;
  //   //iEvent.getByLabel("goodPrimaryVertices",recoPrimaryVertexCollection);
  iEvent.getByToken(vertices_, recoPrimaryVertexCollection);

  //const reco::VertexCollection & vertices = *recoPrimaryVertexCollection.product();

//   // 	edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
//   // 	iEvent.getByType(recoBeamSpotHandle);
//   // 	const reco::BeamSpot bs = *recoBeamSpotHandle;

  //int vtx_counter=0;
  _vtx_N = recoPrimaryVertexCollection->size();

} // end of FillVertices

// =============================================================================================
void Ntuplizer::FillElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
//=============================================================================================
{

  edm::Handle<reco::GsfElectronCollection> electronsCol;
  iEvent.getByToken(electrons_, electronsCol);

  Handle<reco::VertexCollection> thePrimaryVertexColl;
  iEvent.getByToken(vertices_ ,thePrimaryVertexColl);

  Vertex dummy;
  const Vertex *pv = &dummy;
  if (thePrimaryVertexColl->size() != 0) {
    pv = &*thePrimaryVertexColl->begin();
  } else { // create a dummy PV
    Vertex::Error e;
    e(0, 0) = 0.0015 * 0.0015;
    e(1, 1) = 0.0015 * 0.0015;
    e(2, 2) = 15. * 15.;
    Vertex::Point p(0, 0, 0);
    dummy = Vertex(p, e, 0, 0, 0);
  }

  edm::ESHandle<TransientTrackBuilder> builder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
  TransientTrackBuilder thebuilder = *(builder.product());

  //load the conversion collection
  edm::Handle<reco::ConversionCollection> conversions_h;
  iEvent.getByToken(conversions_, conversions_h);

  // get the beam spot
  edm::Handle<reco::BeamSpot> beamspot_h;
  iEvent.getByToken(beamspot_, beamspot_h);
  const reco::BeamSpot &beamSpot = *(beamspot_h.product());

  // get the value map for eiD

//  TClonesArray & electrons = *m_electrons;
  int counter = 0;
  ele_N = electronsCol->size();

  //cout << "ele N = " << ele_N << endl;

  edm::Handle< reco::PFCandidateCollection > pfHandle;
  iEvent.getByToken(pfcandidates_, pfHandle);
  const reco::PFCandidateCollection* pfCandColl = pfHandle.product();

  // -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  //  Loop on electrons
  // -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  for (reco::GsfElectronCollection::const_iterator ielectrons=electronsCol->begin(); ielectrons!=electronsCol->end();++ielectrons) {
//    if(counter>49) continue;

//    setMomentum(myvector, ielectrons->p4());
//    new (electrons[counter]) TLorentzVector(myvector);
    ele_eta.push_back(ielectrons->eta());
    ele_phi.push_back(ielectrons->phi());
    ele_pt.push_back(ielectrons->pt());

    ele_echarge.push_back(ielectrons->charge());
    //
    ele_he.push_back(ielectrons->hcalOverEcal()); //hadronicOverEm() ;
    ele_hebc.push_back(ielectrons-> hcalOverEcalBc());

    // TrackCluster Matching
    ele_eseedpout.push_back(ielectrons->eSeedClusterOverPout());
    ele_ep.push_back(ielectrons->eSuperClusterOverP());
    ele_eseedp.push_back(ielectrons->eSeedClusterOverP());
    ele_eelepout.push_back(ielectrons->eEleClusterOverPout());
    //
    ele_pin_mode.push_back(ielectrons->trackMomentumAtVtx().R() );

    //cout << "pin = " <<  ielectrons->trackMomentumAtVtx().R() << " gsf p = " <<  ielectrons->gsfTrack()->p() << " P from E/p = " << 1./(ielectrons->eSuperClusterOverP()/ielectrons->ecalEnergy()) << endl;

    ele_pout_mode.push_back(ielectrons->trackMomentumOut().R());
    ele_pTin_mode.push_back(ielectrons->trackMomentumAtVtx().Rho());
    ele_pTout_mode.push_back(ielectrons->trackMomentumOut().Rho());
    //
    ele_deltaetaseed.push_back(ielectrons->deltaEtaSeedClusterTrackAtCalo());
    ele_deltaphiseed.push_back(ielectrons->deltaPhiSeedClusterTrackAtCalo());
    ele_deltaetaele.push_back(ielectrons->deltaEtaEleClusterTrackAtCalo());
    ele_deltaphiele.push_back(ielectrons->deltaPhiEleClusterTrackAtCalo());
    ele_deltaetain.push_back(ielectrons->deltaEtaSuperClusterTrackAtVtx());
    ele_deltaphiin.push_back(ielectrons->deltaPhiSuperClusterTrackAtVtx());

    // Shower Shape
    ele_sigmaietaieta.push_back((ielectrons->showerShape()).sigmaIetaIeta); //  ielectrons->
    ele_sigmaetaeta.push_back((ielectrons->showerShape()).sigmaEtaEta); //ielectrons->sigmaEtaEta() ;
    ele_sigmaiphiiphi.push_back((ielectrons->showerShape()).sigmaIphiIphi);
    ele_e15.push_back((ielectrons->showerShape()).e1x5); // ;ielectrons->e1x5() ;
    ele_e25max.push_back((ielectrons->showerShape()).e2x5Max); // ;ielectrons->e2x5Max() ;
    ele_e55.push_back((ielectrons->showerShape()).e5x5); // ;ielectrons->e5x5() ;
    ele_r9.push_back((ielectrons->showerShape()).r9);
    //ele_e1.push_back(    FIXME ;
    //ele_e33.push_back(   FIXME ;
    //
    // Old-Style ShowerShape
    edm::InputTag  reducedBarrelRecHitCollection("reducedEcalRecHitsEB");
    edm::InputTag  reducedEndcapRecHitCollection("reducedEcalRecHitsEE");
    //noZS::
    // noZS::EcalClusterLazyTools lazyToolsNoZS(iEvent, iSetup, reducedBarrelRecHitCollection, reducedEndcapRecHitCollection);

    // const auto & seedCluster = ielectrons->superCluster()->seed());
    // std::vector<float> vCov = lazyToolsNoZS.localCovariances(*seedCluster);
    // std::vector<float> Cov  = lazyToolsNoZS.covariances(*seedCluster);

    ele_oldsigmaetaeta.push_back(ielectrons->full5x5_sigmaEtaEta());    //( !edm::isNotFinite(Cov[0]) ) ? sqrt(Cov[0]) : 0;
    ele_oldsigmaietaieta.push_back(ielectrons->full5x5_sigmaIetaIeta());   //( !edm::isNotFinite(vCov[0]) ) ? sqrt(vCov[0]) : 0;
    ele_oldsigmaiphiiphi.push_back(ielectrons->full5x5_sigmaIphiIphi());   //( !edm::isNotFinite(vCov[2]) ) ? sqrt(vCov[2]) : 0;
    //ele_oldsigmaietaiphi.push_back(vCov[1]; // this is missing in the struct
    ele_oldr9.push_back(ielectrons->full5x5_r9());  //lazyToolsNoZS.e3x3(*seedCluster) / ielectrons->superCluster()->rawEnergy() ;
    ele_olde15.push_back(ielectrons->full5x5_e1x5()); //lazyToolsNoZS.e1x5(*seedCluster);
    ele_olde25max.push_back(ielectrons->full5x5_e2x5Max()); //lazyToolsNoZS.e2x5Max(*seedCluster);
    ele_olde55.push_back(ielectrons->full5x5_e5x5());       // lazyToolsNoZS.e5x5(*seedCluster);
    ele_oldhe.push_back(ielectrons->full5x5_hcalOverEcal());
    ele_oldhebc.push_back( ielectrons->full5x5_hcalOverEcalBc());
    // hcal stuff is not filled
    //electron.full5x5_setShowerShape(ss);
    //electron.full5x5_setSigmaIetaIphi(sigmaIetaIphi);

    // E/P combination
    ele_ecalE.push_back(ielectrons->ecalEnergy());
    ele_ecalErr.push_back(ielectrons->ecalEnergyError());
    ele_trackErr.push_back(ielectrons->trackMomentumError());
    ele_combErr.push_back(ielectrons->p4Error(GsfElectron::P4_COMBINATION));
    ele_PFcombErr.push_back(ielectrons->p4Error(GsfElectron::P4_PFLOW_COMBINATION));
    //cout << "Errors (ecal/track/p4comb/PFprcomb) :" <<ele_ecalErr.push_back( <<" " << ele_trackErr.push_back(<<" "<< ele_combErr.push_back(<<" "<< ele_PFcombErr.push_back( <<endl;

    //regression
    //ele_ecalRegressionEnergy.push_back(  = ielectrons->ecalRegressionEnergy());
    //ele_ecalRegressionError.push_back(ielectrons->ecalRegressionError());
    // ele_ecalTrackRegressionEnergy.push_back(ielectrons->ecalTrackRegressionEnergy());
    //     ele_ecalTrackRegressionError.push_back(ielectrons->ecalTrackRegressionError());
    //     ele_ecalScale.push_back(ielectrons->ecalScale());
    //     ele_ecalSmear.push_back(ielectrons->ecalSmear());
    //     ele_ecalRegressionScale.push_back(ielectrons->ecalRegressionScale());
    //     ele_ecalRegressionSmear.push_back(ielectrons->ecalRegressionSmear());
    //     ele_ecalTrackRegressionScale.push_back(ielectrons->ecalTrackRegressionScale());
    //     ele_ecalTrackRegressionSmear.push_back(ielectrons->ecalTrackRegressionSmear());

    // 		cout << "REGRESSION " << ele_ecalRegressionEnergy.push_back(<< " " << ele_ecalRegressionError.push_back( << endl;
    //
    ele_mva.push_back( 0); //ielectrons->mva() ;
    //
    if (ielectrons->isEB()) ele_isbarrel.push_back(1);
    else  ele_isbarrel.push_back(0);
    if (ielectrons->isEE()) ele_isendcap.push_back(1);
    else  ele_isendcap.push_back(0);
    if (ielectrons->isEBEtaGap()) ele_isEBetaGap.push_back(1);
    if (ielectrons->isEBPhiGap()) ele_isEBphiGap.push_back(1);
    if (ielectrons->isEEDeeGap()) ele_isEEdeeGap.push_back(1);
    if (ielectrons->isEERingGap()) ele_isEEringGap.push_back(1);
    if (ielectrons->ecalDrivenSeed()) ele_isecalDriven.push_back(1);
    if (ielectrons->trackerDrivenSeed()) ele_istrackerDriven.push_back(1);
    //

    // -----------------------------------------------------------------
    // Tracking Variables
    // -----------------------------------------------------------------
    ele_lost_hits.push_back(ielectrons->gsfTrack()->lost()); //numberOfLostHits();
    ele_valid_hits.push_back(ielectrons->gsfTrack()->found()); //numberOfValidHits() ;
    ele_gsfchi2.push_back(ielectrons->gsfTrack()->normalizedChi2());

    //double hits = ielectrons->gsfTrack()->hitPattern().trackerLayersWithMeasurement();
    //cout << "hits = " << hits << " valid = " << ele_valid_hits.push_back( << endl;
    ele_gsfhits.push_back(ielectrons->gsfTrack()->hitPattern().trackerLayersWithMeasurement());

    bool validKF=false;
    reco::TrackRef myTrackRef = ielectrons->closestCtfTrackRef();
    validKF = myTrackRef.isNonnull();
    double kfchi2 = validKF ? myTrackRef->normalizedChi2() : 0 ;
    int kfhits = validKF ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1;
    ele_kfchi2.push_back(kfchi2);  //ielectrons->track()->normalizedChi2() : 0 ;
    ele_kfhits.push_back(kfhits);    //ielectrons->track()->hitPattern().trackerLayersWithMeasurement() : 0 ;
    //
    // 		ele_dxyB.push_back(ielectrons->gsfTrack()->dxy(bs.position()) ;
    ele_dxy.push_back(ielectrons->gsfTrack()->dxy());
    // 		ele_dzB.push_back(  = ielectrons->gsfTrack()->dz(bs.position()) ;
    ele_dz.push_back(ielectrons->gsfTrack()->dz());
    // 		ele_dszB.push_back(ielectrons->gsfTrack()->dsz(bs.position()) ;
    ele_dsz.push_back(ielectrons->gsfTrack()->dsz());
    ele_dzPV.push_back(ielectrons->gsfTrack()->dz(pv->position()));
    //
    // Isolation variables
    // ele_tkSumPt_dr03.push_back(      ielectrons->dr03TkSumPt() ;
    //     ele_ecalRecHitSumEt_dr03.push_back(      = ielectrons->dr03EcalRecHitSumEt() ;
    //     ele_hcalDepth1TowerSumEt_dr03.push_back(ielectrons->dr03HcalDepth1TowerSumEt() ;
    //     ele_hcalDepth2TowerSumEt_dr03.push_back(ielectrons->dr03HcalDepth2TowerSumEt() ;
    //     ele_tkSumPt_dr04.push_back(      ielectrons->dr04TkSumPt() ;
    //     ele_ecalRecHitSumEt_dr04.push_back(      = ielectrons->dr04EcalRecHitSumEt() ;
    //     ele_hcalDepth1TowerSumEt_dr04.push_back(ielectrons->dr04HcalDepth1TowerSumEt() ;
    //     ele_hcalDepth2TowerSumEt_dr04.push_back(ielectrons->dr04HcalDepth2TowerSumEt() ;
    //
    // Custom HCAL
    //  m_calotowers = new edm::Handle<CaloTowerCollection>() ;
    //     if (!iEvent.getByLabel("towerMaker",*m_calotowers)) //hcalTowers_
    //       { edm::LogError("ElectronHcalHelper::readEvent")<<"failed to get the hcal towers of label "<<hcalTowers_ ; }

    //     edm::Ref<pat::ElectronCollection> electronEdmRef(electronsCol, counter);
    //     double egHcalIsoConeSizeOutSmall = 0.3;
    //     double egHcalIsoConeSizeIn       = 0.0, egHcalIsoPtMin=0.0;
    //     int egHcalDepth1 = 1;
    //     int egHcalDepth2 = 2;
    //hadDepth1Isolation03_  = new EgammaTowerIsolation(egHcalIsoConeSizeOutSmall,egHcalIsoConeSizeIn,egHcalIsoPtMin,egHcalDepth1,m_calotowers->product()) ;
    //hadDepth2Isolation03_  = new EgammaTowerIsolation(egHcalIsoConeSizeOutSmall,egHcalIsoConeSizeIn,egHcalIsoPtMin,egHcalDepth2,m_calotowers->product()) ;
    //double hcalDepth1TowerSumEt03 = hadDepth1Isolation03_->getTowerEtSum(&(*electronEdmRef)); //ielectrons)); //electronRef));
    //double hcalDepth2TowerSumEt03 = hadDepth2Isolation03_->getTowerEtSum(&(*electronEdmRef)); // electronEdmRefielectrons)); //electronRef));
    //ele_HCALFullConeSum.push_back(  = hcalDepth1TowerSumEt03+hcalDepth2TowerSumEt03;
    //&((*EleHandle)[i])
    //
    // Conversion Rejection
    ele_conv_dcot.push_back(ielectrons->convDist()); //userFloat("dcot");
    ele_conv_dist.push_back(ielectrons->convDcot()); //ielectrons->userFloat("dist");
    ele_conv_radius.push_back(ielectrons->convRadius());

    // FIXME: Always returns 0 :(
    //cout << " dcot = " << ielectrons->userFloat("dcot") << " dist = " << ielectrons->userFloat("dist") << endl;

    ele_expected_inner_hits.push_back(ielectrons->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS));
    //ielectrons->gsfTrack()->trackerExpectedHitsInner().numberOfHits());


    bool vtxFitConversion = ConversionTools::hasMatchedConversion(*ielectrons, conversions_h, beamSpot.position());
    ele_vtxconv.push_back(vtxFitConversion);

    //

    // -----------------------------------------------------------------
    // PF Isolation for electrons in the cone 0.4
    // -----------------------------------------------------------------
    ele_pfChargedHadIso.push_back( (ielectrons->pfIsolationVariables()).sumChargedHadronPt); //chargedHadronIso());
    ele_pfNeutralHadIso.push_back(  (ielectrons->pfIsolationVariables()).sumNeutralHadronEt); //neutralHadronIso());
    ele_pfPhotonIso.push_back(  (ielectrons->pfIsolationVariables()).sumPhotonEt); //photonIso());
    //
    ele_pfChargedIso.push_back( (ielectrons->pfIsolationVariables()).sumChargedParticlePt);
    ele_pfSumPUIso.push_back((ielectrons->pfIsolationVariables()).sumPUPt);
    //ele_pfCombRelIso.push_back(      = LeptonIsoHelper::combRelIsoPF(lepton_setup, lepton_setup, _PU_Elerho, *ielectrons);
    //FIXME these two for applying the beta corrections
    //ele_pfChargedHadPUIso.push_back(  = ielectrons->chargedAllIso());
    //ele_pfChargedHadPUIso.push_back(  = ielectrons->puChargedHadronIso());

    // -----------------------------------------------------------------
    // SIP3D
    // -----------------------------------------------------------------
    //default values for IP3D
    double ele_ip3D     = -999.0;
    double ele_ip3D_err = 999.; //fMVAVar_ip3dSig = 0.0;
    double d0_corr = 999;
    double d0_err   = 999;

    if (ielectrons->gsfTrack().isNonnull()) {
      const double gsfsign   = ( (-ielectrons->gsfTrack()->dxy(pv->position()))   >=0 ) ? 1. : -1.;

      const reco::TransientTrack &tt = thebuilder.build(ielectrons->gsfTrack()); //transientTrackBuilder.build(ielectrons->gsfTrack());
      const std::pair <bool,Measurement1D> &ip3dpv =  IPTools::absoluteImpactParameter3D(tt, *pv);
      std::pair<bool,Measurement1D> result   = IPTools::absoluteTransverseImpactParameter(tt, * pv);

      if (ip3dpv.first) {
	double ip3d = gsfsign*ip3dpv.second.value();
	double ip3derr = ip3dpv.second.error();

	ele_ip3D     = ip3d;
	ele_ip3D_err = ip3derr;
	// fMVAVar_ip3dSig = ip3d/ip3derr;
      } // if ip3dpv.first

      if(result.first) {
	d0_corr = result.second.value();
	d0_err   = result.second.error();
      } // if d0
    } // if gsftrack.isNonnull

    ele_IP.push_back(ele_ip3D); //fabs(ielectrons->dB(pat::Electron::PV3D));
    ele_IPError.push_back(ele_ip3D_err); //ielectrons->edB(pat::Electron::PV3D);
    double ele_sip = -999; if(ele_ip3D_err!=0) ele_sip = ele_ip3D / ele_ip3D_err;
    ele_SIP.push_back(ele_sip); //ele_IP.push_back(/ele_IPError.push_back(;

    ele_d0.push_back(d0_corr);
    ele_d0err.push_back(d0_err);

    // -----------------------------------------------------------------
    // Get SuperCluster Informations
    // -----------------------------------------------------------------
    //cout << " SuperCluster "<< endl;
    //	if(ielectrons->ecalDrivenSeed()) {
    reco::SuperClusterRef sclRef = ielectrons->superCluster();
    //cout << " SuperClusterRef" << endl;
    ////math::XYZPoint sclPos        = ielectrons->superClusterPosition();
    //      cout << " pflow" << endl;

    //if (!ielectrons->ecalDrivenSeed() && ielectrons->trackerDrivenSeed())
    //sclRef = ielectrons->pflowSuperCluster();

    double R  = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y() +sclRef->z()*sclRef->z());
    double Rt = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y());
    ele_sclRawE.push_back(sclRef->rawEnergy());

    ele_sclE.push_back(sclRef->energy());
    // 		ele_sclE.push_back(     = ielectrons->correctedEcalEnergy();  //for 5XY
    ele_sclEt.push_back(sclRef->energy()*(Rt/R));
    ele_sclEta.push_back(sclRef->eta());
    ele_sclPhi.push_back(sclRef->phi());
    ele_sclNclus.push_back(sclRef->clustersSize());

    //cout << " etawidth = " << sclRef->etaWidth() << endl;

    ele_sclphiwidth.push_back(sclRef->phiWidth());
    ele_scletawidth.push_back(sclRef->etaWidth());

    // cout << " E = " << ielectrons->ecalEnergy() << " E corr = " << ielectrons->correctedEcalEnergy()  << " SCLE = " <<  sclRef->energy() << endl;

    // --------------
    // Sub-Clusters
    // --------------
    int countersub = 0;
    //eSubClusters_ = 0.;
    // Store subclusters
    reco::CaloCluster_iterator itscl  = sclRef->clustersBegin();
    reco::CaloCluster_iterator itsclE = sclRef->clustersEnd();

    //cout << " sub clusters" << endl;
    if(counter<50) {
        for(; itscl < itsclE ; ++itscl) {
            bool isseed = false;
            if((*itscl)==ielectrons->superCluster()->seed()) isseed=true; // continue; // skip seed cluster
            //theBasicClusters_.push_back(&(**itscl));
            //eSubClusters_ += (*itscl)->energy();
            ele_sclsubE[counter][countersub]      = (*itscl)->energy();
            //ele_sclsubE.at(counter).push_back( (*itscl)->energy() );
            ele_sclsubEta[counter][countersub]    = (*itscl)->eta();
            ele_sclsubPhi[counter][countersub]    = (*itscl)->phi();
            ele_sclsubisseed[counter][countersub] = isseed;
            //N subclusters?->sclNclus... to be checked
            countersub++;
        }
    }
    // sort subclusters

    //	} // if ECAL driven

    //cout << " Etawidth = " <<  ele_scletawidth.push_back( << endl;
    //cout << "" << endl;

    // -----------------------------------------------------------------
    // Get PreShower Informations
    // -----------------------------------------------------------------
    //cout << " preshower" << endl;
    ele_psE.push_back(sclRef->preshowerEnergy());
    //T_Elec_PreShowerOverRaw->push_back((eleIt->superCluster()->rawEnergy()>0) ? eleIt->superCluster()->preshowerEnergy() / eleIt->superCluster()->rawEnergy() : -1)

    // -----------------------------------------------------------------
    //fbrem
    // -----------------------------------------------------------------
    ele_fbrem.push_back( ielectrons->fbrem());
    ele_trackfbrem.push_back(ielectrons->trackFbrem());
    ele_SCfbrem.push_back(ielectrons->superClusterFbrem());
    ele_pfSCfbrem.push_back(0); //ielectrons->pfSuperClusterFbrem()); // Should be identical to the previous one...
    ele_eClass.push_back(ielectrons->classification());
    ele_nbrem.push_back(ielectrons->numberOfBrems());

    // FOR MVA...
    //  ele_mvafbrem.push_back(ielectrons->fbrem());
    //     ele_mvadetain.push_back(ielectrons->deltaEtaSuperClusterTrackAtVtx());
    //     ele_mvadphiin.push_back(ielectrons->deltaPhiSuperClusterTrackAtVtx());
    //     ele_mvasieie.push_back(ielectrons->sigmaIetaIeta());
    //     ele_mvahoe.push_back(ielectrons->hcalOverEcal());
    //     ele_mvaeop.push_back(ielectrons->eSuperClusterOverP());
    //     ele_mvae1x5e5x5.push_back((ielectrons->e5x5()) !=0. ? ielectrons->e1x5()/ielectrons->e5x5() : -1. ;
    //     ele_mvaeleopout.push_back(ielectrons->eEleClusterOverPout());
    //     bool validKF= (ielectrons->track().isNonnull());
    //     ele_mvakfchi2.push_back(validKF ? ielectrons->track()->normalizedChi2() : 0 ;
    //     ele_mvakfhits.push_back(validKF ? ielectrons->track()->hitPattern().trackerLayersWithMeasurement() : 0 ;
    //     ele_mvamishits.push_back(ielectrons->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits());
    //     ele_mvadist.push_back(ielectrons->convDist());
    //     ele_mvadcot.push_back(ielectrons->convDcot());
    //     ele_mvaeta.push_back(ielectrons->eta());
    //     ele_mvapt.push_back(ielectrons->pt());
    //     ele_mvaecalseed.push_back(ielectrons->ecalDrivenSeed());

    // -----------------------------------------------------------------
    // Electron ID electronsCol
    // -----------------------------------------------------------------
    edm::Ref<reco::GsfElectronCollection> electronRef(electronsCol, counter); //electronsCollection,i); i++; //reference to the electron
    //cout << "MVA = " << mapMVA[electronRef] << endl;
    //ele_mvaphys14.push_back(mapMVA_phys14[electronRef]) ;
//    ele_mvaphys14fix.push_back(mapMVA_phys14fix[electronRef]);

    //T_Elec_MVAoutput->push_back(mapMVA[electronRef]);


    // -----------------------------------------------------------------
    // PF Candidates around Electrons
    // -----------------------------------------------------------------
     //here loop over the PF candidates around the electrons
    int numberOfPFAroundEle = 0;
    if (counter>49) continue;
    for(unsigned i=0; i<pfCandColl->size(); i++) {
      const reco::PFCandidate& pfc = (*pfCandColl)[i];
      float dR = deltaR(ielectrons->eta(), ielectrons->phi(), pfc.momentum().Eta(), pfc.momentum().Phi());
      int pfID = pfc.particleId(); // 0/X:undefined, 1/h:  charged hadrons, 2/e:electron, 3/mu:muon, 4/gamma:photon, 5/h0:neutral hadron, 6 and 7: HF
      if ( (dR < 0.4) ) //&& (pfID == 1 || pfID == 4 || pfID == 5) )
      {

	float pf_E = pfc.energy();
	float pf_pt =  pfc.pt();
	float pf_eta = pfc.eta(); //momentum().Eta();
	float pf_phi = pfc.phi(); //momentum().Phi();
	float pf_mva_nog = pfc.mva_nothing_gamma();
	float pf_mva_epi = pfc.mva_e_pi();
	//
	//cout << " eta = " << pf_eta << " eta moi = " << pfc.eta() << endl;
	//cout << " energy moi = " << pfc.energy() << endl;
	//
	float dz_pf =  -999.;
	float dxy_pf =   -999.;
	float vx_pf =  -999.;
	float vy_pf =  -999.;
	float vz_pf =  -999.;

	if (pfID == 1) { //|| pfID == 2 || pfID == 3) { // charged hadrons || electrons || muons
	  dz_pf =  pfc.trackRef()->dz();
	  //float dz_ele =  ielectrons->gsfTrack()->dz();
	  dxy_pf =  pfc.trackRef()->dxy();
	  //float dxy_ele =  ielectrons->gsfTrack()->dxy();
	  vx_pf = pfc.vx();
	  vy_pf = pfc.vy();
	  vz_pf = pfc.vz();
	} //

	ele_pf_id[counter][numberOfPFAroundEle] = pfID; //[.push_back(pfID);
	ele_pf_E[counter][numberOfPFAroundEle] =pf_E;
	ele_pf_eta[counter][numberOfPFAroundEle] =pf_eta; // a.push_back(pf_eta);
	ele_pf_phi[counter][numberOfPFAroundEle] = pf_phi; //push_back(pf_phi);
	ele_pf_pt[counter][numberOfPFAroundEle] = pf_pt; //push_back(pf_pt);
	ele_pf_dz[counter][numberOfPFAroundEle] = dz_pf; //.push_back(dz_pf);
	ele_pf_dxy[counter][numberOfPFAroundEle] = dxy_pf; // y.push_back(dxy_pf);
	ele_pf_vx[counter][numberOfPFAroundEle] = vx_pf; //push_back(vx_pf);
	ele_pf_vy[counter][numberOfPFAroundEle] = vy_pf; //push_back(vy_pf);
	ele_pf_vz[counter][numberOfPFAroundEle] = vz_pf; //push_back(vz_pf);
	ele_pf_mva_nog[counter][numberOfPFAroundEle] = pf_mva_nog; //push_back(pf_mva_nog);
	ele_pf_mva_epi[counter][numberOfPFAroundEle] =  pf_mva_epi; //push_back(pf_mva_epi);

	++numberOfPFAroundEle;
      } // if dr04
    } // for loop on PFs
    ele_pf_number.push_back(numberOfPFAroundEle);
    //cout << "PF = " << numberOfPFAroundEle << endl;




    ++counter;
  } // for loop on gsfelectrons

  if(counter>49) { ele_N = 50; cout << "Number of electrons>49, electrons_N set to 50" << endl;}

} // end of FillElectrons

// ====================================================================================
void Ntuplizer::FillMET (const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{

  // caloMET object (negative vector sum of calorimeter towers)
  //edm::Handle< edm::View<reco::CaloMET> > caloMEThandle;
  //iEvent.getByLabel("met", caloMEThandle);

  // MET object that corrects the basic calorimeter MET for muons
  // edm::Handle< edm::View<reco::CaloMET> > muCorrMEThandle;
  //   iEvent.getByLabel("corMetGlobalMuons", muCorrMEThandle);

  // MET object that corrects the basic calorimeter MET for muons and tracks
  //edm::Handle< edm::View<reco::MET> > tcMEThandle;
  //iEvent.getByLabel("tcMet", tcMEThandle);

  // MET object built as the (negative) vector sum of all particles (PFCandidates) reconstructed in the event
  // 	edm::Handle< edm::View<pat::MET> > pfMEThandle;
  // 	iEvent.getByLabel("patMETs", pfMEThandle);

  //edm::Handle< edm::View<cmg::BaseMET> > pfMEThandle;
 // edm::Handle< edm::View<reco::PFMET> > pfMEThandle;
 // iEvent.getByLabel("pfMet", pfMEThandle);
//
//
 // // PFMET
 // _met_pf_et  = (pfMEThandle->front() ).et();
 // _met_pf_px  = (pfMEThandle->front() ).px();
 // _met_pf_py  = (pfMEThandle->front() ).py();
 // _met_pf_phi = (pfMEThandle->front() ).phi();
 // _met_pf_set = (pfMEThandle->front() ).sumEt();
 // _met_pf_sig = (pfMEThandle->front() ).mEtSig();

} // end of Fill MET



// ====================================================================================
void Ntuplizer::FillTruth(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
    Handle<std::vector<reco::GenParticle>> genParticlesHandle;
    iEvent.getByToken(genparticles_, genParticlesHandle);
    for (std::vector<reco::GenParticle>::const_iterator it_p = genParticlesHandle->begin();
         it_p != genParticlesHandle->end(); ++it_p) {
      gen_eta_.push_back(it_p->eta());
      gen_phi_.push_back(it_p->phi());
      gen_pt_.push_back(it_p->pt());
      gen_energy_.push_back(it_p->energy());
      gen_charge_.push_back(it_p->charge());
      gen_pdgid_.push_back(it_p->pdgId());
      gen_status_.push_back(it_p->status());

      std::vector<int> daughters(it_p->daughterRefVector().size(), 0);
      for (unsigned j = 0; j < it_p->daughterRefVector().size(); ++j) {
        daughters[j] = static_cast<int>(it_p->daughterRefVector().at(j).key());
      }
      gen_daughters_.push_back(daughters);
    }
} // end of FillTruth




// ------------ method called once each job just after ending the event loop  ------------
void
Ntuplizer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------

void
Ntuplizer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}


// ------------ method called when ending the processing of a run  ------------

void
Ntuplizer::endRun(edm::Run const&, edm::EventSetup const&)
{
}


// ------------ method called when starting to processes a luminosity block  ------------
/*
void
Ntuplizer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
Ntuplizer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/


 // ====================================================================================================
void Ntuplizer::setMomentum(TLorentzVector &myvector, const LorentzVector & mom)
// ====================================================================================================
{

  myvector.SetPx (mom.Px());
  myvector.SetPy (mom.Py());
  myvector.SetPz (mom.Pz());
  myvector.SetE (mom.E());

}


 // ====================================================================================================
void Ntuplizer::Init()
// ====================================================================================================
{

  _PU_N = 0;

  _vtx_N = 0;

  ele_N = 0;

  _met_pf_et  = 0.;
  _met_pf_px  = 0.;
  _met_pf_py  = 0.;
  _met_pf_phi = 0.;
  _met_pf_set = 0.;
  _met_pf_sig = 0.;

  //  for (int i = 0 ; i < 50 ; ++i) {
    //electrons
    ele_echarge.clear();
    ele_eta.clear();
    ele_phi.clear();
    ele_pt.clear();
    ele_he.clear();
    ele_hebc.clear();
    ele_oldhe.clear();
    ele_oldhebc.clear();
    ele_eseedpout.clear();
    ele_ep.clear();
    ele_eseedp.clear();
    ele_eelepout.clear();
    ele_deltaetaseed.clear();
    ele_deltaetaele.clear();
    ele_deltaphiseed.clear();
    ele_deltaphiele.clear();
    ele_deltaetain.clear();
    ele_deltaphiin.clear();
    ele_sigmaietaieta.clear();
    ele_sigmaiphiiphi.clear();
    ele_sigmaetaeta.clear();
    ele_e15.clear();
    ele_e25max.clear();
    ele_e55.clear();
    ele_e1.clear();
    ele_r9.clear();
    //
    ele_oldsigmaetaeta.clear(); // [i]   = 0;
    ele_oldsigmaietaieta.clear();
    ele_oldsigmaiphiiphi.clear();
    ele_oldsigmaietaiphi.clear();
    ele_oldr9.clear(); //           = 0;
    ele_olde15.clear(); //       [i]           = 0;
    ele_olde25max.clear(); //       [i]        = 0;
    ele_olde55.clear(); //       [i]           = 0;
    //
    //ele_e33.clear();
    //ele_e2overe9.clear();
    ele_pin_mode.clear();
    ele_pout_mode.clear();
    ele_pTin_mode.clear();
    ele_pTout_mode.clear();
    //
    ele_fbrem.clear();
    ele_SCfbrem.clear();
    ele_pfSCfbrem.clear();
    ele_trackfbrem.clear();
    ele_nbrem.clear();
    //
    ele_mva.clear();
    ele_isbarrel.clear();
    ele_isendcap.clear();
    ele_isEBetaGap.clear();
    ele_isEBphiGap.clear();
    ele_isEEdeeGap.clear();
    ele_isEEringGap.clear();
    ele_isecalDriven.clear();
    ele_istrackerDriven.clear();
    ele_eClass.clear();
    ele_valid_hits.clear();
    ele_lost_hits.clear();
    ele_gsfchi2.clear();
    ele_dxyB.clear();
    ele_dxy.clear();
    ele_dzB.clear();
    ele_dz.clear();
    ele_dszB.clear();
    ele_dsz.clear();
    //  ele_tkSumPt_dr03.clear();
    //     ele_ecalRecHitSumEt_dr03.clear();
    //     ele_hcalDepth1TowerSumEt_dr03.clear();
    //     ele_hcalDepth2TowerSumEt_dr03.clear();
    //     ele_tkSumPt_dr04.clear();
    //     ele_ecalRecHitSumEt_dr04.clear();
    //     ele_hcalDepth1TowerSumEt_dr04.clear();
    //     ele_hcalDepth2TowerSumEt_dr04.clear();
    ele_conv_dcot.clear();
    ele_conv_dist.clear();
    ele_conv_radius.clear();
    ele_expected_inner_hits.clear();
    ele_vtxconv.clear();
    //
    ele_pfChargedHadIso.clear(); // [i]   = 0;
    ele_pfNeutralHadIso.clear(); // [[i]   = 0;
    ele_pfPhotonIso.clear(); // [[i]       = 0;

    ele_pfChargedIso.clear();
    ele_pfSumPUIso.clear(); // [[i]   = 0;
    //ele_pfChargedHadPUIso.clear();
    //ele_pfCombRelIso.clear();

    ele_pf_number.clear();

    ele_IP.clear();
    ele_IPError.clear();
    ele_SIP.clear();

    ele_dzPV.clear();
    ele_d0.clear();
    ele_d0err.clear();

    //
    ele_sclRawE.clear(); // [[i]=0;
    ele_sclE.clear(); // [[i]=0;
    ele_sclEt.clear(); // [[i]=0;
    ele_sclEta.clear(); // [[i]=0;
    ele_sclPhi.clear(); // [[i]=0;
    ele_sclNclus.clear(); // [[i]=0;
    ele_sclphiwidth.clear(); // [[i]=0;
    ele_scletawidth.clear(); // [[i]=0;

    ele_ecalE.clear(); // [[i]=0;
    ele_ecalErr.clear(); // [[i]=0;
    ele_trackErr.clear(); // [[i]=0;
    ele_combErr.clear(); // [[i]=0;
    ele_PFcombErr.clear(); // [[i]=0;

    ele_ecalRegressionEnergy.clear(); // [[i]  = 0;
    ele_ecalRegressionError.clear();
    ele_ecalTrackRegressionEnergy.clear(); // [[i]  = 0;
    ele_ecalTrackRegressionError.clear(); // [[i]  = 0;
    ele_ecalScale.clear(); // [[i]  = 0;
    ele_ecalSmear.clear(); // [[i]  = 0;
    ele_ecalRegressionScale.clear(); // [[i]  = 0;
    ele_ecalRegressionSmear.clear(); // [[i]  = 0;
    ele_ecalTrackRegressionScale.clear(); // [[i]  = 0;
    ele_ecalTrackRegressionSmear.clear(); // [[i]  = 0;

    //  ele_mvafbrem[i]=0;
    //     ele_mvadetain[i]=0;
    //     ele_mvadphiin[i]=0;
    //     ele_mvasieie[i]=0;
    //     ele_mvahoe[i]=0;
    //     ele_mvaeop[i]=0;
    //     ele_mvae1x5e5x5[i]=0;
    //     ele_mvaeleopout[i]=0;
    ele_kfchi2.clear(); // [i]=0;
    ele_kfhits.clear(); //[i]=-1;
    ele_gsfhits.clear(); //[i] = -1;
    //     ele_mvamishits[i]=0;
    //     ele_mvadist[i]=0;
    //     ele_mvadcot[i]=0;
    //     ele_mvaeta[i]=0;
    //     ele_mvapt[i]=0;
    //     ele_mvaecalseed[i]=0;

    //ele_sclphiwidth[i]= 0;
    //ele_scletawidth[i]= 0;
    ele_psE.clear(); //[i] = 0.;

    ele_mvaphys14.clear(); //[i] = 0.;
    ele_mvaphys14fix.clear(); //[i] = 0.;

    //  } // for loop on electrons
    //ele_sclsubE.clear();

    // should work at some point with vector of vector or something similar...
    for(int i=0;i<50;i++) {

      for(int j=0;j<20;j++) {
	ele_sclsubE[i][j] = 0; //counter][countersub];
	ele_sclsubEta[i][j] = 0;
	ele_sclsubPhi[i][j] = 0;
	ele_sclsubisseed[i][j] = 0;
      } // for loop on subclusters

      for(int j=0;j<100;j++) {
	ele_pf_id[i][j] = 0;

	ele_pf_E[i][j] = 0;
	ele_pf_eta[i][j] = 0;
	ele_pf_phi[i][j] = 0;
	ele_pf_pt[i][j] = 0;
	ele_pf_dz[i][j] = 0;
	ele_pf_dxy[i][j] = 0;
	ele_pf_vx[i][j] = 0;
	ele_pf_vy[i][j] = 0;
	ele_pf_vz[i][j] = 0;
	ele_pf_mva_nog[i][j] = 0;
	ele_pf_mva_epi[i][j] = 0;
      } // for loop on pf candidates

    }  // for loop on electrons

  // Gen Informations
  ////////////////////
  // reco::GenParticles
  //
  gen_eta_.clear();
  gen_phi_.clear();
  gen_pt_.clear();
  gen_energy_.clear();
  gen_charge_.clear();
  gen_pdgid_.clear();
  gen_status_.clear();
  gen_daughters_.clear();

}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Ntuplizer);
