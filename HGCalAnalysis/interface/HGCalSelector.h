//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Apr 22 23:25:14 2017 by ROOT version 6.06/09
// from TTree Events/
// found on file: /eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/partGun_amartell_PDGid22_nPart1_E100_900pre2_20170117/GSD/partGun_PDGid22_x50_E100.0To100.0_GSD_1.root
//////////////////////////////////////////////////////////

#ifndef HGCalSelector_h
#define HGCalSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/Provenance/interface/StoredProductProvenance.h"
#include "DataFormats/Provenance/interface/Hash.h"
#include "DataFormats/Common/interface/Wrapper.h"

#include "RecoNtuples/HGCalAnalysis/interface/AEvent.h"
#include "RecoNtuples/HGCalAnalysis/interface/AObData.h"

#include <vector>

class HGCalSelector : public TSelector {
    public :
        //TTreeReader     fReader;  //!the tree reader
        TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

        HGCalSelector(TTree * /*tree*/ =0) { }
        virtual ~HGCalSelector() { }
        virtual Int_t   Version() const { return 2; }
        virtual void    Begin(TTree *tree);
        virtual void    SlaveBegin(TTree *tree);
        virtual void    Init(TTree *tree);
        virtual Bool_t  Notify();
        virtual Bool_t  Process(Long64_t entry);
        virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
        virtual void    SetOption(const char *option) { fOption = option; }
        virtual void    SetObject(TObject *obj) { fObject = obj; }
        virtual void    SetInputList(TList *input) { fInput = input; }
        virtual TList  *GetOutputList() const { return fOutput; }
        virtual void    SlaveTerminate();
        virtual void    Terminate();
        int main();

        AEvent                    *event; // not working, don't know why
        AGenPartCollection        *genParticles;
        ARecHitCollection         *recHits;
        ACluster2dCollection      *clusters2D;
        AMultiClusterCollection   *multiClusters;
        ASimClusterCollection     *simClusters;
        APFClusterCollection      *pfClusters;
        ACaloParticleCollection   *caloParticles;

        TBranch *b_event; // not
        TBranch *b_genParticles;
        TBranch *b_recHits;
        TBranch *b_clusters2D;
        TBranch *b_multiClusters;
        TBranch *b_simClusters;
        TBranch *b_pfClusters;
        TBranch *b_caloParticles;

        ClassDef(HGCalSelector,0);
};

#endif

#ifdef HGCalSelector_cxx
void HGCalSelector::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the reader is initialized.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    //fReader.SetTree(tree);
    fChain = tree;

    event         = 0;
    genParticles  = 0;
    recHits       = 0;
    clusters2D    = 0;
    multiClusters = 0;
    simClusters   = 0;
    pfClusters    = 0;
    caloParticles = 0;

    fChain->SetBranchAddress("event", &event, &b_event);
    fChain->SetBranchAddress("genParticles", &genParticles, &b_genParticles);
    fChain->SetBranchAddress("recHits", &recHits, &b_recHits);
    fChain->SetBranchAddress("clusters2D", &clusters2D, &b_clusters2D);
    fChain->SetBranchAddress("multiClusters", &multiClusters, &b_multiClusters);
    fChain->SetBranchAddress("simClusters", &simClusters, &b_simClusters);
    fChain->SetBranchAddress("pfClusters", &pfClusters, &b_pfClusters);
    fChain->SetBranchAddress("caloParticles", &caloParticles, &b_caloParticles);

}

Bool_t HGCalSelector::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}


#endif // #ifdef HGCalSelector_cxx
