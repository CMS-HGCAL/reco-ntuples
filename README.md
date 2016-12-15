# reco-ntuples
Home of the Ntuplizer for the HGCAL reconstruction software studies

Ntuple content definitions can be found at [Definitions.md](Definitions.md).

This temporary version is based on CMSSW_9_0_X. Once github.com/cms-sw/cmssw/pull/16997 is merged, we can use a prerelease and the setup will be cleaner.

```
cmsrel CMSSW_9_0_X_2016-12-07-2300
cd CMSSW_9_0_X_2016-12-07-2300/src
cmsenv
git cms-merge-topic lgray:hgcal_cluster_speed
git cms-merge-topic edjtscott:hgcal_multiclustering_realspacecone
git checkout --theirs RecoLocalCalo/HGCalRecAlgos/src/HGCalDepthPreClusterer.cc
git add -u 
git commit -m "conflict fixed"
git checkout -b topic_${USER}
git clone git@github.com:CMS-HGCAL/reco-ntuples.git RecoNtuples
scram b -j9
cd RecoNtuples/HGCalAnalysis/test
cmsRun twogamma_pt5_eta2_nosmear_calib.py
```

The input file needs to be step3 (i.e. RECO).
