# reco-ntuples
Home of the Ntuplizer for the HGCAL reconstruction software studies

Ntuple content definitions can be found at [Definitions.md](Definitions.md).

based on CMSSW_9_0_X:

```
cmsrel CMSSW_9_0_X_2016-12-07-2300
cd CMSSW_9_0_X_2016-12-07-2300/src
cmsenv
git cms-merge-topic lgray:hgcal_cluster_speed
git cms-merge-topic edjtscott:hgcal_multiclustering_realspacecone
git checkout -b topic_${USER}
git clone git@github.com:CMS-HGCAL/reco-ntuples.git RecoNtuples
scram b -j9
cd RecoNtuples/HGCalAnalysis/test
cmsRun twogamma_pt5_eta2_nosmear_calib.py
```

The input file needs to be step3 (i.e. RECO).
