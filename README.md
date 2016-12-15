# reco-ntuples
Home of the Ntuplizer for the HGCAL reconstruction software studies

Ntuple content definitions can be found at [Definitions.md](Definitions.md).

This version is based on CMSSW_9_0_X. Further updates to the recipe will occur once github.com/cms-sw/cmssw/pull/16997 is merged into CMSSW.

```
cmsrel CMSSW_9_0_0_pre2
cd CMSSW_9_0_pre2/src
cmsenv
git cms-merge-topic edjtscott:hgcal_multiclustering_realspacecone
git checkout -b topic_${USER}
git clone git@github.com:CMS-HGCAL/reco-ntuples.git RecoNtuples
scram b -j9
cd RecoNtuples/HGCalAnalysis/test
cmsRun twogamma_pt5_eta2_nosmear_calib.py
```

The input file needs to be step3 (i.e. RECO).
