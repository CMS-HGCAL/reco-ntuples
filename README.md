# reco-ntuples
Home of the Ntuplizer for the HGCAL reconstruction software studies

Ntuple content definitions can be found at [Definitions.md](Definitions.md).

This version is based on CMSSW_9_0_X.

```
cmsrel CMSSW_9_0_0_pre2
cd CMSSW_9_0_0_pre2/src
cmsenv
git cms-merge-topic CMS-HGCAL:CMSSW_9_0_pre2
# for latest additions for multiclustering with kdtree, also merge:
git cms-merge-topic edjtscott:hgcal_multiclustering_sensordependent_piondev
git clone git@github.com:CMS-HGCAL/reco-ntuples.git RecoNtuples
cd RecoNtuples
git checkout -b topic_${USER}
cd ../
scram b -j4
```

The input file needs to be step3 (i.e. RECO). Example configs are provided in [HGCalAnalysis/test](HGCalAnalysis/test).
