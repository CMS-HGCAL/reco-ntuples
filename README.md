# reco-ntuples
Home of the Ntuplizer for the HGCAL reconstruction software studies

Ntuple content definitions can be found at [Definitions.md](Definitions.md).

This version is based on >= CMSSW_9_3_X, instructions below should provide you with the version for the TDR sample production. Please check if there are later `CMSSW_9_3_X` versions available to profit from bugfixes.

```
cmsrel CMSSW_9_3_2
cd CMSSW_9_3_2/src
cmsenv
git clone git@github.com:CMS-HGCAL/reco-ntuples.git RecoNtuples
cd RecoNtuples
git checkout -b topic_${USER}
cd ../
scram b -j4
```

The input file needs to be step3 (i.e. RECO). Example configs are provided in [HGCalAnalysis/test](HGCalAnalysis/test).
