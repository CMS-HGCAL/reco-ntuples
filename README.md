# reco-ntuples
Home of the Ntuplizer for the HGCAL reconstruction software studies

Ntuple content definitions can be found at [Definitions.md](Definitions.md).

This version is based on >= CMSSW_9_2_1, which contains https://github.com/cms-sw/cmssw/pull/18236

```
cmsrel CMSSW_9_2_1
cd CMSSW_9_2_1/src
cmsenv
git clone git@github.com:CMS-HGCAL/reco-ntuples.git RecoNtuples
cd RecoNtuples
git checkout -b topic_${USER}
cd ../
scram b -j4
```

The input file needs to be step3 (i.e. RECO). Example configs are provided in [HGCalAnalysis/test](HGCalAnalysis/test).
