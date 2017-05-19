# reco-ntuples
Home of the Ntuplizer for the HGCAL reconstruction software studies

Ntuple content definitions can be found at [Definitions.md](Definitions.md).

This version is based on CMSSW_9_1_X.

```
cmsrel CMSSW_9_1_0_pre3
cd CMSSW_9_1_0_pre3/src
cmsenv
# the following PR updating the clustering is currently under review: https://github.com/cms-sw/cmssw/pull/18236
# if you'd like to use it, please use CMSSW_9_2_0 instead of CMSSW_9_1_0_pre3 for now
# git cms-merge-topic CMS-HGCAL:91X
git clone git@github.com:CMS-HGCAL/reco-ntuples.git RecoNtuples
cd RecoNtuples
git checkout -b topic_${USER}
cd ../
scram b -j4
```

The input file needs to be step3 (i.e. RECO). Example configs are provided in [HGCalAnalysis/test](HGCalAnalysis/test).
