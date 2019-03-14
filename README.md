# reco-ntuples

Home of the Ntuplizer for the HGCAL reconstruction software studies

Ntuple content definitions can be found at [Definitions.md](Definitions.md).

This version is based on >= CMSSW_10_6_0_pre2, instructions below should provide you with the version for the TDR sample production. Please check if there are later `CMSSW_10_6_X` versions available to profit from bugfixes.

```shell
cmsrel CMSSW_10_6_0_pre2
cd CMSSW_10_6_0_pre2/src
cmsenv
git clone git@github.com:CMS-HGCAL/reco-ntuples.git RecoNtuples
cd RecoNtuples
git checkout -b topic_${USER}
cd ../
scram b -j4
```

The input file needs to be step3 (i.e. RECO). Example configs are provided in [HGCalAnalysis/test](HGCalAnalysis/test).

Mind that depending on your RECO input file, you need to set `inputTag_HGCalMultiCluster` in the config part for the `EDAnalyzer` of `HGCalAnalysis` (i.e. the ntupliser) to either `hgcalMultiClusters` (newer releases) or `hgcalLayerClusters` (older releases).

For older versions, please check the [releases tab](https://github.com/CMS-HGCAL/reco-ntuples/releases).
