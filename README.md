# reco-ntuples

Home of the Ntuplizer for the HGCAL reconstruction software studies

Ntuple content definitions can be found at [Definitions.md](Definitions.md).

This version is based on >= `CMSSW_11_0_0_pre4`. Please check if there are later `CMSSW_11_X_Y` versions available to profit from bugfixes.

```shell
cmsrel CMSSW_11_0_0_pre4
cd CMSSW_11_0_0_pre4/src
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
