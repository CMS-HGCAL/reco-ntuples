# reco-ntuples
Home of the Ntuplizer for the HGCAL reconstruction software studies

Ntuple content definitions can be found at [Definitions.md](Definitions.md).

to be added on top of CMS-HGCAL:hgcal_clustering_development_810, based on CMSSW_8_1_0_pre15:

```
cmsrel CMSSW_8_1_0_pre15
cd CMSSW_8_1_0_pre15/src
cmsenv
git cms-merge-topic CMS-HGCAL:hgcal_clustering_development_810
git checkout CMS-HGCAL/hgcal_clustering_development_810
git checkout -b topic_${USER}
git clone git@github.com:CMS-HGCAL/reco-ntuples.git RecoNtuples
scram b -j9
cd RecoNtuples/HGCalAnalysis/test
cmsRun twogamma_pt5_eta2_nosmear_calib.py
```

The input file needs to be step3 (i.e. RECO).
