# reco-ntuples
Home of the Ntuplizer for the HGCAL reconstruction software studies

to be added on top of CMS-HGCAL:emilio_plus_sharing:

```
cmsrel CMSSW_8_1_0_pre8
cd CMSSW_8_1_0_pre8/src
cmsenv
git cms-merge-topic CMS-HGCAL:emilio_plus_sharing
git clone git@github.com:CMS-HGCAL/reco-ntuples.git RecoNtuples
scram b -j9
cd RecoNtuples/HGCalAnalysis/test
cmsRun twogamma_pt5_eta2_nosmear_calib.py
```

The input file needs to be step3.
