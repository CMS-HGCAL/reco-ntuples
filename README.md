# reco-ntuples (Electron ID branch)
Home of the Ntuplizer for the HGCAL reconstruction software studies

Ntuple content definitions can be found at [Definitions.md](Definitions.md).

This version is based on >= CMSSW_9_3_2, instructions below should provide you with the version for the TDR sample production.

```bash
cmsrel CMSSW_9_3_2                                                  
cd CMSSW_9_3_2/src                                                        
cmsenv                                             
git clone -b topic_eleID_helper git@github.com:CMS-HGCAL/reco-ntuples.git RecoNtuples
git clone git@github.com:CMS-HGCAL/EgammaTools.git EgammaTools
scram b -j 4
cd RecoNtuples/HGCalAnalysis/test
mkdir ntuples #important to create, otherwise change output file name in config
cmsRun ele15_config_EleID.py   
```

The input file needs to be step3 (i.e. RECO). Example configs are provided in [HGCalAnalysis/test](HGCalAnalysis/test).
