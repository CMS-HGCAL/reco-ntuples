# reco-ntuples
Home of the Ntuplizer for the HGCAL reconstruction software studies

Ntuple content definitions can be found at [Definitions.md](Definitions.md).

This version is based on CMSSW_9_0_X. The recommended procedure for production is to install the [reco-prodtools](https://github.com/CMS-HGCAL/reco-prodtools) repo and follow the instructions there. NB you may need to use the pending PR version for latest updates.

```
cmsrel CMSSW_9_0_0_pre2
cd CMSSW_9_0_pre2/src
cmsenv
git cms-merge-topic edjtscott:hgcal_multiclustering_sensordependent_piondev
git checkout -b topic_${USER}
git clone git@github.com:CMS-HGCAL/reco-ntuples.git RecoNtuples
scram b -j9
```

The input file needs to be step3 (i.e. RECO). You can generate this using the reco-prodtools repo mentioned above.
