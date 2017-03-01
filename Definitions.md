# Ntuple content definitions

| branch name | collection name | collection type | definition  |
| ------------- | ------------- | ----- | ----- |
| event | - | - | general event info |
| particles | `mix::MergedTrackTruth` | `std::vector<TrackingParticle>` | truth level tracks/particles |
| simcluster | `mix:MergedCaloTruth` | `std::vector<SimCluster>` | Geant particle and its associated hits (DetIds) in the HGCal |
| pfcluster | `particleFlowClusterHGCal` | `std::vector<reco::PFCluster>` | mapping of the SimCluster DetIds to the reconstructed hits |
| cluster2d | `imagingClusterHGCal` | `reco::CaloClusterCollection` | reconstructed layer (2D) clusters |
| multicluster | - | `reco::CaloClusterCollection` | reconstructed multi-cluster (3D) built from layer clusters |
| rechits_raw | `HGCalRecHit::HGCEERecHits` <br> `HGCalRecHit::HGCHEFRecHits` <br> `HGCalRecHit::HGCHEBRecHits` | `HGCRecHitCollection` | all reconstructed calorimeter hits |
| rechits | `HGCalRecHit::HGCEERecHits` <br> `HGCalRecHit::HGCHEFRecHits` <br> `HGCalRecHit::HGCHEBRecHits` | `HGCRecHitCollection` | reconstructed calorimeter hits associated to layer clusters |
| caloparticles | `mix:MergedCaloTruth` | `std::vector<CaloParticle>` | |
| tracks | `generalTracks` | `std::vector<reco::Track>` | tracks passing highPurity selection |
