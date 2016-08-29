# Ntuple content definitions

| branch name | collection name | collection type | definition  |
| ------------- | ------------- | ----- | ----- |
| event | - | - | general event info |
| particles | `mix::MergedTrackTruth` | `std::vector<TrackingParticle>` | truth level tracks/particles |
| simcluster | `mix:MergedCaloTruth` | `std::vector<SimCluster>` | Geant particle and its associated hits (DetIds) in the HGCal |
| pfcluster | `mix:MergedCaloTruth` | `std::vector<SimCluster>` | mapping of the SimCluster DetIds to the reconstructed hits |
| cluster2d | `imagingClusterHGCal` | `reco::CaloClusterCollection` | reconstructed 2D clusters |
| multicluster | - | `reco::CaloClusterCollection` | reconstructed 3D cluster built from 2D clusters
| rechits | `HGCalRecHit::HGCEERecHits` <br> `HGCalRecHit::HGCHEFRecHits` | `HGCRecHitCollection` | reconstructed calorimeter hits |
| caloparticles | `mix:MergedCaloTruth` | `std::vector<CaloParticle>` | |
