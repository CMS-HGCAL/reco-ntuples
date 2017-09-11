# Ntuple content definitions

| branch name | collection name | collection type | definition  |
| ------------- | ------------- | ----- | ----- |
| event/run/lumi | - | - | general event info |
| vtx | `mix:MergedCaloTruth` | `std::vector<TrackingVertex>` | general event info |
| genpart | `g4SimHits` | `std::vector<SimTrack>` | truth level tracks/particles |
| simcluster | `mix:MergedCaloTruth` | `std::vector<SimCluster>` | Geant particle and its associated hits (DetIds) in the HGCal |
| pfcluster | `particleFlowClusterHGCal` | `std::vector<reco::PFCluster>` | mapping of the SimCluster DetIds to the reconstructed hits |
| cluster2d | `imagingClusterHGCal` | `reco::CaloClusterCollection` | reconstructed layer (2D) clusters - those that are associated to a multicluster have `cluster2d_multicluster >= 0`, which is the index of the `multiclus` in the ntuple |
| multiclus | - | `reco::CaloClusterCollection` | reconstructed multi-cluster (3D) built from layer clusters |
| rechit | `HGCalRecHit::HGCEERecHits` <br> `HGCalRecHit::HGCHEFRecHits` <br> `HGCalRecHit::HGCHEBRecHits` | `HGCRecHitCollection` | all reconstructed calorimeter hits - those that are associated to layer clusters have `rechit_cluster2d >= 0`, which is the index of the `cluster2d` in the ntuple |
| calopart | `mix:MergedCaloTruth` | `std::vector<CaloParticle>` | |
| track | `generalTracks` | `std::vector<reco::Track>` | tracks passing highPurity selection |

Mind that when running on *centrally produced/RelVal samples*, some collections used are _different_:

| branch name | collection name | collection type | definition  |
| ------------- | ------------- | ----- | ----- |
| vtx | `g4SimHits` | `std::vector<SimVertex>` | general event info |
