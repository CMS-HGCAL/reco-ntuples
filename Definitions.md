# Ntuple content definitions

| branch name | collection name | collection type | definition  |
| ------------- | ------------- | ----- | ----- |
| event/run/lumi | - | - | general event info |
| vtx | `g4SimHits` | `std::vector<SimVertex>` | primary vertex position |
| genpart | `g4SimHits` | `std::vector<SimTrack>` | truth level tracks/particles and information related to their extrapolation towards HGCAL. In particular `reachedEE==2` indicates that the particles reached HGCAL while `reachedEE==1` is for barrel calorimeter and `reachedEE==0` is for the other cases |
| simcluster | `mix:MergedCaloTruth` | `std::vector<SimCluster>` | Geant particle and its associated hits (DetIds) in the HGCal |
| pfcluster | `particleFlowClusterHGCal` | `std::vector<reco::PFCluster>` | mapping of the SimCluster DetIds to the reconstructed hits |
| cluster2d | `hgcalLayerClusters` | `reco::CaloClusterCollection` | reconstructed layer (2D) clusters - those that are associated to a multicluster have `cluster2d_multicluster >= 0`, which is the index of the `multiclus` in the ntuple |
| multiclus | `hgcalLayerClusters` | `reco::HGCalMultiCluster` | reconstructed multi-cluster (3D) built from layer clusters |
| rechit | `HGCalRecHit::HGCEERecHits` <br> `HGCalRecHit::HGCHEFRecHits` <br> `HGCalRecHit::HGCHEBRecHits` | `HGCRecHitCollection` | all reconstructed calorimeter hits - those that are associated to layer clusters have `rechit_cluster2d >= 0`, which is the index of the `cluster2d` in the ntuple |
| calopart | `mix:MergedCaloTruth` | `std::vector<CaloParticle>` | Every CaloParticle is linked to the first stable particle originating from the cascade of particles that left hits in the calorimeters. This stable particle is not included as a SimCluster (unless it itself left hits in the calorimeters). |
| track | `generalTracks` | `std::vector<reco::Track>` | tracks passing highPurity selection |
| gunparticle | `offlinePrimaryVertices` | `reco::Vertex` | `id`, `energy`, `pt`, `eta` and `phi` of gun particles associated to their corresponding vertex |
