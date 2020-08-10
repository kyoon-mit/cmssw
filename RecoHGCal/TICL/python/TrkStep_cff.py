import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.TICLSeedingRegions_cff import ticlSeedingTrk
from RecoHGCal.TICL.ticlLayerTileProducer_cfi import ticlLayerTileProducer as _ticlLayerTileProducer
from RecoHGCal.TICL.trackstersProducer_cfi import trackstersProducer as _trackstersProducer
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer as _filteredLayerClustersProducer
from RecoHGCal.TICL.multiClustersFromTrackstersProducer_cfi import multiClustersFromTrackstersProducer as _multiClustersFromTrackstersProducer

# CLUSTER FILTERING/MASKING

filteredLayerClustersTrk = _filteredLayerClustersProducer.clone(
  clusterFilter = "ClusterFilterByAlgo",
  algo_number = 8,
  iteration_label = "TRK"
)

# CA - PATTERN RECOGNITION

ticlTrackstersTrk = _trackstersProducer.clone(
  filtered_mask = cms.InputTag("filteredLayerClustersTrk", "Trk"),
  seeding_regions = "ticlSeedingTrk",
  filter_on_categories = [2, 4], # filter muons and charged hadrons
  pid_threshold = 0.0,
  missing_layers = 3,
  min_clusters_per_ntuplet = 10,
  min_cos_theta = 0.866, # ~30 degrees
  min_cos_pointing = 0.798, # ~ 37 degrees
  max_delta_time = -1.,
  algo_verbosity = 2,
  oneTracksterPerTrackSeed = True,
  promoteEmptyRegionToTrackster = True,
  itername = "TRK"
)

# MULTICLUSTERS

ticlMultiClustersFromTrackstersTrk = _multiClustersFromTrackstersProducer.clone(
    Tracksters = "ticlTrackstersTrk"
)

ticlTrkStepTask = cms.Task(ticlSeedingTrk
    ,filteredLayerClustersTrk
    ,ticlTrackstersTrk
    ,ticlMultiClustersFromTrackstersTrk)
    
    
# HFNose

filteredLayerClustersTrkHFNose = filteredLayerClustersTrk.clone(
    LayerClusters = 'hgcalLayerClustersHFNose',
    LayerClustersInputMask = cms.InputTag("hgcalLayerClustersHFNose","InitialLayerClustersMask"),
    iteration_label = "TRKn",
    algo_number = 9
#no tracking mask for EM for now
)

ticlTrackstersTrkHFNose = ticlTrackstersTrk.clone(
    detector = "HFNose",
    layer_clusters = "hgcalLayerClustersHFNose",
    layer_clusters_hfnose_tiles = "ticlLayerTileHFNose",
    original_mask = cms.InputTag("hgcalLayerClustersHFNose","InitialLayerClustersMask"),
    filtered_mask = cms.InputTag("filteredLayerClustersHFNoseEM","EMn"),
    seeding_regions = "ticlSeedingTrkHFNose",
    time_layerclusters = cms.InputTag("hgcalLayerClustersHFNose","timeLayerCluster"),
    min_clusters_per_ntuplet = 6,
    itername = "TRKn"
)

ticlHFNoseTrkStepTask = cms.Task(ticlSeedingGlobalHFNose
                              ,filteredLayerClustersTrkHFNose
                              ,ticlTrackstersTrkHFNose
)

