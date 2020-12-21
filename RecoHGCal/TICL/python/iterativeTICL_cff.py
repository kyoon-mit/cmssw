import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.MIPStep_cff import *
from RecoHGCal.TICL.TrkEMStep_cff import *
from RecoHGCal.TICL.TrkStep_cff import *
from RecoHGCal.TICL.EMStep_cff import *
from RecoHGCal.TICL.HADStep_cff import *
from RecoHGCal.TICL.ticlLayerTileProducer_cfi import ticlLayerTileProducer
from RecoHGCal.TICL.pfTICLProducer_cfi import pfTICLProducer as _pfTICLProducer
from RecoHGCal.TICL.trackstersMergeProducer_cfi import trackstersMergeProducer as _trackstersMergeProducer
from RecoHGCal.TICL.multiClustersFromTrackstersProducer_cfi import multiClustersFromTrackstersProducer as _multiClustersFromTrackstersProducer

ticlLayerTileTask = cms.Task(ticlLayerTileProducer)

ticlTrackstersMerge = _trackstersMergeProducer.clone()
ticlMultiClustersFromTrackstersMerge = _multiClustersFromTrackstersProducer.clone(
    Tracksters = "ticlTrackstersMerge"
)
ticlTracksterMergeTask = cms.Task(ticlTrackstersMerge, ticlMultiClustersFromTrackstersMerge)


pfTICL = _pfTICLProducer.clone()
ticlPFTask = cms.Task(pfTICL)

iterTICLTask = cms.Task(ticlLayerTileTask
    ,ticlTrkEMStepTask
    ,ticlEMStepTask
    ,ticlTrkStepTask
    ,ticlHADStepTask
    ,ticlTracksterMergeTask
    ,ticlPFTask
    )

ticlHFNoseLayerTile = ticlLayerTileProducer.clone(
    detector = 'HFNose'
)

ticlHFNoseLayerTileTask = cms.Task(ticlHFNoseLayerTile)

iterHFNoseTICLTask = cms.Task(ticlHFNoseLayerTileTask
    ,ticlHFNoseTrkEMStepTask
    ,ticlHFNoseEMStepTask
    ,ticlHFNoseTrkStepTask
)

def injectTICLintoPF(process):
    if getattr(process,'particleFlowTmp', None):
      process.particleFlowTmp.src = ['particleFlowTmpBarrel', 'pfTICL']

    return process
