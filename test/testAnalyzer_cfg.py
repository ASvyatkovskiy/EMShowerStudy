import FWCore.ParameterSet.Config as cms

# set up process
process = cms.Process("Showers")

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *

# Load geometry
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('MC_38Y_V8::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("FWCore.MessageService.MessageLogger_cfi")

# Muon Reco
process.load("RecoLocalMuon.Configuration.RecoLocalMuon_cff")
process.load("RecoMuon.Configuration.RecoMuon_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.ttrhbwor.ComputeCoarseLocalPositionFromDisk = True
process.ttrhbwr.ComputeCoarseLocalPositionFromDisk = True

process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.destinations += ['DetailedMessages']
process.MessageLogger.categories   = ['MuonShowerFinder', 'GlobalMuonRefitter']
process.MessageLogger.debugModules += ['globalMuons','tevMuons']


process.MessageLogger.DetailedMessages = cms.untracked.PSet(
    threshold  = cms.untracked.string('DEBUG'),
    default    = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    INFO = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    MuonShowerFinder = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
    GlobalMuonRefitter = cms.untracked.PSet(limit = cms.untracked.int32(-1)),

    )

# Muon Reco
process.load("RecoLocalMuon.Configuration.RecoLocalMuon_cff")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('test.root')
)

process.showers = cms.EDAnalyzer("EMShowerAnalyzer",
      MuonServiceProxy
     , TrackLabel = cms.InputTag('globalMuons')
     , debug = cms.bool(True),
    TrackerRecHitBuilder = cms.string('WithTrackAngle'),
    MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),


    EMShowerParameters = cms.PSet(
    DTRecSegmentLabel = cms.InputTag("dt1DRecHits"),
    CSCRecSegmentLabel = cms.InputTag("csc2DRecHits"),
    RPCRecSegmentLabel = cms.InputTag("rpcRecHits"),
    DT4DRecSegmentLabel = cms.InputTag("dt4DSegments"),
    CSCSegmentLabel = cms.InputTag("cscSegments"),

    TrackerRecHitBuilder = cms.string('WithTrackAngle'),
    MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
    sizeThreshold1 = cms.double(20),
    hitThreshold1 = cms.int32(40),
    sizeThreshold2 = cms.double(15),
    hitThreshold2 = cms.int32(60),
    
 )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source(
     "PoolSource",
     noEventSort = cms.untracked.bool(True),
     duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
     skipEvents = cms.untracked.uint32(0),
     fileNames = cms.untracked.vstring(
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_9.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_8.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_7.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_6.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_50.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_5.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_49.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_47.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_46.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_45.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_44.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_43.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_42.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_41.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_40.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_4.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_39.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_38.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_37.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_36.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_35.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_34.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_33.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_32.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_31.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_30.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_3.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_29.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_28.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_27.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_26.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_25.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_24.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_23.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_22.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_21.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_2.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_19.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_18.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_17.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_16.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_15.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_13.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_12.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_11.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_10.root',
       '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt1000/aeverett/SingleMuPt1000_CMSSW_3_4_1_step1/SingleMuPt1000_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773/step2_1.root'
     ),
     secondaryFileNames = cms.untracked.vstring(
    
#    '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt500/aeverett/SingleMuPt500_CMSSW_3_4_1_step1/SingleMuPt500_CMSSW_3_4_1_step1/69f84e5ef963b0cf6410679d5122df01/SingleMuPt500_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco_9.root'
     )
)

process.p = cms.Path(process.muonrecoComplete*process.showers)
