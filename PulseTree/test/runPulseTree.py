from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register ('runtype',
                  '', # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.string,          # string, int, or float
                  "Type of run to analyze [TP,PhiSym,Laser,Zee]")
options.parseArguments()

process = cms.Process('PulseTree')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string(options.outputFile)
)

process.options = cms.untracked.PSet(
#    SkipEvent = cms.untracked.vstring('ProductNotFound'),
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('reco nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v12', '')

import EventFilter.EcalRawToDigi.EcalUnpackerData_cfi
process.ecalDigis = EventFilter.EcalRawToDigi.EcalUnpackerData_cfi.ecalEBunpacker.clone()
process.ecalDigis.DoRegional = False

process.load('PulseStudies.PulseTree.PulseTree_cfi')
use_raw_dat = False
make_collections = ['digis']

if options.runtype=='TP':
    process.PulseTree.doAverage = cms.untracked.bool(True)
    use_raw_dat = True
    process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(50)
elif options.runtype=='PhiSym':
    process.PulseTree.doAverage = cms.untracked.bool(True)
    process.PulseTree.splitByLumi = cms.untracked.bool(True)
    process.PulseTree.filterBx = cms.untracked.vuint32(60,199,338,477,616,803,942,1093,1232,1371,1510,1697,1836,1987,2126,2265,2404,2591,2730,2869,3008,3147,3286) # for run 282184, leading BX in trains
    process.PulseTree.EBDigiCollection = cms.InputTag("hltEcalPhiSymFilter","phiSymEcalDigisEB")
    process.PulseTree.EEDigiCollection = cms.InputTag("hltEcalPhiSymFilter","phiSymEcalDigisEE")
    make_collections.remove('digis')
elif options.runtype=='Laser':
    process.PulseTree.nPedestalSamples = cms.untracked.uint32(2)
    process.PulseTree.minAmplitudeForAverage = cms.untracked.double(20.)
    use_raw_dat = True
elif options.runtype=='Zee':
    process.PulseTree.saveRecHits = cms.untracked.bool(True)
    process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(50)
    make_collections.remove('digis')
    make_collections.append('rechits')
    process.PulseTree.minRecHitAmplitudeWeights = cms.untracked.double(2000.)

if use_raw_dat:
    process.source = cms.Source("NewEventStreamFileReader", fileNames = cms.untracked.vstring(options.inputFiles))
else:
    process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(options.inputFiles),
                                secondaryFileNames = cms.untracked.vstring()
                                )


process.load('RecoLocalCalo.EcalRecProducers.ecalUncalibRecHit_cfi')
process.load('RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi')
process.ecalUncalibRecHit.EBdigiCollection = process.PulseTree.EBDigiCollection
process.ecalUncalibRecHit.EEdigiCollection = process.PulseTree.EEDigiCollection
process.ecalMultiFitUncalibRecHit.EBdigiCollection = process.PulseTree.EBDigiCollection
process.ecalMultiFitUncalibRecHit.EEdigiCollection = process.PulseTree.EEDigiCollection
process.ecalMultiFitUncalibRecHit.algoPSet.useLumiInfoRunHeader = cms.bool(False)

process.ecalMultiFitUncalibRecHit2 = process.ecalMultiFitUncalibRecHit.clone()
process.ecalMultiFitUncalibRecHit2.algoPSet.gainSwitchFix = cms.bool(True)

process.ecalDigis_step = cms.Path(process.ecalDigis)
process.multifit = cms.Path(process.ecalMultiFitUncalibRecHit * process.ecalMultiFitUncalibRecHit2)
process.weights = cms.Path(process.ecalUncalibRecHit)

process.PulseTree_step = cms.Path(process.PulseTree)
process.endjob_step = cms.EndPath(process.endOfProcess)

process.dump = cms.EDAnalyzer("EventContentAnalyzer")
process.dump_step = cms.Path(process.dump)

if 'digis' in make_collections and 'rechits' in make_collections:
    process.schedule = cms.Schedule(process.ecalDigis_step,process.multifit,process.weights,process.PulseTree_step,process.endjob_step)
elif 'digis' in make_collections:
    process.schedule = cms.Schedule(process.ecalDigis_step,process.PulseTree_step,process.endjob_step)
elif 'rechits' in make_collections:
    process.schedule = cms.Schedule(process.multifit,process.weights,process.PulseTree_step,process.endjob_step)
else:
    process.schedule = cms.Schedule(process.PulseTree_step,process.endjob_step)

