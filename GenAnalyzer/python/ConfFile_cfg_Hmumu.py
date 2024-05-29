import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

# Options and Output Report
process.options = cms.untracked.PSet(
                                     wantSummary = cms.untracked.bool(False)
                                     )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#import FWCore.Utilities.FileUtils as FileUtils
inputTxtFile="HIG-Run3Summer22EEwmLHEGS-01111_k"
#readFiles = cms.untracked.vstring(FileUtils.loadListFromFile('inputTxtFiles/'+inputTxtFile+'.txt'))


process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            #fileNames = readFiles,
                            fileNames = cms.untracked.vstring(
                            'file:/eos/user/q/qguo/vbfhmm/production/HIG-Run3Summer22EEwmLHEGS-01111_test.root',
                            'file:/eos/user/q/qguo/vbfhmm/production/HIG-Run3Summer22EEwmLHEGS-01111_1.root',
                            'file:/eos/user/q/qguo/vbfhmm/production/HIG-Run3Summer22EEwmLHEGS-01111_10.root',
                            'file:/eos/user/q/qguo/vbfhmm/production/HIG-Run3Summer22EEwmLHEGS-01111_11.root',
                            'file:/eos/user/q/qguo/vbfhmm/production/HIG-Run3Summer22EEwmLHEGS-01111_12.root',
                            'file:/eos/user/q/qguo/vbfhmm/production/HIG-Run3Summer22EEwmLHEGS-01111_13.root',
                            'file:/eos/user/q/qguo/vbfhmm/production/HIG-Run3Summer22EEwmLHEGS-01111_14.root',
                            'file:/eos/user/q/qguo/vbfhmm/production/HIG-Run3Summer22EEwmLHEGS-01111_15.root',
                            'file:/eos/user/q/qguo/vbfhmm/production/HIG-Run3Summer22EEwmLHEGS-01111_16.root',
                            'file:/eos/user/q/qguo/vbfhmm/production/HIG-Run3Summer22EEwmLHEGS-01111_17.root',
                            'file:/eos/user/q/qguo/vbfhmm/production/HIG-Run3Summer22EEwmLHEGS-01111_18.root',
                            'file:/eos/user/q/qguo/vbfhmm/production/HIG-Run3Summer22EEwmLHEGS-01111_19.root',
                            'file:/eos/user/q/qguo/vbfhmm/production/HIG-Run3Summer22EEwmLHEGS-01111_20.root',
                            'file:/eos/user/q/qguo/vbfhmm/production/HIG-Run3Summer22EEwmLHEGS-01111_21.root',
                            'file:/eos/user/q/qguo/vbfhmm/production/HIG-Run3Summer22EEwmLHEGS-01111_22.root',
                            ),
                            skipBadFiles = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
                            )

process.demo = cms.EDAnalyzer('GenAnalyzer',
                              Verbose		=	cms.bool(False),
                              #OutPutFileName = cms.string("test.root"),
                              OutPutFileName = cms.string(inputTxtFile+'.root'),
                              #OutPutFileName = cms.string("GF_HH_Benchmark3.root"),
                              #genParticlesInputTag  = cms.InputTag("prunedGenParticles"),	# Uncomment if running on MiniAOD
                              #LHEEventInputTag = cms.InputTag("source"),			# Uncomment if running on MiniAOD
                              genParticlesInputTag  = cms.InputTag("genParticles"),		# Uncomment if running on GEN only sample
                              LHEEventInputTag = cms.InputTag("externalLHEProducer"),		# Uncomment if running on GEN only
                              #ak4GenJetsInputTag = cms.InputTag("ak4GenJets"),
                              #ak8GenJetsInputTag = cms.InputTag("ak8GenJets"),

                              )

#process.GenAnalyzer = cms.EDProducer("GenAnalyzer",
#)

process.p = cms.Path(process.demo)
