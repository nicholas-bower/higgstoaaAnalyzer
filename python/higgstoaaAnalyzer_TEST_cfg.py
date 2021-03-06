# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: test2 -s RAW2DIGI,L1Reco,RECO --datatier RECO --era=Run2_2018 --conditions auto:phase1_2018_realistic --eventcontent RECO --filein file:test.root --no_exec
import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
import os
location = os.getcwd() 
options = VarParsing ('python')
options.setDefault('maxEvents',50)
options.setDefault('inputFiles',"root://cmseos.fnal.gov//store/user/nbower/Events/SUSYGluGluToHToAA_AToBB_AToTauTau_M-12_FilterTauTauTrigger_TuneCP5_13TeV_madgraph_pythia8_AOD/SUSYGluGluToHToAA_AToBB_AToTauTau_M-12_FilterTauTauTrigger_TuneCP5_13TeV_madgraph_pythia8_AOD_1.root")


options.parseArguments()



process = cms.Process('TEST') # ,eras.bParkingOpen

# import of standard configurations


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),  
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)



# Output definition


# Path and EndPath definitions


# Schedule definition

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string(location+'/testOutput.root')
                                   )

process.load('higgstoaaAnalyzer/higgstoaaAnalyzer/higgstoaaAnalyzer_cfi')

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
                                   src = cms.InputTag("genParticles"),                                                                 
                                   printP4 = cms.untracked.bool(False),
                                   printPtEtaPhi = cms.untracked.bool(False),
                                   printVertex = cms.untracked.bool(False),
                                   printStatus = cms.untracked.bool(True),
                                   printIndex = cms.untracked.bool(False),
                                   status = cms.untracked.vint32()
                                   )
#process.p = cms.Path(process.printTree)
process.p = cms.Path(process.simple)
# Customisation from command line


# End adding early deletion

#open('pydump.py','w').write(process.dumpPython())
