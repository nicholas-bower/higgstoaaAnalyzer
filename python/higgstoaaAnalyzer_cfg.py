# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: test2 -s RAW2DIGI,L1Reco,RECO --datatier RECO --era=Run2_2018 --conditions auto:phase1_2018_realistic --eventcontent RECO --filein file:test.root --no_exec
import FWCore.ParameterSet.Config as cms
import os
location = os.getcwd() 
query = "ls *.root"
infiles=os.popen(query).read().split()
fileVec=[]
for line in infiles:
    fileVec.append("file:"+location+"/"+line)
    print "file:"+location+"/"+line
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.setDefault('maxEvents',-1)
options.setDefault('inputFiles',fileVec)


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
process.load('higgstoaaAnalyzer/higgstoaaAnalyzer/higgstoaaAnalyzer_cfi')



# Schedule definition

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string(location+'/Outfile.root')
                                   )


process.p = cms.Path(process.simple)
# Customisation from command line


# End adding early deletion

#open('pydump.py','w').write(process.dumpPython())
