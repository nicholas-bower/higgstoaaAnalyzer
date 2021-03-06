import subprocess
import sys,string,math,os
import ConfigParser
import glob
import numpy as np


filesPerList=5

SampleName = sys.argv [1]
def checkAndMakeDir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

def clearDir(dir):
    for fil in glob.glob(dir+"/*"):
        os.remove(fil)

if "QCD" in SampleName:
    masses = ['50to100', '100to200', '200to300','300to500','500to700','700to1000','1000to1500','1500to2000','2000toInf']
else:
    masses = ['FULL_HT']
for Mass in masses:
    first_step=True

    fileListDir="/uscms_data/d3/nbower/FSU/HtoAA/Analyzer/CMSSW_10_6_26/src/higgstoaaAnalyzer/higgstoaaAnalyzer/fileLists/"+SampleName+"/"+Mass+"/"
    checkAndMakeDir("/uscms_data/d3/nbower/FSU/HtoAA/Analyzer/CMSSW_10_6_26/src/higgstoaaAnalyzer/higgstoaaAnalyzer/fileLists/"+SampleName+"/")

    if first_step==True:
        print(fileListDir)
        checkAndMakeDir(fileListDir)
    rootFileDir="/eos/uscms/store/user/nbower/Events/"+SampleName+"/"+Mass+"/"

    first_step=False

    query = "ls " + rootFileDir
    os.system(query)
    files=os.popen(query).read().split()
    for nf in range(1, len(files)+1):
        filelistIdx=int((nf-1)/filesPerList)

        if nf%filesPerList==1:
            out=open(fileListDir+SampleName+"_"+str(filelistIdx)+".txt","w")
        out.write("root://cmseos.fnal.gov/"+rootFileDir+"/"+files[nf-1]+"\n")
    print (str(len(files)))
