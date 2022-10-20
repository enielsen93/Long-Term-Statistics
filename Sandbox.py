LTSFile = r"C:\Dokumenter\Soholt\Alt3aar.MJL"
splitCount = 4
import numpy as np

with open(LTSFile,'r') as f:
    LTSTxt = f.read().split("\n")
LTSTxt = np.array(LTSTxt)

LTSStartLine = [i for i,a in enumerate(LTSTxt) if r"[SIMULATION_EVENT]" in a]
LTSEndLine = [i for i,a in enumerate(LTSTxt) if r"EndSect  // SIMULATION_EVENT" in a]

#LTSNewTxt = LTSTxt[0:LTSStartLine[0]]
jobsPerLTS = int(np.ceil(float(len(LTSStartLine))/4))
LTSNewTxt = [[] for i in range(splitCount)]
job = 0
for i in np.arange(0,len(LTSStartLine),jobsPerLTS):
    LTSNewTxt[job] = LTSTxt[0:LTSStartLine[0]]
    for j in np.arange(i,min(i+jobsPerLTS,len(LTSStartLine))):
        LTSNewTxt[job] = np.concatenate((LTSNewTxt[job],LTSTxt[LTSStartLine[j]:LTSEndLine[j]+1]))
    LTSNewTxt[job] = np.concatenate((LTSNewTxt[job],LTSTxt[LTSEndLine[-1]+1:len(LTSTxt)]))
    
    with open(LTSFile[0:-4] + "_%d.mjl" % job,'w') as f:
        f.write("\n".join(LTSNewTxt[job]))
    job += 1
