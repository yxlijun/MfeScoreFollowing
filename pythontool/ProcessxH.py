# -*- coding:utf-8 -*-
from __future__ import division     
import numpy as np        

xHMinThresh = 240
maxPitchesInEvent = 5
minDurFrame = 7

class ProcessxH(object):
    def __init__(self,minDurFrame=7):
        super(ProcessxH, self).__init__()
        self.minDurFrame = minDurFrame
        self.nFrameCount = np.zeros(88,dtype=np.int)
        self.xHpeak = np.zeros(88)
        self.xHPre7Peak = np.zeros(88)

        self.sameOnset = list()
        self.newPitches = list()
        self.timeHPeakPair = [[] for i in range(88)]

        self.timePitchesPair = list()

    def ProcessingFrame(self,xH,iFrame,timeResolution):
        self.UpdateFrameCount(xH,iFrame,timeResolution)
        self.GenerateNewPitches(xH,iFrame)

    def UpdateFrameCount(self,xH,iFrame,timeResolution):
        maxH = np.max(xH)
        pianoSum = np.sum(np.greater(xH,xHMinThresh))
        pianoRoll = np.zeros(len(xH))
        if pianoSum>maxPitchesInEvent:
            temp = xH.copy()
            index = np.argsort(-temp)[0:maxPitchesInEvent]
            pianoRoll[index] = 1
        else:
            pianoRoll[np.greater(xH,xHMinThresh)] = 1
        
        for i,ncount in enumerate(self.nFrameCount):
            if pianoRoll[i]==1:
                self.nFrameCount[i]+=1
                if self.nFrameCount[i]==self.minDurFrame:
                    onset = (iFrame-self.minDurFrame+1)*timeResolution
                    timePeak = [onset,self.xHPre7Peak[i]]
                    self.timeHPeakPair[i].append(timePeak)
                if self.nFrameCount[i]>self.minDurFrame:
                    if self.xHpeak[i]<xH[i]:
                        self.xHpeak[i] = xH[i]
                else:
                    if self.xHPre7Peak[i]<=xH[i]:
                        self.xHPre7Peak[i] = xH[i]
            else:
                if self.nFrameCount[i]>self.minDurFrame:
                    offset = iFrame*timeResolution         
                    self.timeHPeakPair[i][-1].append(offset)
                    self.timeHPeakPair[i][-1].append(self.xHpeak[i])

                    self.xHpeak[i] = 0
                    self.xHPre7Peak[i] = 0

                self.nFrameCount[i]  = 0


    def GenerateNewPitches(self,xH,iFrame):
        thisNewPitches = list()
        if not self.sameOnset or (self.sameOnset and self.sameOnset[0]==iFrame):
            newPitchesIndex = np.argwhere(self.nFrameCount==self.minDurFrame)
            for index in newPitchesIndex:
                thisNewPitches.append(index[0]+1)
        else:
            thisNewPitches = []
        
        if not self.sameOnset:
            if thisNewPitches:
                newPitchesCandidate = list()
                candidate = np.where((self.nFrameCount>0) & (self.nFrameCount<self.minDurFrame))[0]
                for index in candidate:
                    newPitchesCandidate.append(index+1)
                if newPitchesCandidate:
                    onsetTemp = list()
                    for i,pitchcandidate in enumerate(newPitchesCandidate):
                        onsetTemp.append(self.nFrameCount[pitchcandidate-1])
                        onsetTemp[i] = iFrame+self.minDurFrame-onsetTemp[i]
                    onsetTemp = set(onsetTemp)
                    self.sameOnset.extend(onsetTemp)
            self.newPitches = thisNewPitches
        else:
            if self.sameOnset[0] == iFrame:
                if thisNewPitches:
                    self.newPitches.extend(thisNewPitches)
                    self.newPitches = sorted(self.newPitches)
                del self.sameOnset[0]
            if np.where(np.array(self.nFrameCount==1)>1)[0]:
                self.sameOnset.append(iFrame+minDurFrame-1)


        if not self.sameOnset:
            pairTemp = list()
            appearance = np.zeros(88)
            for i,pitch in enumerate(self.newPitches):
                appearance[pitch-1]+=1
                index = int(len(self.timeHPeakPair[pitch-1])-appearance[pitch-1])
                onset = self.timeHPeakPair[pitch-1][index][0]
                pairTemp.append(onset)
                pairTemp.append(pitch)
            if pairTemp:
                self.timePitchesPair.append(pairTemp)

         
    def UpdateEventFlag(self,iFrame):
        pitches = list()
        flag = (self.newPitches and not self.sameOnset) or (iFrame==-1)
        if flag:
            pitches = self.newPitches
        return pitches


    def GetTimePitchesPair(self):
        return self.timePitchesPair


def readH(Hpath):
    H = list()
    with open(Hpath) as f:
        data = f.readlines()
        for line in data:
            line = line.strip().split('\t')
            H.append(np.array(line,dtype=np.float32))            
    return H



if __name__ =='__main__':
    Hpath = './H.txt'
    H = readH(Hpath)
    timeResolution = 512/44100
    processh = ProcessxH()
    for i,xH in enumerate(H):
        processh.ProcessingFrame(xH,i,timeResolution)
        pitches = processh.UpdateEventFlag(i)
        if pitches:
            print pitches
    timePitches = processh.GetTimePitchesPair()
    print timePitches