#!/usr/bin/env python

import ROOT
from art.morisot import Morisot
from array import array
import sys,os
import argparse

def main():
  # User controlled arguments
  parser = argparse.ArgumentParser()
  parser.add_argument("--inFileName", type=str, default="", help="The path to the input file from runBumpHunter")
  parser.add_argument("--outPath", type=str, default="", help="The path prefix (directory) where you want the output plots")
  parser.add_argument("--lumi", type=float, default=1, help="Luminosity")
  parser.add_argument("--signalFileName", type=str, default="", help="Signal histogram overlaid on background")

  args = parser.parse_args()
  inFileName = args.inFileName
  outPath = args.outPath
  luminosity = args.lumi*1000
  signalFileName = args.signalFileName
    
  print "==================================="
  print "inputFile       : ",inFileName
  print "outPath       : ",outPath
  print "Luminosity       : ",args.lumi
  print "Signal File       : ",signalFileName
  print "==================================="

  # Get input root file
  inFile = ROOT.TFile.Open(inFileName, "READ")
  if not inFile:
    print inFileName, " doesn't exist."
    return
  # make plots folder i.e. make folder extension
  if not os.path.exists(outPath):
      os.makedirs(outPath)

  # Define necessary quantities.
  Ecm = 13

  # Initialize painter
  myPainter = Morisot()
  myPainter.setColourPalette("Teals")
  #myPainter.setEPS(True)
  myPainter.setLabelType(2) # Sets label type i.e. Internal, Work in progress etc.
  # 0 Just ATLAS    
  # 1 "Preliminary"
  # 2 "Internal"
  # 3 "Simulation Preliminary"
  # 4 "Simulation Internal"
  # 5 "Simulation"
  # 6 "Work in Progress"

  # Retrieve search phase inputs
  basicData = inFile.Get("basicData")
  basicBkg = inFile.Get("basicBkg")
  residualHist = inFile.Get("residualHist")
  logLikelihoodPseudoStatHist = inFile.Get("logLikelihoodStatHistNullCase")
  chi2PseudoStatHist = inFile.Get("chi2StatHistNullCase")
  bumpHunterStatHist = inFile.Get("bumpHunterStatHistNullCase")
  bumpHunterTomographyPlot = inFile.Get('bumpHunterTomographyFromPseudoexperiments')
  bumpHunterStatOfFitToData = inFile.Get('bumpHunterStatOfFitToData')

  logLOfFitToDataVec = inFile.Get('logLOfFitToData')
  chi2OfFitToDataVec = inFile.Get('chi2OfFitToData')
  statOfFitToData = inFile.Get('bumpHunterPLowHigh')
  logLOfFitToData = logLOfFitToDataVec[0]
  logLPVal = logLOfFitToDataVec[1]
  chi2OfFitToData = chi2OfFitToDataVec[0]
  chi2PVal = chi2OfFitToDataVec[1]
  bumpHunterStatFitToData = statOfFitToData[0]
  bumpHunterPVal = bumpHunterStatOfFitToData[1]
  bumpLowEdge = statOfFitToData[1]
  bumpHighEdge = statOfFitToData[2]

  print "logL of fit to data is",logLOfFitToData
  print "logL pvalue is",logLPVal
  print "chi2 of fit to data is",chi2OfFitToData
  print "chi2 pvalue is",chi2PVal
  print "bump hunter stat of fit to data is",bumpHunterStatFitToData
  print "bumpLowEdge, bumpHighEdge are",bumpLowEdge,bumpHighEdge
  print "BumpHunter pvalue is",bumpHunterPVal
  print "which is Z value of",GetZVal(bumpHunterPVal,True)

  # Find range
  # Calculate from fit range
  fitRange = inFile.Get("FitRange")
  firstBin = basicData.FindBin(fitRange[0])-1
  lastBin = basicData.FindBin(fitRange[1])+2
  print "firstbin, lastbin: ",firstBin,lastBin
  print "First bin = ",firstBin,": lower edge at",basicData.GetBinLowEdge(firstBin)
  print "Last bin = ",lastBin,": higher edge at" ,basicData.GetBinLowEdge(lastBin)+basicData.GetBinWidth(lastBin)

  # Convert plots into desired final form
  standardbins = basicData.GetXaxis().GetXbins()
  newbins = []#ROOT.TArrayD(standardbins.GetSize())
  for np in range(standardbins.GetSize()) :
    newbins.append(standardbins[np]/1000)

  # Make never versions of old plots
  newbasicdata = ROOT.TH1D("basicData_TeV","basicData_TeV",len(newbins)-1,array('d',newbins))
  newbasicBkg = ROOT.TH1D("basicBkg_TeV","basicBkg_TeV",len(newbins)-1,array('d',newbins))
  newresidualHist = ROOT.TH1D("residualHist_TeV","residualHist_TeV",len(newbins)-1,array('d',newbins))

  for histnew,histold in [[newbasicdata,basicData],[newbasicBkg,basicBkg], [newresidualHist,residualHist]]:
    for bin in range(histnew.GetNbinsX()+2) :
      histnew.SetBinContent(bin,histold.GetBinContent(bin))
      histnew.SetBinError(bin,histold.GetBinError(bin))
 
  # Significances for Todd
  ToddSignificancesHist = ROOT.TH1D("ToddSignificancesHist","ToddSignificancesHist",100,-5,5)
  for bin in range(0,newresidualHist.GetNbinsX()+1):
    if bin < firstBin: continue
    if bin > lastBin: continue
    residualValue = newresidualHist.GetBinContent(bin)
    #print residualValue
    ToddSignificancesHist.Fill(residualValue)
  myPainter.drawBasicHistogram(ToddSignificancesHist,-1,-1,"Residuals","Entries","{0}/ToddSignificancesHist".format(outPath))

  # Search phase plots
  myPainter.drawDataAndFitOverSignificanceHist(newbasicdata,newbasicBkg,newresidualHist,\
            'm_{jj} [TeV]','Events','Significance','{0}/figure1'.format(outPath),\
            luminosity,13,fitRange[0],fitRange[1],firstBin,lastBin+2,True,\
            bumpLowEdge/1000.0,bumpHighEdge/1000.0,[],True,False,[],True,bumpHunterPVal)
  myPainter.drawDataAndFitOverSignificanceHist(newbasicdata,newbasicBkg,newresidualHist,\
            'm_{jj} [TeV]','Events','Significance','{0}/figure1_nobump'.format(outPath),\
            luminosity,13,fitRange[0],fitRange[1],firstBin,lastBin+2,False,\
            bumpLowEdge,bumpHighEdge,[],True,False,[],True,bumpHunterPVal)
  myPainter.drawPseudoExperimentsWithObservedStat(logLikelihoodPseudoStatHist,\
                              float(logLOfFitToData),logLPVal,0,luminosity,13,\
            'logL statistic','Pseudo-exeperiments',"{0}/logLStatPlot".format(outPath))
  myPainter.drawPseudoExperimentsWithObservedStat(chi2PseudoStatHist,
                              float(chi2OfFitToData),chi2PVal,0,luminosity,13,\
            "#chi^{2}",'Pseudo-exeperiments',"{0}/chi2StatPlot".format(outPath))
  myPainter.drawPseudoExperimentsWithObservedStat(bumpHunterStatHist,
                          float(bumpHunterStatFitToData),bumpHunterPVal,0,luminosity,13,\
            'BumpHunter','Pseudo-exeperiments',"{0}/bumpHunterStatPlot".format(outPath))
  myPainter.drawBumpHunterTomographyPlot(bumpHunterTomographyPlot,"{0}/bumpHunterTomographyPlot".format(outPath))


  ####### Draw Signal Template overlaid on Background###########
  signalFile = ROOT.TFile.Open(signalFileName, "read")
  if not signalFile:
    print signalFileName, " doesn't exist!!!!"
    return

  # setup signal information
  signalTitles = {"QStar": "#it{q}*"}
  signalTypes = ["QStar"]
  signalsMasses = {"QStar":[4000, 5000]}
  signalScalingFactors = {"QStar": 10}
  signalAxes = {"QStar": {"X" : "M_{#it{q}*} [GeV]", "Y": "#sigma #times #it{A} #times BR [pb]"} }

  for signalType in signalTypes:
    print "in signal",signalType

    signalMasses = signalsMasses[signalType]
    signalMassesTeV = signalsMasses[signalType][:]
    for index in range(len(signalMasses)) :
      signalMassesTeV[index] = signalMasses[index]/1000.0
    print signalMassesTeV

    signalPlotsTeV = []
    legendlistTeV = []
    for mass in signalMasses :
      sigplot = signalFile.Get("mjj_Scaled_"+signalType+"{0}_10fb".format(mass))
      sigplot.SetDirectory(0)  

      sigplottev = newbasicdata.Clone()
      sigplottev.SetName("sigplot_{0}_{1}_TeV".format(signalType,mass))
      for bin in range(sigplottev.GetNbinsX()+2) :
        sigplottev.SetBinContent(bin,sigplot.GetBinContent(bin))
        sigplottev.SetBinError(bin,sigplot.GetBinError(bin))

      #luminosity = float(luminosity)
      sigplottev.Scale(luminosity/10000.)
      sigplotforfitplusbkg = sigplottev.Clone()
      sigplotforfitplusbkg.SetDirectory(0)
      sigplotforfitplusbkg.SetName(sigplottev.GetName()+"_forfitplusbkg_TeV")
      sigplotforfitplusbkg.Scale(signalScalingFactors[signalType])
      signalPlotsTeV.append(sigplotforfitplusbkg)
      thistitle = signalTitles[signalType] + ", {0}= {1} TeV".format(signalAxes[signalType]["X"].split("[GeV]")[0].replace("M","m"),mass/1000.0)
      legendlistTeV.append(thistitle)

      extLastBin = lastBin
      for bin in range(sigplotforfitplusbkg.GetNbinsX()+2) :
        if bin > extLastBin and sigplotforfitplusbkg.GetBinContent(bin) > 0.01 :
          extLastBin = bin
        if sigplotforfitplusbkg.GetBinLowEdge(bin) > 1.3*mass/1000.0 :
          continue
        if extLastBin < lastBin :
          extLastBin = lastBin

    UserScaleText = signalTitles[signalType]
    if signalScalingFactors[signalType] == 1 :
      UserScaleText = signalTitles[signalType]
    else :
      UserScaleText = UserScaleText+",  #sigma #times "+str(signalScalingFactors[signalType])

    outputName = outPath+"FancyFigure1_"+signalType
    myPainter.drawDataAndFitWithSignalsOverSignificances(newbasicdata,newbasicBkg, None,\
                   newresidualHist,signalPlotsTeV, None, signalMassesTeV,legendlistTeV,\
                   "m_{jj} [TeV]","Events","","Significance ", outputName,luminosity,\
                   Ecm, firstBin,extLastBin+5,\
                   True, bumpLowEdge/1000,bumpHighEdge/1000,\
                   True,False,True,UserScaleText,True,bumpHunterPVal, True, \
                   fitRange[0], fitRange[1])
    outputName = outPath+"FancyFigure1_"+signalType+"_noBump"
    myPainter.drawDataAndFitWithSignalsOverSignificances(newbasicdata,newbasicBkg, None,\
                   newresidualHist,signalPlotsTeV, None, signalMassesTeV,legendlistTeV,\
                   "m_{jj} [TeV]","Events","","Significance ", outputName,luminosity,\
                   Ecm, firstBin,extLastBin+5,\
                   False, bumpLowEdge/1000,bumpHighEdge/1000,\
                   True,False,True,UserScaleText,True,bumpHunterPVal, True, \
                   fitRange[0], fitRange[1])

  inFile.Close()
  print "Done."

def GetZVal (p, excess) :
  #the function normal_quantile converts a p-value into a significance,
  #i.e. the number of standard deviations corresponding to the right-tail of 
  #a Gaussian
  if excess :
    zval = ROOT.Math.normal_quantile(1-p,1);
  else :
    zval = ROOT.Math.normal_quantile(p,1);

  return zval

def MakeHistoFromStats(statistics) :

  nentries = len(statistics)
  nBins = int(float(nentries)/10.0)

  maxVal = max(statistics)
  minVal = min(statistics)
  axisrange = maxVal - minVal;

  thismin = minVal-0.05*axisrange;
  thismax = maxVal+0.05*axisrange;

  statPlot = ROOT.TH1D("statPlot","",nBins,thismin,thismax)
  for val in range(len(statistics)) :
    statPlot.Fill(statistics[val])

  return statPlot

# when calling this script
if __name__ == "__main__":
    main()
