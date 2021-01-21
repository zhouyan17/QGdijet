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
  parser.add_argument("--overlaidSignal", action='store_true', help="Overlaid Signal on Background")
  parser.add_argument("--signalFileName", type=str, default="", help="Signal histogram overlaid on background")
  parser.add_argument("--drawMCComparison", action='store_true', help="Draw the comparison between data and MC")
  parser.add_argument("--mcFileName", type=str, default="", help="MC File Name")

  args = parser.parse_args()
  inFileName = args.inFileName
  outPath = args.outPath
  luminosity = args.lumi*1000
  overlaidSignal = args.overlaidSignal
  signalFileName = args.signalFileName
  drawMCComparison = args.drawMCComparison
  mcFileName = args.mcFileName
    
  print "==================================="
  print "inputFile       : ",inFileName
  print "outPath       : ",outPath
  print "Luminosity       : ",args.lumi
  print "OverlaidSignal       : ",overlaidSignal
  print "Signal File       : ",signalFileName
  print "drawMCComparison: ", drawMCComparison
  print "mcFileName: ", mcFileName
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
  xAxisLabel = "M_{ZX} [TeV]"
  # Initialize painter
  myPainter = Morisot()
  myPainter.setColourPalette("Teals")
  myPainter.setEPS(False)
  myPainter.setisData(False)
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
  lastBin = basicData.FindBin(fitRange[1])
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
            xAxisLabel,'Events','Significance','{0}/figure1'.format(outPath),\
            luminosity,Ecm,fitRange[0],fitRange[1],firstBin,lastBin,True,\
            bumpLowEdge/1000.0,bumpHighEdge/1000.0,[],True,False,[],True,bumpHunterPVal)
  myPainter.drawDataAndFitOverSignificanceHist(newbasicdata,newbasicBkg,newresidualHist,\
            xAxisLabel,'Events','Significance','{0}/figure1_nologx'.format(outPath),\
            luminosity,Ecm,fitRange[0],fitRange[1],firstBin,lastBin,True,\
            bumpLowEdge/1000.0,bumpHighEdge/1000.0,[],False,False,[],True,bumpHunterPVal)
  myPainter.drawDataAndFitOverSignificanceHist(newbasicdata,newbasicBkg,newresidualHist,\
            xAxisLabel,'Events','Significance','{0}/figure1_nobump'.format(outPath),\
            luminosity,Ecm,fitRange[0],fitRange[1],firstBin,lastBin,False,\
            bumpLowEdge,bumpHighEdge,[],True,False,[],True,bumpHunterPVal)
  myPainter.drawDataAndFitOverSignificanceHist(newbasicdata,newbasicBkg,newresidualHist,\
            xAxisLabel,'Events','Significance','{0}/figure1_nobump_nologx'.format(outPath),\
            luminosity,Ecm,fitRange[0],fitRange[1],firstBin,lastBin,False,\
            bumpLowEdge,bumpHighEdge,[],False,False,[],True,bumpHunterPVal)
  myPainter.drawPseudoExperimentsWithObservedStat(logLikelihoodPseudoStatHist,\
                              float(logLOfFitToData),logLPVal,0,luminosity,Ecm,\
            'logL statistic','Pseudo-exeperiments',"{0}/logLStatPlot".format(outPath))
  myPainter.drawPseudoExperimentsWithObservedStat(chi2PseudoStatHist,
                              float(chi2OfFitToData),chi2PVal,0,luminosity,Ecm,\
            "#chi^{2}",'Pseudo-exeperiments',"{0}/chi2StatPlot".format(outPath))
  myPainter.drawPseudoExperimentsWithObservedStat(bumpHunterStatHist,
                          float(bumpHunterStatFitToData),bumpHunterPVal,0,luminosity,Ecm,\
            'BumpHunter','Pseudo-exeperiments',"{0}/bumpHunterStatPlot".format(outPath))
  myPainter.drawBumpHunterTomographyPlot(bumpHunterTomographyPlot,"{0}/bumpHunterTomographyPlot".format(outPath))


  ####### Draw Signal Template overlaid on Background###########
  if overlaidSignal:
    signalFile = ROOT.TFile.Open(signalFileName, "read")
    if not signalFile:
      print signalFileName, " doesn't exist!!!!"
      return

    # setup signal information
    signalTitles = {"qStar": "#it{q}*"}
    signalTypes = ["qStar"]
    signalsMasses = {"qStar":[4000, 5000]}
    signalScalingFactors = {"qStar": 10}
    signalAxes = {"qStar": {"X" : "M_{#it{q}*} [GeV]", "Y": "#sigma #times #it{A} #times BR [pb]"} }

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
        #sigplot = signalFile.Get("{0}_{1}".format(signalType, mass))
        sigplot = signalFile.Get("h_mjj_{0}".format(mass))
        sigplot.SetDirectory(0)  

        sigplottev = newbasicdata.Clone()
        sigplottev.SetName("sigplot_{0}_{1}_TeV".format(signalType,mass))
        for bins in range(sigplot.GetNbinsX()+2) :
          for bin in range(sigplottev.GetNbinsX()+2) :
            if sigplot.GetBinLowEdge(bins)/1000.==sigplottev.GetBinLowEdge(bin) :
              sigplottev.SetBinContent(bin,sigplot.GetBinContent(bins))
              sigplottev.SetBinError(bin,sigplot.GetBinError(bins))

        sigplotforfitplusbkg = sigplottev.Clone()
        sigplotforfitplusbkg.SetDirectory(0)
        sigplotforfitplusbkg.SetName(sigplottev.GetName()+"_forfitplusbkg_TeV")
        sigplotforfitplusbkg.Scale(signalScalingFactors[signalType])
        signalPlotsTeV.append(sigplotforfitplusbkg)
        massTeV = mass/1000.
        if mass/1000.==int(mass/1000.):
          massTeV = int(massTeV)
        thistitle = signalTitles[signalType] + ", {0}= {1} TeV".format(signalAxes[signalType]["X"].split("[GeV]")[0].replace("M","m"),massTeV)
        legendlistTeV.append(thistitle)

        extLastBin = lastBin
        for bin in range(sigplotforfitplusbkg.GetNbinsX()) :
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
                   xAxisLabel,"Events","","Significance ", outputName,luminosity,\
                   Ecm, firstBin,extLastBin,\
                   True, bumpLowEdge/1000,bumpHighEdge/1000,\
                   True,False,True,UserScaleText,True,bumpHunterPVal, True, \
                   fitRange[0], fitRange[1])
      outputName = outPath+"FancyFigure1_"+signalType+"_nologx"
      myPainter.drawDataAndFitWithSignalsOverSignificances(newbasicdata,newbasicBkg, None,\
                   newresidualHist,signalPlotsTeV, None, signalMassesTeV,legendlistTeV,\
                   xAxisLabel,"Events","","Significance ", outputName,luminosity,\
                   Ecm, firstBin,extLastBin,\
                   True, bumpLowEdge/1000,bumpHighEdge/1000,\
                   False,False,True,UserScaleText,True,bumpHunterPVal, True, \
                   fitRange[0], fitRange[1])
      outputName = outPath+"FancyFigure1_"+signalType+"_noBump"
      myPainter.drawDataAndFitWithSignalsOverSignificances(newbasicdata,newbasicBkg, None,\
                   newresidualHist,signalPlotsTeV, None, signalMassesTeV,legendlistTeV,\
                   xAxisLabel,"Events","","Significance ", outputName,luminosity,\
                   Ecm, firstBin,extLastBin,\
                   False, bumpLowEdge/1000,bumpHighEdge/1000,\
                   True,False,True,UserScaleText,True,bumpHunterPVal, True, \
                   fitRange[0], fitRange[1])
      outputName = outPath+"FancyFigure1_"+signalType+"_noBump_nologx"
      myPainter.drawDataAndFitWithSignalsOverSignificances(newbasicdata,newbasicBkg, None,\
                   newresidualHist,signalPlotsTeV, None, signalMassesTeV,legendlistTeV,\
                   xAxisLabel,"Events","","Significance ", outputName,luminosity,\
                   Ecm, firstBin,extLastBin,\
                   False, bumpLowEdge/1000,bumpHighEdge/1000,\
                   False,False,True,UserScaleText,True,bumpHunterPVal, True, \
                   fitRange[0], fitRange[1])

      ###############################
      # Draw the comparison between data and MC in the bottom panel
      if drawMCComparison:
        mcFile = ROOT.TFile(mcFileName, "read") ;
        if not mcFile:
          print "Can not open: ", mcFileName
          return
        mchist_nominal = mcFile.Get("djet_mjj_nominal")
        mchist_jesup = mcFile.Get("djet_mjj_JES_up")
        mchist_jesdown = mcFile.Get("djet_mjj_JES_down")
        newmchist_nominal=ROOT.TH1D("djet_mjj_nominal_TeV","djet_mjj_nominal_TeV",len(newbins)-1,array('d',newbins))
        newmchist_jesup=ROOT.TH1D("djet_mjj_jesup_TeV","djet_mjj_jesup_TeV",len(newbins)-1,array('d',newbins))
        newmchist_jesdown=ROOT.TH1D("djet_mjj_jesdown_TeV","djet_mjj_jesdown_TeV",len(newbins)-1,array('d',newbins))
        for iBin1 in range(1, newmchist_nominal.GetNbinsX()+1):
          for iBin2 in range(1, mchist_nominal.GetNbinsX()+1):
            if newmchist_nominal.GetBinLowEdge(iBin1)*1000==mchist_nominal.GetBinLowEdge(iBin2):
              newmchist_nominal.SetBinContent(iBin1, mchist_nominal.GetBinContent(iBin2))
              newmchist_nominal.SetBinError(iBin1, mchist_nominal.GetBinError(iBin2))
              continue
        for iBin1 in range(1, newmchist_jesup.GetNbinsX()+1):
          for iBin2 in range(1, mchist_jesup.GetNbinsX()+1):
            if newmchist_jesup.GetBinLowEdge(iBin1)*1000==mchist_jesup.GetBinLowEdge(iBin2):
              newmchist_jesup.SetBinContent(iBin1, mchist_jesup.GetBinContent(iBin2))
              newmchist_jesup.SetBinError(iBin1, mchist_jesup.GetBinError(iBin2))
              continue
        for iBin1 in range(1, newmchist_jesdown.GetNbinsX()+1):
          for iBin2 in range(1, mchist_jesdown.GetNbinsX()+1):
            if newmchist_jesdown.GetBinLowEdge(iBin1)*1000==mchist_jesdown.GetBinLowEdge(iBin2):
              newmchist_jesdown.SetBinContent(iBin1, mchist_jesdown.GetBinContent(iBin2))
              newmchist_jesdown.SetBinError(iBin1, mchist_jesdown.GetBinError(iBin2))
              continue
        tmpRatioHist = newbasicdata.Clone()
        tmpRatioHist.SetMarkerColor(ROOT.kBlack)
        tmpRatioHist.Add(newmchist_nominal,-1)
        tmpRatioHist.Divide(newmchist_nominal)
        ## If data is 0 then there should be no ratio drawn
        for iBin in range(1, tmpRatioHist.GetNbinsX()+1):
          if newbasicdata.GetBinContent(iBin) == 0:
            tmpRatioHist.SetBinContent(iBin, 0)
            tmpRatioHist.SetBinError(iBin, 0)

        UpDownRatioHists = []
        if mchist_jesup.GetEntries() >= 0:
          tmpJESRatioHist = newmchist_jesup
          tmpJESRatioHist.Add( newmchist_nominal, -1. )
          tmpJESRatioHist.Divide( newmchist_nominal )
          tmpJESRatioHist.SetMarkerColorAlpha( ROOT.kBlue,0.15)
          tmpJESRatioHist.SetLineColorAlpha( ROOT.kBlue,0.15)
          tmpJESRatioHist.SetFillColorAlpha( ROOT.kBlue, 0.15)
          tmpJESRatioHist.SetFillStyle(1001)
          UpDownRatioHists.append(tmpJESRatioHist)
        if mchist_jesdown.GetEntries() >= 0:
          tmpJESRatioHist = newmchist_jesdown
          tmpJESRatioHist.Add( newmchist_nominal, -1. )
          tmpJESRatioHist.Divide( newmchist_nominal )
          tmpJESRatioHist.SetMarkerColorAlpha( ROOT.kBlue,0.15)
          tmpJESRatioHist.SetLineColorAlpha( ROOT.kBlue,0.15)
          tmpJESRatioHist.SetFillColorAlpha( ROOT.kBlue, 0.15)
          tmpJESRatioHist.SetFillStyle(1001)
          UpDownRatioHists.append(tmpJESRatioHist)
        outputName = outPath+"FancyFigure1_"+signalType+"_WithMCRatio"
        myPainter.drawDataAndFitWithSignalsOverSignificancesWithMCRatio(newbasicdata,newbasicBkg,None,\
                     newresidualHist, signalPlotsTeV, [], signalMassesTeV,legendlistTeV,\
                     xAxisLabel,"Events","#frac{Data-MC}{MC}","Significance",\
                     outputName,luminosity,Ecm,firstBin,lastBin,True,bumpLowEdge/1000.0,bumpHighEdge/1000.0,\
                     True,False,False, UserScaleText,True,bumpHunterPVal,True,fitRange[0],fitRange[1],\
                     newmchist_nominal,tmpRatioHist,UpDownRatioHists[0],UpDownRatioHists[1])
        outputName = outPath+"FancyFigure1_"+signalType+"_WithMCRatio_nologx"
        myPainter.drawDataAndFitWithSignalsOverSignificancesWithMCRatio(newbasicdata,newbasicBkg,None, \
                     newresidualHist, signalPlotsTeV, [],signalMassesTeV,legendlistTeV, \
                     xAxisLabel,"Events","#frac{Data-MC}{MC}","Significance",\
                     outputName,luminosity,Ecm,firstBin,lastBin,True,bumpLowEdge/1000.0,bumpHighEdge/1000.0,\
                     False,False,False, UserScaleText,True,bumpHunterPVal,True,fitRange[0],fitRange[1],\
                     newmchist_nominal,tmpRatioHist,UpDownRatioHists[0],UpDownRatioHists[1])

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
