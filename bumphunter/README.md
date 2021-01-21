## Author: DengfengZhang(dengdeng.zhang@cern.ch)

## BumpHunter: Search for local excess.

## Package Architecture:
### inc/ Header Files:
### src/ Soure Files


## Run BumpHunter?:
1. Setup ROOT
2. Compile: make
3. Run:
  ./runBumpHunter --inFile inputs/BkgPlusSignalScale100_Seed2000.root --outPath results/GenericZX/ --outFileName BumpHunter_BkgPlusSignalScale10_Zmm_CR_LeadFatJ_ZXmass.root --dataHist Zmm_CR_LeadFatJ_ZXmass --bkgHist bkg_Zmm_CR_LeadFatJ_ZXmass --nPseudoExpBH 1000

## Make BumpHunter Plots:
1. Nivagate to scripts/
2. Setup enviromnet: source setup.sh
3. Run:
* Draw BumpHunter Out
> python scriptResonance/BumpHunter/plotBumpHunter.py --inFileName input/BumpHunter_Resonance.root --outPath plotting/BumpHunter/plots_Resonance/ --lumi 140
* Signal overlaid on Background
> python scriptResonance/BumpHunter/plotBumpHunter.py --inFileName input/BumpHunter_Resonance.root --outPath plotting/BumpHunter/Test/ --lumi 140 --overlaidSignal --signalFileName input/dijetMC_qstar_fullBins.root
* WithMCRatio
> python scriptResonance/BumpHunter/plotBumpHunter.py --inFileName input/BumpHunter_Resonance.root --outPath plotting/BumpHunter/Test/ --lumi 140 --overlaidSignal --signalFileName input/dijetMC_qstar_fullBins.root --mcFileName input/pseudoMC.root --drawMCComparison 

