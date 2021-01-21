#include "ATLASStyle/AtlasStyle.C"
#include "ATLASStyle/AtlasUtils.C"
#include "ATLASStyle/AtlasLabels.C"

void Roo() {
    SetAtlasStyle ();
    TCanvas *c = new TCanvas("c","c",800,600);
    
    using namespace std;
    using namespace RooFit;
    using namespace RooStats;
    
    TFile *f = TFile::Open("dijetTotal/dijetBkg.root");
    RooWorkspace* combDijetWS = (RooWorkspace*) f->Get("combWS");

    
    RooDataSet *m_data=dynamic_cast<RooDataSet*>(combDijetWS->data("combData"));
    RooStats::ModelConfig *mc = (RooStats::ModelConfig*)combDijetWS->obj( (TString) "ModelConfig" );
    RooSimultaneous *m_pdf = dynamic_cast<RooSimultaneous*>(mc->GetPdf()); assert (m_pdf);
    RooAbsCategoryLValue* m_cat = const_cast<RooAbsCategoryLValue*>(&m_pdf->indexCat());
    TList *m_dataList = m_data->split( *m_cat, true );
    RooDataSet* datai = ( RooDataSet* )( m_dataList->At( 0 ) );
    m_cat->setBin(0);
    TString channelname=m_cat->getLabel();
    RooAbsPdf* pdfi = m_pdf->getPdf(channelname);
    RooRealVar *x=dynamic_cast<RooRealVar*>(pdfi->getObservables(datai)->first());
    RooPlot *frame = x->frame();
    TH1D *_dataHist = new TH1D("ttt", "ttt", 8900, 1100, 10000);
    const int obsNBins = _dataHist->GetNbinsX();
    _dataHist = (TH1D*)datai->createHistogram("hpdf", *x);
    for( int ibin = 1 ; ibin <= obsNBins; ibin ++ ){
        _dataHist->SetBinError(ibin, sqrt(_dataHist->GetBinContent(ibin)));
    }
    datai->plotOn(frame, DataError(RooAbsData::Poisson), Binning(obsNBins));
    pdfi->plotOn(frame);
    gPad->SetLogy();
    frame->SetMinimum(1e-1);
    //frame->Draw();
    TH1D* hpdf = new TH1D("hpdf", "hpdf", obsNBins, x->getMin(), x->getMax());
    hpdf = (TH1D*)pdfi->createHistogram("hpdf", *x);
    hpdf->Scale(pdfi->expectedEvents(RooArgSet(*x))/hpdf->Integral());
    hpdf->SetMarkerColor(2);
    hpdf->SetLineColor(2);
    //hpdf->Draw();
    //_dataHist->Draw("same");
    TH1D *basicData = (TH1D*)_dataHist -> Clone();
    TH1D *basicBkg = (TH1D*)hpdf -> Clone();
    basicData->Rebin(45);
    basicBkg->Rebin(45);
    
    TFile *fout = TFile::Open("Roo.root", "recreate");
    hpdf->Write("hpdf");
    _dataHist->Write("_dataHist");
    basicData->Write("basicData");
    basicBkg->Write("basicBkg");
    fout->Close();
}
