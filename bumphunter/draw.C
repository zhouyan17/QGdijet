#include "../ATLASStyle/AtlasStyle.C"
#include "../ATLASStyle/AtlasUtils.C"
#include "../ATLASStyle/AtlasLabels.C"

double sigrange = 4;
double bhv = 10.4211;
double chiv = 100.937;
double llv = 691.916;
double rangel = 2855;
double ranger = 4610;
double yl = 1e4;
double yr = 2e2;
TString bhp = "0.019";

void draw() {
    SetAtlasStyle ();
    TCanvas *c = new TCanvas("c","c",1100,750);
    TFile *f2 = new TFile("results/b.root");
    TH1D* diff = (TH1D*) f2->Get("residualHist");
    TH1D* sig = new TH1D("sig", "sig", 20*sigrange, -sigrange, sigrange);
    for (int bin=0; bin<=diff->GetNbinsX(); bin++) {
        sig->Fill(diff->GetBinContent(bin));
    }
    gPad->SetLogy(0);
    sig->SetFillColor(2);
    sig->SetXTitle("Residuals");
    sig->SetYTitle("Entries");
    sig->Draw("][hist");
    c->SaveAs("../p2.png");
    
    c->Clear();
    gPad->SetLogy(1);
    gPad->SetLogx(1);
    TGraphErrors* error = (TGraphErrors*) f2->Get("bumpHunterTomographyFromPseudoexperiments");
    error->SetLineColor(2);
    error->SetMarkerStyle(0);
    error->SetTitle(";Dijet Mass [GeV];Poisson PVal of Interval");
    error->Draw("AP");
    ATLASLabel(0.6,0.2,"Preliminary");
    c->SaveAs("../p3.png");
    
    TH1D* bump = (TH1D*) f2->Get("bumpHunterStatHistNullCase");
    gPad->SetLogx(0);
    bump->SetFillColor(7);
    bump->SetAxisRange(0.5, 3000, "Y");
    bump->SetXTitle("BumpHunter");
    bump->SetYTitle("Pseudo-exeperiments");
    bump->Draw("][hist");
    TArrow *a1 =new TArrow(bhv,1,bhv,0.5,0.03,">");
    a1->SetLineColor(2);
    a1->SetLineWidth(4);
    a1->Draw();
    myText(   0.2,  0.8, 1, "#sqrt{s}= 13 TeV");
    myBoxText( 0.26, 0.73, 0.05, 7, "Pseudo-exeperiments");
    myMarkerText( 0.26, 0.68, 2, 20, "Value in Data",1);
    ATLASLabel(0.2,0.85,"Preliminary");
    c->SaveAs("../p4.png");
    
    TH1D* chi2 = (TH1D*) f2->Get("chi2StatHistNullCase");
    chi2->SetFillColor(7);
    chi2->SetAxisRange(0.5, 3000, "Y");
    chi2->SetXTitle("#chi^{2}");
    chi2->SetYTitle("Pseudo-exeperiments");
    chi2->Draw("][hist");
    TArrow *a2 =new TArrow(chiv,1,chiv,0.5,0.03,">");
    a2->SetLineColor(2);
    a2->SetLineWidth(4);
    a2->Draw();
    myText(   0.2,  0.8, 1, "#sqrt{s}= 13 TeV");
    myBoxText( 0.26, 0.73, 0.05, 7, "Pseudo-exeperiments");
    myMarkerText( 0.26, 0.68, 2, 20, "Value in Data",1);
    ATLASLabel(0.2,0.85,"Preliminary");
    c->SaveAs("../p5.png");
    
    TH1D* llhd = (TH1D*) f2->Get("logLikelihoodStatHistNullCase");
    llhd->SetFillColor(7);
    llhd->SetAxisRange(0.5, 3000, "Y");
    llhd->SetXTitle("logLikelihood");
    llhd->SetYTitle("Pseudo-exeperiments");
    llhd->Draw("][hist");
    TArrow *a3 =new TArrow(llv,1,llv,0.5,0.03,">");
    a3->SetLineColor(2);
    a3->SetLineWidth(4);
    a3->Draw();
    myText(   0.2,  0.8, 1, "#sqrt{s}= 13 TeV");
    myBoxText( 0.26, 0.73, 0.05, 7, "Pseudo-exeperiments");
    myMarkerText( 0.26, 0.68, 2, 20, "Value in Data",1);
    ATLASLabel(0.2,0.85,"Preliminary");
    c->SaveAs("../p6.png");
    
    TPad *pad1 = new TPad("pad1","pad1",0,0.31,1,1);
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.31);
    pad1 -> SetBottomMargin(0);
    pad1 -> SetLeftMargin(0.12);
    pad1 -> SetLogy(1);
    pad1 -> Draw();
    pad2 -> SetBottomMargin(.4);
    pad2 -> SetTopMargin(0);
    pad2 -> SetLeftMargin(0.12);
    pad2 -> SetGridx(1);
    pad2 -> SetGridy(1);
    pad2 -> Draw();
    TH1D* data = (TH1D*) f2->Get("basicData");
    TH1D* bkg = (TH1D*) f2->Get("basicBkg");
    data->SetAxisRange(1000, 10000, "X");
    data->SetAxisRange(0.03, 3000000000, "Y");
    data->SetMarkerStyle(20);
    data->SetMarkerSize(0.8);
    data->SetYTitle("Events");
    bkg->SetLineColor(2);
    pad1->cd();
    data->Draw("E");
    bkg->Draw("][histsame");
    TArrow *edge1 =new TArrow(rangel,yl,rangel,0.01,0.03,"-");
    TArrow *edge2 =new TArrow(ranger,yr,ranger,0.01,0.03,"-");
    edge1->SetLineColor(4);
    edge2->SetLineColor(4);
    edge1->Draw();
    edge2->Draw();
    myText(   0.2,  0.78, 1, "#sqrt{s}= 13 TeV");
    myBoxText( 0.6, 0.73, 0.05, 2, "Background fit");
    myBoxText( 0.6, 0.66, 0.05, 4, "BumpHunter interval");
    myMarkerText( 0.6, 0.8, 1, 20, "Data",1);
    myText(   0.61,  0.4, 1, "p-value = " + bhp);
    ATLASLabel(0.2,0.85,"Preliminary");
    pad2->cd();
    diff->SetAxisRange(1000, 10000, "X");
    diff->SetAxisRange(-sigrange, sigrange, "Y");
    diff->SetFillColor(2);
    diff->SetXTitle("m_{jj}[GeV]");
    diff->SetYTitle("Significance");
    diff->Draw("][hist");
    TArrow *edge11 =new TArrow(rangel,10,rangel,-sigrange,0.03,"-");
    TArrow *edge22 =new TArrow(ranger,10,ranger,-sigrange,0.03,"-");
    edge11->SetLineColor(4);
    edge22->SetLineColor(4);
    edge11->Draw();
    edge22->Draw();
    
    c->SaveAs("../p1.png");
}


/*
 ----------------------dijetTotal---------------------------
 double sigrange = 4;
 double bhv = 10.414;
 double chiv = 100.937;
 double llv = 691.916;
 double rangel = 2855;
 double ranger = 4610;
 double yl = 1e4;
 double yr = 2e2;
 TString bhp = "0.019";
 ----------------------dijetTagged---------------------------
 double sigrange = 4;
 double bhv = 7.13025;
 double chiv = 57.5287;
 double llv = 403.242;
 double rangel = 3080;
 double ranger = 4925;
 double yl = 4e2;
 double yr = 1;
 TString bhp = "0.089";
*/
