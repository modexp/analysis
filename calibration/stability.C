#ifndef _includes_
#define _includes_
#include <vector>
#include "TDirectory.h"
#include "TH1F.h"
#include "TGraph.h"
#include <iostream>
#endif

#define NUMBER_OF_CHANNELS 8

void stability(string var)
{
//    string var = "peak";
    //
    // var - "peak","rate","resolution"
    //
    TCanvas *c1 = new TCanvas();

    gStyle->SetOptStat(0);
    cout << "plotting routine"<<endl;
    
    vector<TGraph*> _gr;
    
    char hname[128];
    TH1F *hld = (TH1F*)gDirectory->Get("hld");
    hld->Draw();
    hld->GetXaxis()->SetTitle("time(sec)");
    hld->GetYaxis()->SetRangeUser(-0.02,0.02);
    hld->GetYaxis()->SetTitle("(v - v(t=0)) / v(t=0)");
    for(int ich = 0; ich<NUMBER_OF_CHANNELS; ich++){
        sprintf(hname,"%s_ch%d",var.c_str(),ich);
        cout << hname <<endl;
        _gr.push_back((TGraph*)gDirectory->Get(hname));
        _gr[ich]->SetLineColor(ich+1);
        _gr[ich]->SetMarkerColor(ich+1);
        if(ich>1) _gr[ich]->Draw("PL");
    }
    c1->Update();
}
