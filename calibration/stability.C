#ifndef __headers__
#define __headers__
#include <vector>
#endif

#define NUMBER_OF_CHANNELS 8

void stability(string var, string type){
    //
    // var   - "energy", "rate","resolution"
    // type  - "abs", "rel"
    //
    TCanvas *c1 = new TCanvas("c1","c1",800,400);
    //
    //
    //
    gStyle->SetOptStat(0);
    cout << "stability:: plotting routine"<<endl;
    
    char hname[128];
    // draw master 1D histogram to hold the graphs
    hld->Draw();
    hld->GetXaxis()->SetTitle("time(sec)");
    
    //
    if(type == "rel") {
        hld->GetYaxis()->SetRangeUser(-0.05,0.05);
        hld->GetYaxis()->SetTitle("(v-v(t=0))/v(t=0)");
    } else {
        string yname;
        if ( var == "rate" ){
            yname = "Rate (Hz)";
        } else if ( var == "energy" ){
            yname = "Energy (keV)";
        } else if ( var == "resolution") {
            yname = "FWHM/E";
        }
        hld->GetYaxis()->SetTitle(yname.c_str());
    }
    
    int istyle;
    
    Double_t ymax = 0;
    // for now forget about the 2 empty detectors...
    int ioff = 0;
    if (yname == "energy") ioff = 1;
    for(int ich = 2; ich<NUMBER_OF_CHANNELS+ioff; ich++){
        if(ich<NUMBER_OF_CHANNELS){
            sprintf(hname,"%s_ch%d",var.c_str(),ich);
            istyle = 1;
        } else if (ich == NUMBER_OF_CHANNELS){
            sprintf(hname,"temperature");
            istyle = 24;
        }
        cout << hname <<endl;
        
        TGraph *gr_tmp = (TGraph*)gDirectory->Get(hname);
        TGraph *gr;
        
        // make an absolute value plot, or relative value
        if (type == "rel") {
            int n = gr_tmp->GetN();
            Double_t *x = gr_tmp->GetX();
            Double_t *y = gr_tmp->GetY();
            Double_t y0 = y[0];
            for (int i=0; i<n; i++) {
                y[i] = (y[i] - y0) / y0;
            }
            gr = new TGraph(n,x,y);
        } else {
            gr = gr_tmp;
            // find the maximum value of all graphs, so we can set the plot range
            Double_t yy = TMath::MaxElement(n,gr_tmp->GetY());
            if(yy>ymax) ymax = yy;
        }
        
        gr->SetLineColor(ich+1);
        gr->SetMarkerColor(ich+1);
        gr->SetMarkerStyle(istyle);
        gr->Draw("PL");
    }
    hld->GetYaxis()->SetRangeUser(0,ymax*1.2);
    
    c1->Update();
    return;
}
