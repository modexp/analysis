#ifndef __headers__
#define __headers__
#include <vector>
#endif

#define NUMBER_OF_CHANNELS 8

string sources[] = {"none","none","^{44}Ti (511keV)","^{44}Ti (511keV)","^{60}Co (1173keV)","^{60}Co (1173keV)","^{137}Cs (662keV)","^{137}Cs (662keV)"};

TCanvas *c1 = new TCanvas("c1","c1",800,400);

string get_figname(string fname){
  string figname = "";
  figname = fname;
  size_t pos = figname.find_last_of("/")+1;
  figname = figname.substr(pos);
  pos = figname.find_last_of(".");
  figname = figname.substr(0,pos);
  pos = figname.find_first_of("_")+1;
  figname = figname.substr(pos);

  return figname;
}


void stability(string rootfile, string var, string type, bool save_plot){
    //
    // Plot rate, energy, resolution as a function of time.
    //
    // Input: rootfile  - input root filename
    //        var       - "energy", "rate","resolution"
    //        type      - "abs", "rel"
    //        save_plot - save plot to .pdf file (or other format)
    //
    // A.P. Colijn
    //
    TFile *_f = new TFile(rootfile.c_str(),"READONLY");
   
    gStyle->SetOptStat(0);
    cout << "stability:: plotting routine"<<endl;
    
    char hname[128],cmd[128];
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
    if (var == "energy") ioff = 1;
    // make a legend
    TLegend *leg = new TLegend(0.65,0.63,0.95,0.89);
    leg->SetFillStyle(0);

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

        // add legend entry
        if( ich<NUMBER_OF_CHANNELS){
          sprintf(cmd, "ch%d - %s ", ich, sources[ich].c_str());
          leg->AddEntry(gr,cmd,"l");
        }
    }
    if(type == "abs") hld->GetYaxis()->SetRangeUser(0,ymax*1.6);
    
    // draw the legend
    leg->SetBorderSize(0);
    leg->Draw();

    c1->Update();

    if(save_plot){
       string figname = "plots/"+get_figname(rootfile)+"_"+var+"_"+type+".pdf";
       c1->Print(figname.c_str());
       string figname = "plots/"+get_figname(rootfile)+"_"+var+"_"+type+".png";
       c1->Print(figname.c_str());
    } 
    return;
}
