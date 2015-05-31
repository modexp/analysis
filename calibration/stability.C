#ifndef __headers__
#define __headers__
#include <vector>
#endif

#define NUMBER_OF_CHANNELS 8
#define MAX_PEAKS 5

float v0[NUMBER_OF_CHANNELS][MAX_PEAKS];

string sources[] = {"none","none","^{44}Ti (511keV)","^{44}Ti (511keV)","^{60}Co (1173keV)","^{60}Co (1173keV)","^{137}Cs (662keV)","^{137}Cs (662keV)"};

TCanvas *c1 = new TCanvas("c1","c1",800,400);

string get_figname(string fname){
  // extract figure name from the filename
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
    char hname[128],cmd[128],cut[128];
    TFile *_f = new TFile(rootfile.c_str(),"READONLY");

    // add files to the chain .....
//    sprintf(cmd,"%s*.root",rootfile.c_str());
//    run->Add(cmd);
   
    gStyle->SetOptStat(0);
    cout << "stability:: plotting routine"<<endl;
    
    string runname = run->GetTitle();
    cout << "runname = "<<runname<<endl;
    Double_t tmax = runtime->GetVal();   

    TH1F *_hld = new TH1F("hld","hld",1,0.,tmax);
    sprintf(hname,"run: %s",runname.c_str());
    _hld->SetTitle(hname);
    _hld->GetXaxis()->SetTitle("time (sec)");
    _hld->GetYaxis()->SetRangeUser(0.,300);
    _hld->Draw();
    if(type == "rel") {
       _hld->GetYaxis()->SetRangeUser(-0.05,0.05);
       _hld->GetYaxis()->SetTitle("(v - v(t=0)) / v(t=0)");
    } else {
        string yname;
        if ( var == "rate" ){
            yname = "Rate (Hz)";
        } else if ( var == "e" ){
            yname = "Energy (keV)";
        } else if ( var == "res") {
            yname = "FWHM/E";
        }
       _hld->GetYaxis()->SetTitle(yname.c_str());
    }
    
    int istyle;
    
    Double_t ymax = 0;
    // make a legend
    TLegend *leg = new TLegend(0.65,0.63,0.95,0.89);
    leg->SetFillStyle(0);

    // get the values of the parameter you want to plot at t=0.
    // needed when plotting relative values
    Double_t v_tmp;
    Int_t    ich_tmp,ipk_tmp;
    sprintf(cmd,"%s",var.c_str());
    
    ana->SetBranchAddress(cmd,&v_tmp); 
    ana->SetBranchAddress("channel",&ich_tmp); 
    ana->SetBranchAddress("peak",&ipk_tmp); 
    for(int i=MAX_PEAKS*NUMBER_OF_CHANNELS; i>=0; i--){
       ana->GetEntry(i); 
       v0[ich_tmp][ipk_tmp] = v_tmp;
    }

    // loop over all the channels
    for(int ich = 2; ich<NUMBER_OF_CHANNELS; ich++){
        if(type == "abs"){
          sprintf(cmd,"%s:time",var.c_str(),ich);
        } else {
          sprintf(cmd,"(%s-%f)/%f:time",var.c_str(),v0[ich][0],v0[ich][0]);
        }
        //sprintf(cut,"(channel==%d)",ich);
     
        sprintf(cut,"(channel==%d)&&(peak==0)",ich);
        ana->SetMarkerStyle(1);
        ana->SetMarkerColor(ich+1);
        ana->SetLineColor(ich+1);
        ana->Draw(cmd,cut,"lsame");
        
        // add legend entry
        sprintf(cmd, "ch%d - %s ", ich, sources[ich].c_str());
        TLine *l1 = new TLine();
        l1->SetLineColor(ich+1);
        leg->AddEntry(l1,cmd,"l");
    }

    if(type == "abs"){
       sprintf(cmd,"%s",var.c_str());
       Double_t ymax = ana->GetMaximum(cmd);
       _hld->GetYaxis()->SetRangeUser(0.,1.6*ymax);
    }
    
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
