#define analyzer_cxx
#include "analyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TVector.h>
#include <TParameter.h>
#include <iostream>
//#include <vector>
//#include <numeric>
#include <stdio.h>
#include <TMath.h>
#include <TF1.h>

#define NUMBER_OF_CHANNELS 8
#define INDEX_TEMPERATURE 24


// ranges for plotting
const int   nbin = 600;
const float emin = 0.; // in keV
const float emax = 3000.; // in keV
const float adc_max_volt = 2.;
const float base_max_val = 2000;

Double_t fitf(Double_t *v, Double_t *par)
{
    Double_t arg = 0;
    if (par[2] != 0) arg = (v[0] - par[1])/par[2];
    
    Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
    fitval += par[3] + par[4]*v[0] + par[5]*v[0]*v[0];

    return fitval;
}

void analyzer::book_histograms(){
    _f = new TFile(analyzer_file.c_str(),"RECREATE");
    
    char hname[128];
    // book histograms
    for (int ich = 0; ich<NUMBER_OF_CHANNELS; ich++){
        // spectra
        sprintf(hname,"_e_all_ch%1d",ich);
        _e_all.push_back(new TH1F(hname,hname,nbin,emin,emax));
        sprintf(hname,"_e_good_ch%1d",ich);
        _e_good.push_back(new TH1F(hname,hname,nbin,emin,emax));
        sprintf(hname,"_e_err1_ch%1d",ich);
        _e_err1.push_back(new TH1F(hname,hname,nbin,emin,emax));
        sprintf(hname,"_e_err2_ch%1d",ich);
        _e_err2.push_back(new TH1F(hname,hname,nbin,emin,emax));
        
        // temporary histograms for stability measurements
        sprintf(hname,"_pk_tmp%1d",ich);
        _pk_tmp.push_back(new TH1F(hname,hname,nbin,emin,emax));
    }
    
    // temporary histogram for stability measurements
    _T = new TH1F("T","T",1000,10+273.15,40+273.15);
    
}

void analyzer::fill_histograms(){
    // fill all the histograms
    
    channel = channel % 100 ;
    _T->Fill(temp+273.15);
    
    _e_all[channel]->Fill(integral);
    if      (error == 0) {
        _pk_tmp[channel]->Fill(integral);
        _e_good[channel]->Fill(integral);
    }
    else if ((error&0x01)!=0) _e_err1[channel]->Fill(integral);
    else if ((error&0x02)!=0) _e_err2[channel]->Fill(integral);
}

void analyzer::get_interval_data(){
    //
    // analyze the data from a time interval.
    // length of interval set by TIME_INTERVAL which can be set in header file
    //
    cout<<"analyzer::get_interval_data:: time_since_start ="<<time_since_start<<endl;
    
    int huh;
//    TCanvas *c1 = new TCanvas("c1","c1",600,400);
    Double_t bin_width = (emax-emin)/nbin;
    for(int ich=0; ich<NUMBER_OF_CHANNELS; ich++){
        int maxbin = _pk_tmp[ich]->GetMaximumBin();
        Double_t energy = _pk_tmp[ich]->GetBinCenter(maxbin);

        _pk_tmp[ich]->GetXaxis()->SetRangeUser(energy-200,energy+200);
//        _pk_tmp[ich]->Draw();

        TF1 *func = new TF1("fit",fitf,energy-200,energy+200,5);
        func->SetParameters(10000,energy,25);
        func->SetParNames("C","mean","sigma");
        _pk_tmp[ich]->Fit("fit","","",energy-100.,energy+75.);
        
//        c1->Update();
//        cin>>huh;
        
        Double_t peak        = func->GetParameter(0);
        energy               = func->GetParameter(1);
        Double_t sigma       = func->GetParameter(2);
        Double_t resolution  = 0;
        if(energy>0) resolution = 2.355*sigma/energy ;
        
        Double_t rate = TMath::Sqrt(2*TMath::Pi())*sigma*peak / TIME_INTERVAL / bin_width;
        cout <<"get_interval_data:: ich ="<<ich<<" E = "<<energy<<" keV  rate = "<<rate<<" Hz  resolution  = "<<resolution<<" % "<<endl;
        
        var_select.at(ich)->push_back(rate);
        var_select.at(ich+  NUMBER_OF_CHANNELS)->push_back(energy);
        var_select.at(ich+2*NUMBER_OF_CHANNELS)->push_back(resolution);
        
        _pk_tmp[ich]->Reset();
        delete func;
    }
    // temperature
    double tmean = _T->GetMean();
    
    var_select.at(INDEX_TEMPERATURE)->push_back(tmean);
    interval_time.push_back((t0+time_since_start)/2.); // almost right....
    
    _T->Reset();
}

void analyzer::write_histograms(){
    //
    // write a few parameters to file
    //
    TParameter<Double_t> * tstartPar = new TParameter<Double_t>("t0",tstart);
    tstartPar->Write();
    TParameter<Double_t> * tendPar = new TParameter<Double_t>("runtime",time_since_start);
    tendPar->Write();
    //
    // write histograms to the output root file
    //
    for(int ich = 0; ich<NUMBER_OF_CHANNELS; ich++){
        _e_all[ich]->Write();
        _e_good[ich]->Write();
        _e_err1[ich]->Write();
        _e_err2[ich]->Write();
    }
    
    char hname[128];
    
    // historam as placeholder for drawing stability graphs
    TH1F *_holder = new TH1F("hld","hld",1,0,time_since_start);
    string htitle = "run: "+run;
    _holder->SetTitle(htitle.c_str());
    //
    // loop over all the vectors and make the stability graphs
    //
    vector<TGraph*> _gr;
    char gname[256];
    for(int ich=0; ich<3*NUMBER_OF_CHANNELS+1; ich++){
        string type="";
        if(ich<3*NUMBER_OF_CHANNELS){
            if( ich/NUMBER_OF_CHANNELS == 0) {
                type = "rate";
                sprintf(gname,"%s_ch%d",type.c_str(),ich%NUMBER_OF_CHANNELS);
            } else if (ich/NUMBER_OF_CHANNELS == 1) {
                type = "energy";
                sprintf(gname,"%s_ch%d",type.c_str(),ich%NUMBER_OF_CHANNELS);
            } else if (ich/NUMBER_OF_CHANNELS == 2) {
                type = "resolution";
                sprintf(gname,"%s_ch%d",type.c_str(),ich%NUMBER_OF_CHANNELS);
            }
        } else if (ich == INDEX_TEMPERATURE){
            type = "temperature";
            sprintf(gname,"%s",type.c_str());
        }
        cout << " graph for ich = "<<ich<< " variable = "<<type<<endl;
        Int_t n = (Int_t)interval_time.size();
        //double vmean = accumulate(var_select.at(ich)->begin(),var_select.at(ich)->end(),0.0)/var_select.at(ich)->size();
        double v0 = var_select.at(ich)->at(0);
        Float_t x[n];
        Float_t y[n];
        for (int i=0 ; i<n ; i++)
        {
            x[i]=interval_time.at(i);
            y[i]= var_select.at(ich)->at(i);

        }
        _gr.push_back(new TGraph(n,x,y));
        _gr[ich]->SetName(gname);
        _gr[ich]->Write();
    }
    _holder->Write();
}

void analyzer::Loop()
{
    // energy calibration for modulation detectors
    if (fChain == 0) return;
    
    Long64_t nentries = fChain->GetEntriesFast();
    //
    // book histograms
    //
    book_histograms();
    
    // initialize the data vectors for teh time interval analysis
    for(int ich=0; ich<3*NUMBER_OF_CHANNELS+1; ich++) var_select.push_back(new vector<float>);
    
    //
    // start the event loop
    //
    cout<<"Start event loop.... nentries ="<<nentries<<endl;
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        //
        // get entry from the tree
        //
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        //
        // process the time information.
        //
        if (jentry == 0) tstart = time;
        time_since_start = time - tstart;
        if (jentry == 0) t0 = time_since_start;
        
        //
        // fill the monitoring histograms
        //
        fill_histograms();
        //
        // if we exceed the maximum time interval, get all the data recorded 
        // during this time. then reset time for a new interval....
        //
        if(time_since_start - t0 > TIME_INTERVAL) {
           get_interval_data();
           // reset the time for the start of the next interval
           t0 = time_since_start;
        }
        
        if(jentry%500000 == 0) cout<<"Processed "<<jentry<<" events"<<endl;
    }
    
    cout<<"Done event loop...."<<endl;
    
    //
    // write histograms to file and generate the stability graphs
    //
    write_histograms();
    cout<<"Done writing histograms"<<endl;
}
