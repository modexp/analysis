#define gain_cxx
#include "gain.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TVector.h>
#include <TParameter.h>
#include <iostream>
//#include <vector>
#include <numeric>
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
    return fitval;
}

void gain::book_histograms(){
    _f = new TFile(gain_file.c_str(),"RECREATE");
    
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

void gain::fill_histograms(){
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

void gain::get_interval_data(){
    //
    // analyze the data from a time interval.
    // length of interval set by TIME_INTERVAL which can be set in header file
    //
    cout<<"gain::get_interval_data:: time_since_start ="<<time_since_start<<endl;
    
    for(int ich=0; ich<NUMBER_OF_CHANNELS; ich++){
        int maxbin = _pk_tmp[ich]->GetMaximumBin();
        Double_t energy = _pk_tmp[ich]->GetBinCenter(maxbin);
        TF1 *func = new TF1("fit",fitf,energy-50,energy+50,3);
        func->SetParameters(10000,energy,25);
        func->SetParNames("C","mean","sigma");
        _pk_tmp[ich]->Fit("fit");
        
        Double_t peak   = func->GetParameter(0);
        energy          = func->GetParameter(1);
        Double_t sigma  = func->GetParameter(2);
        
        var_select.at(ich)->push_back(peak);
        var_select.at(ich+  NUMBER_OF_CHANNELS)->push_back(energy);
        var_select.at(ich+2*NUMBER_OF_CHANNELS)->push_back(sigma);
        
        _pk_tmp[ich]->Reset();
        delete func;
    }
    // temperature
    double tmean = _T->GetMean();
    
    var_select.at(INDEX_TEMPERATURE)->push_back(tmean);
    interval_time.push_back((t0+time_since_start)/2.); // almost right....
    
    // reset the time for the start of the next interval
    t0 = time_since_start;
    
    _T->Reset();
}

void gain::write_histograms(){
    //
    // write historgams to the output root file
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
        cout << " graph for ich = "<<ich<<endl;
        Int_t n = (Int_t)interval_time.size();
        //double vmean = accumulate(var_select.at(ich)->begin(),var_select.at(ich)->end(),0.0)/var_select.at(ich)->size();
        double v0 = var_select.at(ich)->at(0);
        Float_t x[n];
        Float_t y[n];
        for (int i=0 ; i<n ; i++)
        {
            x[i]=interval_time.at(i);
            y[i]=(var_select.at(ich)->at(i) - v0)/v0;
        }
        _gr.push_back(new TGraph(n,x,y));
        if(ich<3*NUMBER_OF_CHANNELS){
            if( ich/NUMBER_OF_CHANNELS == 0) {
                sprintf(gname,"rate_ch%d",ich%NUMBER_OF_CHANNELS);
            } else if (ich/NUMBER_OF_CHANNELS == 1) {
                sprintf(gname,"peak_ch%d",ich%NUMBER_OF_CHANNELS);
            } else if (ich/NUMBER_OF_CHANNELS == 2) {
                sprintf(gname,"resolution_ch%d",ich%NUMBER_OF_CHANNELS);
            }
        } else {
            sprintf(gname,"temperature");
        }
        _gr[ich]->SetName(gname);
        _gr[ich]->Write();
    }
    _holder->Write();
}

void gain::Loop()
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
        // if we exceed the maximum time interval, get all the data recorded during this time. then reset time for a new interval....
        //
        if(time_since_start - t0 > TIME_INTERVAL) get_interval_data();
        
        if(jentry%500000 == 0) cout<<"Processed "<<jentry<<" events"<<endl;
    }
    
    cout<<"Done event loop...."<<endl;
    
    //
    // write histograms to file and generate the stability graphs
    //
    write_histograms();
    cout<<"Done writing histograms"<<endl;
}
