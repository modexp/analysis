#define gain_cxx
#include "gain.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TVector.h>
#include <TParameter.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <stdio.h>
#include <TMath.h>
#include <TF1.h>

#define NUMBER_OF_CHANNELS 8
#define INDEX_TEMPERATURE 8


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


void gain::Loop()
{
    // energy calibration for modulation detectors
    if (fChain == 0) return;
    
    Long64_t nentries = fChain->GetEntriesFast();
    
    // private code below
    TFile *_f = new TFile(gain_file.c_str(),"RECREATE");
    
    vector<TH1F*> _e_all,_e_good,_e_err1,_e_err2,_pk_tmp;
    TH1F *_T = new TH1F("T","T",1000,10+273.15,40+273.15);
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
    
    vector<float> interval_time;
    vector< vector<float> *> var_select;
    for(int ich=0; ich<NUMBER_OF_CHANNELS+1; ich++) var_select.push_back(new vector<float>);
    //vector<float> peak_value;
    
    cout<<"Start event loop.... nentries ="<<nentries<<endl;
    Long64_t nbytes = 0, nb = 0;
    Double_t t0,tstart,ctime;
    Double_t dt_max = 900.;
    
    
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        if(jentry==0) {
            tstart = time;
            cout << "tstart = "<<tstart<<endl;
        }
        ctime = time - tstart;
        
        if (jentry == 0) {
            t0 = ctime;
        }
        
        // analyze the peak position for your selected detector channel
        channel = channel % 100 ;
        _T->Fill(temp+273.15);
        
        _e_all[channel]->Fill(integral);
        if      (error == 0) {
            _pk_tmp[channel]->Fill(integral);
            _e_good[channel]->Fill(integral);
        }
        else if ((error&0x01)!=0) _e_err1[channel]->Fill(integral);
        else if ((error&0x02)!=0) _e_err2[channel]->Fill(integral);
        
        if(jentry%500000 == 0) {
            cout<<"Processed "<<jentry<<" events"<<endl;
        }
        
        // if we exceed the maximum time interval, we store the peak information
        if(ctime - t0 > dt_max){
            cout<<"SPECIAL:: Processed "<<jentry<<" events  ctime ="<<ctime<<endl;
            
            for(int ich=0; ich<NUMBER_OF_CHANNELS; ich++){
                int maxbin = _pk_tmp[ich]->GetMaximumBin();
                Double_t energy = _pk_tmp[ich]->GetBinCenter(maxbin);
                TF1 *func = new TF1("fit",fitf,energy-100,energy+100,3);
                func->SetParameters(10000,energy,25);
                func->SetParNames("C","mean","sigma");
                _pk_tmp[ich]->Fit("fit");
                energy = func->GetParameter(1);
                var_select.at(ich)->push_back(energy);
                cout << ich <<" "<<energy<<" keV" <<endl;
                _pk_tmp[ich]->Reset();
                delete func;
            }
            // temperature
            double tmean = _T->GetMean();
            cout <<"T = "<<tmean<<endl;
            var_select.at(INDEX_TEMPERATURE)->push_back(tmean);
            interval_time.push_back((t0+ctime)/2.); // almost right....
            
            // reset the time for the start of the next interval
            t0 = ctime;
            
            _T->Reset();
        }
    }
    cout<<"Done event loop...."<<endl;
    for(int ich = 0; ich<NUMBER_OF_CHANNELS; ich++){
        _e_all[ich]->Write();
        _e_good[ich]->Write();
        _e_err1[ich]->Write();
        _e_err2[ich]->Write();
    }
    
    TCanvas *c1 = new TCanvas("c1","c1",800,400);
    c1->Draw();
    // historam as placeholder for drawing
    TH1F *_holder = new TH1F("hld","hld",1,0,ctime);
    sprintf(hname,"Gain vs time run = %s",run.c_str());
    _holder->SetTitle(hname);
    //_holder->Fill(0.01,250);
    _holder->GetXaxis()->SetTitle("time (sec)");
    _holder->GetYaxis()->SetTitle("(value - value(t=0)) / value(t=0)");
    _holder->GetYaxis()->SetRangeUser(-0.05,0.05);
    _holder->Draw();
    c1->Update();
    
    // loop over all the vectors and make teh stability graphs
    vector<TGraph*> _gr;
    char gname[256];
    for(int ich=0; ich<NUMBER_OF_CHANNELS+1; ich++){
        Int_t n = (Int_t)interval_time.size();
        double vmean = accumulate(var_select.at(ich)->begin(),var_select.at(ich)->end(),0.0)/var_select.at(ich)->size();
        double v0 = var_select.at(ich)->at(0);
        Float_t x[n];
        Float_t y[n];
        for (int i=0 ; i<n ; i++)
        {
            x[i]=interval_time.at(i);
            y[i]=(var_select.at(ich)->at(i) - v0)/v0;
        }
        _gr.push_back(new TGraph(n,x,y));
        if(ich<NUMBER_OF_CHANNELS){
            sprintf(gname,"peak_ch%d",ich);
        } else {
            sprintf(gname,"temperature");
            _gr[ich]->SetMarkerStyle(24);
            _gr[ich]->SetLineColor(2);;
            _gr[ich]->SetMarkerColor(2);;
        }
        _gr[ich]->SetName(gname);
        _gr[ich]->Write();
        _gr[ich]->SetLineColor(ich+1);
        if(ich>1 && ich<NUMBER_OF_CHANNELS) {
            _gr[ich]->Draw("L");
        } else if (ich == NUMBER_OF_CHANNELS){
            _gr[ich]->Draw("LP");
        }
        
    }
    
    c1->Update();
    string pdfname = gain_file+".pdf";
    c1->Print(pdfname.c_str());
    //_f->Write();
    //_f->Close();
}
