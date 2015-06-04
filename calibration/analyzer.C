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
#define MAX_INDEX 400


#define MAX_PEAKS 5
float source_energy[NUMBER_OF_CHANNELS][MAX_PEAKS] =
//
// the energy peaks you wish to select should be in this list
// NOTE: the first peak should be the highest in the spectrum (sub-optimal, but handy for finding)
//
{
    {1460.,-1,-1,-1,-1}, // channel0: no source
    {1460.,-1,-1,-1,-1}, // channel1: no source
    {511.,1157.020,511.+1157.020,-1,-1}, // channel2: 44Ti
    {511.,1157.020,511.+1157.020,-1,-1}, // channel3: 44Ti
    {1173.2,1332.5,1173.2+1332.5,-1,-1}, // channel4: 60Co
    {1173.2,1332.5,1173.2+1332.5,-1,-1}, // channel5: 60Co
    {662.,-1,-1,-1,-1}, // channel6: 137Cs
    {662.,-1,-1,-1,-1}  // channel7: 137Cs
};

//
// ranges for plotting
//
const int   nbin0 = 600;
// number of bins for the temporary fit histograms.... channels 0+1 have few entries so wider bins
float nbin[NUMBER_OF_CHANNELS]={nbin0/4,nbin0/4,nbin0,nbin0,nbin0,nbin0,nbin0,nbin0};
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
        _e_all.push_back(new TH1F(hname,hname,nbin0,emin,emax));
        sprintf(hname,"_e_good_ch%1d",ich);
        _e_good.push_back(new TH1F(hname,hname,nbin0,emin,emax));
        sprintf(hname,"_e_err1_ch%1d",ich);
        _e_err1.push_back(new TH1F(hname,hname,nbin0,emin,emax));
        sprintf(hname,"_e_err2_ch%1d",ich);
        _e_err2.push_back(new TH1F(hname,hname,nbin0,emin,emax));
        sprintf(hname,"_e_err4_ch%1d",ich);
        _e_err4.push_back(new TH1F(hname,hname,nbin0,emin,emax));
        // baseline
        sprintf(hname,"_b_good_ch%1d",ich);
        _b_good.push_back(new TH1F(hname,hname,nbin0,0,base_max_val));
        sprintf(hname,"_b_err1_ch%1d",ich);
        _b_err1.push_back(new TH1F(hname,hname,nbin0,0,base_max_val));
        sprintf(hname,"_b_err2_ch%1d",ich);
        _b_err2.push_back(new TH1F(hname,hname,nbin0,0,base_max_val));
        sprintf(hname,"_b_err4_ch%1d",ich);
        _b_err4.push_back(new TH1F(hname,hname,nbin0,0,base_max_val));
        // integral vs peak
        sprintf(hname,"_h_vs_E_good_ch%1d",ich);
        _2d_good.push_back(new TH2F(hname,hname,nbin0,emin,emax,nbin0,0.,adc_max_volt));
        sprintf(hname,"_h_vs_E_err1_ch%1d",ich);
        _2d_err1.push_back(new TH2F(hname,hname,nbin0,emin,emax,nbin0,0.,adc_max_volt));
        sprintf(hname,"_h_vs_E_err2_ch%1d",ich);
        _2d_err2.push_back(new TH2F(hname,hname,nbin0,emin,emax,nbin0,0.,adc_max_volt));
        sprintf(hname,"_h_vs_E_err4_ch%1d",ich);
        _2d_err4.push_back(new TH2F(hname,hname,nbin0,emin,emax,nbin0,0.,adc_max_volt));
        
        // temporary histograms for stability measurements
        sprintf(hname,"_pk_tmp%1d",ich);
        _pk_tmp.push_back(new TH1F(hname,hname,nbin[ich],emin,emax));
    }
    
    // temporary histogram for temperature measurements
    _T = new TH1F("T","T",1000,10+273.15,40+273.15);
    
    // Define tree and branches
    tree = new TTree("ana", "Analyzed spectra");
    
    tree->Branch("t0", &_t_t0, "t0/D");
    tree->Branch("time", &_t_time, "time/D");
    tree->Branch("channel", &_t_chanNum, "channel/I");
    tree->Branch("peak", &_t_peakNum, "peak/I");
    tree->Branch("rate", &_t_rate, "rate/D");
    tree->Branch("e", &_t_energy, "e/D");
    tree->Branch("res", &_t_res, "res/D");
    tree->Branch("temp", &_t_temp, "temp/D");
}

void analyzer::fill_histograms(){
    // fill all the histograms
    
    channel = channel % 100 ;
    _T->Fill(temp+273.15);
    
    _e_all[channel]->Fill(integral);
    if      (error == 0) {
        _pk_tmp[channel]->Fill(integral);
        _e_good[channel]->Fill(integral);
        _b_good[channel]->Fill(baseline);
        
        _2d_good[channel]->Fill(integral,height);
    }
    else if ((error&0x01)!=0) {
        _e_err1[channel]->Fill(integral);
        _b_err1[channel]->Fill(baseline);
        
        _2d_err1[channel]->Fill(integral,height);
    }
    else if ((error&0x02)!=0) {
        _e_err2[channel]->Fill(integral);
        _b_err2[channel]->Fill(baseline);
        
        _2d_err2[channel]->Fill(integral,height);
    }
    else if ((error&0x04)!=0) {
    // temporary error handling: err4 is not a real error yet
        _e_err4[channel]->Fill(integral);
        _b_err4[channel]->Fill(baseline);
        
        _2d_err4[channel]->Fill(integral,height);
    }
}

void analyzer::get_interval_data(){
    //
    // analyze the data from a time interval.
    // length of interval set by TIME_INTERVAL which can be set in header file
    //
    
    // time
    _t_t0    = tstart;
    _t_time  = (t0+time_since_start)/2.;
    
    // temperature
    _t_temp  = _T->GetMean();
    _T->Reset();
    
    cout<<"analyzer::get_interval_data:: time_since_start ="<<time_since_start<<endl;
    
    //    int huh;
    TCanvas *c1 = new TCanvas("c1","c1",600,400);
    int huh;
    for(int ich=0; ich<NUMBER_OF_CHANNELS; ich++){
        //
        // find all the selected energy peaks
        //
        int      maxbin;
        double   maxval;
        Double_t e_start, e0, demin, demax;
        for (int ipeak=0; ipeak<MAX_PEAKS; ipeak++){
            if(source_energy[ich][ipeak] >0){
                //
                // find the fit starting values
                //
                
                //
                // first peak is special.... we use the GetMaximumBin() method in order
                // to find this peak even if there is a shift in gain!
                //
    
                
                if (ipeak != 0 ) {
                    // get the position where the peak should be... according to the first fit
                    e_start = e0*source_energy[ich][ipeak] / source_energy[ich][0];
                    _pk_tmp[ich]->GetXaxis()->SetRangeUser(e_start-100,e_start+100);
                } else {
                    // special care for channel 0 & channel 1: these tend to have a high background at low energy!
                    if (ich == 0 || ich ==1){
                        e_start = source_energy[ich][0];
                        _pk_tmp[ich]->GetXaxis()->SetRangeUser(e_start-100,e_start+100);
                    }
                }
                
                maxbin  = _pk_tmp[ich]->GetMaximumBin();
                maxval  = _pk_tmp[ich]->GetBinContent(maxbin);
                e_start = _pk_tmp[ich]->GetBinCenter(maxbin);
                
                
                _pk_tmp[ich]->GetXaxis()->SetRangeUser(0.,3000.);
                _pk_tmp[ich]->Draw();
                
                //
                // fit a Gauss + background to a photopeak
                //
                TF1 *func = new TF1("fit",fitf,e_start-200,e_start+200,5);
                Double_t res_start = 0.06/2.35*sqrt(662.)*sqrt(e_start);
                func->SetParameters(maxval,e_start,res_start);
                func->SetParNames("C","mean","sigma");
                
                Double_t demin = 100;
                Double_t demax = 100;

                if(ich == 4 || ich == 5) { // 60Co has some peaks close to each other
                    if(ipeak == 0) demax = 75;
                    if(ipeak == 1) demin = 75;
                } else if (ich == 0 || ich == 1){
                    demin = 200;
                    demax = 200;
                }
                
                
                Double_t e_low  = e_start - demin;
                Double_t e_high = e_start + demax;

                _pk_tmp[ich]->Fit("fit","Q","",e_low,e_high);
                
                Double_t peak        = func->GetParameter(0);
                _t_energy            = func->GetParameter(1);
                if(ipeak ==0) e0 = _t_energy;
                Double_t sigma       = func->GetParameter(2);
                _t_res = 0;
                if(_t_energy>0) _t_res = 2.355*sigma/_t_energy ;
                
                Double_t bin_width = (emax-emin)/nbin[ich];
                _t_rate = TMath::Sqrt(2*TMath::Pi())*sigma*peak / TIME_INTERVAL / bin_width;
                cout <<"get_interval_data:: ich ="<<ich<<" ipeak = "<<ipeak<<" E = "<<_t_energy<<" keV  rate = "<<_t_rate<<" Hz  resolution  = "<<_t_res<<" % "<<endl;
                
                _t_chanNum = (Int_t)ich;
                _t_peakNum = (Int_t)ipeak;
                
                tree->Fill();
                c1->Update();
//                                cin>>huh;
                
                delete func;
            }
        } // loop over peaks
        _pk_tmp[ich]->Reset(); // reset the histogram
    } // loop over channels
}

void analyzer::write_histograms(){
    //
    // write a few parameters to file
    //
    TNamed *Parameter = new TNamed("run",run.c_str());
    Parameter->Write();
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
        _e_err4[ich]->Write();
        
        _b_good[ich]->Write();
        _b_err1[ich]->Write();
        _b_err2[ich]->Write();
        _b_err4[ich]->Write();
        
        _2d_good[ich]->Write();
        _2d_err1[ich]->Write();
        _2d_err2[ich]->Write();
        _2d_err4[ich]->Write();
    }
    
    char hname[128];
    
    // historam as placeholder for drawing stability graphs
    TH1F *_holder = new TH1F("hld","hld",1,0,time_since_start);
    string htitle = "run: "+run;
    _holder->SetTitle(htitle.c_str());
    //_holder->Write();
    
    tree->Write();
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
