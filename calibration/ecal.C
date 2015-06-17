#define ecal_cxx
#include "ecal.h"
#include <TF1.h>
#include <TMath.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TParameter.h>
#include <iostream>
#include <vector>
#include <stdio.h>

/*---------------------------------------------------------------------------------------------------*/
// (1) You need to identify where the first peak is in the histograms of the integrals.
// (2) Set range_low and range_high to find the range in which the peak with cal_energy should be
const double range_low[NUMBER_OF_CHANNELS] ={0.16e-6,0.14e-6,0.05e-6,0.05e-6,0.,0.,0.,0.};
const double range_high[NUMBER_OF_CHANNELS]={1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6};

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
    {661.7,-1,-1,-1,-1}, // channel6: 137Cs
    {661.7,-1,-1,-1,-1}  // channel7: 137Cs
};

/*---------------------------------------------------------------------------------------------------*/
Double_t fitf_gauss(Double_t *v, Double_t *par)
{
    Double_t arg = 0;
    if (par[2] != 0) arg = (v[0] - par[1])/par[2];
    
    Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
    
    fitval += par[3]+par[4]*v[0];
    
    return fitval;
}
/*---------------------------------------------------------------------------------------------------*/
void ecal::book_histograms(){
    //
    // book histograms for energy calibration
    //
    cout <<"ecal::book_histograms"<<endl;
    _f = new TFile(calFile.c_str(),"RECREATE");
    
    char tmp[100];
    for (int ich = 0; ich <NUMBER_OF_CHANNELS; ich++){
        sprintf(tmp,"integral_ch%02d",ich);
        _integral.push_back(new TH1F(tmp,tmp,1000,0.,1e-6));
        sprintf(tmp,"energy_ch%02d",ich);
        _energy.push_back(new TH1F(tmp,tmp,1000,0.,3000.));
        sprintf(tmp,"energy_all_ch%02d",ich);
        _energy_all.push_back(new TH1F(tmp,tmp,1000,0.,3000.));
        
    }
    cout <<"ecal::book_histograms ... done"<<endl;
}
/*---------------------------------------------------------------------------------------------------*/
void ecal::fill_histograms(int ilevel){
    
    Long64_t nentries = fChain->GetEntriesFast();
    cout<<"Start calibration loop.... nentries ="<<nentries<<" LEVEL = "<<ilevel<<endl;
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // fill the hisotgrams with 'good' events
        channel = channel % 100;
        
        if        (ilevel == BEFORE_CALIBRATION){
            if(error == 0) _integral[channel]->Fill(integral);
        } else if (ilevel == AFTER_CALIBRATION){
            if(error == 0) _energy[channel]->Fill(ccal[channel][0]+integral*ccal[channel][1]);
            _energy_all[channel]->Fill(ccal[channel][0]+integral*ccal[channel][1]);
        }
        if(jentry%100000 == 0) cout<<"Processed "<<jentry<<" events"<<endl;
        
    }
    cout<<"Done calibration loop...."<<endl;
    
}
/*---------------------------------------------------------------------------------------------------*/
void ecal::Loop()
{
    //
    // energy calibration for modulation detectors
    //
    if (fChain == 0) return;
    //
    // book the histograms needed for the energy calibration
    //
    book_histograms();
    //
    // loop over the events and fill histograms: before calibration
    //
    fill_histograms(BEFORE_CALIBRATION);
    
    TCanvas *c1 = new TCanvas("c1","c1",600,400);
    int huh;
    char parname[100];
    for(int ich=0; ich<NUMBER_OF_CHANNELS; ich++){
        cout <<"CHANNEL = "<<ich<<endl;
        // find the number of peaks to calibrate on.....
        int npeak = 0;
        for (int i=0; i<MAX_PEAKS; i++){
            if(source_energy[ich][i]>0) npeak++;
        }
        
        if(npeak == 0){ // no calibration possible ..... E = 1.0*integral.
            ccal[ich][0] = 0;
            ccal[ich][1] = 1;
            ccal[ich][2] = 0;
        }
        
        // STEP1: find the highest peak
        // define the proper range for searching the calibration constant
        _integral[ich]->GetXaxis()->SetRangeUser(range_low[ich],range_high[ich]);
        // find the bin with maximum value in the range
        int ibin      = _integral[ich]->GetMaximumBin();
        double val    = _integral[ich]->GetBinCenter(ibin);
        double maxval = _integral[ich]->GetBinContent(ibin);

        
        Double_t ee[MAX_PEAKS];
        Double_t dee[MAX_PEAKS];
        Double_t area[MAX_PEAKS];
        Double_t darea[MAX_PEAKS];
        Double_t mean_integrals[MAX_PEAKS];
        
        Double_t vlow,vhigh;
        for(int ipeak = 0; ipeak<npeak; ipeak++){
            //
            // fit ranges: if you are the first peak, then we use as start value for
            // the mean the maximum bin as found above. otherwise we scale the
            // starting point from the energy peak we found before in the first fit
            //
            if (ipeak != 0) val = area[0]*source_energy[ich][ipeak]/source_energy[ich][0];
            vlow  = val - 0.035e-6;
            vhigh = val + 0.035e-6;
            
            if(channel == 4 || channel ==5){
                if(ipeak == 0 ) vhigh = val+0.015e-6;
                if(ipeak == 1 ) vlow  = val-0.015e-6;
            }
            
            //
            // fit a Gauss around the desired location in the spectrum
            //

            // starting value for the height
            _integral[ich]->GetXaxis()->SetRangeUser(vlow,vhigh);
            ibin   = _integral[ich]->GetMaximumBin();
            maxval = _integral[ich]->GetBinContent(ibin);
            _integral[ich]->GetXaxis()->SetRangeUser(0.,1.);

            
            // do the fit
            TF1 *func = new TF1("fit",fitf_gauss,vlow,vhigh,4);
            func->SetParameters(maxval,val,0.01e-7);
            func->SetParNames("C","mean","sigma");
            _integral[ich]->Fit("fit","Q","",vlow,vhigh);
            
            // get the parameters
            ee[ipeak]     = source_energy[ich][ipeak];
            dee[ipeak]    = 0.;
            area[ipeak]   = func->GetParameter(1);
            darea[ipeak]  = func->GetParError(1);
        }
        
        //
        // * calculate the calibration parameters.
        // * distinct case when fit is possible (npeak > 1)and when it is not (npeak == 1).
        //
        if (npeak == 1){ // just one point... so the energy calibration is a simple proportionality
            ccal[ich][0] = 0;
            ccal[ich][1] = source_energy[ich][0] / area[0];
            ccal[ich][2] = 0;
        } else if(npeak > 1){ // more than one point... now we fit a straight line
            //
            // Fit a 1st order polynomial to the integral vs energy (inverse calibration). Done
            // like this because I have the errors on the area, and the energy points are
            // infinitely precise - just the photo peak energy values.
            //
            TGraphErrors *gcal = new TGraphErrors(npeak,ee,area,0,darea);
            gcal->Fit("pol1");
            TF1 *ff = gcal->GetFunction("pol1");
            // get the parameters.... but remember I have fitted teh inverse calibration here
            // so now I need to invert as well....
            //   integral = v0+v1*E  -> E = -v0/v1 + 1/v1 * integral
            Double_t v0 = ff->GetParameter(0);
            Double_t v1 = ff->GetParameter(1);
            // construct the calibration constants
            ccal[ich][0] = -v0/v1;
            ccal[ich][1] = 1/v1;
            ccal[ich][2] = 0;
        }
        
        for(int ipar=0; ipar<MAX_PARAMETERS; ipar++){
            sprintf(parname,"cal_ch%02d_c%i",ich,ipar);
            TParameter <double> *p1 = new TParameter<double>(parname,ccal[ich][ipar]);
            p1->Write();
        }
    }
    
    //
    // loop over the events and fill histograms: post-calibration
    //
    fill_histograms(AFTER_CALIBRATION);
    
    
    _f->Write();
    _f->Close();
}
