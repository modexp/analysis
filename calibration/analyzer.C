#define analyzer_cxx
#include "analyzer.h"
// RooFit include files
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooExtendPdf.h>
//
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

using namespace RooFit;
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


/*----------------------------------------------------------------------------------------------------*/
Double_t fitf(Double_t *v, Double_t *par)
{
    Double_t arg = 0;
    if (par[2] != 0) arg = (v[0] - par[1])/par[2];
    
    Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
    fitval += par[3] + par[4]*v[0] + par[5]*v[0]*v[0];
    
    return fitval;
}
/*----------------------------------------------------------------------------------------------------*/
void analyzer::fit_spectrum(int ichannel){
    
//  TCanvas *c2 = new TCanvas("c2","c2",600,400);
  cout <<"analyzer::fit_spectrum  channel = "<<ichannel<<endl;

  // real
  RooRealVar E("E (keV)","E (keV)",emin,emax);

  // the BG function from a histogram
  string mc_file="";
  if       (ichannel == 2 || ichannel == 3){
      mc_file = "/user/z37/Modulation/analysis/RooTest/MC_ti44_modulation.root";
  } else if(ichannel == 4 || ichannel == 5){
      mc_file = "/user/z37/Modulation/analysis/RooTest/MC_co60_modulation.root";
  } else if(ichannel == 6 || ichannel == 7){
      mc_file = "/user/z37/Modulation/analysis/RooTest/MC_cs137_modulation.root";
  }

  TFile *f_mc = new TFile(mc_file.c_str(),"READONLY");
  TH1* h_bg  = (TH1*)f_mc->Get("h2");

  RooDataHist mc1("mc1","mc1",RooArgList(E),h_bg);
  f_mc->Close();
  RooHistPdf bg("bg","bg",E,mc1,0);

  // first Gauss for first photo-peak ....
  Double_t Eval = source_energy[ichannel][0];
  RooRealVar mean1("mean1","mean of gaussian 1",Eval,Eval-50,Eval+50);
  RooRealVar sigma1("sigma1","width of gaussians",25,5.,50.) ;
  RooRealVar g1frac("g1frac","fraction of gauss1",0.2,0.0,1.0) ;
  RooGaussian gauss1("gauss1","gaussian PDF",E,mean1,sigma1) ;

  // second Gauss....
  Eval = source_energy[ichannel][1];
  RooRealVar mean2("mean2","mean of gaussian 2",Eval,Eval-50,Eval+50);
  RooRealVar sigma2("sigma2","width of gaussians",25,5.,100.) ;
  RooRealVar g2frac("g2frac","fraction of gauss2",0.2,0.,1.0) ;
  RooGaussian gauss2("gauss2","gaussian PDF",E,mean2,sigma2) ;
  // third Gauss
  Eval = source_energy[ichannel][2];
  RooRealVar mean3("mean3","mean of gaussian 2",Eval,Eval-50,Eval+50);
  RooRealVar sigma3("sigma3","width of gaussians",25,5.,100.) ;
  RooRealVar g3frac("g3frac","fraction of gauss3",0.05,0.0,1.0) ;
  RooGaussian gauss3("gauss3","gaussian PDF",E,mean3,sigma3) ;

  // normalization
  RooRealVar Norm("Norm","Normalization",1e6,0.,1e12);

  // Gaussian functions for the extra photo-peaks + BG
  double fit_range[2];
  RooAddPdf *sum;
  if       (ichannel == 2 || ichannel ==3){
    // fit range
    fit_range[0] = 400;
    fit_range[1] = 1800;
    // pdf
    sum = new RooAddPdf("sum","g1+g2+g3+bg",RooArgList(gauss1,gauss2,gauss3,bg),RooArgList(g1frac,g2frac,g3frac));
  } else if(ichannel == 4 || ichannel ==5){
    // fit range
    fit_range[0] = 800;
    fit_range[1] = 2800;
    // pdf
    sum = new RooAddPdf("sum","g1+g2+g3+bg",RooArgList(gauss1,gauss2,gauss3,bg),RooArgList(g1frac,g2frac,g3frac));
  } else if(ichannel == 6 || ichannel == 7){
    // fit range
    fit_range[0] = 400;
    fit_range[1] = 1000;
    // pdf
    sum = new RooAddPdf("sum","g1+bg",RooArgList(gauss1,bg),RooArgList(g1frac));
  }
  E.setRange("signalRange",emin,emax);//fit_range[0],fit_range[1]);
  RooExtendPdf esum("esum","extended pdf with Norm",*sum,Norm,"signalRange");
    
  // get the data
  RooDataHist data("data","data",RooArgList(E),(TH1*)_pk_tmp[ichannel]);
    
  //
  // plot the data with the fittted function
  //
  // // RooPlot *Eframe = E.frame();
  // // data.plotOn(Eframe);

  //
  // fit the pdf to the data
  //
  RooFitResult *fr = esum.fitTo(data,Extended(kTRUE),Range(fit_range[0],fit_range[1]),Save());
  // // esum.plotOn(Eframe);
  // // esum.plotOn(Eframe,Components(gauss1),LineColor(2),LineWidth(2));
  // // if(ichannel == 2 || ichannel ==3 || ichannel ==4 || ichannel == 5){
  // //   esum.plotOn(Eframe,Components(gauss2),LineColor(2),LineWidth(2));
  // //   esum.plotOn(Eframe,Components(gauss3),LineColor(2),LineWidth(2));
  // // }
  // // esum.plotOn(Eframe,Components(bg),LineColor(kGreen),LineWidth(2));

  // // Eframe->Draw();
    
  //
  // process the variables to rates
  //
  processFitData(Norm,g1frac,mean1,sigma1,ichannel,0);
  if(ichannel == 2 || ichannel ==3 || ichannel ==4 || ichannel == 5){
    processFitData(Norm,g2frac,mean2,sigma2,ichannel,1);
    processFitData(Norm,g3frac,mean3,sigma3,ichannel,2);
  }
  //c2->Update();
    
  //int huh;
  //cin>>huh;
  //
  //
  //
  delete sum;
  delete fr;
}
  
void analyzer::processFitData(RooRealVar N, RooRealVar f, RooRealVar E, RooRealVar sig, int ichannel, int ipeak){
  //
  // process the fit data in order to get the rate with errors etc into the ntuple
  //
  Double_t E1   = E.getValV();
  Double_t R1   = f.getValV()*N.getValV()/TIME_INTERVAL;
  // correlations between fraction and N1 are very small (?)
  Double_t dR1  = sqrt(pow(N.getValV()*f.getError(),2)+pow(f.getValV()*N.getError(),2))/TIME_INTERVAL;
  Double_t res  = 2.355*sig.getValV()/E1;
  cout <<ichannel<<" "<<ipeak<<" "<<E1<<" "<<res<<" R = "<<R1<<" +- "<<dR1<<endl;
  addTreeEntry(E1,R1,dR1,res,ichannel,ipeak);
}

/*----------------------------------------------------------------------------------------------------*/
void analyzer::addTreeEntry(Double_t E, Double_t R, Double_t dR, Double_t res, Int_t ich, Int_t ipk){
    //
    // fill the tree with the fit results. 
    //
    // the slow data are processed elsewhere, but are also entering this tree:)
    //
    _t_energy = E;
    _t_rate   = R;
    _t_drate  = dR;
    _t_res    = res;
    _t_chanNum = ich;
    _t_peakNum = ipk;

    // this should be the only place where the fill command is called
    tree->Fill();
}
/*----------------------------------------------------------------------------------------------------*/
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
    //_T = new TH1F("T","T",1000,10+273.15,40+273.15);
    
    // Define tree and branches
    tree = new TTree("ana", "Analyzed spectra");
    
    tree->Branch("t0", &_t_t0, "t0/D");
    tree->Branch("time", &_t_time, "time/D");
    tree->Branch("channel", &_t_chanNum, "channel/I");
    tree->Branch("peak", &_t_peakNum, "peak/I");
    tree->Branch("rate", &_t_rate, "rate/D");
    tree->Branch("drate", &_t_drate, "drate/D");
    tree->Branch("e", &_t_energy, "e/D");
    tree->Branch("res", &_t_res, "res/D");
    tree->Branch("temp", &_t_temp, "temp/D");
    tree->Branch("pres", &_t_pres, "pres/D");
    tree->Branch("bx", &_t_bx, "bx/D");
    tree->Branch("by", &_t_by, "by/D");
    tree->Branch("bz", &_t_bz, "bz/D");
    tree->Branch("btot", &_t_btot, "btot/D");
    tree->Branch("humid", &_t_humid, "humid/D");
}
/*----------------------------------------------------------------------------------------------------*/
void analyzer::fill_histograms(){
    // fill all the histograms
    
    channel = channel % 100 ;
    //_T->Fill(temp+273.15);
    
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
/*----------------------------------------------------------------------------------------------------*/
void analyzer::get_interval_data(){
    //
    // analyze the data from a time interval.
    // length of interval set by TIME_INTERVAL which can be set in header file
    //
    
    // time
    _t_t0    = tstart;
    _t_time  = (t0+time_since_start)/2.;
    
    cout<<"analyzer::get_interval_data:: time_since_start ="<<time_since_start<<endl;
    
    for(int ich=0; ich<NUMBER_OF_CHANNELS; ich++){
        if(ich>1) {
            //
            // if we have a nice MC background model we use it to fit the spectrum
            //
            fit_spectrum(ich);
        } else {
            //
            // if we deal with the background detectors we use a linear fit + gauss to fit the signal
            //
            fit_spectrum_simple(ich);
        }
        _pk_tmp[ich]->Reset(); // reset the histogram
    } // loop over channels
}
/*----------------------------------------------------------------------------------------------------*/
void analyzer::fit_spectrum_simple(int ichannel){
    //    int huh;
//    TCanvas *c1 = new TCanvas("c1","c1",600,400);
//    int huh;
    //
    // find all the selected energy peaks
    //
    int      maxbin;
    double   maxval;
    Double_t e_start, e0, demin, demax;
    for (int ipeak=0; ipeak<MAX_PEAKS; ipeak++){
        if(source_energy[ichannel][ipeak] >0){
            //
            // find the fit starting values
            //
            
            //
            // first peak is special.... we use the GetMaximumBin() method in order
            // to find this peak even if there is a shift in gain!
            //
            
            
            if (ipeak != 0 ) {
                // get the position where the peak should be... according to the first fit
                e_start = e0*source_energy[ichannel][ipeak] / source_energy[ichannel][0];
                _pk_tmp[ichannel]->GetXaxis()->SetRangeUser(e_start-100,e_start+100);
            } else {
                // special care for channel 0 & channel 1: these tend to have a high background at low energy!
                if (ichannel == 0 || ichannel ==1){
                    e_start = source_energy[ichannel][0];
                    _pk_tmp[ichannel]->GetXaxis()->SetRangeUser(e_start-100,e_start+100);
                }
            }
            
            maxbin  = _pk_tmp[ichannel]->GetMaximumBin();
            maxval  = _pk_tmp[ichannel]->GetBinContent(maxbin);
            e_start = _pk_tmp[ichannel]->GetBinCenter(maxbin);
            
            _pk_tmp[ichannel]->GetXaxis()->SetRangeUser(0.,3000.);
 //           _pk_tmp[ichannel]->Draw();
            
            //
            // fit a Gauss + background to a photopeak
            //
            TF1 *func = new TF1("fit",fitf,e_start-200,e_start+200,5);
            Double_t res_start = 0.06/2.35*sqrt(662.)*sqrt(e_start);
            func->SetParameters(maxval,e_start,res_start);
            func->SetParNames("C","mean","sigma");
            
            Double_t demin = 100;
            Double_t demax = 100;
            
            if(ichannel == 4 || ichannel == 5) { // 60Co has some peaks close to each other
                if(ipeak == 0) demax = 75;
                if(ipeak == 1) demin = 75;
            } else if (ichannel == 0 || ichannel == 1){
                demin = 200;
                demax = 200;
            }
            
            
            Double_t e_low  = e_start - demin;
            Double_t e_high = e_start + demax;
            
            _pk_tmp[ichannel]->Fit("fit","Q","",e_low,e_high);
            
            Double_t peak        = func->GetParameter(0);
            _t_energy            = func->GetParameter(1);
            if(ipeak ==0) e0 = _t_energy;
            Double_t sigma       = func->GetParameter(2);
            _t_res = 0;
            if(_t_energy>0) _t_res = 2.355*sigma/_t_energy ;
            
            Double_t bin_width = (emax-emin)/nbin[ichannel];
            _t_rate = TMath::Sqrt(2*TMath::Pi())*sigma*peak / TIME_INTERVAL / bin_width;
            cout <<"get_interval_data:: ich ="<<ichannel<<" ipeak = "<<ipeak
                 <<" E = "<<_t_energy<<" keV  rate = "<<_t_rate<<" Hz  resolution  = "<<_t_res<<" % "<<endl;
            
            // fille the output tree.....
            addTreeEntry(_t_energy,_t_rate,1.0,_t_res,ichannel,ipeak);
//            c1->Update();
            
            delete func;
        } // loop over peaks
    }
}
/*----------------------------------------------------------------------------------------------------*/
void analyzer::write_histograms(){
    //
    // write a few parameters to file
    //
    _f->cd();
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
/*----------------------------------------------------------------------------------------------------*/
void analyzer::reset_interval_data(){
    
    n_interval = 0;
    
    _t_temp = 0;
    _t_pres = 0;
    _t_bx = 0;
    _t_by = 0;
    _t_bz = 0;
    _t_btot = 0;
    _t_humid = 0;
    
}
/*----------------------------------------------------------------------------------------------------*/
void analyzer::add_interval_data(){
    
    n_interval++;
    
    _t_temp += temp;
    _t_pres += pres;
    _t_bx += bx;
    _t_by += by;
    _t_bz += bz;
    _t_btot += btot;
    _t_humid += humid;
    
}
/*----------------------------------------------------------------------------------------------------*/
void analyzer::calculate_interval_data(){
    if(n_interval>0){
        _t_temp /= n_interval;
        _t_pres /= n_interval;
        _t_bx /= n_interval;
        _t_by /= n_interval;
        _t_bz /= n_interval;
        _t_btot /= n_interval;
        _t_humid /= n_interval;
    }
}
/*----------------------------------------------------------------------------------------------------*/
//
// MAIN:: Loop routine
//
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
        if (jentry == 0) {
            tstart = time;
            // reset all averages
            reset_interval_data();
        }
        time_since_start = time - tstart;
        if (jentry == 0) t0 = time_since_start;
        
        //
        // fill the monitoring histograms
        //
        fill_histograms();
        //
        // add the data for the slow control average calculations
        //
        add_interval_data();
        //
        // if we exceed the maximum time interval, get all the data recorded
        // during this time. then reset time for a new interval....
        //
        if(time_since_start - t0 > TIME_INTERVAL) {
            // calculate the slow control avarages
            calculate_interval_data();
            // fitting of the peaks
            get_interval_data();
            // reset the time for the start of the next interval
            t0 = time_since_start;
            // reset the interval data for calculating averages
            reset_interval_data();
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
