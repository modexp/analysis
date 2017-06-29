#define analyzer_cxx
/*---------------------------------------------------------------------------------------------------*/
//
// analyzer.C Routine to analyze spectra and calculate the rate of the sources
//
// Usage:
//  prompt> #include "analyzer.C"
//  prompt> analyzer ana(<directory_with_energy_calibrated_rootfiles>,<analysis_output_rootfile>)
//  prompt> ana.Loop()
//
// To inspect the fit results it is possible to plot the fit results. In analyzer.h  set:
// #define PLOT_ON_SCREEN 1
//
// To set the time interval in which a spectrum is calculated, in analyzer.h set:
// #define TIME_INTERVAL <interval_in_seconds>
//
// A.P. Colijn
/*---------------------------------------------------------------------------------------------------*/
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
#include <TMatrixDSym.h>
#include <iostream>
//#include <vector>
//#include <numeric>
#include <stdio.h>
#include <TMath.h>
#include <TF1.h>

#define INDEX_TEMPERATURE 24
#define MAX_INDEX 400

using namespace RooFit;
#define MAX_PEAKS 5

/*---------------------------------------------------------------------------------------------------*/
float source_energy[NUMBER_OF_SOURCES][MAX_PEAKS] =
//
// the energy peaks you wish to select for the analysis should be in this list
//
{
    {1460.,-1,-1,-1,-1},                 // ID0: Background
    {511.,1157.020,511.+1157.020,-1,-1}, // ID1: Ti44
    {1173.2,1332.5,1173.2+1332.5,-1,-1}, // ID2: Co60
    {661.7,-1,-1,-1,-1},                 // ID3: CS137
    {-1,-1,-1,-1,-1},                    // ID4: MN54
    {1460.,-1,-1,-1,-1}                  // ID5: K40
};

/*---------------------------------------------------------------------------------------------------*/

//
// ranges for plotting
//
const int   nbin0 = 600;
// number of bins for the temporary fit histograms.... channels 0+1 have few entries so wider bins
float nbin[NUMBER_OF_CHANNELS]={nbin0,nbin0,nbin0,nbin0,nbin0,nbin0,nbin0,nbin0};
const float emin = 0.; // in keV
const float emax = 3000.; // in keV
const float adc_max_volt = 2.;
const float base_max_val = 2000;

/*----------------------------------------------------------------------------------------------------*/

//
// Gaussian function + 2nd order polynomial for simple rate fitting
//
Double_t fitf(Double_t *v, Double_t *par)
{
    Double_t arg = 0;
    if (par[2] != 0) arg = (v[0] - par[1])/par[2];
    
    Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
    fitval += par[3] + par[4]*v[0] + par[5]*v[0]*v[0];
    
    return fitval;
}
/*----------------------------------------------------------------------------------------------------*/


int plot_counter = 0; //Temporary fix by Joran (jorana@nikhef.nl)
void analyzer::fit_spectrum(int ichannel, double *fit_range){
    //
    // RooFit based spectrum fitter
    //
    cout <<"analyzer::fit_spectrum  channel = "<<ichannel<<endl;
    
    // get the source ID
    int id = source_id[ichannel];
    // identify the peaks in our specified energy range and store their peak ID.
    //
    int nselect = 0;
    int peak_id[MAX_PEAKS];
    for(int ipeak = 0; ipeak<MAX_PEAKS; ipeak++){
        double epeak = source_energy[id][ipeak];
        if ( (epeak > fit_range[0]) && (epeak < fit_range[1])){ // yes the peak is in range
            peak_id[nselect] = ipeak;
            nselect++;
        }
    }
    cout <<"analyzer::fit_spectrum found " <<nselect<<" peaks from "<<fit_range[0]
    <<" keV < E < "<<fit_range[1]<<" keV"<<endl;
    
    //
    // no peaks have been found..... leave fitting routine
    //
    if(nselect == 0) return;
    
    //
    // spectrum is a function of the energy
    //
    RooRealVar E("E","E (keV)",emin,emax);
    //
    // get the data from a TH1 root histogram
    //
    RooDataHist data("data","data",RooArgList(E),(TH1*)_pk_tmp[ichannel]);
    

    //
    // the MC template for each of the sources obtained from a GEANT4 simulation
    //
    string mc_file="";
    if       (id == TI44){
        mc_file = "MC_ti44_modulation.root";
    } else if(id == CO60){
        mc_file = "MC_co60_modulation.root";
    } else if(id == CS137){
        mc_file = "MC_cs137_modulation.root";
    } else if(id == MN54){
        mc_file = "MC_mn54_modulation.root";
    } else if(id == K40){
        mc_file = "MC_k40_modulation.root";
    } else {
        cout <<"analyzer::fit_spectrum BAD source identifier"<<endl;
    }
    
    cout <<"analyzer::fit_spectrum channel = "<<ichannel<<" source_id = "<<id<<" MC template ="<<mc_file<<endl;
    
    TFile *f_mc = new TFile(mc_file.c_str(),"READONLY");
    TH1* h_bg  = (TH1*)f_mc->Get("h2"); 

    //
    // construct the MC pdf
    //
    RooDataHist mc1("mc1","mc1",RooArgList(E),h_bg);
    delete h_bg;
    f_mc->Close();
    RooHistPdf bg("bg","bg",E,mc1,0);

    //=====================================================================================================================
    // Include the measurement the background rate into the fitting procedure. Added by Joran (jorang@xs4all.nl)
    // The BG template should be custom made 
    // http://www.physics.purdue.edu/darkmatters/doku.php?id=modulation:daq:bgtemplatefit

    string ch0_file = "ch_0.root";

    cout <<"analyzer::fit_spectrum channel = "<<ichannel<<" source_id = "<<id<<" BG template ="<<ch0_file<<endl;

    TFile *f_ch0 = new TFile(ch0_file.c_str(),"READONLY");
    
    TH1* h_time  = (TH1*)f_ch0->Get("time");
    double time_template = h_time->GetBinContent(1);

    //
    // construct the background pdf
    //
    TH1* h_ch0  = (TH1*)f_ch0->Get("ch0");  
    RooDataHist bg_data("bg_data","bg_data",RooArgList(E), h_ch0);
    delete h_time; delete h_ch0;
    f_ch0->Close();


    // We could also open the other  background channel but one will be sufficient for now.
    // string ch1_file = "ch_1.root";

    // cout <<"fit_spectrum:: channel = "<<ichannel<<" source_id = "<<id<<" BG template ="<<ch1_file<<endl;

    // TFile *f_ch1 = new TFile(ch1_file.c_str(),"READONLY");
    // TH1* h_ch1   = (TH1*)f_ch1->Get("ch1");  
    

    // RooDataHist mc3("mc3","mc3",RooArgList(E), h_ch1);
    // f_ch1->Close();

    RooHistPdf bg_ch0("bg_ch0","bg_ch0", E, bg_data, 0);
    
    double data_tot_events = 0; double bg_data_tot_events = 0; double mc3_tot_events = 0;
    for (int binnum = fit_range[0]/((emax-emin)/nbin0); binnum < fit_range[1]/((emax-emin)/nbin0); binnum++){
        bg_data.get(binnum) ; double events_bg_data = bg_data.weight() / time_template; // double devents_bg_data = bg_data.weightError() / time_template;
        bg_data_tot_events += events_bg_data;
        data.get(binnum); double events_data = data.weight() / delta_t    ; // double devents_data = data.weightError() / delta_t;
        data_tot_events += events_data;
    }

    std::vector<RooRealVar*> bg_pdf_frac;    
    double bg_frac = (bg_data_tot_events) / (data_tot_events);

    // // A simple background fraction can also be used (as below) but we if we fit the spectrum over a large energy range it is better to not constrain too tightly
    // bg_pdf_frac.push_back(new RooRealVar("bg_pdf_frac","bg_pdf_frac",0.07,0.001,0.50));
    bg_pdf_frac.push_back(new RooRealVar("bg_pdf_frac","bg_pdf_frac", bg_frac, 0.0, 1.0));

    // 
    //=====================================================================================================================
    //
    // define the Gaussians for the photo peaks
    //
    cout <<"analyzer::fit_spectrum Define photo peak Gaussians"<<endl;
    std::vector<RooRealVar*> pk_mean;
    std::vector<RooRealVar*> pk_sigma;
    std::vector<RooRealVar*> pk_frac;
    std::vector<RooGaussian*> pk_gaus;
    
    std::vector<RooRealVar*> pk_frac_tail;
    std::vector<RooRealVar*> pk_sigma_tail;
    std::vector<RooGaussian*> pk_gaus_tail;
    
    
    char vname[128];
    for (int isel=0; isel<nselect; isel++){
        int ipeak = peak_id[isel];
        double epeak = source_energy[id][ipeak];
        cout <<"analyzer::fit_spectrum:: ichannel = "<<ichannel <<" ipeak = "<<ipeak<<" epeak = "<<epeak<<endl;
        sprintf(vname,"mean%i",isel);
        pk_mean.push_back(new RooRealVar(vname,vname,epeak,epeak-50,epeak+50));
        sprintf(vname,"sigma%i",isel);
        pk_sigma.push_back(new RooRealVar(vname,vname,25,5,100));
        sprintf(vname,"frac%i",isel);
        pk_frac.push_back(new RooRealVar(vname,vname,0.3- 0.1 * isel, 0.0,1.0)); // change 0.3 1 - expected Compton - expected background
        
        sprintf(vname,"gaus%i",isel);
        pk_gaus.push_back(new RooGaussian(vname,vname,E,*pk_mean[isel],*pk_sigma[isel]));
        //
        sprintf(vname,"frac_tail%i",isel);
        pk_frac_tail.push_back(new RooRealVar(vname,vname,0.01,0.0,0.1));
        sprintf(vname,"sigma_tail%i",isel);
        pk_sigma_tail.push_back(new RooRealVar(vname,vname,50,5,100));
        sprintf(vname,"gaus_tail%i",isel);
        pk_gaus_tail.push_back(new RooGaussian(vname,vname,E,*pk_mean[isel],*pk_sigma_tail[isel]));
        
    }
    cout <<"analyzer::fit_spectrum Define photo peak Gaussians ---- DONE"<<endl;
    //
    // normalization for the pdf will be a separately fitted parameter
    //
    RooRealVar Norm("Norm","Normalization",1e6,0.,1e12);    
    
    //
    // compose the joined pdf for background + signal
    //
    
    cout <<"analyzer::fit_spectrum Compose the combined pdf"<<endl;
    
    RooAddPdf *sum;
    // Changed by Joran (jorana@nikhef.nl) for including ch0 and ch1 as BG
    if(FIT_BG_TEMPLATE){
        if (nselect==1)  sum = new RooAddPdf("sum","g1+bg_ch0+bg",      RooArgList(*pk_gaus[0], bg_ch0, bg), RooArgList(*pk_frac[0], *bg_pdf_frac[0])); 
        if (nselect==2)  sum = new RooAddPdf("sum","g1+g2+bg_ch0+bg",   RooArgList(*pk_gaus[0],*pk_gaus[1], bg_ch0,bg),RooArgList(*pk_frac[0],*pk_frac[1], *bg_pdf_frac[0]));
        if (nselect==3)  sum = new RooAddPdf("sum","g1+g2+g3+bg_ch0+bg",RooArgList(*pk_gaus[0],*pk_gaus[1], *pk_gaus[2], bg_ch0 ,bg),RooArgList(*pk_frac[0],*pk_frac[1],*pk_frac[2], *bg_pdf_frac[0]));
    } else{
        if (nselect==1)  sum = new RooAddPdf("sum","g1+bg",RooArgList(*pk_gaus[0],bg),RooArgList(*pk_frac[0]));
        if (nselect==2)  sum = new RooAddPdf("sum","g1+g2+bg",RooArgList(*pk_gaus[0],*pk_gaus[1],bg),RooArgList(*pk_frac[0],*pk_frac[1]));
        if (nselect==3)  sum = new RooAddPdf("sum","g1+g2+g3+bg",RooArgList(*pk_gaus[0],*pk_gaus[1],*pk_gaus[2],bg),RooArgList(*pk_frac[0],*pk_frac[1],*pk_frac[2]));
    }
    
    cout <<"analyzer::fit_spectrum Compose the combined pdf ---- DONE "<<endl;
    
    E.setRange("signalRange",emin,emax);//fit_range[0],fit_range[1]);
    RooExtendPdf esum("esum","extended pdf with Norm",*sum,Norm,"signalRange"); 
    
    cout <<"analyzer::fit_spectrum Compose the combined pdf ---- extended DONE "<<endl;
    
    //
    // fit the pdf to the data
    //
    
    fr = esum.fitTo(data,Extended(kTRUE),Range(fit_range[0],fit_range[1]),Save());
    fr->Print();
    //
    // process the fitted variables and store in the output tree
    //
    //============================
    // Added lines below to get a chi2 without plotting
    if       (PLOT_ON_SCREEN == 0){
        for (int id=0; id<nselect; id++){
            //
            // covariance matrix elements for calculation of error on rate
            //
            RooPlot *Eres_frame = E.frame();
            Eres_frame->GetXaxis()->SetRangeUser(fit_range[0],fit_range[1]);
                       
                    
            //
            // plot the data and pdfs. use the plot range as found before
            //
            data.plotOn(Eres_frame);
            esum.plotOn(Eres_frame);

            // Get the residuals 
            cout <<"analyzer::fit_spectrum chi2 for channel = "<<ichannel<< " chi2 = "<< Eres_frame->chiSquare() <<endl;
                        
            Double_t chindf2 = Eres_frame->chiSquare() ;
            if(IGNORE_SMALL_DATASET) {
                if (int(TIME_INTERVAL)==int(delta_t)) processFitData(Norm,*pk_frac[id],*pk_mean[id],*pk_sigma[id], ichannel, peak_id[id], chindf2, *bg_pdf_frac[0]);
            }
            else {
                processFitData(Norm,*pk_frac[id],*pk_mean[id],*pk_sigma[id], ichannel, peak_id[id], chindf2, *bg_pdf_frac[0]);
            }
            delete Eres_frame;            
        }
    }
    //============================
    // for (int id=0; id<nselect; id++){
    //     //
    //     // covariance matrix elements for calculation of error on rate
    //     //
    //     processFitData(Norm,*pk_frac[id],*pk_mean[id],*pk_sigma[id], ichannel,peak_id[id]);
    // }
    
    //
    // draw the fit results..... not for batch running. Set PLOT_ON_SCREEN variable to 0 in analyzer.h file
    //
    if(PLOT_ON_SCREEN){
        TCanvas *c1 = new TCanvas("c1","c1",1800,1200);
        
        //
        // plot the data with the fitted function
        //
        RooPlot *Eframe = E.frame();
        Eframe->SetTitle("");
        Eframe->GetXaxis()->SetRangeUser(fit_range[0],fit_range[1]);
                
        //
        // find maximum data value in plot range
        //
        _pk_tmp[ichannel]->GetXaxis()->SetRangeUser(fit_range[0],fit_range[1]);
        Int_t maxbin  = _pk_tmp[ichannel]->GetMaximumBin();
        Double_t maxval  = _pk_tmp[ichannel]->GetBinContent(maxbin);
        
        //
        // plot the data and pdfs. use the plot range as found before
        //
        data.plotOn(Eframe);
        esum.plotOn(Eframe);

        // Get the residuals by using residHist() and hpull ()
        cout <<"chi2 for channel = "<<ichannel<<endl;
        Double_t chindf = Eframe->chiSquare() ;
        cout << chindf <<endl;


        RooHist *hresid =  Eframe->residHist();
        RooHist* hpull = Eframe->pullHist();

        // pull distribution
        char cmd[128];
        sprintf(cmd, "Pull distribution. Channel = %i Fit range = %i - %i", int(ichannel), int(fit_range[0]), int(fit_range[1]));
        TH1F *pull_dis = new TH1F(cmd, "pull", int(50), -10 * sqrt(HOURS), 10 * sqrt(HOURS));
        
        // Fill the pull distribution form the entries in the pull histogram
        double res, dres;
        for(int i=0; i<(fit_range[1]-fit_range[0])/((emax-emin)/nbin0); i++){             
            res = hpull->GetY()[i];
            dres= hpull->GetErrorY(i);
            pull_dis->Fill(res,dres);
        }
        
        pull_dis->SetTitle(cmd);
        pull_dis->GetYaxis()->SetTitle("frequency");
        pull_dis->GetXaxis()->SetTitle("Residuals/error");
        pull_dis->SetLineColor(1);
        pull_dis->SetMarkerColor(1);
        pull_dis->SetMarkerStyle(8);


        // Create a new frames to draw the residuals and the pull distribution
        RooPlot* res_frame = E.frame(Title("Residuals")) ;
        RooPlot* hpull_frame = E.frame(Title("Residuals / error")) ;
        
        // res_frame->SetTitle("");
        res_frame->addPlotable(hresid,"P") ;
        res_frame->GetYaxis()->SetTitle("Residuals (events)");
        res_frame->GetXaxis()->SetRangeUser(fit_range[0], fit_range[1]);
        
        // hpull_frame->SetTitle("");   
        hpull_frame->addPlotable(hpull,"P") ;
        hpull_frame->GetYaxis()->SetTitle("Residuals");     
        hpull_frame->GetXaxis()->SetRangeUser(fit_range[0], fit_range[1]);

        RooPlot* spec_full_frame = E.frame(Title("Whole spectrum")) ;
        spec_full_frame->SetTitle("");
        spec_full_frame->GetXaxis()->SetRangeUser(emin,emax);
        data.plotOn(spec_full_frame);
        esum.plotOn(spec_full_frame);
      
        // Show the chi2 and the run name on the first frame
        string runstr = run.substr (5,8);
        TPaveLabel *t0 = new TPaveLabel(0.1,1.0,0.9,0.9, Form("%s", runstr.c_str()), "brNDC"); 
        TPaveLabel *t1 = new TPaveLabel(0.1,0.88,0.3,0.80, Form("#chi^{2} = %.2f", chindf), "brNDC"); 
        // Double_t frac = *pk_frac[id].getValV();
        // Double_t R1   = Norm*frac/delta_t;
        // R1 -> paramOn(Eframe,Layout(0.35,0.88,0.85-id*0.15));
        // TPaveLabel *t3 = new TPaveLabel(0.1,0.88,0.3,0.80, Form("#chi^{2} = %.3f", R1), "brNDC"); 
        // Eframe -> addObject(t0);
        Eframe -> addObject(t1);
        spec_full_frame -> addObject(t0);
        // spec_full_frame -> addObject(t1);

        for (int id=0; id<nselect; id++){
            esum.plotOn(Eframe,Components(*pk_gaus[id]),LineColor(2),LineWidth(2));
            pk_gaus[id]->paramOn(Eframe,Layout(0.55,0.88,0.85-id*0.15));
            esum.plotOn(spec_full_frame,Components(*pk_gaus[id]),LineColor(2),LineWidth(2));

            Double_t frac =  pk_frac[id]->getVal();

            Double_t R1   = (Norm.getValV() )*frac/delta_t;
            TPaveLabel *t2 = new TPaveLabel(0.65,0.88-id*0.1,0.88,0.80-id*0.1, Form("R_{%i} = %.2f Hz", id, R1), "brNDC"); 
            spec_full_frame -> addObject(t2);
        }

        Double_t bg_frac_fit=  bg_pdf_frac[0] -> getVal();
        Double_t bg_R    = (Norm.getValV() )* bg_frac_fit / delta_t;
        Double_t bg_R_exp= (Norm.getValV() )* bg_frac     / delta_t;
        TPaveLabel *t3 = new TPaveLabel(0.65,0.88-nselect*0.1, 0.88, 0.80-nselect*0.1, Form("BG Rate fit = %.1f, expect = %.1f Hz", bg_R, bg_R_exp), "brNDC"); 
        spec_full_frame -> addObject(t3);
        
        // For the HV measurements I wanted to add a label to the plot might be useful for other purposes too
        // if (ichannel == 2){
        //     TPaveLabel *t4 = new TPaveLabel(0.1,0.88-0.1,0.3, 0.88,  Form("HV = %.2f V", hv2), "brNDC");
        //     spec_full_frame -> addObject(t4);   
        // }else if (ichannel == 3){ 
        //     TPaveLabel *t4 = new TPaveLabel(0.1,0.88-0.1,0.3, 0.88,  Form("HV = %.2f V", hv3), "brNDC"); 
        //     spec_full_frame -> addObject(t4);      
        // }else if (ichannel == 4){ 
        //     TPaveLabel *t4 = new TPaveLabel(0.1,0.88-0.1,0.3, 0.88,  Form("HV = %.2f V", hv4), "brNDC");
        //     spec_full_frame -> addObject(t4);      
        // }else if (ichannel == 5){ 
        //     TPaveLabel *t4 = new TPaveLabel(0.1,0.88-0.1,0.3, 0.88,  Form("HV = %.2f V", hv5), "brNDC");
        //     spec_full_frame -> addObject(t4);
        // }else if (ichannel == 6){ 
        //     TPaveLabel *t4 = new TPaveLabel(0.1,0.88-0.1,0.3, 0.88,  Form("HV = %.2f V", hv6), "brNDC");
        //     spec_full_frame -> addObject(t4);      
        // }else if (ichannel == 7){ 
        //     TPaveLabel *t4 = new TPaveLabel(0.1,0.88-0.1,0.3, 0.88,  Form("HV = %.2f V", hv7), "brNDC");
        //     spec_full_frame -> addObject(t4);       
        // }
        
        esum.plotOn(Eframe,Components(bg),LineColor(kGreen),LineWidth(2));
        esum.plotOn(spec_full_frame,Components(bg),LineColor(kGreen),LineWidth(2));

        if(FIT_BG_TEMPLATE){
            esum.plotOn(Eframe,Components(bg_ch0),LineColor(kCyan),LineWidth(2));
            esum.plotOn(spec_full_frame,Components(bg_ch0),LineColor(kCyan),LineWidth(2));
        }

        // I save the plots produced somewhere, for that i also use a counter (to keep track of the number of plots). This is a temporary fix.
        plot_counter = plot_counter +1; //     

        char plotpath[128];
        
        sprintf(plotpath, "/data/xenon/mod_admin/data/anaplots/chan%i/E%i%sch%i_E%i_%i_num%i.png", int(ichannel),int (fit_range[0]), run.c_str(), int(ichannel),int (fit_range[0]), int (fit_range[1]), int(plot_counter)); 
        
        c1->Divide(2,2) ;
        c1->cd(1) ; gPad->SetLogy();Eframe->SetMaximum(100000*HOURS); Eframe->SetMinimum(1); gPad->SetBottomMargin(0.15) ; Eframe->GetYaxis()->SetTitleOffset(1.1) ; Eframe->Draw() ;
        c1->cd(2) ; gPad->SetLogy();spec_full_frame->SetMaximum(100000*HOURS); spec_full_frame->SetMinimum(1); gPad->SetBottomMargin(0.15) ; spec_full_frame->GetYaxis()->SetTitleOffset(1.1) ; spec_full_frame->Draw() ;
        c1->cd(3) ; gPad->SetBottomMargin(0.15) ; hpull_frame->SetMaximum(10*sqrt(HOURS)); hpull_frame->SetMinimum(-10*sqrt(HOURS)); hpull_frame->GetYaxis()->SetTitleOffset(1.1) ; hpull_frame->Draw() ;        
        c1->cd(4) ; pull_dis->Fit("gaus", "Q"); gStyle->SetOptFit(0110); pull_dis->SetMinimum(0); pull_dis->SetMaximum((fit_range[1]-fit_range[0])/(((emax-emin)/nbin0)*10));pull_dis->Draw("E1")  ;
 
        if(IGNORE_SMALL_DATASET) { 
            if (int(TIME_INTERVAL)==int(delta_t)) c1->SaveAs(plotpath);
        } else {
            c1->SaveAs(plotpath);
        }
        
        c1->Close();
        
        delete pull_dis; delete Eframe; delete res_frame; delete hpull_frame; delete spec_full_frame;
        
        for (int id=0; id<nselect; id++){         
            if(IGNORE_SMALL_DATASET) {
                if (int(TIME_INTERVAL)==int(delta_t)) processFitData(Norm,*pk_frac[id],*pk_mean[id],*pk_sigma[id], ichannel, peak_id[id], chindf, *bg_pdf_frac[0]);
            } else { 
                processFitData(Norm,*pk_frac[id],*pk_mean[id],*pk_sigma[id], ichannel, peak_id[id], chindf, *bg_pdf_frac[0]);
            }
        }        
        
    }
    //
    // cleanup
    //
    delete sum;
    delete fr;
    
}
// =======================================================================================================================
// =======================================================================================================================
// Lets do some fitting here for channels 0 and 1

void analyzer::fit_spectrum_background(int ichannel, double *fit_range){
    //
    // RooFit based spectrum fitter
    //
    cout <<"analyzer::fit_spectrum_background channel = "<<ichannel<<endl;
    
    // get the source ID
    int id = source_id[ichannel];
        
    //
    // spectrum is a function of the energy
    //
    RooRealVar E("E","E (keV)",emin,emax);
    
    
    //=====================================================================================================================
    // Include the measurement of channel 0 into the fitting procedure. Added by Joran (jorana@nikhef.nl)
    // Loading the file    
    string ch0_file = "ch_0.root";
   
    TFile *f_ch0 = new TFile(ch0_file.c_str(),"READONLY");
    TH1* h_ch0  = (TH1*)f_ch0->Get("ch0");  

    //
    // construct the background pdf
    //
    RooDataHist bg_data("bg_data","bg_data",RooArgList(E), h_ch0);
    delete h_ch0;
    f_ch0->Close();
    RooHistPdf bg_ch0("bg_ch0","bg_ch0", E, bg_data, 0);
    // Introduce the relative weight of this PDF that is fitted
    std::vector<RooRealVar*> bg_pdf_frac;
    
    //
    // normalization for the pdf will be a sperately fitted parameter
    //
    RooRealVar Norm("Norm","Normalization",1e6,0.,1e12);
    
    
    //
    // compose the joined pdf for background + signal
    //
    RooAddPdf *sum;
    sum = new RooAddPdf("sum","bg_pdf_frac",RooArgList(bg_ch0),RooArgList());
    
    E.setRange("signalRange",emin,emax);//fit_range[0],fit_range[1]);
    RooExtendPdf esum("esum","extended pdf with Norm",*sum,Norm,"signalRange"); 
    
    cout <<"analyzer::fit_spectrum_background Compose the combined pdf ---- extended DONE "<<endl;
    //
    // get the data from a TH1 root histogram
    //
    RooDataHist data("data","data",RooArgList(E),(TH1*)_pk_tmp[ichannel]);


    // fit the pdf to the data
    //
    fr = esum.fitTo(data,Extended(kTRUE),Range(fit_range[0],fit_range[1]),Save());
    fr->Print();
    //
    // process the fitted variables and store in the output tree
    //
    //============================
    if       (PLOT_ON_SCREEN == 0){
        //
        // Draw the data and fit on two frames and calculate the chi^2
        //
        RooPlot *Eres_frame = E.frame();
        Eres_frame->GetXaxis()->SetRangeUser(fit_range[0],fit_range[1]);
                                    
        data.plotOn(Eres_frame);
        esum.plotOn(Eres_frame);

        // Get the residuals 
        cout <<"analyzer::fit_spectrum_background chi2 for channel = "<<ichannel<<" ch2 = " << Eres_frame->chiSquare() <<endl;
        
        Double_t chindf2 = Eres_frame->chiSquare() ;
        if(IGNORE_SMALL_DATASET) {
            if (int(TIME_INTERVAL)==int(delta_t)) processFitData_BackGround(Norm, ichannel, chindf2);
        } else {
            processFitData_BackGround(Norm, ichannel, chindf2);
        }
    }
    // ============================
    
    if(PLOT_ON_SCREEN){
        TCanvas *c1 = new TCanvas("c1","c1",1800,1200);
        
        //
        // plot the data with the fittted function
        //
        RooPlot *Eframe = E.frame();
        Eframe->SetTitle("");
        Eframe->GetXaxis()->SetRangeUser(fit_range[0],fit_range[1]);
                
        //
        // find maximum data value in plot range
        //
        _pk_tmp[ichannel]->GetXaxis()->SetRangeUser(fit_range[0],fit_range[1]);
        Int_t maxbin  = _pk_tmp[ichannel]->GetMaximumBin();
        Double_t maxval  = _pk_tmp[ichannel]->GetBinContent(maxbin);
        
        //
        // plot the data and pdfs. use the plot range as found before
        //
        data.plotOn(Eframe);
        esum.plotOn(Eframe);
        // Get the residuals by using residHist() and hpull ()
        Double_t chindf = Eframe->chiSquare() ;
        cout << chindf <<endl;
        cout << Eframe->chiSquare() <<endl;
        RooHist *hresid =  Eframe->residHist();
        RooHist* hpull = Eframe->pullHist();


        char cmd[128];
        sprintf(cmd, "Pull distribution. Channel = %i Fit range = %i - %i", int(ichannel), int(fit_range[0]), int(fit_range[1]));
        TH1F *pull_dis = new TH1F(cmd,"pull",int(sqrt(HOURS)*50),-10*sqrt(HOURS),sqrt(HOURS)*10);
        double res, dres;
        for(int i=0; i<(fit_range[1]-fit_range[0])/((emax-emin)/nbin0); i++){
            res = hpull->GetY()[i];
            dres = hpull->GetErrorY(i);
            pull_dis->Fill(res,dres);
        }
        
        pull_dis->SetTitle(cmd);
        pull_dis->GetYaxis()->SetTitle("frequency");
        pull_dis->GetXaxis()->SetTitle("pull");
        pull_dis->SetLineColor(1);
        pull_dis->SetMarkerColor(1);
        pull_dis->SetMarkerStyle(8);
        

        // Create a new frames to draw the residuals and the pull distribution
        RooPlot* res_frame = E.frame(Title("Residual Distribution")) ;
        RooPlot* hpull_frame = E.frame(Title("Residuals / error")) ;
        RooPlot* spec_full_frame = E.frame(Title("Whole spectrum")) ;
        
        // res_frame->SetTitle("");
        res_frame->addPlotable(hresid,"P") ;
        res_frame->GetYaxis()->SetTitle("Residual (events)");
        res_frame->GetXaxis()->SetRangeUser(fit_range[0],fit_range[1]);
        
        // hpull_frame->SetTitle("");    
        hpull_frame->addPlotable(hpull,"P") ;
        hpull_frame->GetYaxis()->SetTitle("Residuals");      
        hpull_frame->GetXaxis()->SetRangeUser(fit_range[0],fit_range[1]);
        
        // chindf -> paramOn(Eframe, Layout(0.55,0.08,0.85-id*0.15));
        TPaveLabel *t0 = new TPaveLabel(0.1,1.0,0.9,0.9, Form("%s", run.c_str()), "brNDC"); 
        TPaveLabel *t1 = new TPaveLabel(0.1,0.88,0.3,0.80, Form("#chi^{2} = %.2f", chindf), "brNDC"); 
        Eframe -> addObject(t0);
        Eframe -> addObject(t1);
        Double_t R1   = (Norm.getValV() )/delta_t;
        TPaveLabel *t2 = new TPaveLabel(0.65,0.88,0.88,0.80, Form("R = %.2f Hz", R1), "brNDC"); 
        spec_full_frame -> addObject(t2);

        // For the HV measurements I wanted to add something to the plot might be usefull for other purposes too
        // if (ichannel == 0) {
        //     TPaveLabel *t4 = new TPaveLabel(0.1,0.88-0.1,0.3, 0.88, Form("HV = %.2f V", hv0), "brNDC");
        //     spec_full_frame -> addObject(t4);
        // } else if (ichannel == 1) {
        //     TPaveLabel *t4 = new TPaveLabel(0.1,0.88-0.1,0.3, 0.88, Form("HV = %.2f V", hv1), "brNDC"); 
        //     spec_full_frame -> addObject(t4);
        // }

        spec_full_frame->GetXaxis()->SetRangeUser(emin,emax);
        data.plotOn(spec_full_frame);
        esum.plotOn(spec_full_frame);
        esum.plotOn(spec_full_frame,Components(bg_ch0),LineColor(kCyan),LineWidth(2));
        esum.plotOn(Eframe,Components(bg_ch0),LineColor(kCyan),LineWidth(2));
     
        // I save the plots produced somewhere, for that i also use a counter (to keep track of the number of plots). This is a temporary fix.
        plot_counter = plot_counter +1; //
      
        char plotpath[128];
        sprintf(plotpath, "/data/xenon/mod_admin/data/anaplots/chan%i/E%i%sch%i_E%i_%i_num%i.png",int(ichannel),int (fit_range[0]), run.c_str(), int(ichannel),int (fit_range[0]), int (fit_range[1]), int(plot_counter)); 

        c1->Divide(2,2) ;
        c1->cd(1) ; gPad->SetLogy();Eframe->SetMaximum(HOURS*1000); Eframe->SetMinimum(1); gPad->SetBottomMargin(0.15) ; Eframe->GetYaxis()->SetTitleOffset(1.1) ; Eframe->Draw() ;
        c1->cd(3) ; gPad->SetBottomMargin(0.15) ; hpull_frame->SetMaximum(sqrt(HOURS)*10); hpull_frame->SetMinimum(-10*sqrt(HOURS)); hpull_frame->GetYaxis()->SetTitleOffset(1.1) ;        hpull_frame->Draw() ;
        c1->cd(2) ; gPad->SetLogy();spec_full_frame->SetMaximum(HOURS*1000); spec_full_frame->SetMinimum(1); gPad->SetBottomMargin(0.15) ; spec_full_frame->GetYaxis()->SetTitleOffset(1.1) ; spec_full_frame->Draw() ;
        c1->cd(4) ; pull_dis->Fit("gaus", "Q"); gStyle->SetOptFit(0110); pull_dis->SetMinimum(0); pull_dis->SetMaximum((fit_range[1]-fit_range[0])/(((emax-emin)/nbin0)*5));pull_dis->Draw("E1")  ;

        if(IGNORE_SMALL_DATASET) {
            if (int(TIME_INTERVAL)==int(delta_t)) c1->SaveAs(plotpath);
        } else {
            c1->SaveAs(plotpath);
        }

        c1->Close();
        delete pull_dis; delete Eframe; delete res_frame; delete hpull_frame; delete spec_full_frame;

        //int huh;
        //cin>>huh;

        if(IGNORE_SMALL_DATASET) {
            if (int(TIME_INTERVAL)==int(delta_t)) processFitData_BackGround(Norm, ichannel, chindf);
        } else {
            processFitData_BackGround(Norm, ichannel, chindf);
        }
    }
    //
    // cleanup
    //
    delete sum;
    delete fr;
}



// =======================================================================================================================
// =======================================================================================================================
/*----------------------------------------------------------------------------------------------------*/
void analyzer::fit_spectrum(int ichannel){
    //
    // RooFit based spectrum fitter
    //
    //
    cout <<"analyzer::fit_spectrum  channel = "<<ichannel<<endl;
    
    // get the souce id
    int id = source_id[ichannel];
    
    //
    // spectrum is a function of the energy
    //
    RooRealVar E("E","E (keV)",emin,emax);
    
    //
    // the background template for each of the sources obtained from a GEANT4 simulation
    //
    string mc_file="";
    if       (id == TI44){
        mc_file = "/user/z37/Modulation/analysis/calibration/MC_ti44_modulation.root";
    } else if(id == CO60){
        mc_file = "/user/z37/Modulation/analysis/calibration/MC_co60_modulation.root";
    } else if(id == CS137){
        mc_file = "/user/z37/Modulation/analysis/calibration/MC_cs137_modulation.root";
    } else if(id == MN54){
        mc_file = "/user/z37/Modulation/analysis/calibration/MC_mn54_modulation.root";
    } else if(id == K40){
        mc_file = "/user/z37/Modulation/analysis/calibration/MC_k40_modulation.root";
    } else {
        cout <<"fit_spectrum:: BAD source identifier"<<endl;
    }
    
    cout <<"fit_spectrum:: channel = "<<ichannel<<" source_id = "<<id<<" MC template ="<<mc_file<<endl;
    
    TFile *f_mc = new TFile(mc_file.c_str(),"READONLY");
    TH1* h_bg  = (TH1*)f_mc->Get("h2");
    //
    // construct the background pdf
    //
    RooDataHist mc1("mc1","mc1",RooArgList(E),h_bg);
    f_mc->Close();
    RooHistPdf bg("bg","bg",E,mc1,0);
    
    //
    // first Gauss for first photo-peak ....
    //
    Double_t Eval = source_energy[id][0];
    RooRealVar mean1("mean1","mean of gaussian 1",Eval,Eval-50,Eval+50);
    RooRealVar sigma1("sigma1","width of gaussians",25,5.,50.) ;
    RooRealVar g1frac("g1frac","fraction of gauss1",0.2,0.0,1.0) ;
    RooGaussian gauss1("gauss1","gaussian PDF",E,mean1,sigma1) ;
    // second Gauss....
    Eval = source_energy[id][1];
    RooRealVar mean2("mean2","mean of gaussian 2",Eval,Eval-50,Eval+50);
    RooRealVar sigma2("sigma2","width of gaussians",25,5.,100.) ;
    RooRealVar g2frac("g2frac","fraction of gauss2",0.2,0.,1.0) ;
    RooGaussian gauss2("gauss2","gaussian PDF",E,mean2,sigma2) ;
    // third Gauss
    Eval = source_energy[id][2];
    RooRealVar mean3("mean3","mean of gaussian 3",Eval,Eval-50,Eval+50);
    RooRealVar sigma3("sigma3","width of gaussians",25,5.,100.) ;
    RooRealVar g3frac("g3frac","fraction of gauss3",0.05,0.0,1.0) ;
    RooGaussian gauss3("gauss3","gaussian PDF",E,mean3,sigma3) ;
    
    //
    // normalization for the pdf will be a sperately fitted parameter
    //
    RooRealVar Norm("Norm","Normalization",1e6,0.,1e12);
    
    //
    // compose the joined pdf for background + signal
    //
    double fit_range[2];
    RooAddPdf *sum;
    if       (ichannel == 2 || ichannel == 3 ){
        // fit range
        fit_range[0] = 400;
        fit_range[1] = 1800;
        // pdf
        sum = new RooAddPdf("sum","g1+g2+g3+bg",RooArgList(gauss1,gauss2,gauss3,bg),RooArgList(g1frac,g2frac,g3frac));
    } else if(ichannel == 4 || ichannel == 5 ){
        // fit range
        fit_range[0] = 800;
        fit_range[1] = 2800;
        // pdf
        sum = new RooAddPdf("sum","g1+g2+g3+bg",RooArgList(gauss1,gauss2,gauss3,bg),RooArgList(g1frac,g2frac,g3frac));
    } else if(ichannel == 6 || ichannel == 7 ){
        // fit range
        fit_range[0] = 400;
        fit_range[1] = 1000;
        // pdf
        sum = new RooAddPdf("sum","g1+bg",RooArgList(gauss1,bg),RooArgList(g1frac));
    }
    E.setRange("signalRange",emin,emax);//fit_range[0],fit_range[1]);
    RooExtendPdf esum("esum","extended pdf with Norm",*sum,Norm,"signalRange");
    
    //
    // get the data from a TH1 root histogram
    //
    RooDataHist data("data","data",RooArgList(E),(TH1*)_pk_tmp[ichannel]);
    //
    // fit the pdf to the data
    //
    //RooFitResult *fr = esum.fitTo(data,Extended(kTRUE),Range(fit_range[0],fit_range[1]),Save());
    fr = esum.fitTo(data,Extended(kTRUE),Range(fit_range[0],fit_range[1]),Save());
    //
    // process the variables to rates
    //
    processFitData(Norm,g1frac,mean1,sigma1,ichannel,0, 99, g1frac ); //last 99 is chindf and last g1frac should be ignored
    if(ichannel == 2 || ichannel ==3 || ichannel ==4 || ichannel == 5){
        processFitData(Norm,g2frac,mean2,sigma2,ichannel,1, 99, g1frac);//last 99 is chindf and last g1frac should be ignored
        processFitData(Norm,g3frac,mean3,sigma3,ichannel,2, 99, g1frac);//last 99 is chindf and last g1frac should be ignored
    }
    //
    // cleanup
    //
    delete sum;
    delete fr;
}
/*----------------------------------------------------------------------------------------------------*/
void analyzer::processFitData_BackGround(RooRealVar N, int ichannel, Double_t chi2ndfs){
    //
    // process the fit data in order to get the rate with errors etc into the ntuple
    //
    cout<<"analyzer::processFitData_BackGround"<<endl;
    Double_t Norm = N.getValV();
    Double_t chi2s = chi2ndfs;
    Double_t R1   = Norm/delta_t;
    
    Double_t dR1 = Norm;
    dR1 = sqrt(dR1)/delta_t;
    
    cout <<ichannel<<" "<<" "<<" "<<" R = "<<R1<<" +- "<<dR1<<endl;
    
    addTreeEntry(0,R1,dR1,0,ichannel,0, 1, chi2s, 0.0, 1.0); 
}
/*----------------------------------------------------------------------------------------------------*/
void analyzer::processFitData(RooRealVar N, RooRealVar f, RooRealVar E, RooRealVar sig, int ichannel, int ipeak, Double_t chi2ndfs, RooRealVar bg_f){
    //
    // process the fit data in order to get the rate with errors etc into the ntuple
    //
    Double_t E1   = E.getValV();
    Double_t Norm = N.getValV();
    Double_t frac = f.getValV();
    
    Double_t chi2s = chi2ndfs;
    //    Double_t R1   = Norm*frac/TIME_INTERVAL;
    Double_t R1   = Norm*frac/delta_t;
    //
    // Get covariance matrix elements. We calculate Rate = frac*Norm / delta_t, so
    //
    // dRate = sqrt(frac**2*cov(0,0) + 2*frac*Norm*cov(idx,0) + Norm**2*cov(idx,idx))/delta_t
    //
    // The idx should point to the right element in the covariance matrix:
    //     frac0 -> idx=1
    //     frac1 -> idx=2
    //     frac2 -> idx=3
    //     etc etc
    //
    // Be careful to change the order of the elements in the covariance matrix: maybe there is a way to reference
    // by name as well.....
    //
    int idx = 0;
    string fName = f.GetName();

    //    cout <<">>"<<fName<<"<<"<<endl;

    if        (fName == "frac0") {
        idx = 1;
    } else if (fName == "frac1"){
        idx = 2;
    } else if (fName == "frac2"){
        idx = 3;
    }

    Double_t bg_R  = 0;
    Double_t bg_dR = 1; 
    // Calculate the error on the background rate
    if (FIT_BG_TEMPLATE){
        Double_t bg_frac = bg_f.getValV();
        bg_R = Norm * bg_frac / delta_t ;
        bg_dR = Norm * Norm * covariance(1, 1);
        bg_dR += 2 * Norm * bg_frac * covariance(1, 0);
        bg_dR +=  bg_frac * bg_frac * covariance(0, 0);
        bg_dR = sqrt(bg_dR) / delta_t;
        // cout << "got here"<<bg_frac<<bg_R<<" +/- "<<bg_dR<<endl;
        // Since there is one extra item in the covariance matrix (the BG fraction) the proper covariance is one index lower.
        idx = idx + 1; 
    }

    //   cout <<" cov00 = "<<covariance(0,0)<<" cov01 = "<<covariance(idx,0)<<" cov11 = "<<covariance(idx,idx)<<endl;
    Double_t dR1 = Norm*Norm*covariance(idx,idx);
    dR1 += 2*Norm*frac*covariance(idx,0);
    dR1 +=   frac*frac*covariance(0,0);
    dR1 = sqrt(dR1)/delta_t;
    //
    // calculate error on the rate
    //


    Double_t res  = 2.355*sig.getValV()/E1;
    Double_t fractry1 = frac; 
    cout <<"analyzer::processFitData -- channel "<<ichannel<<", pk "<<ipeak<<", E "<<E1<<" keV, res= "<<res<<" R = "<<R1<<" +- "<<dR1<<"Hz, R_BG = "<<bg_R<<"+/-"<< bg_dR<<endl;
 
    addTreeEntry(E1, R1, dR1, res, ichannel, ipeak, fractry1, chi2s, bg_R, bg_dR);
    cout << "analyzer::processFitData -- Done" << endl;
    
}

/*----------------------------------------------------------------------------------------------------*/
void analyzer::addTreeEntry(Double_t E, Double_t R, Double_t dR, Double_t res, Int_t ich, Int_t ipk, Double_t fractry, Double_t chindfin, Double_t bg_rate, Double_t bg_drate){
    //
    // fill the tree with the fit results.
    //
    // the slow data are processed elsewhere, but are also entering this tree:)
    //
    _t_energy  = E;
    _t_rate    = R;
    _t_bg_rate = bg_rate;
    _t_drate   = dR;
    _t_bg_drate= bg_drate;
    _t_res     = res;
    _t_chanNum = ich;
    _t_peakNum = ipk;
    _t_fractry = fractry;
    _t_chi2    = chindfin;

    // To calculate the total number of errors/time I keep a list of these.
    tot_tot_err[ich] = _e_err1[ich] -> GetEntries();
    tot_tot_err[ich] = _e_err2[ich] -> GetEntries();
    tot_tot_err[ich] = _e_err4[ich] -> GetEntries();

    _t_err_rate = (tot_tot_err[ich] - glob_tot_err[ich])/delta_t;
    glob_tot_err[ich] = tot_tot_err[ich];
    cout<<"analyzer::addTreeEntry\n\n analyzer::addTreeEntry rate of error" <<_t_err_rate<<endl;
    
    // this should be the only place where the fill command is called
    tree->Fill();
    
}
/*----------------------------------------------------------------------------------------------------*/
double analyzer::covariance(int i, int j){
    TMatrixDSym cov = fr->covarianceMatrix();
    double cc = cov[i][j];
    return cc;
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
        
        tot_tot_err.push_back(0);
        glob_tot_err.push_back(0);
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
    tree->Branch("bgrate", &_t_bg_rate, "bg_rate/D");
    tree->Branch("bgdrate", &_t_bg_drate, "bg_drate/D");
    tree->Branch("e", &_t_energy, "e/D");
    tree->Branch("res", &_t_res, "res/D");
    tree->Branch("temp", &_t_temp, "temp/D");
    tree->Branch("pres", &_t_pres, "pres/D");
    tree->Branch("bx", &_t_bx, "bx/D");
    tree->Branch("by", &_t_by, "by/D");
    tree->Branch("bz", &_t_bz, "bz/D");
    tree->Branch("btot", &_t_btot, "btot/D");
    tree->Branch("humid", &_t_humid, "humid/D");

    tree->Branch("frac", &_t_fractry, "frac/D"); //joranadd
    tree->Branch("chi2ndf", &_t_chi2, "chi2/D"); //joranadd
    tree->Branch("tot_error", &_t_err_rate, "tot_error/D");

    tree->Branch("hv0", &_t_hv0, "hv0/D");
    tree->Branch("hv1", &_t_hv1, "hv1/D");
    tree->Branch("hv2", &_t_hv2, "hv2/D");
    tree->Branch("hv3", &_t_hv3, "hv3/D");
    tree->Branch("hv4", &_t_hv4, "hv4/D");
    tree->Branch("hv5", &_t_hv5, "hv5/D");
    tree->Branch("hv6", &_t_hv6, "hv6/D");
    tree->Branch("hv7", &_t_hv7, "hv7/D");
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
        // _pk_tmp[channel]->Fill(integral);
        _e_err1[channel]->Fill(integral);
        _b_err1[channel]->Fill(baseline);
        
        _2d_err1[channel]->Fill(integral,height);
    }
    else if ((error&0x02)!=0) {
        // _pk_tmp[channel]->Fill(integral);
        _e_err2[channel]->Fill(integral);
        _b_err2[channel]->Fill(baseline);
        
        _2d_err2[channel]->Fill(integral,height);
    }
    else if ((error&0x04)!=0) {
        // _pk_tmp[channel]->Fill(integral);
        // temporary error handling: err4 is not a real error yet
        _e_err4[channel]->Fill(integral);
        _b_err4[channel]->Fill(baseline);
        
        _2d_err4[channel]->Fill(integral,height);
        // _pk_tmp[channel]->Fill(integral);
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
    
    double range[2] = {0,3000};
    for(int ich=0; ich<NUMBER_OF_CHANNELS; ich++){
        
        // only analyze active channels
        if(channel_active[ich]){
            
            //
            // if we have a nice MC background model we use it to fit the spectrum
            // if we don't have a nice model we will do a simple fit instead
            //
            
            // which is the source?
            int id = source_id[ich];
            if        (id == TI44) {
                range[0] = 400; range[1] = 2800;
                fit_spectrum(ich, range);
            } else if (id == CO60 ) {
                range[0] = 900; range[1] = 3000;
                fit_spectrum(ich, range);
            } else if (id == CS137 ) {
                range[0] = 500; range[1] = 1200;
                fit_spectrum(ich, range);
            } else if (id == MN54) {
                range[0] = 0; range[1] = 2000;
                fit_spectrum(ich, range);
            } else if (id == K40) {
                range[0] = 1300; range[1] = 1600;
                fit_spectrum(ich, range);
            } else if        (id == BACKGROUND) {
                range[0] = 600; range[1] = 3000;
                fit_spectrum_background(ich, range);
            }
        }
        _pk_tmp[ich]->Reset(); // reset the histogram
    } // loop over channels
}
/*----------------------------------------------------------------------------------------------------*/
void analyzer::fit_spectrum_simple(int ichannel){
    //    int huh;
    //    TCanvas *c1 = new TCanvas("c1","c1",600,400);
    //    int huh;
    
    // source id
    int id = source_id[ichannel];
    
    //
    // find all the selected energy peaks
    //
    int      maxbin;
    double   maxval;
    Double_t e_start, e0, demin, demax;
    for (int ipeak=0; ipeak<MAX_PEAKS; ipeak++){
        if(source_energy[id][ipeak] >0){
            //
            // find the fit starting values
            //
            
            //
            // first peak is special.... we use the GetMaximumBin() method in order
            // to find this peak even if there is a shift in gain!
            //
            
            
            if (ipeak != 0 ) {
                // get the position where the peak should be... according to the first fit
                e_start = e0*source_energy[id][ipeak] / source_energy[id][0];
                _pk_tmp[ichannel]->GetXaxis()->SetRangeUser(e_start-100,e_start+100);
            } else {
                // special care for channel 0 & channel 1: these tend to have a high background at low energy!
                if (ichannel == 0 || ichannel ==1){
                    e_start = source_energy[id][0];
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
            _t_rate = TMath::Sqrt(2*TMath::Pi())*sigma*peak / delta_t / bin_width;
            cout <<"get_interval_data:: ich ="<<ichannel<<" ipeak = "<<ipeak
            <<" E = "<<_t_energy<<" keV  rate = "<<_t_rate<<" Hz  resolution  = "<<_t_res<<" % "<<endl;
            
            // fille the output tree.....
            addTreeEntry(_t_energy,_t_rate,1.0,_t_res,ichannel,ipeak, _t_fractry, 0.0, 0.0, 1.0); //0 for chi2ndf, 0.0 for R_bg and 1.0 for dR_bg
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
    
    //
    _t_hv0 = 0;
    _t_hv1 = 0;
    _t_hv2 = 0;
    _t_hv3 = 0;
    _t_hv4 = 0;
    _t_hv5 = 0;
    _t_hv6 = 0;
    _t_hv7 = 0;
    
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
    
    //
    _t_hv0 += hv0;
    _t_hv1 += hv1;
    _t_hv2 += hv2;
    _t_hv3 += hv3;
    _t_hv4 += hv4;
    _t_hv5 += hv5;
    _t_hv6 += hv6;
    _t_hv7 += hv7;
    
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

        // 
        _t_hv0 /= n_interval;
        _t_hv1 /= n_interval;
        _t_hv2 /= n_interval;
        _t_hv3 /= n_interval;
        _t_hv4 /= n_interval;
        _t_hv5 /= n_interval;
        _t_hv6 /= n_interval;
        _t_hv7 /= n_interval;
    }
}
/*---------------------------------------------------------------------------------------------------*/
void analyzer::get_source_id()
{
    cout <<"analyzer::get_source_id"<<endl;
    
    char channel_name[100];
    
    // get the name of the first file in the data chain
    TFile * _f_tmp = fChain->GetFile();
    _f_tmp->cd("info/active");
    
    TNamed *isActive;
    for(int ichannel=0; ichannel<NUMBER_OF_CHANNELS; ichannel++){
        
        sprintf(channel_name,"channel_%i",ichannel);
        gDirectory->GetObject(channel_name,isActive);
        string active = isActive->GetTitle();
        if(active == "On"){
            channel_active[ichannel] = kTRUE;
        } else {
            channel_active[ichannel] = kFALSE;
        }
        
    }
    
    // retrieve the source information
    _f_tmp->cd("info/source");
    
    TNamed *sourceName;
    for(int ichannel=0; ichannel<NUMBER_OF_CHANNELS; ichannel++){
        sprintf(channel_name,"channel_%i",ichannel);
        gDirectory->GetObject(channel_name,sourceName);
        string source = sourceName->GetTitle();
        
        if(channel_active[ichannel]){
            
            cout <<"ecal::get_source_id  channel = "<<ichannel<<" source = "<<source<<endl;
            if(source == "Background"){
                source_id[ichannel] = BACKGROUND;
            } else if ( source == "Ti-44"){
                source_id[ichannel] = TI44;
            } else if ( source == "Co-60"){
                source_id[ichannel] = CO60;
            } else if ( source == "Cs-137"){
                source_id[ichannel] = CS137;
            } else if ( source == "Mn-54"){
                source_id[ichannel] = MN54;
            } else if ( source == "K-40"){
                source_id[ichannel] = K40;
            } else {
                // source_id[ichannel] = TI44; // Do need to change Joran
                cout <<"ecal::get_source_id() Unidentified source ..... TERMINATE"<<endl;
                exit(-1);
            }
        } else {// if the channel is inactive, just set it to BACKGROUND (does not matter)
            source_id[ichannel] = BACKGROUND;
        }
    }
    cout <<"analyzer::get_source_id ... done"<<endl;
    
}
/*----------------------------------------------------------------------------------------------------*/
//
// MAIN:: Loop routine
//
void analyzer::Loop()
{
    if (fChain == 0) return;
    
    //
    // look in the first data file of the chain to see what sources are present
    //
    get_source_id();
    
    
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
    Bool_t last_event = false;
    
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        //
        // get entry from the tree
        //
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) last_event = true;
        
        if(!last_event){
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
        }
        delta_t = time_since_start - t0;
        if((delta_t > TIME_INTERVAL) || last_event) {
            // calculate the slow control averages
            calculate_interval_data();
            // fitting of the peaks
            get_interval_data();
            // reset the time for the start of the next interval
            t0 = time_since_start;
            // reset the interval data for calculating averages
            reset_interval_data();
        }
        
        if(jentry%500000 == 0) cout<<"Processed "<<jentry<<" events"<<endl;
        
        if(last_event) break;
    }
    
    cout<<"Done event loop...."<<endl;
    
    //
    // write histograms to file and generate the stability graphs
    //
    write_histograms();
    cout<<"Done writing histograms"<<endl;
}
