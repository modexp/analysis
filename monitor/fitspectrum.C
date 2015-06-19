#include <iostream>
#include <string>

gSystem->Load("libRooFit") ;
using namespace RooFit ;

const double emin = 0.;
const double emax = 2999.;

#define NUMBER_OF_CHANNELS 8
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

void fitspectrum(string data_file, string mc_file, int ichannel){
  //
  // fit spectrum to a MC - BG function + Gaussian photopeaks
  //
  cout <<"Welcome to the Modulation spectrum fitter"<<endl;
  RooRealVar E("E (keV)","E (keV)",emin,emax);

  // the BG function from a histogram
  TFile *f_mc = new TFile(mc_file.c_str(),"READONLY");
  TH1* h_bg  = (TH1*)f_mc->Get("h2");
  RooDataHist mc1("mc1","mc1",RooArgList(E),h_bg);
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
  RooRealVar Norm("Norm","Number of events in 1st Gauss",1e6,0.,1e12);

  // Gaussian functions for the extra photo-peaks + BG
  double fit_range[2];
  if       (ichannel == 2 || ichannel ==3){ 
    // fit range
    fit_range[0] = 400;
    fit_range[1] = 1800;
    // pdf
    RooAddPdf sum("sum","g1+g2+g3+bg",RooArgList(gauss1,gauss2,gauss3,bg),RooArgList(g1frac,g2frac,g3frac));
  } else if(ichannel == 4 || ichannel ==5){ 
    // fit range
    fit_range[0] = 800;
    fit_range[1] = 2800;
    // pdf
    RooAddPdf sum("sum","g1+g2+g3+bg",RooArgList(gauss1,gauss2,gauss3,bg),RooArgList(g1frac,g2frac,g3frac));
  } else if(ichannel == 6 || ichannel == 7){ 
    // fit range
    fit_range[0] = 400;
    fit_range[1] = 1000;
    // pdf
    RooAddPdf sum("sum","g1+bg",RooArgList(gauss1,bg),RooArgList(g1frac));
  }
  E.setRange("signalRange",emin,emax);//fit_range[0],fit_range[1]);
  RooExtendPdf esum("esum","extended pdf with Norm",sum,Norm,"signalRange");


  //
  // get the data
  //
  TFile *f_data = new TFile(data_file.c_str(),"READONLY");
  char hname[128];
  sprintf(hname,"_e_good_ch%i",ichannel);
  TH1* h_data  = (TH1*)f_data->Get(hname);
  RooDataHist data("data","data",RooArgList(E),h_data);

  Double_t dt = runtime->GetVal();
  //
  // plot the data with the fittted function
  //  
  RooPlot *Eframe = E.frame();
  data.plotOn(Eframe);

  RooFitResult *fr = esum.fitTo(data,Extended(kTRUE),Range(fit_range[0],fit_range[1]),Save());
  fr->Print();
  g1frac.Print();
  cout << " g1 frac = "<<g1frac.getValV()<<endl;
  esum.plotOn(Eframe);
  esum.plotOn(Eframe,Components(gauss1),LineColor(2),LineWidth(2));
  if(ichannel == 2 || ichannel ==3 || ichannel ==4 || ichannel == 5){
    esum.plotOn(Eframe,Components(gauss2),LineColor(2),LineWidth(2));
    esum.plotOn(Eframe,Components(gauss3),LineColor(2),LineWidth(2));
  }
  esum.plotOn(Eframe,Components(bg),LineColor(kGreen),LineWidth(2));

  Eframe->Draw();

  Double_t R1   = g1frac.getValV()*Norm.getValV()/dt;
  cout << " R1 = "<<R1<<endl;

  cout <<"Bye bye ..."<<endl;
}
