#define lifetime_cxx
#include "lifetime.h"
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>

#define LN2 0.69314718055994530942

const int nmax = 100000;

void lifetime::Life(int channel_sel, int peak_sel, string type, bool save)
{
   TCanvas *c1 = new TCanvas("c1","c1",600,400);
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   Double_t t[nmax],R[nmax],dR[nmax],tstart;
   Int_t n = 0,icounter=0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry == 0) tstart = t0;
      // if (Cut(ientry) < 0) continue;
      cout << jentry<<endl;
      if(channel == channel_sel && peak == peak_sel ){
           cout << jentry<< " " <<time << " " << rate <<endl;

           t[n] = t0+time-tstart;
           R[n] = rate;
           dR[n] = drate;//sqrt(rate*900)/900;
         
           n++; 
      }
   }

   //
   // make a TGraphErrors object and fit an exponential
   //
   TGraphErrors *g1 = new TGraphErrors(n,t,R,0,dR);
   TF1 *f1 = new TF1("myfunc","[0]*exp(-x*0.6931471805599/[1]/3600/24/365)");//,0.,4e6);
   //TF1 *f1 = new TF1("myfunc","[0]*(1-x*[1])",0.,1000e6);
   f1->SetParameters(10,10);
   g1->Fit("myfunc");

   TH1F *_pull = new TH1F("pull","pull",50,-5,5);
   //
   // plot the results
   //
   char cmd[128];
   gStyle->SetOptFit(111);
   if(type == "life"){
     g1->SetMarkerStyle(24);
     g1->Draw("AP");

     sprintf(cmd,"Rate as a function of time. Channel = %i Photo-peak = %i",channel_sel,peak_sel);
     g1->SetTitle(cmd);
     g1->GetXaxis()->SetTitle("time (sec)");
     g1->GetYaxis()->SetTitle("rate (Hz)");

   } else if (type == "pull"){
     double res;
     for(int i=0; i<n; i++){
       res = (R[i] - f1->Eval(t[i]) )/dR[i];
       _pull->Fill(res);
     }
     sprintf(cmd,"Pull distribution. Channel = %i Photo-peak = %i",channel_sel,peak_sel);
     _pull->SetTitle(cmd);
     _pull->GetXaxis()->SetTitle("pull");
     _pull->SetLineColor(4);
     _pull->SetMarkerColor(4);
     _pull->SetMarkerStyle(20);
     _pull->Fit("gaus");
     _pull->Draw("pe");
   }

   if(save){
     sprintf(cmd,"t12_ch%i_pk%i.png",channel_sel,peak_sel);
     c1->Print(cmd);
   }
  
}
