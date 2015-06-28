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

void lifetime::Loop(int channel_sel, int peak_sel)
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
   cout << " n = "<<n<<endl;
   TGraphErrors *g1 = new TGraphErrors(n,t,R,0,dR);
   TF1 *f1 = new TF1("myfunc","[0]*exp(-x*0.6931471805599/[1]/3600/24/365)",0.,4e6);
   //TF1 *f1 = new TF1("myfunc","[0]*(1-x*[1])",0.,1000e6);
   f1->SetParameters(10,10);

   gStyle->SetOptFit(111);
   g1->SetMarkerStyle(24);
//   g1->Fit("myfunc","","",300e3, 3000e6);
   g1->Fit("myfunc");
   g1->Draw("AP");

   char cmd[128];
   sprintf(cmd,"Rate as a function of time. Channel = %i Photo-peak = %i",channel_sel,peak_sel);
   g1->SetTitle(cmd);
   g1->GetXaxis()->SetTitle("time (sec)");
   g1->GetYaxis()->SetTitle("rate (Hz)");


   sprintf(cmd,"t12_ch%i_pk%i.png",channel_sel,peak_sel);
   c1->Print(cmd);
  
}
