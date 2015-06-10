#define lifetime_cxx
#include "lifetime.h"
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>

const int nmax = 100000;

void lifetime::Loop(int channel_sel, int peak_sel)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   Double_t t[nmax],R[nmax],dR[nmax],tstart;
   Int_t n = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry == 0) tstart = t0;
      // if (Cut(ientry) < 0) continue;
      if(channel == channel_sel && peak == peak_sel){
         cout << jentry<< " " <<time << " " << rate <<endl;

         t[n] = t0+time-tstart;
         R[n] = rate;
         dR[n] = sqrt(rate*900)/900;
         
         n++; 
      }
   }
   cout << " n = "<<n<<endl;
   TGraphErrors *g1 = new TGraphErrors(n,t,R,0,dR);
   TF1 *f1 = new TF1("myfunc","[0]*exp(-x*[1])",0.,2e6);
   f1->SetParameters(32.5,4e-9);

   gStyle->SetOptFit(111);
   g1->SetMarkerStyle(24);
   g1->Fit("myfunc","","",400e3, 2000e3);
   g1->Draw("AP");
}
