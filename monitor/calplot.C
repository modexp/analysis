#include <string>

float chi2_cut[8] = {1.8,2.5,10.,2.5,10.0,2.5,3.5,5.5};

void calplot(string cal_file, string plot_type, int ich){
   char cmd[128],cut[128];
   TChain *run = new TChain("cal");
   sprintf(cmd,"%s*.root",cal_file.c_str());
   cout <<"cmd >"<<cmd<<"<<"<<endl;
   run->Add(cmd);
//   TFile *_f = new TFile(cal_file.c_str(),"READONLY");

   if        (plot_type == "chi2"){
      run->SetMarkerStyle(24);
      run->SetMarkerColor(1);
      sprintf(cmd,"chi2[%i]:cal_tmin",ich);
      run->Draw(cmd,"","prof");
   } else if (plot_type == "trend"){
      sprintf(cmd,"c0[%i]+c1[%i]*0.2e-6:cal_tmin",ich,ich);
      sprintf(cut,"chi2[%i]>%f",ich,chi2_cut[ich]);
      run->SetMarkerStyle(24);
      run->SetMarkerColor(1);
      run->Draw(cmd);
      run->SetMarkerStyle(24);
      run->SetMarkerColor(2);
      run->Draw(cmd,cut,"same");
   }
}
