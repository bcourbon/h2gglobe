#include "functions.h"
#include "THStack.h"
#include "setTDRStyle.C"

//options:
// - category:  event category: "0" is all cats combined, "1" to "6" correspond to cats "0" to "5" (i.e. inclusive and dijet)
// - var:  which variable to plot: for possible options see runMassPlot.sh
// - blind:  omit data on mass plots for signal region
// - breakdown:  set to true for full breakdown of background MC contributions
// - omitZpeak:  set to true to restrict range on mass plots to 100-180 (as used for public plot)
// - cleanHighWeightEvents:  remove the high weight QCD40_PF events and scale up the GJet40_PF events by ((N_GJet40_PF+N_QCD40_PF)/N_GJet40_PF) where N_X is the histogram integral for sample X (as used for public plot)

//public plot https://twiki.cern.ch/twiki/pub/CMSPublic/Hig13001TWiki/mass.png uses default settings except omitZpeak=true, cleanHighWeightEvents=true

//void massPlot(TString category="0", TString var="all_mass", TString analysis="massfit", bool blind=false, bool breakdown=false, bool omitZpeak=false, bool cleanHighWeightEvents=false) {
int main(int argc, char *argv[])
{ 

        for(int iarg = 0 ; iarg < argc; iarg++)
        {
                cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
        }

        if( argc == 1 ) 
        {
                cerr << "arguments should be passed : fileHistos, filesig, blind, doBands, baseline, useBinnedData" <<endl; 
                return 1;
        }    

	string output="0";
	string input="blah.root";
	TString category="0";
	string var="all_mass";
	TString analysis="massfit";
	bool blind=false;
	bool breakdown=false;
	bool omitZpeak=false;
	bool cleanHighWeightEvents=false;

	if( argc > 1 ) output = argv[1];
	if( argc > 2 ) input = argv[2];
	if( argc > 3 ) category = argv[3];
        if( argc > 4 ) var = argv[4];
	if( argc > 5 ) analysis = argv[5];
        if( argc > 6 ) 
        {
                std::stringstream ss ( argv[6] );
                ss >> blind;
        }	
	if( argc > 7 ) 
        {
                std::stringstream ss ( argv[7] );
                ss >> breakdown;
        }
	if( argc > 8 ) 
        {
                std::stringstream ss ( argv[8] );
                ss >> omitZpeak;
        }
	if( argc > 9 ) 
        {
                std::stringstream ss ( argv[9] );
                ss >> cleanHighWeightEvents;
        }	

	system(Form("mkdir -p %s",output.c_str() ) );

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetTitleFont(42,"XYZ");

  gStyle->SetMarkerSize(0.8);
  gStyle->SetTitleYOffset(1.5);

  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);

  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.05,"XYZ");

  gStyle->SetTitleSize(0.06,"XYZ");
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.24);

  //TFile *f = TFile::Open("root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/histograms_moriond_david/histograms_CMS-HGG_"+analysis+".root");
  TFile *f = TFile::Open(Form("%s",input.c_str())); //FIXME	
  //f->ls();
	cout << endl << "coucou 2" << endl;

  TH1F *hist_mass_sig;
  TH1F *hist_mass_data;
  TH1F *hist_mass_born;
  //TH1F *hist_mass_box;
  TH1F *hist_mass_gjet_pf;
  TH1F *hist_mass_gjet20_pf;
  TH1F *hist_mass_qcd_pf;
  TH1F *hist_mass_qcd30_pf;
  TH1F *hist_mass_qcd_ff;
  TH1F *hist_mass_qcd30_ff;
  TH1F *hist_mass_dy;
  TH1F *hist_mass_bkg;
  THStack* hist_mass_bkg_stack;
  TH1F * hist_mass_qcd40_pf;
  TH1F * hist_mass_pp;
  TH1F * hist_mass_pf;
  TH1F * hist_mass_ff;
  TH1F * hist_mass_bkg_stack_err;
/*
	cout << endl << "coucou 3" << endl;
  hist_mass_data = (TH1F*)(f->Get(var+"_cat"+category+"_Data"))->Clone();
	cout << endl << "coucou 4" << endl;
  hist_mass_sig = (TH1F*)(f->Get(var+"_cat"+category+"_ggh_m125_8TeV"))->Clone();
	cout << endl << "coucou 5" << endl;
  hist_mass_sig->Add((TH1F*)(f->Get(var+"_cat"+category+"_vbf_m125_8TeV")));
	cout << endl << "coucou 6" << endl;
  hist_mass_sig->Add((TH1F*)(f->Get(var+"_cat"+category+"_wzh_m125_8TeV")));
	cout << endl << "coucou 7" << endl;
  hist_mass_sig->Add((TH1F*)(f->Get(var+"_cat"+category+"_tth_m125_8TeV")));
	cout << endl << "coucou 8" << endl;
  hist_mass_born = (TH1F*)(f->Get(var+"_cat"+category+"_diphojet_sherpa_8TeV"))->Clone(); //BOX + BORN
	cout << endl << "coucou 9" << endl;
  //hist_mass_box = (TH1F*)(f->Get(var+"_cat"+category+"_DiPhotonBox25"))->Clone();
	cout << endl << "coucou 10" << endl;
  hist_mass_gjet_pf = (TH1F*)(f->Get(var+"_cat"+category+"_gjet_40_8TeV_pf"))->Clone();
	cout << endl << "coucou 11" << endl;
  hist_mass_gjet20_pf = (TH1F*)(f->Get(var+"_cat"+category+"_gjet_20_8TeV_pf"))->Clone();
	cout << endl << "coucou 12" << endl;  

  if (cleanHighWeightEvents) {
    hist_mass_qcd40_pf = (TH1F*)(f->Get(var+"_cat"+category+"_qcd_40_8TeV_pf"))->Clone();
    float n_qcd40_pf = hist_mass_qcd40_pf->Integral();
    float n_gjet40_pf = hist_mass_gjet_pf->Integral();
    cout << n_qcd40_pf << " " << n_gjet40_pf << endl;
    hist_mass_gjet_pf->Scale((n_gjet40_pf+n_qcd40_pf)/n_gjet40_pf);
  }
  hist_mass_gjet_pf->Add(hist_mass_gjet20_pf);

  if (cleanHighWeightEvents) {
    hist_mass_qcd_pf = (TH1F*)(f->Get(var+"_cat"+category+"_qcd_30_8TeV_pf"))->Clone();
  } else {
    hist_mass_qcd_pf = (TH1F*)(f->Get(var+"_cat"+category+"_qcd_40_8TeV_pf"))->Clone();
    hist_mass_qcd30_pf = (TH1F*)(f->Get(var+"_cat"+category+"_qcd_30_8TeV_pf"))->Clone();
    hist_mass_qcd_pf->Add(hist_mass_qcd30_pf);
  }

  hist_mass_qcd_ff = (TH1F*)(f->Get(var+"_cat"+category+"_qcd_40_8TeV_ff"))->Clone();
  hist_mass_qcd30_ff = (TH1F*)(f->Get(var+"_cat"+category+"_qcd_30_8TeV_ff"))->Clone();
  hist_mass_qcd_ff->Add(hist_mass_qcd30_ff);

	cout << endl << "coucou 13" << endl;
*/
  hist_mass_dy = (TH1F*)(f->Get(var+"_cat"+category+"_DYJetsToLL"))->Clone(); //FIXME

  if (var=="all_mass" || var=="pt") {
    hist_mass_dy->Rebin(2);
  }

  float lumiSF=1.;
  hist_mass_dy->Scale(lumiSF);

  //fix DY normalization
  //hist_mass_dy->Scale(2.089*3532.8/(2950.*1.15));
  //hist_mass_dy->Scale((2950.*1.15)/(3532.8));

//   cout << hist_mass_born->Integral() << ", ";
//   cout << hist_mass_box->Integral() << ", ";
//   cout << hist_mass_gjet_pf->Integral() << ", ";
//   cout << hist_mass_qcd_pf->Integral() << ", ";
//   cout << hist_mass_qcd_ff->Integral() << ", ";
//   cout << hist_mass_dy->Integral() << endl;

  hist_mass_dy->SetFillColor(kMagenta-10);
  hist_mass_dy->SetLineWidth(2.5);

  hist_mass_bkg_stack = new THStack("hist_mass_bkg_stack","Background");
  hist_mass_bkg_stack->Add(hist_mass_dy);


  TCanvas *c_mgg = new TCanvas("c_mgg","Mgg",1000,700);
  c_mgg->SetFillColor(0);

  hist_mass_bkg_stack->Draw("hist");
  if (var=="all_mass") {
    if (omitZpeak) {
      hist_mass_bkg_stack->GetXaxis()->SetRangeUser(100.,179.);
    } else {
      //hist_mass_bkg_stack->GetXaxis()->SetRangeUser(60.,190.); //FIXME
      hist_mass_bkg_stack->GetXaxis()->SetLimits(80.,150.);
    }
  }
  //hist_mass_bkg_stack->GetXaxis()->SetRangeUser(80.,120.);
  //hist_mass_bkg_stack_err->Draw("same,e2");
  //hist_mass_bkg_stack->GetYaxis()->SetTitleSize(0.05);
  hist_mass_bkg_stack->GetYaxis()->SetTitle("Events / 2 GeV");
  //hist_mass_bkg_stack->GetXaxis()->SetTitleOffset(0.9);
  //hist_mass_bkg_stack->GetXaxis()->SetTitleSize(0.05);
  if (var=="all_mass") {
    hist_mass_bkg_stack->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
  } else if (var=="pho1_pt") {
    hist_mass_bkg_stack->GetXaxis()->SetTitle("lead photon p_{T} (GeV)");
  } else if (var=="pho2_pt") {
    hist_mass_bkg_stack->GetXaxis()->SetTitle("sublead photon p_{T} (GeV)");
  } else if (var=="pho1_r9") {
    hist_mass_bkg_stack->GetXaxis()->SetTitle("lead photon R_{9} (GeV)");
  } else if (var=="pho2_r9") {
    hist_mass_bkg_stack->GetXaxis()->SetTitle("sublead photon R_{9} (GeV)");
  } else if (var=="pho_r9") {
    hist_mass_bkg_stack->GetXaxis()->SetTitle("photon R_{9} (GeV)");
  } else if (var=="pho1_eta") {
    hist_mass_bkg_stack->GetXaxis()->SetTitle("lead photon #eta (GeV)");
  } else if (var=="pho2_eta") {
    hist_mass_bkg_stack->GetXaxis()->SetTitle("sublead photon #eta (GeV)");
  } else if (var=="eta") {
    hist_mass_bkg_stack->GetXaxis()->SetTitle("di-photon #eta (GeV)");
  } else if (var=="pt") {
    hist_mass_bkg_stack->GetXaxis()->SetTitle("di-photon p_{T} (GeV)");
    hist_mass_bkg_stack->GetXaxis()->SetRangeUser(0.,120.);
  } else {
    cout << "unknown variable: " << var << endl;
  }
  //hist_mass_bkg_stack->GetXaxis()->SetTitle("lead photon #pt");
  //hist_mass_bkg_stack->GetXaxis()->SetTitle("lead photon p_{T}");
  //hist_mass_bkg_stack->GetXaxis()->SetTitle("lead photon R_{9}");
  //hist_mass_sig->Draw("hist,same");
  //if (var!="all_mass") hist_mass_data->Draw("same,e");
  if (var=="all_mass"&& blind) {
    for (int ibin=16; ibin<36; ibin++) {
      hist_mass_data->SetBinContent(ibin,-999.);
    }
  }
  //hist_mass_data->Draw("same,e");

  TLegend *leg_mass;
  if (var=="pho1_r9" || var=="pho2_r9" || var=="pho_r9") {
    leg_mass = new TLegend(.3,.48,.72,.78);
  }else {
    leg_mass = new TLegend(.7,.8,.92,.9);
  }
  leg_mass->SetBorderSize(0);
  leg_mass->SetFillColor(10);
  leg_mass->SetTextSize(.03);
  if (breakdown) {
    //leg_mass->AddEntry(hist_mass_data,"Data","P");
    //leg_mass->AddEntry(hist_mass_sig,"H#rightarrow#gamma#gamma (125 GeV) #times5","F");
    leg_mass->AddEntry(hist_mass_dy,"DYee+Z","F");
    //leg_mass->AddEntry(hist_mass_qcd_ff,"QCD fake-fake","F");
    //leg_mass->AddEntry(hist_mass_qcd_pf,"QCD prompt-fake","F");
    //leg_mass->AddEntry(hist_mass_gjet_pf,"GJet prompt-fake","F");
    //leg_mass->AddEntry(hist_mass_born,"Sherpa DiPhotonJets","F");
    //leg_mass->AddEntry(hist_mass_box,"Pythia Box","F");
  } else {
    //leg_mass->AddEntry(hist_mass_data,"Data","P");
    //leg_mass->AddEntry(hist_mass_pp,"2 prompt #gamma","F");
    //leg_mass->AddEntry(hist_mass_pf,"1 prompt #gamma 1 fake #gamma","F");
    //leg_mass->AddEntry(hist_mass_ff,"2 fake #gamma","F");
    leg_mass->AddEntry(hist_mass_dy,"Drell-Yan","F");
    //leg_mass->AddEntry(hist_mass_sig,"H#rightarrow#gamma#gamma (125 GeV) #times5","F");
  }
  if (var!="eta") leg_mass->Draw();

  TLatex * txt = new TLatex();
  txt->SetNDC();
  txt->SetTextSize(0.05);
  txt->SetTextAlign(12);
  txt->DrawLatex(0.41,0.85,"#scale[0.8]{#splitline{CMS Private}{#sqrt{s} = 8 TeV L = 19.7 fb^{-1}}}");
  //txt->DrawLatex(0.26,0.82,"CMS preliminary");
  //txt->DrawLatex(0.26,0.75,"#sqrt{s} = 7 TeV L = 1.66 fb^{-1}");

  //c_mgg->SaveAs("mass.gif");
  //c_mgg->SaveAs("mass.png");

  //hist_mass_bkg_stack->GetYaxis()->SetRangeUser(0,30000);
	//hist_mass_bkg_stack->SetMaximum(13000); //FIXME
	hist_mass_bkg_stack->SetMaximum(hist_mass_bkg_stack->GetMaximum() * 1.3);

  string outcat="";
  string outvar=var;
  if (category=="1") outcat = "_cat0";
  if (category=="2") outcat = "_cat1";
  if (category=="3") outcat = "_cat2";
  if (category=="4") outcat = "_cat3";
  if (category=="5") outcat = "_cat4";
  if (category=="6") outcat = "_cat5";
  if (category=="7") outcat = "_cat6";
  if (category=="8") outcat = "_cat7";
  if (category=="9") outcat = "_cat8";
  if (category=="10") outcat = "_cat9";
  if (category=="11") outcat = "_cat10";
  if (category=="12") outcat = "_cat11";
  if (category=="13") outcat = "_cat12";
  if (category=="14") outcat = "_cat13"; 	
  if (category=="15") outcat = "_cat14";
  if (category=="16") outcat = "_cat15";
  if (var=="all_mass") outvar = "mass";
  //c_mgg->SaveAs("/afs/cern.ch/user/f/futyand/www/hgg/moriond_preapproval/dataMC/"+analysis+"/"+outvar+"/"+outvar+outcat+".png");
  c_mgg->Print(Form("%s/%s%s.C",output.c_str(),outvar.c_str(),outcat.c_str()));
  c_mgg->Print(Form("%s/%s%s.png",output.c_str(),outvar.c_str(),outcat.c_str()));
  c_mgg->Print(Form("%s/%s%s.pdf",output.c_str(),outvar.c_str(),outcat.c_str()));
}
