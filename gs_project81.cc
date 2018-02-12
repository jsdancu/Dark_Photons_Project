//23/01/2018

//Reconstructing the Dark Photon's braching ratio vs it mass
//Reconstructing the Dark Photon's width vs its mass and kinetic mixing

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
//#include "Pythia8Plugins/CoupSM.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"
#include "TVectorD.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TRandom2.h"
#include "cmath"
#include <iostream>
#include <iomanip> 
#include <vector>
#include <fstream>
#include <stdio.h>

using namespace Pythia8;
//using namespace CoupSM;

double lepton_width(double m_A, double m_l, double epsilon2){

	double alpha_EM = 1.0/137.0;
	/*HepMC::Pythia8ToHepMC ToHepMC;
	Pythia pythia;
	CoupSM& coup = pythia.CoupSM;

	pythia.init(); 
	
	double alpha_EM = coup.alphaEM(pow(m_A, 2.0));*/

	double Gamma_l = epsilon2 * alpha_EM /3.0 * m_A * (1.0 + 2.0 * pow(m_l, 2.0) / pow(m_A, 2.0)) * sqrt(1.0 - 4.0 * pow(m_l, 2.0) / pow(m_A, 2.0));

	return Gamma_l;

}

double lepton_width1(double m_A, double m_l){

	/*HepMC::Pythia8ToHepMC ToHepMC;
	Pythia pythia;
	CoupSM& coup = pythia.CoupSM;

	pythia.init(); 
	double alpha_EM = coup.alphaEM(pow(m_A, 2.0));*/

	double alpha_EM = 1.0/137.0;

	double Gamma_l = alpha_EM /3.0 * m_A * (1.0 + 2.0 * pow(m_l, 2.0) / pow(m_A, 2.0)) * sqrt(1.0 - 4.0 * pow(m_l, 2.0) / pow(m_A, 2.0));

	return Gamma_l;

}

double hadrons_width(double R_mu, double Gamma_l){

	double Gamma_h = Gamma_l * R_mu;

	return Gamma_h;

}

int main() {

	//HepMC::Pythia8ToHepMC ToHepMC;

	//Open R-data file
	//TFile  *file_in = new TFile("rdata.txt", "READ");

	double x, y, u1, u2;
	std::vector<double> m_A, R_mu, v1, v2;

	for(int i=0;i<14;i++){

		m_A.push_back(0.22 + i * 0.01);
		R_mu.push_back(0.0);
		v1.push_back(0.0);
		v2.push_back(0.0);

	}


	std::ifstream file("rdata.txt");
    	if(file.is_open()){

		file >> x >> y >> u1 >> u2;
		m_A.push_back(x);
		R_mu.push_back(y);
		v1.push_back(u1);
		v2.push_back(u2);

		// Begin reading data
		while(file >> x >> y >> u1 >> u2){
	
			m_A.push_back(x);
			R_mu.push_back(y);
			v1.push_back(u1);
			v2.push_back(u2);

		}

    	}

	double mA[12], Rmu[12], w1[12], w2[12], wx1[12], wx2[12];

	for(int i=0;i<12;i++){

		mA[i] = m_A[i+14];
		Rmu[i] = R_mu[i+14];
		w1[i] = v1[i+14];
		w2[i] = abs(v2[i+14]);
		wx1[i] = 0.0;
		wx2[i] = 0.0;

	}

	TCanvas *c0=new TCanvas("c0","",600,600);

	TF1 *fit = new TF1("fit", "[0]*exp([1]*x)", 0.36, 0.46);

	TGraphAsymmErrors *gre = new TGraphAsymmErrors(12, mA, Rmu, wx2, wx1, w2, w1);
	gre->SetTitle("Fitting the 0.36-0.46 GeV region of R_{#mu}");
	gre->GetXaxis()->SetTitle("m_{A'} (GeV)");
	gre->GetYaxis()->SetTitle("R_{#mu}");
   	gre->Draw("ap");
   	gre->Fit("fit","EX0");

	c0->Print("project81_fit.pdf","pdf");	

	double p0 = fit->GetParameter(0);
	double p1 = fit->GetParameter(1);

	for(int i=6;i<14;i++){

		R_mu[i] = p0 * exp(p1 * m_A[i]);

	}

	/*Pythia pythia;
	ParticleData& pdt = pythia.particleData;

	pythia.init(); 

	double m_mu = pdt.m0(13);//muon mass
	double m_e = pdt.m0(11);//electron mass
	double m_tau = pdt.m0(15);//tau mass
	*/
	

	double m_mu = 0.10566;//muon mass
	double m_e = 5.110e-04;//electron mass
	double m_tau = 1.77682;//tau mass	

	double epsilon2 = 1e-11; 
	double eps2, G, G_mu, G_e, G_h;

	//double Gamma_mu, Gamma_mu1;
	double Gamma_inv = 0.0;

	int N = m_A.size();
	double Gamma[N], Gamma1[N], Br1[N], Gamma_mu[N], Gamma_mu1[N], Br2[N], Gamma_e[N], Gamma_e1[N], Br3[N], Gamma_h[N], Gamma_h1[N];

	TH2D *lifetime1 = new TH2D("lifetime1","Diagram showing the lifetime of Dark Photon based on its mass and #epsilon^{2}", 200, 0.0, 200, 1000, 1e-11, 1e-8);
	lifetime1 -> GetXaxis()-> SetTitle("m_{A'} (GeV)");
	lifetime1 -> GetYaxis()-> SetTitle("#epsilon^{2}");

	for(int i=0;i<N;i++){

		Gamma_mu[i] = lepton_width(m_A[i], m_mu, epsilon2);
		Gamma_e[i] = lepton_width(m_A[i], m_e, epsilon2);
		Gamma_h[i] = hadrons_width(R_mu[i], Gamma_mu[i]);

		Gamma_mu1[i] = lepton_width1(m_A[i], m_mu);
		Gamma_e1[i] = lepton_width1(m_A[i], m_e);
		Gamma_h1[i] = hadrons_width(R_mu[i], Gamma_mu1[i]);

		if(m_A[i] < 2.0*0.28){

			Gamma[i] = Gamma_e[i] + Gamma_mu[i] + Gamma_inv;
			Gamma1[i] = Gamma_e1[i] + Gamma_mu1[i] + Gamma_inv;

		}
		else if(m_A[i] < 2.0*m_tau){

			Gamma[i] = Gamma_e[i] + Gamma_mu[i] + Gamma_h[i] + Gamma_inv;
			Gamma1[i] = Gamma_e1[i] + Gamma_mu1[i] + Gamma_h1[i] + Gamma_inv;

		}
		else if(m_A[i] > 2.0*m_tau){

			Gamma[i] = Gamma_e[i] + Gamma_mu[i] + lepton_width(m_A[i], m_tau, epsilon2) + Gamma_h[i] + Gamma_inv;
			Gamma1[i] = Gamma_e1[i] + Gamma_mu1[i] + lepton_width1(m_A[i], m_tau) + Gamma_h1[i] + Gamma_inv;

		}


		Br1[i] = Gamma_mu[i]/Gamma[i];
		Br2[i] = Gamma_e[i]/Gamma[i];
		Br3[i] = Gamma_h[i]/Gamma[i];

		for(int j=0;j<1000;j++){

			eps2 = epsilon2 + j;

			G_mu = lepton_width(m_A[i], m_mu, eps2);
			G_e = lepton_width(m_A[i], m_e, eps2);
			G_h = hadrons_width(R_mu[i], G_mu);

			if(m_A[i] < 2.0*m_tau){

				G = G_e + G_mu + G_h + Gamma_inv;

			}
			else{
		
				G = G_e + G_mu + lepton_width(m_A[i], m_tau, eps2) + G_h + Gamma_inv;

			}

			lifetime1 -> SetBinContent(m_A[i], eps2, 1.0/G);

		}

	}

	TCanvas *c1=new TCanvas("c1","",600,600);

	lifetime1 -> Draw("COLZ");

	c1->Modified();
	c1->Update();

	c1->Print("project81_lifetime1.pdf","pdf");

	
  	TH1D*lifetime2=new TH1D("h1","Diagram showing the lifetime of Dark Photon based on its mass and #epsilon^{2}",200, 0.3, 200);
	lifetime2 = lifetime1->ProjectionY("#epsilon^{2} projection");
	lifetime2->SetTitle("Diagram showing the lifetime of Dark Photon based on its mass and #epsilon^{2}");
  	lifetime2->GetXaxis()->SetTitle("#epsilon^{2}");
  	lifetime2->GetYaxis()->SetTitle("arbitrary units");
	lifetime2->GetYaxis()->SetTitleOffset(0.8);

	c1->Modified();
	c1->Update();

	c1->Print("project81_lifetime2.pdf","pdf");	


	gPad->SetLogy();

	TGraph *check = new TGraph(N, &m_A[0], &R_mu[0]);
	check -> SetMarkerColor(kBlue);
	check -> SetLineColor(kBlue);
	check -> SetTitle("R_{#mu} vs m_{A'}");
    	check -> GetXaxis()-> SetTitle("m_{A'} (GeV)");
	check -> GetYaxis()-> SetTitle("R_{#mu}");
	check -> GetYaxis()->SetTitleOffset(1.1);

	check->Draw("A*");
	check -> Fit("gaus");

	c1->Modified();
	c1->Update();

	c1->Print("project81_RmuvsmA.pdf","pdf");

	gPad->SetLogx();

	TGraph *width1 = new TGraph(N, &m_A[0], Gamma);
	width1 -> SetMarkerColor(kBlue);
	width1 -> SetLineColor(kBlue);

	TGraph *width_mu1 = new TGraph(N, &m_A[0], Gamma_mu);
	width_mu1 -> SetMarkerColor(kRed);
	width_mu1 -> SetLineColor(kRed);

	TGraph *width_e1 = new TGraph(N, &m_A[0], Gamma_e);
	width_e1 -> SetMarkerColor(kGreen);
	width_e1 -> SetLineColor(kGreen);

	TGraph *width_h1 = new TGraph(N, &m_A[0], Gamma_h);
	width_h1 -> SetMarkerColor(kOrange);
	width_h1 -> SetLineColor(kOrange);

	TMultiGraph *mgr = new TMultiGraph();
	mgr->SetTitle("Dark Photon width distribution based on its mass (#epsilon^{2} = 10^{-11})");
    	mgr->Add(width1,"pl");
    	mgr->Add(width_mu1,"pl");
    	mgr->Add(width_e1,"pl");
    	mgr->Add(width_h1,"pl");
	mgr->Draw("ap");
	mgr->GetXaxis()->SetTitle("m_{A'} (GeV)");
	mgr->GetXaxis()->SetTitleOffset(1);
  	mgr->GetYaxis()->SetTitle("Width (GeV^{-1})");
	mgr->GetYaxis()->SetTitleOffset(1.5);   

	TLegend *legend1 = new TLegend(0.1,0.7,0.4,0.9);
	TLegendEntry *l11 = legend1->AddEntry("width1","#Gamma_{total}","l");	
	l11->SetLineColor(kBlue);
	TLegendEntry *l12 = legend1->AddEntry("width_mu1","#Gamma_{A' #rightarrow #mu^{+} #mu^{-}}","l");
	l12->SetLineColor(kRed);
	TLegendEntry *l13 = legend1->AddEntry("width_e1","#Gamma_{A' #rightarrow e^{+} e^{-}}","l");
	l13->SetLineColor(kGreen);
	TLegendEntry *l14 = legend1->AddEntry("width_h1","#Gamma_{A' #rightarrow hadrons}","l");
	l14->SetLineColor(kOrange);
	legend1->Draw();

	mgr->GetYaxis()->SetRangeUser(1e-16, 1e-8);

	c1->Modified();
	c1->Update();
	c1->Print("project81_width1.pdf","pdf");


	TGraph *width2 = new TGraph(N, &m_A[0], Gamma1);
	width2 -> SetMarkerColor(kBlue);
	width2 -> SetLineColor(kBlue);

	TGraph *width_mu2 = new TGraph(N, &m_A[0], Gamma_mu1);
	width_mu2 -> SetMarkerColor(kRed);
	width_mu2 -> SetLineColor(kRed);

	TGraph *width_e2 = new TGraph(N, &m_A[0], Gamma_e1);
	width_e2 -> SetMarkerColor(kGreen);
	width_e2 -> SetLineColor(kGreen);

	TGraph *width_h2 = new TGraph(N, &m_A[0], Gamma_h1);
	width_h2 -> SetMarkerColor(kOrange);
	width_h2 -> SetLineColor(kOrange);

	TMultiGraph *mgr2 = new TMultiGraph();
	mgr2->SetTitle("Dark Photon width/#epsilon^{2} distribution based on its mass");
    	mgr2->Add(width2,"pl");
    	mgr2->Add(width_mu2,"pl");
    	mgr2->Add(width_e2,"pl");
    	mgr2->Add(width_h2,"pl");
	mgr2->Draw("ap");
	mgr2->GetXaxis()->SetTitle("m_{A'} (GeV)");
	mgr2->GetXaxis()->SetTitleOffset(1);
  	mgr2->GetYaxis()->SetTitle("Width/#epsilon^{2} (GeV^{-1})");
	mgr2->GetYaxis()->SetTitleOffset(1.5);   

	TLegend *legend2 = new TLegend(0.1,0.7,0.4,0.9);
	TLegendEntry *l21 = legend2->AddEntry("width2","#Gamma_{total}","l");	
	l21->SetLineColor(kBlue);
	TLegendEntry *l22 = legend2->AddEntry("width_mu2","#Gamma_{A' #rightarrow #mu^{+} #mu^{-}}","l");
	l22->SetLineColor(kRed);
	TLegendEntry *l23 = legend2->AddEntry("width_e2","#Gamma_{A' #rightarrow e^{+} e^{-}}","l");
	l23->SetLineColor(kGreen);
	TLegendEntry *l24 = legend2->AddEntry("width_h2","#Gamma_{A' #rightarrow hadrons}","l");
	l24->SetLineColor(kOrange);
	legend2->Draw();

	mgr2->GetYaxis()->SetRangeUser(1e-5, 1e3);

	c1->Modified();
	c1->Update();
	c1->Print("project81_width2.pdf","pdf");


	TGraph *branching1 = new TGraph(N, &m_A[0], Br1);
	branching1 -> SetMarkerColor(kBlue);
	branching1 -> SetLineColor(kBlue);

	TGraph *branching2 = new TGraph(N, &m_A[0], Br2);
	branching2 -> SetMarkerColor(kRed);
	branching2 -> SetLineColor(kRed);

	TGraph *branching3 = new TGraph(N, &m_A[0], Br3);
	branching3 -> SetMarkerColor(kGreen);
	branching3 -> SetLineColor(kGreen);

	TMultiGraph *mgr3 = new TMultiGraph();
	mgr3->SetTitle("Dark Photon branching ratios vs its mass");
    	mgr3->Add(branching1,"pl");
    	mgr3->Add(branching2,"pl");
    	mgr3->Add(branching3,"pl");
	mgr3->Draw("ap");
	mgr3->GetXaxis()->SetTitle("m_{A'} (GeV)");
	mgr3->GetXaxis()->SetTitleOffset(1);
  	mgr3->GetYaxis()->SetTitle("Br(A')");
	mgr3->GetYaxis()->SetTitleOffset(1);   

	TLegend *legend3 = new TLegend(0.1,0.3,0.4,0.5);
	TLegendEntry *l31 = legend3->AddEntry("branching1","Br(A' #rightarrow #mu^{-} #mu^{+})","l");	
	l31->SetLineColor(kBlue);
	TLegendEntry *l32 = legend3->AddEntry("branching2","Br(A' #rightarrow e^{-} e^{+})","l");
	l32->SetLineColor(kRed);
	TLegendEntry *l33 = legend3->AddEntry("branching3","Br(A' #rightarrow hadrons)","l");
	l33->SetLineColor(kGreen);
	legend3->Draw();

	mgr3->GetXaxis()->SetRangeUser(0.3, 1.6);

	c1->Modified();
	c1->Update();
	c1->Print("project81_Abranching1.pdf","pdf");

	//gPad->SetLogx();

	mgr3->GetXaxis()->SetRangeUser(0.3, 200);
	//branching->Draw("AP");

	c1->Modified();
	c1->Update();
	c1->Print("project81_Abranching.pdf","pdf");

	//delete file_in;
	file.close();

	// Done.
	return 0;

}
