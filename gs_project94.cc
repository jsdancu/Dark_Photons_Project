//27/02/2018

//Script evaluating the muon misID rate in the LHCb detector

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

double average(TTree *T3){

	double pl_avg = 0.0;

	double id_var, pz_var;

   	T3->SetBranchAddress("id",&id_var);
	T3->SetBranchAddress("pz",&pz_var);

	//Loop through the entries of the tree
	Long64_t nentries = T3->GetEntries();
	for (Long64_t i=0;i<nentries;i++){

	      	T3->GetEntry(i);

		if((id_var == 211) || (id_var == -211)){

			pl_avg = pl_avg + abs(pz_var);

		}

	}

	pl_avg = pl_avg/nentries;

	return pl_avg;

}

double error(double dy, double misid, double misid_avg, double err_misid_avg){

	double err = sqrt((pow(dy / 1e3, 2) + pow(misid * err_misid_avg, 2)) / misid_avg);

	return err;

}

int main() {

	// Open the input TFile 
	TFile  *file_in = new TFile("gs_project34.root", "READ");

	//Get misID tree from file
	TTree *T3 = (TTree*)file_in->Get("T3");

	double pl_avg = average(T3);

	std::cout<<"Average longitudinal momentum in MC data: "<<pl_avg<<std::endl;

	std::vector<double> pl, misid, dxm, dxp, dym, dyp;
	double x, x1, x2, x3, x4, x5, x6;

	//std::ifstream file("misidexport.csv");
	std::ifstream file("misid.dat");
    	if(file.is_open()){

		file >> x1 >> x2 >> x3 >> x4 >> x5 >> x6;
		pl.push_back(x1/1e3);
		misid.push_back(x2);
		dxm.push_back(x3/1e3);
		dxp.push_back(x4/1e3);
		dym.push_back(x5);
		dyp.push_back(x6);

		// Begin reading data
		while(file >> x1 >> x2 >> x3 >> x4 >> x5 >> x6){
	
			pl.push_back(x1/1e3);
			misid.push_back(x2);
			dxm.push_back(x3/1e3);
			dxp.push_back(x4/1e3);
			dym.push_back(x5);
			dyp.push_back(x6);

		}

    	}

	int n = pl.size();

	TCanvas *c0=new TCanvas("c0","",600,600);

	gPad->SetLogy();
	gPad->SetLogx();

	//TF1 *fit1 = new TF1("fit1", "1.0 - exp(-[0]/x) + [1] * x + [2]", 10, 300);
	TF1 *fit1 = new TF1("fit1", "1.0 - exp(-[0]/x) + [1] * x + [2]", 10, 250);

	TGraphAsymmErrors *gre1 = new TGraphAsymmErrors(n, &pl[0], &misid[0], &dxm[0], &dxp[0], &dym[0], &dyp[0]);
	gre1->SetTitle("Graph showing the #mu misID rate in the LHCb detector");
	gre1->GetXaxis()->SetTitle("p_{L} (MeV)");
	gre1->GetYaxis()->SetTitle("Track misID probability");
	gre1->GetYaxis()->SetTitleOffset(1.4);
   	gre1->Draw("ap");
   	gre1->Fit("fit1","EX0");

	c0->Print("project94_fit1.pdf","pdf");

	double alpha1 = fit1->GetParameter(0);
	double beta1 = fit1->GetParameter(1);	
	double gamma1 = fit1->GetParameter(2);
	double err_alpha1 = fit1->GetParError(0);
	double err_beta1 = fit1->GetParError(1);	
	double err_gamma1 = fit1->GetParError(2);

	double misid_avg = 1.0 - exp(-alpha1/pl_avg) + beta1 * pl_avg + gamma1;
	double err_misid_avg = sqrt(pow(err_alpha1 * exp(-alpha1/pl_avg) / pl_avg, 2) + pow(pl_avg * err_beta1, 2) + pow( err_gamma1, 2));

	for(int i=0;i<n;i++){

		dym[i] = error(dym[i], misid[i], misid_avg, err_misid_avg);
		dyp[i] = error(dyp[i], misid[i], misid_avg, err_misid_avg);
		misid[i] = misid[i] /(1e3 * misid_avg);

	}

	//TF1 *fit = new TF1("fit", "1.0 - exp(-[0]/x) + [1] * x + [2]", 10, 300);
	TF1 *fit = new TF1("fit", "1.0 - exp(-[0]/x) + [1] * x + [2]", 10, 250);

	TGraphAsymmErrors *gre = new TGraphAsymmErrors(n, &pl[0], &misid[0], &dxm[0], &dxp[0], &dym[0], &dyp[0]);
	gre->SetTitle("Graph showing the #mu misID rate in the LHCb detector after scaling");
	gre->GetXaxis()->SetTitle("p_{L} (MeV)");
	gre->GetYaxis()->SetTitle("Track misID probability");
	gre->GetYaxis()->SetTitleOffset(1.4);
   	gre->Draw("ap");
   	gre->Fit("fit","EX0");

	c0->Print("project94_fit.pdf","pdf");

	double alpha = fit->GetParameter(0);
	double beta = fit->GetParameter(1);	
	double gamma = fit->GetParameter(2);

	std::cout<<"alpha: "<<alpha<<std::endl;
	std::cout<<"beta: "<<beta<<std::endl;
	std::cout<<"gamma: "<<gamma<<std::endl;

	FILE  *file_out = fopen("results_out94.txt", "w");	
	fprintf(file_out, "%f \n", alpha);
	fprintf(file_out, "%f \n", beta);
	fprintf(file_out, "%f \n", gamma);

	file.close();
	fclose(file_out);

	// Done.
	return 0;

}
