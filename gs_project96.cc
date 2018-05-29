//05/03/2018
//Program reading in the three trees generated in gs_project24.cc and gs_project26.cc and exploring the efficiency/feasibility of cuts on photon pt and p
//Program inspired from gs_project85.cc

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"
#include "cmath"
#include <iostream>
#include <vector>
#include <fstream>
#include <stdio.h>



using namespace std;

//#include "lhcbStyle.C"

struct vect{
	
	std::vector<double> index, id, energy, mass, px, py, pz, mother1, mother2, motherid1, motherid2;
	
};

struct vect1{
	
	std::vector<int> S, SMB, B, misID;
	
};

struct vect2{
	
	int S[80][100]; int SMB[80][100]; int B[80][100]; int misID[80][100];
	
};

double invmass(vect *v, Long64_t i, Long64_t j, Long64_t k){

	double inv_mass=sqrt(pow(v->energy[i]+v->energy[j]+v->energy[k],2)-pow(v->px[i]+v->px[j]+v->px[k],2)-pow(v->py[i]+v->py[j]+v->py[k],2)-pow(v->pz[i]+v->pz[j]+v->pz[k],2));

	return inv_mass;

}

double invmass2(vect *v, Long64_t i, Long64_t j){

	double inv_mass=sqrt(pow(v->energy[i]+v->energy[j],2)-pow(v->px[i]+v->px[j],2)-pow(v->py[i]+v->py[j],2)-pow(v->pz[i]+v->pz[j],2));

	return inv_mass;

}

double p(vect *v, Long64_t i){

	double p1=sqrt(pow(v->px[i],2)+pow(v->py[i],2)+pow(v->pz[i],2));

	return p1;

}

double pt(vect *v, Long64_t i){

	double pt1=sqrt(pow(v->px[i],2)+pow(v->py[i],2));

	return pt1;

}

double eta(vect *v, Long64_t i){

	double eta=0.5 * log((p(v, i) + abs(v->pz[i])) / (p(v, i) - abs(v->pz[i])));

	return eta;

}

double misID_rate(double pl, double alpha, double beta, double gamma){

	double misid = 1.0 - exp(-alpha/pl) + beta * pl + gamma;

	return misid;

}

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

double Gamma_A(double epsilon2, double massA){

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

	TF1 *fit = new TF1("fit", "[0]*exp([1]*x)", 0.36, 0.46);

	TGraphAsymmErrors *gre = new TGraphAsymmErrors(12, mA, Rmu, wx2, wx1, w2, w1);
   	gre->Fit("fit","EX0");

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

	double Gamma_inv = 0.0;

	int N = m_A.size();
	double Gamma[N], Gamma_mu[N], Gamma_e[N], Gamma_h[N];

	double G;

	for(int i=0;i<N;i++){

		Gamma_mu[i] = lepton_width(m_A[i], m_mu, epsilon2);
		Gamma_e[i] = lepton_width(m_A[i], m_e, epsilon2);
		Gamma_h[i] = hadrons_width(R_mu[i], Gamma_mu[i]);

		if(m_A[i] < 2.0*0.28){

			Gamma[i] = Gamma_e[i] + Gamma_mu[i] + Gamma_inv;

		}
		else if(m_A[i] < 2.0*m_tau){

			Gamma[i] = Gamma_e[i] + Gamma_mu[i] + Gamma_h[i] + Gamma_inv;;

		}
		else if(m_A[i] > 2.0*m_tau){

			Gamma[i] = Gamma_e[i] + Gamma_mu[i] + lepton_width(m_A[i], m_tau, epsilon2) + Gamma_h[i] + Gamma_inv;

		}

		if(massA == m_A[i]){G = Gamma[i];}

	}

	return G;

}

double scaling(double epsilon2, double massA, double width){

	double resolution;

	if(massA < 1.0){resolution = 4e-3;}
	else{resolution = 4e-3 * massA;}

	double F = pow(epsilon2, 2) * TMath::Pi()/8.0 * pow(massA, 2)/(width * resolution);

	return F;

}

double metric(double S, double SMB, double B, double misID){

	//double sigma = S/sqrt(B + misID);
	double sigma = S/sqrt(SMB + B + misID);

	return sigma;

}

//void analyze_event1(FILE *file_out1, vect *v, vect1 *vmu, vect2 *vgamma, int column, int row, double mu_ptcut, double mu_pcut, double eta_cut1, double eta_cut2, double gamma_ptcut0, double gamma_pcut0, double gamma_ptcut_step, double gamma_pcut_step, double massA, double dmassA){
void analyze_event1(FILE *file_out1, vect *v, vect1 *vmu, double mu_ptcut, double mu_pcut, double eta_cut1, double eta_cut2, double massA, double dmassA){

	double inv_mass, mu_invmass, pt1, pt2, pt3, p1, p2, p3, eta1, eta2, eta3;

	Long64_t nentries = v->index.size();

	vect1* N = new vect1;
	for(int i=0;i<3;i++){
		N->S.push_back(0);
		N->SMB.push_back(0);
		N->B.push_back(0);
		N->misID.push_back(0);
	}

	double m_eta = 0.54785;
	double dm_eta = 0.02;

	for(Long64_t i=0;i<nentries;i++){ 

		//if entry is muon
		if(v->id[i] == 13){

			pt1 = pt(v, i);
			p1 = p(v, i);
			eta1 = eta(v, i);

				//Loop through the event and find an anti-muon
				for(Long64_t j=0;j<nentries;j++){

					//if entry is anti-muon
					if(v->id[j] == -13){

						pt2 = pt(v, j);
						p2 = p(v, j);
						eta2 = eta(v, j);
								

									//Reconstructing invariant mass of all muon-antimuon-photon combinations and di-muons
									mu_invmass = invmass2(v, i, j);

									if((v->mother1[i] == v->mother1[j]) && (v->mother2[i] == v->mother2[j]) && (v->motherid1[i] == 221)){									

										++N->SMB[0];

										if((pt1>mu_ptcut) && (p1>mu_pcut) && (eta1>eta_cut1) && (eta1<eta_cut2) && (pt2>mu_ptcut) &&  (p2>mu_pcut) && (eta2>eta_cut1) && (eta2<eta_cut2)){ 

											++N->SMB[1];

											if((mu_invmass < massA+dmassA) && (mu_invmass> massA-dmassA)){

												++N->SMB[2];

											}

										}

									}
									else{

										++N->B[0];

										if((pt1>mu_ptcut) && (p1>mu_pcut) && (eta1>eta_cut1) && (eta1<eta_cut2) && (pt2>mu_ptcut) &&  (p2>mu_pcut) && (eta2>eta_cut1) && (eta2<eta_cut2)){

											++N->B[1];

											if((mu_invmass < massA+dmassA) && (mu_invmass> massA-dmassA)){

												++N->B[2];

											}

										}

									}
					

					}

				}

		}

	}

	for(int i=0;i<3;i++){
		vmu->B[i] = vmu->B[i] + N->B[i];
		vmu->SMB[i] = vmu->SMB[i] + N->SMB[i];
	}

	delete N;

}

//void analyze_event2(FILE *file_out1, vect *v, vect1 *vmu, vect2 *vgamma, int column, int row, double mu_ptcut, double mu_pcut, double eta_cut1, double eta_cut2, double gamma_ptcut0, double gamma_pcut0, double gamma_ptcut_step, double gamma_pcut_step, double massA, double dmassA){
void analyze_event2(FILE *file_out1, vect *v, vect1 *vmu, double mu_ptcut, double mu_pcut, double eta_cut1, double eta_cut2, double massA, double dmassA){

	double inv_mass, mu_invmass, pt1, pt2, pt3, p1, p2, p3, eta1, eta2, eta3;

	Long64_t nentries = v->index.size();

	vect1* N = new vect1;
	for(int i=0;i<3;i++){
		N->S.push_back(0);
		N->B.push_back(0);
		N->misID.push_back(0);
	}

	double m_eta = 0.54785;
	double dm_eta = 0.02;

	for(Long64_t i=0;i<nentries;i++){ 

		//if entry is muon
		if(v->id[i] == 13){

			pt1 = pt(v, i);
			p1 = p(v, i);
			eta1 = eta(v, i);

			//if((pt1>mu_ptcut) && (p1>mu_pcut) && (eta1>eta_cut1) && (eta1<eta_cut2)){

				//Loop through the event and find an anti-muon
				for(Long64_t j=0;j<nentries;j++){

					//if entry is anti-muon
					if(v->id[j] == -13){

						pt2 = pt(v, j);
						p2 = p(v, j);
						eta2 = eta(v, j);

									//Reconstructing invariant mass of all muon-antimuon-photon combinations and di-muons
									mu_invmass = invmass2(v, i, j);

									if((v->mother1[i] == v->mother1[j]) && (v->mother2[i] == v->mother2[j]) && (v->motherid1[i] == 221)){

										++N->S[0];

										if((pt1>mu_ptcut) && (p1>mu_pcut) && (eta1>eta_cut1) && (eta1<eta_cut2) && (pt2>mu_ptcut) &&  (p2>mu_pcut) && (eta2>eta_cut1) && (eta2<eta_cut2)){

											++N->S[1];

											if((mu_invmass < massA+dmassA) && (mu_invmass> massA-dmassA)){

												++N->S[2];										

											}

										}

									}
									else{

										++N->B[0];

										if((pt1>mu_ptcut) && (p1>mu_pcut) && (eta1>eta_cut1) && (eta1<eta_cut2) && (pt2>mu_ptcut) &&  (p2>mu_pcut) && (eta2>eta_cut1) && (eta2<eta_cut2)){

											++N->B[1];

											if((mu_invmass < massA+dmassA) && (mu_invmass> massA-dmassA)){

												++N->B[2];												

											}

										}

									}

						}	

					}


		}

	}

	for(int i=0;i<3;i++){
		vmu->S[i] = vmu->S[i] + N->S[i];
		vmu->B[i] = vmu->B[i] + N->B[i];
	}

	delete N;

}

void analyze_event_pi(FILE *file_out1, vect *v, vect1 *vmu, double mu_ptcut, double mu_pcut, double eta_cut1, double eta_cut2, double massA, double dmassA, double alpha, double beta, double gamma){

	double inv_mass, pion_invmass, pt_pp, pt_pn, p_pp, p_pn, eta_pp, eta_pn, pt3, p3, eta3;
	double misid1, misid2;

	Long64_t nentries = v->index.size();

	vect1* N = new vect1;
	for(int i=0;i<3;i++){
		N->S.push_back(0);
		N->B.push_back(0);
		N->misID.push_back(0);
	}

	double m_eta = 0.54785;
	double dm_eta = 0.02;

	for(Long64_t i=0;i<nentries;i++){ 

		//if entry is pi+
		if(v->id[i] == 211){

			pt_pp = pt(v, i);
			p_pp = p(v, i);
			eta_pp = eta(v, i);

			if((pt_pp>mu_ptcut) && (p_pp>mu_pcut) && (eta_pp>eta_cut1) && (eta_pp<eta_cut2)){

				//Loop through the event and find a pi-
				for(Long64_t j=0;j<nentries;j++){

					//if entry is anti-muon
					if(v->id[j] == -211){

						pt_pn = pt(v, j);
						p_pn = p(v, j);
						eta_pn = eta(v, j);

						if((pt_pn>mu_ptcut) && (p_pn>mu_pcut) && (eta_pn>eta_cut1) && (eta_pn<eta_cut2)){
					
									//check that pion pair is not a K^0_S (K short) resonance 
									if((v->motherid1[i] != 310) && (v->motherid1[j] != 310)){

										//Reconstructing invariant mass of all pi+/- & gamma combinations
										pion_invmass = invmass2(v, i, j);
								
										misid1 = misID_rate(abs(v->pz[i]), alpha, beta, gamma);
										misid2 = misID_rate(abs(v->pz[j]), alpha, beta, gamma);

										++N->misID[0];
								
										//if((pion_invmass < m_eta+dm_eta) && (pion_invmass > m_eta-dm_eta)){

											++N->misID[1];

											if((pion_invmass < massA+dmassA) && (pion_invmass> massA-dmassA)){

												N->misID[2] = N->misID[2] + misid1 * misid2; 																						
											}

										//}

									}
								

						}

					}

				}

			}

		}

	}

	for(int i=0;i<3;i++){
		vmu->misID[i] = vmu->misID[i] + N->misID[i];
	}

	delete N;

}

int main() {
/*
	//Loading lhcbStyle for plotting
	gROOT->ProcessLine(".L lhcbstyle.C");
	lhcbStyle();

	// Example of adding stat box - turned off by default in plots
	gStyle->SetOptStat("emr");  // show only nent - e , mean - m , rms - r
*/

	// Open the input TFile 
	TFile  *file_in1 = new TFile("gs_project24.root", "READ");
	TFile  *file_in2 = new TFile("/disk/moose/general/user72/friday0203/gs_project26.root", "READ");

	// Open the output TFile 
	TFile  *file_out = new TFile("gs_project96.root", "recreate");

	FILE  *file_out1 = fopen("results_out96.txt", "w");
	fprintf(file_out1, "Table showing signal significance of signal muons passing various cuts: \n");
	fprintf(file_out1, "S(gamma cuts(pt+p+eta))   B_EM(gamma cuts(pt+p+eta))   B(gamma cuts(pt+p+eta))   misID(gamma cuts(pt+p+eta))   S/B_EM(gamma cuts(pt+p+eta))  S/sqrt(B+misID)(gamma cuts(pt+p+eta))  \n");

	double alpha, beta, gamma;

	std::ifstream file("results_out94.txt");
    	if(file.is_open()){

		file >> alpha >> beta >> gamma;

	}

	//Get trees from file
	TTree *T1 = (TTree*)file_in1->Get("T1");
	TTree *T2 = (TTree*)file_in1->Get("T2");
	TTree *T3 = (TTree*)file_in1->Get("T3");

	TTree *T4 = (TTree*)file_in2->Get("T1");
	TTree *T5 = (TTree*)file_in2->Get("T2");
	
	double index_var, id_var, energy_var, mass_var, px_var, py_var, pz_var, mother1_var, motherid2_var, motherid1_var, mother2_var;

	//Get branches of the trees
	T1->SetBranchAddress("index",&index_var);
   	T1->SetBranchAddress("id",&id_var);
	T1->SetBranchAddress("energy",&energy_var);
   	T1->SetBranchAddress("mass",&mass_var);
	T1->SetBranchAddress("px",&px_var);
   	T1->SetBranchAddress("py",&py_var);
	T1->SetBranchAddress("pz",&pz_var);
	T1->SetBranchAddress("mother1",&mother1_var);
	T1->SetBranchAddress("mother2",&mother2_var);
	T1->SetBranchAddress("motherid1",&motherid1_var);
	T1->SetBranchAddress("motherid2",&motherid2_var);

	T2->SetBranchAddress("index",&index_var);
   	T2->SetBranchAddress("id",&id_var);
	T2->SetBranchAddress("energy",&energy_var);
   	T2->SetBranchAddress("mass",&mass_var);
	T2->SetBranchAddress("px",&px_var);
   	T2->SetBranchAddress("py",&py_var);
	T2->SetBranchAddress("pz",&pz_var);
	T2->SetBranchAddress("mother1",&mother1_var);
	T2->SetBranchAddress("mother2",&mother2_var);
	T2->SetBranchAddress("motherid1",&motherid1_var);
	T2->SetBranchAddress("motherid2",&motherid2_var);

	T3->SetBranchAddress("index",&index_var);
   	T3->SetBranchAddress("id",&id_var);
	T3->SetBranchAddress("energy",&energy_var);
   	T3->SetBranchAddress("mass",&mass_var);
	T3->SetBranchAddress("px",&px_var);
   	T3->SetBranchAddress("py",&py_var);
	T3->SetBranchAddress("pz",&pz_var);
	T3->SetBranchAddress("mother1",&mother1_var);
	T3->SetBranchAddress("mother2",&mother2_var);
	T3->SetBranchAddress("motherid1",&motherid1_var);
	T3->SetBranchAddress("motherid2",&motherid2_var);


	T4->SetBranchAddress("index",&index_var);
   	T4->SetBranchAddress("id",&id_var);
	T4->SetBranchAddress("energy",&energy_var);
   	T4->SetBranchAddress("mass",&mass_var);
	T4->SetBranchAddress("px",&px_var);
   	T4->SetBranchAddress("py",&py_var);
	T4->SetBranchAddress("pz",&pz_var);
	T4->SetBranchAddress("mother1",&mother1_var);
	T4->SetBranchAddress("mother2",&mother2_var);
	T4->SetBranchAddress("motherid1",&motherid1_var);
	T4->SetBranchAddress("motherid2",&motherid2_var);

	T5->SetBranchAddress("index",&index_var);
   	T5->SetBranchAddress("id",&id_var);
	T5->SetBranchAddress("energy",&energy_var);
   	T5->SetBranchAddress("mass",&mass_var);
	T5->SetBranchAddress("px",&px_var);
   	T5->SetBranchAddress("py",&py_var);
	T5->SetBranchAddress("pz",&pz_var);
	T5->SetBranchAddress("mother1",&mother1_var);
	T5->SetBranchAddress("mother2",&mother2_var);
	T5->SetBranchAddress("motherid1",&motherid1_var);
	T5->SetBranchAddress("motherid2",&motherid2_var);

	

	vect* v= new vect;//declares a vector in which we are going to store all the particles from the same event
	double prev_index = 0.0;//saving the event index of the previous entry of the tree

	vect1* vmu = new vect1;
	for(int i=0;i<3;i++){
		vmu->S.push_back(0);
		vmu->SMB.push_back(0);
		vmu->B.push_back(0);
		vmu->misID.push_back(0);
	}

	/*
	vect2* vgamma = new vect2;
	int column = 80;
	int row = 100;
	for (int i=0;i<column;i++){
		for (int j=0;j<row;j++){
			vgamma->S[i][j] = 0;
			vgamma->SMB[i][j] = 0;
			vgamma->B[i][j] = 0;
			vgamma->misID[i][j] = 0;
		}
	}
	*/


	double mu_ptcut = 0.5;
	double mu_pcut = 10.0;
	double eta_cut1 = 2.0;
	double eta_cut2 = 4.5;
	//double gamma_ptcut = 5.0;
	//double gamma_pcut = 10.0;
	//double gamma_ptcut_step = 0.01;
	//double gamma_pcut_step = 0.1;


/*
	TH2D *ptvsp = new TH2D("ptvsp","S/#sqrt{B_{EM}+B+misID} based on p_{t} vs p cuts", column, 0.0, gamma_ptcut_step*column, row, 0.0, gamma_pcut_step*row);
	ptvsp -> GetXaxis()-> SetTitle("p_{t} (GeV)");
	ptvsp -> GetYaxis()-> SetTitle("p (GeV)");

	TH2D *SBem_ptvsp = new TH2D("SBem_ptvsp","S/B_{EM} based on p_{t} vs p cuts", column, 0.0, gamma_ptcut_step*column, row, 0.0, gamma_pcut_step*row);
	SBem_ptvsp -> GetXaxis()-> SetTitle("p_{t} (GeV)");
	SBem_ptvsp -> GetYaxis()-> SetTitle("p (GeV)");

	TH2D *S_ptvsp = new TH2D("S_ptvsp","S_{cut}/S_{total} based on p_{t} vs p cuts", column, 0.0, gamma_ptcut_step*column, row, 0.0, gamma_pcut_step*row);
	S_ptvsp -> GetXaxis()-> SetTitle("p_{t} (GeV)");
	S_ptvsp -> GetYaxis()-> SetTitle("p (GeV)");

	TH2D *Bem_ptvsp = new TH2D("Bem_ptvsp","B_{EM cut}/B_{EM total} based on p_{t} vs p cuts", column, 0.0, gamma_ptcut_step*column, row, 0.0, gamma_pcut_step*row);
	Bem_ptvsp -> GetXaxis()-> SetTitle("p_{t} (GeV)");
	Bem_ptvsp -> GetYaxis()-> SetTitle("p (GeV)");

	TH2D *B_ptvsp = new TH2D("B_ptvsp","B_{cut}/B_{total} based on p_{t} vs p cuts", column, 0.0, gamma_ptcut_step*column, row, 0.0, gamma_pcut_step*row);
	B_ptvsp -> GetXaxis()-> SetTitle("p_{t} (GeV)");
	B_ptvsp -> GetYaxis()-> SetTitle("p (GeV)");
*/
	double epsilon2 = 1e-11;//kinetic mixing coupling
	double massA = 0.27;//dark photon mass
	double dmassA = 0.02;//half the mass window around the A' mass

	//Loop through the entries of the tree
	Long64_t nentries = T2->GetEntries();
std::cout<<"Total number of entries: "<<nentries<<std::endl;
   	for (Long64_t i=0;i<nentries;i++){

	      	T2->GetEntry(i);

		//checks if the entry is the very last one in the tree
		if(i+1==nentries){

			//add last particle to the event vector
			v->index.push_back(index_var);
			v->id.push_back(id_var);
			v->energy.push_back(energy_var);
			v->mass.push_back(mass_var);
			v->px.push_back(px_var);
			v->py.push_back(py_var);
			v->pz.push_back(pz_var);
			v->mother1.push_back(mother1_var);
			v->mother2.push_back(mother2_var);
			v->motherid1.push_back(motherid1_var);
			v->motherid2.push_back(motherid2_var);

			
			//analyze_event1(file_out1, v, vmu, vgamma, column, row, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut0, gamma_pcut0, gamma_ptcut_step, gamma_pcut_step, massA, dmassA);
			analyze_event1(file_out1, v, vmu, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, massA, dmassA);

			//clear vector
			v->index.clear();
			v->id.clear();
			v->energy.clear();
			v->mass.clear();
			v->px.clear();
			v->py.clear();
			v->pz.clear();
			v->mother1.clear();
			v->mother2.clear();
			v->motherid1.clear();
			v->motherid2.clear();

		}

		//checks if the new entry is from the same event as the previous one
		else if(prev_index!=index_var){


			//analyze_event1(file_out1, v, vmu, vgamma, column, row, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut0, gamma_pcut0, gamma_ptcut_step, gamma_pcut_step, massA, dmassA);
			analyze_event1(file_out1, v, vmu, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, massA, dmassA);

			//clear vector
			v->index.clear();
			v->id.clear();
			v->energy.clear();
			v->mass.clear();
			v->px.clear();
			v->py.clear();
			v->pz.clear();
			v->mother1.clear();
			v->mother2.clear();
			v->motherid1.clear();
			v->motherid2.clear();

			//add new particle to the event vector
			v->index.push_back(index_var);
			v->id.push_back(id_var);
			v->energy.push_back(energy_var);
			v->mass.push_back(mass_var);
			v->px.push_back(px_var);
			v->py.push_back(py_var);
			v->pz.push_back(pz_var);
			v->mother1.push_back(mother1_var);
			v->mother2.push_back(mother2_var);
			v->motherid1.push_back(motherid1_var);
			v->motherid2.push_back(motherid2_var);

			prev_index = index_var;//setting the index for current event

		}
		else{

			//add new particle to the event vector
			v->index.push_back(index_var);
			v->id.push_back(id_var);
			v->energy.push_back(energy_var);
			v->mass.push_back(mass_var);
			v->px.push_back(px_var);
			v->py.push_back(py_var);
			v->pz.push_back(pz_var);
			v->mother1.push_back(mother1_var);
			v->mother2.push_back(mother2_var);
			v->motherid1.push_back(motherid1_var);
			v->motherid2.push_back(motherid2_var);

		}
	      	
   	}

	prev_index = 0.0;

	//Loop through the entries of the tree
	nentries = T5->GetEntries();
std::cout<<"Total number of entries: "<<nentries<<std::endl;
   	for (Long64_t i=0;i<nentries;i++){

	      	T5->GetEntry(i);

		//checks if the entry is the very last one in the tree
		if(i+1==nentries){

			//add last particle to the event vector
			v->index.push_back(index_var);
			v->id.push_back(id_var);
			v->energy.push_back(energy_var);
			v->mass.push_back(mass_var);
			v->px.push_back(px_var);
			v->py.push_back(py_var);
			v->pz.push_back(pz_var);
			v->mother1.push_back(mother1_var);
			v->mother2.push_back(mother2_var);
			v->motherid1.push_back(motherid1_var);
			v->motherid2.push_back(motherid2_var);

			
			//analyze_event2(file_out1, v, vmu, vgamma, column, row, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut0, gamma_pcut0, gamma_ptcut_step, gamma_pcut_step, massA, dmassA);
			analyze_event2(file_out1, v, vmu, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, massA, dmassA);

			//clear vector
			v->index.clear();
			v->id.clear();
			v->energy.clear();
			v->mass.clear();
			v->px.clear();
			v->py.clear();
			v->pz.clear();
			v->mother1.clear();
			v->mother2.clear();
			v->motherid1.clear();
			v->motherid2.clear();

		}

		//checks if the new entry is from the same event as the previous one
		else if(prev_index!=index_var){


			//analyze_event2(file_out1, v, vmu, vgamma, column, row, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut0, gamma_pcut0, gamma_ptcut_step, gamma_pcut_step, massA, dmassA);
			analyze_event2(file_out1, v, vmu, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, massA, dmassA);

			//clear vector
			v->index.clear();
			v->id.clear();
			v->energy.clear();
			v->mass.clear();
			v->px.clear();
			v->py.clear();
			v->pz.clear();
			v->mother1.clear();
			v->mother2.clear();
			v->motherid1.clear();
			v->motherid2.clear();

			//add new particle to the event vector
			v->index.push_back(index_var);
			v->id.push_back(id_var);
			v->energy.push_back(energy_var);
			v->mass.push_back(mass_var);
			v->px.push_back(px_var);
			v->py.push_back(py_var);
			v->pz.push_back(pz_var);
			v->mother1.push_back(mother1_var);
			v->mother2.push_back(mother2_var);
			v->motherid1.push_back(motherid1_var);
			v->motherid2.push_back(motherid2_var);

			prev_index = index_var;//setting the index for current event

		}
		else{

			//add new particle to the event vector
			v->index.push_back(index_var);
			v->id.push_back(id_var);
			v->energy.push_back(energy_var);
			v->mass.push_back(mass_var);
			v->px.push_back(px_var);
			v->py.push_back(py_var);
			v->pz.push_back(pz_var);
			v->mother1.push_back(mother1_var);
			v->mother2.push_back(mother2_var);
			v->motherid1.push_back(motherid1_var);
			v->motherid2.push_back(motherid2_var);

		}
	      	
   	}

	prev_index = 0.0;

	//Loop through the entries of tree T3
	nentries = T3->GetEntries();
std::cout<<"Total number of entries T3: "<<nentries<<std::endl;
   	for (Long64_t i=0;i<nentries;i++){

	      	T3->GetEntry(i);

		//if(index_var>=10000){break;}

		//checks if the entry is the very last one in the tree
		if(i+1==nentries){

			//add last particle to the event vector
			v->index.push_back(index_var);
			v->id.push_back(id_var);
			v->energy.push_back(energy_var);
			v->mass.push_back(mass_var);
			v->px.push_back(px_var);
			v->py.push_back(py_var);
			v->pz.push_back(pz_var);
			v->mother1.push_back(mother1_var);
			v->mother2.push_back(mother2_var);
			v->motherid1.push_back(motherid1_var);
			v->motherid2.push_back(motherid2_var);


			//analyze_event_pi(file_out1, v, vmu, vgamma, column, row, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut0, gamma_pcut0, gamma_ptcut_step, gamma_pcut_step, massA, dmassA, alpha, beta, gamma);
			analyze_event_pi(file_out1, v, vmu, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, massA, dmassA, alpha, beta, gamma);

			//clear vector
			v->index.clear();
			v->id.clear();
			v->energy.clear();
			v->mass.clear();
			v->px.clear();
			v->py.clear();
			v->pz.clear();
			v->mother1.clear();
			v->mother2.clear();
			v->motherid1.clear();
			v->motherid2.clear();

		}

		//checks if the new entry is from the same event as the previous one
		else if(prev_index!=index_var){


			//analyze_event_pi(file_out1, v, vmu, vgamma, column, row, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut0, gamma_pcut0, gamma_ptcut_step, gamma_pcut_step, massA, dmassA, alpha, beta, gamma);
			analyze_event_pi(file_out1, v, vmu, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, massA, dmassA, alpha, beta, gamma);

			//clear vector
			v->index.clear();
			v->id.clear();
			v->energy.clear();
			v->mass.clear();
			v->px.clear();
			v->py.clear();
			v->pz.clear();
			v->mother1.clear();
			v->mother2.clear();
			v->motherid1.clear();
			v->motherid2.clear();

			//add new particle to the event vector
			v->index.push_back(index_var);
			v->id.push_back(id_var);
			v->energy.push_back(energy_var);
			v->mass.push_back(mass_var);
			v->px.push_back(px_var);
			v->py.push_back(py_var);
			v->pz.push_back(pz_var);
			v->mother1.push_back(mother1_var);
			v->mother2.push_back(mother2_var);
			v->motherid1.push_back(motherid1_var);
			v->motherid2.push_back(motherid2_var);

			prev_index = index_var;//setting the index for current event

		}
		else{

			//add new particle to the event vector
			v->index.push_back(index_var);
			v->id.push_back(id_var);
			v->energy.push_back(energy_var);
			v->mass.push_back(mass_var);
			v->px.push_back(px_var);
			v->py.push_back(py_var);
			v->pz.push_back(pz_var);
			v->mother1.push_back(mother1_var);
			v->mother2.push_back(mother2_var);
			v->motherid1.push_back(motherid1_var);
			v->motherid2.push_back(motherid2_var);

		}
	      	
   	}

	double sigma, gamma_ptcut, gamma_pcut;
	double width = Gamma_A(epsilon2, massA);//calculating the width of the dark photon

	double F = scaling(epsilon2, massA, width); 
	double Br = 0.0000060;
	double RmisID = 1e-6;

	double S, SMB, B, misID;

	
			S = (double)vmu->S[2] * F * Br * (double)vmu->SMB[0] /(double)vmu->S[0];
			SMB = (double)vmu->SMB[2] * Br;
			B = (double)vmu->B[2] * Br;
			misID = (double)vmu->misID[2];

			sigma = metric(S, SMB, B, misID);

			fprintf(file_out1, "%.5e   %.5e   %.5e   %.5e   %.5e   %.5e \n",  S, SMB, B, misID, S/SMB, sigma);
	


	std::cout<<"Signal entries (no cuts): "<<vmu->S[0]<<std::endl;
	std::cout<<"SM Background entries (no cuts): "<<vmu->SMB[0]<<std::endl;
	std::cout<<"Background entries (no cuts): "<<vmu->B[0]<<std::endl;
	std::cout<<"MisID pion entries (mu cuts + mass cuts): "<<vmu->misID[0]<<std::endl;

	fprintf(file_out1, "Candidate type    Total entries    After eta mass cut    After A mass cut: \n");
	fprintf(file_out1, "S    			   %d    		%d   			%d \n", vmu->S[0], vmu->S[1], vmu->S[2]);
	fprintf(file_out1, "SMB    		   %d    		%d    			%d \n", vmu->SMB[0], vmu->SMB[1], vmu->SMB[2]);
	fprintf(file_out1, "B    			   %d    		%d    			%d \n", vmu->B[0], vmu->B[1], vmu->B[2]);
	fprintf(file_out1, "misID    		   %d    		%d    			%d \n", vmu->misID[0], vmu->misID[1], vmu->misID[2]);
	

	file_out->Write();
	file_out->Close();	
        file_in1->Close();
	file_in2->Close();
	delete file_out;
	file.close();

	fclose(file_out1);
	delete v;
	delete vmu;
	// Done.
	return 0;
}


