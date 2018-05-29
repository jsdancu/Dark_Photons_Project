//13/02/2018
//Program reading in the three trees generated in gs_project92.cc, gs_project34.cc and gs_project62.cc and exploring the efficiency/feasibility of cuts on photon pt and p

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

#include "Pythia8/Pythia.h"
using namespace Pythia8;

using namespace std;

//#include "lhcbStyle.C"

struct vect{
	
	std::vector<double> index, id, energy, mass, px, py, pz, x, y, z, mother1, mother2, motherid1, motherid2;
	
};

struct vect1{
	
	std::vector<int> S, B, EMB, misID;
	
};

struct vect2{
	
	//int S1[80][100]; int S2[80][100]; int B[80][100]; int EMB[80][100]; int misID[80][100];
	int S1[16][20]; int S2[16][20]; int B[16][20]; int EMB[16][20]; double misID[16][20];
	
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

double rt(vect *v, Long64_t i){

	double rt1=sqrt(pow(v->x[i],2)+pow(v->y[i],2));

	return rt1;

}

double p_1(double px, double py, double pz){

	double p1=sqrt(pow(px,2)+pow(py,2)+pow(pz,2));

	return p1;

}

double pt_1(double px, double py){

	double pt1=sqrt(pow(px,2)+pow(py,2));

	return pt1;

}

double eta_1(double px, double py, double pz){

	double eta=0.5 * log((p_1(px, py, pz) + abs(pz)) / (p_1(px, py, pz) - abs(pz)));

	return eta;

}

double rt_1(double x, double y){

	double rt1=sqrt(pow(x,2)+pow(x,2));

	return rt1;

}

double rtg(vect *v, Long64_t i){

	double rt1=sqrt(pow(v->x[i],2)+pow(v->y[i],2));

	return rt1;

}

double misEnergy(double p){

	double m0= 0.10566;
	double energy=sqrt(pow(p, 2) + pow(m0, 2));
	return energy;

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

//double metric(double S, double EMB, double B, double misID){
double metric(double S, double misID){

	double sigma = S/sqrt(misID);
	//double sigma = S/sqrt(EMB + B + misID);

	return sigma;

}

void analyze_event1(FILE *file_out1, vect *v, vect1 *vmu, vect2 *vgamma, int column, int row, double mu_ptcut, double mu_pcut, double eta_cut1, double eta_cut2, double gamma_ptcut0, double gamma_pcut0, double gamma_ptcut_step, double gamma_pcut_step, double massA, double dmassA, double fiducial_min, double fiducial_max){

	double inv_mass, mu_invmass, pt1, pt2, pt3, p1, p2, p3, eta1, eta2, eta3, rt1, rt2, rt3, gamma_ptcut, gamma_pcut;

	Long64_t nentries = v->index.size();

	vect1* N = new vect1;
	for(int i=0;i<4;i++){

		N->EMB.push_back(0);
		N->B.push_back(0);

	}

	vect2* M = new vect2;
	for (int i=0;i<column;i++){
		for (int j=0;j<row;j++){

			M->EMB[i][j] = 0;
			M->B[i][j] = 0;
	
		}
	}

	double m_eta = 0.54785;
	double dm_eta = 0.02;

	for(Long64_t i=0;i<nentries;i++){ 

		//if entry is muon
		if(v->id[i] == 13){

			pt1 = pt(v, i);
			p1 = p(v, i);
			eta1 = eta(v, i);
			rt1 = rt(v, i);

			//if((pt1>mu_ptcut) && (p1>mu_pcut) && (eta1>eta_cut1) && (eta1<eta_cut2)){

				//Loop through the event and find an anti-muon
				for(Long64_t j=0;j<nentries;j++){

					//if entry is anti-muon
					if(v->id[j] == -13){

						pt2 = pt(v, j);
						p2 = p(v, j);
						eta2 = eta(v, j);
						rt2 = rt(v, j);

						//if((pt2>mu_ptcut) &&  (p2>mu_pcut) && (eta2>eta_cut1) && (eta2<eta_cut2)){

							for(Long64_t k=0;k<nentries;k++){

								//if entry is photon
								if(v->id[k] == 22){

									pt3 = pt(v, k);							
									p3 = p(v, k);
									eta3 = eta(v, k);
									rt3 = rt(v, k);	

									//Reconstructing invariant mass of all muon-antimuon-photon combinations and di-muons
									inv_mass = invmass(v, i, j, k);
									mu_invmass = invmass2(v, i, j);

									if((v->mother1[i] == v->mother1[j]) && (v->mother1[i] == v->mother1[k]) && (v->mother2[i] == v->mother2[j]) && (v->mother2[i] == v->mother2[k]) && (v->motherid1[i] == 221)){									

										++N->EMB[0];

										if((pt1>mu_ptcut) && (p1>mu_pcut) && (eta1>eta_cut1) && (eta1<eta_cut2) && (pt2>mu_ptcut) &&  (p2>mu_pcut) && (eta2>eta_cut1) && (eta2<eta_cut2) && (mu_invmass < massA+dmassA) && (mu_invmass> massA-dmassA)){ 

											++N->EMB[1];
											
											//if((rt1 > fiducial_min) && (rt1 < fiducial_max) && (rt2 > fiducial_min) && (rt2 < fiducial_max)){

												++N->EMB[2];

												if((inv_mass < m_eta+dm_eta) && (inv_mass> m_eta-dm_eta)){

													++N->EMB[3];

													if((eta3>eta_cut1) && (eta3<eta_cut2)){

														for(int p=0; p<column; p++){

															gamma_ptcut = gamma_ptcut0 + gamma_ptcut_step * p;

															if(pt3>gamma_ptcut){

																for(int q=0; q<row; q++){

																	gamma_pcut = gamma_pcut0 + gamma_pcut_step * q;
																	if(p3>gamma_pcut){++M->EMB[p][q];}

																}

															}												
							
														}

													}

												}

											//}

										}

									}
									else{

										++N->B[0];

										if((pt1>mu_ptcut) && (p1>mu_pcut) && (eta1>eta_cut1) && (eta1<eta_cut2) && (pt2>mu_ptcut) &&  (p2>mu_pcut) && (eta2>eta_cut1) && (eta2<eta_cut2) && (mu_invmass < massA+dmassA) && (mu_invmass> massA-dmassA)){

											++N->B[1];

											//if((rt1 > fiducial_min) && (rt1 < fiducial_max) && (rt2 > fiducial_min) && (rt2 < fiducial_max)){

												++N->B[2];

												if((inv_mass < m_eta+dm_eta) && (inv_mass> m_eta-dm_eta)){

													++N->B[3];

													if((eta3>eta_cut1) && (eta3<eta_cut2)){

														for(int p=0; p<column; p++){

															gamma_ptcut = gamma_ptcut0 + gamma_ptcut_step * p;

															if(pt3>gamma_ptcut){

																for(int q=0; q<row; q++){

																	gamma_pcut = gamma_pcut0 + gamma_pcut_step * q;
																	if(p3>gamma_pcut){++M->B[p][q];}

																}

															}												
							
														}

													}

												}

											//}

										}

									}

								}	

							}

					}

				}

		}

	}

	for(int i=0;i<4;i++){
		vmu->B[i] = vmu->B[i] + N->B[i];
		vmu->EMB[i] = vmu->EMB[i] + N->EMB[i];
	}

	for(int i=0;i<column;i++){ 
		for(int j=0;j<row;j++){ 
	
			vgamma->EMB[i][j] = vgamma->EMB[i][j] + M->EMB[i][j];
/*			if (vgamma->EMB[i][j] < 0) {
			  std::cout << "--------------------------------------\n";
			  std::cout << vgamma->EMB[i][j] << std::endl;
			  std::cout << M->EMB[i][j] << std::endl;
			}*/
			vgamma->B[i][j] = vgamma->B[i][j] + M->B[i][j];
		}
	}

	delete N;
	delete M;

}

void analyze_event2(FILE *file_out1, vect *v, vect1 *vmu, vect2 *vgamma, int column, int row, double mu_ptcut, double mu_pcut, double eta_cut1, double eta_cut2, double gamma_ptcut0, double gamma_pcut0, double gamma_ptcut_step, double gamma_pcut_step, double massA, double dmassA, double fiducial_min, double fiducial_max){

	double inv_mass, mu_invmass, pt1, pt2, pt3, p1, p2, p3, eta1, eta2, eta3, rt1, rt2, rt3, gamma_ptcut, gamma_pcut;

	Long64_t nentries = v->index.size();

	vect1* N = new vect1;
	for(int i=0;i<4;i++){
		N->S.push_back(0);
	}

	vect2* M = new vect2;
	for (int i=0;i<column;i++){
		for (int j=0;j<row;j++){
			M->S1[i][j] = 0;
			M->S2[i][j] = 0;
		}
	}

	double m_eta = 0.54785;
	double dm_eta = 0.02;

	for(Long64_t i=0;i<nentries;i++){ 

		//if entry is muon
		if(v->id[i] == 13){

			pt1 = pt(v, i);
			p1 = p(v, i);
			eta1 = eta(v, i);
			rt1 = rt(v, i);
//std::cout<<"rt1: "<<rt1<<std::endl;
			//if((pt1>mu_ptcut) && (p1>mu_pcut) && (eta1>eta_cut1) && (eta1<eta_cut2)){

				//Loop through the event and find an anti-muon
				for(Long64_t j=0;j<nentries;j++){

					//if entry is anti-muon
					if(v->id[j] == -13){

						pt2 = pt(v, j);
						p2 = p(v, j);
						eta2 = eta(v, j);
						rt2 = rt(v, j);
//std::cout<<"rt2: "<<rt2<<std::endl;
						//if((pt2>mu_ptcut) &&  (p2>mu_pcut) && (eta2>eta_cut1) && (eta2<eta_cut2)){

							for(Long64_t k=0;k<nentries;k++){

								//if entry is photon
								if(v->id[k] == 22){

									pt3 = pt(v, k);							
									p3 = p(v, k);
									eta3 = eta(v, k);	
									rt3 = rtg(v, k);
	
									//Reconstructing invariant mass of all muon-antimuon-photon combinations and di-muons
									inv_mass = invmass(v, i, j, k);
									mu_invmass = invmass2(v, i, j);

									if((v->mother1[i] == v->mother1[j]) && (v->mother1[i] == v->mother1[k]) && (v->mother2[i] == v->mother2[j]) && (v->mother2[i] == v->mother2[k]) && (v->motherid1[i] == 221)){

										++N->S[0];

										if((pt1>mu_ptcut) && (p1>mu_pcut) && (eta1>eta_cut1) && (eta1<eta_cut2) && (pt2>mu_ptcut) &&  (p2>mu_pcut) && (eta2>eta_cut1) && (eta2<eta_cut2) && (mu_invmass < massA+dmassA) && (mu_invmass> massA-dmassA)){

											++N->S[1];

											if((rt1 > fiducial_min) && (rt1 < fiducial_max) && (rt2 > fiducial_min) && (rt2 < fiducial_max)){++N->S[2];}											
												
												if((inv_mass < m_eta+dm_eta) && (inv_mass> m_eta-dm_eta)){

													++N->S[3];

													if((eta3>eta_cut1) && (eta3<eta_cut2)){

														for(int p=0; p<column; p++){

															gamma_ptcut = gamma_ptcut0 + gamma_ptcut_step * p;

															if(pt3>gamma_ptcut){

																for(int q=0; q<row; q++){

																	gamma_pcut = gamma_pcut0 + gamma_pcut_step * q;
																	if(p3>gamma_pcut){

																		++M->S1[p][q];

																		if((rt1 > fiducial_min) && (rt1 < fiducial_max) && (rt2 > fiducial_min) && (rt2 < fiducial_max)){++M->S2[p][q];}
															
																	}

																}

															}												
							
														}

													}

												//}

											}

										}

									}									

								}	

							}

					}

				}

		}

	}

	for(int i=0;i<4;i++){
		vmu->S[i] = vmu->S[i] + N->S[i];
	}

	for(int i=0;i<column;i++){ 
		for(int j=0;j<row;j++){ 
			vgamma->S1[i][j] = vgamma->S1[i][j] + M->S1[i][j];
			vgamma->S2[i][j] = vgamma->S2[i][j] + M->S2[i][j];
		}
	}

	delete N;
	delete M;

}

void analyze_event_K(FILE *file_out1, vect *v, vect *v1, vect1 *vmu, vect2 *vgamma, int column, int row, double mu_ptcut, double mu_pcut, double eta_cut1, double eta_cut2, double gamma_ptcut0, double gamma_pcut0, double gamma_ptcut_step, double gamma_pcut_step, double massA, double dmassA, double fiducial_min, double fiducial_max, double alpha, double beta, double gamma, double redecay, TH1D *KS_number_event, Pythia &pythia1){

	//Create event
	Event *event1 = &pythia1.event;

	double inv_mass, pion_invmass,  pt1, pt2, pt3, p1, p2, p3, eta1, eta2, eta3, rt1, rt2, rt3, gamma_ptcut, gamma_pcut;
	double misid1, misid2, misEnergy1, misEnergy2;

	Long64_t KS_number = 0;

	Long64_t nentries = v->index.size();
	Long64_t nentries1 = v1->index.size();

	vect1* N = new vect1;
	for(int i=0;i<4;i++){
		N->B.push_back(0);
		N->misID.push_back(0);
	}

	vect2* M = new vect2;
	for (int i=0;i<column;i++){
		for (int j=0;j<row;j++){
			M->B[i][j] = 0;
			M->misID[i][j] = 0;
		}
	}

	double m_eta = 0.54785;
	double dm_eta = 0.02;

	for(Long64_t iEntries=0;iEntries<nentries;iEntries++){ 

		//if entry is KS
		if(v->id[iEntries] == 310){

			++KS_number;

			for(int jEntries = 0; jEntries <redecay; ++jEntries){

				pythia1.event.reset(); 
				pythia1.event.append(v->id[iEntries], 23, 0, 0, 0, 0, 0, 0, v->px[iEntries], v->py[iEntries], v->pz[iEntries], 0, pythia1.particleData.mSel(310)); 
				pythia1.event[1].e(pythia1.event[1].eCalc());
				pythia1.event[1].vProd(v->x[iEntries], v->y[iEntries], v->z[iEntries], 0.0);
				pythia1.event[1].tau(pythia1.event[1].tau0() * pythia1.rndm.exp() );

				//Generate one event. Skip if error.
				if (!pythia1.next()) continue;

				//if (i < 2) {pythia1.info.list(); pythia1.event.list();}

				for (Long64_t i = 0; i < pythia1.event.size(); ++i) {

					if((pythia1.event[i].isFinal()) && (pythia1.event[i].id() == 211)){

						pt1 = pt_1(pythia1.event[i].px(), pythia1.event[i].py());
						p1 = p_1(pythia1.event[i].px(), pythia1.event[i].py(), pythia1.event[i].pz());
						eta1 = eta_1(pythia1.event[i].px(), pythia1.event[i].py(), pythia1.event[i].pz());
						rt1 = rt_1(pythia1.event[i].xProd(), pythia1.event[i].yProd());

						misEnergy1 = misEnergy(p1);

						misid1 = misID_rate(abs(pythia1.event[i].pz()), alpha, beta, gamma);

						if((pt1>mu_ptcut) && (p1>mu_pcut) && (eta1>eta_cut1) && (eta1<eta_cut2)){

							for (Long64_t j = 0; j < pythia1.event.size(); ++j) {

								if((pythia1.event[j].isFinal()) && (pythia1.event[j].id() == -211)){

									pt2 = pt_1(pythia1.event[j].px(), pythia1.event[j].py());
									p2 = p_1(pythia1.event[j].px(), pythia1.event[j].py(), pythia1.event[j].pz());
									eta2 = eta_1(pythia1.event[j].px(), pythia1.event[j].py(), pythia1.event[j].pz());
									rt2 = rt_1(pythia1.event[j].xProd(), pythia1.event[j].yProd());

									misEnergy2 = misEnergy(p2);

									misid2 = misID_rate(abs(pythia1.event[j].pz()), alpha, beta, gamma);

									if((pt2>mu_ptcut) && (p2>mu_pcut) && (eta2>eta_cut1) && (eta2<eta_cut2)){

										for(Long64_t k=0;k<nentries1;k++){

											//if entry is photon
											if(v1->id[k] == 22){

													//Reconstructing invariant mass of all pi+/- & gamma combinations
													inv_mass = sqrt(pow(misEnergy1+misEnergy2+v1->energy[k], 2) - pow(pythia1.event[i].px()+pythia1.event[j].px()+v1->px[k], 2) - pow(pythia1.event[i].py()+pythia1.event[j].py()+v1->py[k], 2) - pow(pythia1.event[i].pz()+pythia1.event[j].pz()+v1->pz[k], 2));
													pion_invmass = sqrt(pow(misEnergy1+misEnergy2, 2) - pow(pythia1.event[i].px()+pythia1.event[j].px(), 2) - pow(pythia1.event[i].py()+pythia1.event[j].py(), 2) - pow(pythia1.event[i].pz()+pythia1.event[j].pz(), 2));

													pt3 = pt(v1, k);
													p3 = p(v1, k);
													eta3 = eta(v1, k);
													rt3 = rt(v1, k);													

													++N->misID[0];
																		
													if((pion_invmass < massA+dmassA) && (pion_invmass> massA-dmassA)){

														++N->misID[1];

														if((rt1 > fiducial_min) && (rt1 < fiducial_max) && (rt2 > fiducial_min) && (rt2 < fiducial_max)){

															++N->misID[2];
												
															if((inv_mass < m_eta+dm_eta) && (inv_mass > m_eta-dm_eta)){

																++N->misID[3];
											
																if((eta3>eta_cut1) && (eta3<eta_cut2)){

																	for(int p=0; p<column; p++){

																		gamma_ptcut = gamma_ptcut0 + gamma_ptcut_step * p;

																		if(pt3>gamma_ptcut){

																			for(int q=0; q<row; q++){

																				gamma_pcut = gamma_pcut0 + gamma_pcut_step * q;
																				if(p3>gamma_pcut){M->misID[p][q] = M->misID[p][q] + misid1 * misid2;}

																			}

																		}												
							
																	}

																}

															}

														}

													}

											}

										}

									}

								}

							}

						}

					}

				}

			}

		}

	}

	if (v->index.size() > 0) {

		KS_number_event->SetBinContent((v->index[1])+1, KS_number);

	}

	for(int i=0;i<4;i++){
		vmu->misID[i] = vmu->misID[i] + N->misID[i];
	}

	for(int i=0;i<column;i++){ 
		for(int j=0;j<row;j++){ 
	
			vgamma->misID[i][j] = vgamma->misID[i][j] + M->misID[i][j];
			
		}
	}

	delete N;
	delete M;


}


/*
void analyze_event_K(FILE *file_out1, vect *v, vect1 *vmu, vect2 *vgamma, int column, int row, double mu_ptcut, double mu_pcut, double eta_cut1, double eta_cut2, double gamma_ptcut0, double gamma_pcut0, double gamma_ptcut_step, double gamma_pcut_step, double massA, double dmassA, double fiducial_min, double fiducial_max, double alpha, double beta, double gamma){

	double inv_mass, pion_invmass, pt_pp, pt_pn, p_pp, p_pn, eta_pp, eta_pn, pt3, p3, eta3, rt1, rt2, rt3, gamma_ptcut, gamma_pcut;
	double misid1, misid2;

	Long64_t nentries = v->index.size();

	vect1* N = new vect1;
	for(int i=0;i<4;i++){
		N->B.push_back(0);
		N->misID.push_back(0);
	}

	vect2* M = new vect2;
	for (int i=0;i<column;i++){
		for (int j=0;j<row;j++){
			M->B[i][j] = 0;
			M->misID[i][j] = 0;
		}
	}

	double m_eta = 0.54785;
	double dm_eta = 0.02;

	for(Long64_t i=0;i<nentries;i++){ 

		//if entry is pi+
		if(v->id[i] == 211){

			pt_pp = pt(v, i);
			p_pp = p(v, i);
			eta_pp = eta(v, i);
			rt1 = rt(v, i);

			if((pt_pp>mu_ptcut) && (p_pp>mu_pcut) && (eta_pp>eta_cut1) && (eta_pp<eta_cut2)){

				//Loop through the event and find a pi-
				for(Long64_t j=0;j<nentries;j++){

					//if entry is anti-muon
					if(v->id[j] == -211){

						pt_pn = pt(v, j);
						p_pn = p(v, j);
						eta_pn = eta(v, j);
						rt2 = rt(v, j);

						if((pt_pn>mu_ptcut) && (p_pn>mu_pcut) && (eta_pn>eta_cut1) && (eta_pn<eta_cut2)){

							for(Long64_t k=0;k<nentries;k++){

								//if entry is photon
								if(v->id[k] == 22){

										//Reconstructing invariant mass of all pi+/- & gamma combinations
										inv_mass = invmass(v, i, j, k);
										pion_invmass = invmass2(v, i, j);

										pt3 = pt(v, k);
										p3 = p(v, k);
										eta3 = eta(v, k);
										rt3 = rt(v, k);

										misid1 = misID_rate(abs(v->pz[i]), alpha, beta, gamma);
										misid2 = misID_rate(abs(v->pz[j]), alpha, beta, gamma);

										++N->misID[0];
																		
										if((pion_invmass < massA+dmassA) && (pion_invmass> massA-dmassA)){

											++N->misID[1];

											if((rt1 > fiducial_min) && (rt1 < fiducial_max) && (rt2 > fiducial_min) && (rt2 < fiducial_max)){

												++N->misID[2];
												
												if((inv_mass < m_eta+dm_eta) && (inv_mass > m_eta-dm_eta)){

													++N->misID[3];
											
													if((eta3>eta_cut1) && (eta3<eta_cut2)){

														for(int p=0; p<column; p++){

															gamma_ptcut = gamma_ptcut0 + gamma_ptcut_step * p;

															if(pt3>gamma_ptcut){

																for(int q=0; q<row; q++){

																	gamma_pcut = gamma_pcut0 + gamma_pcut_step * q;
																	if(p3>gamma_pcut){M->misID[p][q] = M->misID[p][q] + misid1 * misid2;}

																}

															}												
							
														}

													}

												}

											}

										}

								}

							}

						}

					}

				}

			}

		}

	}

	for(int i=0;i<4;i++){
		vmu->misID[i] = vmu->misID[i] + N->misID[i];
	}

	for(int i=0;i<column;i++){ 
		for(int j=0;j<row;j++){ 
	
			vgamma->misID[i][j] = vgamma->misID[i][j] + M->misID[i][j];
			
		}
	}

	delete N;
	delete M;


}
*/
int main() {
/*
	//Loading lhcbStyle for plotting
	gROOT->ProcessLine(".L lhcbstyle.C");
	lhcbStyle();

	// Example of adding stat box - turned off by default in plots
	gStyle->SetOptStat("emr");  // show only nent - e , mean - m , rms - r
*/

        // Create Pythia instance.
	Pythia pythia1;

	pythia1.readString("Random:setSeed = on");
    	pythia1.readString("Random:seed = 0");
	pythia1.readString("ProcessLevel:all = off");
	pythia1.readString("ParticleDecays:FSRinDecays = off");
	pythia1.readString("Print:quiet = on");

	pythia1.readString("310:onMode = off");
	pythia1.readString("310:onIfMatch = 211 -211");
	pythia1.readString("310:mWidth = 0.1");
	pythia1.readString("310:mMin = 0");
	pythia1.readString("310:mMax = -1");

	pythia1.init();


	// Open the input TFile 
	TFile  *file_in1 = new TFile("/disk/moose/general/user72/gs_project92.root", "READ");
	TFile  *file_in2 = new TFile("/disk/moose/general/user72/gs_project64.root", "READ");
	TFile  *file_in3 = new TFile("/disk/moose/general/user72/gs_project36.root", "READ");

	// Open the output TFile 
	TFile  *file_out = new TFile("gs_project84.root", "recreate");

	double alpha, beta, gamma;

	std::ifstream file("results_out94.txt");
    	if(file.is_open()){

		file >> alpha >> beta >> gamma;

	}

	FILE  *file_out1 = fopen("results_out84.txt", "w");
	fprintf(file_out1, "Table showing signal significance of signal muons passing various cuts: \n");
	fprintf(file_out1, "gamma(pt cut)  gamma(p cut)   S(gamma cuts(pt+p+eta))   B(gamma cuts(pt+p+eta))   S/sqrt(B)(gamma cuts(pt+p+eta))  \n");

	//Get trees from file
	TTree *T1 = (TTree*)file_in1->Get("T1");
	TTree *T2 = (TTree*)file_in1->Get("T2");

	TTree *T4 = (TTree*)file_in2->Get("T1");
	TTree *T5 = (TTree*)file_in2->Get("T2");

	TTree *T6 = (TTree*)file_in3->Get("T2");
	
	double index_var, id_var, energy_var, mass_var, px_var, py_var, pz_var, x_var, y_var, z_var, mother1_var, motherid2_var, motherid1_var, mother2_var, index_var1, id_var1, energy_var1, mass_var1, px_var1, py_var1, pz_var1, x_var1, y_var1, z_var1, mother1_var1, motherid2_var1, motherid1_var1, mother2_var1;

	//Get branches of the trees
	T1->SetBranchAddress("index",&index_var);
   	T1->SetBranchAddress("id",&id_var);
	T1->SetBranchAddress("energy",&energy_var);
   	T1->SetBranchAddress("mass",&mass_var);
	T1->SetBranchAddress("px",&px_var);
   	T1->SetBranchAddress("py",&py_var);
	T1->SetBranchAddress("pz",&pz_var);
	T1->SetBranchAddress("x",&x_var);
   	T1->SetBranchAddress("y",&y_var);
	T1->SetBranchAddress("z",&z_var);
	T1->SetBranchAddress("mother1",&mother1_var);
	T1->SetBranchAddress("mother2",&mother2_var);
	T1->SetBranchAddress("motherid1",&motherid1_var);
	T1->SetBranchAddress("motherid2",&motherid2_var);

	T2->SetBranchAddress("index",&index_var1);
   	T2->SetBranchAddress("id",&id_var1);
	T2->SetBranchAddress("energy",&energy_var1);
   	T2->SetBranchAddress("mass",&mass_var1);
	T2->SetBranchAddress("px",&px_var1);
   	T2->SetBranchAddress("py",&py_var1);
	T2->SetBranchAddress("pz",&pz_var1);
	T2->SetBranchAddress("x",&x_var1);
   	T2->SetBranchAddress("y",&y_var1);
	T2->SetBranchAddress("z",&z_var1);
	T2->SetBranchAddress("mother1",&mother1_var1);
	T2->SetBranchAddress("mother2",&mother2_var1);
	T2->SetBranchAddress("motherid1",&motherid1_var1);
	T2->SetBranchAddress("motherid2",&motherid2_var1);


	T4->SetBranchAddress("index",&index_var);
   	T4->SetBranchAddress("id",&id_var);
	T4->SetBranchAddress("energy",&energy_var);
   	T4->SetBranchAddress("mass",&mass_var);
	T4->SetBranchAddress("px",&px_var);
   	T4->SetBranchAddress("py",&py_var);
	T4->SetBranchAddress("pz",&pz_var);
	T4->SetBranchAddress("x",&x_var);
   	T4->SetBranchAddress("y",&y_var);
	T4->SetBranchAddress("z",&z_var);
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
	T5->SetBranchAddress("x",&x_var);
   	T5->SetBranchAddress("y",&y_var);
	T5->SetBranchAddress("z",&z_var);
	T5->SetBranchAddress("mother1",&mother1_var);
	T5->SetBranchAddress("mother2",&mother2_var);
	T5->SetBranchAddress("motherid1",&motherid1_var);
	T5->SetBranchAddress("motherid2",&motherid2_var);

	T6->SetBranchAddress("index",&index_var);
   	T6->SetBranchAddress("id",&id_var);
	T6->SetBranchAddress("energy",&energy_var);
   	T6->SetBranchAddress("mass",&mass_var);
	T6->SetBranchAddress("px",&px_var);
   	T6->SetBranchAddress("py",&py_var);
	T6->SetBranchAddress("pz",&pz_var);
	//T6->SetBranchAddress("x",&x_var);
   	//T6->SetBranchAddress("y",&y_var);
	//T6->SetBranchAddress("z",&z_var);
	T6->SetBranchAddress("mother1",&mother1_var);
	T6->SetBranchAddress("mother2",&mother2_var);
	T6->SetBranchAddress("motherid1",&motherid1_var);
	T6->SetBranchAddress("motherid2",&motherid2_var);

	

	vect* v= new vect;//declares a vector in which we are going to store all the particles from the same event
	double prev_index = 0.0;//saving the event index of the previous entry of the tree

	vect1* vmu = new vect1;
	for(int i=0;i<4;i++){
		vmu->S.push_back(0);
		vmu->B.push_back(0);
		vmu->EMB.push_back(0);
		vmu->misID.push_back(0);
	}
	vect2* vgamma = new vect2;
	int column = 16;
	int row = 20;
	for (int i=0;i<column;i++){
		for (int j=0;j<row;j++){
			vgamma->EMB[i][j] = 0;
			vgamma->S1[i][j] = 0;
			vgamma->S2[i][j] = 0;
			vgamma->B[i][j] = 0;
			vgamma->misID[i][j] = 0;
		}
	}
	


	double mu_ptcut = 0.5;
	double mu_pcut = 10.0;
	double eta_cut1 = 2.0;
	double eta_cut2 = 4.5;
	double gamma_ptcut0 = 0.0;
	double gamma_pcut0 = 0.0;
	double gamma_ptcut_step = 0.05;
	double gamma_pcut_step = 0.5;

	double fiducial_min = 6.0;
	double fiducial_max = 22.0;

	double total_KS_events = 100000.0;

	TH2D *ptvsp = new TH2D("ptvsp","S/#sqrt{B_{K_{S}}} based on p_{t} vs p cuts", column, 0.0, gamma_ptcut_step*column, row, 0.0, gamma_pcut_step*row);
	ptvsp -> GetXaxis()-> SetTitle("p_{t} (GeV)");
	ptvsp -> GetYaxis()-> SetTitle("p (GeV)");

	TH2D *S_ptvsp = new TH2D("S_ptvsp","S_{cut}/S_{total} based on p_{t} vs p cuts", column, 0.0, gamma_ptcut_step*column, row, 0.0, gamma_pcut_step*row);
	S_ptvsp -> GetXaxis()-> SetTitle("p_{t} (GeV)");
	S_ptvsp -> GetYaxis()-> SetTitle("p (GeV)");

	TH2D *B_ptvsp = new TH2D("B_ptvsp","B_{K_{S} cut}/B_{K_{S} total} based on p_{t} vs p cuts", column, 0.0, gamma_ptcut_step*column, row, 0.0, gamma_pcut_step*row);
	B_ptvsp -> GetXaxis()-> SetTitle("p_{t} (GeV)");
	B_ptvsp -> GetYaxis()-> SetTitle("p (GeV)");

	TH1D *KS_number_event = new TH1D("KS_number_event","K_{S} number per event", total_KS_events, 0.0, total_KS_events);
    	KS_number_event -> GetXaxis()-> SetTitle("event index");
	KS_number_event -> GetYaxis()-> SetTitle("number of K_{S}");

	double epsilon2 = 1e-9;//kinetic mixing coupling
	double massA = 0.35;//dark photon mass
	double dmassA = 0.005;//half the mass window around the A' mass

	//EM background: mu and gamma

	//Loop through the entries of the tree
	Long64_t nentries = T6->GetEntries();
std::cout<<"Total number of entries: "<<nentries<<std::endl;
   	for (Long64_t i=0;i<nentries;i++){

	      	T6->GetEntry(i);

		if(index_var>=1000000){break;}

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
			v->x.push_back(x_var);
			v->y.push_back(y_var);
			v->z.push_back(z_var);
			v->mother1.push_back(mother1_var);
			v->mother2.push_back(mother2_var);
			v->motherid1.push_back(motherid1_var);
			v->motherid2.push_back(motherid2_var);
			
			analyze_event1(file_out1, v, vmu, vgamma, column, row, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut0, gamma_pcut0, gamma_ptcut_step, gamma_pcut_step, massA, dmassA, fiducial_min, fiducial_max);

			//clear vector
			v->index.clear();
			v->id.clear();
			v->energy.clear();
			v->mass.clear();
			v->px.clear();
			v->py.clear();
			v->pz.clear();
			v->x.clear();
			v->y.clear();
			v->z.clear();
			v->mother1.clear();
			v->mother2.clear();
			v->motherid1.clear();
			v->motherid2.clear();

		}
		//checks if the new entry is from the same event as the previous one
		else if(prev_index!=index_var){


			analyze_event1(file_out1, v, vmu, vgamma, column, row, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut0, gamma_pcut0, gamma_ptcut_step, gamma_pcut_step, massA, dmassA, fiducial_min, fiducial_max);

			//clear vector
			v->index.clear();
			v->id.clear();
			v->energy.clear();
			v->mass.clear();
			v->px.clear();
			v->py.clear();
			v->pz.clear();
			v->x.clear();
			v->y.clear();
			v->z.clear();
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
			v->x.push_back(x_var);
			v->y.push_back(y_var);
			v->z.push_back(z_var);
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
			v->x.push_back(x_var);
			v->y.push_back(y_var);
			v->z.push_back(z_var);
			v->mother1.push_back(mother1_var);
			v->mother2.push_back(mother2_var);
			v->motherid1.push_back(motherid1_var);
			v->motherid2.push_back(motherid2_var);

		}

   	}

	prev_index = 0;

	//A' signal: mu and gamma

	//Loop through the entries of the tree
	nentries = T5->GetEntries();
std::cout<<"Total number of entries: "<<nentries<<std::endl;
   	for (Long64_t i=0;i<nentries;i++){

	      	T5->GetEntry(i);

		if(index_var>=1000000){break;}

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
			v->x.push_back(x_var);
			v->y.push_back(y_var);
			v->z.push_back(z_var);
			v->mother1.push_back(mother1_var);
			v->mother2.push_back(mother2_var);
			v->motherid1.push_back(motherid1_var);
			v->motherid2.push_back(motherid2_var);

			
			analyze_event2(file_out1, v, vmu, vgamma, column, row, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut0, gamma_pcut0, gamma_ptcut_step, gamma_pcut_step, massA, dmassA, fiducial_min, fiducial_max);

			//clear vector
			v->index.clear();
			v->id.clear();
			v->energy.clear();
			v->mass.clear();
			v->px.clear();
			v->py.clear();
			v->pz.clear();
			v->x.clear();
			v->y.clear();
			v->z.clear();
			v->mother1.clear();
			v->mother2.clear();
			v->motherid1.clear();
			v->motherid2.clear();

		}
		//checks if the new entry is from the same event as the previous one
		else if(prev_index!=index_var){


			analyze_event2(file_out1, v, vmu, vgamma, column, row, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut0, gamma_pcut0, gamma_ptcut_step, gamma_pcut_step, massA, dmassA, fiducial_min, fiducial_max);

			//clear vector
			v->index.clear();
			v->id.clear();
			v->energy.clear();
			v->mass.clear();
			v->px.clear();
			v->py.clear();
			v->pz.clear();
			v->x.clear();
			v->y.clear();
			v->z.clear();
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
			v->x.push_back(x_var);
			v->y.push_back(y_var);
			v->z.push_back(z_var);
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
			v->x.push_back(x_var);
			v->y.push_back(y_var);
			v->z.push_back(z_var);
			v->mother1.push_back(mother1_var);
			v->mother2.push_back(mother2_var);
			v->motherid1.push_back(motherid1_var);
			v->motherid2.push_back(motherid2_var);

		}

   	}

	prev_index = 0.0;

	//K_S tree
	vect* v1= new vect;//declares a vector in which we are going to store all the gamma from the same event
	Long64_t nentries1 = T2->GetEntries();
	Long64_t index_T2 = 0.0;
	double redecay = 100.0;

	//Loop through the entries of tree T2
	//nentries = T2->GetEntries();
	nentries = T1->GetEntries();
std::cout<<"Total number of entries T1: "<<nentries<<std::endl;
   	for (Long64_t i=0;i<nentries;i++){

	      	//T2->GetEntry(i);
		T1->GetEntry(i);

		if(index_var>=10000){break;}

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
			v->x.push_back(x_var);
			v->y.push_back(y_var);
			v->z.push_back(z_var);
			v->mother1.push_back(mother1_var);
			v->mother2.push_back(mother2_var);
			v->motherid1.push_back(motherid1_var);
			v->motherid2.push_back(motherid2_var);

			for (Long64_t j = index_T2 ;j<nentries1;j++){

				T2->GetEntry(j);

				if(index_var == index_var1){

					if(id_var1 == 22){

						v1->index.push_back(index_var1);
						v1->id.push_back(id_var1);
						v1->energy.push_back(energy_var1);
						v1->mass.push_back(mass_var1);
						v1->px.push_back(px_var1);
						v1->py.push_back(py_var1);
						v1->pz.push_back(pz_var1);
						v1->x.push_back(x_var1);
						v1->y.push_back(y_var1);
						v1->z.push_back(z_var1);
						v1->mother1.push_back(mother1_var1);
						v1->mother2.push_back(mother2_var1);
						v1->motherid1.push_back(motherid1_var1);
						v1->motherid2.push_back(motherid2_var1);

					}

				}
				else if(index_var1>index_var){break;}

			}

			analyze_event_K(file_out1, v, v1, vmu, vgamma, column, row, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut0, gamma_pcut0, gamma_ptcut_step, gamma_pcut_step, massA, dmassA, fiducial_min, fiducial_max, alpha, beta, gamma, redecay, KS_number_event, pythia1);

			//clear vector
			v->index.clear();
			v->id.clear();
			v->energy.clear();
			v->mass.clear();
			v->px.clear();
			v->py.clear();
			v->pz.clear();
			v->x.clear();
			v->y.clear();
			v->z.clear();
			v->mother1.clear();
			v->mother2.clear();
			v->motherid1.clear();
			v->motherid2.clear();

		}

		//checks if the new entry is from the same event as the previous one
		else if(prev_index!=index_var){

			for (Long64_t j = index_T2;j<nentries1;j++){

				T2->GetEntry(j);

				if(index_var == index_var1){

					if(id_var1 == 22){

						v1->index.push_back(index_var1);
						v1->id.push_back(id_var1);
						v1->energy.push_back(energy_var1);
						v1->mass.push_back(mass_var1);
						v1->px.push_back(px_var1);
						v1->py.push_back(py_var1);
						v1->pz.push_back(pz_var1);
						v1->x.push_back(x_var1);
						v1->y.push_back(y_var1);
						v1->z.push_back(z_var1);
						v1->mother1.push_back(mother1_var1);
						v1->mother2.push_back(mother2_var1);
						v1->motherid1.push_back(motherid1_var1);
						v1->motherid2.push_back(motherid2_var1);

					}

				}
				else if(index_var1>index_var){index_T2 = j; break;}

			}

			analyze_event_K(file_out1, v, v1, vmu, vgamma, column, row, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut0, gamma_pcut0, gamma_ptcut_step, gamma_pcut_step, massA, dmassA, fiducial_min, fiducial_max, alpha, beta, gamma, redecay, KS_number_event, pythia1);

			//clear vector
			v->index.clear();
			v->id.clear();
			v->energy.clear();
			v->mass.clear();
			v->px.clear();
			v->py.clear();
			v->pz.clear();
			v->x.clear();
			v->y.clear();
			v->z.clear();
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
			v->x.push_back(x_var);
			v->y.push_back(y_var);
			v->z.push_back(z_var);
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
			v->x.push_back(x_var);
			v->y.push_back(y_var);
			v->z.push_back(z_var);
			v->mother1.push_back(mother1_var);
			v->mother2.push_back(mother2_var);
			v->motherid1.push_back(motherid1_var);
			v->motherid2.push_back(motherid2_var);

		}
	      	
   	}

	double sigma, gamma_ptcut, gamma_pcut;
	double width = Gamma_A(epsilon2, massA);//calculating the width of the dark photon

	double F = scaling(epsilon2, massA, width); 
	double Br = 0.0003100;
	double RKS = 0.6922420;

	double S, EMB, B, misID;

	for(int i=0;i<column;i++){ 

		gamma_ptcut = gamma_ptcut0 + i * gamma_ptcut_step;

		for(int j=0;j<row;j++){ 

			gamma_pcut = gamma_pcut0 + j * gamma_pcut_step;

			// this is the displaced case
                        // remove DP signal from B
			S = ((double)vgamma->EMB[i][j] + (double)vgamma->B[i][j]) * Br * F * ((double)vgamma->S2[i][j] /  (double)vgamma->S1[i][j]);
			/*
			if (S < 0) {
			  std::cout << "=============" << endl;
			  std::cout << (double)vgamma->EMB[i][j] << std::endl;
			  std::cout << (double)vgamma->B[i][j] << std::endl;
			  std::cout << Br << std::endl;
			  std::cout << F << std::endl;
			  std::cout << (double)vgamma->S2[i][j] << std::endl;
			  std::cout << (double)vgamma->S1[i][j] << std::endl;
			}
			*/
			B = (double)vgamma->misID[i][j]* RKS ;// / redecay;

			sigma = S / sqrt(B);
	
			//S = (double)vgamma->S[i][j] * F * Br  * (double)vmu->EMB[0] /(double)vmu->S[0];
			//B = (double)vgamma->B[i][j] * Br;
			//misID = (double)vgamma->misID[i][j]* RKS;

			//sigma = metric(S, misID);


			ptvsp->SetBinContent(i+1, j+1, sigma);
			S_ptvsp->SetBinContent(i+1, j+1, (double)vgamma->S2[i][j]/(double)vmu->S[0]);
			B_ptvsp->SetBinContent(i+1, j+1, (double)vgamma->B[i][j]/(double)vmu->B[0]);

			//fprintf(file_out1, "%.5e   %.5e   %.5e   %.5e   %.5e   %.5e  \n", gamma_ptcut, gamma_pcut, S, B, misID, sigma);
			fprintf(file_out1, "%.5e   %.5e   %.5e   %.5e   %.5e  \n", gamma_ptcut, gamma_pcut, S, B, sigma);
	
		}

	}
 
	TCanvas *c1=new TCanvas("c1","",600,600);

	c1->SetRightMargin(0.13);

	gPad->SetLogz();

	ptvsp->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project84_ptvsp.pdf","pdf");

	S_ptvsp->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project84_S_ptvsp.pdf","pdf");

	B_ptvsp->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project84_B_ptvsp.pdf","pdf");	

	KS_number_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project84_KS_number_event.pdf","pdf");

	std::cout<<"Signal entries (no cuts): "<<vmu->S[0]<<std::endl;
	std::cout<<"Background entries (no cuts): "<<vmu->B[0]<<std::endl;
	std::cout<<"MisID pion entries (mu cuts + mass cuts): "<<vmu->misID[0]<<std::endl;

	fprintf(file_out1, "Candidate type    Total entries    After A' mass cut   After transverse flight distance cut     After eta mass cut     Signal significance after mu + A' _fiducial cut: \n");
	fprintf(file_out1, "S    			   %d    			%d   			%d 				%d				%e\n", vmu->S[0], vmu->S[1], vmu->S[2], vmu->S[3], (double)vmu->S[1]/(double)vmu->S[0]);
	fprintf(file_out1, "EMB    			   %d    			%d    			%d 				%d				%e\n", vmu->EMB[0], vmu->EMB[1], vmu->EMB[2], vmu->EMB[3], (double)vmu->EMB[1]/(double)vmu->EMB[0]);
	fprintf(file_out1, "B    			   %d    			%d    			%d 				%d				%e\n", vmu->B[0], vmu->B[1], vmu->B[2], vmu->B[3], (double)vmu->B[1]/(double)vmu->B[0]);
	fprintf(file_out1, "misID    		   %d    			%d    			%d 				%d				%e\n", vmu->misID[0], vmu->misID[1], vmu->misID[2], vmu->misID[3], (double)vmu->misID[1]/(double)vmu->misID[0]);

	//S = (double)vmu->S[2] * F * Br * (double)vmu->EMB[0] /(double)vmu->S[0];
	//misID = (double)vmu->misID[2] * RKS;

	S = ((double)vmu->EMB[1] + (double)vmu->B[1]) * Br * F * (double)vmu->S[2]/(double)vmu->S[1];
	B = (double)vmu->misID[2]* RKS ;// / redecay;
	fprintf(file_out1, "Figure of merit after mu + A' + fiducial cut:  %e \n", S / sqrt(B));
	

	file_out->Write();
	file_out->Close();	
        file_in1->Close();
	file_in2->Close();
	delete file_out;

	fclose(file_out1);
	delete v;
	delete v1;
	delete vmu;
	delete vgamma;
	// Done.
	return 0;
}


