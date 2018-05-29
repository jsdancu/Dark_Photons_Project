//05/01/2018
//Program reading in the three trees generated in gs_project34.cc and performing analysis on the reconstructed invariant mass of eta

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
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

struct vect2{
	
	std::vector<double> n;
	
};

double invmass(vect *v, Long64_t i, Long64_t j, Long64_t k){

	//std::cout<<"i = "<<i<<"    "<<"j = "<<j<<"    "<<"k = "<<k<<"    "<< "energy i = "<<v->energy[i]<<"    "<<"energy j = "<<v->energy[j]<<"energy k = "<<v->energy[k]<<std::endl;

	double inv_mass=sqrt(pow(v->energy[i]+v->energy[j]+v->energy[k],2)-pow(v->px[i]+v->px[j]+v->px[k],2)-pow(v->py[i]+v->py[j]+v->py[k],2)-pow(v->pz[i]+v->pz[j]+v->pz[k],2));

	//std::cout<<"invariant mass = "<<setprecision(15)<<inv_mass<<std::endl;

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


//vect2 analyze_event(FILE *file_out1, vect *v, vect2 *vn, double mu_ptcut, double mu_pcut, double eta_cut1, double eta_cut2, double gamma_ptcut, double gamma_pcut, TH1D *eta_invmass, TH1D *muon_invmass, TH1D *muongamma_invmass_mother, TH1D *muon_invmass_mother, TH1D *eta_invmass_mothers, TH1D *muon_invmass_cuts_mother, TH1D *muon_invmass_cuts_mother2, TH1D *eta_invmass_cuts_mother, TH1D *eta_invmass_cuts_mother2, TH1D *muon_invmass_c3_mother, TH1D *eta_invmass_c3_mother, TH1D *muon_invmass_cuts_bkg, TH1D *muon_invmass_cuts_bkg2, TH1D *eta_invmass_cuts_bkg, TH1D *eta_invmass_cuts_bkg2, TH1D *muon_invmass_c3_bkg, TH1D *eta_invmass_c3_bkg, TH1D *muon_pt_mother, TH1D *muon_pt_bkg, TH1D *muon_ptcut_mother, TH1D *muon_ptetacut_mother, TH1D *muon_ptetamasscut_mother, TH1D *muon_p_mother, TH1D *muon_p_bkg, TH1D *muon_pcut_mother, TH1D *muon_petacut_mother, TH1D *muon_petamasscut_mother, TH2D *muon_ptvsp_mother, TH2D *muon_ptvsp_cuts_mother, TH2D *muon_ptvsp_c1_mother, TH2D *muon_ptvsp_c2_mother, TH2D *muon_ptvsp_c3_mother, TH2D *gamma_ptvsp_mother, TH2D *gamma_ptvsp_cuts_mother, TH2D *gamma_ptvsp_c1_mother, TH2D *gamma_ptvsp_c2_mother, TH2D *gamma_ptvsp_c3_mother, TH2D *muon_ptvseta_mother, TH2D *muon_ptvseta_cuts_mother, TH2D *muon_ptvseta_c1_mother, TH2D *muon_ptvseta_c2_mother, TH2D *muon_ptvseta_c3_mother, TH2D *gamma_ptvseta_mother, TH2D *gamma_ptvseta_cuts_mother, TH2D *gamma_ptvseta_c1_mother, TH2D *gamma_ptvseta_c2_mother, TH2D *gamma_ptvseta_c3_mother, TH1D *muon_eta_mother, TH1D *gamma_pt_mother, TH1D *gamma_pt_bkg, TH1D *gamma_pt_c3_bkg, TH1D *gamma_pt_c4_bkg, TH1D *gamma_pt_c5_bkg, TH1D *gamma_ptcut_mother, TH1D *gamma_ptetacut_mother, TH1D *gamma_ptetamasscut_mother, TH1D *gamma_ptallcut_mother, TH1D *gamma_p_mother, TH1D *gamma_p_bkg, TH1D *gamma_p_c3_bkg, TH1D *gamma_p_c4_bkg, TH1D *gamma_p_c5_bkg, TH1D *gamma_pcut_mother, TH1D *gamma_petacut_mother, TH1D *gamma_petamasscut_mother, TH1D *gamma_pallcut_mother, TH1D *gamma_eta_mother, TH1D *mu_number_event, TH1D *gamma_number_event, TH1D *total_number_event, TH1D *mugamma_combnumber_event, TH1D *gamma_fake){

void analyze_event(FILE *file_out1, vect *v, vect2 *vn, double mu_ptcut, double mu_pcut, double eta_cut1, double eta_cut2, double gamma_ptcut, double gamma_pcut, TH1D *eta_invmass, TH1D *muon_invmass, TH1D *muongamma_invmass_mother, TH1D *muon_invmass_mother, TH1D *eta_invmass_mothers, TH1D *muon_invmass_cuts_mother, TH1D *muon_invmass_cuts_mother2, TH1D *eta_invmass_cuts_mother, TH1D *eta_invmass_cuts_mother2, TH1D *muon_invmass_c3_mother, TH1D *eta_invmass_c3_mother, TH1D *muon_invmass_cuts_sigbkg, TH1D *muon_invmass_cuts_sigbkg2, TH1D *eta_invmass_cuts_sigbkg, TH1D *eta_invmass_cuts_sigbkg2, TH1D *muon_invmass_c3_sigbkg, TH1D *eta_invmass_c3_sigbkg, TH1D *muon_invmass_cuts_bkg, TH1D *muon_invmass_cuts_bkg2, TH1D *eta_invmass_cuts_bkg, TH1D *eta_invmass_cuts_bkg2, TH1D *muon_invmass_c3_bkg, TH1D *eta_invmass_c3_bkg, TH1D *muon_pt_mother, TH1D *muon_pt_bkg, TH1D *muon_ptcut_mother, TH1D *muon_ptetacut_mother, TH1D *muon_ptetamasscut_mother, TH1D *muon_p_mother, TH1D *muon_p_bkg, TH1D *muon_pcut_mother, TH1D *muon_petacut_mother, TH1D *muon_petamasscut_mother, TH1D *muon_eta_mother, TH1D *gamma_pt_mother, TH1D *gamma_pt_bkg, TH1D *gamma_pt_c3_bkg, TH1D *gamma_pt_c4_bkg, TH1D *gamma_pt_c5_bkg, TH1D *gamma_ptcut_mother, TH1D *gamma_ptetacut_mother, TH1D *gamma_ptetamasscut_mother, TH1D *gamma_ptallcut_mother, TH1D *gamma_p_mother, TH1D *gamma_p_bkg, TH1D *gamma_p_c3_bkg, TH1D *gamma_p_c4_bkg, TH1D *gamma_p_c5_bkg, TH1D *gamma_pcut_mother, TH1D *gamma_petacut_mother, TH1D *gamma_petamasscut_mother, TH1D *gamma_pallcut_mother, TH1D *gamma_eta_mother, TH1D *mu_number_event, TH1D *gamma_number_event, TH1D *total_number_event, TH1D *mugamma_combnumber_event, TH1D *gamma_fake){

	double inv_mass, mu_invmass, pt1, pt2, pt3, p1, p2, p3, eta1, eta2, eta3;

	Long64_t nentries = v->index.size();

	Long64_t mu_number = 0;
	Long64_t antimu_number = 0;
	Long64_t gamma_number = 0;
	Long64_t mugamma_combnumber = 0;

	vect2* N = new vect2;
	N->n = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	double m_eta = 0.54785;
	double dm_eta = 0.02;

	bool antimu = false;
	bool gamma = false;

	for(Long64_t i=0;i<nentries;i++){ 

		//if entry is muon
		if(v->id[i] == 13){

			++mu_number;

			//Loop through the event and find an anti-muon
			for(Long64_t j=0;j<nentries;j++){

				//if entry is anti-muon
				if(v->id[j] == -13){

					if(antimu == false){
						++antimu_number;
					}

					for(Long64_t k=0;k<nentries;k++){

						//if entry is photon
						if(v->id[k] == 22){

							++mugamma_combnumber;

							if(gamma == false){
								++gamma_number;
							}

							//Reconstructing invariant mass of all muon-antimuon-photon combinations and di-muons
							inv_mass = invmass(v, i, j, k);
							mu_invmass = invmass2(v, i, j);

							//Plot histogram of all muon-antimuon-photon combinations and di-muons
							eta_invmass->Fill(inv_mass);
							muon_invmass->Fill(mu_invmass);

							pt1 = pt(v, i);
							pt2 = pt(v, j);
							pt3 = pt(v, k);
							p1 = p(v, i);
							p2 = p(v, j);
							p3 = p(v, k);
							eta1 = eta(v, i);
							eta2 = eta(v, j);
							eta3 = eta(v, k);	


							if((pt1>mu_ptcut) && (pt2>mu_ptcut) && (p1>mu_pcut) && (p2>mu_pcut) && (eta1>eta_cut1) && (eta1<eta_cut2) && (eta2>eta_cut1) && (eta2<eta_cut2)){

								eta_invmass_c3_sigbkg->Fill(inv_mass);
								muon_invmass_c3_sigbkg->Fill(mu_invmass);

								/*if(!((v->mother1[i] == v->mother1[j]) && (v->mother1[i] == v->mother1[k]) && (v->mother2[i] == v->mother2[j]) && (v->mother2[i] == v->mother2[k]) && (v->motherid1[i] == 221))){
									eta_invmass_c3_bkg->Fill(inv_mass);
									muon_invmass_c3_bkg->Fill(mu_invmass);
								}*/

								if((inv_mass < m_eta+dm_eta) && (inv_mass> m_eta-dm_eta)){

									eta_invmass_cuts_sigbkg->Fill(inv_mass);
									muon_invmass_cuts_sigbkg->Fill(mu_invmass);
		
									if(!((v->mother1[i] == v->mother1[j]) && (v->mother1[i] == v->mother1[k]) && (v->mother2[i] == v->mother2[j]) && (v->mother2[i] == v->mother2[k]) && (v->motherid1[i] == 221))){
										eta_invmass_cuts_bkg->Fill(inv_mass);
										muon_invmass_cuts_bkg->Fill(mu_invmass);
									}

									if((pt3>gamma_ptcut) && (p3>gamma_pcut) && (eta3>eta_cut1) && (eta3<eta_cut2)){

										eta_invmass_cuts_sigbkg2->Fill(inv_mass);
										muon_invmass_cuts_sigbkg2->Fill(mu_invmass);

										if(!((v->mother1[i] == v->mother1[j]) && (v->mother1[i] == v->mother1[k]) && (v->mother2[i] == v->mother2[j]) && (v->mother2[i] == v->mother2[k]) && (v->motherid1[i] == 221))){
											eta_invmass_cuts_bkg2->Fill(inv_mass);
											muon_invmass_cuts_bkg2->Fill(mu_invmass);
										}

										if((v->motherid1[k] != 221) || ((v->motherid1[k] == 221) && (v->mother1[i] != v->mother1[k]) && (v->mother1[j] != v->mother1[k]))){

											gamma_fake->Fill(v->motherid1[k]);
										

										}

									}

								}

							}

					

							if((v->mother1[i] == v->mother1[j]) && (v->mother1[i] == v->mother1[k]) && (v->mother2[i] == v->mother2[j]) && (v->mother2[i] == v->mother2[k])){

								//If the muon-antimuon pair come from the same mother particle plot them in a separate histogram
								muongamma_invmass_mother->Fill(inv_mass);

								if(v->motherid1[i] == 221){

									muon_pt_mother->Fill(pt1);
									muon_pt_mother->Fill(pt2);
									gamma_pt_mother->Fill(pt3);
									muon_p_mother->Fill(p1);
									muon_p_mother->Fill(p2);
									gamma_p_mother->Fill(p3);
									muon_eta_mother->Fill(eta1);
									muon_eta_mother->Fill(eta2);
									gamma_eta_mother->Fill(eta3);

									++N->n[0];

									muon_invmass_mother->Fill(mu_invmass);
									eta_invmass_mothers->Fill(inv_mass);

									//Plot pt vs p 2D histogram before cuts
									/*muon_ptvsp_mother->Fill(p1, pt1);
									muon_ptvsp_mother->Fill(p2, pt2);
									gamma_ptvsp_mother->Fill(p3, pt3);
									muon_ptvseta_mother->Fill(pt1, eta1);
									muon_ptvseta_mother->Fill(pt2, eta2);
									gamma_ptvseta_mother->Fill(pt3, eta3);*/

									if((pt1>mu_ptcut) && (pt2>mu_ptcut)){

										++N->n[1];

										/*muon_ptvsp_c1_mother->Fill(p1, pt1);
										muon_ptvsp_c1_mother->Fill(p2, pt2);
										gamma_ptvsp_c1_mother->Fill(p3, pt3);
										muon_ptvseta_c1_mother->Fill(pt1, eta1);
										muon_ptvseta_c1_mother->Fill(pt2, eta2);
										gamma_ptvseta_c1_mother->Fill(pt3, eta3);*/

										if((p1>mu_pcut) && (p2>mu_pcut)){

											++N->n[2];

											muon_ptcut_mother->Fill(pt1);
											muon_ptcut_mother->Fill(pt2);
											gamma_ptcut_mother->Fill(pt3);
											muon_pcut_mother->Fill(p1);
											muon_pcut_mother->Fill(p2);
											gamma_pcut_mother->Fill(p3);

											/*muon_ptvsp_c2_mother->Fill(p1, pt1);
											muon_ptvsp_c2_mother->Fill(p2, pt2);
											gamma_ptvsp_c2_mother->Fill(p3, pt3);
											muon_ptvseta_c2_mother->Fill(pt1, eta1);
											muon_ptvseta_c2_mother->Fill(pt2, eta2);
											gamma_ptvseta_c2_mother->Fill(pt3, eta3);*/

											if((eta1>eta_cut1) && (eta1<eta_cut2) && (eta2>eta_cut1) && (eta2<eta_cut2))
											{

												++N->n[3];

												muon_ptetacut_mother->Fill(pt1);
												muon_ptetacut_mother->Fill(pt2);
												gamma_ptetacut_mother->Fill(pt3);
												muon_petacut_mother->Fill(p1);
												muon_petacut_mother->Fill(p2);
												gamma_petacut_mother->Fill(p3);

												/*muon_ptvsp_c3_mother->Fill(p1, pt1);
												muon_ptvsp_c3_mother->Fill(p2, pt2);
												gamma_ptvsp_c3_mother->Fill(p3, pt3);
												muon_ptvseta_c3_mother->Fill(pt1, eta1);
												muon_ptvseta_c3_mother->Fill(pt2, eta2);
												gamma_ptvseta_c3_mother->Fill(pt3, eta3);*/

												muon_invmass_c3_mother->Fill(mu_invmass);
												eta_invmass_c3_mother->Fill(inv_mass);

												if((inv_mass < m_eta+dm_eta) && (inv_mass> m_eta-dm_eta)){

													++N->n[4];

													muon_ptetamasscut_mother->Fill(pt1);
													muon_ptetamasscut_mother->Fill(pt2);
													gamma_ptetamasscut_mother->Fill(pt3);
													muon_petamasscut_mother->Fill(p1);
													muon_petamasscut_mother->Fill(p2);
													gamma_petamasscut_mother->Fill(p3);

													//Plot pt vs p 2D histogram after all cuts
													/*muon_ptvsp_cuts_mother->Fill(p1, pt1);
													muon_ptvsp_cuts_mother->Fill(p2, pt2);
													gamma_ptvsp_cuts_mother->Fill(p3, pt3);
													muon_ptvseta_cuts_mother->Fill(pt1, eta1);
													muon_ptvseta_cuts_mother->Fill(pt2, eta2);
													gamma_ptvseta_cuts_mother->Fill(pt3, eta3);*/

													muon_invmass_cuts_mother->Fill(mu_invmass);
													eta_invmass_cuts_mother->Fill(inv_mass);

												
													if((pt3>gamma_ptcut) && (p3>gamma_pcut) && (eta3>eta_cut1) && (eta3<eta_cut2)){

														++N->n[5];

														eta_invmass_cuts_mother2->Fill(inv_mass);
														muon_invmass_cuts_mother2->Fill(mu_invmass);

														gamma_ptallcut_mother->Fill(pt3);
														gamma_pallcut_mother->Fill(p3);

													}

													
												}

											}

										}

									}

								}

							}
							else if(v->mother1[i] != v->mother1[j]){

								muon_pt_bkg->Fill(pt1);

								muon_pt_bkg->Fill(pt2);

								muon_p_bkg->Fill(p1);

								muon_p_bkg->Fill(p2);

							}
							//else if((v->mother1[i] != v->mother1[j]) && (v->mother1[i] != v->mother1[k]) && (v->mother1[k] != v->mother1[j]) && (v->motherid1[k] != 221)){
							else if((v->motherid1[k] != 221) || ((v->motherid1[k] == 221) && (v->mother1[i] != v->mother1[k]) && (v->mother1[j] != v->mother1[k]))){

								gamma_pt_bkg->Fill(pt3);
								gamma_p_bkg->Fill(p3);

								if((pt1>mu_ptcut) && (pt2>mu_ptcut) && (p1>mu_pcut) && (p2>mu_pcut) && (eta1>eta_cut1) && (eta1<eta_cut2) && (eta2>eta_cut1) && (eta2<eta_cut2)){

									gamma_pt_c3_bkg->Fill(pt3);
									gamma_p_c3_bkg->Fill(p3);

									if((inv_mass < m_eta+dm_eta) && (inv_mass> m_eta-dm_eta)){

										gamma_pt_c4_bkg->Fill(pt3);
										gamma_p_c4_bkg->Fill(p3);

										if((pt3>gamma_ptcut) && (p3>gamma_pcut) && (eta3>eta_cut1) && (eta3<eta_cut2)){

											gamma_pt_c5_bkg->Fill(pt3);
											gamma_p_c5_bkg->Fill(p3);


										}

									}

								}

							}

						}


					}
					gamma = true;

				}

			}
			antimu = true;
		}

	}

	//std::cout<<"muon number = "<<mu_number<<setw(5)<<"anti muon = "<<antimu_number<<std::endl;
	double total_mu_antimu = mu_number+antimu_number;
	double total_number = mu_number+antimu_number+gamma_number;
	if(v->index.size() > 0) {

		mu_number_event->SetBinContent((v->index[1])+1,total_mu_antimu);
		gamma_number_event->SetBinContent((v->index[1])+1,gamma_number);
		total_number_event->SetBinContent((v->index[1])+1,total_number);
		mugamma_combnumber_event->SetBinContent((v->index[1])+1,mugamma_combnumber);

	}
	
	//if (v->index.size() > 0) std::cout<<v->index[0]<<std::endl;
	//std::cout<<v << "    " << v->index.size() <<std::endl;

	if(N->n[0]!=0.0){fprintf(file_out1, "%10f %10f %10f %10f %10f \n", N->n[1]/N->n[0], N->n[2]/N->n[0], N->n[3]/N->n[0], N->n[4]/N->n[0], N->n[5]/N->n[0]);}
	else{fprintf(file_out1, "%10f %10f %10f %10f %10f \n", N->n[1], N->n[2], N->n[3], N->n[4], N->n[5]);}

	for(Long64_t i=0;i<6;i++){ 
	
		vn->n[i] = vn->n[i] + N->n[i];
	
	}

	delete N;

}

//void analyze_event_pi(FILE *file_out2, vect *v, vect2 *vp, double mu_ptcut, double mu_pcut, double eta_cut1, double eta_cut2, double gamma_ptcut, double gamma_pcut, TH1D *pigamma_invmass, TH1D *pigamma_invmass_c3, TH1D *pigamma_invmass_cuts, TH1D *pigamma_invmass_cuts2, TH1D *pi_invmass, TH1D *pi_invmass_c3, TH1D *pi_invmass_cuts, TH1D *pi_invmass_cuts2, TH1D *pi_pt, TH1D *pi_p, TH1D *pi_eta, TH1D *pi_ptcut, TH1D *pi_pcut, TH1D *pi_ptmasscut, TH1D *pi_pmasscut, TH2D *pi_ptvsp, TH2D *pi_ptvsp_cuts, TH1D *pi_number_event, TH1D *pi_combnumber_event, TH1D *total_pigamma_number_event, TH1D *pigamma_combnumber_event){

void analyze_event_pi(FILE *file_out2, vect *v, vect2 *vp, double alpha, double beta, double gamma, double mu_ptcut, double mu_pcut, double eta_cut1, double eta_cut2, double gamma_ptcut, double gamma_pcut, TH1D *pigamma_invmass, TH1D *pigamma_invmass_c3, TH1D *pigamma_invmass_cuts, TH1D *pigamma_invmass_cuts2, TH1D *pi_invmass, TH1D *pi_invmass_c3, TH1D *pi_invmass_cuts, TH1D *pi_invmass_cuts2, TH1D *pi_pt, TH1D *pi_p, TH1D *pi_eta, TH1D *pi_ptcut, TH1D *pi_pcut, TH1D *pi_ptmasscut, TH1D *pi_pmasscut, TH1D *pi_number_event, TH1D *pi_combnumber_event, TH1D *total_pigamma_number_event, TH1D *pigamma_combnumber_event){

	double inv_mass, pion_invmass, pt_pp, pt_pn, p_pp, p_pn, eta_pp, eta_pn, pt3, p3, eta3;
	double misid1, misid2;

	Long64_t nentries = v->index.size();

	Long64_t pi_number_p = 0;
	Long64_t pi_number_n = 0;
	Long64_t gamma_number = 0;
	Long64_t pigamma_combnumber = 0;
	Long64_t pi_combnumber = 0;

	vect2* N = new vect2;
	N->n = {0.0, 0.0, 0.0, 0.0};

	double m_eta = 0.54785;
	double dm_eta = 0.02;

	bool pi_n = false;
	bool gamma1 = false;

	for(Long64_t i=0;i<nentries;i++){ 

		//if entry is pi+
		if(v->id[i] == 211){

			++pi_number_p;

			pt_pp = pt(v, i);
			pi_pt->Fill(pt_pp);
			p_pp = p(v, i);
			pi_p->Fill(p_pp);
			eta_pp = eta(v, i);
			pi_eta->Fill(eta_pp);

			if((pt_pp>mu_ptcut) && (p_pp>mu_pcut) && (eta_pp>eta_cut1) && (eta_pp<eta_cut2)){

				//Loop through the event and find an anti-muon
				for(Long64_t j=0;j<nentries;j++){

					//if entry is anti-muon
					if(v->id[j] == -211){

						++pi_combnumber;

						if(pi_n == false){
							++pi_number_n;
						}

						pt_pn = pt(v, j);
						pi_pt->Fill(pt_pn);
						p_pn = p(v, j);
						pi_p->Fill(p_pn);
						eta_pn = eta(v, j);
						pi_eta->Fill(eta_pn);

						if((pt_pn>mu_ptcut) && (p_pn>mu_pcut) && (eta_pn>eta_cut1) && (eta_pn<eta_cut2)){

							for(Long64_t k=0;k<nentries;k++){

								//if entry is photon
								if(v->id[k] == 22){

									++pigamma_combnumber;

									if(gamma1 == false){
										++gamma_number;
									}

									//check that pion pair is not a K^0_S (K short) resonance 
									if((v->motherid1[i] != 310) && (v->motherid1[j] != 310)){

										//Calculating the misID rate of the "muon" candidates
										misid1 = misID_rate(abs(v->pz[i]), alpha, beta, gamma);
										misid2 = misID_rate(abs(v->pz[j]), alpha, beta, gamma);

										//Reconstructing invariant mass of all pi+/- combinations
										inv_mass = invmass(v, i, j, k);
										pion_invmass = invmass2(v, i, j);

										//Plot histogram of all pi+/- combinations
										pigamma_invmass->Fill(inv_mass, misid1 * misid2);
										pi_invmass->Fill(pion_invmass, misid1 * misid2);	

										pt3 = pt(v, k);
										p3 = p(v, k);
										eta3 = eta(v, k);

										//Plot pt vs p 2D histogram before cuts
										/*pi_ptvsp->Fill(p_pp, pt_pp);
										pi_ptvsp->Fill(p_pn, pt_pn);*/

										++N->n[1];

										pi_ptcut->Fill(pt_pp);
										pi_ptcut->Fill(pt_pn);
										pi_pcut->Fill(p_pp);
										pi_pcut->Fill(p_pn);

										pigamma_invmass_c3->Fill(inv_mass, misid1 * misid2);
										pi_invmass_c3->Fill(pion_invmass, misid1 * misid2);

										if((inv_mass < m_eta+dm_eta) && (inv_mass > m_eta-dm_eta)){

											++N->n[2];

											pi_ptmasscut->Fill(pt_pp);
											pi_ptmasscut->Fill(pt_pn);
											pi_pmasscut->Fill(p_pp);
											pi_pmasscut->Fill(p_pn);

											//Plot pt vs p 2D histogram after cuts
											/*pi_ptvsp_cuts->Fill(p_pp, pt_pp);
											pi_ptvsp_cuts->Fill(p_pn, pt_pn);*/

											pigamma_invmass_cuts->Fill(inv_mass, misid1 * misid2);
											pi_invmass_cuts->Fill(pion_invmass, misid1 * misid2);

											if((pt3>gamma_ptcut) && (p3>gamma_pcut) && (eta3>eta_cut1) && (eta3<eta_cut2)){

												++N->n[3];

												pigamma_invmass_cuts2->Fill(inv_mass, misid1 * misid2);
												pi_invmass_cuts2->Fill(pion_invmass, misid1 * misid2);

											}


										}

									}

								}

							}
							gamma1=true;

						}

					}

				}
				pi_n = true;

			}

		}

	}

	double total_pi_pn = pi_number_p+pi_number_n;
	double total_number = pi_number_p+pi_number_n+gamma_number;
	if (v->index.size() > 0) {

		pi_number_event->SetBinContent((v->index[1])+1,total_pi_pn);
		pi_combnumber_event->SetBinContent((v->index[1])+1, pi_combnumber);
		total_pigamma_number_event->SetBinContent((v->index[1])+1,total_number);
		pigamma_combnumber_event->SetBinContent((v->index[1])+1,pigamma_combnumber);

	}
	
	N->n[0] = pigamma_combnumber;

	if(N->n[0]!=0.0){fprintf(file_out2, "%10f %10f %10f \n", N->n[1]/N->n[0], N->n[2]/N->n[0], N->n[3]/N->n[0]);}
	else{fprintf(file_out2, "%10f %10f %10f \n", N->n[1], N->n[2], N->n[3]);}

	for(Long64_t i=0;i<6;i++){ 
	
		vp->n[i] = vp->n[i] + N->n[i];
	
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
	TFile  *file_in = new TFile("gs_project34.root", "READ");

	// Open the output TFile 
	TFile  *file_out = new TFile("gs_project54.root", "recreate");

	FILE  *file_out1 = fopen("mugamma_out54.txt", "w");
	fprintf(file_out1, "Table showing signal significance of signal muons passing various cuts: \n");
	fprintf(file_out1, "pt     p    eta    mass \n");

	FILE  *file_out2 = fopen("pi_out54.txt", "w");
	fprintf(file_out2, "Table showing signal significance of pions passing various cuts: \n");
	fprintf(file_out2, "pt    p     eta     mass \n");

	double alpha, beta, gamma;

	std::ifstream file("results_out94.txt");
    	if(file.is_open()){

		file >> alpha >> beta >> gamma;

	}


	//Get trees from file
	TTree *T1 = (TTree*)file_in->Get("T1");
	TTree *T2 = (TTree*)file_in->Get("T2");
	TTree *T3 = (TTree*)file_in->Get("T3");

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

	Int_t Nbins1 = 250;
	Double_t Edges1[Nbins1] = {0.2};

	for(Long64_t i=1;i<Nbins1;i++){

		Edges1[i] = Edges1[i-1] * (1.0 + 1.6/100.0);

	}

	Int_t Nbins2 = 100;
	Double_t Edges2[Nbins2] = {0.2};

	for(Long64_t i=1;i<Nbins2;i++){

		Edges2[i] = Edges2[i-1] * (1.0 + 1.6/100.0);

	}



	TH1D *eta_invmass = new TH1D("eta_invmass","Reconstructed eta invariant mass from muon pairs and #gamma", Nbins1 - 1, Edges1);
    	eta_invmass -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *muon_invmass = new TH1D("muon_invmass","Reconstructed di-muon invariant mass (sig+bkg)", Nbins1 - 1, Edges1);
    	muon_invmass -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *muongamma_invmass_mother = new TH1D("muongamma_invmass_mother","Reconstructed muon pair + #gamma invariant mass with same mother", Nbins1 - 1, Edges1);
    	muongamma_invmass_mother -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *muon_invmass_mother = new TH1D("muon_invmass_mother","Reconstructed di-muon invariant mass (sig)", Nbins2 - 1, Edges2);
    	muon_invmass_mother -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_mother = new TH1D("eta_invmass_mother","Reconstructed invariant mass from muon pairs + #gamma (sig)",  Nbins2 - 1, Edges2);
    	eta_invmass_mother -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *muon_invmass_cuts_mother = new TH1D("muon_invmass_cuts_mother","Reconstructed di-muon invariant mass after all cuts (sig)", Nbins2 - 1, Edges2);
    	muon_invmass_cuts_mother -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_cuts_mother = new TH1D("eta_invmass_cuts_mother","Reconstructed invariant mass from muon pairs + #gamma after all cuts (sig)",  Nbins2 - 1, Edges2);
    	eta_invmass_cuts_mother -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *muon_invmass_cuts_mother2 = new TH1D("muon_invmass_cuts_mother2","Reconstructed di-muon invariant mass after all cuts including #gamma (sig)", Nbins2 - 1, Edges2);
    	muon_invmass_cuts_mother2 -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_cuts_mother2 = new TH1D("eta_invmass_cuts_mother2","Reconstructed invariant mass from muon pairs + #gamma after all cuts including #gamma (sig)",  Nbins2 - 1, Edges2);
    	eta_invmass_cuts_mother2 -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *muon_invmass_c3_mother = new TH1D("muon_invmass_c3_mother","Reconstructed muon pair invariant mass p_t, p & pseudorapidity cuts (sig)", Nbins2 - 1, Edges2);
    	muon_invmass_c3_mother -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_c3_mother = new TH1D("eta_invmass_c3_mother","Reconstructed invariant mass from muon pairs + #gamma p_t, p & pseudorapidity cuts (sig)",  Nbins2 - 1, Edges2);
    	eta_invmass_c3_mother -> GetXaxis()-> SetTitle("m (GeV)");


	TH1D *muon_invmass_cuts_sigbkg = new TH1D("muon_invmass_cuts_sigbkg","Reconstructed di-muon invariant mass after all cuts (sig+bkg)", Nbins2 - 1, Edges2);
    	muon_invmass_cuts_sigbkg -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_cuts_sigbkg = new TH1D("eta_invmass_cuts_sigbkg","Reconstructed invariant mass from muon pairs + #gamma after all cuts (sig+bkg)",  Nbins2 - 1, Edges2);
    	eta_invmass_cuts_sigbkg -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *muon_invmass_cuts_sigbkg2 = new TH1D("muon_invmass_cuts_sigbkg2","Reconstructed di-muon invariant mass after all cuts including #gamma (sig+bkg)", Nbins2 - 1, Edges2);
    	muon_invmass_cuts_sigbkg2 -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_cuts_sigbkg2 = new TH1D("eta_invmass_cuts_sigbkg2","Reconstructed invariant mass from muon pairs + #gamma after all cuts including #gamma (sig+bkg)",  Nbins2 - 1, Edges2);
    	eta_invmass_cuts_sigbkg2 -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *muon_invmass_c3_sigbkg = new TH1D("muon_invmass_c3_sigbkg","Reconstructed muon pair invariant mass p_t, p & pseudorapidity cuts (sig+bkg)", Nbins2 - 1, Edges2);
    	muon_invmass_cuts_sigbkg -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_c3_sigbkg = new TH1D("eta_invmass_c3_sigbkg","Reconstructed invariant mass from muon pairs + #gamma p_t, p & pseudorapidity cuts (sig+bkg)",  Nbins2 - 1, Edges2);
    	eta_invmass_cuts_sigbkg -> GetXaxis()-> SetTitle("m (GeV)");



	TH1D *muon_invmass_cuts_bkg = new TH1D("muon_invmass_cuts_bkg","Reconstructed di-muon invariant mass after all cuts (bkg)", Nbins2 - 1, Edges2);
    	muon_invmass_cuts_bkg -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_cuts_bkg = new TH1D("eta_invmass_cuts_bkg","Reconstructed invariant mass from muon pairs + #gamma after all cuts (bkg)",  Nbins2 - 1, Edges2);
    	eta_invmass_cuts_bkg -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *muon_invmass_cuts_bkg2 = new TH1D("muon_invmass_cuts_bkg2","Reconstructed di-muon invariant mass after all cuts including #gamma (bkg)", Nbins2 - 1, Edges2);
    	muon_invmass_cuts_bkg2 -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_cuts_bkg2 = new TH1D("eta_invmass_cuts_bkg2","Reconstructed invariant mass from muon pairs + #gamma after all cuts including #gamma (bkg)",  Nbins2 - 1, Edges2);
    	eta_invmass_cuts_bkg2 -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *muon_invmass_c3_bkg = new TH1D("muon_invmass_c3_bkg","Reconstructed muon pair invariant mass p_t, p & pseudorapidity cuts (bkg)", Nbins2 - 1, Edges2);
    	muon_invmass_cuts_bkg -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_c3_bkg = new TH1D("eta_invmass_c3_bkg","Reconstructed invariant mass from muon pairs + #gamma p_t, p & pseudorapidity cuts (bkg)",  Nbins2 - 1, Edges2);
    	eta_invmass_cuts_bkg -> GetXaxis()-> SetTitle("m (GeV)");



	TH1D *muon_pt_mother = new TH1D("muon_pt_mother","#mu^{-} and #mu^{+} p_{t} distribution with same mother (#eta)", 50, 0.0, 5.0);
    	muon_pt_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *muon_pt_bkg = new TH1D("muon_pt_bkg","#mu^{-} and #mu^{+} p_{t} distribution for background", 50, 0.0, 5.0);
    	muon_pt_bkg -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *muon_ptcut_mother = new TH1D("muon_ptcut_mother","#mu^{-} and #mu^{+} p_{t} distribution after p_{t} & p cuts", 50, 0.0, 5.0);
    	muon_ptcut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *muon_ptetacut_mother = new TH1D("muon_ptetacut_mother","#mu^{-} and #mu^{+} p_{t} distribution after p_{t}, p & pseudorapidity cuts", 50, 0.0, 5.0);
    	muon_ptetacut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *muon_ptetamasscut_mother = new TH1D("muon_ptetamasscut_mother","#mu^{-} and #mu^{+} p_{t} distribution after all cuts", 50, 0.0, 5.0);
    	muon_ptetamasscut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *muon_p_mother = new TH1D("muon_p_mother","#mu^{-} and #mu^{+} p distribution with same mother (#eta)", 60, 0.0, 30.0);
    	muon_p_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *muon_p_bkg = new TH1D("muon_p_bkg","#mu^{-} and #mu^{+} p distribution for background", 60, 0.0, 30.0);
    	muon_p_bkg -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *muon_pcut_mother = new TH1D("muon_pcut_mother","#mu^{-} and #mu^{+} p distribution after p_{t} & p cuts", 60, 0.0, 30.0);
    	muon_pcut_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *muon_petacut_mother = new TH1D("muon_petacut_mother","#mu^{-} and #mu^{+} p distribution after p_{t}, p & pseudorapidity cuts", 60, 0.0, 30.0);
    	muon_petacut_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *muon_petamasscut_mother = new TH1D("muon_petamasscut_mother","#mu^{-} and #mu^{+} p distribution after all cuts", 60, 0.0, 30.0);
    	muon_petamasscut_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *muon_eta_mother = new TH1D("muon_eta_mother","#mu^{-} and #mu^{+} pseudorapidity distribution with same mother (#eta)", 45, 0.0, 4.5);
    	muon_eta_mother -> GetXaxis()-> SetTitle("#eta");




	TH1D *gamma_pt_mother = new TH1D("gamma_pt_mother","#gamma p_{t} distribution with same mother (#eta)", 50, 0.0, 5.0);
    	gamma_pt_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *gamma_pt_bkg = new TH1D("gamma_pt_bkg","#gamma p_{t} distribution for background", 50, 0.0, 5.0);
    	gamma_pt_bkg -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *gamma_pt_c3_bkg = new TH1D("gamma_pt_c3_bkg","#gamma p_{t} distribution for background after p_{t}, p & pseudorapidity cuts", 50, 0.0, 5.0);
    	gamma_pt_c3_bkg -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *gamma_pt_c4_bkg = new TH1D("gamma_pt_c4_bkg","#gamma p_{t} distribution for background after mass cuts", 50, 0.0, 5.0);
    	gamma_pt_c4_bkg -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *gamma_pt_c5_bkg = new TH1D("gamma_pt_c5_bkg","#gamma p_{t} distribution for background after all cuts including #gamma", 50, 0.0, 5.0);
    	gamma_pt_c5_bkg -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *gamma_ptcut_mother = new TH1D("gamma_ptcut_mother","#gamma p_{t} distribution after p_{t} & p cuts", 50, 0.0, 5.0);
    	gamma_ptcut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *gamma_ptetacut_mother = new TH1D("gamma_ptetacut_mother","#gamma p_{t} distribution after p_{t}, p & pseudorapidity cuts", 50, 0.0, 5.0);
    	gamma_ptetacut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *gamma_ptetamasscut_mother = new TH1D("gamma_ptetamasscut_mother","#gamma p_{t} distribution after all cuts", 50, 0.0, 5.0);
    	gamma_ptetamasscut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *gamma_ptallcut_mother = new TH1D("gamma_ptallcut_mother","#gamma p_{t} distribution after all cuts including #gamma", 50, 0.0, 5.0);
    	gamma_ptallcut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *gamma_p_mother = new TH1D("gamma_p_mother","#gamma p distribution with same mother (#eta)", 60, 0.0, 30.0);
    	gamma_p_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *gamma_p_bkg = new TH1D("gamma_p_bkg","#gamma p distribution for background", 60, 0.0, 30.0);
    	gamma_p_bkg -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *gamma_p_c3_bkg = new TH1D("gamma_p_c3_bkg","#gamma p distribution for background after p_{t}, p & pseudorapidity cuts", 60, 0.0, 30.0);
    	gamma_p_c3_bkg -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *gamma_p_c4_bkg = new TH1D("gamma_p_c4_bkg","#gamma p distribution for background after mass cuts", 60, 0.0, 30.0);
    	gamma_p_c4_bkg -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *gamma_p_c5_bkg = new TH1D("gamma_p_c5_bkg","#gamma p distribution for background after all cuts including #gamma", 60, 0.0, 30.0);
    	gamma_p_c5_bkg -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *gamma_pcut_mother = new TH1D("gamma_pcut_mother","#gamma p distribution after p_{t} & p cuts", 60, 0.0, 30.0);
    	gamma_pcut_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *gamma_petacut_mother = new TH1D("gamma_petacut_mother","#gamma p distribution after p_{t}, p & pseudorapidity cuts", 60, 0.0, 30.0);
    	gamma_petacut_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *gamma_petamasscut_mother = new TH1D("gamma_petamasscut_mother","#gamma p distribution after all cuts", 60, 0.0, 30.0);
    	gamma_petamasscut_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *gamma_pallcut_mother = new TH1D("gamma_pallcut_mother","#gamma p distribution after all cuts including #gamma", 60, 0.0, 30.0);
    	gamma_pallcut_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *gamma_eta_mother = new TH1D("gamma_eta_mother","#gamma pseudorapidity distribution with same mother (#eta)", 45, 0.0, 4.5);
    	gamma_eta_mother -> GetXaxis()-> SetTitle("#eta");




/*
	TH2D *muon_ptvsp_mother = new TH2D("muon_ptvsp_mother","signal #mu^{#pm} p vs p_t distributions before cuts", 60, 0.0, 30.0, 50, 0.0, 5.0);
	muon_ptvsp_mother -> GetXaxis()-> SetTitle("p (GeV)");
	muon_ptvsp_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *muon_ptvsp_cuts_mother = new TH2D("muon_ptvsp_cuts_mother","signal #mu^{#pm} p vs p_t distributions after cuts", 60, 0.0, 30.0, 50, 0.0, 5.0);
	muon_ptvsp_cuts_mother -> GetXaxis()-> SetTitle("p (GeV)");
	muon_ptvsp_cuts_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *muon_ptvsp_c1_mother = new TH2D("muon_ptvsp_c1_mother","signal #mu^{#pm} p vs p_t distributions after p_t cuts", 60, 0.0, 30.0, 50, 0.0, 5.0);
	muon_ptvsp_c1_mother -> GetXaxis()-> SetTitle("p (GeV)");
	muon_ptvsp_c1_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *muon_ptvsp_c2_mother = new TH2D("muon_ptvsp_c2_mother","signal #mu^{#pm} p vs p_t distributions after p_t & p cuts", 60, 0.0, 30.0, 50, 0.0, 5.0);
	muon_ptvsp_c2_mother -> GetXaxis()-> SetTitle("p (GeV)");
	muon_ptvsp_c2_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *muon_ptvsp_c3_mother = new TH2D("muon_ptvsp_c3_mother","signal #mu^{#pm} p vs p_t distributions after p_t, p & pseudorapidity cuts", 60, 0.0, 30.0, 50, 0.0, 5.0);
	muon_ptvsp_c3_mother -> GetXaxis()-> SetTitle("p (GeV)");
	muon_ptvsp_c3_mother -> GetYaxis()-> SetTitle("p_t (GeV)");



	TH2D *muon_ptvseta_mother = new TH2D("muon_ptvseta_mother","signal #mu^{#pm} p_t vs pseudorapidity distributions before cuts", 50, 0.0, 5.0, 45, 0.0, 4.5);
	muon_ptvseta_mother -> GetXaxis()-> SetTitle("p_t (GeV)");
	muon_ptvseta_mother -> GetYaxis()-> SetTitle("#eta");

	TH2D *muon_ptvseta_cuts_mother = new TH2D("muon_ptvseta_cuts_mother","signal #mu^{#pm} p_t vs pseudorapidity distributions after all cuts", 50, 0.0, 5.0, 45, 0.0, 4.5);
	muon_ptvseta_cuts_mother -> GetXaxis()-> SetTitle("p_t (GeV)");
	muon_ptvseta_cuts_mother -> GetYaxis()-> SetTitle("#eta");

	TH2D *muon_ptvseta_c1_mother = new TH2D("muon_ptvseta_c1_mother","signal #mu^{#pm} p_t vs pseudorapidity distributions after p_t cuts", 50, 0.0, 5.0, 45, 0.0, 4.5);
	muon_ptvseta_c1_mother -> GetXaxis()-> SetTitle("p_t (GeV)");
	muon_ptvseta_c1_mother -> GetYaxis()-> SetTitle("#eta");

	TH2D *muon_ptvseta_c2_mother = new TH2D("muon_ptvseta_c2_mother","signal #mu^{#pm} p_t vs pseudorapidity distributions after p_t & p cuts", 50, 0.0, 5.0, 45, 0.0, 4.5);
	muon_ptvseta_c2_mother -> GetXaxis()-> SetTitle("p_t (GeV)");
	muon_ptvseta_c2_mother -> GetYaxis()-> SetTitle("#eta");

	TH2D *muon_ptvseta_c3_mother = new TH2D("muon_ptvseta_c3_mother","signal #mu^{#pm} p_t vs pseudorapidity distributions after p_t, p & pseudorapidity cuts", 50, 0.0, 5.0, 45, 0.0, 4.5);
	muon_ptvseta_c3_mother -> GetXaxis()-> SetTitle("p_t (GeV)");
	muon_ptvseta_c3_mother -> GetYaxis()-> SetTitle("#eta");



	TH2D *gamma_ptvsp_mother = new TH2D("gamma_ptvsp_mother","signal #gamma p vs p_t distributions before cuts", 60, 0.0, 30.0, 50, 0.0, 5.0);
	gamma_ptvsp_mother -> GetXaxis()-> SetTitle("p (GeV)");
	gamma_ptvsp_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *gamma_ptvsp_cuts_mother = new TH2D("gamma_ptvsp_cuts_mother","signal #gamma p vs p_t distributions after cuts", 60, 0.0, 30.0, 50, 0.0, 5.0);
	gamma_ptvsp_cuts_mother -> GetXaxis()-> SetTitle("p (GeV)");
	gamma_ptvsp_cuts_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *gamma_ptvsp_c1_mother = new TH2D("gamma_ptvsp_c1_mother","signal #gamma p vs p_t distributions after p_t cuts", 60, 0.0, 30.0, 50, 0.0, 5.0);
	gamma_ptvsp_c1_mother -> GetXaxis()-> SetTitle("p (GeV)");
	gamma_ptvsp_c1_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *gamma_ptvsp_c2_mother = new TH2D("gamma_ptvsp_c2_mother","signal #gamma p vs p_t distributions after p_t & p cuts", 60, 0.0, 30.0, 50, 0.0, 5.0);
	gamma_ptvsp_c2_mother -> GetXaxis()-> SetTitle("p (GeV)");
	gamma_ptvsp_c2_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *gamma_ptvsp_c3_mother = new TH2D("gamma_ptvsp_c3_mother","signal #gamma p vs p_t distributions after p_t, p & pseudorapidity cuts", 60, 0.0, 30.0, 50, 0.0, 5.0);
	gamma_ptvsp_c3_mother -> GetXaxis()-> SetTitle("p (GeV)");
	gamma_ptvsp_c3_mother -> GetYaxis()-> SetTitle("p_t (GeV)");



	TH2D *gamma_ptvseta_mother = new TH2D("gamma_ptvseta_mother","signal #gamma p_t vs pseudorapidity distributions before cuts", 50, 0.0, 5.0, 45, 0.0, 4.5);
	gamma_ptvseta_mother -> GetXaxis()-> SetTitle("p_t (GeV)");
	gamma_ptvseta_mother -> GetYaxis()-> SetTitle("#eta");

	TH2D *gamma_ptvseta_cuts_mother = new TH2D("gamma_ptvseta_cuts_mother","signal #gamma p_t vs pseudorapidity distributions after all cuts", 50, 0.0, 5.0, 45, 0.0, 4.5);
	gamma_ptvseta_cuts_mother -> GetXaxis()-> SetTitle("p_t (GeV)");
	gamma_ptvseta_cuts_mother -> GetYaxis()-> SetTitle("#eta");

	TH2D *gamma_ptvseta_c1_mother = new TH2D("gamma_ptvseta_c1_mother","signal #gamma p_t vs pseudorapidity distributions after p_t cuts", 50, 0.0, 5.0, 45, 0.0, 4.5);
	gamma_ptvseta_c1_mother -> GetXaxis()-> SetTitle("p_t (GeV)");
	gamma_ptvseta_c1_mother -> GetYaxis()-> SetTitle("#eta");

	TH2D *gamma_ptvseta_c2_mother = new TH2D("gamma_ptvseta_c2_mother","signal #gamma p_t vs pseudorapidity distributions after p_t & p cuts", 50, 0.0, 5.0, 45, 0.0, 4.5);
	gamma_ptvseta_c2_mother -> GetXaxis()-> SetTitle("p_t (GeV)");
	gamma_ptvseta_c2_mother -> GetYaxis()-> SetTitle("#eta");

	TH2D *gamma_ptvseta_c3_mother = new TH2D("gamma_ptvseta_c3_mother","signal #gamma p_t vs pseudorapidity distributions after p_t, p & pseudorapidity cuts", 50, 0.0, 5.0, 45, 0.0, 4.5);
	gamma_ptvseta_c3_mother -> GetXaxis()-> SetTitle("p_t (GeV)");
	gamma_ptvseta_c3_mother -> GetYaxis()-> SetTitle("#eta");
*/


	TH1D *mu_number_event = new TH1D("mu_number_event","Muon-antimuon number per event", 100000, 0.0, 100000.0);
    	mu_number_event -> GetXaxis()-> SetTitle("event index");
	mu_number_event -> GetYaxis()-> SetTitle("number of #mu^{#pm}");

	TH1D *gamma_number_event = new TH1D("gamma_number_event","#gamma number per event", 100000, 0.0, 100000.0);
    	gamma_number_event -> GetXaxis()-> SetTitle("event index");
	gamma_number_event -> GetYaxis()-> SetTitle("number of #gamma");

	TH1D *total_number_event = new TH1D("total_number_event","Total number of particles per event", 100000, 0.0, 100000.0);
    	total_number_event -> GetXaxis()-> SetTitle("event index");
	total_number_event -> GetYaxis()-> SetTitle("number of particles");

	TH1D *mugamma_combnumber_event = new TH1D("mugamma_combnumber_event","#mu^{#pm} & #gamma number of combinations per event", 100000, 0.0, 100000.0);
    	mugamma_combnumber_event -> GetXaxis()-> SetTitle("event index");
	mugamma_combnumber_event -> GetYaxis()-> SetTitle("number of combinations");

	TH1D *gamma_fake = new TH1D("gamma_fake","Parent ID of background #gamma passing the cuts", 2000, -1000.0, 1000.0);
    	gamma_fake -> GetXaxis()-> SetTitle("#gamma parent ID");


	TH1D *pigamma_invmass = new TH1D("pigamma_invmass","Reconstructed #eta invariant mass from misID #pi^{+}+#pi^{-}+#gamma", Nbins1 - 1, Edges1);
    	pigamma_invmass -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *pigamma_invmass_c3 = new TH1D("pigamma_invmass_c3","Reconstructed #eta invariant mass from misID #pi^{+}+#pi^{-}+#gamma after p_t, p & #eta cuts", Nbins2 - 1, Edges2);
    	pigamma_invmass_c3 -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *pigamma_invmass_cuts = new TH1D("pigamma_invmass_cuts","Reconstructed #eta invariant mass from misID #pi^{+}+#pi^{-}+#gamma after all cuts", Nbins2 - 1, Edges2);
    	pigamma_invmass_cuts -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *pigamma_invmass_cuts2 = new TH1D("pigamma_invmass_cuts2","Reconstructed #eta invariant mass from misID #pi^{+}+#pi^{-}+#gamma after all cuts including #gamma", Nbins2 - 1, Edges2);
    	pigamma_invmass_cuts2 -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *pi_invmass = new TH1D("pi_invmass","Reconstructed invariant mass from misID #pi^{#pm}", Nbins1 - 1, Edges1);
    	pi_invmass -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *pi_invmass_c3 = new TH1D("pi_invmass_c3","Reconstructed invariant mass from misID #pi^{#pm} after p_t, p & pseudorapidity cuts", Nbins2 - 1, Edges2);
    	pi_invmass_c3 -> GetXaxis()-> SetTitle("m (GeV)");
	
	TH1D *pi_invmass_cuts = new TH1D("pi_invmass_cuts","Reconstructed invariant mass from misID #pi^{#pm} after all cuts", Nbins2 - 1, Edges2);
    	pi_invmass_cuts -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *pi_invmass_cuts2 = new TH1D("pi_invmass_cuts2","Reconstructed invariant mass from misID #pi^{#pm} after all cuts including #gamma", Nbins2 - 1, Edges2);
    	pi_invmass_cuts2 -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *pi_pt = new TH1D("pi_pt","#pi^{-} and #pi^{+} p_{t} distribution", 50, 0.0, 5.0);
    	pi_pt -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *pi_ptcut = new TH1D("pi_ptcut","#pi^{-} and #pi^{+} p_{t} distribution after cuts", 50, 0.0, 5.0);
    	pi_ptcut -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *pi_ptmasscut = new TH1D("pi_ptmasscut","#pi^{-} and #pi^{+} p_{t} distribution after cuts including mass", 50, 0.0, 5.0);
    	pi_ptmasscut -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *pi_p = new TH1D("pi_p","#pi^{-} and #pi^{+} p distribution", 60, 0.0, 30.0);
    	pi_p -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *pi_pcut = new TH1D("pi_pcut","#pi^{-} and #pi^{+} p distribution after cuts", 60, 0.0, 30.0);
    	pi_pcut -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *pi_pmasscut = new TH1D("pi_pmasscut","#pi^{-} and #pi^{+} p distribution after cuts including mass", 60, 0.0, 30.0);
    	pi_pmasscut -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *pi_eta = new TH1D("pi_eta","#pi^{-} and #pi^{+} pseudorapidity distribution", 45, 0.0, 4.5);
    	pi_eta -> GetXaxis()-> SetTitle("#eta");



/*
	TH2D *pi_ptvsp = new TH2D("pi_ptvsp","#pi^{#pm} p vs p_t distributions before cuts", 60, 0.0, 30.0, 50, 0.0, 5.0);
	pi_ptvsp -> GetXaxis()-> SetTitle("p (GeV)");
	pi_ptvsp -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *pi_ptvsp_cuts = new TH2D("pi_ptvsp_cuts","#pi^{#pm} p vs p_t distributions after cuts", 60, 0.0, 30.0, 50, 0.0, 5.0);
	pi_ptvsp_cuts -> GetXaxis()-> SetTitle("p (GeV)");
	pi_ptvsp_cuts -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *pi_ptvsp_c1 = new TH2D("pi_ptvsp_c1","#pi^{#pm} p vs p_t distributions after p_t cuts", 60, 0.0, 30.0, 50, 0.0, 5.0);
	pi_ptvsp_c1 -> GetXaxis()-> SetTitle("p (GeV)");
	pi_ptvsp_c1 -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *pi_ptvsp_c2 = new TH2D("pi_ptvsp_c2","#pi^{#pm} p vs p_t distributions after p_t & p cuts", 60, 0.0, 30.0, 50, 0.0, 5.0);
	pi_ptvsp_c2 -> GetXaxis()-> SetTitle("p (GeV)");
	pi_ptvsp_c2 -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *pi_ptvsp_c3 = new TH2D("pi_ptvsp_c3","#pi^{#pm} p vs p_t distributions after p_t, p & pseudorapidity cuts", 60, 0.0, 30.0, 50, 0.0, 5.0);
	pi_ptvsp_c3 -> GetXaxis()-> SetTitle("p (GeV)");
	pi_ptvsp_c3 -> GetYaxis()-> SetTitle("p_t (GeV)");



	TH2D *pi_ptvseta = new TH2D("pi_ptvseta","#pi^{#pm} p_t vs pseudorapidity distributions before cuts", 50, 0.0, 5.0, 45, 0.0, 4.5);
	pi_ptvseta -> GetXaxis()-> SetTitle("p_t (GeV)");
	pi_ptvseta -> GetYaxis()-> SetTitle("#eta");

	TH2D *pi_ptvseta_cuts = new TH2D("pi_ptvseta_cuts","#pi^{#pm} p_t vs pseudorapidity distributions after all cuts", 50, 0.0, 5.0, 45, 0.0, 4.5);
	pi_ptvseta_cuts -> GetXaxis()-> SetTitle("p_t (GeV)");
	pi_ptvseta_cuts -> GetYaxis()-> SetTitle("#eta");

	TH2D *pi_ptvseta_c1 = new TH2D("pi_ptvseta_c1","#pi^{#pm} p_t vs pseudorapidity distributions after p_t cuts", 50, 0.0, 5.0, 45, 0.0, 4.5);
	pi_ptvseta_c1 -> GetXaxis()-> SetTitle("p_t (GeV)");
	pi_ptvseta_c1 -> GetYaxis()-> SetTitle("#eta");

	TH2D *pi_ptvseta_c2 = new TH2D("pi_ptvseta_c2","#pi^{#pm} p_t vs pseudorapidity distributions after p_t & p cuts", 50, 0.0, 5.0, 45, 0.0, 4.5);
	pi_ptvseta_c2 -> GetXaxis()-> SetTitle("p_t (GeV)");
	pi_ptvseta_c2 -> GetYaxis()-> SetTitle("#eta");

	TH2D *pi_ptvseta_c3 = new TH2D("pi_ptvseta_c3","#pi^{#pm} p_t vs pseudorapidity distributions after p_t, p & pseudorapidity cuts", 50, 0.0, 5.0, 45, 0.0, 4.5);
	pi_ptvseta_c3 -> GetXaxis()-> SetTitle("p_t (GeV)");
	pi_ptvseta_c3 -> GetYaxis()-> SetTitle("#eta");
*/


	TH1D *pi_number_event = new TH1D("pi_number_event","#pi^{#pm} number per event", 1000, 0.0, 1000.0);
    	pi_number_event -> GetXaxis()-> SetTitle("event index");
	pi_number_event -> GetYaxis()-> SetTitle("number of #pi^{#pm}");

	TH1D *pi_combnumber_event = new TH1D("pi_combnumber_event","#pi^{#pm} number of combinations  per event", 1000, 0.0, 1000.0);
    	pi_combnumber_event -> GetXaxis()-> SetTitle("event index");
	pi_combnumber_event -> GetYaxis()-> SetTitle("number of combinations");

	TH1D *total_pigamma_number_event = new TH1D("total_pigamma_number_event","#pi^{+}+#pi^{-}+#gamma number per event", 1000, 0.0, 1000.0);
    	total_pigamma_number_event -> GetXaxis()-> SetTitle("event index");
	total_pigamma_number_event -> GetYaxis()-> SetTitle("#pi^{+}+#pi^{-}+#gamma");

	TH1D *pigamma_combnumber_event = new TH1D("pigamma_combnumber_event","#pi^{+}+#pi^{-}+#gamma number of combinations  per event", 1000, 0.0, 1000.0);
    	pigamma_combnumber_event -> GetXaxis()-> SetTitle("event index");
	pigamma_combnumber_event -> GetYaxis()-> SetTitle("number of combinations");

	vect* v= new vect;//declares a vector in which we are going to store all the particles from the same event
	double prev_index = 0.0;//saving the event index of the previous entry of the tree

	vect2* vn = new vect2;
	vn->n = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	double mu_ptcut = 0.5;
	double mu_pcut = 10.0;
	double eta_cut1 = 2.0;
	double eta_cut2 = 4.5;
	double gamma_ptcut = 0.5;
	double gamma_pcut = 10.0;

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

			//analyze_event(file_out1, v, vn, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut, gamma_pcut, eta_invmass, muon_invmass, muongamma_invmass_mother, muon_invmass_mother, eta_invmass_mother, muon_invmass_cuts_mother, muon_invmass_cuts_mother2, eta_invmass_cuts_mother, eta_invmass_cuts_mother2, muon_invmass_c3_mother, eta_invmass_c3_mother, muon_invmass_cuts_bkg, muon_invmass_cuts_bkg2, eta_invmass_cuts_bkg, eta_invmass_cuts_bkg2, muon_invmass_c3_bkg, eta_invmass_c3_bkg, muon_pt_mother, muon_pt_bkg, muon_ptcut_mother, muon_ptetacut_mother, muon_ptetamasscut_mother, muon_p_mother, muon_p_bkg, muon_pcut_mother, muon_petacut_mother, muon_petamasscut_mother, muon_ptvsp_mother, muon_ptvsp_cuts_mother, muon_ptvsp_c1_mother, muon_ptvsp_c1_mother, muon_ptvsp_c2_mother, gamma_ptvsp_mother, gamma_ptvsp_cuts_mother, gamma_ptvsp_c1_mother, gamma_ptvsp_c2_mother, gamma_ptvsp_c3_mother, muon_ptvseta_mother, muon_ptvseta_cuts_mother, muon_ptvseta_c1_mother, muon_ptvseta_c2_mother, muon_ptvseta_c3_mother, gamma_ptvseta_mother, gamma_ptvseta_cuts_mother, gamma_ptvseta_c1_mother, gamma_ptvseta_c2_mother, gamma_ptvseta_c3_mother, muon_eta_mother, gamma_pt_mother, gamma_pt_bkg, gamma_pt_c3_bkg, gamma_pt_c4_bkg, gamma_pt_c5_bkg, gamma_ptcut_mother, gamma_ptetacut_mother, gamma_ptetamasscut_mother, gamma_ptallcut_mother, gamma_p_mother, gamma_p_bkg, gamma_p_c3_bkg, gamma_p_c4_bkg, gamma_p_c5_bkg, gamma_pcut_mother, gamma_petacut_mother, gamma_petamasscut_mother, gamma_pallcut_mother, gamma_eta_mother, mu_number_event, gamma_number_event, total_number_event, mugamma_combnumber_event, gamma_fake);

			analyze_event(file_out1, v, vn, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut, gamma_pcut, eta_invmass, muon_invmass, muongamma_invmass_mother, muon_invmass_mother, eta_invmass_mother, muon_invmass_cuts_mother, muon_invmass_cuts_mother2, eta_invmass_cuts_mother, eta_invmass_cuts_mother2, muon_invmass_c3_mother, eta_invmass_c3_mother, muon_invmass_cuts_sigbkg, muon_invmass_cuts_sigbkg2, eta_invmass_cuts_sigbkg, eta_invmass_cuts_sigbkg2, muon_invmass_c3_sigbkg, eta_invmass_c3_sigbkg, muon_invmass_cuts_bkg, muon_invmass_cuts_bkg2, eta_invmass_cuts_bkg, eta_invmass_cuts_bkg2, muon_invmass_c3_bkg, eta_invmass_c3_bkg, muon_pt_mother, muon_pt_bkg, muon_ptcut_mother, muon_ptetacut_mother, muon_ptetamasscut_mother, muon_p_mother, muon_p_bkg, muon_pcut_mother, muon_petacut_mother, muon_petamasscut_mother, muon_eta_mother, gamma_pt_mother, gamma_pt_bkg, gamma_pt_c3_bkg, gamma_pt_c4_bkg, gamma_pt_c5_bkg, gamma_ptcut_mother, gamma_ptetacut_mother, gamma_ptetamasscut_mother, gamma_ptallcut_mother, gamma_p_mother, gamma_p_bkg, gamma_p_c3_bkg, gamma_p_c4_bkg, gamma_p_c5_bkg, gamma_pcut_mother, gamma_petacut_mother, gamma_petamasscut_mother, gamma_pallcut_mother, gamma_eta_mother, mu_number_event, gamma_number_event, total_number_event, mugamma_combnumber_event, gamma_fake);

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

			//analyze_event(file_out1, v, vn, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut, gamma_pcut, eta_invmass, muon_invmass, muongamma_invmass_mother, muon_invmass_mother, eta_invmass_mother, muon_invmass_cuts_mother, muon_invmass_cuts_mother2, eta_invmass_cuts_mother, eta_invmass_cuts_mother2, muon_invmass_c3_mother, eta_invmass_c3_mother, muon_invmass_cuts_bkg, muon_invmass_cuts_bkg2, eta_invmass_cuts_bkg, eta_invmass_cuts_bkg2, muon_invmass_c3_bkg, eta_invmass_c3_bkg, muon_pt_mother, muon_pt_bkg, muon_ptcut_mother, muon_ptetacut_mother, muon_ptetamasscut_mother, muon_p_mother, muon_p_bkg, muon_pcut_mother, muon_petacut_mother, muon_petamasscut_mother, muon_ptvsp_mother, muon_ptvsp_cuts_mother, muon_ptvsp_c1_mother, muon_ptvsp_c1_mother, muon_ptvsp_c2_mother, gamma_ptvsp_mother, gamma_ptvsp_cuts_mother, gamma_ptvsp_c1_mother, gamma_ptvsp_c2_mother, gamma_ptvsp_c3_mother, muon_ptvseta_mother, muon_ptvseta_cuts_mother, muon_ptvseta_c1_mother, muon_ptvseta_c2_mother, muon_ptvseta_c3_mother, gamma_ptvseta_mother, gamma_ptvseta_cuts_mother, gamma_ptvseta_c1_mother, gamma_ptvseta_c2_mother, gamma_ptvseta_c3_mother, muon_eta_mother, gamma_pt_mother, gamma_pt_bkg, gamma_pt_c3_bkg, gamma_pt_c4_bkg, gamma_pt_c5_bkg, gamma_ptcut_mother, gamma_ptetacut_mother, gamma_ptetamasscut_mother, gamma_ptallcut_mother, gamma_p_mother, gamma_p_bkg, gamma_p_c3_bkg, gamma_p_c4_bkg, gamma_p_c5_bkg, gamma_pcut_mother, gamma_petacut_mother, gamma_petamasscut_mother, gamma_pallcut_mother, gamma_eta_mother, mu_number_event, gamma_number_event, total_number_event, mugamma_combnumber_event, gamma_fake);

			analyze_event(file_out1, v, vn, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut, gamma_pcut, eta_invmass, muon_invmass, muongamma_invmass_mother, muon_invmass_mother, eta_invmass_mother, muon_invmass_cuts_mother, muon_invmass_cuts_mother2, eta_invmass_cuts_mother, eta_invmass_cuts_mother2, muon_invmass_c3_mother, eta_invmass_c3_mother, muon_invmass_cuts_sigbkg, muon_invmass_cuts_sigbkg2, eta_invmass_cuts_sigbkg, eta_invmass_cuts_sigbkg2, muon_invmass_c3_sigbkg, eta_invmass_c3_sigbkg, muon_invmass_cuts_bkg, muon_invmass_cuts_bkg2, eta_invmass_cuts_bkg, eta_invmass_cuts_bkg2, muon_invmass_c3_bkg, eta_invmass_c3_bkg, muon_pt_mother, muon_pt_bkg, muon_ptcut_mother, muon_ptetacut_mother, muon_ptetamasscut_mother, muon_p_mother, muon_p_bkg, muon_pcut_mother, muon_petacut_mother, muon_petamasscut_mother, muon_eta_mother, gamma_pt_mother, gamma_pt_bkg, gamma_pt_c3_bkg, gamma_pt_c4_bkg, gamma_pt_c5_bkg, gamma_ptcut_mother, gamma_ptetacut_mother, gamma_ptetamasscut_mother, gamma_ptallcut_mother, gamma_p_mother, gamma_p_bkg, gamma_p_c3_bkg, gamma_p_c4_bkg, gamma_p_c5_bkg, gamma_pcut_mother, gamma_petacut_mother, gamma_petamasscut_mother, gamma_pallcut_mother, gamma_eta_mother, mu_number_event, gamma_number_event, total_number_event, mugamma_combnumber_event, gamma_fake);

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

	vect2* vp = new vect2;
	vp->n = {0.0, 0.0, 0.0, 0.0};

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

			//analyze_event_pi(file_out2, v, vp, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut, gamma_pcut, pigamma_invmass, pigamma_invmass_c3, pigamma_invmass_cuts, pigamma_invmass_cuts2, pi_invmass, pi_invmass_c3, pi_invmass_cuts, pi_invmass_cuts2, pi_pt, pi_p, pi_eta, pi_ptcut, pi_pcut, pi_ptmasscut, pi_pmasscut, pi_ptvsp, pi_ptvsp_cuts, pi_number_event, pi_combnumber_event, total_pigamma_number_event, pigamma_combnumber_event);

			analyze_event_pi(file_out2, v, vp, alpha, beta, gamma, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut, gamma_pcut, pigamma_invmass, pigamma_invmass_c3, pigamma_invmass_cuts, pigamma_invmass_cuts2, pi_invmass, pi_invmass_c3, pi_invmass_cuts, pi_invmass_cuts2, pi_pt, pi_p, pi_eta, pi_ptcut, pi_pcut, pi_ptmasscut, pi_pmasscut, pi_number_event, pi_combnumber_event, total_pigamma_number_event, pigamma_combnumber_event);

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

			//analyze_event_pi(file_out2, v, vp, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut, gamma_pcut, pigamma_invmass, pigamma_invmass_c3, pigamma_invmass_cuts, pigamma_invmass_cuts2, pi_invmass, pi_invmass_c3, pi_invmass_cuts, pi_invmass_cuts2, pi_pt, pi_p, pi_eta, pi_ptcut, pi_pcut, pi_ptmasscut, pi_pmasscut, pi_ptvsp, pi_ptvsp_cuts, pi_number_event, pi_combnumber_event, total_pigamma_number_event, pigamma_combnumber_event);

			analyze_event_pi(file_out2, v, vp, alpha, beta, gamma, mu_ptcut, mu_pcut, eta_cut1, eta_cut2, gamma_ptcut, gamma_pcut, pigamma_invmass, pigamma_invmass_c3, pigamma_invmass_cuts, pigamma_invmass_cuts2, pi_invmass, pi_invmass_c3, pi_invmass_cuts, pi_invmass_cuts2, pi_pt, pi_p, pi_eta, pi_ptcut, pi_pcut, pi_ptmasscut, pi_pmasscut, pi_number_event, pi_combnumber_event, total_pigamma_number_event, pigamma_combnumber_event);

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
   
	TCanvas *c1=new TCanvas("c1","",600,600);

	nentries = muon_pt_mother->GetEntries();
	muon_pt_mother->Scale(1.0 / nentries, "width");
	muon_pt_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_pt_mother1.pdf","pdf");

	nentries = muon_pt_bkg->GetEntries();
	muon_pt_bkg->Scale(1.0 / nentries, "width");
	muon_pt_bkg->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_pt_bkg1.pdf","pdf");

	nentries = muon_ptcut_mother->GetEntries();
	muon_ptcut_mother->Scale(1.0 / nentries, "width");
	muon_ptcut_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_ptcut_mother1.pdf","pdf");

	nentries = muon_ptetacut_mother->GetEntries();
	muon_ptetacut_mother->Scale(1.0 / nentries, "width");
	muon_ptetacut_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_ptetacut_mother1.pdf","pdf");

	nentries = muon_ptetamasscut_mother->GetEntries();
	muon_ptetamasscut_mother->Scale(1.0 / nentries, "width");
	muon_ptetamasscut_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_ptetamasscut_mother1.pdf","pdf");

	nentries = muon_p_mother->GetEntries();
	muon_p_mother->Scale(1.0 / nentries, "width");
	muon_p_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_p_mother1.pdf","pdf");

	nentries = muon_p_bkg->GetEntries();
	muon_p_bkg->Scale(1.0 / nentries, "width");
	muon_p_bkg->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_p_bkg1.pdf","pdf");

	nentries = muon_pcut_mother->GetEntries();
	muon_pcut_mother->Scale(1.0 / nentries, "width");
	muon_pcut_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_pcut_mother1.pdf","pdf");

	nentries = muon_petacut_mother->GetEntries();
	muon_petacut_mother->Scale(1.0 / nentries, "width");
	muon_petacut_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_petacut_mother1.pdf","pdf");

	nentries = muon_petamasscut_mother->GetEntries();
	muon_petamasscut_mother->Scale(1.0 / nentries, "width");
	muon_petamasscut_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_petamasscut_mother1.pdf","pdf");

/*
	//Scaling???
	nentries = muon_ptvsp_mother->GetEntries();
	muon_ptvsp_mother->Scale(1.0 / nentries, "width");
	muon_ptvsp_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_ptvsp_mother1.pdf","pdf");


	//Scaling???
	nentries = muon_ptvsp_cuts_mother->GetEntries();
	muon_ptvsp_cuts_mother->Scale(1.0 / nentries, "width");
	muon_ptvsp_cuts_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_ptvsp_cuts_mother1.pdf","pdf");

	nentries = muon_ptvsp_c1_mother->GetEntries();
	muon_ptvsp_c1_mother->Scale(1.0 / nentries, "width");
	muon_ptvsp_c1_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_ptvsp_c1_mother1.pdf","pdf");

	nentries = muon_ptvsp_c2_mother->GetEntries();
	muon_ptvsp_c2_mother->Scale(1.0 / nentries, "width");
	muon_ptvsp_c2_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_ptvsp_c2_mother1.pdf","pdf");

	nentries = muon_ptvsp_c3_mother->GetEntries();
	muon_ptvsp_c3_mother->Scale(1.0 / nentries, "width");
	muon_ptvsp_c3_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_ptvsp_c3_mother1.pdf","pdf");


	nentries = muon_ptvseta_mother->GetEntries();
	muon_ptvseta_mother->Scale(1.0 / nentries, "width");
	muon_ptvseta_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_ptvseta_mother1.pdf","pdf");

	nentries = muon_ptvseta_cuts_mother->GetEntries();
	muon_ptvseta_cuts_mother->Scale(1.0 / nentries, "width");
	muon_ptvseta_cuts_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_ptvseta_cuts_mother1.pdf","pdf");

	nentries = muon_ptvseta_c1_mother->GetEntries();
	muon_ptvseta_c1_mother->Scale(1.0 / nentries, "width");
	muon_ptvseta_c1_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_ptvseta_c1_mother1.pdf","pdf");

	nentries = muon_ptvseta_c2_mother->GetEntries();
	muon_ptvseta_c2_mother->Scale(1.0 / nentries, "width");
	muon_ptvseta_c2_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_ptvseta_c2_mother1.pdf","pdf");

	nentries = muon_ptvseta_c3_mother->GetEntries();
	muon_ptvseta_c3_mother->Scale(1.0 / nentries, "width");
	muon_ptvseta_c3_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_ptvseta_c3_mother1.pdf","pdf");
*/

	nentries = muon_eta_mother->GetEntries();
	muon_eta_mother->Scale(1.0 / nentries, "width");
	muon_eta_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_pseudorap_mother1.pdf","pdf");





	nentries = gamma_pt_mother->GetEntries();
	gamma_pt_mother->Scale(1.0 / nentries, "width");
	gamma_pt_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_pt_mother1.pdf","pdf");

	nentries = gamma_pt_bkg->GetEntries();
	gamma_pt_bkg->Scale(1.0 / nentries, "width");
	gamma_pt_bkg->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_pt_bkg1.pdf","pdf");

	nentries = gamma_pt_c3_bkg->GetEntries();
	gamma_pt_c3_bkg->Scale(1.0 / nentries, "width");
	gamma_pt_c3_bkg->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_pt_c3_bkg1.pdf","pdf");

	nentries = gamma_pt_c4_bkg->GetEntries();
	gamma_pt_c4_bkg->Scale(1.0 / nentries, "width");
	gamma_pt_c4_bkg->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_pt_c4_bkg1.pdf","pdf");

	nentries = gamma_pt_c5_bkg->GetEntries();
	gamma_pt_c5_bkg->Scale(1.0 / nentries, "width");
	gamma_pt_c5_bkg->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_pt_c5_bkg1.pdf","pdf");

	nentries = gamma_ptcut_mother->GetEntries();
	gamma_ptcut_mother->Scale(1.0 / nentries, "width");
	gamma_ptcut_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_ptcut_mother1.pdf","pdf");

	nentries = gamma_ptetacut_mother->GetEntries();
	gamma_ptetacut_mother->Scale(1.0 / nentries, "width");
	gamma_ptetacut_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_ptetacut_mother1.pdf","pdf");

	nentries = gamma_ptetamasscut_mother->GetEntries();
	gamma_ptetamasscut_mother->Scale(1.0 / nentries, "width");
	gamma_ptetamasscut_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_ptetamasscut_mother1.pdf","pdf");

	nentries = gamma_ptallcut_mother->GetEntries();
	gamma_ptallcut_mother->Scale(1.0 / nentries, "width");
	gamma_ptallcut_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_ptallcut_mother1.pdf","pdf");

	nentries = gamma_p_mother->GetEntries();
	gamma_p_mother->Scale(1.0 / nentries, "width");
	gamma_p_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_p_mother1.pdf","pdf");

	nentries = gamma_p_bkg->GetEntries();
	gamma_p_bkg->Scale(1.0 / nentries, "width");
	gamma_p_bkg->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_p_bkg1.pdf","pdf");

	nentries = gamma_p_c3_bkg->GetEntries();
	gamma_p_c3_bkg->Scale(1.0 / nentries, "width");
	gamma_p_c3_bkg->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_p_c3_bkg1.pdf","pdf");

	nentries = gamma_p_c4_bkg->GetEntries();
	gamma_p_c4_bkg->Scale(1.0 / nentries, "width");
	gamma_p_c4_bkg->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_p_c4_bkg1.pdf","pdf");

	nentries = gamma_p_c5_bkg->GetEntries();
	gamma_p_c5_bkg->Scale(1.0 / nentries, "width");
	gamma_p_c5_bkg->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_p_c5_bkg1.pdf","pdf");

	nentries = gamma_pcut_mother->GetEntries();
	gamma_pcut_mother->Scale(1.0 / nentries, "width");
	gamma_pcut_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_pcut_mother1.pdf","pdf");

	nentries = gamma_petacut_mother->GetEntries();
	gamma_petacut_mother->Scale(1.0 / nentries, "width");
	gamma_petacut_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_petacut_mother1.pdf","pdf");

	nentries = gamma_petamasscut_mother->GetEntries();
	gamma_petamasscut_mother->Scale(1.0 / nentries, "width");
	gamma_petamasscut_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_petamasscut_mother1.pdf","pdf");

	nentries = gamma_pallcut_mother->GetEntries();
	gamma_pallcut_mother->Scale(1.0 / nentries, "width");
	gamma_pallcut_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_pallcut_mother1.pdf","pdf");

/*
	//Scaling???
	nentries = gamma_ptvsp_mother->GetEntries();
	gamma_ptvsp_mother->Scale(1.0 / nentries, "width");
	gamma_ptvsp_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_ptvsp_mother1.pdf","pdf");

	nentries = gamma_ptvsp_cuts_mother->GetEntries();
	gamma_ptvsp_cuts_mother->Scale(1.0 / nentries, "width");
	gamma_ptvsp_cuts_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_ptvsp_cuts_mother1.pdf","pdf");

	nentries = gamma_ptvsp_c1_mother->GetEntries();
	gamma_ptvsp_c1_mother->Scale(1.0 / nentries, "width");
	gamma_ptvsp_c1_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_ptvsp_c1_mother1.pdf","pdf");

	nentries = gamma_ptvsp_c2_mother->GetEntries();
	gamma_ptvsp_c2_mother->Scale(1.0 / nentries, "width");
	gamma_ptvsp_c2_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_ptvsp_c2_mother1.pdf","pdf");

	nentries = gamma_ptvsp_c3_mother->GetEntries();
	gamma_ptvsp_c3_mother->Scale(1.0 / nentries, "width");
	gamma_ptvsp_c3_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_ptvsp_c3_mother1.pdf","pdf");


	nentries = gamma_ptvseta_mother->GetEntries();
	gamma_ptvseta_mother->Scale(1.0 / nentries, "width");
	gamma_ptvseta_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_ptvseta_mother1.pdf","pdf");

	nentries = gamma_ptvseta_cuts_mother->GetEntries();
	gamma_ptvseta_cuts_mother->Scale(1.0 / nentries, "width");
	gamma_ptvseta_cuts_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_ptvseta_cuts_mother1.pdf","pdf");

	nentries = gamma_ptvseta_c1_mother->GetEntries();
	gamma_ptvseta_c1_mother->Scale(1.0 / nentries, "width");
	gamma_ptvseta_c1_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_ptvseta_c1_mother1.pdf","pdf");

	nentries = gamma_ptvseta_c2_mother->GetEntries();
	gamma_ptvseta_c2_mother->Scale(1.0 / nentries, "width");
	gamma_ptvseta_c2_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_ptvseta_c2_mother1.pdf","pdf");

	nentries = gamma_ptvseta_c3_mother->GetEntries();
	gamma_ptvseta_c3_mother->Scale(1.0 / nentries, "width");
	gamma_ptvseta_c3_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_ptvseta_c3_mother1.pdf","pdf");
*/

	nentries = gamma_eta_mother->GetEntries();
	gamma_eta_mother->Scale(1.0 / nentries, "width");
	gamma_eta_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_pseudorap_mother1.pdf","pdf");






	mu_number_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_number_event.pdf","pdf");

	mugamma_combnumber_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_mugamma_combnumber_event.pdf","pdf");

	nentries = pi_pt->GetEntries();
	pi_pt->Scale(1.0 / nentries, "width");
	pi_pt->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_pt1.pdf","pdf");

	nentries = pi_p->GetEntries();
	pi_p->Scale(1.0 / nentries, "width");
	pi_p->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_p1.pdf","pdf");

	nentries = pi_eta->GetEntries();
	pi_eta->Scale(1.0 / nentries, "width");
	pi_eta->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_eta1.pdf","pdf");

	nentries = pi_ptcut->GetEntries();
	pi_ptcut->Scale(1.0 / nentries, "width");
	pi_ptcut->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_ptcut1.pdf","pdf");

	nentries = pi_pcut->GetEntries();
	pi_pcut->Scale(1.0 / nentries, "width");
	pi_pcut->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_pcut1.pdf","pdf");

	nentries = pi_ptmasscut->GetEntries();
	pi_ptmasscut->Scale(1.0 / nentries, "width");
	pi_ptmasscut->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_ptmasscut1.pdf","pdf");

	nentries = pi_pmasscut->GetEntries();
	pi_pmasscut->Scale(1.0 / nentries, "width");
	pi_pmasscut->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_pmasscut1.pdf","pdf");
/*
	//Scaling???
	nentries = pi_ptvsp->GetEntries();
	pi_ptvsp->Scale(1.0 / nentries, "width");
	pi_ptvsp->Draw("COLZ");
	
	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_ptvsp1.pdf","pdf");

	//Scaling???
	nentries = pi_ptvsp_cuts->GetEntries();
	pi_ptvsp_cuts->Scale(1.0 / nentries, "width");
	pi_ptvsp_cuts->Draw("COLZ");
	
	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_ptvsp_cuts1.pdf","pdf");
*/
	pi_number_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_number_event1.pdf","pdf");

	pi_combnumber_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_combnumber_event1.pdf","pdf");

	total_pigamma_number_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_total_pigamma_number_event1.pdf","pdf");

	pigamma_combnumber_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_pigamma_combnumber_event1.pdf","pdf");

	gamma_fake->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_fake1.pdf","pdf");



	gPad->SetLogy();

	double Br2 = 0.0003100;
	double sigma = 100.304 *1e12;//Pythia event cross section converted to fb from mb
	double L = 50.0;//RUN4 luminosity
	double total = 1e5;

	nentries = eta_invmass->GetEntries();
	eta_invmass->Scale(Br2 * sigma * L/(2.0 * total), "width");
	eta_invmass->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etainvmass1.pdf","pdf");


	nentries = muon_invmass->GetEntries();
	muon_invmass->Scale(Br2 * sigma * L/(2.0 * total), "width");
	muon_invmass->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_muoninvmass1.pdf","pdf");

	nentries = muon_invmass_mother->GetEntries();
	muon_invmass_mother->Scale(Br2 * sigma * L/(2.0 * total), "width");
	muon_invmass_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_muoninvmass_mother1.pdf","pdf");

	nentries = eta_invmass_mother->GetEntries();
	eta_invmass_mother->Scale(Br2 * sigma * L/(2.0 * total), "width");
	eta_invmass_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etainvmass_mother1.pdf","pdf");

	nentries = muon_invmass_cuts_mother->GetEntries();
	muon_invmass_cuts_mother->Scale(Br2 * sigma * L/(2.0 * total), "width");
	muon_invmass_cuts_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_muoninvmass_cuts_mother1.pdf","pdf");
///
	nentries = eta_invmass_cuts_mother->GetEntries();
	eta_invmass_cuts_mother->Scale(Br2 * sigma * L/(2.0 * total), "width");
	eta_invmass_cuts_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etainvmass_cuts_mother1.pdf","pdf");

	nentries = muon_invmass_cuts_mother2->GetEntries();
	muon_invmass_cuts_mother2->Scale(Br2 * sigma * L/(2.0 * total), "width");
	muon_invmass_cuts_mother2->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_muoninvmass_cuts_mother2.pdf","pdf");
///
	nentries = eta_invmass_cuts_mother2->GetEntries();
	eta_invmass_cuts_mother2->Scale(Br2 * sigma * L/(2.0 * total), "width");
	eta_invmass_cuts_mother2->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etainvmass_cuts_mother2.pdf","pdf");

	nentries = muon_invmass_c3_mother->GetEntries();
	muon_invmass_c3_mother->Scale(Br2 * sigma * L/(2.0 * total), "width");
	muon_invmass_c3_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_muoninvmass_c3_mother1.pdf","pdf");

	nentries = eta_invmass_c3_mother->GetEntries();
	eta_invmass_c3_mother->Scale(Br2 * sigma * L/(2.0 * total), "width");
	eta_invmass_c3_mother->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etainvmass_c3_mother1.pdf","pdf");


	nentries = muon_invmass_cuts_sigbkg->GetEntries();
	muon_invmass_cuts_sigbkg->Scale(Br2 * sigma * L/(2.0 * total), "width");
	muon_invmass_cuts_sigbkg->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_muoninvmass_cuts_sigbkg1.pdf","pdf");

	nentries = eta_invmass_cuts_sigbkg->GetEntries();;
	eta_invmass_cuts_sigbkg->Scale(Br2 * sigma * L/(2.0 * total), "width");
	eta_invmass_cuts_sigbkg->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etainvmass_cuts_sigbkg1.pdf","pdf");

	nentries = muon_invmass_cuts_sigbkg2->GetEntries();
	muon_invmass_cuts_sigbkg2->Scale(Br2 * sigma * L/(2.0 * total), "width");
	muon_invmass_cuts_sigbkg2->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_muoninvmass_cuts_sigbkg2.pdf","pdf");

	nentries = eta_invmass_cuts_sigbkg2->GetEntries();;
	eta_invmass_cuts_sigbkg2->Scale(Br2 * sigma * L/(2.0 * total), "width");
	eta_invmass_cuts_sigbkg2->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etainvmass_cuts_sigbkg2.pdf","pdf");

	nentries = muon_invmass_c3_sigbkg->GetEntries();
	muon_invmass_c3_sigbkg->Scale(Br2 * sigma * L/(2.0 * total), "width");
	muon_invmass_c3_sigbkg->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_muoninvmass_c3_sigbkg1.pdf","pdf");

	nentries = eta_invmass_c3_sigbkg->GetEntries();
	eta_invmass_c3_sigbkg->Scale(Br2 * sigma * L/(2.0 * total), "width");
	eta_invmass_c3_sigbkg->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etainvmass_c3_sigbkg1.pdf","pdf");



	nentries = muon_invmass_cuts_bkg->GetEntries();
	muon_invmass_cuts_bkg->Scale(Br2 * sigma * L/(2.0 * total), "width");
	muon_invmass_cuts_bkg->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_muoninvmass_cuts_bkg1.pdf","pdf");

	nentries = eta_invmass_cuts_bkg->GetEntries();;
	eta_invmass_cuts_bkg->Scale(Br2 * sigma * L/(2.0 * total), "width");
	eta_invmass_cuts_bkg->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etainvmass_cuts_bkg1.pdf","pdf");

	nentries = muon_invmass_cuts_bkg2->GetEntries();
	muon_invmass_cuts_bkg2->Scale(Br2 * sigma * L/(2.0 * total), "width");
	muon_invmass_cuts_bkg2->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_muoninvmass_cuts_bkg2.pdf","pdf");

	nentries = eta_invmass_cuts_bkg2->GetEntries();;
	eta_invmass_cuts_bkg2->Scale(Br2 * sigma * L/(2.0 * total), "width");
	eta_invmass_cuts_bkg2->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etainvmass_cuts_bkg2.pdf","pdf");

	nentries = muon_invmass_c3_bkg->GetEntries();
	muon_invmass_c3_bkg->Scale(Br2 * sigma * L/(2.0 * total), "width");
	muon_invmass_c3_bkg->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_muoninvmass_c3_bkg1.pdf","pdf");

	nentries = eta_invmass_c3_bkg->GetEntries();
	eta_invmass_c3_bkg->Scale(Br2 * sigma * L/(2.0 * total), "width");
	eta_invmass_c3_bkg->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etainvmass_c3_bkg1.pdf","pdf");

	double misID = 1e-6;

	nentries = pigamma_invmass->GetEntries();
	pigamma_invmass->Scale(sigma * L/(2.0 * total), "width");
	pigamma_invmass->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_pigamma_invmass1.pdf","pdf");

	nentries = pigamma_invmass_c3->GetEntries();
	pigamma_invmass_c3->Scale(sigma * L/(2.0 * total), "width");
	pigamma_invmass_c3->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_pigamma_invmass_c31.pdf","pdf");

	nentries = pigamma_invmass_cuts->GetEntries();
	pigamma_invmass_cuts->Scale(sigma * L/(2.0 * total), "width");
	pigamma_invmass_cuts->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_pigamma_invmass_cuts1.pdf","pdf");

	nentries = pigamma_invmass_cuts2->GetEntries();
	pigamma_invmass_cuts2->Scale(sigma * L/(2.0 * total), "width");
	pigamma_invmass_cuts2->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_pigamma_invmass_cuts2.pdf","pdf");

	nentries = pi_invmass->GetEntries();
	pi_invmass->Scale(sigma * L/(2.0 * total), "width");
	pi_invmass->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_invmass1.pdf","pdf");

	nentries = pi_invmass_c3->GetEntries();
	pi_invmass_c3->Scale(sigma * L/(2.0 * total), "width");
	pi_invmass_c3->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_invmass_c31.pdf","pdf");

	nentries = pi_invmass_cuts->GetEntries();
	pi_invmass_cuts->Scale(sigma * L/(2.0 * total), "width");
	pi_invmass_cuts->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_invmass_cuts1.pdf","pdf");

	nentries = pi_invmass_cuts2->GetEntries();
	pi_invmass_cuts2->Scale(sigma * L/(2.0 * total), "width");
	pi_invmass_cuts2->Draw("h");

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_invmass_cuts2.pdf","pdf");


	muon_pt_bkg->SetLineColor(kYellow);
	muon_pt_bkg->SetFillStyle(1001);
	muon_pt_bkg->Draw("h");
	muon_ptcut_mother->SetLineColor(kGreen);
	muon_ptcut_mother->SetFillStyle(1001);
	muon_ptcut_mother->Draw("same h");
	muon_ptetacut_mother->SetLineColor(kRed);
	muon_ptetacut_mother->SetFillStyle(1001);
	muon_ptetacut_mother->Draw("same h");
	muon_ptetamasscut_mother->SetLineColor(kMagenta);
	muon_ptetamasscut_mother->SetFillStyle(1001);
	muon_ptetamasscut_mother->Draw("same h");
	muon_pt_mother->SetLineColor(kBlue);
	muon_pt_mother->SetFillStyle(1001);
	muon_pt_mother->Draw("same h");

	TLegend *legend3 = new TLegend(0.5,0.55,0.9,0.75);	
	TLegendEntry *leg32 = legend3->AddEntry("muon_pt_bkg","p_{t} distribution for background","f");
  	leg32->SetFillColor(kYellow);
	TLegendEntry *leg3 = legend3->AddEntry("muon_pt_mother","p_{t} distribution before cuts (sig)","f");
  	leg3->SetFillColor(kBlue);
	TLegendEntry *leg33 = legend3->AddEntry("muon_ptcut_mother","p_{t} distribution after p_{t} & p cuts (sig)","f");
  	leg33->SetFillColor(kGreen);
	TLegendEntry *leg31 = legend3->AddEntry("muon_ptetacut_mother","p_{t} distribution after p_{t}, p & #eta cuts (sig)","f");
  	leg31->SetFillColor(kRed);
	TLegendEntry *leg34 = legend3->AddEntry("muon_ptetamasscut_mother","p_{t} distribution after all cuts (sig)","f");
  	leg34->SetFillColor(kMagenta);
	legend3->Draw("same");

	muon_pt_bkg->SetTitle("#mu^{#pm} p_{t} distributions");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_pt_mother12.pdf","pdf");


	gamma_pt_mother->SetLineColor(kBlue);
	gamma_pt_mother->SetFillStyle(1001);
	gamma_pt_mother->Draw("h");
	gamma_ptetacut_mother->SetLineColor(kRed);
	gamma_ptetacut_mother->SetFillStyle(1001);
	gamma_ptetacut_mother->Draw("same h");
	gamma_ptetamasscut_mother->SetLineColor(kGreen);
	gamma_ptetamasscut_mother->SetFillStyle(1001);
	gamma_ptetamasscut_mother->Draw("same h");
	gamma_ptallcut_mother->SetLineColor(kCyan);
	gamma_ptallcut_mother->SetFillStyle(1001);
	gamma_ptallcut_mother->Draw("same h");
	

	TLegend *legend8 = new TLegend(0.3,0.7,0.7,0.9);	
	TLegendEntry *leg82 = legend8->AddEntry("gamma_pt_mother","p_{t} distribution before cuts","f");
  	leg82->SetFillColor(kBlue);
	TLegendEntry *leg84 = legend8->AddEntry("gamma_ptetacut_mother","p_{t} distribution after p_{t}, p & #eta cuts","f");
  	leg84->SetFillColor(kRed);
	TLegendEntry *leg85 = legend8->AddEntry("gamma_ptetamasscut_mother","p_{t} distribution after mass cut","f");
  	leg85->SetFillColor(kGreen);
	TLegendEntry *leg86 = legend8->AddEntry("gamma_ptallcut_mother","p_{t} distribution after all cuts including #gamma","f");
  	leg86->SetFillColor(kCyan);
	legend8->Draw("same");

	gamma_pt_mother->SetTitle("#gamma p_{t} distributions for signal");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_pt_mother12.pdf","pdf");

	gamma_pt_bkg->SetLineColor(kBlue);
	gamma_pt_bkg->SetFillStyle(1001);
	gamma_pt_bkg->Draw("h");
	gamma_pt_c3_bkg->SetLineColor(kRed);
	gamma_pt_c3_bkg->SetFillStyle(1001);
	gamma_pt_c3_bkg->Draw("same h");
	gamma_pt_c4_bkg->SetLineColor(kGreen);
	gamma_pt_c4_bkg->SetFillStyle(1001);
	gamma_pt_c4_bkg->Draw("same h");
	gamma_pt_c5_bkg->SetLineColor(kCyan);
	gamma_pt_c5_bkg->SetFillStyle(1001);
	gamma_pt_c5_bkg->Draw("same h");
	

	TLegend *legend26 = new TLegend(0.5,0.5,0.9,0.7);	
	TLegendEntry *leg262 = legend26->AddEntry("gamma_pt_bkg","p_{t} distribution before cuts","f");
  	leg262->SetFillColor(kBlue);
	TLegendEntry *leg264 = legend26->AddEntry("gamma_pt_c3_bkg","p_{t} distribution after p_{t}, p & #eta cuts","f");
  	leg264->SetFillColor(kRed);
	TLegendEntry *leg265 = legend26->AddEntry("gamma_pt_c4_bkg","p_{t} distribution after mass cut","f");
  	leg265->SetFillColor(kGreen);
	TLegendEntry *leg266 = legend26->AddEntry("gamma_pt_c5_bkg","p_{t} distribution after all cuts including #gamma","f");
  	leg266->SetFillColor(kCyan);
	legend26->Draw("same");

	gamma_pt_bkg->SetTitle("#gamma p_{t} distributions for background");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_pt_bkg12.pdf","pdf");



	pi_pt->SetLineColor(kBlue);
	pi_pt->SetFillStyle(1001);
	pi_pt->Draw("h");
	pi_ptcut->SetLineColor(kRed);
	pi_ptcut->SetFillStyle(1001);
	pi_ptcut->Draw("same h");
	pi_ptmasscut->SetLineColor(kGreen);
	pi_ptmasscut->SetFillStyle(1001);
	pi_ptmasscut->Draw("same h");

	TLegend *legend4 = new TLegend(0.3,0.7,0.7,0.9);
	TLegendEntry *leg41 = legend4->AddEntry("pi_pt","#pi^{#pm} p_{t} distribution before cuts","f");
  	leg41->SetFillColor(kBlue);	
	TLegendEntry *leg42 = legend4->AddEntry("pi_ptcut","#pi^{#pm} p_{t} distribution after p_{t}, p & #eta cuts","f");
  	leg42->SetFillColor(kRed);
	TLegendEntry *leg43 = legend4->AddEntry("pi_ptmasscut","#pi^{#pm} p_{t} distribution after all cuts","f");
  	leg43->SetFillColor(kGreen);
	legend4->Draw("same");

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_pt12.pdf","pdf");

	muon_p_bkg->SetLineColor(kYellow);
	muon_p_bkg->SetFillStyle(1001);
	muon_p_bkg->Draw("h");
	muon_pcut_mother->SetLineColor(kGreen);
	muon_pcut_mother->SetFillStyle(1001);
	muon_pcut_mother->Draw("same h");
	muon_petacut_mother->SetLineColor(kRed);
	muon_petacut_mother->SetFillStyle(1001);
	muon_petacut_mother->Draw("same h");
	muon_petamasscut_mother->SetLineColor(kMagenta);
	muon_petamasscut_mother->SetFillStyle(1001);
	muon_petamasscut_mother->Draw("same h");
	muon_p_mother->SetLineColor(kBlue);
	muon_p_mother->SetFillStyle(1001);
	muon_p_mother->Draw("same h");

	TLegend *legend6 = new TLegend(0.3,0.7,0.7,0.9);
	TLegendEntry *leg5 = legend6->AddEntry("muon_p_mother","p distribution before cuts (sig)","f");
  	leg3->SetFillColor(kBlue);	
	TLegendEntry *leg62 = legend6->AddEntry("muon_p_bkg","p distribution for background","f");
  	leg62->SetFillColor(kYellow);
	TLegendEntry *leg6 = legend6->AddEntry("muon_pcut_mother","p distribution p_{t} & p cuts","f");
  	leg6->SetFillColor(kGreen);
	TLegendEntry *leg61 = legend6->AddEntry("muon_petacut_mother","p distribution after p_{t}, p & #eta cuts","f");
  	leg61->SetFillColor(kRed);
	TLegendEntry *leg63 = legend6->AddEntry("muon_petamasscut_mother","p distribution after all cuts","f");
  	leg63->SetFillColor(kMagenta);
	legend6->Draw("same");

	muon_p_bkg->SetTitle("#mu^{#pm} p distributions");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_p_mother12.pdf","pdf");


	gamma_p_mother->SetLineColor(kBlue);
	gamma_p_mother->SetFillStyle(1001);
	gamma_p_mother->Draw("h");	
	gamma_petacut_mother->SetLineColor(kRed);
	gamma_petacut_mother->SetFillStyle(1001);
	gamma_petacut_mother->Draw("same h");
	gamma_petamasscut_mother->SetLineColor(kGreen);
	gamma_petamasscut_mother->SetFillStyle(1001);
	gamma_petamasscut_mother->Draw("same h");
	gamma_pallcut_mother->SetLineColor(kCyan);
	gamma_pallcut_mother->SetFillStyle(1001);
	gamma_pallcut_mother->Draw("same h");
	

	TLegend *legend9 = new TLegend(0.3,0.7,0.7,0.9);	
	TLegendEntry *leg92 = legend9->AddEntry("gamma_p_mother","p distribution before cuts","f");
  	leg92->SetFillColor(kBlue);
	TLegendEntry *leg94 = legend9->AddEntry("gamma_petacut_mother","p distribution after p_{t}, p & #eta cuts","f");
  	leg94->SetFillColor(kRed);
	TLegendEntry *leg95 = legend9->AddEntry("gamma_petamasscut_mother","p distribution after mass cut","f");
  	leg95->SetFillColor(kGreen);
	TLegendEntry *leg96 = legend9->AddEntry("gamma_pallcut_mother","p distribution after all cuts including #gamma","f");
  	leg96->SetFillColor(kCyan);
	legend9->Draw("same");

	gamma_p_mother->SetTitle("#gamma p distributions for signal");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_p_mother12.pdf","pdf");


	gamma_p_bkg->SetLineColor(kBlue);
	gamma_p_bkg->SetFillStyle(1001);
	gamma_p_bkg->Draw("h");
	gamma_p_c3_bkg->SetLineColor(kRed);
	gamma_p_c3_bkg->SetFillStyle(1001);
	gamma_p_c3_bkg->Draw("same h");
	gamma_p_c4_bkg->SetLineColor(kGreen);
	gamma_p_c4_bkg->SetFillStyle(1001);
	gamma_p_c4_bkg->Draw("same h");
	gamma_p_c5_bkg->SetLineColor(kCyan);
	gamma_p_c5_bkg->SetFillStyle(1001);
	gamma_p_c5_bkg->Draw("same h");
	

	TLegend *legend25 = new TLegend(0.3,0.7,0.7,0.9);	
	TLegendEntry *leg252 = legend25->AddEntry("gamma_p_bkg","p distribution before cuts ","f");
  	leg252->SetFillColor(kBlue);
	TLegendEntry *leg254 = legend25->AddEntry("gamma_p_c3_bkg","p distribution after p_{t}, p & #eta cuts ","f");
  	leg254->SetFillColor(kRed);
	TLegendEntry *leg255 = legend25->AddEntry("gamma_p_c4_bkg","p distribution after mass cut","f");
  	leg255->SetFillColor(kGreen);
	TLegendEntry *leg256 = legend25->AddEntry("gamma_p_c5_bkg","p distribution after all cuts including #gamma","f");
  	leg256->SetFillColor(kCyan);
	legend25->Draw("same");

	gamma_p_bkg->SetTitle("#gamma p distributions for background");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_p_bkg12.pdf","pdf");



	pi_p->SetLineColor(kBlue);
	pi_p->SetFillStyle(1001);
	pi_p->Draw("h");
	pi_pcut->SetLineColor(kRed);
	pi_pcut->SetFillStyle(1001);
	pi_pcut->Draw("same h");
	pi_pmasscut->SetLineColor(kGreen);
	pi_pmasscut->SetFillStyle(1001);
	pi_pmasscut->Draw("same h");

	TLegend *legend5 = new TLegend(0.3,0.7,0.7,0.9);
	TLegendEntry *leg51 = legend5->AddEntry("pi_p","#pi^{#pm} p distribution before cuts","f");
  	leg51->SetFillColor(kBlue);	
	TLegendEntry *leg52 = legend5->AddEntry("pi_pcut","#pi^{#pm} p distribution after cuts","f");
  	leg52->SetFillColor(kRed);
	TLegendEntry *leg53 = legend5->AddEntry("pi_pmasscut","#pi^{#pm} p distribution after cuts including mass","f");
  	leg53->SetFillColor(kGreen);
	legend5->Draw("same");

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_p12.pdf","pdf");


	

	eta_invmass->SetLineColor(kBlue);
	eta_invmass->SetFillStyle(1001);
	eta_invmass->Draw("h");
	eta_invmass_mother->SetLineColor(kGreen);
	eta_invmass_mother->SetFillStyle(1001);
	eta_invmass_mother->Draw("same h");

	TLegend *legend1 = new TLegend(0.7,0.5,0.9,0.7);
	TLegendEntry *leg1 = legend1->AddEntry("eta_invmass","#eta invariant mass (sig+bkg)","f");
  	leg1->SetFillColor(kBlue);	
	TLegendEntry *leg2 = legend1->AddEntry("eta_invmass_mother","#eta invariant mass (sig)","f");
  	leg2->SetFillColor(kGreen);
	legend1->Draw("same");

	//gPad->SetLogx();

	c1->Modified();
	c1->Update();
	c1->Print("project54_etainvmass12.pdf","pdf");


	eta_invmass->SetLineColor(kBlue);
	eta_invmass->SetFillStyle(1001);
	eta_invmass->Draw("h");
	eta_invmass_mother->SetLineColor(kGreen);
	eta_invmass_mother->SetFillStyle(1001);
	eta_invmass_mother->Draw("same h");
	pigamma_invmass->SetLineColor(kOrange);
	pigamma_invmass->SetFillStyle(1001);
	pigamma_invmass->Draw("same h");

	TLegend *legend7 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg70 = legend7->AddEntry("eta_invmass","#eta invariant mass (sig+bkg)","f");
  	leg70->SetFillColor(kBlue);	
	TLegendEntry *leg71 = legend7->AddEntry("eta_invmass_mother","#eta invariant mass (sig)","f");
  	leg71->SetFillColor(kGreen);
	TLegendEntry *leg72 = legend7->AddEntry("pigamma_invmass","misID 2#pi+#gamma invariant mass","f");
  	leg71->SetFillColor(kOrange);
	legend7->Draw("same");

	//eta_invmass->SetAxisRange(1e6, 1e12,"Y");
	eta_invmass->SetTitle("Reconstructed #eta invariant mass");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etapigamma_invmass12.pdf","pdf");


	eta_invmass_cuts_mother->SetLineColor(kBlue);
	eta_invmass_cuts_mother->SetFillStyle(1001);
	eta_invmass_cuts_mother->Draw("h");
	eta_invmass_cuts_mother2->SetLineColor(kRed);
	eta_invmass_cuts_mother2->SetFillStyle(1001);
	eta_invmass_cuts_mother2->Draw("same h");
	pigamma_invmass_cuts->SetLineColor(kGreen);
	pigamma_invmass_cuts->SetFillStyle(1001);
	pigamma_invmass_cuts->Draw("same h");
	pigamma_invmass_cuts2->SetLineColor(kOrange);
	pigamma_invmass_cuts2->SetFillStyle(1001);
	pigamma_invmass_cuts2->Draw("same h");

	TLegend *legend17 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg170 = legend17->AddEntry("eta_invmass_cuts_mother","#eta invariant mass (sig)","f");
  	leg170->SetFillColor(kBlue);
	TLegendEntry *leg172 = legend17->AddEntry("eta_invmass_cuts_mother2","#eta invariant mass including #gamma (sig)","f");
  	leg172->SetFillColor(kRed);	
	TLegendEntry *leg171 = legend17->AddEntry("pigamma_invmass_cuts","misID 2#pi+#gamma invariant mass","f");
  	leg171->SetFillColor(kGreen);
	TLegendEntry *leg173 = legend17->AddEntry("pigamma_invmass_cuts2","misID 2#pi+#gamma invariant mass including #gamma","f");
  	leg173->SetFillColor(kOrange);
	legend17->Draw("same");

	//eta_invmass_cuts_mother->SetAxisRange(1e-6, 1e2,"Y");
	eta_invmass_cuts_mother->SetTitle("Reconstructed #eta invariant mass after all cuts");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etapigamma_invmass_cuts12.pdf","pdf");


	eta_invmass_c3_mother->SetLineColor(kBlue);
	eta_invmass_c3_mother->SetFillStyle(1001);
	eta_invmass_c3_mother->Draw("h");
	pigamma_invmass_c3->SetLineColor(kRed);
	pigamma_invmass_c3->SetFillStyle(1001);
	pigamma_invmass_c3->Draw("same h");

	TLegend *legend18 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg180 = legend18->AddEntry("eta_invmass_c3_mother","#eta invariant mass (sig)","f");
  	leg180->SetFillColor(kBlue);	
	TLegendEntry *leg181 = legend18->AddEntry("pigamma_invmass_c3","misID 2#pi+#gamma invariant mass","f");
  	leg181->SetFillColor(kRed);
	legend18->Draw("same");

	//eta_invmass_c3_mother->SetAxisRange(1e1, 1e7,"Y");
	eta_invmass_c3_mother->SetTitle("Reconstructed #eta invariant mass after p_{t}, p and #eta cuts");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etapigamma_invmass_c312.pdf","pdf");




	muon_invmass->SetLineColor(kBlue);
	muon_invmass->SetFillStyle(1001);
	muon_invmass->Draw("h");
	muon_invmass_mother->SetLineColor(kGreen);
	muon_invmass_mother->SetFillStyle(1001);
	muon_invmass_mother->Draw("same h");
	pi_invmass->SetLineColor(kOrange);
	pi_invmass->SetFillStyle(1001);
	pi_invmass->Draw("same h");

	TLegend *legend10 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg100 = legend10->AddEntry("muon_invmass","di-muon invariant mass (sig+bkg)","f");
  	leg100->SetFillColor(kBlue);	
	TLegendEntry *leg101 = legend10->AddEntry("muon_invmass_mother","di-muon invariant mass (sig)","f");
  	leg101->SetFillColor(kGreen);
	TLegendEntry *leg102 = legend10->AddEntry("pi_invmass","misID #pi invariant mass","f");
  	leg101->SetFillColor(kOrange);
	legend10->Draw("same");

	//muon_invmass->SetAxisRange(1e4, 1e10,"Y");
	muon_invmass->SetAxisRange(0, 2,"X");
	muon_invmass->SetTitle("Reconstructed di-muon invariant mass");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mupi_invmass12.pdf","pdf");


	muon_invmass_cuts_mother->SetLineColor(kBlue);
	muon_invmass_cuts_mother->SetFillStyle(1001);
	muon_invmass_cuts_mother->Draw("h");
	muon_invmass_cuts_mother2->SetLineColor(kRed);
	muon_invmass_cuts_mother2->SetFillStyle(1001);
	muon_invmass_cuts_mother2->Draw("same h");
	pi_invmass_cuts->SetLineColor(kGreen);
	pi_invmass_cuts->SetFillStyle(1001);
	pi_invmass_cuts->Draw("same h");
	pi_invmass_cuts2->SetLineColor(kOrange);
	pi_invmass_cuts2->SetFillStyle(1001);
	pi_invmass_cuts2->Draw("same h");

	TLegend *legend11 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg110 = legend11->AddEntry("muon_invmass_cuts_mother","di-muon invariant mass (sig)","f");
  	leg110->SetFillColor(kBlue);
	TLegendEntry *leg112 = legend11->AddEntry("muon_invmass_cuts_mother2","di-muon invariant mass including #gamma (sig)","f");
  	leg112->SetFillColor(kRed);	
	TLegendEntry *leg111 = legend11->AddEntry("pi_invmass_cuts","misID #pi invariant mass","f");
  	leg111->SetFillColor(kGreen);
	TLegendEntry *leg113 = legend11->AddEntry("pi_invmass_cuts2","misID #pi invariant mass including #gamma","f");
  	leg113->SetFillColor(kOrange);
	legend11->Draw("same");

	//muon_invmass_cuts_mother->SetAxisRange(1e-5, 1e1,"Y");
	muon_invmass_cuts_mother->SetTitle("Reconstructed di-muon invariant mass after all cuts (sig)");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mupi_invmass_cuts12.pdf","pdf");


	muon_invmass_c3_mother->SetLineColor(kBlue);
	muon_invmass_c3_mother->SetFillStyle(1001);
	muon_invmass_c3_mother->Draw("h");
	pi_invmass_c3->SetLineColor(kRed);
	pi_invmass_c3->SetFillStyle(1001);
	pi_invmass_c3->Draw("same h");

	TLegend *legend12 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg120 = legend12->AddEntry("muon_invmass_c3_mother","di-muon invariant mass (sig)","f");
  	leg120->SetFillColor(kBlue);	
	TLegendEntry *leg121 = legend12->AddEntry("pi_invmass_c3","misID #pi invariant mass","f");
  	leg121->SetFillColor(kRed);
	legend12->Draw("same");

	//muon_invmass_c3_mother->SetAxisRange(1e1, 1e7,"Y");
	muon_invmass_c3_mother->SetTitle("Reconstructed di-muon invariant mass after p_{t}, p and #eta cuts");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mupi_invmass_c312.pdf","pdf");


/*
	eta_invmass_cuts_sigbkg->SetLineColor(kBlue);
	eta_invmass_cuts_sigbkg->SetFillStyle(1001);
	eta_invmass_cuts_sigbkg->Draw("h");
	eta_invmass_cuts_sigbkg2->SetLineColor(kRed);
	eta_invmass_cuts_sigbkg2->SetFillStyle(1001);
	eta_invmass_cuts_sigbkg2->Draw("same h");
	pigamma_invmass_cuts->SetLineColor(kGreen);
	pigamma_invmass_cuts->SetFillStyle(1001);
	pigamma_invmass_cuts->Draw("same h");
	pigamma_invmass_cuts2->SetLineColor(kOrange);
	pigamma_invmass_cuts2->SetFillStyle(1001);
	pigamma_invmass_cuts2->Draw("same h");

	TLegend *legend13 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg130 = legend13->AddEntry("eta_invmass_cuts_sigbkg","#eta invariant mass (sig+bkg)","f");
  	leg130->SetFillColor(kBlue);
	TLegendEntry *leg132 = legend13->AddEntry("eta_invmass_cuts_sigbkg2","#eta invariant mass including #gamma (sig+bkg)","f");
  	leg132->SetFillColor(kRed);	
	TLegendEntry *leg131 = legend13->AddEntry("pigamma_invmass_cuts","misID 2#pi+#gamma invariant mass","f");
  	leg131->SetFillColor(kGreen);
	TLegendEntry *leg133 = legend13->AddEntry("pigamma_invmass_cuts2","misID 2#pi+#gamma invariant mass including #gamma","f");
  	leg133->SetFillColor(kOrange);
	legend13->Draw("same");

	eta_invmass_cuts_sigbkg->SetAxisRange(1e-5, 1e2,"Y");
	eta_invmass_cuts_sigbkg->SetTitle("Reconstructed #eta invariant mass after all cuts (sig+bkg)");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etapigamma_invmass_cuts_sigbkg12.pdf","pdf");



	eta_invmass_cuts_bkg->SetLineColor(kBlue);
	eta_invmass_cuts_bkg->SetFillStyle(1001);
	eta_invmass_cuts_bkg->Draw("h");
	eta_invmass_cuts_bkg2->SetLineColor(kRed);
	eta_invmass_cuts_bkg2->SetFillStyle(1001);
	eta_invmass_cuts_bkg2->Draw("same h");
	pigamma_invmass_cuts->SetLineColor(kGreen);
	pigamma_invmass_cuts->SetFillStyle(1001);
	pigamma_invmass_cuts->Draw("same h");
	pigamma_invmass_cuts2->SetLineColor(kOrange);
	pigamma_invmass_cuts2->SetFillStyle(1001);
	pigamma_invmass_cuts2->Draw("same h");

	TLegend *legend23 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg230 = legend23->AddEntry("eta_invmass_cuts_bkg","#eta invariant mass (bkg)","f");
  	leg230->SetFillColor(kBlue);
	TLegendEntry *leg232 = legend23->AddEntry("eta_invmass_cuts_bkg2","#eta invariant mass including #gamma (bkg)","f");
  	leg232->SetFillColor(kRed);	
	TLegendEntry *leg231 = legend23->AddEntry("pigamma_invmass_cuts","misID 2#pi+#gamma invariant mass","f");
  	leg231->SetFillColor(kGreen);
	TLegendEntry *leg233 = legend23->AddEntry("pigamma_invmass_cuts2","misID 2#pi+#gamma invariant mass including #gamma","f");
  	leg233->SetFillColor(kOrange);
	legend23->Draw("same");

	eta_invmass_cuts_bkg->SetAxisRange(1e-5, 1e2,"Y");
	eta_invmass_cuts_bkg->SetTitle("Reconstructed #eta invariant mass after all cuts (bkg)");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etapigamma_invmass_cuts_bkg12.pdf","pdf");



	eta_invmass_c3_sigbkg->SetLineColor(kBlue);
	eta_invmass_c3_sigbkg->SetFillStyle(1001);
	eta_invmass_c3_sigbkg->Draw("h");
	pigamma_invmass_c3->SetLineColor(kRed);
	pigamma_invmass_c3->SetFillStyle(1001);
	pigamma_invmass_c3->Draw("same h");

	TLegend *legend14 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg140 = legend14->AddEntry("eta_invmass_c3_sigbkg","#eta invariant mass (sig+bkg)","f");
  	leg140->SetFillColor(kBlue);	
	TLegendEntry *leg141 = legend14->AddEntry("pigamma_invmass_c3","misID 2#pi+#gamma invariant mass","f");
  	leg141->SetFillColor(kRed);
	legend14->Draw("same");

	//eta_invmass_c3_bkg->SetAxisRange(1e6, 1e12,"Y");
	eta_invmass_c3_sigbkg->SetTitle("Reconstructed #eta invariant mass after p_{t}, p and #eta cuts (sig+bkg)");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etapigamma_invmass_c3_sigbkg12.pdf","pdf");*/


	muon_invmass_cuts_sigbkg->SetLineColor(kBlue);
	muon_invmass_cuts_sigbkg->SetFillStyle(1001);
	muon_invmass_cuts_sigbkg->Draw("h");
	muon_invmass_cuts_sigbkg2->SetLineColor(kRed);
	muon_invmass_cuts_sigbkg2->SetFillStyle(1001);
	muon_invmass_cuts_sigbkg2->Draw("same h");
	pi_invmass_cuts->SetLineColor(kGreen);
	pi_invmass_cuts->SetFillStyle(1001);
	pi_invmass_cuts->Draw("same h");
	pi_invmass_cuts2->SetLineColor(kOrange);
	pi_invmass_cuts2->SetFillStyle(1001);
	pi_invmass_cuts2->Draw("same h");

	TLegend *legend15 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg150 = legend15->AddEntry("muon_invmass_cuts_sigbkg","di-muon invariant mass (sig+bkg)","f");
  	leg150->SetFillColor(kBlue);	
	TLegendEntry *leg152 = legend15->AddEntry("muon_invmass_cuts_sigbkg2","di-muon invariant mass including #gamma (sig+bkg)","f");
  	leg152->SetFillColor(kRed);	
	TLegendEntry *leg151 = legend15->AddEntry("pi_invmass_cuts","misID #pi invariant mass","f");
  	leg151->SetFillColor(kGreen);
	TLegendEntry *leg153 = legend15->AddEntry("pi_invmass_cuts2","misID #pi invariant mass including #gamma","f");
  	leg153->SetFillColor(kOrange);
	legend15->Draw("same");

	//muon_invmass_cuts_sigbkg->SetAxisRange(1e-5, 1e1,"Y");
	muon_invmass_cuts_sigbkg->SetTitle("Reconstructed di-muon invariant mass after all cuts (sig+bkg)");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mupi_invmass_cuts_sigbkg12.pdf","pdf");


	muon_invmass_cuts_bkg->SetLineColor(kBlue);
	muon_invmass_cuts_bkg->SetFillStyle(1001);
	muon_invmass_cuts_bkg->Draw("h");
	muon_invmass_cuts_bkg2->SetLineColor(kRed);
	muon_invmass_cuts_bkg2->SetFillStyle(1001);
	muon_invmass_cuts_bkg2->Draw("same h");
	pi_invmass_cuts->SetLineColor(kGreen);
	pi_invmass_cuts->SetFillStyle(1001);
	pi_invmass_cuts->Draw("same h");
	pi_invmass_cuts2->SetLineColor(kOrange);
	pi_invmass_cuts2->SetFillStyle(1001);
	pi_invmass_cuts2->Draw("same h");

	TLegend *legend27 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg270 = legend27->AddEntry("muon_invmass_cuts_bkg","di-muon invariant mass (bkg)","f");
  	leg270->SetFillColor(kBlue);	
	TLegendEntry *leg272 = legend27->AddEntry("muon_invmass_cuts_bkg2","di-muon invariant mass including #gamma (bkg)","f");
  	leg272->SetFillColor(kRed);	
	TLegendEntry *leg271 = legend27->AddEntry("pi_invmass_cuts","misID #pi invariant mass","f");
  	leg271->SetFillColor(kGreen);
	TLegendEntry *leg273 = legend27->AddEntry("pi_invmass_cuts2","misID #pi invariant mass including #gamma","f");
  	leg273->SetFillColor(kOrange);
	legend27->Draw("same");

	//muon_invmass_cuts_bkg->SetAxisRange(1e-5, 1e1,"Y");
	muon_invmass_cuts_bkg->SetTitle("Reconstructed di-muon invariant mass after all cuts (bkg)");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mupi_invmass_cuts_bkg123.pdf","pdf");


	muon_invmass_c3_sigbkg->SetLineColor(kBlue);
	muon_invmass_c3_sigbkg->SetFillStyle(1001);
	muon_invmass_c3_sigbkg->Draw("h");
	pi_invmass_c3->SetLineColor(kRed);
	pi_invmass_c3->SetFillStyle(1001);
	pi_invmass_c3->Draw("same h");

	TLegend *legend16 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg160 = legend16->AddEntry("muon_invmass_c3_sigbkg","di-muon invariant mass (sig+bkg)","f");
  	leg160->SetFillColor(kBlue);	
	TLegendEntry *leg161 = legend16->AddEntry("pi_invmass_c3","misID #pi invariant mass","f");
  	leg161->SetFillColor(kRed);
	legend16->Draw("same");

	//muon_invmass_c3_bkg->SetAxisRange(1e6, 1e12,"Y");
	muon_invmass_c3_sigbkg->SetTitle("Reconstructed di-muon invariant mass after p_{t}, p and #eta cuts (sig+bkg)");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mupi_invmass_c3_sig12.pdf","pdf");



	muon_invmass_c3_mother->SetLineColor(kBlue);
	muon_invmass_c3_mother->SetFillStyle(1001);
	muon_invmass_c3_mother->Draw("h");
	muon_invmass_c3_sigbkg->SetLineColor(kGreen);
	muon_invmass_c3_sigbkg->SetFillStyle(1001);
	muon_invmass_c3_sigbkg->Draw("same h");
	pi_invmass_c3->SetLineColor(kRed);
	pi_invmass_c3->SetFillStyle(1001);
	pi_invmass_c3->Draw("same h");

	TLegend *legend20 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg200 = legend20->AddEntry("muon_invmass_c3_mother","di-muon invariant mass (sig)","f");
  	leg200->SetFillColor(kBlue);	
	TLegendEntry *leg201 = legend20->AddEntry("muon_invmass_c3_sigbkg","di-muon invariant mass (sig+bkg)","f");
  	leg201->SetFillColor(kGreen);	
	TLegendEntry *leg202 = legend20->AddEntry("pi_invmass_c3","misID #pi invariant mass","f");
  	leg202->SetFillColor(kRed);
	legend20->Draw("same");

	//muon_invmass_c3_mother->SetAxisRange(1e-5, 1e3,"Y");
	muon_invmass_c3_mother->SetTitle("Reconstructed di-muon invariant mass after p_{t}, p and #eta cuts (sig & sig+bkg)");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mupi_invmass_c3_sig123.pdf","pdf");


	muon_invmass_cuts_mother->SetLineColor(kBlue);
	muon_invmass_cuts_mother->SetFillStyle(1001);
	muon_invmass_cuts_mother->Draw("h");
	muon_invmass_cuts_sigbkg->SetLineColor(kGreen);
	muon_invmass_cuts_sigbkg->SetFillStyle(1001);
	muon_invmass_cuts_sigbkg->Draw("same h");
	pi_invmass_cuts->SetLineColor(kRed);
	pi_invmass_cuts->SetFillStyle(1001);
	pi_invmass_cuts->Draw("same h");

	TLegend *legend21 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg210 = legend21->AddEntry("muon_invmass_cuts_mother","di-muon invariant mass (sig)","f");
  	leg210->SetFillColor(kBlue);	
	TLegendEntry *leg211 = legend21->AddEntry("muon_invmass_cuts_sigbkg","di-muon invariant mass (sig+bkg)","f");
  	leg211->SetFillColor(kGreen);	
	TLegendEntry *leg212 = legend21->AddEntry("pi_invmass_cuts","misID #pi invariant mass","f");
  	leg212->SetFillColor(kRed);
	legend21->Draw("same");

	//muon_invmass_cuts_mother->SetAxisRange(1e-6, 1e1,"Y");
	muon_invmass_cuts_mother->SetTitle("Reconstructed di-muon invariant mass after all cuts (sig & sig+bkg)");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mupi_invmass_cuts_sig12.pdf","pdf");


	muon_invmass_cuts_mother->SetLineColor(kBlue);
	muon_invmass_cuts_mother->SetFillStyle(1001);
	muon_invmass_cuts_mother->Draw("h");
	muon_invmass_cuts_bkg->SetLineColor(kGreen);
	muon_invmass_cuts_bkg->SetFillStyle(1001);
	muon_invmass_cuts_bkg->Draw("same h");
	pi_invmass_cuts->SetLineColor(kRed);
	pi_invmass_cuts->SetFillStyle(1001);
	pi_invmass_cuts->Draw("same h");

	TLegend *legend24 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg240 = legend24->AddEntry("muon_invmass_cuts_mother","di-muon invariant mass (sig)","f");
  	leg240->SetFillColor(kBlue);	
	TLegendEntry *leg241 = legend24->AddEntry("muon_invmass_cuts_bkg","di-muon invariant mass (bkg)","f");
  	leg241->SetFillColor(kGreen);	
	TLegendEntry *leg242 = legend24->AddEntry("pi_invmass_cuts","misID #pi invariant mass","f");
  	leg242->SetFillColor(kRed);
	legend24->Draw("same");

	//muon_invmass_cuts_mother->SetAxisRange(1e-6, 1e1,"Y");
	muon_invmass_cuts_mother->SetTitle("Reconstructed di-muon invariant mass after all cuts (sig & bkg)");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mupi_invmass_cuts_bkg12.pdf","pdf");


	muon_invmass_cuts_sigbkg2->SetLineColor(kGreen);
	muon_invmass_cuts_sigbkg2->SetFillStyle(1001);
	muon_invmass_cuts_sigbkg2->Draw("h");
	muon_invmass_cuts_mother2->SetLineColor(kBlue);
	muon_invmass_cuts_mother2->SetFillStyle(1001);
	muon_invmass_cuts_mother2->Draw("same h");
	pi_invmass_cuts2->SetLineColor(kRed);
	pi_invmass_cuts2->SetFillStyle(1001);
	pi_invmass_cuts2->Draw("same h");

	TLegend *legend22 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg220 = legend22->AddEntry("muon_invmass_cuts_mother2","di-muon invariant mass (sig)","f");
  	leg220->SetFillColor(kBlue);	
	TLegendEntry *leg221 = legend22->AddEntry("muon_invmass_cuts_sigbkg2","di-muon invariant mass (sig+bkg)","f");
  	leg221->SetFillColor(kGreen);	
	TLegendEntry *leg222 = legend22->AddEntry("pi_invmass_cuts2","misID #pi invariant mass","f");
  	leg222->SetFillColor(kRed);
	legend22->Draw("same");

	//muon_invmass_cuts_sigbkg2->SetAxisRange(1e-6, 1e1,"Y");
	muon_invmass_cuts_sigbkg2->SetTitle("Reconstructed di-muon invariant mass after all cuts including #gamma (sig & sig+bkg)");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mupi_invmass_cuts2_sigbkg12.pdf","pdf");


	muon_invmass_cuts_bkg2->SetLineColor(kGreen);
	muon_invmass_cuts_bkg2->SetFillStyle(1001);
	muon_invmass_cuts_bkg2->Draw("h");
	muon_invmass_cuts_mother2->SetLineColor(kBlue);
	muon_invmass_cuts_mother2->SetFillStyle(1001);
	muon_invmass_cuts_mother2->Draw("same h");
	pi_invmass_cuts2->SetLineColor(kRed);
	pi_invmass_cuts2->SetFillStyle(1001);
	pi_invmass_cuts2->Draw("same h");

	TLegend *legend28 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg280 = legend28->AddEntry("muon_invmass_cuts_mother2","di-muon invariant mass (sig)","f");
  	leg280->SetFillColor(kBlue);	
	TLegendEntry *leg281 = legend28->AddEntry("muon_invmass_cuts_bkg2","di-muon invariant mass (bkg)","f");
  	leg281->SetFillColor(kGreen);	
	TLegendEntry *leg282 = legend28->AddEntry("pi_invmass_cuts2","misID #pi invariant mass","f");
  	leg282->SetFillColor(kRed);
	legend28->Draw("same");

	//muon_invmass_cuts_bkg2->SetAxisRange(1e-6, 1e1,"Y");
	muon_invmass_cuts_bkg2->SetTitle("Reconstructed di-muon invariant mass after all cuts including #gamma (sig & bkg)");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mupi_invmass_cuts2_bkg12.pdf","pdf");



	std::cout<<"Parent ID of background photons that pass the cuts: "<<std::endl;
	for(int i=1;i<=2000;i++){

		if(gamma_fake->GetBinContent(i)!=0.0){std::cout<<gamma_fake->GetXaxis()->GetBinLowEdge(i)<<std::endl;}

	}


	std::cout<<"Total signal entries: "<<vn->n[0]<<std::endl;
	std::cout<<"After pt cuts: "<<vn->n[1]<<std::endl;
	std::cout<<"After pt & p cuts: "<<vn->n[2]<<std::endl;
	std::cout<<"After pt, p & pseudorapidity cuts: "<<vn->n[3]<<std::endl;
	std::cout<<"After mass cut: "<<vn->n[4]<<std::endl;
	std::cout<<"After all cuts (including photons): "<<vn->n[5]<<std::endl;

	std::cout<<""<<std::endl;

	std::cout<<"Total misID pion entries: "<<vp->n[0]<<std::endl;
	std::cout<<"After pt, p & pseudorapidity cuts: "<<vp->n[1]<<std::endl;
	std::cout<<"After mass cut: "<<vp->n[2]<<std::endl;
	std::cout<<"After all cuts (including photons): "<<vp->n[3]<<std::endl;

	std::cout<<""<<std::endl;

	std::cout<<"Signal significance of signal entries: "<<std::endl;
	std::cout<<"After pt cuts: "<<vn->n[1]/vn->n[0]<<std::endl;
	std::cout<<"After pt & p cuts: "<<vn->n[2]/vn->n[0]<<std::endl;
	std::cout<<"After pt, p & pseudorapidity cuts: "<<vn->n[3]/vn->n[0]<<std::endl;
	std::cout<<"After mass cut: "<<vn->n[4]/vn->n[0]<<std::endl;
	std::cout<<"After all cuts (including photons): "<<vn->n[5]/vn->n[0]<<std::endl;

	std::cout<<""<<std::endl;
/*
	std::cout<<"Signal significance of misID pion entries: "<<std::endl;
	std::cout<<"After pt, p & pseudorapidity cuts: "<<vp->n[1]/vp->n[0]<<std::endl;
	std::cout<<"After mass cut: "<<vp->n[2]/vp->n[0]<<std::endl;
	std::cout<<"After all cuts (including photons): "<<vp->n[3]/vp->n[0]<<std::endl;*/

	FILE  *file_out3 = fopen("signal_background_out54.txt", "w");

	fprintf(file_out3, "Total signal entries:  %f \n", vn->n[0]);
	fprintf(file_out3, "After pt cuts: %f \n", vn->n[1]);
	fprintf(file_out3, "After pt & p cuts:  %f \n", vn->n[2]);
	fprintf(file_out3, "After pt, p & pseudorapidity cuts:  %f \n", vn->n[3]);
	fprintf(file_out3, "After mass cut:  %f \n ", vn->n[4]);
	fprintf(file_out3, "After all cuts (including photons):  %f \n \n", vn->n[5]);

	fprintf(file_out3, "Total misID pion entries:  %f \n", vp->n[0]);
	fprintf(file_out3, "After pt, p & pseudorapidity cuts:  %f \n", vp->n[1]);
	fprintf(file_out3, "After mass cut:  %f \n", vp->n[2]);
	fprintf(file_out3, "After all cuts (including photons):  %f \n \n", vp->n[3]);

	fprintf(file_out3, "Signal significance of signal entries: \n");
	fprintf(file_out3, "After pt cuts: %e \n", vn->n[1]/vn->n[0]);
	fprintf(file_out3, "After pt & p cuts:  %e \n", vn->n[2]/vn->n[0]);
	fprintf(file_out3, "After pt, p & pseudorapidity cuts:  %e \n", vn->n[3]/vn->n[0]);
	fprintf(file_out3, "After mass cut:  %e \n ", vn->n[4]/vn->n[0]);
	fprintf(file_out3, "After all cuts (including photons):  %e \n \n", vn->n[5]/vn->n[0]);
/*
	fprintf(file_out3, "Signal significance of misID pion entries: \n");
	fprintf(file_out3, "After pt, p & pseudorapidity cuts:  %e \n", vp->n[1]/vp->n[0]);
	fprintf(file_out3, "After mass cut:  %e \n", vp->n[2]/vp->n[0]);
	fprintf(file_out3, "After all cuts (including photons):  %e \n", vp->n[3]/vp->n[0]);*/

//std::cout<<"here1"<<std::endl;

	file_out->Write();
//std::cout<<"here1a"<<std::endl;
	file_out->Close();	
//std::cout<<"here1b"<<std::endl;
//std::cout<<"here2"<<std::endl;
        //file_in->Close();
	//delete file_out;
//std::cout<<"here3"<<std::endl;
	fclose(file_out1);
	fclose(file_out2);
	fclose(file_out3);
	file.close();
	delete v;
	delete vn;
	delete vp;
	// Done.
	return 0;
}


