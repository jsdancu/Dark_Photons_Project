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

void analyze_event(FILE *file_out1, vect *v, TH1D *eta_invmass, TH1D *muon_invmass, TH1D *muongamma_invmass_mother, TH1D *muon_invmass_mother, TH1D *eta_invmass_mothers, TH1D *muon_invmass_cuts_mother, TH1D *eta_invmass_cuts_mother, TH1D *muon_invmass_c3_mother, TH1D *eta_invmass_c3_mother, TH1D *muon_invmass_cuts_bkg, TH1D *eta_invmass_cuts_bkg, TH1D *muon_invmass_c3_bkg, TH1D *eta_invmass_c3_bkg, TH1D *muon_pt_mother, TH1D *muon_pt_bkg, TH1D *muon_ptcut_mother, TH1D *muon_ptetacut_mother, TH1D *muon_ptetamasscut_mother, TH1D *muon_p_mother, TH1D *muon_p_bkg, TH1D *muon_pcut_mother, TH1D *muon_petacut_mother, TH1D *muon_petamasscut_mother, TH2D *muon_ptvsp_mother, TH2D *muon_ptvsp_cuts_mother, TH2D *muon_ptvsp_c1_mother, TH2D *muon_ptvsp_c2_mother, TH2D *muon_ptvsp_c3_mother, TH2D *gamma_ptvsp_mother, TH2D *gamma_ptvsp_cuts_mother, TH2D *gamma_ptvsp_c1_mother, TH2D *gamma_ptvsp_c2_mother, TH2D *gamma_ptvsp_c3_mother, TH2D *muon_ptvseta_mother, TH2D *muon_ptvseta_cuts_mother, TH2D *muon_ptvseta_c1_mother, TH2D *muon_ptvseta_c2_mother, TH2D *muon_ptvseta_c3_mother, TH2D *gamma_ptvseta_mother, TH2D *gamma_ptvseta_cuts_mother, TH2D *gamma_ptvseta_c1_mother, TH2D *gamma_ptvseta_c2_mother, TH2D *gamma_ptvseta_c3_mother, TH1D *muon_eta_mother, TH1D *gamma_pt_mother, TH1D *gamma_pt_bkg, TH1D *gamma_ptcut_mother, TH1D *gamma_ptetacut_mother, TH1D *gamma_ptetamasscut_mother, TH1D *gamma_p_mother, TH1D *gamma_p_bkg, TH1D *gamma_pcut_mother, TH1D *gamma_petacut_mother, TH1D *gamma_petamasscut_mother, TH1D *gamma_eta_mother, TH1D *mu_number_event, TH1D *gamma_number_event, TH1D *total_number_event, TH1D *mugamma_combnumber_event){

	double inv_mass, mu_invmass, pt1, pt2, pt3, p1, p2, p3, eta1, eta2, eta3;

	Long64_t nentries = v->index.size();

	Long64_t mu_number = 0;
	Long64_t antimu_number = 0;
	Long64_t gamma_number = 0;
	Long64_t mugamma_combnumber = 0;

	Long64_t Ncut0 = 0.0;
	Long64_t Ncut1 = 0.0;
	Long64_t Ncut2 = 0.0;
	Long64_t Ncut3 = 0.0;
	Long64_t Ncut4 = 0.0;

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


							if((pt1>0.5) && (pt2>0.5) && (p1>10.0) && (p2>10.0) && (eta1>2.0) && (eta1<4.5) && (eta2>2.0) && (eta2<4.5)){

								eta_invmass_c3_bkg->Fill(inv_mass);
								muon_invmass_c3_bkg->Fill(mu_invmass);

								if((inv_mass < m_eta+dm_eta) && (inv_mass> m_eta-dm_eta)){

									eta_invmass_cuts_bkg->Fill(inv_mass);
									muon_invmass_cuts_bkg->Fill(mu_invmass);

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

									++Ncut0;

									muon_invmass_mother->Fill(mu_invmass);
									eta_invmass_mothers->Fill(inv_mass);

									//Plot pt vs p 2D histogram before cuts
									muon_ptvsp_mother->Fill(p1, pt1);
									muon_ptvsp_mother->Fill(p2, pt2);
									gamma_ptvsp_mother->Fill(p3, pt3);
									muon_ptvseta_mother->Fill(pt1, eta1);
									muon_ptvseta_mother->Fill(pt2, eta2);
									gamma_ptvseta_mother->Fill(pt3, eta3);

									if((pt1>0.5) && (pt2>0.5)){

										++Ncut1;

										muon_ptvsp_c1_mother->Fill(p1, pt1);
										muon_ptvsp_c1_mother->Fill(p2, pt2);
										gamma_ptvsp_c1_mother->Fill(p3, pt3);
										muon_ptvseta_c1_mother->Fill(pt1, eta1);
										muon_ptvseta_c1_mother->Fill(pt2, eta2);
										gamma_ptvseta_c1_mother->Fill(pt3, eta3);

										if((p1>10.0) && (p2>10.0)){

											++Ncut2;

											muon_ptcut_mother->Fill(pt1);
											muon_ptcut_mother->Fill(pt2);
											gamma_ptcut_mother->Fill(pt3);
											muon_pcut_mother->Fill(p1);
											muon_pcut_mother->Fill(p2);
											gamma_pcut_mother->Fill(p3);

											muon_ptvsp_c2_mother->Fill(p1, pt1);
											muon_ptvsp_c2_mother->Fill(p2, pt2);
											gamma_ptvsp_c2_mother->Fill(p3, pt3);
											muon_ptvseta_c2_mother->Fill(pt1, eta1);
											muon_ptvseta_c2_mother->Fill(pt2, eta2);
											gamma_ptvseta_c2_mother->Fill(pt3, eta3);

											if((eta1>2.0) && (eta1<4.5) && (eta2>2.0) && (eta2<4.5))
											{

												++Ncut3;

												muon_ptetacut_mother->Fill(pt1);
												muon_ptetacut_mother->Fill(pt2);
												gamma_ptetacut_mother->Fill(pt3);
												muon_petacut_mother->Fill(p1);
												muon_petacut_mother->Fill(p2);
												gamma_petacut_mother->Fill(p3);

												muon_ptvsp_c3_mother->Fill(p1, pt1);
												muon_ptvsp_c3_mother->Fill(p2, pt2);
												gamma_ptvsp_c3_mother->Fill(p3, pt3);
												muon_ptvseta_c3_mother->Fill(pt1, eta1);
												muon_ptvseta_c3_mother->Fill(pt2, eta2);
												gamma_ptvseta_c3_mother->Fill(pt3, eta3);

												muon_invmass_c3_mother->Fill(mu_invmass);
												eta_invmass_c3_mother->Fill(inv_mass);

												if((inv_mass < m_eta+dm_eta) && (inv_mass> m_eta-dm_eta)){

													++Ncut4;

													muon_ptetamasscut_mother->Fill(pt1);
													muon_ptetamasscut_mother->Fill(pt2);
													gamma_ptetamasscut_mother->Fill(pt3);
													muon_petamasscut_mother->Fill(p1);
													muon_petamasscut_mother->Fill(p2);
													gamma_petamasscut_mother->Fill(p3);

													//Plot pt vs p 2D histogram after all cuts
													muon_ptvsp_cuts_mother->Fill(p1, pt1);
													muon_ptvsp_cuts_mother->Fill(p2, pt2);
													gamma_ptvsp_cuts_mother->Fill(p3, pt3);
													muon_ptvseta_cuts_mother->Fill(pt1, eta1);
													muon_ptvseta_cuts_mother->Fill(pt2, eta2);
													gamma_ptvseta_cuts_mother->Fill(pt3, eta3);

													muon_invmass_cuts_mother->Fill(mu_invmass);
													eta_invmass_cuts_mother->Fill(inv_mass);
													
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
							else if(v->motherid1[k] != 221){

								gamma_pt_bkg->Fill(pt3);
								gamma_p_bkg->Fill(p3);

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

	if(Ncut0!=0.0){fprintf(file_out1, "%10f %10f %10f %10f \n", (double)Ncut1/Ncut0, (double)Ncut2/Ncut0, (double)Ncut3/Ncut0, (double)Ncut4/Ncut0);}
	else{fprintf(file_out1, "%10f %10f %10f %10f \n", (double)Ncut1, (double)Ncut2, (double)Ncut3,(double)Ncut4);}


}

void analyze_event_pi(FILE *file_out2, vect *v, TH1D *pigamma_invmass, TH1D *pigamma_invmass_c3, TH1D *pigamma_invmass_cuts, TH1D *pi_invmass, TH1D *pi_invmass_c3, TH1D *pi_invmass_cuts, TH1D *pi_pt, TH1D *pi_p, TH1D *pi_eta, TH1D *pi_ptcut, TH1D *pi_pcut, TH1D *pi_ptmasscut, TH1D *pi_pmasscut, TH2D *pi_ptvsp, TH2D *pi_ptvsp_cuts, TH1D *pi_number_event, TH1D *pi_combnumber_event, TH1D *total_pigamma_number_event, TH1D *pigamma_combnumber_event){

	double inv_mass, pion_invmass, pt_pp, pt_pn, p_pp, p_pn, eta_pp, eta_pn;

	Long64_t nentries = v->index.size();

	Long64_t pi_number_p = 0;
	Long64_t pi_number_n = 0;
	Long64_t gamma_number = 0;
	Long64_t pigamma_combnumber = 0;
	Long64_t pi_combnumber = 0;

	Long64_t Ncut0 = 0;
	Long64_t Ncut1 = 0;
	Long64_t Ncut2 = 0;
	Long64_t Ncut3 = 0;
	Long64_t Ncut4 = 0;

	double m_eta = 0.54785;
	double dm_eta = 0.02;

	bool pi_n = false;
	bool gamma = false;

	for(Long64_t i=0;i<nentries;i++){ 

		//if entry is pi+
		if(v->id[i] == 211){

			++pi_number_p;

			//Loop through the event and find an anti-muon
			for(Long64_t j=0;j<nentries;j++){

				//if entry is anti-muon
				if(v->id[j] == -211){

					++pi_combnumber;

					if(pi_n == false){
						++pi_number_n;
					}

					for(Long64_t k=0;k<nentries;k++){

						//if entry is photon
						if(v->id[k] == 22){

							++pigamma_combnumber;

							if(gamma == false){
								++gamma_number;
							}

							//Reconstructing invariant mass of all pi+/- combinations
							inv_mass = invmass(v, i, j, k);
							pion_invmass = invmass2(v, i, j);

							//Plot histogram of all pi+/- combinations
							pigamma_invmass->Fill(inv_mass);
							pi_invmass->Fill(pion_invmass);

							pt_pp = pt(v, i);
							pi_pt->Fill(pt_pp);

							pt_pn = pt(v, j);
							pi_pt->Fill(pt_pn);

							p_pp = p(v, i);
							pi_p->Fill(p_pp);

							p_pn = p(v, j);
							pi_p->Fill(p_pn);

							eta_pp = eta(v, i);
							pi_eta->Fill(eta_pp);

							eta_pn = eta(v, j);
							pi_eta->Fill(eta_pn);

							//Plot pt vs p 2D histogram before cuts
							pi_ptvsp->Fill(p_pp, pt_pp);
							pi_ptvsp->Fill(p_pn, pt_pn);

							++Ncut0;

							if((pt_pp>0.5) && (pt_pn>0.5)){

								++Ncut1;

								if((p_pp>10.0) && (p_pn>10.0)){

									++Ncut2;

									if((eta_pp>2.0) && (eta_pp<4.5) && (eta_pn>2.0) && (eta_pn<4.5)){

										++Ncut3;

										pi_ptcut->Fill(pt_pp);
										pi_ptcut->Fill(pt_pn);
										pi_pcut->Fill(p_pp);
										pi_pcut->Fill(p_pn);

										pigamma_invmass_c3->Fill(inv_mass);
										pi_invmass_c3->Fill(pion_invmass);

										if((inv_mass < m_eta+dm_eta) && (inv_mass > m_eta-dm_eta)){

											++Ncut3;

											pi_ptmasscut->Fill(pt_pp);
											pi_ptmasscut->Fill(pt_pn);
											pi_pmasscut->Fill(p_pp);
											pi_pmasscut->Fill(p_pn);

											//Plot pt vs p 2D histogram after cuts
											pi_ptvsp_cuts->Fill(p_pp, pt_pp);
											pi_ptvsp_cuts->Fill(p_pn, pt_pn);

											pigamma_invmass_cuts->Fill(inv_mass);
											pi_invmass_cuts->Fill(pion_invmass);

										}

									}

								}

							}

						}

					}
					gamma=true;

				}

			}
			pi_n = true;
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
	

	if(Ncut0!=0.0){fprintf(file_out2, "%10f %10f %10f %10f \n", (double)Ncut1/Ncut0, (double)Ncut2/Ncut0, (double)Ncut3/Ncut0, (double)Ncut4/Ncut0);}
	else{fprintf(file_out2, "%10f %10f %10f %10f \n", (double)Ncut1, (double)Ncut2, (double)Ncut3,(double)Ncut4);}


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

	TH1D *muon_invmass_c3_mother = new TH1D("muon_invmass_c3_mother","Reconstructed muon pair invariant mass p_t, p & pseudorapidity cuts (sig)", Nbins2 - 1, Edges2);
    	muon_invmass_cuts_mother -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_c3_mother = new TH1D("eta_invmass_c3_mother","Reconstructed invariant mass from muon pairs + #gamma p_t, p & pseudorapidity cuts (sig)",  Nbins2 - 1, Edges2);
    	eta_invmass_cuts_mother -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *muon_invmass_cuts_bkg = new TH1D("muon_invmass_cuts_bkg","Reconstructed di-muon invariant mass after all cuts (bkg)", Nbins2 - 1, Edges2);
    	muon_invmass_cuts_bkg -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_cuts_bkg = new TH1D("eta_invmass_cuts_bkg","Reconstructed invariant mass from muon pairs + #gamma after all cuts (bkg)",  Nbins2 - 1, Edges2);
    	eta_invmass_cuts_bkg -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *muon_invmass_c3_bkg = new TH1D("muon_invmass_c3_bkg","Reconstructed muon pair invariant mass p_t, p & pseudorapidity cuts (bkg)", Nbins2 - 1, Edges2);
    	muon_invmass_cuts_bkg -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_c3_bkg = new TH1D("eta_invmass_c3_bkg","Reconstructed invariant mass from muon pairs + #gamma p_t, p & pseudorapidity cuts (bkg)",  Nbins2 - 1, Edges2);
    	eta_invmass_cuts_bkg -> GetXaxis()-> SetTitle("m (GeV)");



	TH1D *muon_pt_mother = new TH1D("muon_pt_mother","#mu^{-} and #mu^{+} p_{t} distribution with same mother (#eta)", 500, 0.0, 5.0);
    	muon_pt_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *muon_pt_bkg = new TH1D("muon_pt_bkg","#mu^{-} and #mu^{+} p_{t} distribution for background", 500, 0.0, 5.0);
    	muon_pt_bkg -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *muon_ptcut_mother = new TH1D("muon_ptcut_mother","#mu^{-} and #mu^{+} p_{t} distribution after p_{t} & p cuts", 500, 0.0, 5.0);
    	muon_ptcut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *muon_ptetacut_mother = new TH1D("muon_ptetacut_mother","#mu^{-} and #mu^{+} p_{t} distribution after p_{t}, p & pseudorapidity cuts", 500, 0.0, 5.0);
    	muon_ptetacut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *muon_ptetamasscut_mother = new TH1D("muon_ptetamasscut_mother","#mu^{-} and #mu^{+} p_{t} distribution after all cuts", 500, 0.0, 5.0);
    	muon_ptetamasscut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *muon_p_mother = new TH1D("muon_p_mother","#mu^{-} and #mu^{+} p distribution with same mother (#eta)", 10000, 0.0, 100.0);
    	muon_p_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *muon_p_bkg = new TH1D("muon_p_bkg","#mu^{-} and #mu^{+} p distribution for background", 10000, 0.0, 100.0);
    	muon_p_bkg -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *muon_pcut_mother = new TH1D("muon_pcut_mother","#mu^{-} and #mu^{+} p distribution after p_{t} & p cuts", 10000, 0.0, 100.0);
    	muon_pcut_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *muon_petacut_mother = new TH1D("muon_petacut_mother","#mu^{-} and #mu^{+} p distribution after p_{t}, p & pseudorapidity cuts", 10000, 0.0, 100.0);
    	muon_petacut_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *muon_petamasscut_mother = new TH1D("muon_petamasscut_mother","#mu^{-} and #mu^{+} p distribution after all cuts", 10000, 0.0, 100.0);
    	muon_petamasscut_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *muon_eta_mother = new TH1D("muon_eta_mother","#mu^{-} and #mu^{+} pseudorapidity distribution with same mother (#eta)", 45, 0.0, 4.5);
    	muon_eta_mother -> GetXaxis()-> SetTitle("#eta");




	TH1D *gamma_pt_mother = new TH1D("gamma_pt_mother","#gamma p_{t} distribution with same mother (#eta)", 500, 0.0, 5.0);
    	gamma_pt_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *gamma_pt_bkg = new TH1D("gamma_pt_bkg","#gamma p_{t} distribution for background", 500, 0.0, 5.0);
    	gamma_pt_bkg -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *gamma_ptcut_mother = new TH1D("gamma_ptcut_mother","#gamma p_{t} distribution after p_{t} & p cuts", 500, 0.0, 5.0);
    	gamma_ptcut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *gamma_ptetacut_mother = new TH1D("gamma_ptetacut_mother","#gamma p_{t} distribution after p_{t}, p & pseudorapidity cuts", 500, 0.0, 5.0);
    	gamma_ptetacut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *gamma_ptetamasscut_mother = new TH1D("gamma_ptetamasscut_mother","#gamma p_{t} distribution after all cuts", 500, 0.0, 5.0);
    	gamma_ptetamasscut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *gamma_p_mother = new TH1D("gamma_p_mother","#gamma p distribution with same mother (#eta)", 10000, 0.0, 100.0);
    	gamma_p_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *gamma_p_bkg = new TH1D("gamma_p_bkg","#gamma p distribution for background", 10000, 0.0, 100.0);
    	gamma_p_bkg -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *gamma_pcut_mother = new TH1D("gamma_pcut_mother","#gamma p distribution after p_{t} & p cuts", 10000, 0.0, 100.0);
    	gamma_pcut_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *gamma_petacut_mother = new TH1D("gamma_petacut_mother","#gamma p distribution after p_{t}, p & pseudorapidity cuts", 10000, 0.0, 100.0);
    	gamma_petacut_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *gamma_petamasscut_mother = new TH1D("gamma_petamasscut_mother","#gamma p distribution after all cuts", 10000, 0.0, 100.0);
    	gamma_petamasscut_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *gamma_eta_mother = new TH1D("gamma_eta_mother","#gamma pseudorapidity distribution with same mother (#eta)", 45, 0.0, 4.5);
    	gamma_eta_mother -> GetXaxis()-> SetTitle("#eta");





	TH2D *muon_ptvsp_mother = new TH2D("muon_ptvsp_mother","signal #mu^{#pm} p vs p_t distributions before cuts", 1000, 0.0, 100.0, 50, 0.0, 5.0);
	muon_ptvsp_mother -> GetXaxis()-> SetTitle("p (GeV)");
	muon_ptvsp_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *muon_ptvsp_cuts_mother = new TH2D("muon_ptvsp_cuts_mother","signal #mu^{#pm} p vs p_t distributions after cuts", 1000, 0.0, 100.0, 50, 0.0, 5.0);
	muon_ptvsp_cuts_mother -> GetXaxis()-> SetTitle("p (GeV)");
	muon_ptvsp_cuts_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *muon_ptvsp_c1_mother = new TH2D("muon_ptvsp_c1_mother","signal #mu^{#pm} p vs p_t distributions after p_t cuts", 1000, 0.0, 100.0, 50, 0.0, 5.0);
	muon_ptvsp_c1_mother -> GetXaxis()-> SetTitle("p (GeV)");
	muon_ptvsp_c1_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *muon_ptvsp_c2_mother = new TH2D("muon_ptvsp_c2_mother","signal #mu^{#pm} p vs p_t distributions after p_t & p cuts", 1000, 0.0, 100.0, 50, 0.0, 5.0);
	muon_ptvsp_c2_mother -> GetXaxis()-> SetTitle("p (GeV)");
	muon_ptvsp_c2_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *muon_ptvsp_c3_mother = new TH2D("muon_ptvsp_c3_mother","signal #mu^{#pm} p vs p_t distributions after p_t, p & pseudorapidity cuts", 1000, 0.0, 100.0, 50, 0.0, 5.0);
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



	TH2D *gamma_ptvsp_mother = new TH2D("gamma_ptvsp_mother","signal #gamma p vs p_t distributions before cuts", 1000, 0.0, 100.0, 50, 0.0, 5.0);
	gamma_ptvsp_mother -> GetXaxis()-> SetTitle("p (GeV)");
	gamma_ptvsp_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *gamma_ptvsp_cuts_mother = new TH2D("gamma_ptvsp_cuts_mother","signal #gamma p vs p_t distributions after cuts", 1000, 0.0, 100.0, 50, 0.0, 5.0);
	gamma_ptvsp_cuts_mother -> GetXaxis()-> SetTitle("p (GeV)");
	gamma_ptvsp_cuts_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *gamma_ptvsp_c1_mother = new TH2D("gamma_ptvsp_c1_mother","signal #gamma p vs p_t distributions after p_t cuts", 1000, 0.0, 100.0, 50, 0.0, 5.0);
	gamma_ptvsp_c1_mother -> GetXaxis()-> SetTitle("p (GeV)");
	gamma_ptvsp_c1_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *gamma_ptvsp_c2_mother = new TH2D("gamma_ptvsp_c2_mother","signal #gamma p vs p_t distributions after p_t & p cuts", 1000, 0.0, 100.0, 50, 0.0, 5.0);
	gamma_ptvsp_c2_mother -> GetXaxis()-> SetTitle("p (GeV)");
	gamma_ptvsp_c2_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *gamma_ptvsp_c3_mother = new TH2D("gamma_ptvsp_c3_mother","signal #gamma p vs p_t distributions after p_t, p & pseudorapidity cuts", 1000, 0.0, 100.0, 50, 0.0, 5.0);
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




	TH1D *pigamma_invmass = new TH1D("pigamma_invmass","Reconstructed #eta invariant mass from misID #pi^{+}+#pi^{-}+#gamma", Nbins1 - 1, Edges1);
    	pigamma_invmass -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *pigamma_invmass_c3 = new TH1D("pigamma_invmass_c3","Reconstructed #eta invariant mass from misID #pi^{+}+#pi^{-}+#gamma after p_t, p & #eta cuts", Nbins2 - 1, Edges2);
    	pigamma_invmass_c3 -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *pigamma_invmass_cuts = new TH1D("pigamma_invmass_cuts","Reconstructed #eta invariant mass from misID #pi^{+}+#pi^{-}+#gamma after all cuts", Nbins2 - 1, Edges2);
    	pigamma_invmass_cuts -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *pi_invmass = new TH1D("pi_invmass","Reconstructed #eta invariant mass from misID #pi^{#pm}", Nbins1 - 1, Edges1);
    	pi_invmass -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *pi_invmass_c3 = new TH1D("pi_invmass_c3","Reconstructed #eta invariant mass from misID #pi^{#pm} after p_t, p & pseudorapidity cuts", Nbins2 - 1, Edges2);
    	pi_invmass_c3 -> GetXaxis()-> SetTitle("m (GeV)");
	
	TH1D *pi_invmass_cuts = new TH1D("pi_invmass_cuts","Reconstructed #eta invariant mass from misID #pi^{#pm} after all cuts", Nbins2 - 1, Edges2);
    	pi_invmass_cuts -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *pi_pt = new TH1D("pi_pt","#pi^{-} and #pi^{+} p_{t} distribution", 500, 0.0, 5.0);
    	pi_pt -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *pi_ptcut = new TH1D("pi_ptcut","#pi^{-} and #pi^{+} p_{t} distribution after cuts", 500, 0.0, 5.0);
    	pi_ptcut -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *pi_ptmasscut = new TH1D("pi_ptmasscut","#pi^{-} and #pi^{+} p_{t} distribution after cuts including mass", 500, 0.0, 5.0);
    	pi_ptmasscut -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *pi_p = new TH1D("pi_p","#pi^{-} and #pi^{+} p distribution", 1000, 0.0, 100.0);
    	pi_p -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *pi_pcut = new TH1D("pi_pcut","#pi^{-} and #pi^{+} p distribution after cuts", 1000, 0.0, 100.0);
    	pi_pcut -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *pi_pmasscut = new TH1D("pi_pmasscut","#pi^{-} and #pi^{+} p distribution after cuts including mass", 1000, 0.0, 100.0);
    	pi_pmasscut -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *pi_eta = new TH1D("pi_eta","#pi^{-} and #pi^{+} pseudorapidity distribution", 45, 0.0, 4.5);
    	pi_eta -> GetXaxis()-> SetTitle("#eta");




	TH2D *pi_ptvsp = new TH2D("pi_ptvsp","#pi^{#pm} p vs p_t distributions before cuts", 1000, 0.0, 100.0, 50, 0.0, 5.0);
	pi_ptvsp -> GetXaxis()-> SetTitle("p (GeV)");
	pi_ptvsp -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *pi_ptvsp_cuts = new TH2D("pi_ptvsp_cuts","#pi^{#pm} p vs p_t distributions after cuts", 1000, 0.0, 100.0, 50, 0.0, 5.0);
	pi_ptvsp_cuts -> GetXaxis()-> SetTitle("p (GeV)");
	pi_ptvsp_cuts -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *pi_ptvsp_c1 = new TH2D("pi_ptvsp_c1","#pi^{#pm} p vs p_t distributions after p_t cuts", 1000, 0.0, 100.0, 50, 0.0, 5.0);
	pi_ptvsp_c1 -> GetXaxis()-> SetTitle("p (GeV)");
	pi_ptvsp_c1 -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *pi_ptvsp_c2 = new TH2D("pi_ptvsp_c2","#pi^{#pm} p vs p_t distributions after p_t & p cuts", 1000, 0.0, 100.0, 50, 0.0, 5.0);
	pi_ptvsp_c2 -> GetXaxis()-> SetTitle("p (GeV)");
	pi_ptvsp_c2 -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *pi_ptvsp_c3 = new TH2D("pi_ptvsp_c3","#pi^{#pm} p vs p_t distributions after p_t, p & pseudorapidity cuts", 1000, 0.0, 100.0, 50, 0.0, 5.0);
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

			analyze_event(file_out1, v, eta_invmass, muon_invmass, muongamma_invmass_mother, muon_invmass_mother, eta_invmass_mother, muon_invmass_cuts_mother, eta_invmass_cuts_mother, muon_invmass_c3_mother, eta_invmass_c3_mother, muon_invmass_cuts_bkg, eta_invmass_cuts_bkg, muon_invmass_c3_bkg, eta_invmass_c3_bkg, muon_pt_mother, muon_pt_bkg, muon_ptcut_mother, muon_ptetacut_mother, muon_ptetamasscut_mother, muon_p_mother, muon_p_bkg, muon_pcut_mother, muon_petacut_mother, muon_petamasscut_mother, muon_ptvsp_mother, muon_ptvsp_cuts_mother, muon_ptvsp_c1_mother, muon_ptvsp_c1_mother, muon_ptvsp_c2_mother, gamma_ptvsp_mother, gamma_ptvsp_cuts_mother, gamma_ptvsp_c1_mother, gamma_ptvsp_c2_mother, gamma_ptvsp_c3_mother, muon_ptvseta_mother, muon_ptvseta_cuts_mother, muon_ptvseta_c1_mother, muon_ptvseta_c2_mother, muon_ptvseta_c3_mother, gamma_ptvseta_mother, gamma_ptvseta_cuts_mother, gamma_ptvseta_c1_mother, gamma_ptvseta_c2_mother, gamma_ptvseta_c3_mother, muon_eta_mother, gamma_pt_mother, gamma_pt_bkg, gamma_ptcut_mother, gamma_ptetacut_mother, gamma_ptetamasscut_mother, gamma_p_mother, gamma_p_bkg, gamma_pcut_mother, gamma_petacut_mother, gamma_petamasscut_mother, gamma_eta_mother, mu_number_event, gamma_number_event, total_number_event, mugamma_combnumber_event);

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

			analyze_event(file_out1, v, eta_invmass, muon_invmass, muongamma_invmass_mother, muon_invmass_mother, eta_invmass_mother, muon_invmass_cuts_mother, eta_invmass_cuts_mother, muon_invmass_c3_mother, eta_invmass_c3_mother, muon_invmass_cuts_bkg, eta_invmass_cuts_bkg, muon_invmass_c3_bkg, eta_invmass_c3_bkg, muon_pt_mother, muon_pt_bkg, muon_ptcut_mother, muon_ptetacut_mother, muon_ptetamasscut_mother, muon_p_mother, muon_p_bkg, muon_pcut_mother, muon_petacut_mother, muon_petamasscut_mother, muon_ptvsp_mother, muon_ptvsp_cuts_mother, muon_ptvsp_c1_mother, muon_ptvsp_c1_mother, muon_ptvsp_c2_mother, gamma_ptvsp_mother, gamma_ptvsp_cuts_mother, gamma_ptvsp_c1_mother, gamma_ptvsp_c2_mother, gamma_ptvsp_c3_mother, muon_ptvseta_mother, muon_ptvseta_cuts_mother, muon_ptvseta_c1_mother, muon_ptvseta_c2_mother, muon_ptvseta_c3_mother, gamma_ptvseta_mother, gamma_ptvseta_cuts_mother, gamma_ptvseta_c1_mother, gamma_ptvseta_c2_mother, gamma_ptvseta_c3_mother, muon_eta_mother, gamma_pt_mother, gamma_pt_bkg, gamma_ptcut_mother, gamma_ptetacut_mother, gamma_ptetamasscut_mother, gamma_p_mother, gamma_p_bkg, gamma_pcut_mother, gamma_petacut_mother, gamma_petamasscut_mother, gamma_eta_mother, mu_number_event, gamma_number_event, total_number_event, mugamma_combnumber_event);

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

	//Loop through the entries of tree T3
	nentries = T3->GetEntries();
std::cout<<"Total number of entries T3: "<<nentries<<std::endl;
   	for (Long64_t i=0;i<nentries;i++){

	      	T3->GetEntry(i);

		if(index_var>=1000){break;}

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

			analyze_event_pi(file_out2, v, pigamma_invmass, pigamma_invmass_c3, pigamma_invmass_cuts, pi_invmass, pi_invmass_c3, pi_invmass_cuts, pi_pt, pi_p, pi_eta, pi_ptcut, pi_pcut, pi_ptmasscut, pi_pmasscut, pi_ptvsp, pi_ptvsp_cuts, pi_number_event, pi_combnumber_event, total_pigamma_number_event, pigamma_combnumber_event);

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

			analyze_event_pi(file_out2, v, pigamma_invmass, pigamma_invmass_c3, pigamma_invmass_cuts, pi_invmass, pi_invmass_c3, pi_invmass_cuts, pi_pt, pi_p, pi_eta, pi_ptcut, pi_pcut, pi_ptmasscut, pi_pmasscut, pi_ptvsp, pi_ptvsp_cuts, pi_number_event, pi_combnumber_event, total_pigamma_number_event, pigamma_combnumber_event);

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
	muon_pt_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_pt_mother1.pdf","pdf");

	nentries = muon_pt_bkg->GetEntries();
	muon_pt_bkg->Scale(1.0 / nentries, "width");
	muon_pt_bkg->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_pt_bkg1.pdf","pdf");

	nentries = muon_ptcut_mother->GetEntries();
	muon_ptcut_mother->Scale(1.0 / nentries, "width");
	muon_ptcut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_ptcut_mother1.pdf","pdf");

	nentries = muon_ptetacut_mother->GetEntries();
	muon_ptetacut_mother->Scale(1.0 / nentries, "width");
	muon_ptetacut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_ptetacut_mother1.pdf","pdf");

	nentries = muon_ptetamasscut_mother->GetEntries();
	muon_ptetamasscut_mother->Scale(1.0 / nentries, "width");
	muon_ptetamasscut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_ptetamasscut_mother1.pdf","pdf");

	nentries = muon_p_mother->GetEntries();
	muon_p_mother->Scale(1.0 / nentries, "width");
	muon_p_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_p_mother1.pdf","pdf");

	nentries = muon_p_bkg->GetEntries();
	muon_p_bkg->Scale(1.0 / nentries, "width");
	muon_p_bkg->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_p_bkg1.pdf","pdf");

	nentries = muon_pcut_mother->GetEntries();
	muon_pcut_mother->Scale(1.0 / nentries, "width");
	muon_pcut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_pcut_mother1.pdf","pdf");

	nentries = muon_petacut_mother->GetEntries();
	muon_petacut_mother->Scale(1.0 / nentries, "width");
	muon_petacut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_petacut_mother1.pdf","pdf");

	nentries = muon_petamasscut_mother->GetEntries();
	muon_petamasscut_mother->Scale(1.0 / nentries, "width");
	muon_petamasscut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_petamasscut_mother1.pdf","pdf");


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


	nentries = muon_eta_mother->GetEntries();
	muon_eta_mother->Scale(1.0 / nentries, "width");
	muon_eta_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_pseudorap_mother1.pdf","pdf");





	nentries = gamma_pt_mother->GetEntries();
	gamma_pt_mother->Scale(1.0 / nentries, "width");
	gamma_pt_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_pt_mother1.pdf","pdf");

	nentries = gamma_pt_bkg->GetEntries();
	gamma_pt_bkg->Scale(1.0 / nentries, "width");
	gamma_pt_bkg->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_pt_bkg1.pdf","pdf");

	nentries = gamma_ptcut_mother->GetEntries();
	gamma_ptcut_mother->Scale(1.0 / nentries, "width");
	gamma_ptcut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_ptcut_mother1.pdf","pdf");

	nentries = gamma_ptetacut_mother->GetEntries();
	gamma_ptetacut_mother->Scale(1.0 / nentries, "width");
	gamma_ptetacut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_ptetacut_mother1.pdf","pdf");

	nentries = gamma_ptetamasscut_mother->GetEntries();
	gamma_ptetamasscut_mother->Scale(1.0 / nentries, "width");
	gamma_ptetamasscut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_ptetamasscut_mother1.pdf","pdf");

	nentries = gamma_p_mother->GetEntries();
	gamma_p_mother->Scale(1.0 / nentries, "width");
	gamma_p_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_p_mother1.pdf","pdf");

	nentries = gamma_p_bkg->GetEntries();
	gamma_p_bkg->Scale(1.0 / nentries, "width");
	gamma_p_bkg->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_p_bkg1.pdf","pdf");

	nentries = gamma_pcut_mother->GetEntries();
	gamma_pcut_mother->Scale(1.0 / nentries, "width");
	gamma_pcut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_pcut_mother1.pdf","pdf");

	nentries = gamma_petacut_mother->GetEntries();
	gamma_petacut_mother->Scale(1.0 / nentries, "width");
	gamma_petacut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_petacut_mother1.pdf","pdf");

	nentries = gamma_petamasscut_mother->GetEntries();
	gamma_petamasscut_mother->Scale(1.0 / nentries, "width");
	gamma_petamasscut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_petamasscut_mother1.pdf","pdf");


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


	nentries = gamma_eta_mother->GetEntries();
	gamma_eta_mother->Scale(1.0 / nentries, "width");
	gamma_eta_mother->Draw();

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
	pi_pt->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_pt1.pdf","pdf");

	nentries = pi_p->GetEntries();
	pi_p->Scale(1.0 / nentries, "width");
	pi_p->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_p1.pdf","pdf");

	nentries = pi_eta->GetEntries();
	pi_eta->Scale(1.0 / nentries, "width");
	pi_eta->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_eta1.pdf","pdf");

	nentries = pi_ptcut->GetEntries();
	pi_ptcut->Scale(1.0 / nentries, "width");
	pi_ptcut->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_ptcut1.pdf","pdf");

	nentries = pi_pcut->GetEntries();
	pi_pcut->Scale(1.0 / nentries, "width");
	pi_pcut->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_pcut1.pdf","pdf");

	nentries = pi_ptmasscut->GetEntries();
	pi_ptmasscut->Scale(1.0 / nentries, "width");
	pi_ptmasscut->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_ptmasscut1.pdf","pdf");

	nentries = pi_pmasscut->GetEntries();
	pi_pmasscut->Scale(1.0 / nentries, "width");
	pi_pmasscut->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_pmasscut1.pdf","pdf");

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


	gPad->SetLogy();

	double Br2 = 0.0003100;

	nentries = eta_invmass->GetEntries();
	//eta_invmass->Scale(1.0 / nentries, "width");
	eta_invmass->Scale(Br2*nentries, "width");
	eta_invmass->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_etainvmass1.pdf","pdf");


	nentries = muon_invmass->GetEntries();
	muon_invmass->Scale(Br2*nentries, "width");
	muon_invmass->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_muoninvmass1.pdf","pdf");

	nentries = muon_invmass_mother->GetEntries();
	//muon_invmass_mother->Scale(1.0 / nentries, "width");
	muon_invmass_mother->Scale(Br2*nentries, "width");
	muon_invmass_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_muoninvmass_mother1.pdf","pdf");

	nentries = eta_invmass_mother->GetEntries();
	//eta_invmass_mother->Scale(1.0 / nentries, "width");
	eta_invmass_mother->Scale(Br2*nentries, "width");
	eta_invmass_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_etainvmass_mother1.pdf","pdf");

	nentries = muon_invmass_cuts_mother->GetEntries();
	muon_invmass_cuts_mother->Scale(Br2*nentries, "width");
	muon_invmass_cuts_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_muoninvmass_cuts_mother1.pdf","pdf");
///
	nentries = eta_invmass_cuts_mother->GetEntries();
	eta_invmass_cuts_mother->Scale(Br2*nentries, "width");
	eta_invmass_cuts_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_etainvmass_cuts_mother1.pdf","pdf");

	nentries = muon_invmass_c3_mother->GetEntries();
	muon_invmass_c3_mother->Scale(Br2*nentries, "width");
	muon_invmass_c3_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_muoninvmass_c3_mother1.pdf","pdf");

	nentries = eta_invmass_c3_mother->GetEntries();
	eta_invmass_c3_mother->Scale(Br2*nentries, "width");
	eta_invmass_c3_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_etainvmass_c3_mother1.pdf","pdf");

	nentries = muon_invmass_cuts_bkg->GetEntries();
	muon_invmass_cuts_bkg->Scale(Br2*nentries, "width");
	muon_invmass_cuts_bkg->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_muoninvmass_cuts_bkg1.pdf","pdf");

	nentries = eta_invmass_cuts_bkg->GetEntries();
	eta_invmass_cuts_bkg->Scale(Br2*nentries, "width");
	eta_invmass_cuts_bkg->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_etainvmass_cuts_bkg1.pdf","pdf");

	nentries = muon_invmass_c3_bkg->GetEntries();
	muon_invmass_c3_bkg->Scale(Br2*nentries, "width");
	muon_invmass_c3_bkg->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_muoninvmass_c3_bkg1.pdf","pdf");

	nentries = eta_invmass_c3_bkg->GetEntries();
	eta_invmass_c3_bkg->Scale(Br2*nentries, "width");
	eta_invmass_c3_bkg->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_etainvmass_c3_bkg1.pdf","pdf");

	double misID = 1e-6;

	nentries = pigamma_invmass->GetEntries();
	pigamma_invmass->Scale(misID*nentries, "width");
	pigamma_invmass->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_pigamma_invmass1.pdf","pdf");

	nentries = pigamma_invmass_c3->GetEntries();
	pigamma_invmass_c3->Scale(misID*nentries, "width");
	pigamma_invmass_c3->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_pigamma_invmass_c31.pdf","pdf");

	nentries = pigamma_invmass_cuts->GetEntries();
	pigamma_invmass_cuts->Scale(misID*nentries, "width");
	pigamma_invmass_cuts->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_pigamma_invmass_cuts1.pdf","pdf");

	nentries = pi_invmass->GetEntries();
	pi_invmass->Scale(misID*nentries, "width");
	pi_invmass->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_invmass1.pdf","pdf");

	nentries = pi_invmass_c3->GetEntries();
	pi_invmass_c3->Scale(misID*nentries, "width");
	pi_invmass_c3->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_invmass_c31.pdf","pdf");

	nentries = pi_invmass_cuts->GetEntries();
	pi_invmass_cuts->Scale(misID*nentries, "width");
	pi_invmass_cuts->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project54_pi_invmass_cuts1.pdf","pdf");


	muon_pt_bkg->SetLineColor(kYellow);
	muon_pt_bkg->SetFillStyle(1001);
	muon_pt_bkg->Draw();
	muon_ptcut_mother->SetLineColor(kGreen);
	muon_ptcut_mother->SetFillStyle(1001);
	muon_ptcut_mother->Draw("same");
	muon_ptetacut_mother->SetLineColor(kRed);
	muon_ptetacut_mother->SetFillStyle(1001);
	muon_ptetacut_mother->Draw("same");
	muon_ptetamasscut_mother->SetLineColor(kMagenta);
	muon_ptetamasscut_mother->SetFillStyle(1001);
	muon_ptetamasscut_mother->Draw("same");
	muon_pt_mother->SetLineColor(kBlue);
	muon_pt_mother->SetFillStyle(1001);
	muon_pt_mother->Draw("same");

	TLegend *legend3 = new TLegend(0.5,0.5,0.9,0.7);	
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


	gamma_pt_bkg->SetLineColor(kYellow);
	gamma_pt_bkg->SetFillStyle(1001);
	gamma_pt_bkg->Draw();
	gamma_ptcut_mother->SetLineColor(kGreen);
	gamma_ptcut_mother->SetFillStyle(1001);
	gamma_ptcut_mother->Draw("same");
	gamma_ptetacut_mother->SetLineColor(kRed);
	gamma_ptetacut_mother->SetFillStyle(1001);
	gamma_ptetacut_mother->Draw("same");
	gamma_ptetamasscut_mother->SetLineColor(kMagenta);
	gamma_ptetamasscut_mother->SetFillStyle(1001);
	gamma_ptetamasscut_mother->Draw("same");
	gamma_pt_mother->SetLineColor(kBlue);
	gamma_pt_mother->SetFillStyle(1001);
	gamma_pt_mother->Draw("same");

	TLegend *legend8 = new TLegend(0.5,0.5,0.9,0.7);	
	TLegendEntry *leg81 = legend8->AddEntry("gamma_pt_bkg","p_{t} distribution for background","f");
  	leg81->SetFillColor(kYellow);
	TLegendEntry *leg82 = legend8->AddEntry("gamma_pt_mother","p_{t} distribution before cuts (sig)","f");
  	leg82->SetFillColor(kBlue);
	TLegendEntry *leg83 = legend8->AddEntry("gamma_ptcut_mother","p_{t} distribution after p_{t} & p cuts (sig)","f");
  	leg83->SetFillColor(kGreen);
	TLegendEntry *leg84 = legend8->AddEntry("gamma_ptetacut_mother","p_{t} distribution after p_{t}, p & #eta cuts (sig)","f");
  	leg84->SetFillColor(kRed);
	TLegendEntry *leg85 = legend8->AddEntry("gamma_ptetamasscut_mother","p_{t} distribution after all cuts (sig)","f");
  	leg85->SetFillColor(kMagenta);
	legend8->Draw("same");

	gamma_pt_bkg->SetTitle("#gamma p_{t} distributions");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_pt_mother12.pdf","pdf");



	pi_pt->SetLineColor(kBlue);
	pi_pt->SetFillStyle(1001);
	pi_pt->Draw();
	pi_ptcut->SetLineColor(kRed);
	pi_ptcut->SetFillStyle(1001);
	pi_ptcut->Draw("same");
	pi_ptmasscut->SetLineColor(kGreen);
	pi_ptmasscut->SetFillStyle(1001);
	pi_ptmasscut->Draw("same");

	TLegend *legend4 = new TLegend(0.5,0.5,0.9,0.7);
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
	muon_p_bkg->Draw();
	muon_pcut_mother->SetLineColor(kGreen);
	muon_pcut_mother->SetFillStyle(1001);
	muon_pcut_mother->Draw("same");
	muon_petacut_mother->SetLineColor(kRed);
	muon_petacut_mother->SetFillStyle(1001);
	muon_petacut_mother->Draw("same");
	muon_petamasscut_mother->SetLineColor(kMagenta);
	muon_petamasscut_mother->SetFillStyle(1001);
	muon_petamasscut_mother->Draw("same");
	muon_p_mother->SetLineColor(kBlue);
	muon_p_mother->SetFillStyle(1001);
	muon_p_mother->Draw("same");

	TLegend *legend6 = new TLegend(0.5,0.5,0.9,0.7);
	TLegendEntry *leg5 = legend6->AddEntry("muon_p_mother","p distribution before cuts","f");
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

	c1->Modified();
	c1->Update();
	c1->Print("project54_mu_p_mother12.pdf","pdf");


	gamma_p_bkg->SetLineColor(kYellow);
	gamma_p_bkg->SetFillStyle(1001);
	gamma_p_bkg->Draw();
	gamma_pcut_mother->SetLineColor(kGreen);
	gamma_pcut_mother->SetFillStyle(1001);
	gamma_pcut_mother->Draw("same");
	gamma_petacut_mother->SetLineColor(kRed);
	gamma_petacut_mother->SetFillStyle(1001);
	gamma_petacut_mother->Draw("same");
	gamma_petamasscut_mother->SetLineColor(kMagenta);
	gamma_petamasscut_mother->SetFillStyle(1001);
	gamma_petamasscut_mother->Draw("same");
	gamma_p_mother->SetLineColor(kBlue);
	gamma_p_mother->SetFillStyle(1001);
	gamma_p_mother->Draw("same");

	TLegend *legend9 = new TLegend(0.5,0.5,0.9,0.7);	
	TLegendEntry *leg91 = legend9->AddEntry("gamma_p_bkg","p distribution for background","f");
  	leg91->SetFillColor(kYellow);
	TLegendEntry *leg92 = legend9->AddEntry("gamma_p_mother","p distribution before cuts (sig)","f");
  	leg92->SetFillColor(kBlue);
	TLegendEntry *leg93 = legend9->AddEntry("gamma_pcut_mother","p distribution after p_{t} & p cuts (sig)","f");
  	leg93->SetFillColor(kGreen);
	TLegendEntry *leg94 = legend9->AddEntry("gamma_petacut_mother","p distribution after p_{t}, p & #eta cuts (sig)","f");
  	leg94->SetFillColor(kRed);
	TLegendEntry *leg95 = legend9->AddEntry("gamma_petamasscut_mother","p distribution after all cuts (sig)","f");
  	leg95->SetFillColor(kMagenta);
	legend9->Draw("same");

	gamma_p_bkg->SetTitle("#gamma p distributions");

	c1->Modified();
	c1->Update();
	c1->Print("project54_gamma_p_mother12.pdf","pdf");



	pi_p->SetLineColor(kBlue);
	pi_p->SetFillStyle(1001);
	pi_p->Draw();
	pi_pcut->SetLineColor(kRed);
	pi_pcut->SetFillStyle(1001);
	pi_pcut->Draw("same");
	pi_pmasscut->SetLineColor(kGreen);
	pi_pmasscut->SetFillStyle(1001);
	pi_pmasscut->Draw("same");

	TLegend *legend5 = new TLegend(0.5,0.5,0.9,0.7);
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
	eta_invmass->Draw();
	eta_invmass_mother->SetLineColor(kGreen);
	eta_invmass_mother->SetFillStyle(1001);
	eta_invmass_mother->Draw("same");

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
	eta_invmass->Draw();
	eta_invmass_mother->SetLineColor(kGreen);
	eta_invmass_mother->SetFillStyle(1001);
	eta_invmass_mother->Draw("same");
	pigamma_invmass->SetLineColor(kOrange);
	pigamma_invmass->SetFillStyle(1001);
	pigamma_invmass->Draw("same");

	TLegend *legend7 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg70 = legend7->AddEntry("eta_invmass","#eta invariant mass (sig+bkg)","f");
  	leg70->SetFillColor(kBlue);	
	TLegendEntry *leg71 = legend7->AddEntry("eta_invmass_mother","#eta invariant mass (sig)","f");
  	leg71->SetFillColor(kGreen);
	TLegendEntry *leg72 = legend7->AddEntry("pigamma_invmass","misID 2#pi+#gamma invariant mass","f");
  	leg71->SetFillColor(kOrange);
	legend7->Draw("same");

	eta_invmass->SetAxisRange(1e6, 1e12,"Y");
	eta_invmass->SetTitle("Reconstructed #eta invariant mass");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etapigamma_invmass12.pdf","pdf");


	eta_invmass_cuts_mother->SetLineColor(kBlue);
	eta_invmass_cuts_mother->SetFillStyle(1001);
	eta_invmass_cuts_mother->Draw();
	pigamma_invmass_cuts->SetLineColor(kRed);
	pigamma_invmass_cuts->SetFillStyle(1001);
	pigamma_invmass_cuts->Draw("same");

	TLegend *legend17 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg170 = legend17->AddEntry("eta_invmass_cuts_mother","#eta invariant mass (sig)","f");
  	leg170->SetFillColor(kBlue);	
	TLegendEntry *leg171 = legend17->AddEntry("pigamma_invmass_cuts","misID 2#pi+#gamma invariant mass","f");
  	leg171->SetFillColor(kRed);
	legend17->Draw("same");

	eta_invmass_cuts_mother->SetAxisRange(1e0, 1e4,"Y");
	eta_invmass_cuts_mother->SetTitle("Reconstructed #eta invariant mass after all cuts");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etapigamma_invmass_cuts12.pdf","pdf");


	eta_invmass_c3_mother->SetLineColor(kBlue);
	eta_invmass_c3_mother->SetFillStyle(1001);
	eta_invmass_c3_mother->Draw();
	pigamma_invmass_c3->SetLineColor(kRed);
	pigamma_invmass_c3->SetFillStyle(1001);
	pigamma_invmass_c3->Draw("same");

	TLegend *legend18 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg180 = legend18->AddEntry("eta_invmass_c3_mother","#eta invariant mass (sig)","f");
  	leg180->SetFillColor(kBlue);	
	TLegendEntry *leg181 = legend18->AddEntry("pigamma_invmass_c3","misID 2#pi+#gamma invariant mass","f");
  	leg181->SetFillColor(kRed);
	legend18->Draw("same");

	eta_invmass_c3_mother->SetAxisRange(1e1, 1e7,"Y");
	eta_invmass_c3_mother->SetTitle("Reconstructed #eta invariant mass after p_{t}, p and #eta cuts");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etapigamma_invmass_c312.pdf","pdf");




	muon_invmass->SetLineColor(kBlue);
	muon_invmass->SetFillStyle(1001);
	muon_invmass->Draw();
	muon_invmass_mother->SetLineColor(kGreen);
	muon_invmass_mother->SetFillStyle(1001);
	muon_invmass_mother->Draw("same");
	pi_invmass->SetLineColor(kOrange);
	pi_invmass->SetFillStyle(1001);
	pi_invmass->Draw("same");

	TLegend *legend10 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg100 = legend10->AddEntry("muon_invmass","di-muon invariant mass (sig+bkg)","f");
  	leg100->SetFillColor(kBlue);	
	TLegendEntry *leg101 = legend10->AddEntry("muon_invmass_mother","di-muon invariant mass (sig)","f");
  	leg101->SetFillColor(kGreen);
	TLegendEntry *leg102 = legend10->AddEntry("pi_invmass","misID #pi invariant mass","f");
  	leg101->SetFillColor(kOrange);
	legend10->Draw("same");

	muon_invmass->SetAxisRange(1e6, 1e12,"Y");
	muon_invmass->SetAxisRange(0, 2,"X");
	muon_invmass->SetTitle("Reconstructed di-muon invariant mass");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mupi_invmass12.pdf","pdf");


	muon_invmass_cuts_mother->SetLineColor(kBlue);
	muon_invmass_cuts_mother->SetFillStyle(1001);
	muon_invmass_cuts_mother->Draw();
	pi_invmass_cuts->SetLineColor(kRed);
	pi_invmass_cuts->SetFillStyle(1001);
	pi_invmass_cuts->Draw("same");

	TLegend *legend11 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg110 = legend11->AddEntry("muon_invmass_cuts_mother","di-muon invariant mass (sig)","f");
  	leg110->SetFillColor(kBlue);	
	TLegendEntry *leg111 = legend11->AddEntry("pi_invmass_cuts","misID #pi invariant mass","f");
  	leg111->SetFillColor(kRed);
	legend11->Draw("same");

	muon_invmass_cuts_mother->SetAxisRange(1e-1, 1e3,"Y");
	muon_invmass_cuts_mother->SetTitle("Reconstructed di-muon invariant mass after all cuts");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mupi_invmass_cuts12.pdf","pdf");


	muon_invmass_c3_mother->SetLineColor(kBlue);
	muon_invmass_c3_mother->SetFillStyle(1001);
	muon_invmass_c3_mother->Draw();
	pi_invmass_c3->SetLineColor(kRed);
	pi_invmass_c3->SetFillStyle(1001);
	pi_invmass_c3->Draw("same");

	TLegend *legend12 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg120 = legend12->AddEntry("muon_invmass_c3_mother","di-muon invariant mass (sig)","f");
  	leg120->SetFillColor(kBlue);	
	TLegendEntry *leg121 = legend12->AddEntry("pi_invmass_c3","misID #pi invariant mass","f");
  	leg121->SetFillColor(kRed);
	legend12->Draw("same");

	muon_invmass_c3_mother->SetAxisRange(1e1, 1e7,"Y");
	muon_invmass_c3_mother->SetTitle("Reconstructed di-muon invariant mass after p_{t}, p and #eta cuts");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mupi_invmass_c312.pdf","pdf");



	eta_invmass_cuts_bkg->SetLineColor(kBlue);
	eta_invmass_cuts_bkg->SetFillStyle(1001);
	eta_invmass_cuts_bkg->Draw();
	pigamma_invmass_cuts->SetLineColor(kRed);
	pigamma_invmass_cuts->SetFillStyle(1001);
	pigamma_invmass_cuts->Draw("same");

	TLegend *legend13 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg130 = legend13->AddEntry("eta_invmass_cuts_bkg","#eta invariant mass (bkg)","f");
  	leg130->SetFillColor(kBlue);	
	TLegendEntry *leg131 = legend13->AddEntry("pigamma_invmass_cuts","misID 2#pi+#gamma invariant mass","f");
  	leg131->SetFillColor(kRed);
	legend13->Draw("same");

	eta_invmass_cuts_bkg->SetAxisRange(1e0, 1e4,"Y");
	eta_invmass_cuts_bkg->SetTitle("Reconstructed #eta invariant mass after all cuts (bkg)");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etapigamma_invmass_cuts_bkg12.pdf","pdf");


	eta_invmass_c3_bkg->SetLineColor(kBlue);
	eta_invmass_c3_bkg->SetFillStyle(1001);
	eta_invmass_c3_bkg->Draw();
	pigamma_invmass_c3->SetLineColor(kRed);
	pigamma_invmass_c3->SetFillStyle(1001);
	pigamma_invmass_c3->Draw("same");

	TLegend *legend14 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg140 = legend14->AddEntry("eta_invmass_c3_bkg","#eta invariant mass (bkg)","f");
  	leg140->SetFillColor(kBlue);	
	TLegendEntry *leg141 = legend14->AddEntry("pigamma_invmass_c3","misID 2#pi+#gamma invariant mass","f");
  	leg141->SetFillColor(kRed);
	legend14->Draw("same");

	//eta_invmass_c3_bkg->SetAxisRange(1e6, 1e12,"Y");
	eta_invmass_c3_bkg->SetTitle("Reconstructed #eta invariant mass after p_{t}, p and #eta cuts (bkg)");

	c1->Modified();
	c1->Update();
	c1->Print("project54_etapigamma_invmass_c3_bkg12.pdf","pdf");


	muon_invmass_cuts_bkg->SetLineColor(kBlue);
	muon_invmass_cuts_bkg->SetFillStyle(1001);
	muon_invmass_cuts_bkg->Draw();
	pi_invmass_cuts->SetLineColor(kRed);
	pi_invmass_cuts->SetFillStyle(1001);
	pi_invmass_cuts->Draw("same");

	TLegend *legend15 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg150 = legend15->AddEntry("muon_invmass_cuts_bkg","di-muon invariant mass (bkg)","f");
  	leg150->SetFillColor(kBlue);	
	TLegendEntry *leg151 = legend15->AddEntry("pi_invmass_cuts","misID #pi invariant mass","f");
  	leg151->SetFillColor(kRed);
	legend15->Draw("same");

	muon_invmass_cuts_bkg->SetAxisRange(1e0, 1e4,"Y");
	muon_invmass_cuts_bkg->SetTitle("Reconstructed di-muon invariant mass after all cuts (bkg)");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mupi_invmass_cuts_bkg12.pdf","pdf");


	muon_invmass_c3_bkg->SetLineColor(kBlue);
	muon_invmass_c3_bkg->SetFillStyle(1001);
	muon_invmass_c3_bkg->Draw();
	pi_invmass_c3->SetLineColor(kRed);
	pi_invmass_c3->SetFillStyle(1001);
	pi_invmass_c3->Draw("same");

	TLegend *legend16 = new TLegend(0.5,0.2,0.9,0.4);
	TLegendEntry *leg160 = legend16->AddEntry("muon_invmass_c3_bkg","di-muon invariant mass (bkg)","f");
  	leg160->SetFillColor(kBlue);	
	TLegendEntry *leg161 = legend16->AddEntry("pi_invmass_c3","misID #pi invariant mass","f");
  	leg161->SetFillColor(kRed);
	legend16->Draw("same");

	//muon_invmass_c3_bkg->SetAxisRange(1e6, 1e12,"Y");
	muon_invmass_c3_bkg->SetTitle("Reconstructed di-muon invariant mass after p_{t}, p and #eta cuts (bkg)");

	c1->Modified();
	c1->Update();
	c1->Print("project54_mupi_invmass_c3_bkg12.pdf","pdf");



	file_out->Write();

	delete file_in;
	delete file_out;
	fclose(file_out1);
	fclose(file_out2);

	// Done.
	return 0;
}


