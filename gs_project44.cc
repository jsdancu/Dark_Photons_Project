//27/11/2017
//Program reading in the three trees generated in gs_project24.cc and performing analysis on the reconstructed invariant mass of eta

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

double invmass(vect *v, Long64_t i, Long64_t j){

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

void analyze_event(FILE *file_out1, vect *v, TH1D *eta_invmass, TH1D *muon_invmass_mother, TH1D *eta_invmass_mothers, TH1D *muon_pt_mother, TH1D *muon_pt_bkg, TH1D *muon_ptcut_mother, TH1D *muon_ptetacut_mother, TH1D *muon_ptetamasscut_mother, TH1D *muon_p_mother, TH1D *muon_p_bkg, TH1D *muon_pcut_mother, TH1D *muon_petacut_mother, TH1D *muon_petamasscut_mother, TH2D *muon_ptvsp_mother, TH2D *muon_ptvsp_cuts_mother, TH2D *muon_ptvsp_c1_mother, TH2D *muon_ptvsp_c2_mother, TH2D *muon_ptvsp_c3_mother, TH2D *muon_ptvseta_mother, TH2D *muon_ptvseta_cuts_mother, TH2D *muon_ptvseta_c1_mother, TH2D *muon_ptvseta_c2_mother, TH2D *muon_ptvseta_c3_mother, TH1D *muon_eta_mother, TH1D *mu_number_event, TH1D *mu_combnumber_event){

	double inv_mass;

	Long64_t nentries = v->index.size();

	Long64_t mu_number = 0;
	Long64_t antimu_number = 0;
	Long64_t mu_combnumber = 0;

	Long64_t Ncut0 = 0.0;
	Long64_t Ncut1 = 0.0;
	Long64_t Ncut2 = 0.0;
	Long64_t Ncut3 = 0.0;
	Long64_t Ncut4 = 0.0;

	double m_eta = 0.54785;
	double dm_eta = 0.02;

	bool antimu = false;

	for(Long64_t i=0;i<nentries;i++){ 

		//if entry is muon
		if(v->id[i] == 13){

			++mu_number;

			//Loop through the event and find an anti-muon
			for(Long64_t j=0;j<nentries;j++){

				//if entry is anti-muon
				if(v->id[j] == -13){

					++mu_combnumber;

					if(antimu == false){
						++antimu_number;
					}

					//Reconstructing invariant mass of all muon-anti-muon combinations
					inv_mass = invmass(v, i, j);

					//Plot histogram of all muon-anti-muon combinations
					eta_invmass->Fill(inv_mass);

					double pt1 = pt(v, i);
					double pt2 = pt(v, j);
					double p1 = p(v, i);
					double p2 = p(v, j);
					double eta1 = eta(v, i);
					double eta2 = eta(v, j);

					if((v->mother1[i] == v->mother1[j]) && (v->mother2[i] == v->mother2[j])){

						//If the muon-antimuon pair come from the same mother particle plot them in a separate histogram
						muon_invmass_mother->Fill(inv_mass);

						if(v->motherid1[i] == 221){

							muon_pt_mother->Fill(pt1);
							muon_pt_mother->Fill(pt2);
							muon_p_mother->Fill(p1);
							muon_p_mother->Fill(p2);
							muon_eta_mother->Fill(eta1);
							muon_eta_mother->Fill(eta2);

							++Ncut0;

							eta_invmass_mothers->Fill(inv_mass);

							//Plot pt vs p 2D histogram before cuts
							muon_ptvsp_mother->Fill(p1, pt1);
							muon_ptvsp_mother->Fill(p2, pt2);
							muon_ptvseta_mother->Fill(pt1, eta1);
							muon_ptvseta_mother->Fill(pt2, eta2);

							if((pt1>0.5) && (pt2>0.5)){

								++Ncut1;

								muon_ptvsp_c1_mother->Fill(p1, pt1);
								muon_ptvsp_c1_mother->Fill(p2, pt2);
								muon_ptvseta_c1_mother->Fill(pt1, eta1);
								muon_ptvseta_c1_mother->Fill(pt2, eta2);

								if((p1>10.0) && (p2>10.0)){

									++Ncut2;

									muon_ptcut_mother->Fill(pt1);
									muon_ptcut_mother->Fill(pt2);
									muon_pcut_mother->Fill(p1);
									muon_pcut_mother->Fill(p2);

									muon_ptvsp_c2_mother->Fill(p1, pt1);
									muon_ptvsp_c2_mother->Fill(p2, pt2);
									muon_ptvseta_c2_mother->Fill(pt1, eta1);
									muon_ptvseta_c2_mother->Fill(pt2, eta2);

									if((eta1>2.0) && (eta1<4.5) && (eta2>2.0) && (eta2<4.5))
									{

										++Ncut3;

										muon_ptetacut_mother->Fill(pt1);
										muon_ptetacut_mother->Fill(pt2);
										muon_petacut_mother->Fill(p1);
										muon_petacut_mother->Fill(p2);

										muon_ptvsp_c3_mother->Fill(p1, pt1);
										muon_ptvsp_c3_mother->Fill(p2, pt2);
										muon_ptvseta_c3_mother->Fill(pt1, eta1);
										muon_ptvseta_c3_mother->Fill(pt2, eta2);

										if((inv_mass < m_eta+dm_eta) && (inv_mass> m_eta-dm_eta)){

											++Ncut4;

											muon_ptetamasscut_mother->Fill(pt1);
											muon_ptetamasscut_mother->Fill(pt2);
											muon_petamasscut_mother->Fill(p1);
											muon_petamasscut_mother->Fill(p2);

											//Plot pt vs p 2D histogram after all cuts
											muon_ptvsp_cuts_mother->Fill(p1, pt1);
											muon_ptvsp_cuts_mother->Fill(p2, pt2);
											muon_ptvseta_cuts_mother->Fill(pt1, eta1);
											muon_ptvseta_cuts_mother->Fill(pt2, eta2);

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

				}

			}
			antimu = true;
		}

	}

	//std::cout<<"muon number = "<<mu_number<<setw(5)<<"anti muon = "<<antimu_number<<std::endl;
	double total_mu_antimu = mu_number+antimu_number;
	if (v->index.size() > 0) {

		mu_number_event->SetBinContent((v->index[1])+1,total_mu_antimu);
		mu_combnumber_event->SetBinContent((v->index[1])+1,mu_combnumber);

	}
	//if (v->index.size() > 0) std::cout<<v->index[0]<<std::endl;
	//std::cout<<v << "    " << v->index.size() <<std::endl;

	if(Ncut0!=0.0){fprintf(file_out1, "%10f %10f %10f %10f \n", (double)Ncut1/Ncut0, (double)Ncut2/Ncut0, (double)Ncut3/Ncut0, (double)Ncut4/Ncut0);}
	else{fprintf(file_out1, "%10f %10f %10f %10f \n", (double)Ncut1, (double)Ncut2, (double)Ncut3,(double)Ncut4);}


}

void analyze_event_pi(FILE *file_out2, vect *v, TH1D *pi_invmass, TH1D *pi_pt, TH1D *pi_p, TH1D *pi_eta, TH1D *pi_ptcut, TH1D *pi_pcut, TH1D *pi_pcut1, TH1D *pi_pcut2, TH1D *pi_ptmasscut, TH1D *pi_pmasscut, TH2D *pi_ptvsp, TH2D *pi_ptvsp_cuts, TH2D *pi_ptvsp_c1, TH2D *pi_ptvsp_c2, TH2D *pi_ptvsp_c3, TH2D *pi_ptvseta, TH2D *pi_ptvseta_cuts, TH2D *pi_ptvseta_c1, TH2D *pi_ptvseta_c2, TH2D *pi_ptvseta_c3, TH1D *pi_number_event, TH1D *pi_combnumber_event){

	double inv_mass;

	Long64_t nentries = v->index.size();

	Long64_t pi_number_p = 0;
	Long64_t pi_number_n = 0;
	Long64_t pi_combnumber = 0;

	Long64_t Ncut0 = 0;
	Long64_t Ncut1 = 0;
	Long64_t Ncut2 = 0;
	Long64_t Ncut3 = 0;
	Long64_t Ncut4 = 0;

	double m_eta = 0.54785;
	double dm_eta = 0.02;

	bool pi_n = false;

	for(Long64_t i=0;i<nentries;i++){ 

		//if entry is pi+
		if(v->id[i] == 211){

			++pi_number_p;

			//Loop through the event and find a pi-
			for(Long64_t j=0;j<nentries;j++){

				//if entry is pi-
				if(v->id[j] == -211){

					++pi_combnumber;

					if(pi_n == false){
						++pi_number_n;
					}

					//Reconstructing invariant mass of all pi+/- combinations
					inv_mass = invmass(v, i, j);

					//Plot histogram of all pi+/- combinations
					pi_invmass->Fill(inv_mass);

					double pt_pp = pt(v, i);
					pi_pt->Fill(pt_pp);

					double pt_pn = pt(v, j);
					pi_pt->Fill(pt_pn);

					double p_pp = p(v, i);
					pi_p->Fill(p_pp);

					double p_pn = p(v, j);
					pi_p->Fill(p_pn);

					double eta_pp = eta(v, i);
					pi_eta->Fill(eta_pp);

					double eta_pn = eta(v, j);
					pi_eta->Fill(eta_pn);

					//Plot pt vs p 2D histogram before cuts
					pi_ptvsp->Fill(p_pp, pt_pp);
					pi_ptvsp->Fill(p_pn, pt_pn);
					pi_ptvseta->Fill(pt_pp, eta_pp);
					pi_ptvseta->Fill(pt_pn, eta_pn);

					++Ncut0;

					if((pt_pp>0.5) && (pt_pn>0.5)){

						++Ncut1;

						pi_pcut1->Fill(p_pp);
						pi_pcut1->Fill(p_pn);

						pi_ptvsp_c1->Fill(p_pp, pt_pp);
						pi_ptvsp_c1->Fill(p_pn, pt_pn);
						pi_ptvseta_c1->Fill(pt_pp, eta_pp);
						pi_ptvseta_c1->Fill(pt_pn, eta_pn);

						if((p_pp>10.0) && (p_pn>10.0)){

							++Ncut2;

							pi_pcut2->Fill(p_pp);
							pi_pcut2->Fill(p_pn);

							pi_ptvsp_c2->Fill(p_pp, pt_pp);
							pi_ptvsp_c2->Fill(p_pn, pt_pn);
							pi_ptvseta_c2->Fill(pt_pp, eta_pp);
							pi_ptvseta_c2->Fill(pt_pn, eta_pn);

							if((eta_pp>2.0) && (eta_pp<4.5) && (eta_pn>2.0) && (eta_pn<4.5)){

								++Ncut3;

								pi_ptcut->Fill(pt_pp);
								pi_ptcut->Fill(pt_pn);
								pi_pcut->Fill(p_pp);
								pi_pcut->Fill(p_pn);

								pi_ptvsp_c3->Fill(p_pp, pt_pp);
								pi_ptvsp_c3->Fill(p_pn, pt_pn);
								pi_ptvseta_c3->Fill(pt_pp, eta_pp);
								pi_ptvseta_c3->Fill(pt_pn, eta_pn);

								if((inv_mass < m_eta+dm_eta) && (inv_mass > m_eta-dm_eta)){

									++Ncut4;

									pi_ptmasscut->Fill(pt_pp);
									pi_ptmasscut->Fill(pt_pn);
									pi_pmasscut->Fill(p_pp);
									pi_pmasscut->Fill(p_pn);

									//Plot pt vs p 2D histogram after cuts
									pi_ptvsp_cuts->Fill(p_pp, pt_pp);
									pi_ptvsp_cuts->Fill(p_pn, pt_pn);
									pi_ptvseta_cuts->Fill(pt_pp, eta_pp);
									pi_ptvseta_cuts->Fill(pt_pn, eta_pn);

								}

							}

						}

					}


				}

			}
			pi_n = true;
		}

	}

	//std::cout<<"muon number = "<<mu_number<<setw(5)<<"anti muon = "<<antimu_number<<std::endl;
	double total_pi_pn = pi_number_p+pi_number_n;
	if (v->index.size() > 0) {

		pi_number_event->SetBinContent((v->index[1])+1,total_pi_pn);
		pi_combnumber_event->SetBinContent((v->index[1])+1, pi_combnumber);

	}
	//if (v->index.size() > 0) std::cout<<v->index[0]<<std::endl;
	//std::cout<<v << "    " << v->index.size() <<std::endl;

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
	TFile  *file_in = new TFile("gs_project24.root", "READ");

	// Open the output TFile 
	TFile  *file_out = new TFile("gs_project44.root", "recreate");

	FILE  *file_out1 = fopen("mu_out.txt", "w");
	fprintf(file_out1, "Table showing signal significance of signal muons passing various cuts: \n");
	fprintf(file_out1, "pt     p    eta    mass \n");

	FILE  *file_out2 = fopen("pi_out.txt", "w");
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

	TH1D *eta_invmass = new TH1D("eta_invmass","Reconstructed #eta invariant mass from #mu pairs", Nbins1 - 1, Edges1);
	//TH1D *eta_invmass = new TH1D("eta_invmass","Reconstructed eta invariant mass from muon pairs", 3000, 0.0, 30.0);
    	eta_invmass -> GetXaxis()-> SetTitle("m (GeV)");

	//TH1D *muon_invmass_mother = new TH1D("muon_invmass_mother","Reconstructed muon pair invariant mass with same mother", 3000, 0.0, 30.0);
	TH1D *muon_invmass_mother = new TH1D("muon_invmass_mother","Reconstructed #mu pair invariant mass with same mother", Nbins2 - 1, Edges2);
    	muon_invmass_mother -> GetXaxis()-> SetTitle("m (GeV)");

	//TH1D *eta_invmass_mother = new TH1D("eta_invmass_mother","Reconstructed eta invariant mass from muon pairs with same mother", 100, 0.0, 1.0);
	TH1D *eta_invmass_mother = new TH1D("eta_invmass_mother","Reconstructed #eta invariant mass from #mu pairs with same mother (#eta)",  Nbins2 - 1, Edges2);
    	eta_invmass_mother -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *muon_pt_mother = new TH1D("muon_pt_mother","#mu^{-} and #mu^{+} p_{t} distribution with same mother (#eta)", 500, 0.0, 5.0);
    	muon_pt_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *muon_pt_bkg = new TH1D("muon_pt_bkg","#mu^{-} and #mu^{+} p_{t} distribution for background", 500, 0.0, 5.0);
    	muon_pt_bkg -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *muon_ptcut_mother = new TH1D("muon_ptcut_mother","#mu^{-} and #mu^{+} p_{t} distribution after cuts (with same mother (#eta))", 500, 0.0, 5.0);
    	muon_ptcut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *muon_ptetacut_mother = new TH1D("muon_ptetacut_mother","#mu^{-} and #mu^{+} p_{t} distribution after p_t, p, pseudorapidity cut(with same mother (#eta))", 500, 0.0, 5.0);
    	muon_ptetacut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *muon_ptetamasscut_mother = new TH1D("muon_ptetamasscut_mother","#mu^{-} and #mu^{+} p_{t} distribution after all cuts (with same mother (#eta))", 500, 0.0, 5.0);
    	muon_ptetamasscut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *muon_p_mother = new TH1D("muon_p_mother","#mu^{-} and #mu^{+} p distribution with same mother (#eta)", 10000, 0.0, 100.0);
    	muon_p_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *muon_p_bkg = new TH1D("muon_p_bkg","#mu^{-} and #mu^{+} p distribution for background", 10000, 0.0, 100.0);
    	muon_p_bkg -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *muon_pcut_mother = new TH1D("muon_pcut_mother","#mu^{-} and #mu^{+} p distribution after cuts (with same mother (#eta))", 10000, 0.0, 100.0);
    	muon_pcut_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *muon_petacut_mother = new TH1D("muon_petacut_mother","#mu^{-} and #mu^{+} p distribution after p_t, p, pseudorapidity cut(with same mother (#eta))", 10000, 0.0, 100.0);
    	muon_petacut_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *muon_petamasscut_mother = new TH1D("muon_petamasscut_mother","#mu^{-} and #mu^{+} p distribution after all cuts (with same mother (#eta))", 10000, 0.0, 100.0);
    	muon_petamasscut_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *muon_eta_mother = new TH1D("muon_eta_mother","#mu^{-} and #mu^{+} pseudorapidity distribution with same mother (#eta)", 45, 0.0, 4.5);
    	muon_eta_mother -> GetXaxis()-> SetTitle("#eta");




	TH2D *muon_ptvsp_mother = new TH2D("muon_ptvsp_mother","signal #mu^{#pm} p vs p_t distributions before cuts", 10000, 0.0, 100.0, 50, 0.0, 5.0);
	muon_ptvsp_mother -> GetXaxis()-> SetTitle("p (GeV)");
	muon_ptvsp_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *muon_ptvsp_cuts_mother = new TH2D("muon_ptvsp_cuts_mother","signal #mu^{#pm} p vs p_t distributions after all cuts", 10000, 0.0, 100.0, 50, 0.0, 5.0);
	muon_ptvsp_cuts_mother -> GetXaxis()-> SetTitle("p (GeV)");
	muon_ptvsp_cuts_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *muon_ptvsp_c1_mother = new TH2D("muon_ptvsp_c1_mother","signal #mu^{#pm} p vs p_t distributions after p_t cuts", 10000, 0.0, 100.0, 50, 0.0, 5.0);
	muon_ptvsp_c1_mother -> GetXaxis()-> SetTitle("p (GeV)");
	muon_ptvsp_c1_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *muon_ptvsp_c2_mother = new TH2D("muon_ptvsp_c2_mother","signal #mu^{#pm} p vs p_t distributions after p_t & p cuts", 10000, 0.0, 100.0, 50, 0.0, 5.0);
	muon_ptvsp_c2_mother -> GetXaxis()-> SetTitle("p (GeV)");
	muon_ptvsp_c2_mother -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *muon_ptvsp_c3_mother = new TH2D("muon_ptvsp_c3_mother","signal #mu^{#pm} p vs p_t distributions after p_t, p & pseudorapidity cuts", 10000, 0.0, 100.0, 50, 0.0, 5.0);
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



	TH1D *mu_number_event = new TH1D("mu_number_event","Muon-antimuon number per event", 100000, 0.0, 100000.0);
    	mu_number_event -> GetXaxis()-> SetTitle("event index");
	mu_number_event -> GetYaxis()-> SetTitle("number of #mu^{#pm}");

	TH1D *mu_combnumber_event = new TH1D("mu_combnumber_event","Muon-antimuon number of combinations per event", 100000, 0.0, 100000.0);
    	mu_combnumber_event -> GetXaxis()-> SetTitle("event index");
	mu_combnumber_event -> GetYaxis()-> SetTitle("number of combinations");

	//TH1D *pi_invmass = new TH1D("pi_invmass","Reconstructed #eta invariant mass from misID #pi^{#pm}", 300, 0.0, 30.0);
	TH1D *pi_invmass = new TH1D("pi_invmass","Reconstructed #eta invariant mass from misID #pi^{#pm}", Nbins1 - 1, Edges1);
    	pi_invmass -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *pi_pt = new TH1D("pi_pt","#pi^{-} and #pi^{+} p_{t} distribution", 500, 0.0, 5.0);
    	pi_pt -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *pi_ptcut = new TH1D("pi_ptcut","#pi^{-} and #pi^{+} p_{t} distribution after cuts", 500, 0.0, 5.0);
    	pi_ptcut -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *pi_ptmasscut = new TH1D("pi_ptmasscut","#pi^{-} and #pi^{+} p_{t} distribution after cuts including mass", 500, 0.0, 5.0);
    	pi_ptmasscut -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *pi_p = new TH1D("pi_p","#pi^{-} and #pi^{+} p distribution", 1000, 10.0, 20.0);
    	pi_p -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *pi_pcut = new TH1D("pi_pcut","#pi^{-} and #pi^{+} p distribution after p_t, p & pseudorapidity cuts", 1000, 10.0, 20.0);
    	pi_pcut -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *pi_pcut1 = new TH1D("pi_pcut1","#pi^{-} and #pi^{+} p distribution after p_t cut", 1000, 10.0, 20.0);
    	pi_pcut1 -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *pi_pcut2 = new TH1D("pi_pcut2","#pi^{-} and #pi^{+} p distribution after p_t & p cuts", 1000, 10.0, 20.0);
    	pi_pcut2 -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *pi_pmasscut = new TH1D("pi_pmasscut","#pi^{-} and #pi^{+} p distribution after all cuts", 1000, 10.0, 20.0);
    	pi_pmasscut -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *pi_eta = new TH1D("pi_eta","#pi^{-} and #pi^{+} pseudorapidity distribution", 45, 0.0, 4.5);
    	pi_eta -> GetXaxis()-> SetTitle("#eta");



	TH2D *pi_ptvsp = new TH2D("pi_ptvsp","#pi^{#pm} p vs p_t distributions before cuts", 1000, 10.0, 20.0, 50, 0.0, 5.0);
	pi_ptvsp -> GetXaxis()-> SetTitle("p (GeV)");
	pi_ptvsp -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *pi_ptvsp_cuts = new TH2D("pi_ptvsp_cuts","#pi^{#pm} p vs p_t distributions after all cuts", 1000, 10.0, 20.0, 50, 0.0, 5.0);
	pi_ptvsp_cuts -> GetXaxis()-> SetTitle("p (GeV)");
	pi_ptvsp_cuts -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *pi_ptvsp_c1 = new TH2D("pi_ptvsp_c1","#pi^{#pm} p vs p_t distributions after p_t cuts", 1000, 10.0, 20.0, 50, 0.0, 5.0);
	pi_ptvsp_c1 -> GetXaxis()-> SetTitle("p (GeV)");
	pi_ptvsp_c1 -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *pi_ptvsp_c2 = new TH2D("pi_ptvsp_c2","#pi^{#pm} p vs p_t distributions after p_t & p cuts", 1000, 10.0, 20.0, 50, 0.0, 5.0);
	pi_ptvsp_c2 -> GetXaxis()-> SetTitle("p (GeV)");
	pi_ptvsp_c2 -> GetYaxis()-> SetTitle("p_t (GeV)");

	TH2D *pi_ptvsp_c3 = new TH2D("pi_ptvsp_c3","#pi^{#pm} p vs p_t distributions after p_t, p & pseudorapidity cuts", 1000, 10.0, 20.0, 50, 0.0, 5.0);
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




	TH1D *pi_number_event = new TH1D("pi_number_event","#pi^{#pm} number per event", 100000, 0.0, 100000.0);
    	pi_number_event -> GetXaxis()-> SetTitle("event index");
	pi_number_event -> GetYaxis()-> SetTitle("number of #pi^{#pm}");

	TH1D *pi_combnumber_event = new TH1D("pi_combnumber_event","#pi^{#pm} number of combinations  per event", 100000, 0.0, 100000.0);
    	pi_combnumber_event -> GetXaxis()-> SetTitle("event index");
	pi_combnumber_event -> GetYaxis()-> SetTitle("number of combinations");

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

			analyze_event(file_out1, v, eta_invmass, muon_invmass_mother, eta_invmass_mother, muon_pt_mother, muon_pt_bkg, muon_ptcut_mother, muon_ptetacut_mother, muon_ptetamasscut_mother, muon_p_mother, muon_p_bkg, muon_pcut_mother, muon_petacut_mother, muon_petamasscut_mother, muon_ptvsp_mother, muon_ptvsp_cuts_mother, muon_ptvsp_c1_mother, muon_ptvsp_c2_mother, muon_ptvsp_c3_mother, muon_ptvseta_mother, muon_ptvseta_cuts_mother, muon_ptvseta_c1_mother, muon_ptvseta_c2_mother, muon_ptvseta_c3_mother, muon_eta_mother, mu_number_event, mu_combnumber_event);

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

			analyze_event(file_out1, v, eta_invmass, muon_invmass_mother, eta_invmass_mother, muon_pt_mother, muon_pt_bkg, muon_ptcut_mother, muon_ptetacut_mother, muon_ptetamasscut_mother, muon_p_mother, muon_p_bkg, muon_pcut_mother, muon_petacut_mother, muon_petamasscut_mother, muon_ptvsp_mother, muon_ptvsp_cuts_mother, muon_ptvsp_c1_mother, muon_ptvsp_c2_mother, muon_ptvsp_c3_mother, muon_ptvseta_mother, muon_ptvseta_cuts_mother, muon_ptvseta_c1_mother, muon_ptvseta_c2_mother, muon_ptvseta_c3_mother, muon_eta_mother, mu_number_event, mu_combnumber_event);

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

			analyze_event_pi(file_out2, v, pi_invmass, pi_pt, pi_p, pi_eta, pi_ptcut, pi_pcut, pi_pcut1, pi_pcut2, pi_ptmasscut, pi_pmasscut, pi_ptvsp, pi_ptvsp_cuts, pi_ptvsp_c1, pi_ptvsp_c2, pi_ptvsp_c3, pi_ptvseta, pi_ptvseta_cuts, pi_ptvseta_c1, pi_ptvseta_c2, pi_ptvseta_c3, pi_number_event, pi_combnumber_event);

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

			analyze_event_pi(file_out2, v, pi_invmass, pi_pt, pi_p, pi_eta, pi_ptcut, pi_pcut, pi_pcut1, pi_pcut2, pi_ptmasscut, pi_pmasscut, pi_ptvsp, pi_ptvsp_cuts, pi_ptvsp_c1, pi_ptvsp_c2, pi_ptvsp_c3, pi_ptvseta, pi_ptvseta_cuts, pi_ptvseta_c1, pi_ptvseta_c2, pi_ptvseta_c3, pi_number_event, pi_combnumber_event);

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
	c1->Print("project44_mu_pt_mother1.pdf","pdf");

	nentries = muon_pt_bkg->GetEntries();
	muon_pt_bkg->Scale(1.0 / nentries, "width");
	muon_pt_bkg->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_pt_bkg1.pdf","pdf");

	nentries = muon_ptcut_mother->GetEntries();
	muon_ptcut_mother->Scale(1.0 / nentries, "width");
	muon_ptcut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_ptcut_mother1.pdf","pdf");

	nentries = muon_ptetacut_mother->GetEntries();
	muon_ptetacut_mother->Scale(1.0 / nentries, "width");
	muon_ptetacut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_ptetacut_mother1.pdf","pdf");

	nentries = muon_ptetamasscut_mother->GetEntries();
	muon_ptetamasscut_mother->Scale(1.0 / nentries, "width");
	muon_ptetamasscut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_ptetamasscut_mother1.pdf","pdf");

	nentries = muon_p_mother->GetEntries();
	muon_p_mother->Scale(1.0 / nentries, "width");
	muon_p_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_p_mother1.pdf","pdf");

	nentries = muon_p_bkg->GetEntries();
	muon_p_bkg->Scale(1.0 / nentries, "width");
	muon_p_bkg->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_p_bkg1.pdf","pdf");

	nentries = muon_pcut_mother->GetEntries();
	muon_pcut_mother->Scale(1.0 / nentries, "width");
	muon_pcut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_pcut_mother1.pdf","pdf");

	nentries = muon_petacut_mother->GetEntries();
	muon_petacut_mother->Scale(1.0 / nentries, "width");
	muon_petacut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_petacut_mother1.pdf","pdf");

	nentries = muon_petamasscut_mother->GetEntries();
	muon_petamasscut_mother->Scale(1.0 / nentries, "width");
	muon_petamasscut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_petamasscut_mother1.pdf","pdf");


	//Scaling???
	nentries = muon_ptvsp_mother->GetEntries();
	muon_ptvsp_mother->Scale(1.0 / nentries, "width");
	muon_ptvsp_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_ptvsp_mother1.pdf","pdf");

	//Scaling???
	nentries = muon_ptvsp_cuts_mother->GetEntries();
	muon_ptvsp_cuts_mother->Scale(1.0 / nentries, "width");
	muon_ptvsp_cuts_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_ptvsp_cuts_mother1.pdf","pdf");

	//Scaling???
	nentries = muon_ptvsp_c1_mother->GetEntries();
	muon_ptvsp_c1_mother->Scale(1.0 / nentries, "width");
	muon_ptvsp_c1_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_ptvsp_c1_mother1.pdf","pdf");

	//Scaling???
	nentries = muon_ptvsp_c2_mother->GetEntries();
	muon_ptvsp_c2_mother->Scale(1.0 / nentries, "width");
	muon_ptvsp_c2_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_ptvsp_c2_mother1.pdf","pdf");

	//Scaling???
	nentries = muon_ptvsp_c3_mother->GetEntries();
	muon_ptvsp_c3_mother->Scale(1.0 / nentries, "width");
	muon_ptvsp_c3_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_ptvsp_c3_mother1.pdf","pdf");

	//Scaling???
	nentries = muon_ptvseta_mother->GetEntries();
	muon_ptvseta_mother->Scale(1.0 / nentries, "width");
	muon_ptvseta_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_ptvseta_mother1.pdf","pdf");

	//Scaling???
	nentries = muon_ptvseta_cuts_mother->GetEntries();
	muon_ptvseta_cuts_mother->Scale(1.0 / nentries, "width");
	muon_ptvseta_cuts_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_ptvseta_cuts_mother1.pdf","pdf");

	//Scaling???
	nentries = muon_ptvseta_c1_mother->GetEntries();
	muon_ptvseta_c1_mother->Scale(1.0 / nentries, "width");
	muon_ptvseta_c1_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_ptvseta_c1_mother1.pdf","pdf");

	//Scaling???
	nentries = muon_ptvseta_c2_mother->GetEntries();
	muon_ptvseta_c2_mother->Scale(1.0 / nentries, "width");
	muon_ptvseta_c2_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_ptvseta_c2_mother1.pdf","pdf");

	//Scaling???
	nentries = muon_ptvseta_c3_mother->GetEntries();
	muon_ptvseta_c3_mother->Scale(1.0 / nentries, "width");
	muon_ptvseta_c3_mother->Draw("COLZ");

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_ptvseta_c3_mother1.pdf","pdf");



	nentries = muon_eta_mother->GetEntries();
	muon_eta_mother->Scale(1.0 / nentries, "width");
	muon_eta_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_pseudorap_mother1.pdf","pdf");

	mu_number_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_number_event.pdf","pdf");

	mu_combnumber_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_combnumber_event.pdf","pdf");

	nentries = pi_pt->GetEntries();
	pi_pt->Scale(1.0 / nentries, "width");
	pi_pt->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_pt1.pdf","pdf");

	nentries = pi_p->GetEntries();
	pi_p->Scale(1.0 / nentries, "width");
	pi_p->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_p1.pdf","pdf");

	nentries = pi_eta->GetEntries();
	pi_eta->Scale(1.0 / nentries, "width");
	pi_eta->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_eta1.pdf","pdf");

	nentries = pi_ptcut->GetEntries();
	pi_ptcut->Scale(1.0 / nentries, "width");
	pi_ptcut->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_ptcut1.pdf","pdf");

	nentries = pi_pcut->GetEntries();
	pi_pcut->Scale(1.0 / nentries, "width");
	pi_pcut->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_pcut1.pdf","pdf");

	nentries = pi_pcut1->GetEntries();
	pi_pcut1->Scale(1.0 / nentries, "width");
	pi_pcut1->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_pcut11.pdf","pdf");

	nentries = pi_pcut2->GetEntries();
	pi_pcut2->Scale(1.0 / nentries, "width");
	pi_pcut2->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_pcut21.pdf","pdf");

	nentries = pi_ptmasscut->GetEntries();
	pi_ptmasscut->Scale(1.0 / nentries, "width");
	pi_ptmasscut->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_ptmasscut1.pdf","pdf");

	nentries = pi_pmasscut->GetEntries();
	pi_pmasscut->Scale(1.0 / nentries, "width");
	pi_pmasscut->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_pmasscut1.pdf","pdf");

	//Scaling???
	nentries = pi_ptvsp->GetEntries();
	pi_ptvsp->Scale(1.0 / nentries, "width");
	pi_ptvsp->Draw("COLZ");
	
	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_ptvsp1.pdf","pdf");

	//Scaling???
	nentries = pi_ptvsp_cuts->GetEntries();
	pi_ptvsp_cuts->Scale(1.0 / nentries, "width");
	pi_ptvsp_cuts->Draw("COLZ");
	
	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_ptvsp_cuts1.pdf","pdf");

	//Scaling???
	nentries = pi_ptvsp_c1->GetEntries();
	pi_ptvsp_c1->Scale(1.0 / nentries, "width");
	pi_ptvsp_c1->Draw("COLZ");
	
	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_ptvsp_c11.pdf","pdf");

	//Scaling???
	nentries = pi_ptvsp_c2->GetEntries();
	pi_ptvsp_c2->Scale(1.0 / nentries, "width");
	pi_ptvsp_c2->Draw("COLZ");
	
	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_ptvsp_c21.pdf","pdf");

	//Scaling???
	nentries = pi_ptvsp_c3->GetEntries();
	pi_ptvsp_c3->Scale(1.0 / nentries, "width");
	pi_ptvsp_c3->Draw("COLZ");
	
	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_ptvsp_c31.pdf","pdf");

	//Scaling???
	nentries = pi_ptvseta->GetEntries();
	pi_ptvseta->Scale(1.0 / nentries, "width");
	pi_ptvseta->Draw("COLZ");
	
	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_ptvseta1.pdf","pdf");

	//Scaling???
	nentries = pi_ptvseta_cuts->GetEntries();
	pi_ptvseta_cuts->Scale(1.0 / nentries, "width");
	pi_ptvseta_cuts->Draw("COLZ");
	
	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_ptvseta_cuts1.pdf","pdf");

	//Scaling???
	nentries = pi_ptvseta_c1->GetEntries();
	pi_ptvseta_c1->Scale(1.0 / nentries, "width");
	pi_ptvseta_c1->Draw("COLZ");
	
	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_ptvseta_c11.pdf","pdf");

	//Scaling???
	nentries = pi_ptvseta_c2->GetEntries();
	pi_ptvseta_c2->Scale(1.0 / nentries, "width");
	pi_ptvseta_c2->Draw("COLZ");
	
	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_ptvseta_c21.pdf","pdf");

	//Scaling???
	nentries = pi_ptvseta_c3->GetEntries();
	pi_ptvseta_c3->Scale(1.0 / nentries, "width");
	pi_ptvseta_c3->Draw("COLZ");
	
	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_ptvseta_c31.pdf","pdf");

	pi_number_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_number_event1.pdf","pdf");

	pi_combnumber_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_combnumber_event1.pdf","pdf");


	gPad->SetLogy();

	double Br1 = 0.0000060;

	//nentries = eta_invmass->GetEntries();
	//eta_invmass->Scale(1.0 / nentries, "width");
	eta_invmass->Scale(Br1, "width");
	eta_invmass->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_etainvmass1.pdf","pdf");

	//nentries = muon_invmass_mother->GetEntries();
	//muon_invmass_mother->Scale(1.0 / nentries, "width");
	muon_invmass_mother->Scale(Br1, "width");
	muon_invmass_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_muoninvmass_mother1.pdf","pdf");

	//nentries = eta_invmass_mother->GetEntries();
	//eta_invmass_mother->Scale(1.0 / nentries, "width");
	eta_invmass_mother->Scale(Br1, "width");
	eta_invmass_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_etainvmass_mother1.pdf","pdf");

	double misID = 1e-6;

	//nentries = pi_invmass->GetEntries();
	//pi_invmass->Scale(1.0 / nentries, "width");
	pi_invmass->Scale(misID, "width");
	pi_invmass->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_invmass1.pdf","pdf");

	//muon_pt_bkg->SetFillColor(kYellow);
	muon_pt_bkg->SetLineColor(kYellow);
	muon_pt_bkg->SetFillStyle(1001);
	muon_pt_bkg->Draw();
	//muon_ptcut_mother->SetFillColor(kGreen);
	muon_ptcut_mother->SetLineColor(kGreen);
	muon_ptcut_mother->SetFillStyle(1001);
	muon_ptcut_mother->Draw("same");
	//muon_ptetacut_mother->SetFillColor(kRed);
	muon_ptetacut_mother->SetLineColor(kRed);
	muon_ptetacut_mother->SetFillStyle(1001);
	muon_ptetacut_mother->Draw("same");
	muon_ptetamasscut_mother->SetLineColor(kMagenta);
	muon_ptetamasscut_mother->SetFillStyle(1001);
	muon_ptetamasscut_mother->Draw("same");
	//muon_pt_mother->SetFillColor(kBlue);
	muon_pt_mother->SetLineColor(kBlue);
	muon_pt_mother->SetFillStyle(1001);
	muon_pt_mother->Draw("same");

	TLegend *legend3 = new TLegend(0.5,0.5,0.9,0.7);	
	TLegendEntry *leg3 = legend3->AddEntry("muon_pt_mother","p_{t} distribution before cuts","f");
  	leg3->SetFillColor(kBlue);
	TLegendEntry *leg32 = legend3->AddEntry("muon_pt_bkg","p_{t} distribution for background","f");
  	leg32->SetFillColor(kYellow);
	TLegendEntry *leg33 = legend3->AddEntry("muon_ptcut_mother","p_{t} distribution after cuts","f");
  	leg33->SetFillColor(kGreen);
	TLegendEntry *leg31 = legend3->AddEntry("muon_ptetacut_mother","p_{t} distribution after cuts including pseudorapidity","f");
  	leg31->SetFillColor(kRed);
	TLegendEntry *leg34 = legend3->AddEntry("muon_ptetamasscut_mother","p_{t} distribution after cuts including pseudorapidity and mass","f");
  	leg34->SetFillColor(kMagenta);
	legend3->Draw("same");

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_pt_mother12.pdf","pdf");

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
	TLegendEntry *leg42 = legend4->AddEntry("pi_ptcut","#pi^{#pm} p_{t} distribution after cuts","f");
  	leg42->SetFillColor(kRed);
	TLegendEntry *leg43 = legend4->AddEntry("pi_ptmasscut","#pi^{#pm} p_{t} distribution after cuts including mass","f");
  	leg43->SetFillColor(kGreen);
	legend4->Draw("same");

	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_pt12.pdf","pdf");

	//muon_p_bkg->SetFillColor(kYellow);
	muon_p_bkg->SetLineColor(kYellow);
	muon_p_bkg->SetFillStyle(1001);
	muon_p_bkg->Draw();
	//muon_pcut_mother->SetFillColor(kGreen);
	muon_pcut_mother->SetLineColor(kGreen);
	muon_pcut_mother->SetFillStyle(1001);
	muon_pcut_mother->Draw("same");
	//muon_petacut_mother->SetFillColor(kRed);
	muon_petacut_mother->SetLineColor(kRed);
	muon_petacut_mother->SetFillStyle(1001);
	muon_petacut_mother->Draw("same");
	muon_petamasscut_mother->SetLineColor(kMagenta);
	muon_petamasscut_mother->SetFillStyle(1001);
	muon_petamasscut_mother->Draw("same");
	//muon_p_mother->SetFillColor(kBlue);
	muon_p_mother->SetLineColor(kBlue);
	muon_p_mother->SetFillStyle(1001);
	muon_p_mother->Draw("same");

	TLegend *legend6 = new TLegend(0.5,0.5,0.9,0.7);
	TLegendEntry *leg5 = legend6->AddEntry("muon_p_mother","p distribution before cuts","f");
  	leg3->SetFillColor(kBlue);	
	TLegendEntry *leg62 = legend6->AddEntry("muon_p_bkg","p distribution for background","f");
  	leg62->SetFillColor(kYellow);
	TLegendEntry *leg6 = legend6->AddEntry("muon_pcut_mother","p distribution after cuts","f");
  	leg6->SetFillColor(kGreen);
	TLegendEntry *leg61 = legend6->AddEntry("muon_petacut_mother","p distribution after cuts including pseudorapidity","f");
  	leg61->SetFillColor(kRed);
	TLegendEntry *leg63 = legend6->AddEntry("muon_petamasscut_mother","p distribution after cuts including pseudorapidity and mass","f");
  	leg63->SetFillColor(kMagenta);
	legend6->Draw("same");

	c1->Modified();
	c1->Update();
	c1->Print("project44_mu_p_mother12.pdf","pdf");

	pi_p->SetLineColor(kBlue);
	pi_p->SetFillStyle(1001);
	pi_p->Draw();
	pi_pcut1->SetLineColor(kMagenta);
	pi_pcut1->SetFillStyle(1001);
	pi_pcut1->Draw("same");
	pi_pcut2->SetLineColor(kOrange);
	pi_pcut2->SetFillStyle(1001);
	pi_pcut2->Draw("same");
	pi_pcut->SetLineColor(kRed);
	pi_pcut->SetFillStyle(1001);
	pi_pcut->Draw("same");
	pi_pmasscut->SetLineColor(kGreen);
	pi_pmasscut->SetFillStyle(1001);
	pi_pmasscut->Draw("same");

	TLegend *legend5 = new TLegend(0.1,0.7,0.5,0.9);
	TLegendEntry *leg51 = legend5->AddEntry("pi_p","#pi^{#pm} p distribution before cuts","f");
  	leg51->SetFillColor(kBlue);
	TLegendEntry *leg54 = legend5->AddEntry("pi_pcut1","#pi^{#pm} p distribution after p_{t} cut","f");
  	leg54->SetFillColor(kMagenta);
	TLegendEntry *leg55 = legend5->AddEntry("pi_pcut2","#pi^{#pm} p distribution after p_{t} & p cuts","f");
  	leg55->SetFillColor(kOrange);	
	TLegendEntry *leg52 = legend5->AddEntry("pi_pcut","#pi^{#pm} p distribution after p_{t}, p & #eta cuts","f");
  	leg52->SetFillColor(kRed);
	TLegendEntry *leg53 = legend5->AddEntry("pi_pmasscut","#pi^{#pm} p distribution after all cuts","f");
  	leg53->SetFillColor(kGreen);
	legend5->Draw("same");

	pi_p->SetAxisRange(1e-3, 1e0,"Y");

	c1->Modified();
	c1->Update();
	c1->Print("project44_pi_p12.pdf","pdf");

	eta_invmass_mother->SetLineColor(kBlue);
	eta_invmass_mother->SetFillColor(kWhite);
	eta_invmass_mother->SetFillStyle(1001);
	eta_invmass_mother->Draw();
	pi_invmass->SetLineColor(kOrange);
	pi_invmass->SetFillStyle(1001);
	pi_invmass->Draw("same");

	TLegend *legend7 = new TLegend(0.7,0.5,0.9,0.7);	
	TLegendEntry *leg71 = legend7->AddEntry("eta_invmass_mother","#eta invariant mass (sig)","f");
  	leg71->SetFillColor(kBlue);
	TLegendEntry *leg72 = legend7->AddEntry("pi_invmass","misID #pi invariant mass","f");
  	leg72->SetFillColor(kOrange);
	legend7->Draw("same");

	c1->Modified();
	c1->Update();
	c1->Print("project44_etapi_invmass12.pdf","pdf");

	eta_invmass->SetFillColor(kBlue);
	eta_invmass->SetFillStyle(1001);
	eta_invmass->Draw();
	eta_invmass_mother->SetFillColor(kGreen);
	eta_invmass_mother->SetFillStyle(1001);
	eta_invmass_mother->Draw("same");

	TLegend *legend1 = new TLegend(0.7,0.5,0.9,0.7);
	TLegendEntry *leg1 = legend1->AddEntry("eta_invmass","#eta invariant mass (sig+bkg)","f");
  	leg1->SetFillColor(kBlue);	
	TLegendEntry *leg2 = legend1->AddEntry("eta_invmass_mother","#eta invariant mass (sig)","f");
  	leg2->SetFillColor(kGreen);
	legend1->Draw("same");

	gPad->SetLogx();

	c1->Modified();
	c1->Update();
	c1->Print("project44_etainvmass12.pdf","pdf");


	file_out->Write();

	delete file_in;
	delete file_out;
	fclose(file_out1);
	fclose(file_out2);

	delete v;

	// Done.
	return 0;
}


