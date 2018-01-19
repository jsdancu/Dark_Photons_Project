//06/01/2018
//Program reading in the three trees generated in gs_project24.cc and gs_project32.cc and performing analysis on the reconstructed invariant mass of eta
//Program created based on gs_project44.cc and gs_project52.cc

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

double invmass1(vect *v, Long64_t i, Long64_t j){

	double inv_mass=sqrt(pow(v->energy[i]+v->energy[j],2)-pow(v->px[i]+v->px[j],2)-pow(v->py[i]+v->py[j],2)-pow(v->pz[i]+v->pz[j],2));

	return inv_mass;

}

double invmass2(vect *v, Long64_t i, Long64_t j, Long64_t k){

	double inv_mass=sqrt(pow(v->energy[i]+v->energy[j]+v->energy[k],2)-pow(v->px[i]+v->px[j]+v->px[k],2)-pow(v->py[i]+v->py[j]+v->py[k],2)-pow(v->pz[i]+v->pz[j]+v->pz[k],2));

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

//void analyze_event1(vect *v, TTree *T1, TH1D *eta_invmass_mother1, TH1D *pt_mother1, TH1D *p_mother1, TH1D *eta_invmass_cuts_mother1, TH1D *eta_massdiff_mother1){
void analyze_event1(vect *v,  TH1D *eta_invmass_mother1, TH1D *pt_mother1, TH1D *p_mother1, TH1D *eta_invmass_cuts_mother1, TH1D *eta_invmass_c3_mother1, TH1D *eta_invmass_bkg1, TH1D *eta_invmass_cuts_bkg1, TH1D *eta_invmass_c3_bkg1, TH1D *eta_massdiff_mother1, TH1D *mu_combnumber_event1, TH1D *mu_combnumber_cuts_event1){

	double inv_mass, pt1, pt2, p1, p2, eta1, eta2;

	double m_eta = 0.54785;
	double dm_eta = 0.02;

	int mu_combnumber = 0;
	int mu_combnumber_cuts = 0;
/*
	double index_var1, id_var1, energy_var1, mass_var1, px_var1, py_var1, pz_var1, mother1_var1, motherid2_var1, motherid1_var1, mother2_var1;
	double massdiff;

	//Get branches of the trees
	T1->SetBranchAddress("index",&index_var1);
   	T1->SetBranchAddress("id",&id_var1);
	T1->SetBranchAddress("energy",&energy_var1);
   	T1->SetBranchAddress("mass",&mass_var1);
	T1->SetBranchAddress("px",&px_var1);
   	T1->SetBranchAddress("py",&py_var1);
	T1->SetBranchAddress("pz",&pz_var1);
	T1->SetBranchAddress("mother1",&mother1_var1);
	T1->SetBranchAddress("mother2",&mother2_var1);
	T1->SetBranchAddress("motherid1",&motherid1_var1);
	T1->SetBranchAddress("motherid2",&motherid2_var1);*/

	Long64_t nentries = v->index.size();

	for(Long64_t i=0;i<nentries;i++){ 

		//if entry is muon
		if(v->id[i] == 13){

			//Loop through the event and find an anti-muon
			for(Long64_t j=0;j<nentries;j++){

				//if entry is anti-muon
				if(v->id[j] == -13){

					++mu_combnumber;

					//Reconstructing invariant mass of all muon-anti-muon combinations
					inv_mass = invmass1(v, i, j);

					pt1 = pt(v, i);
					pt2 = pt(v, j);
					p1 = p(v, i);
					p2 = p(v, j);
					eta1 = eta(v, i);
					eta2 = eta(v, j);


					//Building the sig+bkg histograms
					eta_invmass_bkg1->Fill(inv_mass);

					if((pt1>0.5) && (pt2>0.5) && (p1>10.0) && (p2>10.0) && (eta1>2.0) && (eta1<4.5) && (eta2>2.0) && (eta2<4.5)){

						eta_invmass_c3_bkg1->Fill(inv_mass);

						if((inv_mass < m_eta+dm_eta) && (inv_mass> m_eta-dm_eta)){

							eta_invmass_cuts_bkg1->Fill(inv_mass);

						}

					}


					//If muons are signal
					if((v->mother1[i] == v->mother1[j]) && (v->mother2[i] == v->mother2[j]) && (v->motherid1[i] == 221)){

						eta_invmass_mother1->Fill(inv_mass);

						pt_mother1->Fill(pt1);
						pt_mother1->Fill(pt2);
						p_mother1->Fill(p1);
						p_mother1->Fill(p2);

/*
						T1->GetEntry(0);
						int k = 1;
						//searches for the mother etas
						while(index_var1 < v->index[i]){

							T1->GetEntry(k);
							k++;

						}

						if(index_var1==v->index[i]){

							massdiff = abs(mass_var1-inv_mass);
							eta_massdiff_mother1->Fill(massdiff);

						}*/


						if((pt1>0.5) && (pt2>0.5) && (p1>10.0) && (p2>10.0) && (eta1>2.0) && (eta1<4.5) && (eta2>2.0) && (eta2<4.5)){

							eta_invmass_c3_mother1->Fill(inv_mass);

							if((inv_mass < m_eta+dm_eta) && (inv_mass> m_eta-dm_eta)){

								eta_invmass_cuts_mother1->Fill(inv_mass);

								++mu_combnumber_cuts;

							}

						}

					}

				}

			}
			
		}

	}

	if (v->index.size() > 0) {

		mu_combnumber_event1->SetBinContent((v->index[1])+1,mu_combnumber);
		mu_combnumber_cuts_event1->SetBinContent((v->index[1])+1,mu_combnumber_cuts);

	}

}

void analyze_event_pi(vect *v, TH1D *pi_invmass, TH1D *pi_pt, TH1D *pi_p, TH1D *pi_invmass_cuts, TH1D *pi_invmass_c3, TH1D *pi_combnumber_event, TH1D *pi_combnumber_cuts_event){

	double inv_mass, pt1, pt2, p1, p2, eta1, eta2;

	double m_eta = 0.54785;
	double dm_eta = 0.02;

	int pi_combnumber = 0;
	int pi_combnumber_cuts = 0;

	Long64_t nentries = v->index.size();

	for(Long64_t i=0;i<nentries;i++){ 

		//if entry is pi+
		if(v->id[i] == 211){

			//Loop through the event and find a pi-
			for(Long64_t j=0;j<nentries;j++){

				//if entry is pi-
				if(v->id[j] == -211){

					++pi_combnumber;

					//Reconstructing invariant mass of all pi+/- combinations
					inv_mass = invmass1(v, i, j);

					//Plot histogram of all pi+/- combinations
					pi_invmass->Fill(inv_mass);

					pt1 = pt(v, i);
					pi_pt->Fill(pt1);

					pt2 = pt(v, j);
					pi_pt->Fill(pt2);

					p1 = p(v, i);
					pi_p->Fill(p1);

					p2 = p(v, j);
					pi_p->Fill(p2);

					eta1 = eta(v, i);

					eta2 = eta(v, j);


					if((pt1>0.5) && (pt2>0.5) && (p1>10.0) && (p2>10.0) && (eta1>2.0) && (eta1<4.5) && (eta2>2.0) && (eta2<4.5)){

						pi_invmass_c3->Fill(inv_mass);

						if((inv_mass < m_eta+dm_eta) && (inv_mass> m_eta-dm_eta)){

							pi_invmass_cuts->Fill(inv_mass);

							++pi_combnumber_cuts;

						}

					}

				}

			}
		
		}

	}

	if (v->index.size() > 0) {

		pi_combnumber_event->SetBinContent((v->index[1])+1,pi_combnumber);
		pi_combnumber_cuts_event->SetBinContent((v->index[1])+1,pi_combnumber_cuts);

	}

}

void analyze_event2(vect *v, TH1D *eta_invmass_mother2, TH1D *pt_mother2, TH1D *p_mother2, TH1D *eta_invmass_cuts_mother2, TH1D *eta_invmass_c3_mother2, TH1D *dimuons_invmass_mother2, TH1D *dimuons_invmass_cuts_mother2, TH1D *dimuons_invmass_c3_mother2, TH1D *eta_invmass_bkg2, TH1D *eta_invmass_cuts_bkg2, TH1D *eta_invmass_c3_bkg2, TH1D *dimuons_invmass_bkg2, TH1D *dimuons_invmass_cuts_bkg2, TH1D *dimuons_invmass_c3_bkg2, TH1D *eta_massdiff_mother2, TH1D *mugamma_combnumber_event2, TH1D *mugamma_combnumber_cuts_event2){

	double inv_mass, inv_mass2, pt1, pt2, p1, p2, eta1, eta2;

	double m_eta = 0.54785;
	double dm_eta = 0.02;

	int mugamma_combnumber = 0;
	int mugamma_combnumber_cuts = 0;

	Long64_t nentries = v->index.size();

	for(Long64_t i=0;i<nentries;i++){

		//if entry is muon
		if(v->id[i] == 13){

			//Loop through the event and find an anti-muon
			for(Long64_t j=0;j<nentries;j++){

				//if entry is anti-muon
				if(v->id[j] == -13){

					for(Long64_t k=0;k<nentries;k++){

						//if entry is photon
						if(v->id[k] == 22){

							++mugamma_combnumber;

							//Reconstructing invariant mass of all muon-antimuon-photon combinations
							inv_mass = invmass2(v, i, j, k);

							inv_mass2 = invmass1(v, i, j);
							pt1 = pt(v, i);
							pt2 = pt(v, j);
							p1 = p(v, i);
							p2 = p(v, j);
							eta1 = eta(v, i);
							eta2 = eta(v, j);


							//Building sig+bkg histograms
							eta_invmass_bkg2->Fill(inv_mass);
							dimuons_invmass_bkg2->Fill(inv_mass2);

							if((pt1>0.5) && (pt2>0.5) && (p1>10.0) && (p2>10.0) && (eta1>2.0) && (eta1<4.5) && (eta2>2.0) && (eta2<4.5)){

								eta_invmass_c3_bkg2->Fill(inv_mass);
								dimuons_invmass_c3_bkg2->Fill(inv_mass2);

								if((inv_mass < m_eta+dm_eta) && (inv_mass> m_eta-dm_eta)){

									eta_invmass_cuts_bkg2->Fill(inv_mass);
									dimuons_invmass_cuts_bkg2->Fill(inv_mass2);

								}

							}



							if((v->mother1[i] == v->mother1[j]) && (v->mother1[i] == v->mother1[k]) && (v->mother2[i] == v->mother2[j]) && (v->mother2[i] == v->mother2[k]) && (v->motherid1[i] == 221)){

								eta_invmass_mother2->Fill(inv_mass);

								dimuons_invmass_mother2->Fill(inv_mass2);
								pt_mother2->Fill(pt1);
								pt_mother2->Fill(pt2);	
								p_mother2->Fill(p1);	
								p_mother2->Fill(p2);

								


								if((pt1>0.5) && (pt2>0.5) && (p1>10.0) && (p2>10.0) && (eta1>2.0) && (eta1<4.5) && (eta2>2.0) && (eta2<4.5)){

									eta_invmass_c3_mother2->Fill(inv_mass);
									dimuons_invmass_c3_mother2->Fill(inv_mass2);

									if((inv_mass < m_eta+dm_eta) && (inv_mass> m_eta-dm_eta)){

										eta_invmass_cuts_mother2->Fill(inv_mass);
										dimuons_invmass_cuts_mother2->Fill(inv_mass2);

										++mugamma_combnumber_cuts;

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

		mugamma_combnumber_event2->SetBinContent((v->index[1])+1,mugamma_combnumber);
		mugamma_combnumber_cuts_event2->SetBinContent((v->index[1])+1,mugamma_combnumber_cuts);

	}

}

int main() {
/*
	//Loading lhcbStyle for plotting
	gROOT->ProcessLine(".L lhcbstyle.C");
	lhcbStyle();

	// Example of adding stat box - turned off by default in plots
	gStyle->SetOptStat("emr");  // show only nent - e , mean - m , rms - r
*/
	// Open the input TFiles 
	TFile  *file_in1 = new TFile("gs_project24.root", "READ");
	TFile  *file_in2 = new TFile("gs_project32.root", "READ");

	// Open the output TFile 
	TFile  *file_out = new TFile("gs_project4452.root", "recreate");

	//Get trees from file
	TTree *T1 = (TTree*)file_in1->Get("T1");
	TTree *T2 = (TTree*)file_in1->Get("T2");
	TTree *T3 = (TTree*)file_in1->Get("T3");

	TTree *T4 = (TTree*)file_in2->Get("T1");
	TTree *T5 = (TTree*)file_in2->Get("T2");

	double index_var, id_var, energy_var, mass_var, px_var, py_var, pz_var, mother1_var, motherid2_var, motherid1_var, mother2_var, index_var1, id_var1, energy_var1, mass_var1, px_var1, py_var1, pz_var1, mother1_var1, motherid2_var1, motherid1_var1, mother2_var1;

	//Get branches of the trees
	T1->SetBranchAddress("index",&index_var1);
   	T1->SetBranchAddress("id",&id_var1);
	T1->SetBranchAddress("energy",&energy_var1);
   	T1->SetBranchAddress("mass",&mass_var1);
	T1->SetBranchAddress("px",&px_var1);
   	T1->SetBranchAddress("py",&py_var1);
	T1->SetBranchAddress("pz",&pz_var1);
	T1->SetBranchAddress("mother1",&mother1_var1);
	T1->SetBranchAddress("mother2",&mother2_var1);
	T1->SetBranchAddress("motherid1",&motherid1_var1);
	T1->SetBranchAddress("motherid2",&motherid2_var1);

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

	T4->SetBranchAddress("index",&index_var1);
   	T4->SetBranchAddress("id",&id_var1);
	T4->SetBranchAddress("energy",&energy_var1);
   	T4->SetBranchAddress("mass",&mass_var1);
	T4->SetBranchAddress("px",&px_var1);
   	T4->SetBranchAddress("py",&py_var1);
	T4->SetBranchAddress("pz",&pz_var1);
	T4->SetBranchAddress("mother1",&mother1_var1);
	T4->SetBranchAddress("mother2",&mother2_var1);
	T4->SetBranchAddress("motherid1",&motherid1_var1);
	T4->SetBranchAddress("motherid2",&motherid2_var1);

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


	TH1D *eta_invmass_mother1 = new TH1D("eta_invmass_mother1","Reconstructed #eta invariant mass from di-muons",  Nbins2 - 1, Edges2);
    	eta_invmass_mother1 -> GetXaxis()-> SetTitle("m (GeV)");
	
	TH1D *pi_invmass = new TH1D("pi_invmass","Reconstructed #eta invariant mass from misID #pi^{#pm}", Nbins1 - 1, Edges1);
    	pi_invmass -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_mother2 = new TH1D("eta_invmass_mother2","Reconstructed invariant mass from #mu pairs + #gamma ",  Nbins2 - 1, Edges2);
    	eta_invmass_mother2 -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *dimuons_invmass_mother2 = new TH1D("dimuons_invmass_mother2","Di-muon invariant mass distributions from #eta -> #mu #mu #gamma ",  Nbins2 - 1, Edges2);
    	dimuons_invmass_mother2 -> GetXaxis()-> SetTitle("m (GeV)");



	TH1D *eta_invmass_cuts_mother1 = new TH1D("eta_invmass_cuts_mother1","Reconstructed #eta invariant mass from di-muons after all cuts",  Nbins2 - 1, Edges2);
    	eta_invmass_cuts_mother1 -> GetXaxis()-> SetTitle("m (GeV)");
	
	TH1D *pi_invmass_cuts = new TH1D("pi_invmass_cuts","Reconstructed #eta invariant mass from misID #pi^{#pm} after all cuts", Nbins1 - 1, Edges1);
    	pi_invmass_cuts -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_cuts_mother2 = new TH1D("eta_invmass_cuts_mother2","Reconstructed invariant mass from #mu pairs + #gamma after all cuts",  Nbins2 - 1, Edges2);
    	eta_invmass_cuts_mother2 -> GetXaxis()-> SetTitle("m (GeV)");	

	TH1D *dimuons_invmass_cuts_mother2 = new TH1D("dimuons_invmass_cuts_mother2","Di-muon invariant mass distributions from #eta -> #mu #mu #gamma after all cuts",  Nbins2 - 1, Edges2);
    	dimuons_invmass_cuts_mother2 -> GetXaxis()-> SetTitle("m (GeV)");



	TH1D *eta_invmass_c3_mother1 = new TH1D("eta_invmass_c3_mother1","Reconstructed #eta invariant mass from di-muons after p_t, p & pseudorapidity cuts",  Nbins2 - 1, Edges2);
    	eta_invmass_c3_mother1 -> GetXaxis()-> SetTitle("m (GeV)");
	
	TH1D *pi_invmass_c3 = new TH1D("pi_invmass_c3","Reconstructed #eta invariant mass from misID #pi^{#pm} after p_t, p & pseudorapidity cuts", Nbins1 - 1, Edges1);
    	pi_invmass_c3 -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_c3_mother2 = new TH1D("eta_invmass_c3_mother2","Reconstructed invariant mass from #mu pairs + #gamma after p_t, p & pseudorapidity cuts",  Nbins2 - 1, Edges2);
    	eta_invmass_c3_mother2 -> GetXaxis()-> SetTitle("m (GeV)");	

	TH1D *dimuons_invmass_c3_mother2 = new TH1D("dimuons_invmass_c3_mother2","Di-muon invariant mass distributions from #eta -> #mu #mu #gamma after p_t, p & pseudorapidity cuts",  Nbins2 - 1, Edges2);
    	dimuons_invmass_c3_mother2 -> GetXaxis()-> SetTitle("m (GeV)");





	TH1D *eta_invmass_bkg1 = new TH1D("eta_invmass_bkg1","Reconstructed #eta invariant mass from di-muons (bkg)",  Nbins2 - 1, Edges2);
    	eta_invmass_bkg1 -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_bkg2 = new TH1D("eta_invmass_bkg2","Reconstructed invariant mass from #mu pairs + #gamma (bkg)",  Nbins2 - 1, Edges2);
    	eta_invmass_bkg2 -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *dimuons_invmass_bkg2 = new TH1D("dimuons_invmass_bkg2","Di-muon invariant mass distributions from #eta -> #mu #mu #gamma (bkg)",  Nbins2 - 1, Edges2);
    	dimuons_invmass_bkg2 -> GetXaxis()-> SetTitle("m (GeV)");



	TH1D *eta_invmass_cuts_bkg1 = new TH1D("eta_invmass_cuts_bkg1","Reconstructed #eta invariant mass from di-muons after all cuts (bkg)",  Nbins2 - 1, Edges2);
    	eta_invmass_cuts_bkg1 -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_cuts_bkg2 = new TH1D("eta_invmass_cuts_bkg2","Reconstructed invariant mass from #mu pairs + #gamma after all cuts (bkg)",  Nbins2 - 1, Edges2);
    	eta_invmass_cuts_bkg2 -> GetXaxis()-> SetTitle("m (GeV)");	

	TH1D *dimuons_invmass_cuts_bkg2 = new TH1D("dimuons_invmass_cuts_bkg2","Di-muon invariant mass distributions from #eta -> #mu #mu #gamma after all cuts (bkg)",  Nbins2 - 1, Edges2);
    	dimuons_invmass_cuts_bkg2 -> GetXaxis()-> SetTitle("m (GeV)");



	TH1D *eta_invmass_c3_bkg1 = new TH1D("eta_invmass_c3_bkg1","Reconstructed #eta invariant mass from di-muons after p_t, p & pseudorapidity cuts (bkg)",  Nbins2 - 1, Edges2);
    	eta_invmass_c3_bkg1 -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_c3_bkg2 = new TH1D("eta_invmass_c3_bkg2","Reconstructed invariant mass from #mu pairs + #gamma after p_t, p & pseudorapidity cuts (bkg)",  Nbins2 - 1, Edges2);
    	eta_invmass_c3_bkg2 -> GetXaxis()-> SetTitle("m (GeV)");	

	TH1D *dimuons_invmass_c3_bkg2 = new TH1D("dimuons_invmass_c3_bkg2","Di-muon invariant mass distributions from #eta -> #mu #mu #gamma after p_t, p & pseudorapidity cuts (bkg)",  Nbins2 - 1, Edges2);
    	dimuons_invmass_c3_bkg2 -> GetXaxis()-> SetTitle("m (GeV)");





	TH1D *eta_massdiff_mother1 = new TH1D("eta_massdiff_mother1","Mass difference between #eta and #mu+#mu",  100000, 0.0, 0.01);
    	eta_massdiff_mother1 -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_massdiff_mother2 = new TH1D("eta_massdiff_mother2","Mass difference between #eta and #mu+#mu+#gamma",  100000, 0.0, 0.01);
    	eta_massdiff_mother2 -> GetXaxis()-> SetTitle("m (GeV)");


	TH1D *pt_mother1 = new TH1D("pt_mother1","#mu^{#pm} p_{t} distribution", 500, 0.0, 5.0);
    	pt_mother1 -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *pt_mother2 = new TH1D("pt_mother2","#mu^{#pm} p_{t} distribution", 500, 0.0, 5.0);
    	pt_mother2 -> GetXaxis()-> SetTitle("p_{t} (GeV)");




	TH1D *pi_pt = new TH1D("pi_pt","#pi^{#pm} p_{t} distribution", 500, 0.0, 5.0);
    	pi_pt -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *p_mother1 = new TH1D("p_mother1","#mu^{#pm} p distribution", 10000, 0.0, 100.0);
    	p_mother1 -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *p_mother2 = new TH1D("p_mother2","#mu^{#pm} p distribution", 10000, 0.0, 100.0);
    	p_mother2 -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *pi_p = new TH1D("pi_p","#pi^{#pm} p distribution", 10000, 0.0, 100.0);
    	pi_p -> GetXaxis()-> SetTitle("p (GeV)");



	TH1D *mu_combnumber_event1 = new TH1D("mu_combnumber_event1","Muon-antimuon number of combinations per event", 100000, 0.0, 100000.0);
    	mu_combnumber_event1 -> GetXaxis()-> SetTitle("event index");
	mu_combnumber_event1 -> GetYaxis()-> SetTitle("number of combinations");

	TH1D *pi_combnumber_event = new TH1D("pi_combnumber_event","#pi^{#pm} number of combinations per event", 100000, 0.0, 100000.0);
    	pi_combnumber_event -> GetXaxis()-> SetTitle("event index");
	pi_combnumber_event -> GetYaxis()-> SetTitle("number of combinations");

	TH1D *mugamma_combnumber_event2 = new TH1D("mugamma_combnumber_event2","#mu^{#pm}+#gamma number of combinations per event", 100000, 0.0, 100000.0);
    	mugamma_combnumber_event2 -> GetXaxis()-> SetTitle("event index");
	mugamma_combnumber_event2 -> GetYaxis()-> SetTitle("number of combinations");

	TH1D *mu_combnumber_cuts_event1 = new TH1D("mu_combnumber_cuts_event1","Muon-antimuon number of combinations per event after cuts", 100000, 0.0, 100000.0);
    	mu_combnumber_cuts_event1 -> GetXaxis()-> SetTitle("event index");
	mu_combnumber_cuts_event1 -> GetYaxis()-> SetTitle("number of combinations");

	TH1D *pi_combnumber_cuts_event = new TH1D("pi_combnumber_cuts_event","#pi^{#pm} number of combinations per event after cuts", 100000, 0.0, 100000.0);
    	pi_combnumber_cuts_event -> GetXaxis()-> SetTitle("event index");
	pi_combnumber_cuts_event -> GetYaxis()-> SetTitle("number of combinations");

	TH1D *mugamma_combnumber_cuts_event2 = new TH1D("mugamma_combnumber_cuts_event2","#mu^{#pm}+#gamma number of combinations per event after cuts", 100000, 0.0, 100000.0);
    	mugamma_combnumber_cuts_event2 -> GetXaxis()-> SetTitle("event index");
	mugamma_combnumber_cuts_event2 -> GetYaxis()-> SetTitle("number of combinations");

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

			//analyze_event1(v, T1, eta_invmass_mother1, pt_mother1, p_mother1, eta_invmass_cuts_mother1, eta_massdiff_mother1);
			analyze_event1(v, eta_invmass_mother1, pt_mother1, p_mother1, eta_invmass_cuts_mother1, eta_invmass_c3_mother1, eta_invmass_bkg1, eta_invmass_cuts_bkg1, eta_invmass_c3_bkg1, eta_massdiff_mother1, mu_combnumber_event1, mu_combnumber_cuts_event1);

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

			//analyze_event1(v, T1, eta_invmass_mother1, pt_mother1, p_mother1, eta_invmass_cuts_mother1, eta_massdiff_mother1);
			analyze_event1(v,  eta_invmass_mother1, pt_mother1, p_mother1, eta_invmass_cuts_mother1, eta_invmass_c3_mother1, eta_invmass_bkg1, eta_invmass_cuts_bkg1, eta_invmass_c3_bkg1, eta_massdiff_mother1, mu_combnumber_event1, mu_combnumber_cuts_event1);

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

			analyze_event_pi(v, pi_invmass, pi_pt, pi_p, pi_invmass_cuts, pi_invmass_c3, pi_combnumber_event, pi_combnumber_cuts_event);

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

			analyze_event_pi(v, pi_invmass, pi_pt, pi_p, pi_invmass_cuts, pi_invmass_c3, pi_combnumber_event, pi_combnumber_cuts_event);

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

			analyze_event2(v, eta_invmass_mother2, pt_mother2, p_mother2, eta_invmass_cuts_mother2, eta_invmass_c3_mother2, dimuons_invmass_mother2, dimuons_invmass_cuts_mother2, dimuons_invmass_c3_mother2, eta_invmass_bkg2, eta_invmass_cuts_bkg2, eta_invmass_c3_bkg2, dimuons_invmass_bkg2, dimuons_invmass_cuts_bkg2, dimuons_invmass_c3_bkg2, eta_massdiff_mother2, mugamma_combnumber_event2, mugamma_combnumber_cuts_event2);

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

			analyze_event2(v, eta_invmass_mother2, pt_mother2, p_mother2, eta_invmass_cuts_mother2, eta_invmass_c3_mother2, dimuons_invmass_mother2, dimuons_invmass_cuts_mother2, dimuons_invmass_c3_mother2, eta_invmass_bkg2, eta_invmass_cuts_bkg2, eta_invmass_c3_bkg2, dimuons_invmass_bkg2, dimuons_invmass_cuts_bkg2, dimuons_invmass_c3_bkg2, eta_massdiff_mother2, mugamma_combnumber_event2, mugamma_combnumber_cuts_event2);

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

	nentries = pt_mother1->GetEntries();
	pt_mother1->Scale(1.0 / nentries, "width");
	pt_mother1->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_mu_pt_mother1.pdf","pdf");

	nentries = pi_pt->GetEntries();
	pi_pt->Scale(1.0 / nentries, "width");
	pi_pt->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_pi_pt.pdf","pdf");

	nentries = pt_mother2->GetEntries();
	pt_mother2->Scale(1.0 / nentries, "width");
	pt_mother2->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_mu_pt_mother2.pdf","pdf");


	nentries = p_mother1->GetEntries();
	p_mother1->Scale(1.0 / nentries, "width");
	p_mother1->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_mu_p_mother1.pdf","pdf");

	nentries = pi_p->GetEntries();
	pi_p->Scale(1.0 / nentries, "width");
	pi_p->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_pi_p.pdf","pdf");

	nentries = p_mother2->GetEntries();
	p_mother2->Scale(1.0 / nentries, "width");
	p_mother2->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_mu_p_mother2.pdf","pdf");


	mu_combnumber_event1->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_mu_combnumber_event1.pdf","pdf");

	mugamma_combnumber_event2->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_mugamma_combnumber_event2.pdf","pdf");

	pi_combnumber_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_pi_combnumber_event1.pdf","pdf");

	mu_combnumber_cuts_event1->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_mu_combnumber_cuts_event1.pdf","pdf");

	mugamma_combnumber_cuts_event2->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_mugamma_combnumber_cuts_event2.pdf","pdf");

	pi_combnumber_cuts_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_pi_combnumber_cuts_event1.pdf","pdf");


	gPad->SetLogy();


	pt_mother1->SetLineColor(kBlue);
	pt_mother1->SetFillStyle(1001);
	pt_mother1->Draw();
	pt_mother2->SetLineColor(kRed);
	pt_mother2->SetFillStyle(1001);
	pt_mother2->Draw("same");
	pi_pt->SetLineColor(kGreen);
	pi_pt->SetFillStyle(1001);
	pi_pt->Draw("same");

	TLegend *legend2 = new TLegend(0.5,0.5,0.9,0.7);	
	TLegendEntry *leg21 = legend2->AddEntry("pt_mother1","di-muon p_{t} from #mu #mu","f");
  	leg21->SetFillColor(kBlue);
	TLegendEntry *leg22 = legend2->AddEntry("pt_mother2","di-muon p_{t} from #mu #mu #gamma","f");
  	leg22->SetFillColor(kRed);
	TLegendEntry *leg23 = legend2->AddEntry("pi_pt","misID #pi p_{t}","f");
  	leg23->SetFillColor(kGreen);
	legend2->Draw("same");

	pt_mother1->SetTitle("Reconstructed p_{t} distributions");

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etapi_pt123.pdf","pdf");


	p_mother1->SetLineColor(kBlue);
	p_mother1->SetFillStyle(1001);
	p_mother1->Draw();
	p_mother2->SetLineColor(kRed);
	p_mother2->SetFillStyle(1001);
	p_mother2->Draw("same");
	pi_p->SetLineColor(kGreen);
	pi_p->SetFillStyle(1001);
	pi_p->Draw("same");

	TLegend *legend3 = new TLegend(0.5,0.5,0.9,0.7);	
	TLegendEntry *leg31 = legend3->AddEntry("p_mother1","di-muon p from #mu #mu","f");
  	leg31->SetFillColor(kBlue);
	TLegendEntry *leg32 = legend3->AddEntry("p_mother2","di-muon p from #mu #mu #gamma","f");
  	leg32->SetFillColor(kRed);
	TLegendEntry *leg33 = legend3->AddEntry("pi_p","misID #pi p","f");
  	leg33->SetFillColor(kGreen);
	legend3->Draw("same");

	p_mother1->SetTitle("Reconstructed p distributions");

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etapi_p123.pdf","pdf");
	

	double Br1 = 0.0000060;

	//nentries = eta_invmass_mother1->GetEntries();
	//eta_invmass_mother1->Scale(1.0 / nentries, "width");
	//eta_invmass_mother1->Scale(Br1);
	eta_invmass_mother1->Scale(Br1, "width");
	eta_invmass_mother1->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etainvmass_mother1.pdf","pdf");

	double Br2 = 0.0003100;

	//nentries = eta_invmass_mother2->GetEntries();
	//eta_invmass_mother2->Scale(1.0 / nentries, "width");
	//eta_invmass_mother2->Scale(Br2);
	eta_invmass_mother2->Scale(Br2, "width");
	eta_invmass_mother2->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etainvmass_mother2.pdf","pdf");


	//nentries = dimuons_invmass_mother2->GetEntries();
	//dimuons_invmass_mother2->Scale(1.0 / nentries, "width");
	//dimuons_invmass_mother2->Scale(Br2);
	dimuons_invmass_mother2->Scale(Br2, "width");
	dimuons_invmass_mother2->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_dimuon_invmass_mother2.pdf","pdf");

	dimuons_invmass_cuts_mother2->Scale(Br2, "width");
	dimuons_invmass_cuts_mother2->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_dimuon_invmass_cuts_mother2.pdf","pdf");


	double misID = 1e-6;

	//nentries = pi_invmass->GetEntries();
	//pi_invmass->Scale(1.0 / nentries, "width");
	//pi_invmass->Scale(misID);
	pi_invmass->Scale(misID, "width");
	pi_invmass->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_pi_invmass1.pdf","pdf");



	//nentries = eta_invmass_cuts_mother1->GetEntries();
	//eta_invmass_cuts_mother1->Scale(1.0 / nentries, "width");
	//eta_invmass_cuts_mother1->Scale(Br1);
	eta_invmass_cuts_mother1->Scale(Br1, "width");
	eta_invmass_cuts_mother1->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etainvmass_cuts_mother1.pdf","pdf");


	//nentries = eta_invmass_cuts_mother2->GetEntries();
	//eta_invmass_cuts_mother2->Scale(1.0 / nentries, "width");
	//eta_invmass_cuts_mother2->Scale(Br2);
	eta_invmass_cuts_mother2->Scale(Br2, "width");
	eta_invmass_cuts_mother2->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etainvmass_cuts_mother2.pdf","pdf");


	//nentries = pi_invmass_cuts->GetEntries();
	//pi_invmass_cuts->Scale(1.0 / nentries, "width");
	//pi_invmass_cuts->Scale(misID);
	pi_invmass_cuts->Scale(misID, "width");
	pi_invmass_cuts->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_pi_invmass_cuts1.pdf","pdf");



	eta_invmass_c3_mother1->Scale(Br1, "width");
	eta_invmass_c3_mother1->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etainvmass_c3_mother1.pdf","pdf");


	eta_invmass_c3_mother2->Scale(Br2, "width");
	eta_invmass_c3_mother2->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etainvmass_c3_mother2.pdf","pdf");


	pi_invmass_c3->Scale(misID, "width");
	pi_invmass_c3->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_pi_invmass_c31.pdf","pdf");

	
	dimuons_invmass_c3_mother2->Scale(Br2, "width");
	dimuons_invmass_c3_mother2->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_dimuon_invmass_c3_mother2.pdf","pdf");




	eta_invmass_bkg1->Scale(Br1, "width");
	eta_invmass_bkg1->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etainvmass_bkg1.pdf","pdf");

	eta_invmass_bkg2->Scale(Br2, "width");
	eta_invmass_bkg2->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etainvmass_bkg2.pdf","pdf");

	dimuons_invmass_bkg2->Scale(Br2, "width");
	dimuons_invmass_bkg2->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_dimuon_invmass_bkg2.pdf","pdf");

	dimuons_invmass_cuts_bkg2->Scale(Br2, "width");
	dimuons_invmass_cuts_bkg2->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_dimuon_invmass_cuts_bkg2.pdf","pdf");

	eta_invmass_cuts_bkg1->Scale(Br1, "width");
	eta_invmass_cuts_bkg1->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etainvmass_cuts_bkg1.pdf","pdf");

	eta_invmass_cuts_bkg2->Scale(Br2, "width");
	eta_invmass_cuts_bkg2->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etainvmass_cuts_bkg2.pdf","pdf");

	eta_invmass_c3_bkg1->Scale(Br1, "width");
	eta_invmass_c3_bkg1->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etainvmass_c3_bkg1.pdf","pdf");

	eta_invmass_c3_bkg2->Scale(Br2, "width");
	eta_invmass_c3_bkg2->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etainvmass_c3_bkg2.pdf","pdf");
	
	dimuons_invmass_c3_bkg2->Scale(Br2, "width");
	dimuons_invmass_c3_bkg2->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_dimuon_invmass_c3_bkg2.pdf","pdf");

	

	eta_invmass_cuts_mother2->SetLineColor(kRed);
	eta_invmass_cuts_mother2->SetFillStyle(1001);
	eta_invmass_cuts_mother2->Draw();
	eta_invmass_cuts_mother1->SetLineColor(kBlue);
	eta_invmass_cuts_mother1->SetFillStyle(1001);
	eta_invmass_cuts_mother1->Draw("same");
	pi_invmass_cuts->SetLineColor(kGreen);
	pi_invmass_cuts->SetFillStyle(1001);
	pi_invmass_cuts->Draw("same");

	TLegend *legend4 = new TLegend(0.5,0.5,0.9,0.7);	
	TLegendEntry *leg42 = legend4->AddEntry("eta_invmass_cuts_mother2","#eta invariant mass from #mu #mu #gamma","f");
  	leg42->SetFillColor(kRed);
	TLegendEntry *leg41 = legend4->AddEntry("eta_invmass_cuts_mother1","#eta invariant mass from #mu #mu","f");
  	leg41->SetFillColor(kBlue);	
	TLegendEntry *leg43 = legend4->AddEntry("pi_invmass_cuts","misID #pi invariant mass","f");
  	leg43->SetFillColor(kGreen);
	legend4->Draw("same");

	eta_invmass_cuts_mother1->SetTitle("Reconstructed #eta invariant mass after cuts");

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etapi_invmass_cuts123.pdf","pdf");


	dimuons_invmass_mother2->SetLineColor(kRed);
	dimuons_invmass_mother2->SetFillStyle(1001);
	dimuons_invmass_mother2->Draw();
	eta_invmass_mother1->SetLineColor(kBlue);
	eta_invmass_mother1->SetFillStyle(1001);
	eta_invmass_mother1->Draw("same");
	pi_invmass->SetLineColor(kGreen);
	pi_invmass->SetFillStyle(1001);
	pi_invmass->Draw("same");

	TLegend *legend5 = new TLegend(0.5,0.5,0.9,0.7);
	TLegendEntry *leg51 = legend5->AddEntry("dimuons_invmass_mother2","di-muon invariant mass from #mu #mu #gamma","f");
  	leg51->SetFillColor(kRed);	
	TLegendEntry *leg52 = legend5->AddEntry("eta_invmass_mother1","di-muon invariant mass from #mu #mu","f");
  	leg52->SetFillColor(kBlue);	
	TLegendEntry *leg53 = legend5->AddEntry("pi_invmass","misID #pi invariant mass","f");
  	leg53->SetFillColor(kGreen);
	legend5->Draw("same");

	dimuons_invmass_mother2->SetTitle("Reconstructed di-muon invariant mass");
	dimuons_invmass_mother2 -> GetXaxis()-> SetTitle("m_{#mu #mu} (GeV)");

	c1->Modified();
	c1->Update();
	c1->Print("project4452_dimuon_invmass123.pdf","pdf");

	dimuons_invmass_cuts_mother2->SetLineColor(kRed);
	dimuons_invmass_cuts_mother2->SetFillStyle(1001);
	dimuons_invmass_cuts_mother2->Draw();
	eta_invmass_cuts_mother1->SetLineColor(kBlue);
	eta_invmass_cuts_mother1->SetFillStyle(1001);
	eta_invmass_cuts_mother1->Draw("same");
	pi_invmass_cuts->SetLineColor(kGreen);
	pi_invmass_cuts->SetFillStyle(1001);
	pi_invmass_cuts->Draw("same");

	TLegend *legend6 = new TLegend(0.5,0.5,0.9,0.7);	
	TLegendEntry *leg62 = legend6->AddEntry("dimuons_invmass_cuts_mother2","#eta invariant mass from #mu #mu #gamma","f");
  	leg62->SetFillColor(kRed);
	TLegendEntry *leg61 = legend6->AddEntry("eta_invmass_cuts_mother1","#eta invariant mass from #mu #mu","f");
  	leg61->SetFillColor(kBlue);	
	TLegendEntry *leg63 = legend6->AddEntry("pi_invmass_cuts","misID #pi invariant mass","f");
  	leg63->SetFillColor(kGreen);
	legend6->Draw("same");

	eta_invmass_cuts_mother1->SetTitle("Reconstructed di-muon invariant mass after all cuts");
	dimuons_invmass_cuts_mother2 -> GetXaxis()-> SetTitle("m_{#mu #mu} (GeV)");

	c1->Modified();
	c1->Update();
	c1->Print("project4452_dimuon_invmass_cuts123.pdf","pdf");


	dimuons_invmass_c3_mother2->SetLineColor(kRed);
	dimuons_invmass_c3_mother2->SetFillStyle(1001);
	dimuons_invmass_c3_mother2->Draw();
	eta_invmass_c3_mother1->SetLineColor(kBlue);
	eta_invmass_c3_mother1->SetFillStyle(1001);
	eta_invmass_c3_mother1->Draw("same");
	pi_invmass_c3->SetLineColor(kGreen);
	pi_invmass_c3->SetFillStyle(1001);
	pi_invmass_c3->Draw("same");

	TLegend *legend7 = new TLegend(0.5,0.1,0.9,0.3);	
	TLegendEntry *leg72 = legend7->AddEntry("dimuons_invmass_c3_mother2","#eta invariant mass from #mu #mu #gamma","f");
  	leg72->SetFillColor(kRed);
	TLegendEntry *leg71 = legend7->AddEntry("eta_invmass_c3_mother1","#eta invariant mass from #mu #mu","f");
  	leg71->SetFillColor(kBlue);	
	TLegendEntry *leg73 = legend7->AddEntry("pi_invmass_c3","misID #pi invariant mass","f");
  	leg73->SetFillColor(kGreen);
	legend7->Draw("same");

	eta_invmass_c3_mother1->SetTitle("Reconstructed di-muon invariant mass after p_{t}, p & pseudorapidity cuts");
	dimuons_invmass_c3_mother2 -> GetXaxis()-> SetTitle("m_{#mu #mu} (GeV)");

	c1->Modified();
	c1->Update();
	c1->Print("project4452_dimuon_invmass_c3123.pdf","pdf");



	eta_invmass_bkg2->SetLineColor(kRed);
	eta_invmass_bkg2->SetFillStyle(1001);
	eta_invmass_bkg2->Draw();
	eta_invmass_bkg1->SetLineColor(kBlue);
	eta_invmass_bkg1->SetFillStyle(1001);
	eta_invmass_bkg1->Draw("same");
	pi_invmass->SetLineColor(kGreen);
	pi_invmass->SetFillStyle(1001);
	pi_invmass->Draw("same");

	TLegend *legend11 = new TLegend(0.5,0.5,0.9,0.7);
	TLegendEntry *leg112 = legend11->AddEntry("eta_invmass_bkg2","#eta invariant mass from #mu #mu #gamma","f");
  	leg112->SetFillColor(kRed);	
	TLegendEntry *leg111 = legend11->AddEntry("eta_invmass_bkg1","#eta invariant mass from #mu #mu","f");
  	leg111->SetFillColor(kBlue);	
	TLegendEntry *leg113 = legend11->AddEntry("pi_invmass","misID #pi invariant mass","f");
  	leg113->SetFillColor(kGreen);
	legend11->Draw("same");

	eta_invmass_bkg2->SetTitle("Reconstructed #eta invariant mass (bkg)");

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etapi_invmass_bkg123.pdf","pdf");


	eta_invmass_cuts_bkg2->SetLineColor(kRed);
	eta_invmass_cuts_bkg2->SetFillStyle(1001);
	eta_invmass_cuts_bkg2->Draw();
	eta_invmass_cuts_bkg1->SetLineColor(kBlue);
	eta_invmass_cuts_bkg1->SetFillStyle(1001);
	eta_invmass_cuts_bkg1->Draw("same");
	pi_invmass_cuts->SetLineColor(kGreen);
	pi_invmass_cuts->SetFillStyle(1001);
	pi_invmass_cuts->Draw("same");

	TLegend *legend12 = new TLegend(0.5,0.5,0.9,0.7);	
	TLegendEntry *leg122 = legend12->AddEntry("eta_invmass_cuts_bkg2","#eta invariant mass from #mu #mu #gamma","f");
  	leg122->SetFillColor(kRed);
	TLegendEntry *leg121 = legend12->AddEntry("eta_invmass_cuts_bkg1","#eta invariant mass from #mu #mu","f");
  	leg121->SetFillColor(kBlue);	
	TLegendEntry *leg123 = legend12->AddEntry("pi_invmass_cuts","misID #pi invariant mass","f");
  	leg123->SetFillColor(kGreen);
	legend12->Draw("same");

	eta_invmass_cuts_bkg2->SetAxisRange(1e-2, 1e2,"Y");
	eta_invmass_cuts_bkg2->SetTitle("Reconstructed #eta invariant mass after all cuts (bkg)");

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etapi_invmass_bkg_cuts123.pdf","pdf");


	eta_invmass_c3_bkg2->SetLineColor(kRed);
	eta_invmass_c3_bkg2->SetFillStyle(1001);
	eta_invmass_c3_bkg2->Draw();
	eta_invmass_c3_bkg1->SetLineColor(kBlue);
	eta_invmass_c3_bkg1->SetFillStyle(1001);
	eta_invmass_c3_bkg1->Draw("same");
	pi_invmass_c3->SetLineColor(kGreen);
	pi_invmass_c3->SetFillStyle(1001);
	pi_invmass_c3->Draw("same");

	TLegend *legend13 = new TLegend(0.5,0.4,0.9,0.6);	
	TLegendEntry *leg132 = legend13->AddEntry("eta_invmass_c3_bkg2","#eta invariant mass from #mu #mu #gamma","f");
  	leg132->SetFillColor(kRed);
	TLegendEntry *leg131 = legend13->AddEntry("eta_invmass_c3_bkg1","#eta invariant mass from #mu #mu","f");
  	leg131->SetFillColor(kBlue);	
	TLegendEntry *leg133 = legend13->AddEntry("pi_invmass_c3","misID #pi invariant mass","f");
  	leg133->SetFillColor(kGreen);
	legend13->Draw("same");

	eta_invmass_c3_bkg2->SetTitle("Reconstructed #eta invariant mass after p_t, p & pseudorapidity cuts (bkg)");

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etapi_invmass_bkg_c3123.pdf","pdf");


	dimuons_invmass_bkg2->SetLineColor(kRed);
	dimuons_invmass_bkg2->SetFillStyle(1001);
	dimuons_invmass_bkg2->Draw();
	eta_invmass_bkg1->SetLineColor(kBlue);
	eta_invmass_bkg1->SetFillStyle(1001);
	eta_invmass_bkg1->Draw("same");
	pi_invmass->SetLineColor(kGreen);
	pi_invmass->SetFillStyle(1001);
	pi_invmass->Draw("same");

	TLegend *legend8 = new TLegend(0.5,0.5,0.9,0.7);
	TLegendEntry *leg81 = legend8->AddEntry("dimuons_invmass_bkg2","di-muon invariant mass from #mu #mu #gamma","f");
  	leg81->SetFillColor(kRed);	
	TLegendEntry *leg82 = legend8->AddEntry("eta_invmass_bkg1","di-muon invariant mass from #mu #mu","f");
  	leg82->SetFillColor(kBlue);	
	TLegendEntry *leg83 = legend8->AddEntry("pi_invmass","misID #pi invariant mass","f");
  	leg83->SetFillColor(kGreen);
	legend8->Draw("same");

	dimuons_invmass_bkg2->SetTitle("Reconstructed di-muon invariant mass (bkg)");
	dimuons_invmass_bkg2 -> GetXaxis()-> SetTitle("m_{#mu #mu} (GeV)");

	c1->Modified();
	c1->Update();
	c1->Print("project4452_dimuon_invmass_bkg123.pdf","pdf");

	dimuons_invmass_cuts_bkg2->SetLineColor(kRed);
	dimuons_invmass_cuts_bkg2->SetFillStyle(1001);
	dimuons_invmass_cuts_bkg2->Draw();
	eta_invmass_cuts_bkg1->SetLineColor(kBlue);
	eta_invmass_cuts_bkg1->SetFillStyle(1001);
	eta_invmass_cuts_bkg1->Draw("same");
	pi_invmass_cuts->SetLineColor(kGreen);
	pi_invmass_cuts->SetFillStyle(1001);
	pi_invmass_cuts->Draw("same");

	TLegend *legend9 = new TLegend(0.5,0.5,0.9,0.7);	
	TLegendEntry *leg92 = legend9->AddEntry("dimuons_invmass_cuts_bkg2","#eta invariant mass from #mu #mu #gamma","f");
  	leg92->SetFillColor(kRed);
	TLegendEntry *leg91 = legend9->AddEntry("eta_invmass_cuts_bkg1","#eta invariant mass from #mu #mu","f");
  	leg91->SetFillColor(kBlue);	
	TLegendEntry *leg93 = legend9->AddEntry("pi_invmass_cuts","misID #pi invariant mass","f");
  	leg93->SetFillColor(kGreen);
	legend9->Draw("same");

	eta_invmass_cuts_bkg2->SetTitle("Reconstructed di-muon invariant mass after all cuts (bkg)");
	dimuons_invmass_cuts_bkg2 -> GetXaxis()-> SetTitle("m_{#mu #mu} (GeV)");

	c1->Modified();
	c1->Update();
	c1->Print("project4452_dimuon_invmass_bkg_cuts123.pdf","pdf");


	dimuons_invmass_c3_bkg2->SetLineColor(kRed);
	dimuons_invmass_c3_bkg2->SetFillStyle(1001);
	dimuons_invmass_c3_bkg2->Draw();
	eta_invmass_c3_bkg1->SetLineColor(kBlue);
	eta_invmass_c3_bkg1->SetFillStyle(1001);
	eta_invmass_c3_bkg1->Draw("same");
	pi_invmass_c3->SetLineColor(kGreen);
	pi_invmass_c3->SetFillStyle(1001);
	pi_invmass_c3->Draw("same");

	TLegend *legend10 = new TLegend(0.5,0.55,0.9,0.75);	
	TLegendEntry *leg102 = legend10->AddEntry("dimuons_invmass_c3_bkg2","#eta invariant mass from #mu #mu #gamma","f");
  	leg102->SetFillColor(kRed);
	TLegendEntry *leg101 = legend10->AddEntry("eta_invmass_c3_bkg1","#eta invariant mass from #mu #mu","f");
  	leg101->SetFillColor(kBlue);	
	TLegendEntry *leg103 = legend10->AddEntry("pi_invmass_c3","misID #pi invariant mass","f");
  	leg103->SetFillColor(kGreen);
	legend10->Draw("same");

	dimuons_invmass_c3_bkg2->SetAxisRange(1e-2, 1e3,"Y");

	eta_invmass_c3_bkg2->SetTitle("Reconstructed di-muon invariant mass after p_{t}, p & pseudorapidity cuts (bkg)");
	dimuons_invmass_c3_bkg2 -> GetXaxis()-> SetTitle("m_{#mu #mu} (GeV)");

	c1->Modified();
	c1->Update();
	c1->Print("project4452_dimuon_invmass_bkg_c3123.pdf","pdf");





	eta_invmass_mother2->SetLineColor(kRed);
	eta_invmass_mother2->SetFillStyle(1001);
	eta_invmass_mother2->Draw();
	eta_invmass_mother1->SetLineColor(kBlue);
	eta_invmass_mother1->SetFillStyle(1001);
	eta_invmass_mother1->Draw("same");
	pi_invmass->SetLineColor(kGreen);
	pi_invmass->SetFillStyle(1001);
	pi_invmass->Draw("same");

	TLegend *legend1 = new TLegend(0.5,0.5,0.9,0.7);
	TLegendEntry *leg12 = legend1->AddEntry("eta_invmass_mother2","#eta invariant mass from #mu #mu #gamma","f");
  	leg12->SetFillColor(kRed);	
	TLegendEntry *leg11 = legend1->AddEntry("eta_invmass_mother1","#eta invariant mass from #mu #mu","f");
  	leg11->SetFillColor(kBlue);	
	TLegendEntry *leg13 = legend1->AddEntry("pi_invmass","misID #pi invariant mass","f");
  	leg13->SetFillColor(kGreen);
	legend1->Draw("same");

	eta_invmass_mother2->SetTitle("Reconstructed #eta invariant mass");

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etapi_invmass123.pdf","pdf");

	gPad->SetLogx();

	c1->Modified();
	c1->Update();
	c1->Print("project4452_etainvmass123.pdf","pdf");


	file_out->Write();

	delete file_in1;
	delete file_in2;
	delete file_out;

	delete v;

	// Done.
	return 0;
}


