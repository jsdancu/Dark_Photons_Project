//24/11/2017
//Program reading in the three trees generated in gs_project25.cc and performing analysis on the reconstructed invariant mass of eta

#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "cmath"
#include <iostream>
#include <vector>

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

void analyze_event(vect *v, TH1D *eta_invmass, TH1D *muon_invmass_mother, TH1D *eta_invmass_mothers, TH1D *muon_pt_mother, TH1D *muon_pt_bkg, TH1D *muon_ptcut_mother, TH1D *muon_ptetacut_mother, TH1D *muon_p_mother, TH1D *muon_p_bkg, TH1D *muon_pcut_mother, TH1D *muon_petacut_mother, TH1D *muon_eta_mother, TH1D *mu_number_event){

	double inv_mass;

	Long64_t nentries = v->index.size();

	Long64_t mu_number = 0;
	Long64_t antimu_number = 0;

	bool antimu = false;

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

					//Reconstructing invariant mass of all muon-anti-muon combinations
					inv_mass = invmass(v, i, j);

					//Plot histogram of all muon-anti-muon combinations
					eta_invmass->Fill(inv_mass);

					if((v->mother1[i] == v->mother1[j]) && (v->mother2[i] == v->mother2[j])){

						//If the muon-antimuon pair come from the same mother particle plot them in a separate histogram
						muon_invmass_mother->Fill(inv_mass);

						if(v->motherid1[i] == 221){

							eta_invmass_mothers->Fill(inv_mass);

							double pt_mother1 = pt(v, i);
							muon_pt_mother->Fill(pt_mother1);

							double pt_mother2 = pt(v, j);
							muon_pt_mother->Fill(pt_mother2);

							double p_mother1 = p(v, i);
							muon_p_mother->Fill(p_mother1);

							double p_mother2 = p(v, j);
							muon_p_mother->Fill(p_mother2);

							double eta_mother1 = eta(v, i);
							muon_eta_mother->Fill(eta_mother1);

							double eta_mother2 = eta(v, j);
							muon_eta_mother->Fill(eta_mother2);

							if((pt_mother1>0.5) && (p_mother1>10.0) && (pt_mother2>0.5) && (p_mother2>10.0))
							{

								muon_ptcut_mother->Fill(pt_mother1);
								muon_ptcut_mother->Fill(pt_mother2);
								muon_pcut_mother->Fill(p_mother1);
								muon_pcut_mother->Fill(p_mother2);

								if((eta_mother1>2.0) && (eta_mother1<4.5) && (eta_mother2>2.0) && (eta_mother2<4.5))
								{
									muon_ptetacut_mother->Fill(pt_mother1);
									muon_ptetacut_mother->Fill(pt_mother2);
									muon_petacut_mother->Fill(p_mother1);
									muon_petacut_mother->Fill(p_mother2);
								}

							}

						}

					}
					else if(v->mother1[i] != v->mother1[j]){

						double pt_bkg1 = pt(v, i);
						muon_pt_bkg->Fill(pt_bkg1);

						double pt_bkg2 = pt(v, j);
						muon_pt_bkg->Fill(pt_bkg2);

						double p_bkg1 = p(v, i);
						muon_p_bkg->Fill(p_bkg1);

						double p_bkg2 = p(v, j);
						muon_p_bkg->Fill(p_bkg2);

					}

				}

			}
			antimu = true;
		}

	}

	//std::cout<<"muon number = "<<mu_number<<setw(5)<<"anti muon = "<<antimu_number<<std::endl;
	double total_mu_antimu = mu_number+antimu_number;
	if (v->index.size() > 0) {mu_number_event->SetBinContent((v->index[1])+1,total_mu_antimu);}
	//if (v->index.size() > 0) std::cout<<v->index[0]<<std::endl;
	//std::cout<<v << "    " << v->index.size() <<std::endl;

}

void analyze_event_pi(vect *v, TH1D *pi_invmass, TH1D *pi_pt, TH1D *pi_p, TH1D *pi_eta, TH1D *pi_ptcut, TH1D *pi_pcut, TH1D *pi_number_event){

	double inv_mass;

	Long64_t nentries = v->index.size();

	Long64_t pi_number_p = 0;
	Long64_t pi_number_n = 0;

	bool pi_n = false;

	for(Long64_t i=0;i<nentries;i++){ 

		//if entry is pi+
		if(v->id[i] == 211){

			++pi_number_p;

			//Loop through the event and find an anti-muon
			for(Long64_t j=0;j<nentries;j++){

				//if entry is anti-muon
				if(v->id[j] == -211){

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

					if((pt_pp>0.5) && (p_pp>10.0) && (pt_pn>0.5) && (p_pn>10.0) && (eta_pp>2.0) && (eta_pp<4.5) && (eta_pn>2.0) && (eta_pn<4.5))
					{

						pi_ptcut->Fill(pt_pp);
						pi_ptcut->Fill(pt_pn);
						pi_pcut->Fill(p_pp);
						pi_pcut->Fill(p_pn);


					}


				}

			}
			pi_n = true;
		}

	}

	//std::cout<<"muon number = "<<mu_number<<setw(5)<<"anti muon = "<<antimu_number<<std::endl;
	double total_pi_pn = pi_number_p+pi_number_n;
	if (v->index.size() > 0) {pi_number_event->SetBinContent((v->index[1])+1,total_pi_pn);}
	//if (v->index.size() > 0) std::cout<<v->index[0]<<std::endl;
	//std::cout<<v << "    " << v->index.size() <<std::endl;

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
	TFile  *file_in = new TFile("gs_project25.root", "READ");

	// Open the output TFile 
	TFile  *file_out = new TFile("gs_project45.root", "recreate");

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

	TH1D *eta_invmass = new TH1D("eta_invmass","Reconstructed #eta invariant mass from muon pairs", 300, 0.0, 30.0);
    	eta_invmass -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *muon_invmass_mother = new TH1D("muon_invmass_mother","Reconstructed muon pair invariant mass with same mother", 3000, 0.0, 30.0);
    	muon_invmass_mother -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_mother = new TH1D("eta_invmass_mother","Reconstructed #eta invariant mass from muon pairs with same mother", 100, 0.0, 1.0);
    	eta_invmass_mother -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *muon_pt_mother = new TH1D("muon_pt_mother","#mu^{-} and #mu^{+} p_{t} distribution with same mother (#eta)", 500, 0.0, 5.0);
    	muon_pt_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *muon_pt_bkg = new TH1D("muon_pt_bkg","#mu^{-} and #mu^{+} p_{t} distribution for background", 500, 0.0, 5.0);
    	muon_pt_bkg -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *muon_ptcut_mother = new TH1D("muon_ptcut_mother","#mu^{-} and #mu^{+} p_{t} distribution after cuts (with same mother (#eta))", 500, 0.0, 5.0);
    	muon_ptcut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *muon_ptetacut_mother = new TH1D("muon_ptetacut_mother","#mu^{-} and #mu^{+} p_{t} distribution after cuts including pseudorapidity (with same mother (#eta))", 500, 0.0, 5.0);
    	muon_ptetacut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *muon_p_mother = new TH1D("muon_p_mother","#mu^{-} and #mu^{+} p distribution with same mother (#eta)", 10000, 0.0, 100.0);
    	muon_p_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *muon_p_bkg = new TH1D("muon_p_bkg","#mu^{-} and #mu^{+} p distribution for background", 10000, 0.0, 100.0);
    	muon_p_bkg -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *muon_pcut_mother = new TH1D("muon_pcut_mother","#mu^{-} and #mu^{+} p distribution after cuts (with same mother (#eta))", 10000, 0.0, 100.0);
    	muon_pcut_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *muon_petacut_mother = new TH1D("muon_petacut_mother","#mu^{-} and #mu^{+} p distribution after cuts including pseudorapidity (with same mother (#eta))", 10000, 0.0, 100.0);
    	muon_petacut_mother -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *muon_eta_mother = new TH1D("muon_eta_mother","#mu^{-} and #mu^{+} pseudorapidity distribution with same mother (#eta)", 45, 0.0, 4.5);
    	muon_eta_mother -> GetXaxis()-> SetTitle("#eta");
	
	TH1D *mu_number_event = new TH1D("mu_number_event","Muon-antimuon number per event", 10, 0.0, 10.0);
    	mu_number_event -> GetXaxis()-> SetTitle("event index");
	mu_number_event -> GetYaxis()-> SetTitle("number of #mu^{#pm}");

	TH1D *pi_invmass = new TH1D("pi_invmass","Reconstructed #eta invariant mass from misID #pi^{#pm}", 300, 0.0, 30.0);
    	pi_invmass -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *pi_pt = new TH1D("pi_pt","#pi^{-} and #pi^{+} p_{t} distribution", 500, 0.0, 5.0);
    	pi_pt -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *pi_ptcut = new TH1D("pi_ptcut","#pi^{-} and #pi^{+} p_{t} distribution after cuts", 500, 0.0, 5.0);
    	pi_ptcut -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *pi_p = new TH1D("pi_p","#pi^{-} and #pi^{+} p distribution", 1000, 0.0, 100.0);
    	pi_p -> GetXaxis()-> SetTitle("p_{t} (GeV)");

	TH1D *pi_pcut = new TH1D("pi_pcut","#pi^{-} and #pi^{+} p distribution after cuts", 1000, 0.0, 100.0);
    	pi_pcut -> GetXaxis()-> SetTitle("p (GeV)");

	TH1D *pi_eta = new TH1D("pi_eta","#pi^{-} and #pi^{+} pseudorapidity distribution", 45, 0.0, 4.5);
    	pi_eta -> GetXaxis()-> SetTitle("#eta");

	TH1D *pi_number_event = new TH1D("pi_number_event","#pi^{#pm} number per event", 10, 0.0, 10.0);
    	pi_number_event -> GetXaxis()-> SetTitle("event index");
	pi_number_event -> GetYaxis()-> SetTitle("number of #pi^{#pm}");

	vect* v= new vect;//declares a vector in which we are going to store all the particles from the same event
	double prev_index = 0.0;//saving the event index of the previous entry of the tree

	//Loop through the entries of tree T2
	Long64_t nentries = T2->GetEntries();
std::cout<<"Total number of entries T2: "<<nentries<<std::endl;
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

			analyze_event(v, eta_invmass, muon_invmass_mother, eta_invmass_mother, muon_pt_mother, muon_pt_bkg, muon_ptcut_mother, muon_ptetacut_mother, muon_p_mother, muon_p_bkg, muon_pcut_mother, muon_petacut_mother, muon_eta_mother, mu_number_event);

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

			analyze_event(v, eta_invmass, muon_invmass_mother, eta_invmass_mother, muon_pt_mother, muon_pt_bkg, muon_ptcut_mother, muon_ptetacut_mother, muon_p_mother, muon_p_bkg, muon_pcut_mother, muon_petacut_mother, muon_eta_mother, mu_number_event);

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

			analyze_event_pi(v, pi_invmass, pi_pt, pi_p, pi_eta, pi_ptcut, pi_pcut, pi_number_event);

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

			analyze_event_pi(v, pi_invmass, pi_pt, pi_p, pi_eta, pi_ptcut, pi_pcut, pi_number_event);

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

	muon_pt_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_mu_pt_mother1.pdf","pdf");

	muon_pt_bkg->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_mu_pt_bkg1.pdf","pdf");

	muon_ptcut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_mu_ptcut_mother1.pdf","pdf");

	muon_ptetacut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_mu_ptetacut_mother1.pdf","pdf");

	muon_p_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_mu_p_mother1.pdf","pdf");

	muon_p_bkg->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_mu_p_bkg1.pdf","pdf");

	muon_pcut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_mu_pcut_mother1.pdf","pdf");

	muon_petacut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_mu_petacut_mother1.pdf","pdf");

	muon_eta_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_mu_pseudorap_mother1.pdf","pdf");

	mu_number_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_mu_number_event1.pdf","pdf");

	pi_pt->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_pi_pt1.pdf","pdf");

	pi_p->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_pi_p1.pdf","pdf");

	pi_eta->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_pi_eta1.pdf","pdf");

	pi_ptcut->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_pi_ptcut1.pdf","pdf");

	pi_pcut->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_pi_pcut1.pdf","pdf");

	pi_number_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_pi_number_event1.pdf","pdf");

	gPad->SetLogy();

	eta_invmass->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_etainvmass1.pdf","pdf");

	muon_invmass_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_muoninvmass_mother1.pdf","pdf");

	eta_invmass_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_etainvmass_mother1.pdf","pdf");
 
	pi_invmass->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project45_pi_invmass1.pdf","pdf");

	//muon_pt_mother->SetFillColor(kBlue);
	muon_pt_mother->SetLineColor(kBlue);
	muon_pt_mother->SetFillStyle(1001);
	muon_pt_mother->Draw();
	//muon_pt_bkg->SetFillColor(kYellow);
	muon_pt_bkg->SetLineColor(kYellow);
	muon_pt_bkg->SetFillStyle(1001);
	muon_pt_bkg->Draw("same");
	//muon_ptcut_mother->SetFillColor(kGreen);
	muon_ptcut_mother->SetLineColor(kGreen);
	muon_ptcut_mother->SetFillStyle(1001);
	muon_ptcut_mother->Draw("same");
	//muon_ptetacut_mother->SetFillColor(kRed);
	muon_ptetacut_mother->SetLineColor(kRed);
	muon_ptetacut_mother->SetFillStyle(1001);
	muon_ptetacut_mother->Draw("same");

	TLegend *legend2 = new TLegend(0.5,0.5,0.9,0.7);
	TLegendEntry *leg3 = legend2->AddEntry("muon_pt_mother","p_{t} distribution before cuts","f");
  	leg3->SetFillColor(kBlue);	
	TLegendEntry *leg32 = legend2->AddEntry("muon_pt_bkg","p_{t} distribution for background","f");
  	leg32->SetFillColor(kYellow);
	TLegendEntry *leg33 = legend2->AddEntry("muon_ptcut_mother","p_{t} distribution after cuts","f");
  	leg33->SetFillColor(kGreen);
	TLegendEntry *leg31 = legend2->AddEntry("muon_ptetacut_mother","p_{t} distribution after cuts including pseudorapidity","f");
  	leg31->SetFillColor(kRed);
	legend2->Draw("same");

	c1->Modified();
	c1->Update();
	c1->Print("project45_mu_pt_mother12.pdf","pdf");

	pi_pt->SetLineColor(kBlue);
	pi_pt->SetFillStyle(1001);
	pi_pt->Draw();
	pi_ptcut->SetLineColor(kRed);
	pi_ptcut->SetFillStyle(1001);
	pi_ptcut->Draw("same");

	TLegend *legend4 = new TLegend(0.5,0.5,0.9,0.7);
	TLegendEntry *leg41 = legend4->AddEntry("pi_pt","#pi^{#pm} p_{t} distribution before cuts","f");
  	leg41->SetFillColor(kBlue);	
	TLegendEntry *leg42 = legend4->AddEntry("pi_ptcut","#pi^{#pm} p_{t} distribution after cuts","f");
  	leg42->SetFillColor(kRed);
	legend4->Draw("same");

	c1->Modified();
	c1->Update();
	c1->Print("project45_pi_pt12.pdf","pdf");

	//muon_p_mother->SetFillColor(kBlue);
	muon_p_mother->SetLineColor(kBlue);
	muon_p_mother->SetFillStyle(1001);
	muon_p_mother->Draw();
	//muon_p_bkg->SetFillColor(kYellow);
	muon_p_bkg->SetLineColor(kYellow);
	muon_p_bkg->SetFillStyle(1001);
	muon_p_bkg->Draw("same");
	//muon_pcut_mother->SetFillColor(kGreen);
	muon_pcut_mother->SetLineColor(kGreen);
	muon_pcut_mother->SetFillStyle(1001);
	muon_pcut_mother->Draw("same");
	//muon_petacut_mother->SetFillColor(kRed);
	muon_petacut_mother->SetLineColor(kRed);
	muon_petacut_mother->SetFillStyle(1001);
	muon_petacut_mother->Draw("same");

	TLegend *legend3 = new TLegend(0.5,0.5,0.9,0.7);
	TLegendEntry *leg5 = legend3->AddEntry("muon_p_mother","p distribution before cuts","f");
  	leg3->SetLineColor(kBlue);	
	TLegendEntry *leg62 = legend3->AddEntry("muon_p_bkg","p distribution for background","f");
  	leg62->SetLineColor(kYellow);
	TLegendEntry *leg6 = legend3->AddEntry("muon_pcut_mother","p distribution after cuts","f");
  	leg6->SetLineColor(kGreen);
	TLegendEntry *leg61 = legend3->AddEntry("muon_petacut_mother","p distribution after cuts including pseudorapidity","f");
  	leg61->SetLineColor(kRed);
	legend3->Draw("same");

	c1->Modified();
	c1->Update();
	c1->Print("project45_mu_p_mother12.pdf","pdf");

	pi_p->SetLineColor(kBlue);
	pi_p->SetFillStyle(1001);
	pi_p->Draw();
	pi_pcut->SetLineColor(kRed);
	pi_pcut->SetFillStyle(1001);
	pi_pcut->Draw("same");

	TLegend *legend5 = new TLegend(0.5,0.5,0.9,0.7);
	TLegendEntry *leg51 = legend5->AddEntry("pi_p","#pi^{#pm} p distribution before cuts","f");
  	leg51->SetFillColor(kBlue);	
	TLegendEntry *leg52 = legend5->AddEntry("pi_pcut","#pi^{#pm} p distribution after cuts","f");
  	leg52->SetFillColor(kRed);
	legend5->Draw("same");

	c1->Modified();
	c1->Update();
	c1->Print("project45_pi_p12.pdf","pdf");

	//eta_invmass->SetFillColor(kBlue);
	eta_invmass->SetFillStyle(1001);
	eta_invmass->Draw();
	//eta_invmass_mother->SetFillColor(kGreen);
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
	c1->Print("project45_etainvmass12.pdf","pdf");

	file_out->Write();

	delete file_in;
	delete file_out;

	// Done.
	return 0;
}


