//16/11/2017
//Program reading in the two trees generated in gs_project23.cc and performing analysis on the reconstructed invariant mass of eta

//#include "treeio.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "cmath"
#include <iostream>
#include <vector>

struct vect{
	
	std::vector<double> index, id, energy, mass, px, py, pz, mother1, mother2, motherid1, motherid2;
	
};

double invmass(vect *v, Long64_t i, Long64_t j){

	//std::cout<<"i = "<<i<<"    "<<"j = "<<j<<"    "<< "energy i = "<<v->energy[i]<<"    "<<"energy j = "<<v->energy[j]<<std::endl;

	double inv_mass=sqrt(pow(v->energy[i]+v->energy[j],2)-pow(v->px[i]+v->px[j],2)-pow(v->py[i]+v->py[j],2)-pow(v->pz[i]+v->pz[j],2));

	//std::cout<<"invariant mass = "<<inv_mass<<std::endl;

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

void analyze_event(vect *v, TH1D *eta_invmass, TH1D *muon_invmass_mother, TH1D *eta_invmass_mothers, TH1D *muon_pt_mother, TH1D *muon_ptcut_mother, TH1D *mu_number_event){

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

//std::cout<<"i = "<<i<<"    "<<"j = "<<j<<"    "<< "energy i = "<<v->energy[i]<<"    "<<"energy j = "<<v->energy[j]<<std::endl;

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
							//muon_p_mother->Fill(p_mother1);

							if((pt_mother1>0.5) && (p_mother1>10.0))
							{
							muon_ptcut_mother->Fill(pt_mother1);
							}

							double p_mother2 = p(v, j);
							//muon_p_mother->Fill(p_mother2);

							if((pt_mother2>0.5) && (p_mother2>10.0))
							{
							muon_ptcut_mother->Fill(pt_mother2);
							}

						}

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
	std::cout<<v << "    " << v->index.size() <<std::endl;

}

int main() {

	// Open the input TFile 
	TFile  *file_in = new TFile("gs_project23.root", "READ");

	// Open the output TFile 
	TFile  *file_out = new TFile("gs_project43.root", "recreate");

	//Get trees from file
	TTree *T1 = (TTree*)file_in->Get("T1");
	TTree *T2 = (TTree*)file_in->Get("T2");

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


	TH1D *eta_invmass = new TH1D("eta_invmass","Reconstructed eta invariant mass from muon pairs", 100, 0.0, 1.0);
    	eta_invmass -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *muon_invmass_mother = new TH1D("muon_invmass_mother","Reconstructed muon pair invariant mass with same mother", 100, 0.0, 1.0);
    	muon_invmass_mother -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *eta_invmass_mother = new TH1D("eta_invmass_mother","Reconstructed eta invariant mass from muon pairs with same mother", 100, 0.0, 1.0);
    	eta_invmass_mother -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *muon_pt_mother = new TH1D("muon_pt_mother","#mu^{-} and #mu^{+} p_{t} distribution with same mother (#eta)", 100, 0.0, 1.0);
    	muon_pt_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");
	//muon_pt_mother -> GetYaxis()-> SetLog();

	TH1D *muon_ptcut_mother = new TH1D("muon_ptcut_mother","#mu^{-} and #mu^{+} p_{t} distribution after cuts (with same mother (#eta))", 300, 0.0, 3.0);
    	muon_ptcut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");
	//muon_ptcut_mother -> GetYaxis()-> SetLog();

	TH1D *mu_number_event = new TH1D("mu_number_event","Muon-antimuon number per event", 10, 0.0, 10.0);
    	mu_number_event -> GetXaxis()-> SetTitle("event index");
	mu_number_event -> GetYaxis()-> SetTitle("number of #mu^{#pm}");

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

			analyze_event(v, eta_invmass, muon_invmass_mother, eta_invmass_mother, muon_pt_mother, muon_ptcut_mother, mu_number_event);

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

			analyze_event(v, eta_invmass, muon_invmass_mother, eta_invmass_mother, muon_pt_mother, muon_ptcut_mother, mu_number_event);

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

	eta_invmass->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project43_etainvmass1.pdf","pdf");

	muon_invmass_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project43_muoninvmass_mother1.pdf","pdf");

	eta_invmass_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project43_etainvmass_mother1.pdf","pdf");

	//c1->SetLogy();
	muon_pt_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project43_mu_pt_mother1.pdf","pdf");

	muon_ptcut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project43_mu_ptcut_mother1.pdf","pdf");

	mu_number_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project43_mu_number_event.pdf","pdf");

	file_out->Write();

	delete file_in;
	delete file_out;

	// Done.
	return 0;
}


