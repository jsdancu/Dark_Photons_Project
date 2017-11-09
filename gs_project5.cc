//12/10/2017
//Program reading in the tree generated in gs_project2.cc and performing analysis on the reconstructed invariant mass
//credits to Phil for the reading the tree from file part
//credits to Nigel for analysis procedure suggestions
#include "treeio.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
//#include "cmath"
#include <iostream>
#include <vector>

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

void analyze_event(vect *v, TH1F *eta_invmass, TH1F *muongamma_invmass_mother, TH1F *eta_invmass_mother, TH1F *muon_invmass_mother, TH1F *muon_pt_mother, TH1F *muon_ptcut_mother, TH1F *gamma_pt_mother, TH1F *mu_number_event, TH1F *gamma_number_event, TH1F *total_number_event){

	double inv_mass;

	Long64_t nentries = v->index.size();

//std::cout<<"event size = "<<nentries<<std::endl;

	Long64_t mu_number = 0;
	Long64_t antimu_number = 0;
	Long64_t gamma_number = 0;

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

							if(gamma == false){
								++gamma_number;
							}

							//Reconstructing invariant mass of all muon-antimuon-photon combinations
							inv_mass = invmass(v, i, j, k);

							//Plot histogram of all muon-antimuon-photon combinations
							eta_invmass->Fill(inv_mass);

							if((v->mother1[i] == v->mother1[j]) && (v->mother1[i] == v->mother1[k]) && (v->mother2[i] == v->mother2[j]) && (v->mother2[i] == v->mother2[k])){

								//If the muon-antimuon-photon pair come from the same mother particle plot them in a separate histogram

								if((v->mother1[i] == v->mother1[j]) && (v->mother2[i] == v->mother2[j])){

									//If the muon-antimuon pair come from the same mother particle plot them in a separate histogram
									muongamma_invmass_mother->Fill(inv_mass);

									if(v->motherid1[i] == 221){

										eta_invmass_mother->Fill(inv_mass);

										double mu_invmass = invmass2(v, i, j);
										muon_invmass_mother->Fill(mu_invmass);

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

										double pt_mother3 = pt(v, k);
										gamma_pt_mother->Fill(pt_mother3);

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

	//if(int(v->index[0]) % 1000 == 0){std::cout<<"muon number = "<<mu_number<<setw(5)<<"anti muon = "<<antimu_number<<setw(5)<<"photon = "<<gamma_number<<std::endl;}

	if (v->index.size() > 0) 
	{
		double total_mu_antimu = mu_number+antimu_number;
		double total_number = mu_number+antimu_number+gamma_number;
		mu_number_event->SetBinContent((v->index[1])+1,total_mu_antimu);
		gamma_number_event->SetBinContent((v->index[1])+1,gamma_number);
		total_number_event->SetBinContent((v->index[1])+1,total_number);

	}
	//if (v->index.size() > 0) std::cout<<v->index[0]<<std::endl;
	//std::cout<<v << "    " << v->index.size() <<std::endl;

}

int main() {

	// Open the input TFile 
	TFile  *file_in = new TFile("gs_project3.root", "READ");

	// Open the output TFile 
	TFile  *file_out = new TFile("gs_project5.root", "recreate");

	//Get tree from file
	TTree *T = (TTree*)file_in->Get("T");

	double index_var, id_var, energy_var, mass_var, px_var, py_var, pz_var, mother1_var, mother2_var, motherid1_var, motherid2_var;

	//Get branches of the tree
	T->SetBranchAddress("index",&index_var);
   	T->SetBranchAddress("id",&id_var);
	T->SetBranchAddress("energy",&energy_var);
   	T->SetBranchAddress("mass",&mass_var);
	T->SetBranchAddress("px",&px_var);
   	T->SetBranchAddress("py",&py_var);
	T->SetBranchAddress("pz",&pz_var);
	T->SetBranchAddress("mother1",&mother1_var);
	T->SetBranchAddress("mother2",&mother2_var);
	T->SetBranchAddress("motherid1",&motherid1_var);
	T->SetBranchAddress("motherid2",&motherid2_var);

	TH1F *eta_invmass = new TH1F("eta_invmass","Reconstructed eta invariant mass from muon pairs and gamma", 400, 0.0, 40.0);
    	eta_invmass -> GetXaxis()-> SetTitle("m (GeV)");

	TH1F *muongamma_invmass_mother = new TH1F("muongamma_invmass_mother","Reconstructed muon pair + #gamma invariant mass with same mother", 100, 0.0, 1.0);
    	muongamma_invmass_mother -> GetXaxis()-> SetTitle("m (GeV)");

	TH1F *eta_invmass_mother = new TH1F("eta_invmass_mother","Reconstructed invariant mass from muon pairs + #gamma with same mother (#eta)", 100, 0.0, 1.0);
    	eta_invmass_mother -> GetXaxis()-> SetTitle("m (GeV)");

	TH1F *muon_invmass_mother = new TH1F("muon_invmass_mother","Reconstructed muon pair invariant mass with same mother(#eta)", 100, 0.0, 1.0);
    	muon_invmass_mother -> GetXaxis()-> SetTitle("m (GeV)");

	TH1F *muon_pt_mother = new TH1F("muon_pt_mother","#mu^{-} and #mu^{+} p_{t} distribution with same mother (#eta)", 150, 0.0, 1.5);
    	muon_pt_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");
	//muon_pt_mother -> GetYaxis()-> SetLog();

	TH1F *muon_ptcut_mother = new TH1F("muon_ptcut_mother","#mu^{-} and #mu^{+} p_{t} distribution after cuts (with same mother (#eta))", 300, 0.0, 3.0);
    	muon_ptcut_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");
	
	TH1F *gamma_pt_mother = new TH1F("gamma_pt_mother","#gamma p_{t} distribution as #eta decay product", 150, 0.0, 1.5);
    	muon_pt_mother -> GetXaxis()-> SetTitle("p_{t} (GeV)");
	//muon_pt_mother -> GetYaxis()-> SetLog();

	TH1F *mu_number_event = new TH1F("mu_number_event","Muon-antimuon number per event", 100000, 0.0, 100000.0);
    	mu_number_event -> GetXaxis()-> SetTitle("event index");
	mu_number_event -> GetYaxis()-> SetTitle("number of #mu^{#pm}");

	TH1F *gamma_number_event = new TH1F("gamma_number_event","#gamma number per event", 100000, 0.0, 100000.0);
    	gamma_number_event -> GetXaxis()-> SetTitle("event index");
	gamma_number_event -> GetYaxis()-> SetTitle("number of #gamma");

	TH1F *total_number_event = new TH1F("total_number_event","Total number of particles per event", 100000, 0.0, 100000.0);
    	total_number_event -> GetXaxis()-> SetTitle("event index");
	total_number_event -> GetYaxis()-> SetTitle("number of particles");

	vect* v= new vect;//declares a vector in which we are going to store all the particles from the same event
	double prev_index = 0.0;//saving the event index of the previous entry of the tree

	//Loop through the entries of the tree
	Long64_t nentries = T->GetEntries();
std::cout<<"Total number of entries: "<<nentries<<std::endl;
   	for (Long64_t i=0;i<nentries;i++){

	      	T->GetEntry(i);

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

			analyze_event(v, eta_invmass, muongamma_invmass_mother, eta_invmass_mother, muon_invmass_mother, muon_pt_mother, muon_ptcut_mother, gamma_pt_mother, mu_number_event, gamma_number_event, total_number_event);

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

			analyze_event(v, eta_invmass, muongamma_invmass_mother, eta_invmass_mother, muon_invmass_mother, muon_pt_mother, muon_ptcut_mother, gamma_pt_mother, mu_number_event, gamma_number_event, total_number_event);

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
	c1->Print("project5_etainvmass.pdf","pdf");

	muongamma_invmass_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project5_muongamma_invmass_mother.pdf","pdf");

	muon_invmass_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project5_muoninvmass_mother.pdf","pdf");

	eta_invmass_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project5_etainvmass_mother.pdf","pdf");

	muon_pt_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project5_mu_pt_mother.pdf","pdf");

	muon_ptcut_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project5_mu_ptcut_mother.pdf","pdf");


	gamma_pt_mother->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project5_gamma_pt_mother.pdf","pdf");

	mu_number_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project5_mu_number_event.pdf","pdf");

	gamma_number_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project5_gamma_number_event.pdf","pdf");

	total_number_event->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project5_total_number_event.pdf","pdf");

	muon_pt_mother->Draw();
	muon_ptcut_mother->Draw("same");

	TLegend *legend2 = new TLegend(0.7,0.5,0.9,0.7);
	TLegendEntry *leg3 = legend2->AddEntry("muon_pt_mother","P_t distribution before cuts","f");
  	leg3->SetFillColor(kBlue);	
	TLegendEntry *leg4 = legend2->AddEntry("muon_ptcut_mother","P_t distribution after cuts","f");
  	leg4->SetFillColor(kGreen);
	legend2->Draw("same");

	gPad->SetLogy();

	c1->Modified();
	c1->Update();
	c1->Print("project5_mu_pt_mother12.pdf","pdf");

	eta_invmass->Draw();
	eta_invmass_mother->Draw("same");

	TLegend *legend1 = new TLegend(0.7,0.5,0.9,0.7);
	TLegendEntry *leg1 = legend1->AddEntry("eta_invmass","#eta invariant mass (sig+bkg)","f");
  	leg1->SetFillColor(kBlue);	
	TLegendEntry *leg2 = legend2->AddEntry("eta_invmass_mother","#eta invariant mass (sig)","f");
  	leg2->SetFillColor(kGreen);
	legend1->Draw("same");

	gPad->SetLogx();

	c1->Modified();
	c1->Update();
	c1->Print("project5_etainvmass12.pdf","pdf");

	file_out->Write();

	delete file_in;
	delete file_out;

	// Done.
	return 0;
}


