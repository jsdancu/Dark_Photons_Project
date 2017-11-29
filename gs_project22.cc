//16/11/2017
//Program inspired by gs_project23.cc 
//This code: -sets up LHC environment (pp collision, 13 TeV CM energy) and soft QCD interactions;
//		-saves some features of all the etas and its decay products into two separate trees;
//		-saves muons from other processes in the decay products tree;
//		-restricts decays for eta: eta->mu anti-mu for the second instance of pythia
//		-tries to reconstruct invariant mass of eta from decay products by pairing up the decay particles

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom2.h"
#include "Pythia8/SigmaHiggs.h"
#include "cmath"
#include "Pythia8Plugins/HepMC2.h"
#include <vector>

using namespace Pythia8;

double invmass(Vec4 p1, Vec4 p2){

	double inv_mass=(p1 + p2).mCalc();
	return inv_mass;

}

int main() {

	HepMC::Pythia8ToHepMC ToHepMC;

    	// Specify file where HepMC events will be stored.
    	HepMC::IO_GenEvent ascii_io("gs_project22.dat", std::ios::out);

	// Generator. Process selection. LHC initialization. Histogram.
    	Pythia pythia0, pythia1;

	pythia0.readString("Random:setSeed = on");
    	pythia0.readString("Random:seed = 0");
    	pythia0.readString("SoftQCD:all = on");
	pythia0.readString("PhaseSpace:pTHatMin = 20.");    	

	pythia1.readString("Random:setSeed = on");
    	pythia1.readString("Random:seed = 0");
    	pythia1.readString("SoftQCD:all = on");
	pythia1.readString("PhaseSpace:pTHatMin = 20.");
	pythia1.readString("ProcessLevel:all = off");

	//Set eta to decay 
	pythia0.readString("221:onMode = on");

	pythia1.readString("221:onMode = off");
	pythia1.readString("221:onIfMatch = 13 -13");

	//Initilise for p (2212) and  p (2212) collisins at 13 TeV
    	//initialization of LHC environment
    	pythia0.readString("Beams:idA = 2212");
    	pythia0.readString("Beams:idB = 2212");

	pythia1.readString("Beams:idA = 2212");
    	pythia1.readString("Beams:idB = 2212");

    	//Start at 13 TeV
    	double myEcm =13000.;
    	pythia0.settings.parm("Beams:eCM", myEcm);
	pythia1.settings.parm("Beams:eCM", myEcm);

	pythia0.init();
	pythia1.init();

	// Set up the ROOT TFile and TTree.
	TFile *file = TFile::Open("gs_project22.root","recreate");

	//Create event
	Event *event = &pythia0.event;
	Event *event1 = &pythia1.event;

	//Define a histograms into which we accumulate the invariant mass distribution for muon combinations
    	TH1D *eta_invmass = new TH1D("eta_invmass","Reconstructed eta invariant mass from muon pairs", 100, 0.0, 1.0);
    	eta_invmass -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *mu_number_event = new TH1D("mu_number_event","Muon-antimuon number per event without #eta-> #mu^{+} #mu^{-}", 100000, 0.0, 100000.0);
    	mu_number_event -> GetXaxis()-> SetTitle("event index");
	mu_number_event -> GetYaxis()-> SetTitle("number of #mu^{#pm}");

	//Create TTree for etas
	TTree *T1 = new TTree("T1","ev1 Tree");

	double index_var, id_var, energy_var, mass_var, px_var, py_var, pz_var, mother1_var, mother2_var, motherid1_var, motherid2_var;

	TBranch *index1 = T1->Branch("index", &index_var);
	TBranch *id1 = T1->Branch("id", &id_var);
	TBranch *energy1 = T1->Branch("energy", &energy_var);
	TBranch *mass1 = T1->Branch("mass", &mass_var);
	TBranch *px1 = T1->Branch("px", &px_var);
	TBranch *py1 = T1->Branch("py", &py_var);
	TBranch *pz1 = T1->Branch("pz", &pz_var);
	TBranch *mother11 = T1->Branch("mother1", &mother1_var);
	TBranch *mother21 = T1->Branch("mother2", &mother2_var);
	TBranch *motherid11 = T1->Branch("motherid1", &motherid1_var);
	TBranch *motherid21 = T1->Branch("motherid2", &motherid2_var);

	//Creat TTree for eta products
	TTree *T2 = new TTree("T2","ev1 Tree");

	TBranch *index2 = T2->Branch("index", &index_var);
	TBranch *id2 = T2->Branch("id", &id_var);
	TBranch *energy2 = T2->Branch("energy", &energy_var);
	TBranch *mass2 = T2->Branch("mass", &mass_var);
	TBranch *px2 = T2->Branch("px", &px_var);
	TBranch *py2 = T2->Branch("py", &py_var);
	TBranch *pz2 = T2->Branch("pz", &pz_var);
	TBranch *mother12 = T2->Branch("mother1", &mother1_var);
	TBranch *mother22 = T2->Branch("mother2", &mother2_var);
	TBranch *motherid12 = T2->Branch("motherid1", &motherid1_var);
	TBranch *motherid22 = T2->Branch("motherid2", &motherid2_var);

	//How many events shall we generate?
    	Long64_t nEvents=100000;

	double inv_mass;
	const int nPrint=5;

	std::vector<int> v;

	// Begin event loop. Generate event; skip if generation aborted.
	for (Long64_t iEvent = 0; iEvent < nEvents; ++iEvent) {

		//Generate one event. Skip if error.
		if (!pythia0.next()) continue;

		//Print out entire event contents and decays for first five events.
		//if (iEvent < nPrint) {pythia0.info.list(); pythia0.event.list();}

		Long64_t eta_number = 0;
		Long64_t mu_antimu_number = 0;

		//Loop over all particles that have been generated in this event
		for (Long64_t i = 0; i < pythia0.event.size(); ++i) {	

			if(pythia0.event[i].isFinal() && ((pythia0.event[i].id() == 13) || (pythia0.event[i].id() == -13))){

				mu_antimu_number++;

				index_var = iEvent;
				id_var = pythia0.event[i].id();
				energy_var = pythia0.event[i].e();
				mass_var = pythia0.event[i].m();
				px_var = pythia0.event[i].px();
				py_var = pythia0.event[i].py();
				pz_var = pythia0.event[i].pz();
				mother1_var = pythia0.event[i].mother1();
				mother2_var = pythia0.event[i].mother2();
				motherid1_var = pythia0.event[pythia0.event[i].mother1()].id();
				motherid2_var = pythia0.event[pythia0.event[i].mother2()].id();

				T2->Fill();

			}

			//check if particle is eta
			if (pythia0.event[i].id() == 221){

				++eta_number;

				//adding eta to T1
				index_var = iEvent;
				id_var = pythia0.event[i].id();
				energy_var = pythia0.event[i].e();
				mass_var = pythia0.event[i].m();
				px_var = pythia0.event[i].px();
				py_var = pythia0.event[i].py();
				pz_var = pythia0.event[i].pz();
				mother1_var = pythia0.event[i].mother1();
				mother2_var = pythia0.event[i].mother2();
				motherid1_var = pythia0.event[pythia0.event[i].mother1()].id();
				motherid2_var = pythia0.event[pythia0.event[i].mother2()].id();

				T1->Fill();

				v.push_back(i);

			} 

		}

		mu_number_event -> SetBinContent(iEvent, mu_antimu_number);//filling the mu-antimu/event histogram;

		TRandom2 *rndm=new TRandom2(0);
		int n = v.size();
		int x = (int) rndm->Uniform(n);
		//std::cout<< "size of vector: " << n <<std::endl;
		//std::cout<< "random number: " << x << std::endl;

		for (Long64_t i = 0; i < n; ++i) {

			if(i==x){

				pythia1.event.reset(); 
				pythia1.event.append(pythia0.event[v[x]].id(), 23, 0, 0, 0, 0, 0, 0, pythia0.event[v[x]].p(), pythia0.event[v[x]].m()); 

				//Generate one event. Skip if error.
				if (!pythia1.next()) continue;

				//if (i < 2) {pythia1.info.list(); pythia1.event.list();}

				Double_t E = 0.0;
				Double_t px = 0.0;
				Double_t py = 0.0;
				Double_t pz = 0.0;

				for (Long64_t j = 0; j < pythia1.event.size(); ++j) {

					if(pythia1.event[j].isFinal()){

						//adding eta products to T2
						index_var = iEvent;
						id_var = pythia1.event[j].id();
//std::cout<< "id_var" << id_var <<std::endl;
						energy_var = pythia1.event[j].e();
						mass_var = pythia1.event[j].m();
						px_var = pythia1.event[j].px();
						py_var = pythia1.event[j].py();
						pz_var = pythia1.event[j].pz();
						mother1_var = pythia1.event[j].mother1();
						mother2_var = pythia1.event[j].mother2();
						motherid1_var = pythia1.event[pythia1.event[j].mother1()].id();
						motherid2_var = pythia1.event[pythia1.event[j].mother2()].id();

						T2->Fill();

						E = E + energy_var;
						px = px + px_var;
						py = py + py_var;
						pz = pz + pz_var;

					}

				}

				Double_t inv_mass = sqrt(pow(E, 2)- pow(px, 2) - pow(py, 2) - pow(pz, 2));
				eta_invmass -> Fill(inv_mass);

			}
			else{

				double m = pythia0.event[v[i]].daughterList().size();
				for (Long64_t j = 0; j < m; ++j) {

					//adding eta products to T2
					index_var = iEvent;
					id_var = pythia0.event[pythia0.event[v[i]].daughterList()[j]].id();
					energy_var = pythia0.event[pythia0.event[v[i]].daughterList()[j]].e();
					mass_var = pythia0.event[pythia0.event[v[i]].daughterList()[j]].m();
					px_var = pythia0.event[pythia0.event[v[i]].daughterList()[j]].px();
					py_var = pythia0.event[pythia0.event[v[i]].daughterList()[j]].py();
					pz_var = pythia0.event[pythia0.event[v[i]].daughterList()[j]].pz();
					mother1_var = pythia0.event[pythia0.event[v[i]].daughterList()[j]].mother1();
					mother2_var = pythia0.event[pythia0.event[v[i]].daughterList()[j]].mother2();
					motherid1_var = pythia0.event[pythia0.event[pythia0.event[v[i]].daughterList()[j]].mother1()].id();//i.e. 221, but just for checking
					motherid2_var = pythia0.event[pythia0.event[pythia0.event[v[i]].daughterList()[j]].mother2()].id();//i.e. 0

					T2->Fill();

				}

			}

		}

		v.clear();
	
	}// End event loop.


	// Statistics on event generation.
	pythia0.stat();

	//  Print and Write trees.
	T1->Print();
	T1->Write();
	T2->Print();
	T2->Write();

	file->Write();
	delete file;

	// Done.
	return 0;
}
