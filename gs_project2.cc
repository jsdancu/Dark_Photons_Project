//09/10/2017
//Program inspired by gs_project1.cc and main92.cc
//This code: -sets up LHC environment (pp collision, 13 TeV CM energy) and soft QCD interactions;
//		-saves some features of the final particles of each event in n-tuples;
//		-restricts decays for eta: eta->mu anti-mu
//		-tries to reconstruct invariant mass of eta from decay products by pairing up the decay particles

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "Pythia8/SigmaHiggs.h"
#include "cmath"
#include "Pythia8Plugins/HepMC2.h"

using namespace Pythia8;

double invmass(Vec4 p1, Vec4 p2){

	double inv_mass=(p1 + p2).mCalc();
	return inv_mass;

}

int main() {

	HepMC::Pythia8ToHepMC ToHepMC;

    	// Specify file where HepMC events will be stored.
    	HepMC::IO_GenEvent ascii_io("gs_project2.dat", std::ios::out);

	// Generator. Process selection. LHC initialization. Histogram.
    	Pythia pythia;
    	pythia.readString("Random:setSeed = on");
    	pythia.readString("Random:seed = 0");
    	//pythia.readString("HardQCD:all = on");
    	pythia.readString("SoftQCD:all = on");
	pythia.readString("PhaseSpace:pTHatMin = 20.");

	//Set eta to decay to 2 muons and gamma
	//pythia.readString("221:onMode = 1 1 11 22 13 -13");

	//Set eta to decay to 2 muons
	pythia.readString("221:onMode = 0");
	pythia.readString("221:onIfMatch = 13 -13");

	//Initilise for p (2212) and  p (2212) collisins at 13 TeV
    	//initialization of LHC environment
    	pythia.readString("Beams:idA = 2212");
    	pythia.readString("Beams:idB = 2212");

    	//Start at 13 TeV
    	double myEcm =13000.;
    	pythia.settings.parm("Beams:eCM", myEcm);

	pythia.init();

	// Set up the ROOT TFile and TTree.
	TFile *file = TFile::Open("gs_project2.root","recreate");

	//Create event
	Event *event = &pythia.event;

	//Define a histograms into which we accumulate the invariant mass distribution for muon combinations
 
    	TH1F *eta_invmass = new TH1F("eta_invmass","Reconstructed eta invariant mass from muon pairs", 300, 0.0, 30.0);
    	eta_invmass -> GetXaxis()-> SetTitle("m (GeV)");

	TH1F *mu_number_event = new TH1F("mu_number_event","Muon-antimuon number per event", 100000, 0.0, 10000.0);
    	mu_number_event -> GetXaxis()-> SetTitle("event index");
	mu_number_event -> GetYaxis()-> SetTitle("number of #mu^{#pm}");

	//Create TTree
	TTree *T = new TTree("T","ev1 Tree");

	double index_var, id_var, energy_var, mass_var, px_var, py_var, pz_var, mother1_var, mother2_var, motherid1_var, motherid2_var;

	TBranch *index = T->Branch("index", &index_var);
	TBranch *id = T->Branch("id", &id_var);
	TBranch *energy = T->Branch("energy", &energy_var);
	TBranch *mass = T->Branch("mass", &mass_var);
	//TBranch *p = T->Branch("p", &p_var);
	TBranch *px = T->Branch("px", &px_var);
	TBranch *py = T->Branch("py", &py_var);
	TBranch *pz = T->Branch("pz", &pz_var);
	TBranch *mother1 = T->Branch("mother1", &mother1_var);
	TBranch *mother2 = T->Branch("mother2", &mother2_var);
	TBranch *motherid1 = T->Branch("motherid1", &motherid1_var);
	TBranch *motherid2 = T->Branch("motherid2", &motherid2_var);

	//How many events shall we generate?
    	Long64_t nEvents=100000;

	double inv_mass;
	const int nPrint=5;

	Long64_t total_final_state = 0;
	Long64_t total_mu_antimu = 0;

	bool antimu;

	// Begin event loop. Generate event; skip if generation aborted.
	for (Long64_t iEvent = 0; iEvent < nEvents; ++iEvent) {

		//Generate one event. Skip if error.
		if (!pythia.next()) continue;

		//Print out entire event contents and decays for first five events.
		if (iEvent < nPrint) {pythia.info.list(); pythia.event.list();}

		Long64_t final_state = 0;
		Long64_t mu_number = 0;
		Long64_t antimu_number = 0;

		antimu = false; //sets anti-muons status to false, as in they haven't been saved in the tree

		//Loop over all particles that have been generated in this event
		for (Long64_t i = 0; i < pythia.event.size(); ++i) {	

			//but only consider those that are "final state", i.e. ignore
			//intermediate ones that have decayed.
			if (pythia.event[i].isFinal()){

				++final_state;

				//Looking for muon pairs: first find a muon
				if (pythia.event[i].id()==13){

					++mu_number;

					//adding muon to the tree
					index_var = iEvent;
					id_var = pythia.event[i].id();
					energy_var = pythia.event[i].e();
					mass_var = pythia.event[i].m();
					px_var = pythia.event[i].px();
					py_var = pythia.event[i].py();
					pz_var = pythia.event[i].pz();
					mother1_var = pythia.event[i].mother1();
					mother2_var = pythia.event[i].mother2();
					motherid1_var = pythia.event[pythia.event[i].mother1()].id();
					motherid2_var = pythia.event[pythia.event[i].mother2()].id();

					T->Fill();
std::cout<<"iEvent = "<<iEvent<<"    "<<"id = "<<id_var<<"    "<< "energy = "<<energy_var<<std::endl;

				   	for (Long64_t j = 0; j < pythia.event.size(); ++j){

						//If particle is anti-muon
						if (pythia.event[j].id()==-13){

							//check if the anti-muon has not been filled in the tree (prevent multiple saving of the same anti-muon)
							if(antimu == false){

								++antimu_number;

								//adding anti-muon to the tree
								index_var = iEvent;
								id_var = pythia.event[j].id();
								energy_var = pythia.event[j].e();
								mass_var = pythia.event[j].m();
								px_var = pythia.event[j].px();
								py_var = pythia.event[j].py();
								pz_var = pythia.event[j].pz();
								mother1_var = pythia.event[j].mother1();
								mother2_var = pythia.event[j].mother2();
								motherid1_var = pythia.event[pythia.event[j].mother1()].id();
								motherid2_var = pythia.event[pythia.event[j].mother2()].id();

								T->Fill();
std::cout<<"i = "<<i<<"    "<<"j = "<<j<<"    "<<"iEvent = "<<iEvent<<"    "<<"id = "<<id_var<<"    "<< "energy = "<<energy_var<<std::endl;
								
							}

							//Reconstructing invariant mass of all muon-anti-muon combinations
							inv_mass = invmass(pythia.event[i].p(), pythia.event[j].p());

							//Plot histogram of all muon-anti-muon combinations
							eta_invmass->Fill(inv_mass);	

						}

					}
					antimu = true;//all the anti-muons have been saved in the event

				}

			} 

		}

		std::cout<<"final state = "<<final_state<<setw(5)<<"muon number = "<<mu_number<<setw(5)<<"anti muon = "<<antimu_number<<setw(5)<<"difference = "<<final_state - mu_number - antimu_number<<std::endl;

		total_final_state = total_final_state + final_state;
		total_mu_antimu = total_mu_antimu + mu_number + antimu_number;
		mu_number_event->Fill(total_mu_antimu);
	
	}// End event loop.

	std::cout<<"total final state = "<<total_final_state<<setw(5)<<"muon+antimuon number = "<<total_mu_antimu<<std::endl;

	// Statistics on event generation.
	pythia.stat();

	//  Print and Write tree.
	T->Print();
	T->Write();

	file->Write();
	delete file;

	// Done.
	return 0;
}
