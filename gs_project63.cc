//12/02/2018
//Program inspired by gs_project34.cc 
//This code: -sets up LHC environment (pp collision, 13 TeV CM energy) and soft QCD interactions;
//		-saves some features of all the etas and its decay products into two separate trees;
//		-saves muons, photons and misID pions from other processes in the decay products tree;
//		-restricts decays for eta: eta-> gamma A' where A'->mu anti-mu for the second instance of pythia
//		-tries to reconstruct invariant mass of eta from decay products by pairing up the decay particles

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"
#include "Pythia8/SigmaHiggs.h"
#include "cmath"
#include "Pythia8Plugins/HepMC2.h"
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"
#include "TVectorD.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TRandom2.h"

using namespace Pythia8;

int main() {

	HepMC::Pythia8ToHepMC ToHepMC;

    	// Specify file where HepMC events will be stored.
    	HepMC::IO_GenEvent ascii_io("/disk/moose/general/user72/gs_project63.dat", std::ios::out);

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

	//Set the features of the second pythia instance to force the eta to decay via A' and a gamma where the A' decays into a di-muon pair. 
	//The features of A' (this case is Z'0) are determined: its mass and width
	pythia1.particleData.m0(32, 0.27);
	pythia1.particleData.mMin(32, 0);

	double width = 1.191e-14;//calculating the width of the dark photon
	
	pythia1.particleData.mWidth(32, width);

	pythia1.readString("32:isResonance = false");

	pythia1.readString("221:oneChannel = 1 1 0 32 22");
	pythia1.readString("32:oneChannel = 1 1 0 13 -13");

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
	TFile *file = TFile::Open("/disk/moose/general/user72/gs_project63.root","recreate");

	//Create event
	Event *event = &pythia0.event;
	Event *event1 = &pythia1.event;


	//How many events shall we generate?
    	Long64_t nEvents=10;

	double inv_mass;

	std::vector<int> v;

	// Begin event loop. Generate event; skip if generation aborted.
	for (Long64_t iEvent = 0; iEvent < nEvents; ++iEvent) {

		//Generate one event. Skip if error.
		if (!pythia0.next()) continue;

		//Print out entire event contents and decays 
		pythia0.info.list(); pythia0.event.list();

		//Loop over all particles that have been generated in this event
		for (Long64_t i = 0; i < pythia0.event.size(); ++i) {	

			//check if particle is eta
			if (pythia0.event[i].id() == 221){

				v.push_back(i);

			} 

		}

		int n = v.size();

		for (Long64_t i = 0; i < n; ++i) {

				pythia1.event.reset(); 
				pythia1.event.append(pythia0.event[v[i]].id(), 22, 0, 0, 0, 0, 0, 0, pythia0.event[v[i]].p(), pythia0.event[v[i]].m()); 

				//Generate one event. Skip if error.
				if (!pythia1.next()) continue;

				pythia1.info.list(); pythia1.event.list();


		}

		v.clear();
	
	}// End event loop.


	// Statistics on event generation.
	pythia0.stat();

std::cout<<"A' width: "<<width<<std::endl;
std::cout<<"Z'0 is resonance? "<<pythia1.particleData.isResonance(32)<<std::endl;

	file->Write();
	delete file;

	// Done.
	return 0;
}
