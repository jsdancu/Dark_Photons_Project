#include "Pythia8/Pythia.h"
#include "Pythia8/SigmaHiggs.h"

using namespace Pythia8;

/*int main() {

	

	// Generator. Process selection. LHC initialization. Histogram.
    	Pythia pythia0;

	pythia0.readString("Random:setSeed = on");
    	pythia0.readString("Random:seed = 0");
    	pythia0.readString("SoftQCD:all = on");
	pythia0.readString("PhaseSpace:pTHatMin = 20.");    	


	std::cout << "K_S width: " << pythia0.particleData.tau0(23) << std::endl;

}*/

//26/02/2018
//Program inspired by gs_project91.cc  & gs_project62.cc
//This code: -sets up LHC environment (pp collision, 13 TeV CM energy) and soft QCD interactions;
//		-restricts decays for K_S: K_S->pi pi
//		-generates samples for further analysis

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TRandom2.h"
#include "Pythia8/SigmaHiggs.h"
#include "cmath"
//#include "Pythia8Plugins/HepMC2.h"
#include <vector>

using namespace Pythia8;

double invmass(Vec4 p1, Vec4 p2){

	double inv_mass=(p1 + p2).mCalc();
	return inv_mass;

}

double misEnergy(double p2){

	double m0= 0.10566;
	double energy=sqrt(p2 + pow(m0, 2));
	return energy;

}


int main() {

	//HepMC::Pythia8ToHepMC ToHepMC;

    	// Specify file where HepMC events will be stored.
    	//HepMC::IO_GenEvent ascii_io("/disk/moose/general/user72/gs_project92.dat", std::ios::out);

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
		pythia1.readString("ParticleDecays:FSRinDecays = off");

	pythia0.readString("310:onMode = on");

	pythia1.readString("310:onMode = off");
	pythia1.readString("310:onIfMatch = 211 -211");

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
	TFile *file = TFile::Open("/disk/moose/general/user72/gs_project92.root","recreate");

	//Create event
	Event *event = &pythia0.event;
	Event *event1 = &pythia1.event;


	//Create TTree for K_S
	TTree *T1 = new TTree("T1","ev1 Tree");

	double index_var, id_var, energy_var, mass_var, px_var, py_var, pz_var, x_var, y_var, z_var, mother1_var, mother2_var, motherid1_var, motherid2_var;

	TBranch *index1 = T1->Branch("index", &index_var);
	TBranch *id1 = T1->Branch("id", &id_var);
	TBranch *energy1 = T1->Branch("energy", &energy_var);
	TBranch *mass1 = T1->Branch("mass", &mass_var);
	TBranch *px1 = T1->Branch("px", &px_var);
	TBranch *py1 = T1->Branch("py", &py_var);
	TBranch *pz1 = T1->Branch("pz", &pz_var);
	TBranch *x1 = T1->Branch("x", &x_var);
	TBranch *y1 = T1->Branch("y", &y_var);
	TBranch *z1 = T1->Branch("z", &z_var);
	TBranch *mother11 = T1->Branch("mother1", &mother1_var);
	TBranch *mother21 = T1->Branch("mother2", &mother2_var);
	TBranch *motherid11 = T1->Branch("motherid1", &motherid1_var);
	TBranch *motherid21 = T1->Branch("motherid2", &motherid2_var);

	//Creat TTree for K_S products (pi pi)
	TTree *T2 = new TTree("T2","ev1 Tree");

	TBranch *index2 = T2->Branch("index", &index_var);
	TBranch *id2 = T2->Branch("id", &id_var);
	TBranch *energy2 = T2->Branch("energy", &energy_var);
	TBranch *mass2 = T2->Branch("mass", &mass_var);
	TBranch *px2 = T2->Branch("px", &px_var);
	TBranch *py2 = T2->Branch("py", &py_var);
	TBranch *pz2 = T2->Branch("pz", &pz_var);
	TBranch *x2 = T2->Branch("x", &x_var);
	TBranch *y2 = T2->Branch("y", &y_var);
	TBranch *z2 = T2->Branch("z", &z_var);
	TBranch *mother12 = T2->Branch("mother1", &mother1_var);
	TBranch *mother22 = T2->Branch("mother2", &mother2_var);
	TBranch *motherid12 = T2->Branch("motherid1", &motherid1_var);
	TBranch *motherid22 = T2->Branch("motherid2", &motherid2_var);


	//How many events shall we generate?
    	Long64_t nEvents=10;//0000;

	const int nPrint=5;

	std::vector<int> v;
	double angle;
	Vec4 r1, r2;

	// Begin event loop. Generate event; skip if generation aborted.
	for (Long64_t iEvent = 0; iEvent < nEvents; ++iEvent) {

		//Generate one event. Skip if error.
		if (!pythia0.next()) continue;

		//Print out entire event contents and decays for first five events.
		//if (iEvent < nPrint) {pythia0.info.list(); pythia0.event.list();}

		Long64_t KS_number = 0;
		Long64_t gamma_number = 0;
		Long64_t pi_number = 0;

		//Loop over all particles that have been generated in this event
		for (Long64_t i = 0; i < pythia0.event.size(); ++i) {	

			if(pythia0.event[i].id() == 310){

				++KS_number;

				//adding K_S to T1
				index_var = iEvent;
				id_var = pythia0.event[i].id();
				energy_var = pythia0.event[i].e();
				mass_var = pythia0.event[i].m();
				px_var = pythia0.event[i].px();
				py_var = pythia0.event[i].py();
				pz_var = pythia0.event[i].pz();
				x_var = pythia0.event[i].xProd();
				y_var = pythia0.event[i].yProd();
				z_var = pythia0.event[i].zProd();
				mother1_var = pythia0.event[i].mother1();
				mother2_var = pythia0.event[i].mother2();
				motherid1_var = pythia0.event[pythia0.event[i].mother1()].id();
				motherid2_var = pythia0.event[pythia0.event[i].mother2()].id();

				T1->Fill();

				v.push_back(i);

			}
			if (pythia0.event[i].id() == 22){

				gamma_number++;

				index_var = iEvent;
				id_var = pythia0.event[i].id();
				energy_var = pythia0.event[i].e();
				mass_var = pythia0.event[i].m();
				px_var = pythia0.event[i].px();
				py_var = pythia0.event[i].py();
				pz_var = pythia0.event[i].pz();
				x_var = pythia0.event[i].xProd();
				y_var = pythia0.event[i].yProd();
				z_var = pythia0.event[i].zProd();
				mother1_var = pythia0.event[i].mother1();
				mother2_var = pythia0.event[i].mother2();
				motherid1_var = pythia0.event[pythia0.event[i].mother1()].id();
				motherid2_var = pythia0.event[pythia0.event[i].mother2()].id();

				T2->Fill();

			}

		}

		int n = v.size();

		for (Long64_t i = 0; i < n; ++i) {

			//if(pythia0.event[v[i]].id() == 211){

				pythia1.event.reset(); 
				pythia1.event.append(pythia0.event[v[i]].id(), 23, 0, 0, 0, 0, 0, 0, pythia0.event[v[i]].p(), pythia0.event[v[i]].m()); 
				pythia1.event[1].vProd(pythia0.event[v[i]].vProd());
				pythia1.event[1].tau(pythia1.event[1].tau0() * pythia1.rndm.exp() );
				//Generate one event. Skip if error.
				if (!pythia1.next()) continue;

				//if (i < 2) {pythia1.info.list(); pythia1.event.list();}

				for (Long64_t j = 0; j < pythia1.event.size(); ++j) {

					if(pythia1.event[j].isFinal()){

						//adding K_S products to T2
						index_var = iEvent;
						id_var = pythia1.event[j].id();
//std::cout<< "id_var" << id_var <<std::endl;
						//Changing the energy of pions as misID muons
						energy_var = misEnergy(pythia1.event[j].pAbs2());

						mass_var = pythia1.event[j].m();
						px_var = pythia1.event[j].px();
						py_var = pythia1.event[j].py();
						pz_var = pythia1.event[j].pz();

						//Saving the pion production vertex
						x_var = pythia1.event[j].xProd();
						y_var = pythia1.event[j].yProd();
						z_var = pythia1.event[j].zProd();
std::cout<< "i "<<i<<"  id_var" << id_var << "x_var  " << x_var <<  "   y_var" << y_var <<std::endl;
						mother1_var = pythia1.event[j].mother1();
						mother2_var = pythia1.event[j].mother2();
						motherid1_var = pythia1.event[pythia1.event[j].mother1()].id();
						motherid2_var = pythia1.event[pythia1.event[j].mother2()].id();

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
