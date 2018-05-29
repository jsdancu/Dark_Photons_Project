//25/02/2018
//Program inspired by gs_project34.cc  
//This code: -sets up LHC environment (pp collision, 13 TeV CM energy) and soft QCD interactions;
//		-restricts decays for pi: pi->mu nu_mu
//		-measures angular distribution between pi and mu

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

double ang(Vec4 r1, Vec4 r2){

	double scalar = r1.px() * r2.px() + r1.py() * r2.py() + r1.pz() * r2.pz();
std::cout<<"scalar: "<<scalar<<std::endl;
	double abs_r1 = sqrt(pow( r1.px(), 2) + pow( r1.py(), 2) + pow( r1.pz(), 2));
std::cout<<"abs_r1: "<<abs_r1<<std::endl;
	double abs_r2 = sqrt(pow( r2.px(), 2) + pow( r2.py(), 2) + pow( r2.pz(), 2));
std::cout<<"abs_r2: "<<abs_r2<<std::endl;

	double angle = acos(scalar / (abs_r1 * abs_r2));

	return angle;

}

int main() {

	//HepMC::Pythia8ToHepMC ToHepMC;

    	// Specify file where HepMC events will be stored.
    	//HepMC::IO_GenEvent ascii_io("/disk/moose/general/user72/gs_project91.dat", std::ios::out);

  Pythia pythia0("", false), pythia1("", false);

	pythia0.readString("Random:setSeed = on");
    	pythia0.readString("Random:seed = 0");
    	pythia0.readString("SoftQCD:all = on");
	pythia0.readString("PhaseSpace:pTHatMin = 20.");    	
	pythia0.readString("Print:quiet = on");    	

	pythia1.readString("Random:setSeed = on");
    	pythia1.readString("Random:seed = 0");
	pythia1.readString("ProcessLevel:all = off");

	pythia0.readString("211:onMode = off");
	pythia0.readString("-211:onMode = off");

	pythia1.readString("211:mayDecay = true");
	pythia1.readString("211:onMode = off");
	pythia1.readString("211:onIfMatch = -13 14");
	
	pythia1.readString("-211:mayDecay = true");
	pythia1.readString("-211:onMode = off");
	pythia1.readString("-211:onIfMatch = 13 -14");

	//Initilise for p (2212) and  p (2212) collisins at 13 TeV
    	//initialization of LHC environment
    	pythia0.readString("Beams:idA = 2212");
    	pythia0.readString("Beams:idB = 2212");

    	//Start at 13 TeV
    	double myEcm =13000.;
    	pythia0.settings.parm("Beams:eCM", myEcm);

	pythia0.init();
	pythia1.init();

	// Set up the ROOT TFile and TTree.
	TFile *file = TFile::Open("/disk/moose/general/user72/gs_project91.root","recreate");

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

	TH1D *angular_distribution = new TH1D("angular_distribution","Angular distribution between mother #pi^{#pm} and #mu^{#mp}", 50, -0.2, 1.0);
    	angular_distribution -> GetXaxis()-> SetTitle("angle #theta (rad)");
	angular_distribution -> GetYaxis()-> SetTitle("counts");

	//How many events shall we generate?
    	Long64_t nEvents=1000;

	const int nPrint=5;

	std::vector<int> v;
	//friend double angle;
	double angle;
	Vec4 r1, r2;

	// Begin event loop. Generate event; skip if generation aborted.
	for (Long64_t iEvent = 0; iEvent < nEvents; ++iEvent) {

		//Generate one event. Skip if error.
		if (!pythia0.next()) continue;

		//Print out entire event contents and decays for first five events.
		//if (iEvent < nPrint) {pythia0.info.list(); pythia0.event.list();}

		//Loop over all particles that have been generated in this event
		for (Long64_t i = 0; i < pythia0.event.size(); ++i) {	

		  if((pythia0.event[i].id() == 310) || (pythia0.event[i].id() == 310)){

				//adding eta to T1
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

		}

		int n = v.size();

		for (Long64_t i = 0; i < n; ++i) {

				pythia1.event.reset(); 
				pythia1.event.append(pythia0.event[v[i]].id(), 23, 0, 0, 0, 0, 0, 0, pythia0.event[v[i]].p(), pythia0.event[v[i]].m());

				//Generate one event. Skip if error.
				if (!pythia1.next()) continue;
				std::cout << pythia1.event[2].vProd() << "\n";

				//if (i < 2) {pythia1.info.list(); pythia1.event.list();}

				for (Long64_t j = 0; j < pythia1.event.size(); ++j) {

					if(pythia1.event[j].isFinal()){

						//adding eta products to T2
						index_var = iEvent;
						id_var = pythia1.event[j].id();
						energy_var = pythia1.event[j].e();
						mass_var = pythia1.event[j].m();
						px_var = pythia1.event[j].px();
						py_var = pythia1.event[j].py();
						pz_var = pythia1.event[j].pz();
						x_var = pythia1.event[j].xProd();
						y_var = pythia1.event[j].yProd();
						z_var = pythia1.event[j].zProd();
						mother1_var = pythia1.event[j].mother1();
						mother2_var = pythia1.event[j].mother2();
						motherid1_var = pythia1.event[pythia1.event[j].mother1()].id();
						motherid2_var = pythia1.event[pythia1.event[j].mother2()].id();

						T2->Fill();

						if((pythia1.event[j].id() == 13) || (pythia1.event[j].id() == -13)){
						//if((pythia1.event[j].id() == 14) || (pythia1.event[j].id() == -14)){

							r1 = pythia1.event[pythia1.event[j].mother1()].p();
							
std::cout<<"px1: "<<pythia1.event[pythia1.event[j].mother1()].px()<<std::endl;
std::cout<<"py1: "<<pythia1.event[pythia1.event[j].mother1()].py()<<std::endl;
std::cout<<"pz1: "<<pythia1.event[pythia1.event[j].mother1()].pz()<<std::endl;
							r2 = pythia1.event[j].p();
std::cout<<"px2: "<<pythia1.event[j].px()<<std::endl;
std::cout<<"py2: "<<pythia1.event[j].py()<<std::endl;
std::cout<<"pz2: "<<pythia1.event[j].pz()<<std::endl;

							angle = theta(r1, r2);
							angular_distribution->Fill(angle);
							
						
						}

					}

				}

		}

		v.clear();
	
	}// End event loop.


	// Statistics on event generation.
	pythia0.stat();

	TCanvas *c1=new TCanvas("c1","",600,600);

	gPad->SetLogy();

	angular_distribution->Draw();

	c1->Modified();
	c1->Update();
	c1->Print("project91_angular_distribution.pdf","pdf");

	file->Write();
	delete file;

	// Done.
	return 0;
}
