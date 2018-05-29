//28/02/2018
//Program inspired by gs_project92.cc  & the one sent by Phil in Python
//This code: -sets up LHC environment (pp collision, 13 TeV CM energy) and soft QCD interactions;
//		-produces decays of mesons+baryons containing a b-quark into a particle containing a c-quark+muon and further decaying the particle with c-quark into whatever+muon: B->D mu ->X mu mu
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

void doDecays(Particle &bParticle, Long64_t iEvent, TTree *T2, TTree *T3, Pythia &pythia1){
/*
	Pythia pythia1("", false);

	pythia1.readString("Random:setSeed = on");
    	pythia1.readString("Random:seed = 0");
    	pythia1.readString("SoftQCD:all = on");
	pythia1.readString("PhaseSpace:pTHatMin = 20.");
	pythia1.readString("ProcessLevel:all = off");
	pythia1.readString("ParticleDecays:FSRinDecays = off");
	pythia1.readString("Print:quiet = on");

	//Initilise for p (2212) and  p (2212) collisins at 13 TeV
    	//initialization of LHC environment
	pythia1.readString("Beams:idA = 2212");
    	pythia1.readString("Beams:idB = 2212");

    	//Start at 13 TeV
    	double myEcm =13000.;
	pythia1.settings.parm("Beams:eCM", myEcm);

	pythia1.init();
*/
	//Create event
	Event *event1 = &pythia1.event;

	pythia1.event.reset(); 
	pythia1.event.append(bParticle.id(), 23, 0, 0, 0, 0, 0, 0, bParticle.p(), bParticle.m()); 
	pythia1.event[1].vProd(bParticle.vProd());
	pythia1.event[1].tau(bParticle.tau0() * pythia1.rndm.exp() );

	//Generate one event. Skip if error.
	//if (!pythia1.next()) continue;

	//if (i < 2) {pythia1.info.list(); pythia1.event.list();}
	//pythia1.info.list(); pythia1.event.list();

	bool mu1 = false;
	bool mu2 = false;
	bool cquark = false;
	bool candidate = false;
	bool the_one = false;
	std::vector<int> v1, v2, v3;

	double index_var, id_var, energy_var, mass_var, px_var, py_var, pz_var, x_var, y_var, z_var, mother1_var, mother2_var, motherid1_var, motherid2_var, granny1_var, granny2_var, grannyid1_var, grannyid2_var;
				
	for (Long64_t j = 0; j < pythia1.event.size(); ++j) {

		if((pythia1.event[j].id() == 13) || (pythia1.event[j].id() == -13)){ 				

			//std::cout<< "mu from B: "<<pythia1.event[j].id()  <<std::endl;

			//If particle comes from the decayed particle that contains a b-quark
			if(pythia1.event[pythia1.event[j].mother1()].id() == bParticle.id()){
				//std::cout<< "B->mu + whatever: "<<pythia1.event[j].id()  <<std::endl;
				mu1 = true;
				v1.push_back(j);

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

			}
			else{
				//std::cout<< "B->whatever->mu + whatever: "<<pythia1.event[j].id()  <<std::endl;
				mu2 = true;
				v2.push_back(j);
			}

		}

					
		if(abs(pythia1.particleData.heaviestQuark(pythia1.event[j].id())) == 4){ 				

			//std::cout<< "c quark particle: "<<pythia1.event[j].id()  <<std::endl;
			cquark = true;

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

			int m =pythia1.event[j].daughterList().size();
			for (int k = 0; k < m; ++k) {

				//v4.push_back(pythia1.event[pythia1.event[j].daughterList()[k]].id());

				index_var = iEvent;
				id_var = pythia1.event[k].id();
				energy_var = pythia1.event[k].e();
				mass_var = pythia1.event[k].m();
				px_var = pythia1.event[k].px();
				py_var = pythia1.event[k].py();
				pz_var = pythia1.event[k].pz();
				x_var = pythia1.event[k].xProd();
				y_var = pythia1.event[k].yProd();
				z_var = pythia1.event[k].zProd();
				mother1_var = pythia1.event[k].mother1();
				mother2_var = pythia1.event[k].mother2();
				motherid1_var = pythia1.event[pythia1.event[k].mother1()].id();
				motherid2_var = pythia1.event[pythia1.event[k].mother2()].id();
				granny1_var = pythia1.event[pythia1.event[k].mother1()].mother1();
				granny2_var = pythia1.event[pythia1.event[k].mother2()].mother2();
				grannyid1_var = pythia1.event[pythia1.event[pythia1.event[k].mother1()].mother1()].id();
				grannyid2_var = pythia1.event[pythia1.event[pythia1.event[k].mother2()].mother2()].id();

				T3->Fill();

				//std::cout<< "c quark particle daughter: "<<pythia1.event[pythia1.event[j].daughterList()[k]].id() <<std::endl;
				if((pythia1.event[pythia1.event[j].daughterList()[k]].id() == 13) || (pythia1.event[pythia1.event[j].daughterList()[k]].id() == -13)){
					v3.push_back(j);
					
				} 

			}

	
		}
					
	}

	int n1 = v1.size();
	int n2 = v2.size();
	int n3 = v3.size();

	for (Long64_t i = 0; i < n1; ++i) {
		for (Long64_t j = 0; j < n2; ++j) {
			for (Long64_t k = 0; k < n3; ++k) {
				std::cout<< "Candidate event iEvent+ i: " << iEvent<< "  " << i << "   particles: " <<pythia1.particleData.name(pythia1.event[v1[i]].id()) << "   " << pythia1.particleData.name(pythia1.event[v2[j]].id()) << "   " << pythia1.particleData.name(pythia1.event[v3[k]].id())<<std::endl;
			}
		}
	}				

	v1.clear();
	v2.clear();
	v3.clear();

}


int main() {

	//HepMC::Pythia8ToHepMC ToHepMC;

    	// Specify file where HepMC events will be stored.
    	//HepMC::IO_GenEvent ascii_io("/disk/moose/general/user72/gs_project93.dat", std::ios::out);

	Pythia pythia0("", false), pythia1("", false);

	pythia0.readString("Random:setSeed = on");
    	pythia0.readString("Random:seed = 0");
    	pythia0.readString("SoftQCD:all = on");
	pythia0.readString("PhaseSpace:pTHatMin = 20.");    	
	pythia0.readString("Print:quiet = on");

	pythia1.readString("Random:setSeed = on");
    	pythia1.readString("Random:seed = 0");
    	pythia1.readString("SoftQCD:all = on");
	pythia1.readString("PhaseSpace:pTHatMin = 20.");
	pythia1.readString("ProcessLevel:all = off");
	pythia1.readString("ParticleDecays:FSRinDecays = off");
	pythia1.readString("Print:quiet = on");

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
	TFile *file = TFile::Open("/disk/moose/general/user72/gs_project93.root","recreate");

	//Create event
	Event *event = &pythia0.event;
	Event *event1 = &pythia1.event;

	int pid = pythia0.particleData.nextId(1);
	ParticleDataEntry *pde;
	DecayChannel channel;
	std::string prods;
/*
	while(pid != 0){

		pid = pythia0.particleData.nextId(pid);

		 if((abs(pythia0.particleData.heaviestQuark(pid)) == 5) || (abs(pythia0.particleData.heaviestQuark(pid)) == 4)){ 
			
			pde = pythia0.particleData.particleDataEntryPtr(pid);
			std::cout<< pde->name()<<std::endl;
			std::cout<<"---------"<<std::endl;

			for(int i=0; i<pde->sizeChannels(); i++){
					
			    	channel = pde->channel(i);
				
			    	if (!channel.contains(13) || !channel.contains(-13)){

					channel.onMode(0);

				}
			    	else{

					prods = "";

					for(int j=0; j< channel.multiplicity(); j++){

					    prods += pythia0.particleData.name(channel.product(j)) + "   ";

					}

					std::cout<<channel.bRatio() << "   " << prods <<std::endl;

				//}

			}

			std::cout<<"---------"<<std::endl;

		}

	}
*/

	while(pid != 0){

		pid = pythia0.particleData.nextId(pid);

		 if(abs(pythia0.particleData.heaviestQuark(pid)) == 5){ 

			pde = pythia0.particleData.particleDataEntryPtr(pid);

			//Set B particle not to decay 
			pythia0.readString(std::to_string(pid)+":onMode = off");

		}
/*		else if(abs(pythia0.particleData.heaviestQuark(pid)) == 4){

			for(int i=0; i<pde->sizeChannels(); i++){
					
				    channel = pde->channel(i);
				
				    if (!channel.contains(13) || !channel.contains(-13)){

					channel.onMode(0);

				}

			}

		}
*/
	}




	//Create TTree for parent particles containing b-quarks
	TTree *T1 = new TTree("T1","ev1 Tree");

	double index_var, id_var, energy_var, mass_var, px_var, py_var, pz_var, x_var, y_var, z_var, mother1_var, mother2_var, motherid1_var, motherid2_var, granny1_var, granny2_var, grannyid1_var, grannyid2_var;

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

	//Create TTree for decay products (particles with c-quark and muon)
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

	//Create TTree for decay products (particles with light-quark and muon)
	TTree *T3 = new TTree("T3","ev1 Tree");

	TBranch *index3 = T3->Branch("index", &index_var);
	TBranch *id3 = T3->Branch("id", &id_var);
	TBranch *energy3 = T3->Branch("energy", &energy_var);
	TBranch *mass3 = T3->Branch("mass", &mass_var);
	TBranch *px3 = T3->Branch("px", &px_var);
	TBranch *py3 = T3->Branch("py", &py_var);
	TBranch *pz3 = T3->Branch("pz", &pz_var);
	TBranch *x3 = T3->Branch("x", &x_var);
	TBranch *y3 = T3->Branch("y", &y_var);
	TBranch *z3 = T3->Branch("z", &z_var);
	TBranch *mother13 = T3->Branch("mother1", &mother1_var);
	TBranch *mother23 = T3->Branch("mother2", &mother2_var);
	TBranch *motherid13 = T3->Branch("motherid1", &motherid1_var);
	TBranch *motherid23 = T3->Branch("motherid2", &motherid2_var);
	TBranch *granny13 = T3->Branch("granny1", &granny1_var);
	TBranch *granny23 = T3->Branch("granny", &granny2_var);
	TBranch *grannyid13 = T3->Branch("grannyid1", &grannyid1_var);
	TBranch *grannyid23 = T3->Branch("grannyid2", &grannyid2_var);

	//Create TTree for random photons
	TTree *T4 = new TTree("T4","ev1 Tree");

	TBranch *index4 = T4->Branch("index", &index_var);
	TBranch *id4 = T4->Branch("id", &id_var);
	TBranch *energy4 = T4->Branch("energy", &energy_var);
	TBranch *mass4 = T4->Branch("mass", &mass_var);
	TBranch *px4 = T4->Branch("px", &px_var);
	TBranch *py4 = T4->Branch("py", &py_var);
	TBranch *pz4 = T4->Branch("pz", &pz_var);
	TBranch *x4 = T4->Branch("x", &x_var);
	TBranch *y4 = T4->Branch("y", &y_var);
	TBranch *z4 = T4->Branch("z", &z_var);
	TBranch *mother14 = T4->Branch("mother1", &mother1_var);
	TBranch *mother24 = T4->Branch("mother2", &mother2_var);
	TBranch *motherid14 = T4->Branch("motherid1", &motherid1_var);
	TBranch *motherid24 = T4->Branch("motherid2", &motherid2_var);


	//How many events shall we generate?
    	Long64_t nEvents=100;
	int bDecays = 1000;
	int repeat;

	const int nPrint=5;

	std::vector<int> v;
	Vec4 r1, r2;
	Particle bParticle;

	// Begin event loop. Generate event; skip if generation aborted.
	for (Long64_t iEvent = 0; iEvent < nEvents; ++iEvent) {

		//Generate one event. Skip if error.
		if (!pythia0.next()) continue;

		//Print out entire event contents and decays for first five events.
		//if (iEvent < nPrint) {pythia0.info.list(); pythia0.event.list();}

		Long64_t b_number = 0;

		//Loop over all particles that have been generated in this event
		for (Long64_t i = 0; i < pythia0.event.size(); ++i) {	

			 if(abs(pythia0.particleData.heaviestQuark(pythia0.event[i].id())) == 5){ 

				b_number++;

//std::cout<<"pid: "<<pythia0.event[i].id()<<std::endl;


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

				//gamma_number++;

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

				T4->Fill();

			}		

		}	

		int n = v.size();

		TRandom2 *rndm=new TRandom2(0);
		int x = (int) rndm->Uniform(n);

		bool mu1 = false;
		bool mu2 = false;
		bool cquark = false;
		bool candidate = false;
		bool the_one = false;
		std::vector<int> v1, v2, v3, v4, v5;

		for (Long64_t i = 0; i < n; ++i) {

			bParticle = pythia0.event[v[i]];

			if(i==x){repeat = bDecays;}
			else{repeat = 1;}
/*
				for (Long64_t ii = 0; ii < bDecays; ++ii) {

					doDecays(bParticle, iEvent, T2, T3, pythia1);

				}	

			}
			else{

				doDecays(bParticle, iEvent, T2, T3, pythia1);

			}
*/
				for (Long64_t ii = 0; ii < bDecays; ++ii) {
				
					pythia1.event.reset(); 
					pythia1.event.append(pythia0.event[v[i]].id(), 23, 0, 0, 0, 0, 0, 0, pythia0.event[v[i]].p(), pythia0.event[v[i]].m()); 
					pythia1.event[1].vProd(pythia0.event[v[i]].vProd());
					pythia1.event[1].tau(pythia1.event[1].tau0() * pythia1.rndm.exp() );

					//Generate one event. Skip if error.
					if (!pythia1.next()) continue;

					//if (i < 2) {pythia1.info.list(); pythia1.event.list();}
					//pythia1.info.list(); pythia1.event.list();
				
					for (Long64_t j = 0; j < pythia1.event.size(); ++j) {

						if((pythia1.event[j].id() == 13) || (pythia1.event[j].id() == -13)){ 				

							//std::cout<< "mu from B: "<<pythia1.event[j].id()  <<std::endl;

							//If particle comes from the decayed particle that contains a b-quark
							if(pythia1.event[pythia1.event[j].mother1()].id() == pythia0.event[v[i]].id()){
								//std::cout<< "B->mu + whatever: "<<pythia1.event[j].id()  <<std::endl;
								mu1 = true;
								v1.push_back(j);

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

							}
							else{
								//std::cout<< "B->whatever->mu + whatever: "<<pythia1.event[j].id()  <<std::endl;
								mu2 = true;
								v2.push_back(j);
							}

						}

					
						if(abs(pythia0.particleData.heaviestQuark(pythia1.event[j].id())) == 4){ 				

							//std::cout<< "c quark particle: "<<pythia1.event[j].id()  <<std::endl;
							cquark = true;

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

							int m =pythia1.event[j].daughterList().size();
							for (int k = 0; k < m; ++k) {

								//v4.push_back(pythia1.event[pythia1.event[j].daughterList()[k]].id());

								index_var = iEvent;
								id_var = pythia1.event[k].id();
								energy_var = pythia1.event[k].e();
								mass_var = pythia1.event[k].m();
								px_var = pythia1.event[k].px();
								py_var = pythia1.event[k].py();
								pz_var = pythia1.event[k].pz();
								x_var = pythia1.event[k].xProd();
								y_var = pythia1.event[k].yProd();
								z_var = pythia1.event[k].zProd();
								mother1_var = pythia1.event[k].mother1();
								mother2_var = pythia1.event[k].mother2();
								motherid1_var = pythia1.event[pythia1.event[k].mother1()].id();
								motherid2_var = pythia1.event[pythia1.event[k].mother2()].id();
								granny1_var = pythia1.event[pythia1.event[k].mother1()].mother1();
								granny2_var = pythia1.event[pythia1.event[k].mother2()].mother2();
								grannyid1_var = pythia1.event[pythia1.event[pythia1.event[k].mother1()].mother1()].id();
								grannyid2_var = pythia1.event[pythia1.event[pythia1.event[k].mother2()].mother2()].id();

								T3->Fill();

								//std::cout<< "c quark particle daughter: "<<pythia1.event[pythia1.event[j].daughterList()[k]].id() <<std::endl;
								if((pythia1.event[pythia1.event[j].daughterList()[k]].id() == 13) || (pythia1.event[pythia1.event[j].daughterList()[k]].id() == -13)){
									v3.push_back(j);
					
								} 

							}

	/*
							pde = pythia1.particleData.particleDataEntryPtr(pythia1.event[j].id());
							m = v4.size();

							for(int r=0; r<2; r++){
								channel = pde->channel(r);
								prods = "";
								for (int k = 0; k < m; ++k) {
									std::cout<<iEvent <<" " << i << " " << " " << ii << " " << j << " " << r <<  "  " <<  k << std::endl;// " "  << v4[k] <<std::endl;
									std::cout<<v4[k] <<std::endl;
									/*if(channel.contains(v4[k])){
										the_one = true;
									}*/
	/*							}
							}

							v4.clear();
	/*
	//						for(int r=0; r<pde->sizeChannels(); r++){
					
								channel = pde->channel(r);
								prods = "";
				
								//if(m == channel.multiplicity()){
									the_one = false;
									for (int p = 0; p < m; ++p) {
										std::cout<<iEvent <<" " << i << " " << " " << ii << " " << j << " " << r <<  "  " <<  p << " "  << v4[p] <<std::endl;
										if(channel.contains(v4[p])){
											the_one = true;
										}// break;}
										//else{prods = prods + "  " + std::to_string(v4[p]);}
										//prods = prods + "  " + std::to_string(v4[p]);
										//if(!channel.contains(v4[p])){std::cout<< "D channel contains : "<<channel.contains(v4[p])<<std::endl;}
									}
								//}
								//std::cout<< "D products : "<<pythia1.event[j].id() << ":   "<< prods <<std::endl;

							}

							//channel = pde->channel.contains(pythia1.event[pythia1.event[j].daughterList()].id());
	/*
							if(the_one == true){
								std::cout<< "D products +  channel.bRatio(): "<<pythia1.event[j].id() << ":   "<< prods<< "   " << channel.bRatio() <<std::endl;
								v5.push_back(channel.bRatio());
							}

							the_one = false;
	*/
						
						}
					
					}

					int n1 = v1.size();
					int n2 = v2.size();
					int n3 = v3.size();

					for (Long64_t iii = 0; iii < n1; ++iii) {
						for (Long64_t j = 0; j < n2; ++j) {
							for (Long64_t k = 0; k < n3; ++k) {
								std::cout<< "Candidate event iEvent+ i: " << iEvent<< "  " << i << "   particles: " <<pythia0.particleData.name(pythia1.event[v1[iii]].id()) << "   " << pythia0.particleData.name(pythia1.event[v2[j]].id()) << "   " << pythia0.particleData.name(pythia1.event[v3[k]].id())<<std::endl;
							}
						}
					}				

					v1.clear();
					v2.clear();
					v3.clear();
					v5.clear();

				}

			//}
	
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
	T3->Print();
	T3->Write();
	T4->Print();
	T4->Write();

	file->Write();

	delete file;

	// Done.
	return 0;
}
