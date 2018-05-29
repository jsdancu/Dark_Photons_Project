//01/03/2018
//Program inspired by gs_project62.cc 
//This code: -sets up LHC environment (pp collision, 13 TeV CM energy) and soft QCD interactions;
//		-saves some features of all the etas and its decay products into two separate trees;
//		-saves muons, photons and misID pions from other processes in the decay products tree;
//		-restricts decays for eta: eta-> A' where A'->mu anti-mu for the second instance of pythia
//		-tries to reconstruct invariant mass of eta from decay products by pairing up the decay particles
//		-produces 10^5 events

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"
#include "Pythia8/SigmaHiggs.h"
#include "cmath"
//#include "Pythia8Plugins/HepMC2.h"
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

double invmass(Vec4 p1, Vec4 p2){

	double inv_mass=(p1 + p2).mCalc();
	return inv_mass;

}

double misEnergy(double p2){

	double m0= 0.10566;
	double energy=sqrt(p2 + pow(m0, 2));
	return energy;

}

double lepton_width(double m_A, double m_l, double epsilon2){

	double alpha_EM = 1.0/137.0;
	/*HepMC::Pythia8ToHepMC ToHepMC;
	Pythia pythia;
	CoupSM& coup = pythia.CoupSM;

	pythia.init(); 
	
	double alpha_EM = coup.alphaEM(pow(m_A, 2.0));*/

	double Gamma_l = epsilon2 * alpha_EM /3.0 * m_A * (1.0 + 2.0 * pow(m_l, 2.0) / pow(m_A, 2.0)) * sqrt(1.0 - 4.0 * pow(m_l, 2.0) / pow(m_A, 2.0));

	return Gamma_l;

}

double lepton_width1(double m_A, double m_l){

	/*HepMC::Pythia8ToHepMC ToHepMC;
	Pythia pythia;
	CoupSM& coup = pythia.CoupSM;

	pythia.init(); 
	double alpha_EM = coup.alphaEM(pow(m_A, 2.0));*/

	double alpha_EM = 1.0/137.0;

	double Gamma_l = alpha_EM /3.0 * m_A * (1.0 + 2.0 * pow(m_l, 2.0) / pow(m_A, 2.0)) * sqrt(1.0 - 4.0 * pow(m_l, 2.0) / pow(m_A, 2.0));

	return Gamma_l;

}

double hadrons_width(double R_mu, double Gamma_l){

	double Gamma_h = Gamma_l * R_mu;

	return Gamma_h;

}

double Gamma_A(double epsilon2, double massA){

	double x, y, u1, u2;
	std::vector<double> m_A, R_mu, v1, v2;

	for(int i=0;i<14;i++){

		m_A.push_back(0.22 + i * 0.01);
		R_mu.push_back(0.0);
		v1.push_back(0.0);
		v2.push_back(0.0);

	}


	std::ifstream file("rdata.txt");
    	if(file.is_open()){

		file >> x >> y >> u1 >> u2;
		m_A.push_back(x);
		R_mu.push_back(y);
		v1.push_back(u1);
		v2.push_back(u2);

		// Begin reading data
		while(file >> x >> y >> u1 >> u2){
	
			m_A.push_back(x);
			R_mu.push_back(y);
			v1.push_back(u1);
			v2.push_back(u2);

		}

    	}

	double mA[12], Rmu[12], w1[12], w2[12], wx1[12], wx2[12];

	for(int i=0;i<12;i++){

		mA[i] = m_A[i+14];
		Rmu[i] = R_mu[i+14];
		w1[i] = v1[i+14];
		w2[i] = abs(v2[i+14]);
		wx1[i] = 0.0;
		wx2[i] = 0.0;

	}	

	TF1 *fit = new TF1("fit", "[0]*exp([1]*x)", 0.36, 0.46);

	TGraphAsymmErrors *gre = new TGraphAsymmErrors(12, mA, Rmu, wx2, wx1, w2, w1);
   	gre->Fit("fit","EX0");

	double p0 = fit->GetParameter(0);
	double p1 = fit->GetParameter(1);

	for(int i=6;i<14;i++){

		R_mu[i] = p0 * exp(p1 * m_A[i]);

	}

	/*Pythia pythia;
	ParticleData& pdt = pythia.particleData;

	pythia.init(); 

	double m_mu = pdt.m0(13);//muon mass
	double m_e = pdt.m0(11);//electron mass
	double m_tau = pdt.m0(15);//tau mass
	*/
	

	double m_mu = 0.10566;//muon mass
	double m_e = 5.110e-04;//electron mass
	double m_tau = 1.77682;//tau mass	

	double Gamma_inv = 0.0;

	int N = m_A.size();
	double Gamma[N], Gamma_mu[N], Gamma_e[N], Gamma_h[N];

	double G;

	for(int i=0;i<N;i++){

		Gamma_mu[i] = lepton_width(m_A[i], m_mu, epsilon2);
		Gamma_e[i] = lepton_width(m_A[i], m_e, epsilon2);
		Gamma_h[i] = hadrons_width(R_mu[i], Gamma_mu[i]);

		if(m_A[i] < 2.0*0.28){

			Gamma[i] = Gamma_e[i] + Gamma_mu[i] + Gamma_inv;

		}
		else if(m_A[i] < 2.0*m_tau){

			Gamma[i] = Gamma_e[i] + Gamma_mu[i] + Gamma_h[i] + Gamma_inv;;

		}
		else if(m_A[i] > 2.0*m_tau){

			Gamma[i] = Gamma_e[i] + Gamma_mu[i] + lepton_width(m_A[i], m_tau, epsilon2) + Gamma_h[i] + Gamma_inv;

		}

		if(massA == m_A[i]){G = Gamma[i];}

	}

	return G;

}


int main() {

	//Open R-data file
	TFile  *file_in = new TFile("rdata.txt", "READ");

	//HepMC::Pythia8ToHepMC ToHepMC;

    	// Specify file where HepMC events will be stored.
    	//HepMC::IO_GenEvent ascii_io("/disk/moose/general/user72/friday0203/gs_project28.dat", std::ios::out);

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
	pythia1.readString("ParticleDecays:FSRinDecays = off");

	//Set eta to decay 
	pythia0.readString("221:onMode = on");

	double epsilon2 = 1e-7;//kinetic mixing coupling
	double massA = 0.27;//dark photon mass
	double width = Gamma_A(epsilon2, massA);//calculating the width of the dark photon

	//Set the features of the second pythia instance to force the eta to decay via A' and a gamma where the A' decays into a di-muon pair. 
	//The features of A' (this case is Z'0) are determined: its mass and width
	pythia1.particleData.m0(32, massA);
	pythia1.particleData.mMin(32, 0);
	pythia1.particleData.mWidth(32, width);

	pythia1.readString("32:isResonance = false");

	pythia1.readString("221:oneChannel = 1 1 0 32");
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
	TFile *file = TFile::Open("/disk/moose/general/user72/friday0203/gs_project28.root","recreate");

	//Create event
	Event *event = &pythia0.event;
	Event *event1 = &pythia1.event;
/*
	//Define a histograms into which we accumulate the invariant mass distribution for muon combinations
    	TH1D *eta_invmass = new TH1D("eta_invmass","Reconstructed eta invariant mass from #gamma + muon pairs", 100, 0.0, 1.0);
    	eta_invmass -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *mu_invmass = new TH1D("mu_invmass","Reconstructed di-muon invariant mass", 100, 0.0, 2.0);
    	mu_invmass -> GetXaxis()-> SetTitle("m (GeV)");

	TH1D *mu_number_event = new TH1D("mu_number_event","Muon-antimuon number per event without #eta-> #gamma #mu^{+} #mu^{-}", 100000, 0.0, 100000.0);
    	mu_number_event -> GetXaxis()-> SetTitle("event index");
	mu_number_event -> GetYaxis()-> SetTitle("number of #mu^{#pm}");

	TH1D *gamma_number_event = new TH1D("gamma_number_event","#gamma number per event without #eta-> #gamma #mu^{+} #mu^{-}", 100000, 0.0, 100000.0);
    	gamma_number_event -> GetXaxis()-> SetTitle("event index");
	gamma_number_event -> GetYaxis()-> SetTitle("number of #gamma");
*/
	//Create TTree for etas
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

	//Creat TTree for eta products
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
/*
	//Creat TTree for misID pions
	TTree *T3 = new TTree("T3","ev1 Tree");

	TBranch *index3 = T3->Branch("index", &index_var);
	TBranch *id3 = T3->Branch("id", &id_var);
	TBranch *energy3 = T3->Branch("energy", &energy_var);
	TBranch *mass3 = T3->Branch("mass", &mass_var);
	TBranch *px3 = T3->Branch("px", &px_var);
	TBranch *py3 = T3->Branch("py", &py_var);
	TBranch *pz3 = T3->Branch("pz", &pz_var);
	TBranch *mother13 = T3->Branch("mother1", &mother1_var);
	TBranch *mother23 = T3->Branch("mother2", &mother2_var);
	TBranch *motherid13 = T3->Branch("motherid1", &motherid1_var);
	TBranch *motherid23 = T3->Branch("motherid2", &motherid2_var);
*/
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
		if (iEvent < nPrint) {pythia0.info.list(); pythia0.event.list();}
/*
		Long64_t eta_number = 0;
		Long64_t mu_antimu_number = 0;
		Long64_t gamma_number = 0;
		Long64_t pi_number = 0;
*/
		//Loop over all particles that have been generated in this event
		for (Long64_t i = 0; i < pythia0.event.size(); ++i) {	
/*
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

				T2->Fill();

			}

			//saving the pions from the event onto a separate tree with energy as if they were identified as muons
			if (pythia0.event[i].isFinal() && ((pythia0.event[i].id() == 211) || (pythia0.event[i].id() == -211))){

				pi_number++;

				//adding pi to T3
				index_var = iEvent;
				id_var = pythia0.event[i].id();

				//change energy of pion as if it was misidentified as a muon
				energy_var = misEnergy(pythia0.event[i].pAbs2());

				mass_var = pythia0.event[i].m();
				px_var = pythia0.event[i].px();
				py_var = pythia0.event[i].py();
				pz_var = pythia0.event[i].pz();
				mother1_var = pythia0.event[i].mother1();
				mother2_var = pythia0.event[i].mother2();
				motherid1_var = pythia0.event[pythia0.event[i].mother1()].id();
				motherid2_var = pythia0.event[pythia0.event[i].mother2()].id();

				T3->Fill();

			}
*/
			//check if particle is eta
			if (pythia0.event[i].id() == 221){

				//++eta_number;

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

		TRandom2 *rndm=new TRandom2(0);
		int n = v.size();
		int x = (int) rndm->Uniform(n);
		//std::cout<< "size of vector: " << n <<std::endl;
		//std::cout<< "random number: " << x << std::endl;

		for (Long64_t i = 0; i < n; ++i) {

			if(i==x){

				pythia1.event.reset(); 
				pythia1.event.append(pythia0.event[v[x]].id(), 22, 0, 0, 0, 0, 0, 0, pythia0.event[v[x]].p(), pythia0.event[v[x]].m()); 

				//Generate one event. Skip if error.
				if (!pythia1.next()) continue;

				//if (i < 2) {pythia1.info.list(); pythia1.event.list();}
/*
				Double_t E = 0.0;
				Double_t px = 0.0;
				Double_t py = 0.0;
				Double_t pz = 0.0;

				Vec4 p1, p2;
*/
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

						//Saving the muon production vertex
						x_var = pythia1.event[j].xProd();
						y_var = pythia1.event[j].yProd();
						z_var = pythia1.event[j].zProd();

						if((pythia1.event[j].id() == 13) || (pythia1.event[j].id() == -13)){

							//mu_invmass -> Fill(pythia1.event[j].m());
							//++mu_antimu_number;

							//Saving the eta as the mother particle instead of A'
							mother1_var = pythia1.event[pythia1.event[j].mother1()].mother1();
							mother2_var = pythia1.event[pythia1.event[j].mother2()].mother2();
							motherid1_var = pythia1.event[pythia1.event[pythia1.event[j].mother1()].mother1()].id();
							motherid2_var = pythia1.event[pythia1.event[pythia1.event[j].mother2()].mother2()].id();


						}
						else if(pythia1.event[j].id() == 22){
							mother1_var = pythia1.event[j].mother1();
							mother2_var = pythia1.event[j].mother2();
							motherid1_var = pythia1.event[pythia1.event[j].mother1()].id();
							motherid2_var = pythia1.event[pythia1.event[j].mother2()].id();
						}

						T2->Fill();

						/*E = E + energy_var;
						px = px + px_var;
						py = py + py_var;
						pz = pz + pz_var;

						if(pythia1.event[j].id() == 13){p1 = pythia1.event[j].p();}
						else if(pythia1.event[j].id() == -13){p2 = pythia1.event[j].p();}

						if((pythia1.event[j].id() == 13) || (pythia1.event[j].id() == -13)){
							++mu_antimu_number;
						}
						//else if(pythia1.event[j].id() == 22){T3->Fill();}
						*/
					}

				}

				//Double_t inv_mass = sqrt(pow(E, 2)- pow(px, 2) - pow(py, 2) - pow(pz, 2));
				//eta_invmass -> Fill(inv_mass);

				//Double_t inv_mass2 = invmass(p1, p2);
				//mu_invmass -> Fill(inv_mass2);

			}
			/*else{

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

					if(pythia0.event[pythia0.event[v[i]].daughterList()[j]].id() == 22){
						++gamma_number;
						T3->Fill();
					}
					else if((pythia0.event[pythia0.event[v[i]].daughterList()[j]].id() == 211) || (pythia0.event[pythia0.event[v[i]].daughterList()[j]].id() == -211)){
						T3->Fill();
					}

				}

			}*/

		}

		//mu_number_event -> SetBinContent(iEvent+1, mu_antimu_number);//filling the mu-antimu/event histogram;
		//gamma_number_event -> SetBinContent(iEvent+1, gamma_number);//filling the gamma/event histogram;

		v.clear();
	
	}// End event loop.


	// Statistics on event generation.
	pythia0.stat();

	//  Print and Write trees.
	T1->Print();
	T1->Write();
	T2->Print();
	T2->Write();

std::cout<<"A' width: "<<width<<std::endl;

	file->Write();
	delete file;
	delete file_in;

	// Done.
	return 0;
}
