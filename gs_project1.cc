 //gSystem.Load("../libPythia8");
 // File: main01.cc
 // This is a simple test program.
 // It studies a few particle multiplicities within some nominal p_T, energy, angular
 //  acceptances in pp collisions at 13 TeV.
 // Adapted from one of many Pythia example programs, see
 //  /home/nkw/gs2015/pythia8201_sl6/share/Pythia8/examples:

#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "Pythia8/SigmaHiggs.h"
#include "cmath"
#include "Pythia8Plugins/HepMC2.h"

using namespace Pythia8;
int main() {

  bool write_data=true;
  //Initialisations
  // Interface for conversion from Pythia8::Event to HepMC event.
  //if (write_data){
    HepMC::Pythia8ToHepMC ToHepMC;

    // Specify file where HepMC events will be stored.
    HepMC::IO_GenEvent ascii_io("gs_project1.dat", std::ios::out);
  
    // Generator. Process selection. LHC initialization. Histogram.
    Pythia pythia;
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 0");
    pythia.readString("HardQCD:all = on");
    pythia.readString("SoftQCD:all = on");
    // pythia.readString("SLHA:file = ../sps1a.spc");
    //pythia.readString("SLHA:readFrom = 2");
    // pythia.readString("SUSY:all = on");
    // pythia.readString("HiggsSM:all = on");
    //pythia.readString("WeakSingleBoson:all = on");
    pythia.readString("Top:all = on");

    //pythia.readString("Top:ffbar2ttbar(s:gmZ) = on");

    pythia.readString("HiggsSM:ffbar2HZ = on");
    //Restrict Z decay modes to e and mu
    pythia.readString("23:onMode = off");
    pythia.readString("23:onIfAny = 11 13");
    // Set m_H
    pythia.readString("25:m0=126.22");
    pythia.settings.flag("SpaceShower:QEDshowerByL", false);
    pythia.settings.flag("SpaceShower:QEDshowerByQ", false);
    pythia.settings.flag("SpaceShower:QCDshower", false);
    pythia.settings.flag("ISR:all",false);
    pythia.settings.flag("FSR:all",false);



    //Initilise for p (2212) and  p (2212) collisins at 13 TeV
    // LHC initialization at Higgs mass.
    pythia.readString("Beams:idA =  2212");
    pythia.readString("Beams:idB = 2212");

    //Start at 13 TeV
    double myEcm =13000.;
    pythia.settings.parm("Beams:eCM", myEcm);
    pythia.init();

    // Check that Z0 decay channels set correctly.
    // pythia.particleData.listChanged();

    //Open root file to store the data
    TFile *file = TFile::Open("gs_project1.root","recreate");

    //Create event
    Event *event = &pythia.event;

    //Create a TTree (n-tuple) with branch called event
    TTree *T = new TTree("T","ev1 Tree");
    T->Branch("event",&event);

    //Define a few histograms into which we accumulate
    //numbers of particles, momenta, etc., event-by-event.

    TH1F *mult = new TH1F("mult","charged multiplicity", 200, -0.5, 199.5);
    // ptcharged-> GetXaxis()->SetTitle("Multiplicity");
    TH1F *mult_acc = new TH1F("mult_acc","charged multiplicity inside acceptance", 200, -0.5, 199.5);
    //ptchargedbbb-> GetXaxis()->SetTitle("Multiplicity");
    TH1F *multn = new TH1F("multn","neutral multiplicity", 200, -0.5, 199.5);
    TH1F *multn_acc = new TH1F("multn_acc","neutral multiplicity inside acceptance", 200, -0.5, 199.5);

    //ptcharged-> GetXaxis()->SetTitle("Multiplicity");
    TH1F *cthetacharged = new TH1F("cthetacharged","cos(theta) of charged particles", 100, -1., 1.);
    TH1F *cthetacharged_acc = new TH1F("cthetacharged_acc","cos(theta) of charged particles inside acceptance", 100, -1., 1.);

    TH1F *cthetaneutral = new TH1F("cthetaneutral","cos(theta) of neutral particles", 100, -1., 1.);
    TH1F *cthetaneutral_acc = new TH1F("cthetaneutral_acc","cos(theta) of neutral particles inside acceptance", 100, -1., 1.);

    TH1F *ptcharged = new TH1F("ptcharged","pt of charged particles", 100, 0, 100);
    TH1F *ptcharged_acc = new TH1F("ptcharged_acc","pt of charged particles inside acceptance", 100, 0, 10);

    TH1F *eneutral = new TH1F("eneutral","e of neutral particles", 100, 0, 500);
    TH1F *eneutral_acc = new TH1F("eneutral_acc","e of neutral particles inside acceptance", 100, 0, 100);

    TH1F *ptH = new TH1F("ptH","pt of Higgs", 1000, 0, 1000);

    TH1F *diphoton_invmass = new TH1F("diphoton_invmass","Inv. mass of two photons", 100, 120, 130);
    TH1F *diphoton_invmass_b = new TH1F("diphoton_invmass_b","Inv. mass of two photons", 100, 0, 130);
    TH1F *dineutral_invmass = new TH1F("dineutral_invmass","Inv. mass of two neutral particles", 100, 0, 200);

    double costheta=0;
    double momn=0;
    //Define number of events for which we may want to print a lot of information

    const int nPrint=5;

  // Say for example we can detect particles down to cos(theta)<0.9
    double costheta_acceptance=0.9;

  // Reconstruction of charged particles only for pt>0.2 GeV, given the magnetic field
  // neutrals for energy>0.1GeV (is this to low/high?),
    double pt_acceptance=0.2;
    double e_acceptance=0.1;

  //How many events shall we generate?
    int nEvents=100;

  //You will want to read this web page to understand how to get at the event record.
  //http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html

    // Begin event loop.
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {

   	//Generate one event. Skip if error.
      	if (!pythia.next()) continue;

	//Print out entire event contents and decays for first five events.
	if (iEvent < nPrint) {pythia.info.list(); pythia.event.list();}

	// Fill the pythia event into the TTree.  
	// Warning: the files will rapidly become large if all events 
	// are saved. In some cases it may be convenient to do some 
	// processing of events and only save those that appear 
	// interesting for future analyses.
	T->Fill();

	//Find out the number of all final state particles that are neutral/charged and similarly
	//in different regions of the detector "acceptance" (defined by polar angles
	//(relative to beam axis)).
	int nCharged = 0;
	int nCharged_acc = 0;
	int nNeutral = 0;
	int nNeutral_acc = 0;

	//Loop over all particles that have been generated in this event
	for (int i = 0; i < pythia.event.size(); ++i) {

		//a few counters
		double neutrals_invmass=0;

		//but only consider those that are "final state", i.e. ignore
		//intermediate ones that have decayed.
		if (pythia.event[i].isFinal()){ 

			//Momentum of current particle in the loop
			momn=sqrt(pow(pythia.event[i].px(),2)+pow(pythia.event[i].py(),2)+pow(pythia.event[i].pz(),2));

			//cos(theta) (polar angle, relative to incoming beam axis)
			costheta=cos(pythia.event[i].theta());

			//Start by looking only at charged particles
			if (pythia.event[i].isCharged()){

				//increment number of charged, final state particles in this event
				++nCharged;
			  	//Fill polar angle and p_T distributons.
				cthetacharged->Fill(costheta);
				ptcharged->Fill(pythia.event[i].pT());

			    	if (abs(costheta) <costheta_acceptance && pythia.event[i].pT()>pt_acceptance){
		   		      // And additionally, fill separate histos for those particles that satisfy min. p_T and polar angle requirements.
				      ++nCharged_acc;
				      cthetacharged_acc->Fill(costheta);
				      ptcharged_acc->Fill(pythia.event[i].pT());
			    	}
			  }

			  //Next, look at neutral particles (also final state only)
			  else {

				//increment number of charged, final state particles in this event
				++nNeutral;

			  	//Fill polar angle and energy distributons.
				cthetaneutral->Fill(costheta);
				eneutral->Fill(pythia.event[i].e());

			    	if (abs(costheta) <costheta_acceptance && pythia.event[i].e()>e_acceptance){
		   		// And additionally, fill separate histos for those particles that satisfy min. energy and polar angle requirements.
			        ++nNeutral_acc;
			        cthetaneutral_acc->Fill(costheta);
			        eneutral_acc->Fill(pythia.event[i].e());
			   	}

			   	//We may want to look at all pairwise combinations of particles - here we do neutral with neutral only.
			   	//Avoid including a pairwise combination of with the same particle!
			   	for (int j = i+1; j < pythia.event.size(); ++j){

			   		if (pythia.event[j].isFinal() && pythia.event[j].isNeutral()){ 

						// Loop over all pairwise combinations of neutral particles
						neutrals_invmass=(pythia.event[i].p() + pythia.event[j].p()).mCalc();
						// std::cout<<"neutrals_invmass, 1 = "<<neutrals_invmass<<std::endl;
						//Plot histo of all neutral-neutral combinations
						dineutral_invmass->Fill(neutrals_invmass);
						//and the same, but only if both are photons
						if (pythia.event[i].id()==22 && pythia.event[j].id()==22) {
						  diphoton_invmass->Fill(neutrals_invmass);
						  diphoton_invmass_b->Fill(neutrals_invmass);
						}
			      		}
			    	}
			    
			  }


		  	  //Back to considering all final state particles now


		  	  //Keep our eyes open for a Higgs!
			  if (pythia.event[i].id()==25){
			    ptH->Fill(pythia.event[i].pT());
			  }

		} // End of all final state particles only condition

	} // End of particle loop. Statistics. Histogram. Done.
		  


	//Fill histo. with total number of charged and neutral particles in the event and within our nominal acceptance.
	mult->Fill(nCharged);
	mult_acc->Fill(nCharged_acc);
	multn->Fill(nNeutral);
	multn_acc->Fill(nNeutral_acc);
	//Only write out data for use in Delphes if we asked for it.
	//Try to write out as HepMC event file
	// Construct new empty HepMC event and fill it.
	// Units will be as chosen for HepMC build; but can be changed
	// by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
	if (write_data){
		HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
		ToHepMC.fill_next_event( pythia, hepmcevt );
			    
		// Write the HepMC event to file. Done with it.
		ascii_io << hepmcevt;
		delete hepmcevt;

	}


    } //End event loop

    // Statistics on event generation.
    pythia.stat();

    //  Write tree.
    T->Print();
    T->Write();

    cout << mult;
    file->Write();

    delete file;

    return 0;

  }
