#include "Pythia8/Pythia.h"
using namespace Pythia8;


//==========================================================================

int main() {

  // Create minimum bias generator.
  Pythia minbias;
  minbias.readString("SoftQcd:all = on");
  minbias.init();
  Event& mbevent = minbias.event;
  
  // Create the decayer.
  Pythia decayer;
  decayer.readString("ProcessLevel:all = off");
  decayer.readString("221:onMode = off");
  decayer.readString("221:onIfMatch = 22 13 -13");
  decayer.init();
  Event& dcevent = decayer.event;


  // Begin of event loop.
  for (int iEvent = 0; iEvent < 10; ++iEvent) {

    // Find an eta.
    if (!minbias.next()) continue;
    Particle *prt(0);
    for (int iPrt = 0; iPrt < mbevent.size(); ++iPrt) {
      prt = &mbevent[iPrt];
      if (prt->idAbs() == 221) break;
    }
    if (!prt || prt->idAbs() != 221) continue;
    
    // Do the decay.
    dcevent.reset();
    dcevent.append(prt->id(), 23, 0, 0, 0, 0, 0, 0, prt->p(), prt->m());
    
    if (!decayer.next()) continue;
    for (int iDec = 0; iDec < dcevent.size(); ++iDec) {
      Particle *dec = &dcevent[iDec];
      if (dec->isFinal()) cout << dec->id() << "\n";
    }
    
  // End of event loop.
  }

  // Done.
  return 0;
}
