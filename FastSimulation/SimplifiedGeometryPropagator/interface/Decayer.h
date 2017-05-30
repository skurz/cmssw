#ifndef FASTSIM_DECAYER_H
#define FASTSIM_DECAYER_H

#include <memory>
#include <vector>


///////////////////////////////////////////////
// Decayer
//
// Description: Implementation of non-stable particle decays
//
// Author: L. Vanelderen
// Date: 13 May 2014
//
// Revision: Class structure modified to match SimplifiedGeometryPropagator
//           S. Kurz, 29 May 2017
//////////////////////////////////////////////////////////


namespace gen {
  class P8RndmEngine;
}

namespace CLHEP {
  class HepRandomEngine;
}

namespace Pythia8 {
  class Pythia;
}

namespace fastsim
{
    class Particle;
    class Decayer 
    {
    public:
    	
    	Decayer();
    	~Decayer();
    	void decay(const Particle & particle,std::vector<std::unique_ptr<Particle> > & secondaries,CLHEP::HepRandomEngine & engine) const;
    	
        private:
    	
    	std::unique_ptr<Pythia8::Pythia> pythia_; 
    	std::unique_ptr<gen::P8RndmEngine> pythiaRandomEngine_;
    };
}
#endif
