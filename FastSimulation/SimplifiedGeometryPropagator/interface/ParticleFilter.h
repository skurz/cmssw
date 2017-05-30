#ifndef FASTSIM_PARTICLEFILTER
#define FASTSIM_PARTICLEFILTER

#include "DataFormats/Math/interface/LorentzVector.h"
#include <vector>


///////////////////////////////////////////////
// ParticleFilter
//
// Description: (Kinematic) cuts on the particles that are propagated 
//
// Author: Patrick Janot
// Date: 09 Dez 2003
//
// Revision: Class structure modified to match SimplifiedGeometryPropagator
//           S. Kurz, 29 May 2017
//////////////////////////////////////////////////////////


namespace edm
{
    class ParameterSet;
}

namespace fastsim
{
    class Particle;
    class ParticleFilter
    {
    public:
	ParticleFilter(const edm::ParameterSet & cfg);
	bool accepts(const Particle & particle) const;
	bool acceptsEn(const Particle & particle) const;
	bool acceptsVtx(const math::XYZTLorentzVector & originVertexPosition) const;

    private:
	// see constructor for comments
	double chargedPtMin2_, EMin_, protonEMin_;
	double cos2ThetaMax_;
	double vertexRMax2_,vertexZMax_;
	std::vector<int> skipParticles_;
    };
}

#endif
