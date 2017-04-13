#ifndef FASTSIM_PARTCILE_H
#define FASTSIM_PARTICLE_H

#include "DataFormats/Math/interface/LorentzVector.h"

namespace fastsim
{
    class Particle
    {
    public:
	Particle(int pdgId,
		 const math::XYZTLorentzVector & position,
		 const math::XYZTLorentzVector & momentum)
	    : pdgId_(pdgId)
	    , charge_(-999.)
	    , position_(position)
	    , momentum_(momentum)
	    , remainingProperLifeTimeC_(-999.) // lifetime in ct
	    , simTrackIndex_(-1)
	    , simVertexIndex_(-1)
	    , genParticleIndex_(-1)
	{;}
	
	// setters
	void setSimTrackIndex(int index) {simTrackIndex_ = index;}
	void setSimVertexIndex(int index) {simVertexIndex_ = index;}
	void setGenParticleIndex(int index){genParticleIndex_ = index;}
	void setStable(){remainingProperLifeTimeC_ = -1.;}
	void setRemainingProperLifeTimeC(double remainingProperLifeTime){remainingProperLifeTimeC_ = remainingProperLifeTime;}
	void setCharge(double charge){charge_ = charge;}


	// ordinary getters
	int pdgId() const {return pdgId_;}
    double charge() const {return charge_;}
	const math::XYZTLorentzVector & position() const {return position_;}
	const math::XYZTLorentzVector & momentum() const {return momentum_;}
	double remainingProperLifeTimeC() const {return remainingProperLifeTimeC_;}
	int simTrackIndex() const {return simTrackIndex_;}
	int simVertexIndex() const {return simVertexIndex_;}
	int genParticleIndex() const {return genParticleIndex_;}
	bool isStable() const {return remainingProperLifeTimeC_ == -1.;}

	// other
    bool chargeIsSet() const {return charge_!=-999.;}
	bool remainingProperLifeTimeIsSet() const {return remainingProperLifeTimeC_ != -999.;}
	double gamma() const { return momentum().Gamma(); };

	// non-const getters
	math::XYZTLorentzVector & position() {return position_;}
	math::XYZTLorentzVector & momentum() {return momentum_;}

	friend std::ostream& operator << (std::ostream& os , const Particle & particle);

    private:
	const int pdgId_;
	double charge_;
	math::XYZTLorentzVector position_;
	math::XYZTLorentzVector momentum_;
	double remainingProperLifeTimeC_;
	int simTrackIndex_;
	int simVertexIndex_;
	int genParticleIndex_;
    };

    std::ostream& operator << (std::ostream& os , const Particle & particle);

}

#endif
