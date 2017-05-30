#ifndef FASTSIM_PARTCILE_H
#define FASTSIM_PARTICLE_H

#include "DataFormats/Math/interface/LorentzVector.h"


///////////////////////////////////////////////
// Particle
//
// Description: Definition of a generic FastSim Particle which can be propagated through the detector (formerly ParticlePropagator)
//
// Author: L. Vanelderen, S. Kurz
// Date: 29 May 2017
//////////////////////////////////////////////////////////


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
	    , energyDeposit_(0)
	    , isLooper_(false)
	    , motherDeltaR_(-1)
		, motherPpdId_(0)
		, motherSimTrackIndex_(-999)
	{;}
	
	// setters
	void setSimTrackIndex(int index) {simTrackIndex_ = index;}
	void setSimVertexIndex(int index) {simVertexIndex_ = index;}
	void setGenParticleIndex(int index){genParticleIndex_ = index;}
	void setStable(){remainingProperLifeTimeC_ = -1.;}
	void setRemainingProperLifeTimeC(double remainingProperLifeTime){remainingProperLifeTimeC_ = remainingProperLifeTime;}
	void setCharge(double charge){charge_ = charge;}
	void setEnergyDeposit(double energyDeposit){energyDeposit_ = energyDeposit;}
	void setLooper(){isLooper_ = true;}
	void setMotherDeltaR(math::XYZTLorentzVector motherMomentum){
		motherDeltaR_ = (momentum_.Vect().Unit().Cross(motherMomentum.Vect().Unit())).R();
	}
	void setMotherPdgId(int id){motherPpdId_ = id;}
	void setMotherSimTrackIndex(int id){motherSimTrackIndex_ = id;}
	void resetMother(){motherDeltaR_ = -1; motherPpdId_ = 0; motherSimTrackIndex_ = -999;}


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
	double getEnergyDeposit() const {return energyDeposit_;}
	double isLooper() const {return isLooper_;}
	double getMotherDeltaR() const {return motherDeltaR_;}
	int getMotherPdgId() const {return motherPpdId_;}
	int getMotherSimTrackIndex() const {return motherSimTrackIndex_;}

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
	double energyDeposit_;
	bool isLooper_;
	double motherDeltaR_;
	int motherPpdId_;
	int motherSimTrackIndex_;
    };

    std::ostream& operator << (std::ostream& os , const Particle & particle);

}

#endif
