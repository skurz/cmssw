#ifndef FASTSIM_STRAIGHTTRAJECTORY_H
#define FASTSIM_STRAIGHTTRAJECTORY_H

#include "FastSimulation/SimplifiedGeometryPropagator/interface/Trajectory.h"


///////////////////////////////////////////////
// StraightTrajectory
//
// Description: Definition the trajectory of an uncharged particle
//
// Author: L. Vanelderen, S. Kurz
// Date: 29 May 2017
//////////////////////////////////////////////////////////


namespace fastsim
{
    class StraightTrajectory : public Trajectory
    {
    public:
	StraightTrajectory(const Particle & particle) : Trajectory(particle) {;}
	bool crosses(const BarrelSimplifiedGeometry & layer) const override {return true;}
	double nextCrossingTimeC(const BarrelSimplifiedGeometry & layer) const override;
	void move(double deltaTimeC) override;
    };
}

#endif
