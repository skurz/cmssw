
#ifndef FASTSIM_HELIXTRAJECTORY_H
#define FASTSIM_HELIXTRAJECTORY_H

#include "FastSimulation/SimplifiedGeometryPropagator/interface/Trajectory.h"


///////////////////////////////////////////////
// HelixTrajectory
//
// Description: Definition the trajectory of a charged particle
//
// Author: L. Vanelderen, S. Kurz
// Date: 29 May 2017
//////////////////////////////////////////////////////////


namespace fastsim
{
    class HelixTrajectory : public Trajectory
    {
    public:
	HelixTrajectory(const Particle & particle,double magneticFieldZ);
	bool crosses(const BarrelSimplifiedGeometry & layer) const override;
	double nextCrossingTimeC(const BarrelSimplifiedGeometry & layer) const override;
	void move(double deltaTimeC) override;
    private:
	const double radius_;
	const double phi_;
	const double centerX_;
	const double centerY_;
	const double centerR_;
	const double minR_;
	const double maxR_;
	const double phiSpeed_;
    };
}

#endif
