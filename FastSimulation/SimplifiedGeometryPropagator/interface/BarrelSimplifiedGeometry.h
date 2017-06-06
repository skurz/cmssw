#ifndef FASTSIM_BARRELSIMPLIFIEDGEOMETRY_H
#define FASTSIM_BARRELSIMPLIFIEDGEOMETRY_H

#include "FastSimulation/SimplifiedGeometryPropagator/interface/SimplifiedGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "TH1F.h"


///////////////////////////////////////////////
// BarrelSimplifiedGeometry
//
// Description: Implementation of a barrel detector layer
//
// Author: L. Vanelderen, S. Kurz
// Date: 29 May 2017
//////////////////////////////////////////////////////////


namespace fastsim{

    class BarrelSimplifiedGeometry : public SimplifiedGeometry
    {
    public:
	~BarrelSimplifiedGeometry(){};
	
	BarrelSimplifiedGeometry(double radius) :
	    SimplifiedGeometry(radius) {}
	
	BarrelSimplifiedGeometry(BarrelSimplifiedGeometry &&) = default;
	
	const double getRadius() const { return position_; }
	
	const double getThickness(const math::XYZTLorentzVector & position) const override
	{
	    if(!this->isOnSurface(position))
	    {
		throw cms::Exception("fastsim::BarrelSimplifiedGeometry::getThickness") << "position is not on layer's surface";
	    }
	    return thicknessHist_->GetBinContent(thicknessHist_->GetXaxis()->FindBin(fabs(position.Z())));
	}

	const double getThickness(const math::XYZTLorentzVector & position, const math::XYZTLorentzVector & momentum) const override
	{
	    math::XYZVector normVec(position.Px(), position.Py(), 0.);
	    double fabsCosTheta = fabs(momentum.Vect().Dot(normVec)) / (momentum.P() * normVec.R());
	    return getThickness(position) / fabsCosTheta;
	}
	
	const double getMagneticFieldZ(const math::XYZTLorentzVector & position) const override
	{
	    if(!this->isOnSurface(position))
	    {
		throw cms::Exception("fastsim::BarrelSimplifiedGeometry::getMagneticFieldZ") << "position is not on layer's surface";
	    }
	    return magneticFieldHist_->GetBinContent(magneticFieldHist_->GetXaxis()->FindBin(fabs(position.z())));
	}

	bool isForward() const override 
	{ 
	    return false;
	}

	bool isOnSurface(const math::XYZTLorentzVector & position) const override
	{
	    return fabs(position_ - position.Rho()) < epsilonDistanceR_;
	}
    };

}

#endif
