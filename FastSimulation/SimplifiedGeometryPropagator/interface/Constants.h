#ifndef FASTSIM_CONSTANTS_H
#define FASTSIM_CONSTANTS_H


///////////////////////////////////////////////
// Constants
//
// Description: Definition of constants needed for the SimplifiedGeometryPropagator package
//
// Author: L. Vanelderen, S. Kurz
// Date: 29 May 2017
//////////////////////////////////////////////////////////


namespace fastsim
{
    namespace Constants
    {
		static double constexpr speedOfLight = 29.9792458;  // [cm / ns]
		static double constexpr eMass = 0.0005109990615;
		static double constexpr muMass = 0.1056583745;
		static double constexpr epsilonDistance_ =  1.0e-7; // [cm]
		static double constexpr NA = 6.022e+23;
    };
}

#endif
