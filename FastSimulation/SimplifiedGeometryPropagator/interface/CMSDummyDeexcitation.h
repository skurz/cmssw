#ifndef FASTSIM_CMSDUMMYDEEXCITATION_H
#define FASTSIM_CMSDUMMYDEEXCITATION_H

/** 
 * This class is needed as a dummy interface to Geant4  
 * nuclear de-excitation module; no secondary produced
 *
 * \author Vladimir Ivanchenko
 * $Date: 20-Jan-2015
 */ 

#include "G4VPreCompoundModel.hh"
#include "G4ReactionProductVector.hh"

class G4Fragment;
class G4HadFinalState;
class G4HadProjectile;
class G4Nucleus;

namespace fastsim
{
	class CMSDummyDeexcitation : public G4VPreCompoundModel
	{ 
	public:

	  CMSDummyDeexcitation():G4VPreCompoundModel(0, "PRECO") {}; 

	  virtual ~CMSDummyDeexcitation() {};

	  G4HadFinalState* ApplyYourself(const G4HadProjectile&, G4Nucleus&) { return 0; } 

	  G4ReactionProductVector* DeExcite(G4Fragment&) { return new G4ReactionProductVector(); };

	};
}

#endif
