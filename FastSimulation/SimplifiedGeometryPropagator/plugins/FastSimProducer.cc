// system include files
#include <memory>
#include <string>

// framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// data formats
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/LorentzVector.h"

// fastsim
#include "FastSimulation/Utilities/interface/RandomEngineAndDistribution.h"
#include "FastSimulation/SimplifiedGeometryPropagator/interface/Geometry.h"
#include "FastSimulation/SimplifiedGeometryPropagator/interface/SimplifiedGeometry.h"
#include "FastSimulation/SimplifiedGeometryPropagator/interface/Decayer.h"
#include "FastSimulation/SimplifiedGeometryPropagator/interface/LayerNavigator.h"
#include "FastSimulation/SimplifiedGeometryPropagator/interface/Particle.h"
#include "FastSimulation/SimplifiedGeometryPropagator/interface/ParticleFilter.h"
#include "FastSimulation/SimplifiedGeometryPropagator/interface/InteractionModel.h"
#include "FastSimulation/SimplifiedGeometryPropagator/interface/InteractionModelFactory.h"
#include "FastSimulation/SimplifiedGeometryPropagator/interface/ParticleManager.h"
#include "FastSimulation/Particle/interface/ParticleTable.h" // TODO: move this


///////////////////////////////////////////////
// Author: L. Vanelderen, S. Kurz
// Date: 29 May 2017
//////////////////////////////////////////////////////////


//! The core class of the new SimplifiedGeometryPropagator.
/*!
	Coordinates the propagation of all particles, this means it does the following loop:
    1) Get particle from ParticleManager
    2) Call LayerNavigator to move particle to next intersection with layer
    3) Loop over all the interactions and add secondaries to the event
    4) Repeat steps 2), 3) until particle left the tracker, lost all its energy or is about to decay
    5) If particle is about to decay: do decay and add secondaries to the event
    6) Restart from 1) with the next particle
    7) If last particle was propagated add SimTracks, SimVertices, SimHits,... to the event
*/
class FastSimProducer : public edm::stream::EDProducer<> {
	public:
    explicit FastSimProducer(const edm::ParameterSet&);
    ~FastSimProducer(){;}

	private:
	virtual void beginStream(edm::StreamID id);
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream();

    edm::EDGetTokenT<edm::HepMCProduct> genParticlesToken_; //!< Token to get the genParticles
    fastsim::Geometry geometry_; //!< The definition of the detector according to python config
    double beamPipeRadius_; //!< The radius of the beampipe
    double deltaRchargedMother_;  //!< Cut on deltaR for ClosestChargedDaughter algorithm (FastSim tracking)
    fastsim::ParticleFilter particleFilter_;  //!< Decides which particles have to be propagated
    std::unique_ptr<RandomEngineAndDistribution> _randomEngine;  //!< The random engine
    fastsim::Decayer decayer_;  //!< Handles decays of non-stable particles using pythia
    std::vector<std::unique_ptr<fastsim::InteractionModel> > interactionModels_;  //!< All defined interaction models
    std::map<std::string, fastsim::InteractionModel *> interactionModelMap_;  //!< Each interaction model has a unique name
    edm::IOVSyncValue iovSyncValue_;  //!< The iov thing (interval of validity)
    static const std::string MESSAGECATEGORY;  //!< Category of debugging messages ("FastSimulation")
};

const std::string FastSimProducer::MESSAGECATEGORY = "FastSimulation";

FastSimProducer::FastSimProducer(const edm::ParameterSet& iConfig)
    : genParticlesToken_(consumes<edm::HepMCProduct>(iConfig.getParameter<edm::InputTag>("src"))) 
    , geometry_(iConfig.getParameter<edm::ParameterSet>("detectorDefinition"))
    , beamPipeRadius_(iConfig.getParameter<double>("beamPipeRadius"))
    , deltaRchargedMother_(iConfig.getParameter<double>("deltaRchargedMother"))
    , particleFilter_(iConfig.getParameter<edm::ParameterSet>("particleFilter"))
    , _randomEngine(nullptr)
{

    //----------------
    // define interaction models
    //---------------

    const edm::ParameterSet & modelCfgs = iConfig.getParameter<edm::ParameterSet>("interactionModels");
    for(const std::string & modelName : modelCfgs.getParameterNames())
    {
		const edm::ParameterSet & modelCfg = modelCfgs.getParameter<edm::ParameterSet>(modelName);
		std::string modelClassName(modelCfg.getParameter<std::string>("className"));
		// Use plugin-factory to create model
		std::unique_ptr<fastsim::InteractionModel> interactionModel(fastsim::InteractionModelFactory::get()->create(modelClassName, modelName, modelCfg));
		if(!interactionModel.get()){
			throw cms::Exception("FastSimProducer") << "InteractionModel " << modelName << " could not be created" << std::endl;
		}
		// Add model to list
		interactionModels_.push_back(std::move(interactionModel));
		// and create the map
		interactionModelMap_[modelName] = interactionModels_.back().get();
    }

    //----------------
    // register products
    //----------------

    // SimTracks and SimVertices
    produces<edm::SimTrackContainer>();
    produces<edm::SimVertexContainer>();
    // products of interaction models, i.e. simHits
    for(auto & interactionModel : interactionModels_)
    {
		interactionModel->registerProducts(*this);
    }
}

void
FastSimProducer::beginStream(const edm::StreamID id)
{
    _randomEngine = std::make_unique<RandomEngineAndDistribution>(id);
}

void
FastSimProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    LogDebug(MESSAGECATEGORY) << "   produce";

    // do the iov thing
    if(iovSyncValue_!=iSetup.iovSyncValue())
    {
		LogDebug(MESSAGECATEGORY) << "   triggering update of event setup" << std::endl;
		iovSyncValue_=iSetup.iovSyncValue();
		geometry_.update(iSetup, interactionModelMap_);
    }

    // Define containers for SimTracks, SimVertices
    std::unique_ptr<edm::SimTrackContainer> output_simTracks(new edm::SimTrackContainer);
    std::unique_ptr<edm::SimVertexContainer> output_simVertices(new edm::SimVertexContainer);

    // Get the particle data table (in case lifetime or charge of GenParticles not set)
    edm::ESHandle <HepPDT::ParticleDataTable> pdt;
    iSetup.getData(pdt);
    ParticleTable::Sentry ptable(&(*pdt));

    // Get the GenParticle collection
    edm::Handle<edm::HepMCProduct> genParticles;
    iEvent.getByToken(genParticlesToken_, genParticles);

    // Load the ParticleManager which returns the particles that have to be propagated
    // Creates a fastsim::Particle out of a GenParticle/secondary
    fastsim::ParticleManager particleManager(*genParticles->GetEvent()
											,*pdt
											,beamPipeRadius_
											,deltaRchargedMother_
											,particleFilter_
											,output_simTracks
											,output_simVertices);
	
    LogDebug(MESSAGECATEGORY) << "################################"
			      << "\n###############################";    

	// loop over particles
    for(std::unique_ptr<fastsim::Particle> particle = particleManager.nextParticle(*_randomEngine); particle != 0;particle=particleManager.nextParticle(*_randomEngine)) 
    {
    	LogDebug(MESSAGECATEGORY) << "\n   moving NEXT particle: " << *particle;

		// move the particle through the layers
		fastsim::LayerNavigator layerNavigator(geometry_);
		const fastsim::SimplifiedGeometry * layer = 0;

		// moveParticleToNextLayer(..) returns 0 in case that particle decays
		// in this case particle is propagated up to its decay vertex
		while(layerNavigator.moveParticleToNextLayer(*particle,layer))
		{
		    LogDebug(MESSAGECATEGORY) << "   moved to next layer: " << *layer;
			LogDebug(MESSAGECATEGORY) << "   new state: " << *particle;

		    // perform interaction between layer and particle
		    // do only if there is actual material
		    if(layer->getThickness(particle->position(), particle->momentum()) > 1E-10){
		    	int nSecondaries = 0;
		    	// loop on interaction models
			    for(fastsim::InteractionModel * interactionModel : layer->getInteractionModels())
			    {
					LogDebug(MESSAGECATEGORY) << "   interact with " << *interactionModel;
					std::vector<std::unique_ptr<fastsim::Particle> > secondaries;
					interactionModel->interact(*particle,*layer,secondaries,*_randomEngine);
					nSecondaries += secondaries.size();
					particleManager.addSecondaries(particle->position(),particle->simTrackIndex(),secondaries);
			    }

			    // kinematic cuts: particle might e.g. lost all its energy
			    if(!particleFilter_.acceptsEn(*particle))
			    {	
			    	// Add endvertex if particle did not create any secondaries
			    	if(nSecondaries==0) particleManager.addEndVertex(particle.get());
			    	layer = 0;
			    	break;
			    }
			}

		    // temporary: break after 25 ns or if outermost layer hit
		    // no calorimetry simulated yet!
		    if(particle->position().T() > 25 || particle->position().Perp2() > 119.*119.)
		    {
			    layer = 0;
				break;
		    }
		    
		    LogDebug(MESSAGECATEGORY) << "--------------------------------"
					      << "\n-------------------------------";
		}

		// do decays
		if(!particle->isStable() && particle->remainingProperLifeTimeC() < 1E-10)
		{
		    LogDebug(MESSAGECATEGORY) << "Decaying particle...";
		    std::vector<std::unique_ptr<fastsim::Particle> > secondaries;
		    decayer_.decay(*particle,secondaries, _randomEngine->theEngine());
		    LogDebug(MESSAGECATEGORY) << "   decay has " << secondaries.size() << " products";
		    particleManager.addSecondaries(particle->position(), particle->simTrackIndex(),secondaries);
		}
		
		LogDebug(MESSAGECATEGORY) << "################################"
					  << "\n###############################";
    }

    // store simTracks and simVertices
    iEvent.put(particleManager.harvestSimTracks());
    iEvent.put(particleManager.harvestSimVertices());
    // store products of interaction models, i.e. simHits
    for(auto & interactionModel : interactionModels_)
    {
		interactionModel->storeProducts(iEvent);
    }
}

void
FastSimProducer::endStream()
{
	_randomEngine.reset();
}

DEFINE_FWK_MODULE(FastSimProducer);
