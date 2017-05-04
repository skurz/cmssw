import FWCore.ParameterSet.Config as cms

from FastSimulation.Event.ParticleFilter_cfi import  ParticleFilterBlock
from FastSimulation.SimplifiedGeometryPropagator.TrackerMaterial_cfi import TrackerMaterialBlock

fastSimProducer = cms.EDProducer(
    "FastSimProducer",
    src = cms.InputTag("generatorSmeared"),
    particleFilter =  ParticleFilterBlock.ParticleFilter,
    detectorDefinition = TrackerMaterialBlock.TrackerMaterial,
    beamPipeRadius = cms.double(3.),
    deltaRchargedMother = cms.double(0.02), # Maximum angle to associate a charged daughter to a charged mother (mostly done to associate muons to decaying pions)
    interactionModels = cms.PSet(
            pairProduction = cms.PSet(
                className = cms.string("fastsim::PairProduction"),
                photonEnergyCut = cms.double(0.1)
                ),
            nuclearInteraction = cms.PSet(
                className = cms.string("fastsim::NuclearInteraction"),
                distCut = cms.double(0.020),
                hadronEnergy = cms.double(3.5), # the smallest momentum for elastic interactions
                inputFile = cms.string("NuclearInteractionInputFile.txt"), # the file to read the starting interaction in each files (random reproducibility in case of a crash)
                ),
            #nuclearInteractionFTF = cms.PSet(
            #    className = cms.string("fastsim::NuclearInteractionFTF"),
            #    distCut = cms.double(0.020),
            #    bertiniLimit = cms.double(3.5), # upper energy limit for the Bertini cascade 
            #    energyLimit = cms.double(0.1), # Kinetic energy threshold for secondaries 
            #    ),
            bremsstrahlung = cms.PSet(
                className = cms.string("fastsim::Bremsstrahlung"),
                minPhotonEnergy = cms.double(0.1),
                minPhotonEnergyFraction = cms.double(0.005)
                ),
            energyLoss = cms.PSet(
                className = cms.string("fastsim::EnergyLoss"),
                minMomentumCut = cms.double(0.1)
                ),
            multipleScattering = cms.PSet(
                className = cms.string("fastsim::MultipleScattering"),
                minPt = cms.double(0.2)
                ),
            trackerSimHits = cms.PSet(
                className = cms.string("fastsim::TrackerSimHitProducer"),
                minMomentumCut = cms.double(0.1),
                doHitsFromInboundParticles = cms.bool(False), # Track reconstruction not possible for those particles so hits do not have to be simulated
                ),    
            dummyHits = cms.PSet(
                className = cms.string("fastsim::DummyHitProducer")
                ),
            ),
    )
