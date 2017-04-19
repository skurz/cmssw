import FWCore.ParameterSet.Config as cms

from FastSimulation.Event.ParticleFilter_cfi import  ParticleFilterBlock
from FastSimulation.SimplifiedGeometryPropagator.TrackerMaterial_cfi import TrackerMaterialBlock

fastSimProducer = cms.EDProducer(
    "FastSimProducer",
    src = cms.InputTag("generatorSmeared"),
    particleFilter =  ParticleFilterBlock.ParticleFilter,
    detectorDefinition = TrackerMaterialBlock.TrackerMaterial,
    beamPipeRadius = cms.double(3.),
    interactionModels = cms.PSet(
            bremsstrahlung = cms.PSet(
                className = cms.string("fastsim::Bremsstrahlung"),
                minPhotonEnergy = cms.double(0.1),
                minPhotonEnergyFraction = cms.double(0.005)
                ),
            energyLoss = cms.PSet(
                className = cms.string("fastsim::EnergyLoss")
                ),
            trackerSimHits = cms.PSet(
                className = cms.string("fastsim::TrackerSimHitProducer")
                ),    
            dummyHits = cms.PSet(
                className = cms.string("fastsim::DummyHitProducer")
                ),
            ),
    )
