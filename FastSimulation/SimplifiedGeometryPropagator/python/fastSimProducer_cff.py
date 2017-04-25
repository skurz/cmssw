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
            pairProduction = cms.PSet(
                className = cms.string("fastsim::PairProduction"),
                photonEnergyCut = cms.double(0.1)
                ),
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
                minMomentumCut = cms.double(0.1)
                ),    
            dummyHits = cms.PSet(
                className = cms.string("fastsim::DummyHitProducer")
                ),
            ),
    )
