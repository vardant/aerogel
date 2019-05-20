#include "aerogelStackingAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"

aerogelStackingAction::aerogelStackingAction()
: gammaCounter(0)
{}

aerogelStackingAction::~aerogelStackingAction()
{}

G4ClassificationOfNewTrack
aerogelStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
  { // particle is optical photon
    if(aTrack->GetParentID()>0)
    { // particle is secondary
      gammaCounter++;
    }
  }
  return fUrgent;
}

void aerogelStackingAction::NewStage()
{
  G4cout << "Total number of optical photons produced : "
         << gammaCounter << G4endl;
}

void aerogelStackingAction::PrepareNewEvent()
{ gammaCounter = 0; }
