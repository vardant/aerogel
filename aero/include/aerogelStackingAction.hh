#ifndef aerogelStackingAction_H
#define aerogelStackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

class aerogelStackingAction : public G4UserStackingAction
{
  public:
    aerogelStackingAction();
   ~aerogelStackingAction();

  public:
    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    void NewStage();
    void PrepareNewEvent();

  private:
    G4int gammaCounter;
};

#endif

