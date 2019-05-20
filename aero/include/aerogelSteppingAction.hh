#ifndef aerogelSteppingAction_h
#define aerogelSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "aerogelRunAction.hh"
#include "G4Step.hh"
class aerogelDetectorConstruction;
class aerogelEventAction;

class aerogelSteppingAction : public G4UserSteppingAction
{
  public:

    aerogelSteppingAction(aerogelDetectorConstruction*, aerogelEventAction*);
   ~aerogelSteppingAction();

    void UserSteppingAction(const G4Step*);

public:
G4int QE;
G4double QE0;
G4double PhotEn[104];
G4double qe[104];

private:
aerogelDetectorConstruction* detector;
aerogelEventAction* eventaction;
};

#endif
