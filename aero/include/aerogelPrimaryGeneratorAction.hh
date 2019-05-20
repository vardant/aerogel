#ifndef aerogelPrimaryGeneratorAction_h
#define aerogelPrimaryGeneratorAction_h 1
#include "aerogelEventAction.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include <string>
using std::string;

class G4ParticleGun;
class G4Event;
//class aerogelPrimaryGeneratorMessenger;

class aerogelPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    aerogelPrimaryGeneratorAction(aerogelEventAction*);
   ~aerogelPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);

    void SetOptPhotonPolar();
    void SetOptPhotonPolar(G4double);

public:
   G4int Line;
   G4int Ev_number;
   G4double Phi;
   G4double Psi;
   G4double Deltai;
   char Boarder;
string string0;
  G4double ekin;
  G4double Ypos;
  G4double Zpos;
  G4int ID;
  string part;
  G4double mass;
  G4double momentum;

  private:
    G4ParticleGun* particleGun;
    //aerogelPrimaryGeneratorMessenger* gunMessenger;
    aerogelEventAction* eventaction;
};

#endif /*aerogelPrimaryGeneratorAction_h*/
