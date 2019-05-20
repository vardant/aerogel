#include "aerogelPrimaryGeneratorAction.hh"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "CLHEP/Random/RandFlat.h"
#include <sstream>
using std::ifstream;
using namespace CLHEP;

aerogelPrimaryGeneratorAction::aerogelPrimaryGeneratorAction(aerogelEventAction*evt)
  :eventaction(evt)
{
  //part = "opticalphoton";
  //part = "proton";
  //part ="pi-";
  part ="mu-";

  if(part == "mu-"){
    ID = 1;
    ////    mass = 0.105658369; //Gev/c2
  }

  if(part == "pi-"){
    ID = 9;
    ////    mass = 0.13957; //Gev/c2
  }

  if(part == "e-"){
    ID = 3;
    ////    mass = 0.0005109989;  //GeV/c2
  }

  eventaction->SetID(ID);

  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  
}

aerogelPrimaryGeneratorAction::~aerogelPrimaryGeneratorAction()
{
  delete particleGun;
}

void aerogelPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  /*
  //G4double Xpos = -206.*cm; //SHMS
  G4double Xpos = -57.3*cm; //HMS

  momentum = 3.; // Gev/c
  G4int P = G4int(momentum);

  std::ostringstream stm;
  //stm << "random_inputs/p"<<P<<".0";
  stm << string0 << P << ".00";//HMS
  //stm << string0 << P << ".inp";//SHMS

  G4String string1 = stm.str();
  G4cout << "PrimaryGeneratorAction: string1 = " << string1 << G4endl;

  //  Line = RandFlat::shoot(20000);
  G4double L = 20000.*G4UniformRand();
  Line = (G4int)L;


  ifstream indata;
  indata.open(string1);
  if(!indata) { // file couldn't be opened
    G4cerr << "Error: Primary generation data file could not be opened" << G4endl;
    exit(1);
  }


  //  reading from file, HMS
  for(G4int i=0;i<Line;i++){
    indata >> Deltai >> Ypos >> Zpos >> Phi >> Psi;
  }
  */
  /*
  //  reading from file, SHMS
  for(G4int i=0;i<Line;i++){
  indata >> Ypos >> Zpos >> Phi >> Psi >> Deltai;
  }
  */
  /*
  Ypos = -Ypos*cm;
  Zpos =  Zpos*cm;

  momentum = momentum*(1.+Deltai/100.)*GeV;
  ////  ekin = sqrt(mass*mass + momentum*momentum) - mass;

  eventaction->SetY(Ypos);
  eventaction->SetZ(Zpos);
  ////  eventaction->SetMomentum(ekin); ???
  eventaction->SetMomentum(momentum);
  eventaction->SetPhi(-Phi);
  eventaction->SetPsi(Psi);
  */

  //default kinematic
  /*
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(part);
 
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleTime(0.0*ns);
  //  particleGun->SetParticleMomentum(momentum*GeV);
  //  particleGun->SetParticleEnergy(ekin*GeV);
  particleGun->SetParticleEnergy(3.3*GeV);

  //particleGun->SetParticlePosition(G4ThreeVector(Xpos,Ypos,Zpos));
  particleGun->SetParticlePosition(G4ThreeVector(Xpos,0.*cm,0.*cm));
  // particleGun->SetParticlePosition(G4ThreeVector(7.2*cm,0.,-20.*cm));
  //particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,-Phi,Psi));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  //particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  */

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("mu-");
 
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleTime(0.0*ns);
  particleGun->SetParticleMomentum(3.*GeV);

  G4double Xpos = -50.*cm;
  Ypos = (G4UniformRand() - 0.5) * 50.*cm;
  Zpos = (G4UniformRand() - 0.5) * 50.*cm;
  particleGun->SetParticlePosition(G4ThreeVector(Xpos,Ypos,Zpos));

  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));

  particleGun->GeneratePrimaryVertex(anEvent);
}

void aerogelPrimaryGeneratorAction::SetOptPhotonPolar()
{
 G4double angle = G4UniformRand() * 360.0*deg;
 SetOptPhotonPolar(angle);
}

void aerogelPrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
 if (particleGun->GetParticleDefinition()->GetParticleName() != "opticalphoton")
   {
     G4cout << "--> warning from PrimaryGeneratorAction::SetOptPhotonPolar() :"
               "the particleGun is not an opticalphoton" << G4endl;
     return;
   }
     	       
 G4ThreeVector normal (1., 0., 0.);
 G4ThreeVector kphoton = particleGun->GetParticleMomentumDirection();
 G4ThreeVector product = normal.cross(kphoton); 
 G4double modul2       = product*product;
 
 G4ThreeVector e_perpend (0., 0., 1.);
 if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product; 
 G4ThreeVector e_paralle    = e_perpend.cross(kphoton);
 
 G4ThreeVector polar = std::cos(angle)*e_paralle + std::sin(angle)*e_perpend;
 particleGun->SetParticlePolarization(polar);
}
