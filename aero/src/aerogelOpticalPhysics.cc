#include "aerogelOpticalPhysics.hh"
//#include "aerogelOpPhysMessenger.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   

#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpWLS.hh"

#include "G4VPhysicsConstructor.hh"

aerogelOpticalPhysics::aerogelOpticalPhysics(const G4String& name)
  :  G4VPhysicsConstructor(name)
{
  //  new aerogelOpPhysMessenger(this);
}

aerogelOpticalPhysics::~aerogelOpticalPhysics()
{
}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4OpticalPhoton.hh"

void aerogelOpticalPhysics::ConstructParticle()
{
  G4OpticalPhoton::OpticalPhotonDefinition();
}


#include "G4ProcessManager.hh"

void aerogelOpticalPhysics::ConstructProcess()
{
  theScintProcess = new G4Scintillation();
  theCerenkovProcess=new G4Cerenkov();
  theAbsorptionProcess=new G4OpAbsorption();
  theRayleighScattering=new G4OpRayleigh();
  theBoundaryProcess=new G4OpBoundaryProcess();
  theWLSProcess=new G4OpWLS();

  theWLSProcess->UseTimeProfile("delta");
//theWLSProcess->UseTimeProfile("exponential");

  G4ProcessManager * pManager = 0;
  
  pManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
  pManager->AddDiscreteProcess(theAbsorptionProcess);
  pManager->AddDiscreteProcess(theRayleighScattering);
  ////  theBoundaryProcess->SetModel(unified);
  pManager->AddDiscreteProcess(theBoundaryProcess);
  pManager->AddDiscreteProcess(theWLSProcess);
  
  theScintProcess->SetScintillationYieldFactor(1.);
  theScintProcess->SetScintillationExcitationRatio(0.0);
  theScintProcess->SetTrackSecondariesFirst(true);

  auto theParticleIterator = GetParticleIterator();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    pManager = particle->GetProcessManager();
    if(theCerenkovProcess->IsApplicable(*particle)){
      pManager->AddProcess(theCerenkovProcess);
      pManager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
    }
    if(theScintProcess->IsApplicable(*particle)){
      pManager->AddProcess(theScintProcess);
      pManager->SetProcessOrderingToLast(theScintProcess,idxAtRest);
      pManager->SetProcessOrderingToLast(theScintProcess,idxPostStep);
    }
  
  }
}

void aerogelOpticalPhysics::SetScintYieldFactor(G4double yf){
  if(theScintProcess)
    theScintProcess->SetScintillationYieldFactor(yf);
}
