#include "aerogelSteppingAction.hh"
#include "aerogelDetectorConstruction.hh"
#include "aerogelEventAction.hh"
#include "G4Track.hh"
#include "aerogelRunAction.hh"
#include "G4Step.hh"
#include "G4ProcessManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "CLHEP/Random/RandPoisson.h"

#include <math.h>
#include <stdio.h>
#include "G4TouchableHistory.hh"

using namespace CLHEP;
using namespace std;

aerogelSteppingAction::aerogelSteppingAction(aerogelDetectorConstruction*det, aerogelEventAction*evt)
:detector(det), eventaction(evt)
{ }

aerogelSteppingAction::~aerogelSteppingAction()
{ }

void aerogelSteppingAction::UserSteppingAction(const G4Step* aStep)
{

  G4Track* track = aStep->GetTrack();
  G4ParticleDefinition* particle_def = track->GetDefinition();

  G4StepPoint* point1 = aStep->GetPreStepPoint();
  //G4ThreeVector pos1 = point1->GetPosition();

  G4StepPoint* point2 = aStep->GetPostStepPoint();
  G4ThreeVector pos2 = point2->GetPosition();

  G4VPhysicalVolume* prevol = point1->GetPhysicalVolume();
  G4VPhysicalVolume* postvol = point2->GetPhysicalVolume();

  G4double zpos2 = pos2.z();

  // Count generated optical photons.

  if(particle_def != G4OpticalPhoton::OpticalPhotonDefinition()) {

    const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep();

    if (secondaries->size()>0) {

      for (unsigned int i=0; i<secondaries->size(); ++i) {

        if (secondaries->at(i)->GetParentID()>0) {
	  if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()
	     == G4OpticalPhoton::OpticalPhotonDefinition()) {

	    //	    if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
	    //		== "Scintillation")fScintillationCounter++;
	    //	    if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
	    //		== "Cerenkov")fCerenkovCounter++;

	    if (postvol != NULL && postvol->GetName() != "World") {
	      eventaction->AddNumGammaGenerated(1);
	      //	      cout << " ==> Optical photon generated" << endl;
	    }

	  }
        }
      }
    }
  }

  if(particle_def == G4OpticalPhoton::OpticalPhotonDefinition()) {

    auto step_process = (G4VProcess*)point2->GetProcessDefinedStep();
    auto proc_name = step_process->GetProcessName(); 
    if (proc_name == "OpAbsorption") {
      eventaction->AddNumGammaBulkAbsorbed(1);
      //      cout << " ==> Optical photon absorption: prevol = " << prevol->GetName();
      //      cout << ", postvol = " << postvol->GetName() << endl;
    }

    // Kill optical photons trapped in the mesh.
    /*
    if (postvol != NULL && postvol->GetName() == "fiber mesh") {
      double track_length = track->GetTrackLength();
      if (track_length > 10*m) {
	track->SetTrackStatus(fStopAndKill);
	eventaction->AddNumGammaMeshTrapped(1);
	//      G4cout << "   post.vol. = " << postvol->GetName() << G4endl;
	//      G4cout << "   track length = " << track->GetTrackLength() << G4endl;
	//      	G4cout << "*** aerogelSteppingAction::UserSteppingAction: opt. photon killed, track length = "
	//      	       << track_length/m << " m ***" << G4endl;
      }
    }
    */

  }
  
  //PMT hits 

  static G4OpBoundaryProcess* boundary=NULL;
  
  //find the boundary process only once
  if(!boundary){
    G4ProcessManager* pm 
      = track->GetDefinition()->GetProcessManager();
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    G4int i;
    for( i=0;i<nprocesses;i++){
      if((*pv)[i]->GetProcessName()=="OpBoundary"){
        boundary = (G4OpBoundaryProcess*)(*pv)[i];
        break;
      }
    }
  }

  if(particle_def == G4OpticalPhoton::OpticalPhotonDefinition()) {

    //Check to see if the partcile was actually at a boundary

    if(point2->GetStepStatus()==fGeomBoundary){

      //    G4OpBoundaryProcessStatus boundaryStatus=Undefined;
      G4OpBoundaryProcessStatus boundaryStatus = boundary->GetStatus();

      if (boundaryStatus == Absorption) {
	if (postvol != NULL) {
	  if (postvol->GetName() == "cathode_phys")
	    eventaction->AddNumGammaCathodeAbsorbed(1);
	  else if (postvol->GetName() == "fiber mesh") {
	    eventaction->AddNumGammaMeshAbsorbed(1);
	    cout << "Absorption at mesh: " << pos2.x()/cm << " " << pos2.y()/cm << " " << pos2.z()/cm << endl;
	  }
	  else {
	    eventaction->AddNumGammaBoundaryAbsorbed(1);

	    //	cout << "OpPhoton absorption at boundary: prevol = " << prevol->GetName()
	    //	     << "  postvol = ";
	    //	postvol != NULL ? cout << postvol->GetName() : cout << "NULL";
	    //	cout << endl;
	  }
	}
      }

      if (postvol != NULL && postvol->GetName() == "fiber mesh"
	  && boundaryStatus == TotalInternalReflection) {
	track->SetTrackStatus(fStopAndKill);
	eventaction->AddNumGammaMeshTrapped(1);
	//      G4cout << "   post.vol. = " << postvol->GetName() << G4endl;
	//      G4cout << "*** Opt. photon trapped in mesh, killed ***" G4endl;
      }
 
      if (boundaryStatus == Detection) {

	G4int col,row;

	if(zpos2 > 0.){
	  col = 0;
	}
	else{
	  col = 1;
	}

	//G4double row_d = (ypos2 + (detector->AerogelY() - detector->Y0()))/(detector->DeltaY());
	// G4int row= G4int(row_d);
	row = 0;
	eventaction->AddNumGamma(col,row);
	// G4cout << "Detection!!" << G4endl;
      }
    }
  
  }   //optical photon

}
