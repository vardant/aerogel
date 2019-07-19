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
using std::ifstream;

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

  G4double zpos2 = pos2.z();

  if(particle_def == G4OpticalPhoton::OpticalPhotonDefinition()) {

    // Kill optical photons trapped in the mesh.
    G4VPhysicalVolume* postvol = point2->GetPhysicalVolume();
    G4VPhysicalVolume* prevol = point1->GetPhysicalVolume();

    //    if (prevol->GetName() == "World") {
    //    if (postvol->GetName() == "World") {
    //      track->SetTrackStatus(fStopAndKill);
    //      G4cout << "*** aerogelSteppingAction::UserSteppingAction: opt. photon killed in the World volume ***" << G4endl;
    //    }

    if (postvol != NULL && postvol->GetName() == "fiber mesh") {
      //      G4cout << "   post.vol. = " << postvol->GetName() << G4endl;
      //      G4cout << "   track length = " << track->GetTrackLength() << G4endl;
      double track_length = track->GetTrackLength();
      //      if (track_length > 100*m) {
      if (track_length > 10*m) {
	track->SetTrackStatus(fStopAndKill);
      	G4cout << "*** aerogelSteppingAction::UserSteppingAction: opt. photon killed, track length = "
      	       << track_length/m << " m ***" << G4endl;
      }
    }

  }
  
  //PMT hits 

  G4OpBoundaryProcessStatus boundaryStatus=Undefined;

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

    boundaryStatus = boundary->GetStatus();
 
    //Check to see if the partcile was actually at a boundary

    if(point2->GetStepStatus()==fGeomBoundary){

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
