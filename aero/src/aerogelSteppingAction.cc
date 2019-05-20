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
  /*
    G4double h = 4.13566733e-15;  //eV/c
    G4double c = 2.99796e17; //nm/s

    G4int n = 104;
    G4double wl;
    ifstream Indata_qe;
    Indata_qe.open("inputs/qe.dat");
  
    if(!Indata_qe) { // file couldn't be opened
    G4cerr << "Error: Cannot read QEs" << G4endl;
    exit(1);
    }

    // Reading QEs
    for (G4int k = n - 1; k > -1; k-- ) {
    Indata_qe >> wl >> qe[k];
    PhotEn[k] = h*c/wl; //190 - 880 nm		
    }
  */
  /*
    for (G4int k = 0; k < n; k++) {
    G4cout << "phot en = " << PhotEn[k] << "ev qe = " << qe[k] << G4endl; //190 - 700 nm		
    }
  */


  //G4StepPoint* point2 = aStep->GetPostStepPoint();
  //G4ThreeVector pos2 = point2->GetPosition();
  //G4double zpos2 = pos2.z();
  //G4double ypos2 = pos2.y();


  //  G4StepPoint* point1 = aStep->GetPreStepPoint();
  //G4ThreeVector pos1 = point1->GetPosition();

  G4StepPoint* point2 = aStep->GetPostStepPoint();
  G4ThreeVector pos2 = point2->GetPosition();

  //G4double xpos1 = pos1.x();
  //G4double ypos1 = pos1.y();
  //G4double zpos1 = pos1.z();

  //G4double xpos2 = pos2.x();
  //G4double ypos2 = pos2.y();
  G4double zpos2 = pos2.z();


  //const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
  const G4Track* track = aStep->GetTrack();
  //G4double edep = aStep->GetTotalEnergyDeposit();
  //G4double ekin = track->GetKineticEnergy();
  //G4VPhysicalVolume* volume = track->GetVolume();
  //G4VPhysicalVolume* nextvolume = track->GetNextVolume();
  G4ParticleDefinition* particle_def = track->GetDefinition();

  //// Kill optical photons escaping PMT glass.
  //  if(particle_def == G4OpticalPhoton::OpticalPhotonDefinition()){
  //    G4VPhysicalVolume* prevol = point1->GetPhysicalVolume();
  //    G4VPhysicalVolume* postvol = point2->GetPhysicalVolume();
  //    if (prevol->GetName() == "PMT glass" && postvol->GetName() == "World") {
  //      aStep->GetTrack()->SetTrackStatus(fStopAndKill);
  //    }
  //  }

  //G4ThreeVector dir = track->GetMomentumDirection();
  //G4double x_dir = dir.x();
  //G4double y_dir = dir.y();
  //G4double z_dir = dir.z();

  G4int col,row;

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

  if(particle_def == G4OpticalPhoton::OpticalPhotonDefinition()){

    //  G4VPhysicalVolume* prevol = 
    //  point1->GetPhysicalVolume();

    boundaryStatus=boundary->GetStatus();
 
    //Check to see if the partcile was actually at a boundary

    if(point2->GetStepStatus()==fGeomBoundary){

      //  G4VPhysicalVolume* postvol = 
      // 	point2->GetPhysicalVolume();

    //  if (postvol->GetName() == "World") 
      //	aStep->GetTrack()->SetTrackStatus(fStopAndKill);

      if (boundaryStatus == Detection) {

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
  

/*
 if(detector->WindowEquals(volume)) {
   if(nextvolume == detector->GetCathode() && particle_def->GetParticleName() == "opticalphoton"){
  // G4cout << " Photon at cathode!!!!!!!! "  << G4endl;
  // G4cout << " ekin = " << ekin*1.e+6 << "eV " << G4endl;
     if(zpos2 > 0.){
       col = 0;
     }
     else{
       col = 1;
     }
   row = (ypos2 + (detector->AerogelY() - detector->Y0()))/(detector->DeltaY());
     for(G4int k=0; k < n; k++){
       if(ekin*1.e6 > PhotEn[k] && ekin*1.e6 < PhotEn[k+1]){
	 QE0 = qe[k] + ((qe[k+1] - qe[k])*(ekin*1.e6 - PhotEn[k]))/(PhotEn[k+1] - PhotEn[k]);
	 QE = RandPoisson::shoot(QE0);
	//QE = 1;
	 eventaction->AddNumGamma(col,row,QE);
       }
     }
     if(ekin > PhotEn[n -1]*eV){
       QE = 0;
       eventaction->AddNumGamma(col,row,QE);
     }
   }
 }
*/
//eventaction->AddNumGamma(col,row);

}
