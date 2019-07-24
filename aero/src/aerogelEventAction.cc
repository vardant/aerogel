#include "aerogelEventAction.hh"
#include "aerogelRunAction.hh"

#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandFlat.h"

#include "CLHEP/Random/RandFlat.h"
#include <fstream>
 
aerogelEventAction::aerogelEventAction(aerogelRunAction* run)
:runAct(run)
{
}

aerogelEventAction::~aerogelEventAction()
{
}

void aerogelEventAction::BeginOfEventAction(const G4Event* anEvent)
{
  G4int iev = anEvent->GetEventID();
 G4cout<<"EVENT NUMBER "<<iev+1<<"\n";
 nPMT = 0;
 num_tot = 0.;

 for(G4int i=0;i<2;i++){
   for(G4int j=0;j<8;j++){
     num_gamma[i][j] = 0;
   }
 }

 num_gamma_generated = 0;
 num_gamma_bulk_absorbed = 0;
 num_gamma_boundary_absorbed = 0;
 num_gamma_mesh_absorbed = 0;
 num_gamma_mesh_trapped = 0;
 num_gamma_cathode_absorbed = 0;

} 

void aerogelEventAction::EndOfEventAction(const G4Event* )
{
  for(G4int i=0;i<2;i++){
    for(G4int j=0;j<8;j++){
      if( num_gamma[i][j]>0 ){
       	 nPMT = nPMT + 1;
       }
   }
 } 

outFileaerogel = runAct->GetFileaerogel();
fprintf(outFileaerogel,"%d %.7f %.7f %.7f %.7f %.4f %.d\n",nPMT,Z,Y,psi,phi,mom,particleID);
 for(G4int i=0;i<2;i++){
   for(G4int j=0;j<8;j++){
       if( num_gamma[i][j]>0 ){
	 //	 fprintf(outFileaerogel, "%d %d %d\n",num_gamma[i][j],j+1,i+1);
	 num_tot = num_tot + num_gamma[i][j];
//	  G4cout << "NPE in PMT " << i+1 << " " << j+1 << "    ----    " << num_gamma[i][j] << G4endl; 
      }
   }
 } 


 outFileNPE = runAct->GetFileNPE();
 fprintf(outFileNPE, "%5.f %5d %5d %5d %5d %5d %5d\n", num_tot, num_gamma_generated,
	 num_gamma_bulk_absorbed, num_gamma_boundary_absorbed, num_gamma_mesh_absorbed, num_gamma_mesh_trapped,
	 num_gamma_cathode_absorbed);
 fprintf(outFileaerogel, "%0.f\n",num_tot);
 G4cout<<"Number of photoelectrons from Aerogel: "<<num_tot<<"\n";

}
