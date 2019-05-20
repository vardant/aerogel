#include "aerogelRunAction.hh"

#include "G4Run.hh"

#include <fstream>

aerogelRunAction::aerogelRunAction()
{}

aerogelRunAction::~aerogelRunAction()
{}

void aerogelRunAction::BeginOfRunAction(const G4Run* )
{
  photoelectron_numFile = fopen("photEl_num.dat","w");
  aerogel_outFile = fopen("aerogel_out.dat","w");

}

void aerogelRunAction::EndOfRunAction(const G4Run*)
{ 
  fclose(photoelectron_numFile);
  fclose(aerogel_outFile);
}

