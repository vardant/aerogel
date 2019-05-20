#ifndef aerogelRunAction_h
#define aerogelRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <fstream>

class G4Run;

class aerogelRunAction : public G4UserRunAction
{
private:
  FILE* photoelectron_numFile;
  FILE* aerogel_outFile;
 
  public:
    aerogelRunAction();
   ~aerogelRunAction();

  public:
  FILE* GetFileNPE(){return photoelectron_numFile;};
  FILE* GetFileaerogel(){return aerogel_outFile;};

    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
};

#endif





