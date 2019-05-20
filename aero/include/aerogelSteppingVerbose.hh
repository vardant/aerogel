class aerogelSteppingVerbose;

#ifndef aerogelSteppingVerbose_h
#define aerogelSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class aerogelSteppingVerbose : public G4SteppingVerbose
{
 public:   

   aerogelSteppingVerbose();
  ~aerogelSteppingVerbose();

   void StepInfo();
   void TrackingStarted();

};

#endif
