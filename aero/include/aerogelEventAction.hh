#ifndef aerogelEventAction_h
#define aerogelEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class aerogelRunAction;

//class G4Event;

class aerogelEventAction : public G4UserEventAction
{
  public:
    aerogelEventAction(aerogelRunAction*);
   ~aerogelEventAction();

  //void AddNumGamma(G4int col, G4int row, G4int QEff){num_gamma[col][row] = num_gamma[col][row]+QEff + 1;};
  void AddNumGamma(G4int col, G4int row){num_gamma[col][row] = num_gamma[col][row]+1;};

  void AddNumGammaGenerated(int n) {num_gamma_generated += n;};
  void AddNumGammaBulkAbsorbed(int n) {num_gamma_bulk_absorbed += n;};
  void AddNumGammaBoundaryAbsorbed(int n) {num_gamma_boundary_absorbed += n;};
  void AddNumGammaMeshAbsorbed(int n) {num_gamma_mesh_absorbed += n;};
  void AddNumGammaMeshTrapped(int n) {num_gamma_mesh_trapped += n;};
  void AddNumGammaCathodeAbsorbed(int n) {num_gamma_cathode_absorbed += n;};

  void BeginOfEventAction(const G4Event* anEvent);
  void EndOfEventAction(const G4Event* anEvent);

  void SetMomentum (G4double momentum) {mom = momentum;};
  //void SetEnergy (G4double ekin) {en = ekin/1000.;};
  void SetY (G4double Ypos) {Y = Ypos/10.;};
  void SetZ (G4double Zpos) {Z = Zpos/10.;};
  void SetID (G4int ID) {particleID = ID;};		
  void SetPhi (G4double Phi) {phi = Phi;};	
  void SetPsi (G4double Psi) {psi = Psi;};
	
private:

   G4double phi;
   G4double psi;
   G4int  particleID;
   G4double en;
   G4double mom;
   G4double Y;
   G4double Z;

  FILE* outFileNPE;
  FILE* outFileaerogel;

  G4int num_gamma[2][8];
  G4int nPMT;

  G4double num_tot;
  int num_gamma_generated;
  int num_gamma_bulk_absorbed;
  int num_gamma_boundary_absorbed;
  int num_gamma_mesh_absorbed;
  int num_gamma_mesh_trapped;
  int num_gamma_cathode_absorbed;

   aerogelRunAction* runAct;
};

#endif
