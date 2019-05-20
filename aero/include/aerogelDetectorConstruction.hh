#ifndef aerogelDetectorConstruction_h
#define aerogelDetectorConstruction_h 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"
class G4LogicalVolume;
class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"
using std::string;

class aerogelDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    aerogelDetectorConstruction();
   ~aerogelDetectorConstruction();

  public:
    G4VPhysicalVolume* Construct();


const G4VPhysicalVolume* GetWorld() {return expHall_phys;};
const G4VPhysicalVolume* GetCathode() {return photocathode_phys;};
bool WindowEquals(G4VPhysicalVolume *v) {
  bool result=false;
  for (G4int i = 0; i < 8; i++)
  {
	for(G4int k = 0; k < 2; k++){	
    		if (PMT_glass_phys[i][k]==v) result = true;
	}
  }
  return result;
}


  G4double AerogelY() {return hight_wrap;};
  G4double Y0() {return y0;};
  G4double DeltaY() {return d;};

  private:
    G4double expHall_x;
    G4double expHall_y;
    G4double expHall_z;

G4double hight_wrap; // x distance from back edge of the box to the holes center
G4double y0; // y distance from bottom edge to first hole center
G4double d;

G4double PhotEn[104]; //for cathode
G4double qe[104];
G4double RIndex_cathode[104];
G4double RIndex_Aerogel[35];
G4double RIndex_Air[35];
G4double RIndex_glass[35];
G4double Absorption[35];
G4double Scattering[35];
G4double Scattering1[35];
G4double PhotonEnergy[35];
G4MaterialPropertiesTable* MPT_Aerogel;
G4MaterialPropertiesTable* MPT_glass;

G4VPhysicalVolume* photocathode_phys_test;
G4VPhysicalVolume* PMT_glass_phys_test;
G4VPhysicalVolume* millipore_aer_phys;
G4VPhysicalVolume* millipore_dif_phys;
G4VPhysicalVolume* millipore_wrap_phys;
G4VPhysicalVolume* photocathode_phys;
G4VPhysicalVolume* Dif_phys;
G4VPhysicalVolume* Dif_side_phys_right;
G4VPhysicalVolume* Dif_side_phys_left;
G4VPhysicalVolume* Gore_side_phys_right;
G4VPhysicalVolume* Gore_side_phys_left;
G4VPhysicalVolume* Dif_top_phys;
G4VPhysicalVolume* Dif_bottom_phys;
G4VPhysicalVolume* Al_back_phys;
G4VPhysicalVolume* HC_phys;
G4VPhysicalVolume* A_side_phys_right;
G4VPhysicalVolume* A_side_phys_left;
G4VPhysicalVolume* A_top_phys;
G4VPhysicalVolume* A_bottom_phys;
G4VPhysicalVolume* PMT_glass_phys[8][2];
G4VPhysicalVolume* Aerogel_phys;
G4VPhysicalVolume* Aerogel_box_phys;
G4VPhysicalVolume* expHall_phys;
G4VPhysicalVolume* test_phys;
G4VPhysicalVolume* Gore_phys;
G4VPhysicalVolume* Gore_top_phys;
G4VPhysicalVolume* Gore_bot_phys;
G4VPhysicalVolume* cath_black_phys;
G4VPhysicalVolume* fiber_mesh_z_phys;
G4VPhysicalVolume* fiber_mesh_y_phys;
G4LogicalSkinSurface* MilliporeSurface;	
G4LogicalSkinSurface* MilliporeSurface1;
G4LogicalSkinSurface* MilliporeSurface2;
G4LogicalSkinSurface* GoreSurface;
G4LogicalSkinSurface* GoreSurface1;
G4LogicalSkinSurface* GoreSurface2;
};

#endif /*aerogelDetectorConstruction_h*/
