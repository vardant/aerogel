#include "aerogelDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Ellipsoid.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4MultiUnion.hh"
#include "G4UnionSolid.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4PVReplica.hh"
#include "G4NistManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using std::ifstream;
using std::string;

using namespace std;

aerogelDetectorConstruction::aerogelDetectorConstruction()
{
  expHall_x = 2.5*m;
  expHall_y = expHall_z = 0.8*m;
  
}

aerogelDetectorConstruction::~aerogelDetectorConstruction(){;}

G4VPhysicalVolume* aerogelDetectorConstruction::Construct()
{
  const G4double h = 4.13566733e-15;  //eV/c
  const G4double c = 2.99796e+17; //nm/s

  const G4double Aer_depth = 9./2.*cm; //Aerogel thickness

  const G4double wrap_thick = 0.01*cm; //Millipore paper thickness
  d = 5.875*2.54*cm;; // distance between PMTs centers
  const G4double PMT_r = 5.*2.54/2.*cm; //PMT radius
  const G4double PMT_curve_r = 13.*cm; //PMT curve radius
  const G4double PMT_glass_thick = 2.*mm; // PMT glass thickness
  const G4double cathode_r = 11./2.*cm; //photocathode radius

  const G4double hole_r = 5.487*2.54/2.*cm; // radius of holes for PMTs in Al support structure


  //---------SHMS configuration
  const G4double Ax = 4.267*2.54/2.*cm; //Aerogel box depth ----------------------------------------------SHMS
  const G4double Dx = 6.742*2.54/2.*cm;  //Diffusion box depth--------------------------------------------SHMS
  const G4double hight = 44.5*2.54/2.*cm; //Detector hight (without Al support structure and wrapping)----SHMS
  const G4double width = 40.05*2.54/2.*cm; //Detector width (without Al support structure and wrapping)---SHMS
  y0 = 4.625*2.54*cm; // y distance from bottom edge to first hole center------------------SHMS
  const G4double hole_pos_x = 3.25*2.54*cm; // x distance from back edge of the box to the holes center---SHMS
  const G4int n_PMT = 7; //number of PMTs at each side----------------------------------------------------SHMS
  const G4String qe_file = "./inputs/qe_xp4500.dat"; // name of file with qe fot PMTs-----SHMS
  const G4int n_qe = 104; //number of qe values in file for xp4500 PMT-----------------------------------SHMS

/*
  //---------HMS configurarion
//G4double Ax = 9.5/2.*cm; //Aerogel box depth -------------------------------------------------HMS
  G4double Ax = 9.271/2.*cm; //Aerogel box depth -------------------------------------------------HMS 
  G4double Dx = 16.83/2.*cm;  //Diffusion box depth-----------------------------------------------HMS
  //G4double Dx = 16.83/2.*cm;  //Diffusion box depth-----------------------------------------------HMS
  G4double hight = 120./2.*cm; //Detector hight (without Al support structure and wrapping)-----HMS
  G4double width = 70./2.*cm; //Detector width (without Al support structure and wrapping)------HMS
  G4double y0 = 8.*cm; // y distance from bottom edge to first hole center----------------------HMS
//  G4double hole_pos_x = 12.5*cm; // x distance from back edge of the box to the holes center----HMS
  //G4double hole_pos_x = Dx; // x distance from back edge of the box to the holes center----HMS
   G4double hole_pos_x = 16.83/2.*cm;
  G4int n_PMT = 8; //number of PMTs at each side------------------------------------------------HMS
  G4String qe_file = "/w/hallc/shms/simon/geant4/g4work/aerogel/inputs/qe_xp4572.dat"; // name of file with qe fot PMTs------HMS
  G4int n_qe = 79; //number of qe values in file for xp4572 PMT---------------------------------HMS
*/

  //Fiber mesh.
  struct {
    ///    const double fiber_diam = 2.*mm;
    const double fiber_diam = 0.1*mm;
    const double step_x = 1.*cm;     //keep constant
    const double step_y = 0.05*cm;
    const double step_z = 0.05*cm;
  } mesh;

  cout << "aerogelDetectorConstruction::Construct: fiber mesh sizes:" << endl;
  cout << "  fiber_diam = " << mesh.fiber_diam/mm << " mm" << endl;
  cout << "  step_x = " << mesh.step_x/cm << " cm" << endl;
  cout << "  step_y = " << mesh.step_y/cm << " cm" << endl;
  cout << "  step_z = " << mesh.step_z/cm << " cm" << endl;

 
  //	------------- Materials -------------
  G4double a, z, density;
  G4int nelements,natoms;
  
  //Elements
  G4Element* O  = new G4Element("Oxygen"  ,"O" , z= 8., a= 16.00*g/mole);
  G4Element* SI = new G4Element("Silicium", "SI", z=14 , a=28.088*g/mole);
  G4Element* NA = new G4Element("Natrium", "NA", z=11 , a=22.990*g/mole);
  G4Element* CA = new G4Element("Calcium", "Ca", z=20 , a=40.8*g/mole); 
  G4Element* H = new G4Element("Hydrogen", "H", z=1, a=1.008*g/mole);
  G4Element* C = new G4Element("Carbon", "C", z=6, a=12.011*g/mole);
  
  //SiO2
  G4Material* SiO2 = new G4Material ("Si Oxyd", density = 2.634*g/cm3,  nelements = 2);
  SiO2->AddElement(SI, natoms = 1);
  SiO2->AddElement(O, natoms = 2);
  
  //H2O
  G4Material* H2O = new G4Material("Water", 1.000*g/cm3, 2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  
  //Millipore (paper)
  G4Material* Millipore = new G4Material("Millipoe", 1.000*g/cm3, 3);
  Millipore->AddElement(C, natoms=6);
  Millipore->AddElement(H, natoms=10);
  Millipore->AddElement(O, natoms=5);
  
  //Aerogel
  G4Material* Aerog = new G4Material("Aerogel", 0.1100*g/cm3, 3);
  Aerog->AddMaterial(SiO2, 62.5*perCent);
  Aerog->AddMaterial(H2O , 37.4*perCent);
  Aerog->AddElement (C , 0.1*perCent);
  
  //Air
  G4NistManager* man = G4NistManager::Instance();
  G4Material* Air = man->FindOrBuildMaterial("G4_AIR");

  //Argentum
  G4Material* Ag = new G4Material("Ag", z= 47., a= 107.868*g/mole, density= 10.5*g/cm3);

  G4Material* CarbonMaterial = man->FindOrBuildMaterial("G4_C");
  
  G4Material* CopperMaterial = man->FindOrBuildMaterial("G4_Cu");
  
  //glass
  G4Material* glass = new G4Material("glass", density=2.4*g/cm3, nelements=4);
  glass->AddElement(O, 45.98*perCent);
  glass->AddElement(NA, 9.64412*perCent);
  glass->AddElement(SI, 33.6553*perCent);
  glass->AddElement(CA, 10.7205*perCent);
  
  // ------------ Generate & Add Material Properties Table ------------
  const G4int nEntries = 35;
  
  //---Aerogel
  for(int i = 0; i < nEntries; i++){
    RIndex_Aerogel[i] = 1.0304;// SP30
    //	RIndex_Aerogel[i] = 1.015;// SP15
  }
  
  G4double wl;
  ifstream Indata_abs, Indata_scat;
  Indata_abs.open("./inputs/abs_length1.dat");
  Indata_scat.open("./inputs/scat_length.dat");

  if(!Indata_abs || !Indata_scat) { // file couldn't be opened
    G4cerr << "Error: Cannot read attenuation parameters" << G4endl;
    exit(1);
  }
  
  for (G4int k = nEntries - 1; k > -1; k-- ) {
    Indata_abs >> wl >> Absorption[k]; //from Aschenauer for one block, differs from mean absorption lengths values, not used further
    Indata_scat >> wl >> Scattering[k];//from Aschenauer for one block, close to mean values, may be used 
    PhotonEnergy[k] = h*c/(wl*1000.)*eV; //220 - 880 nm
    Absorption[k] *= 2.*cm;
    Scattering[k] *= cm;
    Scattering1[k] = 1.125/0.0094*(wl*wl*wl*wl)*cm; // from fit with C*tmean = 0.0094 mkm4, block thickness = 1.125cm		
  }
  
  G4double PhEnTest[3] = {h*c/880.*eV, h*c/300.*eV, h*c/190.*eV};
  //G4double AbsorptionTest[3] = {30.7*cm, 30.7*cm, 1.*cm};
    G4double AbsorptionTest[3] = {110.*cm, 110.*cm, 1.*cm};
  // G4double AbsorptionTest[3] = {90.*cm, 90.*cm, 1.*cm};

  MPT_Aerogel = new G4MaterialPropertiesTable();
  MPT_Aerogel->AddProperty("RINDEX", PhotonEnergy, RIndex_Aerogel, nEntries);
  MPT_Aerogel->AddProperty("ABSLENGTH", PhEnTest, AbsorptionTest, 3);
  MPT_Aerogel->AddProperty("RAYLEIGH", PhotonEnergy, Scattering1, nEntries);
  Aerog->SetMaterialPropertiesTable(MPT_Aerogel);	
 
  // ----Air
  for(int i = 0; i < nEntries; i++){
    RIndex_Air[i] = 1.000293;
    //  RIndex_Air[i] = 1.0;
  }
  G4MaterialPropertiesTable* MPT_Air = new G4MaterialPropertiesTable();
  MPT_Air->AddProperty("RINDEX", PhotonEnergy, RIndex_Air, nEntries);
  Air->SetMaterialPropertiesTable(MPT_Air);
  
  //---Glass
  for (int i = 0; i < nEntries; i++){
    RIndex_glass[i] = 1.54;
  }
  MPT_glass = new G4MaterialPropertiesTable();
  MPT_glass->AddProperty("RINDEX", PhotonEnergy, RIndex_glass,nEntries);
  glass->SetMaterialPropertiesTable(MPT_glass);

  //Carbon optics MPT.

  const int ndat_c = 92;
  G4double energy_c[ndat_c];
  G4double refind_c[ndat_c];
  G4double extind_c[ndat_c];

  ifstream Indata_c;
  Indata_c.open("./inputs/carbon_optic.dat");

  if(!Indata_c) { // file couldn't be opened
    G4cerr << "Error: Cannot read carbon optics data" << G4endl;
    exit(1);
  }
  
  for (G4int k = ndat_c - 1; k > -1; k-- ) {
    double wl_c;
    Indata_c >> wl_c >> refind_c[k] >> extind_c[k];
    energy_c[k] = h*c/(wl_c*1000.)*eV;
  }

  Indata_c.close();

  G4MaterialPropertiesTable* MPT_Carbon = new G4MaterialPropertiesTable();
  MPT_Carbon->AddProperty("REALRINDEX", energy_c, refind_c, ndat_c);
  MPT_Carbon->AddProperty("IMAGINARYRINDEX", energy_c, extind_c, ndat_c);
  //  MPT_Carbon->DumpTable();
  //  getchar();

  
  //Copper optics MPT.

  const int ndat_cu = 43;
  G4double energy_cu[ndat_cu];
  G4double refind_cu[ndat_cu];
  G4double extind_cu[ndat_cu];

  ifstream Indata_cu;
  Indata_cu.open("./inputs/copper_optic.dat");

  if(!Indata_cu) { // file couldn't be opened
    G4cerr << "Error: Cannot read copper optics data" << G4endl;
    exit(1);
  }
  
  for (G4int k = ndat_cu - 1; k > -1; k-- ) {
    double wl_cu;
    Indata_cu >> wl_cu >> refind_cu[k] >> extind_cu[k];
    energy_cu[k] = h*c/(wl_cu*1000.)*eV;
  }

  Indata_cu.close();

  G4MaterialPropertiesTable* MPT_Copper = new G4MaterialPropertiesTable();
  MPT_Copper->AddProperty("REALRINDEX", energy_cu, refind_cu, ndat_cu);
  MPT_Copper->AddProperty("IMAGINARYRINDEX", energy_cu, extind_cu, ndat_cu);
  MPT_Copper->DumpTable();
  //  getchar();

  
//--------Build detector-------------------// 

  //--------------- The experimental Hall
  G4Box* expHall_box = new G4Box("World",expHall_x,expHall_y,expHall_z);
  G4LogicalVolume* expHall_log = new G4LogicalVolume(expHall_box,Air,"World",0,0,0);
  expHall_phys = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);
  
  //----------PMTs------------//
  //PMT glass log
  G4double zTopCutGlass = -PMT_curve_r*std::cos(std::asin(PMT_r/PMT_curve_r));
  G4Ellipsoid* PMT_glass = new G4Ellipsoid("glass",PMT_curve_r, PMT_curve_r, 
						PMT_curve_r, -PMT_curve_r, zTopCutGlass);  
  G4LogicalVolume *PMT_glass_log = new G4LogicalVolume(PMT_glass, glass, "PMT glass", 0,0,0);
  
  //Photocathode
  G4double cathode_curve_r = PMT_curve_r - PMT_glass_thick;
  G4double zTopCutCathode = -cathode_curve_r*std::cos(std::asin(cathode_r/cathode_curve_r));
  G4Ellipsoid* photocathode = new G4Ellipsoid("cathode",
			cathode_curve_r, cathode_curve_r, cathode_curve_r, -cathode_curve_r, zTopCutCathode);
  G4LogicalVolume *photocathode_log = new G4LogicalVolume(photocathode,Ag,"photocathode test", 0,0,0);
  photocathode_phys = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.*cm),
					photocathode_log, "cathode_phys",  PMT_glass_log, false, 0);
  
  G4double black_thick = 0.01*mm;   //To not let photons out of the detector
  ////  G4Tubs* cath_black = new G4Tubs("cath black",0., PMT_r, black_thick,0.*deg,360.*deg);
  G4Ellipsoid* cath_black = new G4Ellipsoid("glass",PMT_curve_r, PMT_curve_r, 
					    PMT_curve_r, zTopCutGlass - black_thick, zTopCutGlass);  
  G4LogicalVolume *cath_black_log = new G4LogicalVolume(cath_black,Ag,"cath black log", 0,0,0);
  ////  cath_black_phys = new G4PVPlacement(0,G4ThreeVector(0.,0.,zTopCutGlass - black_thick),
  ////				      cath_black_log, "cathode_black_phys",  PMT_glass_log, false, 0);
  cath_black_phys = new G4PVPlacement(0,G4ThreeVector(),
				      cath_black_log, "cathode_black_phys",  PMT_glass_log, false, 0);
    
  //-------Millipore
  hight_wrap = hight + wrap_thick;
  G4double width_wrap = width + wrap_thick;
  G4double depth_wrap = Ax + Dx + wrap_thick;  //Detector depth
  
  ////  G4Box *Millipore_wrap_box = new G4Box("Mil wrap box", depth_wrap, hight_wrap, width_wrap);
  ////  G4Tubs *hole = new G4Tubs("hole for pmt", 0., PMT_r, width_wrap, 0.*deg, 360.*deg);
  ////  G4Tubs *hole = new G4Tubs("hole for pmt", 0., PMT_r, width_wrap+1.*mm, 0.*deg, 360.*deg);
  ////  G4double PMT_pos_x = depth_wrap - hole_pos_x;
  ////  
  ////  G4double PMT_pos_y0 = - hight_wrap + y0;
  ////  G4ThreeVector trans0(PMT_pos_x, PMT_pos_y0, 0.);
  ////  G4SubtractionSolid *Millipore_sub = new G4SubtractionSolid("Millipore wrapping sub",
  ////							 Millipore_wrap_box, hole, 0, trans0);
  ////  
  ////  for(G4int i  = 0; i < n_PMT - 1; i++){
  ////    G4ThreeVector trans_hole(PMT_pos_x, PMT_pos_y0 + d*(i+1), 0.);
  ////    Millipore_sub = new G4SubtractionSolid("Millipore sub", Millipore_sub, hole, 0, trans_hole);
  ////  }

  G4double PMT_pos_x = depth_wrap - hole_pos_x;
  G4double PMT_pos_y0 = - hight_wrap + y0;

  G4MultiUnion *holes = new G4MultiUnion("holes solid");
  G4Tubs *hole = new G4Tubs("hole for pmt", 0., PMT_r, width_wrap+1.*mm, 0.*deg, 360.*deg);

  for (int i=0; i<n_PMT; i++) {
    G4RotationMatrix rotm  = G4RotationMatrix();
    G4ThreeVector trans_hole(PMT_pos_x, PMT_pos_y0 + d*i, 0.);
    G4Transform3D tr = G4Transform3D(rotm,trans_hole);
    holes->AddNode(*hole,tr);
  }

  holes->Voxelize();

  G4Box *Millipore_wrap_box = new G4Box("Mil wrap box", depth_wrap, hight_wrap, width_wrap);
  G4SubtractionSolid *Millipore_sub = new G4SubtractionSolid("Millipore wrapping sub",
							     Millipore_wrap_box, holes, 0, G4ThreeVector());

  G4Box *dif_box = new G4Box("dif box", Dx, hight, width);
  G4ThreeVector trans_dif(depth_wrap - wrap_thick - Dx, 0., 0.);
  Millipore_sub = new G4SubtractionSolid("Millipore sub", Millipore_sub, dif_box, 0, trans_dif);

  G4LogicalVolume *millipore_wrap_log = new G4LogicalVolume(Millipore_sub, Millipore, "millipore dif", 0,0,0);
  millipore_wrap_phys = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.),
					  millipore_wrap_log, "millipore wrap", expHall_log, false, 0); 

  //---------Aerogel box
  G4Box *aerogel_box = new G4Box("aerogel box", Ax, hight, width);
  G4LogicalVolume *aerogel_box_log = new G4LogicalVolume(aerogel_box, Air, "aerogel  log", 0,0,0);
  Aerogel_box_phys = new G4PVPlacement(0, G4ThreeVector(-depth_wrap + wrap_thick + Ax, 0., 0.),
				       aerogel_box_log, "Aerogel box phys", millipore_wrap_log, false, 0);

  //---------Aerogel
  G4double Aer_hight = hight;
  G4double Aer_width = width;
  G4Box *aerogel = new G4Box("aerogel box", Aer_depth, Aer_hight, Aer_width);
  G4LogicalVolume *aerogel_log = new G4LogicalVolume(aerogel, Aerog, "aerogel  log", 0,0,0);

  //Fiber mesh in aerogel.

  //With overlays, but can be visualized.
  /*
  if (mesh.step_x !=0. && mesh.step_y != 0. && mesh.step_z !=0) {

    G4MultiUnion *fiber_mesh_z = new G4MultiUnion("fiber mesh Z");
    G4Tubs *fiber_z = new G4Tubs("fiber Z", 0., mesh.fiber_diam/2., Aer_width - 0.1*mm, 0.*deg, 360.*deg);

    //    double xf = -Aer_depth + 0.1*mm + mesh.fiber_diam/2.;
    //    while (xf <= Aer_depth - 0.1*mm - mesh.fiber_diam/2.) {
    double xf = -Aer_depth + 0.5*cm;
    while (xf < Aer_depth) {

      double yf = -Aer_hight + 0.1*mm + mesh.fiber_diam/2.;
      while (yf <= Aer_hight - 0.1*mm - mesh.fiber_diam/2.) {
	G4RotationMatrix rotm;
	G4ThreeVector trans_f(xf, yf, 0.);
	G4Transform3D tr_f = G4Transform3D(rotm,trans_f);
	fiber_mesh_z->AddNode(*fiber_z,tr_f);
	yf += mesh.step_y;
      }

      xf += mesh.step_x;
    }

    fiber_mesh_z->Voxelize();

    G4MultiUnion *fiber_mesh_y = new G4MultiUnion("fiber mesh Y");
    G4Tubs *fiber_y = new G4Tubs("fiber Y", 0., mesh.fiber_diam/2., Aer_hight - 0.1*mm, 0.*deg, 360.*deg);

    //    xf = -Aer_depth + 0.1*mm + mesh.fiber_diam/2.;
    xf = -Aer_depth + 0.1*mm + mesh.fiber_diam + 0.001*mm;
    while (xf <= Aer_depth - 0.1*mm - mesh.fiber_diam/2.) {

      double zf = -Aer_width + 0.1*mm + mesh.fiber_diam/2.;
      while (zf <= Aer_width - 0.1*mm - mesh.fiber_diam/2.) {
	G4RotationMatrix rotm;
	rotm.rotateX(90.*deg);
	G4ThreeVector trans_f(xf, 0., zf);
	G4Transform3D tr_f = G4Transform3D(rotm,trans_f);
	fiber_mesh_y->AddNode(*fiber_y,tr_f);
	zf += mesh.step_z;
      }

      xf += mesh.step_x;
    }

    fiber_mesh_y->Voxelize();

    //    G4UnionSolid* fiber_mesh = new G4UnionSolid("fibermesh", fiber_mesh_z, fiber_mesh_y);
    //    G4LogicalVolume *fiber_mesh_log = new G4LogicalVolume(fiber_mesh, Air, "fiber mesh log", 0,0,0);
    //    fiber_mesh_z_phys = new G4PVPlacement(0,G4ThreeVector(),fiber_mesh_log, "fiber mesh", aerogel_log, false, 0);

    G4LogicalVolume *fiber_mesh_z_log = new G4LogicalVolume(fiber_mesh_z, Air, "fiber mesh Z log", 0,0,0);
    G4LogicalVolume *fiber_mesh_y_log = new G4LogicalVolume(fiber_mesh_y, Air, "fiber mesh Y log", 0,0,0);

    fiber_mesh_z_phys = new G4PVPlacement(0,G4ThreeVector(),fiber_mesh_z_log, "fiber mesh z", aerogel_log, false, 0);
    fiber_mesh_y_phys = new G4PVPlacement(0,G4ThreeVector(),fiber_mesh_y_log, "fiber mesh y", aerogel_log, false, 0);
  }
  */

  //No overlays, no visualization.
  
  if (mesh.step_x !=0. && mesh.step_y != 0. && mesh.step_z !=0) {

    G4MultiUnion *fiber_mesh = new G4MultiUnion("fiber mesh");
    G4Tubs *fiber_z = new G4Tubs("fiber Z", 0., mesh.fiber_diam/2., Aer_width - 0.1*mm, 0.*deg, 360.*deg);

    //--    double xf = -Aer_depth + 0.1*mm + mesh.fiber_diam/2.;
    //--    while (xf <= Aer_depth - 0.1*mm - mesh.fiber_diam/2.) {
    double xf = -Aer_depth + 0.5*cm;
    while (xf < Aer_depth) {
      
      double yf = -Aer_hight + 0.1*mm + mesh.fiber_diam/2.;
      while (yf <= Aer_hight - 0.1*mm - mesh.fiber_diam/2.) {
      	G4RotationMatrix rotm;
      	G4ThreeVector trans_f(xf, yf, 0.);
      	G4Transform3D tr_f = G4Transform3D(rotm,trans_f);
      	fiber_mesh->AddNode(*fiber_z,tr_f);
      	yf += mesh.step_y;
      }
      
      xf += mesh.step_x;
    }

    G4Tubs *fiber_y = new G4Tubs("fiber Y", 0., mesh.fiber_diam/2., Aer_hight - 0.1*mm, 0.*deg, 360.*deg);

    xf = -Aer_depth + 0.1*mm + mesh.fiber_diam/2.;
    while (xf <= Aer_depth - 0.1*mm - mesh.fiber_diam/2.) {

      double zf = -Aer_width + 0.1*mm + mesh.fiber_diam/2.;
      while (zf <= Aer_width - 0.1*mm - mesh.fiber_diam/2.) {
	G4RotationMatrix rotm;
	rotm.rotateX(90.*deg);
	G4ThreeVector trans_f(xf, 0., zf);
	G4Transform3D tr_f = G4Transform3D(rotm,trans_f);
	fiber_mesh->AddNode(*fiber_y,tr_f);
	zf += mesh.step_z;
      }

      xf += mesh.step_x;
    }

    fiber_mesh->Voxelize();

    //    G4LogicalVolume *fiber_mesh_log = new G4LogicalVolume(fiber_mesh, Air, "fiber mesh log", 0,0,0);
    //    G4LogicalVolume *fiber_mesh_log = new G4LogicalVolume(fiber_mesh, CarbonMaterial, "fiber mesh log", 0,0,0);
    G4LogicalVolume *fiber_mesh_log = new G4LogicalVolume(fiber_mesh, CopperMaterial, "fiber mesh log", 0,0,0);

    G4OpticalSurface* FiberOptSurf = new G4OpticalSurface("FiberOptSurf");
    FiberOptSurf -> SetType(dielectric_metal);
    FiberOptSurf -> SetFinish(polished);
    FiberOptSurf -> SetModel(glisur);
    //    FiberOptSurf-> SetMaterialPropertiesTable(MPT_Carbon); 
    FiberOptSurf-> SetMaterialPropertiesTable(MPT_Copper); 
    //  FiberOptSurf->DumpInfo();

    new G4LogicalSkinSurface("FiberOptSurf", fiber_mesh_log, FiberOptSurf);

    fiber_mesh_z_phys = new G4PVPlacement(0,G4ThreeVector(),fiber_mesh_log, "fiber mesh", aerogel_log, false, 0);
 }

  Aerogel_phys = new G4PVPlacement(0, G4ThreeVector(-Ax + Aer_depth, 0., 0.),
				   aerogel_log, "Aerogel phys", aerogel_box_log, false, 0);

  //---------Gore
  // G4double Gore_hight = 12*3.*2.54/2.*cm;
  // G4double Gore_width = 30.*2.54/2.*cm ;
  G4double Gore_hight = hight;
  G4double Gore_width = width;
  G4double Gore_depth = 1.5*mm;
  G4Box *Gore = new G4Box("gore box", Gore_depth, Gore_hight, Gore_width);
  G4LogicalVolume *gore_log = new G4LogicalVolume(Gore, Millipore, "gore log", 0,0,0);
  Gore_phys = new G4PVPlacement(0, G4ThreeVector(depth_wrap - wrap_thick - Gore_depth, 0., 0.),
  				gore_log , "Gore phys", expHall_log , false, 0);
 
  ////  G4Box *Gore_top_bot = new G4Box("gore box", Dx, Gore_depth, width);
  G4Box *Gore_top_bot = new G4Box("gore top/bot", Dx-Gore_depth, Gore_depth, Gore_width);
  G4LogicalVolume *gore_top_bot_log = new G4LogicalVolume(Gore_top_bot, Millipore, "gore top/bot log", 0,0,0);
  ////  Gore_top_phys = new G4PVPlacement(0, G4ThreeVector(Ax, hight - Gore_depth, 0.),
  Gore_top_phys = new G4PVPlacement(0, G4ThreeVector(Ax-Gore_depth, hight - Gore_depth, 0.),
  				    gore_top_bot_log , "Gore top phys",expHall_log , false, 0);

  ////  Gore_bot_phys = new G4PVPlacement(0, G4ThreeVector(Ax, -hight + Gore_depth, 0.),
  Gore_bot_phys = new G4PVPlacement(0, G4ThreeVector(Ax-Gore_depth, -hight + Gore_depth, 0.),
  				    gore_top_bot_log , "Gore bot phys",expHall_log , false, 0);
 
 
  //Gore side covers with holes for PMTs
  G4Box* gore_side_box = new G4Box("gore side box", Dx, hight, Gore_depth);
  G4Tubs* PMT_hole = new G4Tubs("sub",0., hole_r, Gore_depth, 0.*deg,360.*deg);
  
  G4ThreeVector trans00(Dx - hole_pos_x, PMT_pos_y0, 0.);
  G4SubtractionSolid *Gore_side_sub = new G4SubtractionSolid("Gore side sub", gore_side_box, PMT_hole, 0, trans00);
  
  for(G4int i  = 0; i < n_PMT - 1; i++){
    G4ThreeVector trans01(Dx - hole_pos_x, -hight + (y0 + d*(i+1)), 0.);
    Gore_side_sub = new G4SubtractionSolid("Gore side sub", Gore_side_sub, PMT_hole, 0, trans01);
  }

  G4LogicalVolume *Gore_side_log = new G4LogicalVolume(Gore_side_sub, Millipore, "Gore side cover", 0,0,0);  
  // Gore_side_phys_right = new G4PVPlacement(0,G4ThreeVector(Ax, 0.*cm, -width + Gore_depth),
  //					  Gore_side_log,"Diffusion Box Al Side",expHall_log,false,0);

  // Gore_side_phys_left = new G4PVPlacement(0,G4ThreeVector(Ax, 0.*cm, width - Gore_depth),
  //					  Gore_side_log,"Diffusion Box Al Side",expHall_log,false,0);


  // -------PMT placements
  G4double PMT_pos_z;
  G4RotationMatrix* yRot = new G4RotationMatrix;
  for(G4int k = 0; k < 2; k++){
    if (k == 0){   //left PMTs
      PMT_pos_z = -width_wrap + zTopCutGlass;
      yRot->rotateY(M_PI*rad);
    }
    else { //right PMTs
      PMT_pos_z = width_wrap -zTopCutGlass;
      yRot=0;
    }
    
    for(G4int i  = 0; i < n_PMT; i++){
      G4double PMT_pos_y = PMT_pos_y0 + d*i;
      PMT_glass_phys[i][k] = new G4PVPlacement(yRot,G4ThreeVector(PMT_pos_x, PMT_pos_y, PMT_pos_z),
					       PMT_glass_log,"PMT glass", expHall_log,false,0);   
    }
  }

   //	------------- Surfaces --------------

  const G4int num = 11;    //from Benot
  
  G4double Ephoton[num] = 
    {h*c/900.*eV, h*c/600.*eV, h*c/550.*eV,
     h*c/520.*eV, h*c/500.*eV, h*c/450.*eV, h*c/420.*eV,
     h*c/400.*eV, h*c/365.*eV, h*c/315.*eV, h*c/190.*eV };
  
  /*
  G4double Ephoton_gore[num] = // for 3mm Gore reflector
    {h*c/900.*eV, h*c/600.*eV, h*c/550.*eV,
     h*c/520.*eV, h*c/500.*eV, h*c/450.*eV, h*c/420.*eV,
     h*c/400.*eV, h*c/365.*eV, h*c/250.*eV, h*c/190.*eV };
  */
  
  /*
  G4double Ephoton_gore[6] = // for Gore reflector
    {h*c/900.*eV, h*c/750.*eV, h*c/500.*eV, 
     h*c/400.*eV,  h*c/250.*eV, h*c/190.*eV };
  */

  G4double Ephoton_gore[6] = // for 0.5mm Gore reflector
    {h*c/900.*eV, h*c/750.*eV, h*c/500.*eV, 
     h*c/350.*eV,  h*c/300.*eV, h*c/190.*eV };

  // Original
  G4double MilliporeRefl[num] =   //from Benot, 1 layer
    {0.973, 0.973, 0.969, 
     0.969, 0.969, 0.967, 0.961, 
     0.956, 0.943, 0.854, 0.500,};
  

  
  /*  
  G4double GoreRefl[num] =   //3mm Gore reflector
    {0.99, 0.99, 0.99, 
     0.99, 0.99, 0.99, 0.99, 
     0.99, 0.99, 0.99, 0.500,};
  */

  /*
 G4double GoreRefl[6] =   //1mm Gore reflector
   {0.9, 0.98, 0.99,  
     0.999, 0.999, 0.500,};
  */

 G4double GoreRefl[6] =   //0.5mm Gore reflector
   {0.9, 0.97, 0.98,  
     0.995, 0.99, 0.500,};

/*
G4double MilliporeRefl[num] =   //from Benot, 3 layers
    {0.980, 0.980, 0.981, 
     0.982, 0.981, 0.975, 0.968, 
     0.963, 0.944, 0.888, 0.888,};
*/
  G4double MilliporeSpecLobe[num] = 
    {0.0,0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 
     0.0, 0.0, 0.0, 0.0,};
  
  G4double MilliporeSpecSpike[num] = 
    {0.05, 0.05, 0.05,
     0.05, 0.05, 0.05, 0.05, 
     0.05, 0.125, 0.2, 0.2,};
  
  G4double MilliporeBackScat[num] = 
    {0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 
     0.0, 0.0, 0.0, 0.0,};
  

  // Millipore wrapping
  G4OpticalSurface* OpMilSurface = new G4OpticalSurface("MilliporeSurface");
  OpMilSurface->SetType(dielectric_dielectric);
  OpMilSurface->SetFinish(groundfrontpainted);
  OpMilSurface->SetModel(unified);
  // OpMilSurface->SetSigmaAlpha(1.);

  
  // Gore wrapping
  G4OpticalSurface* OpGoreSurface = new G4OpticalSurface("GoreSurface");
  OpGoreSurface->SetType(dielectric_dielectric);
  OpGoreSurface->SetFinish(groundfrontpainted);
  OpGoreSurface->SetModel(unified);
  // OpMilSurface->SetSigmaAlpha(1.);
  
  G4MaterialPropertiesTable *myST3 = new G4MaterialPropertiesTable();
  myST3->AddProperty("REFLECTIVITY", Ephoton, MilliporeRefl, num);
  myST3->AddProperty("SPECULARLOBECONSTANT", Ephoton, MilliporeSpecLobe, num);
  myST3->AddProperty("SPECULARSPIKECONSTANT", Ephoton, MilliporeSpecSpike, num);
  myST3->AddProperty("BACKSCATTERCONSTANT", Ephoton, MilliporeBackScat, num);
  
  OpMilSurface->SetMaterialPropertiesTable(myST3);
  
  MilliporeSurface = new G4LogicalSkinSurface("millipore surface",millipore_wrap_log, OpMilSurface);

  if(MilliporeSurface->GetLogicalVolume() == millipore_wrap_log) G4cout << "Equal" << G4endl;
  ((G4OpticalSurface*)
   (MilliporeSurface->GetSurface(millipore_wrap_log)->GetSurfaceProperty()))->DumpInfo();

  
  //Gore
  G4MaterialPropertiesTable *myST4 = new G4MaterialPropertiesTable();
  myST4->AddProperty("REFLECTIVITY", Ephoton_gore, GoreRefl, num);
  myST4->AddProperty("SPECULARLOBECONSTANT", Ephoton, MilliporeSpecLobe, num);
  myST4->AddProperty("SPECULARSPIKECONSTANT", Ephoton, MilliporeSpecSpike, num);
  myST4->AddProperty("BACKSCATTERCONSTANT", Ephoton, MilliporeBackScat, num);
  
  OpGoreSurface->SetMaterialPropertiesTable(myST4);
   
  GoreSurface = new G4LogicalSkinSurface("Gore surface",gore_log, OpGoreSurface);
  GoreSurface1 = new G4LogicalSkinSurface("Gore surface",Gore_side_log, OpGoreSurface);
  GoreSurface2 = new G4LogicalSkinSurface("Gore surface",gore_top_bot_log, OpGoreSurface);

  if(GoreSurface->GetLogicalVolume() == gore_log) G4cout << "Equal" << G4endl;
  ((G4OpticalSurface*)
   (GoreSurface->GetSurface(gore_log)->GetSurfaceProperty()))->DumpInfo();

  if(GoreSurface1->GetLogicalVolume() == Gore_side_log) G4cout << "Equal" << G4endl;
  ((G4OpticalSurface*)
   (GoreSurface1->GetSurface(Gore_side_log)->GetSurfaceProperty()))->DumpInfo();

 if(GoreSurface2->GetLogicalVolume() == gore_top_bot_log ) G4cout << "Equal" << G4endl;
  ((G4OpticalSurface*)
   (GoreSurface2->GetSurface(gore_top_bot_log)->GetSurfaceProperty()))->DumpInfo();

  

  //--------Photocathode
  G4double wl_cathode;
  ifstream Indata_qe;
 // Indata_qe.open("inputs/qe_xp4500.dat");
    Indata_qe.open(qe_file);

  if(!Indata_qe) { // file couldn't be opened
    G4cerr << "Error: Cannot read QEs" << G4endl;
    exit(1);
  }

  // Reading QEs
  for (G4int k = n_qe - 1; k > -1; k-- ) {
    Indata_qe >> wl_cathode >> qe[k];
    PhotEn[k] = h*c/wl_cathode*eV; //190 - 700 nm		
  }  
  /*
    for (G4int i = 0; i < n; i++){
    G4cout << "WL = " << PhotEn[i]*1.e6 << "eV  qe = " << qe[i] << G4endl; 
    }
  */

  for (G4int i = 0; i < n_qe; i++) {
    RIndex_cathode[i] = 0.;
  }

  G4OpticalSurface* Cathode_surface = new G4OpticalSurface("Cathode");
  
  Cathode_surface -> SetType(dielectric_metal);
  Cathode_surface -> SetFinish(polished);
  Cathode_surface -> SetModel(glisur);
  
  G4MaterialPropertiesTable* Cathode_MPT = new G4MaterialPropertiesTable();
  Cathode_MPT -> AddProperty("REFLECTIVITY",PhotEn,RIndex_cathode,n_qe);
  Cathode_MPT -> AddProperty("EFFICIENCY",PhotEn, qe, n_qe);
  
  Cathode_surface -> SetMaterialPropertiesTable(Cathode_MPT);
  new G4LogicalSkinSurface("Cathode", photocathode_log, Cathode_surface);
  
  // Vis attributes
  expHall_log->SetVisAttributes (G4VisAttributes::Invisible);
  
  return expHall_phys;
}
