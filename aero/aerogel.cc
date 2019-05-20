#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "G4ios.hh"

#include "aerogelDetectorConstruction.hh"
//#include "aerogelPhysicsList.hh"
#include "aerogelOpticalPhysics.hh"
//#include "LHEP.hh"
////#include "QGSP.hh"
#include "QGSP_BERT.hh"
#include "aerogelPrimaryGeneratorAction.hh"
#include "aerogelRunAction.hh"
#include "aerogelEventAction.hh"
#include "aerogelStackingAction.hh"
#include "aerogelSteppingVerbose.hh"
#include "aerogelSteppingAction.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <error.h>
#include <errno.h>

using namespace std;

int getSeed();

int getSeed()
{
    // open /dev/random
    int fd;
    if ((fd = open("/dev/random", O_RDONLY)) == -1)
    {
        error(1, errno, "Fatal: unable to open /dev/random for reading.");
    }
    // read it
    int32_t value;
    ssize_t len;

    switch (len = read(fd, &value, sizeof(int32_t)))
    {
        case sizeof(int):
            printf("The seed is %i\n", value);
            break;
        case -1:
            error(1, errno, "Fatal: unable to read from /dev/random.");
            break;
        default:
            error(1, 0, "Fatal: %li bytes was read from " 
                  "/dev/random instead of %li.", len, sizeof(int32_t));
    }
    // close file
    if (close(fd) == -1)
    {
        error(0, errno, "Warning: unable to close /dev/random properly.");
    }
    return value;
}


//----Batch mode

/*
int main()
{
   
  // Seed the random number generator manually
  //
  G4long myseed = getSeed();
  CLHEP::HepRandom::setTheSeed(myseed);
  
   
  // Run manager
  //
  G4RunManager* runManager = new G4RunManager;

  // Initialization classes - mandatory
  //
  aerogelDetectorConstruction *detector = new aerogelDetectorConstruction;
  runManager->SetUserInitialization(detector);
  // runManager->SetUserInitialization(new aerogelPhysicsList);
  
  //  LHEP* theLHEP = new LHEP;
  //  theLHEP->RegisterPhysics(new aerogelOpticalPhysics("optical")); 
  //  runManager->SetUserInitialization(theLHEP);	
  


 // QGSP* theQGSP = new QGSP;
 // theQGSP->RegisterPhysics(new aerogelOpticalPhysics("optical")); 
 // runManager->SetUserInitialization(theQGSP);



  QGSP_BERT* theQGSP_BERT = new QGSP_BERT;
  theQGSP_BERT->RegisterPhysics(new aerogelOpticalPhysics("optical")); 
  runManager->SetUserInitialization(theQGSP_BERT);


  // Action classes
  //
  aerogelRunAction* runaction = new aerogelRunAction;
  runManager->SetUserAction(runaction);
  //
  aerogelEventAction* eventaction = new aerogelEventAction(runaction); 
  runManager->SetUserAction(eventaction);
  //
 // G4VUserPrimaryGeneratorAction* gen_action = new aerogelPrimaryGeneratorAction;
  aerogelPrimaryGeneratorAction* gen_action = new aerogelPrimaryGeneratorAction(eventaction);
  runManager->SetUserAction(gen_action);

  aerogelSteppingAction* steppingaction = new aerogelSteppingAction(detector, eventaction);
  runManager->SetUserAction(steppingaction);  

  G4UserStackingAction* stacking_action = new aerogelStackingAction;
  runManager->SetUserAction(stacking_action);
  
  // Initialize G4 kernel
  //
  runManager->Initialize();
    
  //  int numberOfEvent = 400;
  int numberOfEvent = 100;

  runManager->BeamOn(numberOfEvent);

    delete runManager;

  return 0;
}

*/

   
//---Interactive mode

int main(int argc,char** argv)
{
  // are user asking for some help?
  if (argc >= 2 &&
      (!strncmp(argv[1], "-h", 2) ||
       !strncmp(argv[1], "--help", 6)) )
    {
      printf("This program returns veritably random 32-bit signed integer " 
	     "using /dev/random source.\n");
      exit(0);
    }
  // Seed the random number generator manually
  //
  G4long myseed = getSeed();
  CLHEP::HepRandom::setTheSeed(myseed);

  cout << "==> Aerogel: argc = " << argc << endl;
  for (int i=0; i<argc; i++) {
    cout << "==>  argv " << i << " = " << argv[i] << endl;
  }

#ifdef G4UI_USE

  // Detect interactive mode (if no arguments) and define UI session

  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
    cout << "==>  UI session defined." << endl;
  }

#endif
    
  // Run manager
  //
  G4RunManager* runManager = new G4RunManager;

  // Initialization classes - mandatory
  //
  aerogelDetectorConstruction *detector = new aerogelDetectorConstruction;
  runManager->SetUserInitialization(detector);
  // runManager->SetUserInitialization(new HadronPhysicsLHEP);
  //runManager->SetUserInitialization(new aerogelPhysicsList);
  cout << "==>  RunManager initialized with detector." << endl;
  

//  LHEP* theLHEP = new LHEP;
//  theLHEP->RegisterPhysics(new aerogelOpticalPhysics("optical")); 
//  runManager->SetUserInitialization(theLHEP);	
  
////  QGSP* theQGSP = new QGSP;
////  theQGSP->RegisterPhysics(new aerogelOpticalPhysics("optical")); 
////  runManager->SetUserInitialization(theQGSP);

  QGSP_BERT* theQGSP_BERT = new QGSP_BERT;
  theQGSP_BERT->RegisterPhysics(new aerogelOpticalPhysics("optical")); 
  runManager->SetUserInitialization(theQGSP_BERT);
  cout << "==>  RunManager initialized with QGSP_BERT." << endl;

  // Action classes
  //
  aerogelRunAction* runaction = new aerogelRunAction;
  runManager->SetUserAction(runaction);
  //
  aerogelEventAction* eventaction = new aerogelEventAction(runaction); 
  runManager->SetUserAction(eventaction);
  //
  aerogelPrimaryGeneratorAction* gen_action = new aerogelPrimaryGeneratorAction(eventaction);
  runManager->SetUserAction(gen_action);

  aerogelSteppingAction* steppingaction = new aerogelSteppingAction(detector, eventaction);
  runManager->SetUserAction(steppingaction);  

  G4UserStackingAction* stacking_action = new aerogelStackingAction;
  runManager->SetUserAction(stacking_action);

  cout << "==>  Runmanager actions set." << endl;
  
  // Initialize G4 kernel
  //
  runManager->Initialize();
  cout << "==>  RunManager initialized." << endl;
    
  // Get the pointer to the User Interface manager
  //
  /*     

   G4UImanager* UI = G4UImanager::GetUIpointer(); 

   if (argc==1)   // Define UI session for interactive mode
     {
       G4UIsession* session = 0;
#ifdef G4UI_USE_TCSH
       session = new G4UIterminal(new G4UItcsh);      
#else
       session = new G4UIterminal();
#endif    
       UI->ApplyCommand("/control/execute vis.mac"); 
       session->SessionStart();
       delete session;
     }
   
   else         // Batch mode
     {
       G4String command = "/control/execute ";
       G4String fileName = argv[1];
       UI->ApplyCommand(command+fileName);
     }
   */

#ifdef G4VIS_USE
  // visualization manager
  //
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
  G4cout << "==> visManager initialized" << G4endl;
#endif
  
   G4UImanager* UImanager = G4UImanager::GetUIpointer(); 
   cout << "==>  Got UImanager." << endl;

   if ( argc != 1 ) {

     // execute an argument macro file if exist
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UImanager->ApplyCommand(command+fileName);

     cout << "==>  macro file executed." << endl;
   }
   else {

     cout << "==>  argc==1 branch:" << endl;

#ifdef G4UI_USE

#ifdef G4VIS_USE
     cout << "==>  init_vis.mac to be executed." << endl;
     UImanager->ApplyCommand("/control/execute init_vis.mac"); 
     cout << "==>  init_vis.mac executed." << endl;
#else
     UImanager->ApplyCommand("/control/execute init.mac"); 
     cout << "==>  init.mac executed." << endl;
#endif

     if (ui->IsGUI()) {
       ///         UImanager->ApplyCommand("/control/execute gui.mac");
       cout << "==>  ui is GUI." << endl;
     }     

     // start interactive session
     ui->SessionStart();
     cout << "==>  UI session started." << endl;
     delete ui;

#endif

   }   //argc!=1
  
#ifdef G4VIS_USE
   delete visManager;
#endif

   delete runManager;
   // delete verbosity;

   return 0;
}
