#include "PHG4BarrelEcalDetector.h"
#include "PHG4BarrelEcalDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <phool/recoConsts.h>
#include <Geant4/G4Trap.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4ThreeVector.hh>      // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>      // for G4Transform3D
#include <Geant4/G4Types.hh>            // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume
#include <Geant4/G4DisplacedSolid.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4MaterialPropertiesTable.hh>  // for G4MaterialProperties...
#include <Geant4/G4MaterialPropertyVector.hh>   // for G4MaterialPropertyVector

#include <TSystem.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>  // for pair, make_pair

using namespace std;

class G4VSolid;
class PHCompositeNode;

//_______________________________________________________________________
PHG4BarrelEcalDetector::PHG4BarrelEcalDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4BarrelEcalDisplayAction*>(subsys->GetDisplayAction()))
  , m_Params(parameters)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_AbsorberActiveFlag(m_Params->get_int_param("absorberactive"))
  , m_SupportActiveFlag(m_Params->get_int_param("supportactive"))
  , m_TowerLogicNamePrefix("bcalTower")
  , m_SuperDetector("NONE")
{
}
//_______________________________________________________________________
int PHG4BarrelEcalDetector::IsInBarrelEcal(G4VPhysicalVolume* volume) const
{
  G4LogicalVolume* mylogvol = volume->GetLogicalVolume();
  if (m_ActiveFlag)
  {
    if (m_ScintiLogicalVolSet.find(mylogvol) != m_ScintiLogicalVolSet.end())
    {
      return 1;
    }
  }

  if (m_AbsorberActiveFlag)
  {
    if (m_AbsorberLogicalVolSet.find(mylogvol) != m_AbsorberLogicalVolSet.end())
    {
      return -1;
    }
  }

  if (m_SupportActiveFlag)
  {
    if (m_SupportLogicalVolSet.find(mylogvol) != m_SupportLogicalVolSet.end())
    {
      return -2;
    }
  }
  return 0;
}

//_______________________________________________________________________
void PHG4BarrelEcalDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  if (Verbosity() > 0)
  {
    std::cout << "PHG4BarrelEcalDetector: Begin Construction" << std::endl;
  }

  if (m_Params->get_string_param("mapping_file").empty())
  {
    std::cout << "ERROR in PHG4BarrelEcalDetector: No mapping file specified. Abort detector construction." << std::endl;
    std::cout << "Please run set_string_param(\"mapping_file\", std::string filename ) first." << std::endl;
    gSystem->Exit(1);
  }

  ParseParametersFromTable();

 const double Length = becal_length;
 const double max_radius = Radius + tower_length + elec_length + support_length;
 const double pos_x1 = 0*cm;
 const double pos_y1 = 0*cm;
 const double pos_z1 = 0*cm;

  G4Tubs *cylinder_solid = new G4Tubs("BCAL_SOLID",
                                       Radius, max_radius,
                                       Length/ 2.0, 0, 2*M_PI);


  G4Material *cylinder_mat = G4Material::GetMaterial("G4_AIR");
  assert(cylinder_mat);


  G4LogicalVolume *cylinder_logic = new G4LogicalVolume(cylinder_solid, cylinder_mat,
                                       "BCAL_SOLID", 0, 0, 0);

  m_DisplayAction->AddVolume(cylinder_logic, "BCalCylinder");

  //cylinder_physi = 

  std::string name_envelope = m_TowerLogicNamePrefix + "_envelope";

  new G4PVPlacement(0, G4ThreeVector(pos_x1, pos_y1, pos_z1), cylinder_logic, name_envelope,
                                     logicWorld, false, 0, OverlapCheck());


  /* Construct single calorimeter tower */
  G4LogicalVolume* singletower = ConstructTower();

  PlaceTower(cylinder_logic, singletower);

  return;
}


G4LogicalVolume*
PHG4BarrelEcalDetector::ConstructTower()
{

  double lower_side = 2*M_PI*Radius/nTowers_layer;
  double upper_side = 2*M_PI*(Radius+tower_length)/nTowers_layer;


  G4Trap* block_solid = new G4Trap(
                  "solid_tower",
      tower_length/2,                                                 // G4double pDz,
      0,  0,                                                    // G4double pTheta, G4double pPhi,
      lower_side/2, lower_side/2, lower_side/2,                 // G4double pDy1, G4double pDx1, G4double pDx2,
      0,                                                       // G4double pAlp1,
      upper_side/2, lower_side/2,lower_side/2,                                    // G4double pDy2, G4double pDx3, G4double pDx4,
      0                                                         // G4double pAlp2 //
  );

  
  G4double density;
  G4int ncomponents;
  G4Material *sciglass1 = new G4Material("sciglass1", density = 4.22 * g / cm3, ncomponents = 4,  kStateSolid);
  sciglass1->AddElement(G4Element::GetElement("Ba"), 0.3875);
  sciglass1->AddElement(G4Element::GetElement("Gd"), 0.2146);
  sciglass1->AddElement(G4Element::GetElement("Si"), 0.1369);
  sciglass1->AddElement(G4Element::GetElement("O"),  0.2610);

  G4Material* cylinder_mat = G4Material::GetMaterial("sciglass1");
  assert(cylinder_mat);

  G4LogicalVolume* block_logic = new G4LogicalVolume(block_solid, cylinder_mat,
                                                     "solid_tower", 0, 0,
                                                     nullptr);
  m_ScintiLogicalVolSet.insert(block_logic);
  m_DisplayAction->AddVolume(block_logic, "Block");

  return block_logic;
}

int PHG4BarrelEcalDetector::PlaceTower(G4LogicalVolume* sec, G4LogicalVolume* singletower)
{
  /* Loop over all tower positions in vector and place tower */
  for (std::map<std::string, towerposition>::iterator iterator = m_TowerPostionMap.begin(); iterator != m_TowerPostionMap.end(); ++iterator)
  {
    if (Verbosity() > 0)
    {
      std::cout << "PHG4BarrelEcalDetector: Place tower " << iterator->first
                << " idx_j = " << iterator->second.idx_j << ", idx_k = " << iterator->second.idx_k<< std::endl;
                //<< " at x = " << iterator->second.x << " , y = " << iterator->second.y << " , z = " << iterator->second.z << std::endl;
    }
     int copyno = (iterator->second.idx_j << 16) + iterator->second.idx_k;


     G4RotationMatrix becal_rotm;
     becal_rotm.rotateY(iterator->second.roty);
     becal_rotm.rotateZ(iterator->second.rotz);

      new G4PVPlacement(G4Transform3D(becal_rotm, G4ThreeVector(iterator->second.centerx, iterator->second.centery, iterator->second.centerz)),
                      singletower,
                      iterator->first,
                      sec,
                      0, copyno, OverlapCheck());


    
  }
  return 0;
}


int PHG4BarrelEcalDetector::ParseParametersFromTable()
{
  /* Open the datafile, if it won't open return an error */
  std::ifstream istream_mapping;
  istream_mapping.open(m_Params->get_string_param("mapping_file"));
  if (!istream_mapping.is_open())
  {
    std::cout << "ERROR in PHG4BarrelEcalDetector: Failed to open mapping file " << m_Params->get_string_param("mapping_file") << std::endl;
    gSystem->Exit(1);
  }

  /* loop over lines in file */
  std::string line_mapping;
  while (getline(istream_mapping, line_mapping))
  {

    std::istringstream iss(line_mapping);

    unsigned idphi_j, ideta_k;
    G4double cx, cy, cz;
    G4double rot_z, rot_y;
    std::string dummys;
 
    if (!(iss >> dummys >> idphi_j >> ideta_k >> cx >> cy >> cz >> rot_y >> rot_z))
    {
      std::cout << "ERROR in PHG4BarrelEcalDetector: Failed to read line in mapping file " << m_Params->get_string_param("mapping_file") << std::endl;
      gSystem->Exit(1);
    }

    /* Construct unique name for tower */
      /* Mapping file uses cm, this class uses mm for length */
      std::ostringstream towername;
      towername.str("");
      towername << m_TowerLogicNamePrefix << "_j_" << idphi_j << "_k_" << ideta_k;

      /* insert tower into tower map */
      towerposition tower_new;
      tower_new.centerx = cx;
      tower_new.centery = cy;
      tower_new.centerz = cz;
      tower_new.roty = rot_y;
      tower_new.rotz = rot_z;
      tower_new.idx_j = idphi_j;
      tower_new.idx_k = ideta_k;
      m_TowerPostionMap.insert(make_pair(towername.str(), tower_new));

      
    }

  return 0;
}
