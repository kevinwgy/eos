#include <time.h>
#include <VarFcnSG.h>
#include <VarFcnNASG.h>
#include <VarFcnMG.h>
#include <VarFcnJWL.h>
#include <VarFcnANEOSEx1.h>
#include <VarFcnDummy.h>

#include <polylogarithm_function.h>

using std::cout;
using std::endl;

int verbose = 0;

/*************************************
 * Main Function
 ************************************/
int main(int argc, char* argv[])
{

  clock_t start_time = clock(); //for timing purpose only

  print("\033[0;32m==========================================\033[0m\n");
  print("\033[0;32m                 START                    \033[0m\n");
  print("\033[0;32m==========================================\033[0m\n");
  print("\n");

  //! Read user's input file (read the parameters)
  IoData iod(argc, argv);
  verbose = iod.output.verbose;

  //! Finalize IoData (read additional files and check for errors)
  iod.finalize();

  VarFcnANEOSEx1 *vf0 = NULL;

  //! Initialize VarFcn (EOS, etc.) 
  std::vector<VarFcnBase *> vf;
  for(int i=0; i<(int)iod.eqs.materials.dataMap.size(); i++)
    vf.push_back(NULL); //allocate space for the VarFcn pointers

  for(auto it = iod.eqs.materials.dataMap.begin(); it != iod.eqs.materials.dataMap.end(); it++) {
    int matid = it->first;
    if(matid < 0 || matid >= (int)vf.size()) {
      print_error("*** Error: Detected error in the specification of material indices (id = %d).\n", matid);
      exit_mpi();
    }
    if(it->second->eos == MaterialModelData::STIFFENED_GAS) {
      vf[matid] = new VarFcnSG(*it->second);
      print("- Initialized vf[%d]: Stiffened Gas.\n", matid);
    } else if(it->second->eos == MaterialModelData::NOBLE_ABEL_STIFFENED_GAS) {
      vf[matid] = new VarFcnNASG(*it->second);
      print("- Initialized vf[%d]: Noble-Abel Stiffened Gas.\n", matid);
    } else if(it->second->eos == MaterialModelData::MIE_GRUNEISEN) {
      vf[matid] = new VarFcnMG(*it->second);
      print("- Initialized vf[%d]: Mie-Gruneisen.\n", matid);
    } else if(it->second->eos == MaterialModelData::JWL) {
      vf[matid] = new VarFcnJWL(*it->second);
      print("- Initialized vf[%d]: Jones-Wilkins-Lee (JWL).\n", matid);
    } else if(it->second->eos == MaterialModelData::ANEOS_BIRCH_MURNAGHAN_DEBYE) {
      vf[matid] = new VarFcnANEOSEx1(*it->second);
      vf0 = new VarFcnANEOSEx1(*it->second);
      print("- Initialized vf[%d]: ANEOS: Birch-Murnaghan-Debye.\n", matid);
    } else {
      print_error("*** Error: Unable to initialize variable functions (VarFcn) for the specified material model.\n");
      exit_mpi();
    }
  }


  //--------------------------------------------------
  //                  TEST SECTION 
  //--------------------------------------------------

  double rho = 0.00896;
  double p   = 1.0e5;
  double e, T;
/*
  fprintf(stderr,"\n\n");
  fprintf(stderr,"Mie-Gruneisen:\n");
  e = vf[1]->GetInternalEnergyPerUnitMass(rho, p);
  fprintf(stderr,"e = %e.\n", e);
  T = vf[1]->GetTemperature(rho, e);
  fprintf(stderr,"T = %e.\n", T);
*/

  fprintf(stderr,"ANEOS:\n");
  e = vf[2]->GetInternalEnergyPerUnitMass(rho, p);
  fprintf(stderr,"e = %e.\n", e);
  T = vf[2]->GetTemperature(rho, e);
  fprintf(stderr,"T = %e.\n", T);



  double e_hot0  = vf0->ComputeThermalSpecificEnergy(rho,T-1.0e-6*T);
  double e_hot1  = vf0->ComputeThermalSpecificEnergy(rho,T+1.0e-6*T);
  double cv = (e_hot1 - e_hot0)/(2.0e-6*T);
  fprintf(stderr,"cv = %e.\n", cv);


  fprintf(stderr,"Speed of sound = %e.\n", vf0->ComputeSoundSpeed(rho,e));

  T = vf[2]->GetTemperature(rho, e);
  fprintf(stderr,"T = %e.\n", T);
  //--------------------------------------------------
  //               END OF TEST SECTION
  //--------------------------------------------------



  print("\n");
  print("\033[0;32m==========================================\033[0m\n");
  print("\033[0;32m           NORMAL TERMINATION             \033[0m\n");
  print("\033[0;32m==========================================\033[0m\n");
  print("Total Computation Time: %f sec.\n", ((double)(clock()-start_time))/CLOCKS_PER_SEC);
  print("\n");

  for(int i=0; i<(int)vf.size(); i++)
    delete vf[i];

  return 0;
}

//--------------------------------------------------------------

