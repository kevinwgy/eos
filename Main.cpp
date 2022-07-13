#include <time.h>
#include <VarFcnSG.h>
#include <VarFcnMG.h>
#include <VarFcnJWL.h>
#include <VarFcnANEOSEx1.h>
#include <VarFcnDummy.h>

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

  //! Initialize VarFcn (EOS, etc.) 
  std::vector<VarFcnBase *> vf;
  for(int i=0; i<iod.eqs.materials.dataMap.size(); i++)
    vf.push_back(NULL); //allocate space for the VarFcn pointers

  for(auto it = iod.eqs.materials.dataMap.begin(); it != iod.eqs.materials.dataMap.end(); it++) {
    int matid = it->first;
    if(matid < 0 || matid >= vf.size()) {
      print_error("*** Error: Detected error in the specification of material indices (id = %d).\n", matid);
      exit_mpi();
    }
    if(it->second->eos == MaterialModelData::STIFFENED_GAS) {
      vf[matid] = new VarFcnSG(*it->second);
      print("- Initialized vf[%d]: Stiffened Gas.\n", matid);
    } else if(it->second->eos == MaterialModelData::MIE_GRUNEISEN) {
      vf[matid] = new VarFcnMG(*it->second);
      print("- Initialized vf[%d]: Mie-Gruneisen.\n", matid);
    } else if(it->second->eos == MaterialModelData::JWL) {
      vf[matid] = new VarFcnJWL(*it->second);
      print("- Initialized vf[%d]: Jones-Wilkins-Lee (JWL).\n", matid);
    } else if(it->second->eos == MaterialModelData::ANEOS_BIRCH_MURNAGHAN_DEBYE) {
      vf[matid] = new VarFcnANEOSEx1(*it->second);
      print("- Initialized vf[%d]: ANEOS: Birch-Murnaghan-Debye.\n", matid);
    } else {
      print_error("*** Error: Unable to initialize variable functions (VarFcn) for the specified material model.\n");
      exit_mpi();
    }
  }


  //--------------------------------------------------
  //                  TEST SECTION 
  //--------------------------------------------------












  //--------------------------------------------------
  //               END OF TEST SECTION
  //--------------------------------------------------



  print("\n");
  print("\033[0;32m==========================================\033[0m\n");
  print("\033[0;32m           NORMAL TERMINATION             \033[0m\n");
  print("\033[0;32m==========================================\033[0m\n");
  print("Total Computation Time: %f sec.\n", ((double)(clock()-start_time))/CLOCKS_PER_SEC);
  print("\n");

  for(int i=0; i<vf.size(); i++)
    delete vf[i];

  return 0;
}

//--------------------------------------------------------------

