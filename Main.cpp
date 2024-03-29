/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <time.h>
#include <VarFcnSG.h>
#include <VarFcnNASG.h>
#include <VarFcnMG.h>
#include <VarFcnMGExt.h>
#include <VarFcnTillot.h>
#include <VarFcnJWL.h>
#include <VarFcnANEOSEx1.h>
#include <VarFcnDummy.h>
#include <EOSAnalyzer.h>

//#include <polylogarithm_function.h>
//#include<ordinary_differential_equations.h>

#include <chrono>
using namespace std::chrono;

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
    } else if(it->second->eos == MaterialModelData::EXTENDED_MIE_GRUNEISEN) {
      vf[matid] = new VarFcnMGExt(*it->second);
      print("- Initialized vf[%d]: Extended Mie-Gruneisen.\n", matid);
    } else if(it->second->eos == MaterialModelData::TILLOTSON) {
      vf[matid] = new VarFcnTillot(*it->second);
      print("- Initialized vf[%d]: Tillotson.\n", matid);
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
  //   Output EOS Tabulation (If Specified by User)
  //--------------------------------------------------
  if(iod.special_tools.type == SpecialToolsData::EOS_TABULATION) {
    EOSAnalyzer eos_analyzer(iod.special_tools.eos_tabulationMap, vf);
    eos_analyzer.GenerateAllEOSTables();
  }

  //--------------------------------------------------
  //                  TEST SECTION 
  //--------------------------------------------------

  VarFcnBase* vf0 = vf[0];
  VarFcnBase* vf1 = vf[1];

  double rho = 1.0e-3;
  double e   = 0.0;
  double p0 = vf0->GetPressure(rho,e);
  double p1 = vf1->GetPressure(rho,e);
  fprintf(stdout,"p0 = %e, p1 = %e.\n", p0, p1);
  double p   = 1.0e13;
  double e0 = vf0->GetInternalEnergyPerUnitMass(rho,p);
  double e1 = vf1->GetInternalEnergyPerUnitMass(rho,p);
  fprintf(stdout,"e0 = %e, e1 = %e.\n", e0, e1);

  double rho0 = vf0->GetDensity(p,e0);
  double rho1 = vf1->GetDensity(p,e1);
  fprintf(stdout,"rho0 = %e, rho1 = %e.\n", rho0, rho1);

  double dpdrho0 = vf0->GetDpdrho(rho,e);
  double dpdrho1 = vf1->GetDpdrho(rho,e);
  double drho = 1e-8*rho;
  double dpdrho2 = (vf1->GetPressure(rho+drho, e) - vf1->GetPressure(rho-drho, e))/(2.0*drho);
  fprintf(stdout,"dpdrho0 = %e, dpdrho1 = %e, dpdrho2 = %e.\n", dpdrho0, dpdrho1, dpdrho2);

  double Gamma0 = vf0->GetBigGamma(rho,e);
  double Gamma1 = vf1->GetBigGamma(rho,e);
  double de = 1.0;
  double Gamma2 = 1.0/rho*(vf1->GetPressure(rho,e+de)-vf1->GetPressure(rho,e-de))/(2.0*de);
  fprintf(stdout,"Gamma0 = %e, Gamma1 = %e, Gamma2 = %e.\n", Gamma0, Gamma1, Gamma2);
  
/*
  double h = 5e14;
  e0 = vf0->GetInternalEnergyPerUnitMassFromEnthalpy(rho,h);
  e1 = vf1->GetInternalEnergyPerUnitMassFromEnthalpy(rho,h);
  fprintf(stdout,"[from entropy] e0 = %e, e1 = %e.\n", e0, e1);

  double T0, T1;
  auto start = high_resolution_clock::now();
  for(int i=0; i<1000000; i++)
    T0 = vf0->GetTemperature(rho,e);
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  fprintf(stdout,"T0 = %.12e. [Time] %e us.\n", T0, (double)duration.count());

  start = high_resolution_clock::now();
  for(int i=0; i<1000000; i++)
    T1 = vf1->GetTemperature(rho,e);
  stop = high_resolution_clock::now();
  duration = duration_cast<microseconds>(stop - start);
  fprintf(stdout,"T1 = %.12e. [Time] %e us.\n", T1, (double)duration.count());

  e0 = vf0->GetInternalEnergyPerUnitMassFromTemperature(rho,T0);
  e1 = vf1->GetInternalEnergyPerUnitMassFromTemperature(rho,T1);
  fprintf(stdout,"[from T] e0 = %e, e1 = %e.\n", e0, e1);
*/


/*
  double rho = 3.365839e-4;
  double p = -7.082216e6;
  double e = vf0->GetInternalEnergyPerUnitMass(rho, p);
  fprintf(stdout,"e = %e.\n", e);


  double c2 = vf0->ComputeSoundSpeedSquare(rho, e);
  fprintf(stdout,"c2 = %e.\n", c2);
*/

/*
  double rho = 0.998e-3;
  double e = 0.0; //7.0e12;
  fprintf(stdout,"------------------------\n");
  fprintf(stdout,"Input: rho = %e, e = %e.\n", rho, e); 
  fprintf(stdout,"------------------------\n");
  double p = vf0->GetPressure(rho,e);
  fprintf(stdout,"p = %e. (calculated)\n", p);
  double e1 = vf0->GetInternalEnergyPerUnitMass(rho,p);
  fprintf(stdout,"e = %e. (calculated)\n", e1);
  double e0 = vf0->GetReferenceInternalEnergyPerUnitMass();
  fprintf(stdout,"e0 = %e.\n", e0);
  double rho2 = vf0->GetDensity(p,e);
  fprintf(stdout,"rho = %e. (calculated)\n", rho2);
  double dpdrho = vf0->GetDpdrho(rho,e);
  fprintf(stdout,"dpdrho = %e.\n", dpdrho);
  double Gamma = vf0->GetBigGamma(rho,e);
  fprintf(stdout,"Gamma = %e.\n", Gamma);
  double T = vf0->GetTemperature(rho,e);
  fprintf(stdout,"T = %e.\n", T);
  double T0 = vf0->GetReferenceTemperature();
  fprintf(stdout,"T0 = %e.\n", T0);
  double e2 = vf0->GetInternalEnergyPerUnitMassFromTemperature(rho,T);
  fprintf(stdout,"e = %e. (calculated from T)\n", e2);
  double h = e + p/rho;
  fprintf(stdout,"h = %e.\n", h);
  double e3 = vf0->GetInternalEnergyPerUnitMassFromEnthalpy(rho,h);
  fprintf(stdout,"e = %e. (calculated from h)\n", e3);

*/

/*

  auto fun = [&]([[maybe_unused]] double *U, [[maybe_unused]] double t, double *dUdt) {
    dUdt[0] = U[0];  
    dUdt[1] = exp(t);
    return;
  };

  auto state_checker = [&]([[maybe_unused]] double *U) {
    return false;
  };

  lalala test;

  double U0[2] = {1.0, 1.0}, U[2], last_dt(1000.0);
  std::vector<double> Uall;
  std::vector<double> tall;
  //MathTools::runge_kutta_4(fun, 2, 0.0, U0, 0.1, 1.0, U, state_checker, &Uall, &last_dt);
  MathTools::runge_kutta_45(fun, 2, 0.0, U0, 0.1, 1.0, U, 
                            state_checker,
                            1e-9, (int)1e7, 
                            &tall, &Uall);

  fprintf(stderr, "Sol = %e %e. last_dt = %e.\n", U[0], U[1], last_dt);
  //for(int i=0; i<(int)(Uall.size()/2); i++)
  //  fprintf(stderr,"Uall[%d] = %e %e.\n", i, Uall[2*i], Uall[2*i+1]);
  for(int i=0; i<(int)tall.size(); i++)
    fprintf(stderr,"tall[%d] = %e, Uall = %e %e.\n", i, tall[i], Uall[2*i], Uall[2*i+1]);

*/

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

