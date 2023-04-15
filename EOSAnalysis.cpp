#include<EOSAnalysis.h>
#include<fstream>
#include<cstring> //strcmp
#include<cassert>
using std::vector;

//------------------------------------------------------------------------------

EOSAnalysis::EOSAnalysis(ObjectMap<EOSTabulationData> &eos_tabulationMap_,
                         vector<VarFcnBase*>& vf_)
           : eos_tabulationMap(eos_tabulationMap_), vf(vf_) 
{ }

//------------------------------------------------------------------------------

EOSAnalysis::~EOSAnalysis()
{ }

//------------------------------------------------------------------------------

void
EOSAnalysis::GenerateAllEOSTables()
{
  for(auto it = eos_tabulationMap.dataMap.begin(); 
           it != eos_tabulationMap.dataMap.end(); it++)
    GenerateEOSTable(*it->second);
}

//------------------------------------------------------------------------------

void
EOSAnalysis::GenerateEOSTable(EOSTabulationData& tab)
{

  int id = tab.materialid;
  if(id<0 || id>=(int)vf.size()) {
    print_error("*** Error: Cannot generate EOS table for undefined "
                "material id: %d.\n", id);
    exit(-1);
  }

  if(!strcmp(tab.filename,"")) {
    print_error("*** Error: Cannot generate EOS table for"
                   "material (%d). Output file is not specified.\n", id);
    exit(-1);
  }


  FILE *file = fopen(tab.filename, "w");
  if(!file) {
    print_error("*** Error: Cannot write file %s.\n", tab.filename);
    exit(-1);
  } else
    print("- Outputing EOS Tabulation (material id: %d) to %s.\n", id, tab.filename);



  vector<vector<double> > result;
  double x0 = tab.x0, y0 = tab.y0, xmax = tab.xmax, ymax = tab.ymax;
  int Nx = tab.Nx, Ny = tab.Ny;

  if(tab.output == EOSTabulationData::PRESSURE) {

    if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::SPECIFIC_INTERNAL_ENERGY) {
      print(file, "# Tabulating Pressure as a function of ");
      print(file, "Density [%e, %e] and Specific Internal Energy [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double rho, double e) {return vf[id]->GetPressure(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::TEMPERATURE) {
      print(file, "# Tabulating Pressure as a function of ");
      print(file, "Density [%e, %e] and Temperature [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double rho, double T) {
        double e = vf[id]->GetInternalEnergyPerUnitMassFromTemperature(rho,T); 
        return vf[id]->GetPressure(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }
  }

  else if(tab.output == EOSTabulationData::SPECIFIC_INTERNAL_ENERGY) {

    if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::PRESSURE) {
      print(file, "# Tabulating Specific Internal Energy as a function of ");
      print(file, "Density [%e, %e] and Pressure [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double rho, double p) {return vf[id]->GetInternalEnergyPerUnitMass(rho,p);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::TEMPERATURE) {
      print(file, "# Tabulating Specific Internal Energy as a function of ");
      print(file, "Density [%e, %e] and Temperature [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double rho, double T) {return vf[id]->GetInternalEnergyPerUnitMassFromTemperature(rho,T);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }
  }

  else if(tab.output == EOSTabulationData::DENSITY) {
    if(tab.xvar == EOSTabulationData::PRESSURE && tab.yvar == EOSTabulationData::SPECIFIC_INTERNAL_ENERGY) {
      print(file, "# Tabulating Density as a function of ");
      print(file, "Pressure [%e, %e] and Specific Internal Energy [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double p, double e) {return vf[id]->GetDensity(p,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }
  }

  else if(tab.output == EOSTabulationData::DP_DE) {

    if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::SPECIFIC_INTERNAL_ENERGY) {
      print(file, "# Tabulating PressureDerivativeEnergy as a function of ");
      print(file, "Density [%e, %e] and Specific Internal Energy [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double rho, double e) {return rho*vf[id]->GetBigGamma(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::PRESSURE) {
      print(file, "# Tabulating PressureDerivativeEnergy as a function of ");
      print(file, "Density [%e, %e] and Pressure [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double rho, double p) {
        double e = vf[id]->GetInternalEnergyPerUnitMass(rho,p);
        return rho*vf[id]->GetBigGamma(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::TEMPERATURE) {
      print(file, "# Tabulating PressureDerivativeEnergy as a function of ");
      print(file, "Density [%e, %e] and Temperature [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double rho, double T) {
        double e = vf[id]->GetInternalEnergyPerUnitMassFromTemperature(rho,T);
        return rho*vf[id]->GetBigGamma(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }
  }

  else if(tab.output == EOSTabulationData::GRUNEISEN_PARAMETER) {

    if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::SPECIFIC_INTERNAL_ENERGY) {
      print(file, "# Tabulating Gruneisen Parameter as a function of ");
      print(file, "Density [%e, %e] and Specific Internal Energy [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double rho, double e) {return vf[id]->GetBigGamma(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::PRESSURE) {
      print(file, "# Tabulating Gruneisen Parameter as a function of ");
      print(file, "Density [%e, %e] and Pressure [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double rho, double p) {
        double e = vf[id]->GetInternalEnergyPerUnitMass(rho,p);
        return vf[id]->GetBigGamma(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }
    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::TEMPERATURE) {
      print(file, "# Tabulating Gruneisen Parameter as a function of ");
      print(file, "Density [%e, %e] and Temperature [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double rho, double T) {
        double e = vf[id]->GetInternalEnergyPerUnitMassFromTemperature(rho,T);
        return vf[id]->GetBigGamma(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }
  }

  else if(tab.output == EOSTabulationData::DP_DRHO) {

    if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::SPECIFIC_INTERNAL_ENERGY) {
      print(file, "# Tabulating PressureDerivativeDensity as a function of ");
      print(file, "Density [%e, %e] and Specific Internal Energy [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double rho, double e) {return vf[id]->GetDpdrho(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::PRESSURE) {
      print(file, "# Tabulating PressureDerivativeDensity as a function of ");
      print(file, "Density [%e, %e] and Pressure [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double rho, double p) {
        double e = vf[id]->GetInternalEnergyPerUnitMass(rho,p);
        return vf[id]->GetDpdrho(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::TEMPERATURE) {
      print(file, "# Tabulating PressureDerivativeDensity as a function of ");
      print(file, "Density [%e, %e] and Temperature [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double rho, double T) {
        double e = vf[id]->GetInternalEnergyPerUnitMassFromTemperature(rho,T);
        return vf[id]->GetDpdrho(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }
  }

  else if(tab.output == EOSTabulationData::BULK_MODULUS) {

    if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::SPECIFIC_INTERNAL_ENERGY) {
      print(file, "# Tabulating Bulk Modulus as a function of ");
      print(file, "Density [%e, %e] and Specific Internal Energy [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double rho, double e) {return rho*vf[id]->GetDpdrho(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::PRESSURE) {
      print(file, "# Tabulating Bulk Modulus as a function of ");
      print(file, "Density [%e, %e] and Pressure [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double rho, double p) {
        double e = vf[id]->GetInternalEnergyPerUnitMass(rho,p);
        return rho*vf[id]->GetDpdrho(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::TEMPERATURE) {
      print(file, "# Tabulating Bulk Modulus as a function of ");
      print(file, "Density [%e, %e] and Temperature [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double rho, double T) {
        double e = vf[id]->GetInternalEnergyPerUnitMassFromTemperature(rho,T);
        return rho*vf[id]->GetDpdrho(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }
  }

  else if(tab.output == EOSTabulationData::TEMPERATURE) {

    if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::SPECIFIC_INTERNAL_ENERGY) {
      print(file, "# Tabulating Temperature as a function of ");
      print(file, "Density [%e, %e] and Specific Internal Energy [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double rho, double e) {return vf[id]->GetTemperature(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::DENSITY && tab.yvar == EOSTabulationData::PRESSURE) {
      print(file, "# Tabulating Temperature as a function of ");
      print(file, "Density [%e, %e] and Pressure [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double rho, double p) {
        double e = vf[id]->GetInternalEnergyPerUnitMass(rho,p);
        return vf[id]->GetTemperature(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }

    else if(tab.xvar == EOSTabulationData::PRESSURE && tab.yvar == EOSTabulationData::SPECIFIC_INTERNAL_ENERGY) {
      print(file, "# Tabulating Temperature as a function of ");
      print(file, "Pressure [%e, %e] and Specific Internal Energy [%e, %e]\n", x0, xmax, y0, ymax);
      auto fun = [&](double p, double e) {
        double rho = vf[id]->GetDensity(p,e);
        return vf[id]->GetTemperature(rho,e);};
      tabulate2Dfunction_uniform(fun, x0, xmax, Nx, y0, ymax, Ny, result);
    }
  }

  else {
    print_error("*** Error: Unable to create EOS table for the two independent variables specified by user. "
                "(Try switching the two variables.)\n");
    fclose(file);
    delete file;
    exit(-1);
  }

  assert(result.size()>0);


  // Write data to file

  if(result.size()==1) { //1D 
    double vmin, dv;
    if(x0==xmax || Nx==1) {
      print(file,"# First variable fixed at %e.\n", x0);
      vmin = y0;
      dv = Ny==1 ? 0.0 : (ymax-y0)/(Ny-1); 
    } else {
      assert(y0==ymax || Ny==1);
      print(file,"# Second variable fixed at %e.\n", y0);
      vmin = x0;
      dv = Nx==1 ? 0.0 : (xmax-x0)/(Nx-1); 
    }
    for(int i=0; i<(int)result[0].size(); i++)
      print(file, "%16.8e  %16.8e\n", vmin+i*dv, result[0][i]);

    fclose(file);
    delete file;

    return;
  }

  
  //2D
  assert(result.size()>1);
  print(file,"# Unstructured data\n");
  double dx = Nx<=1 ? 0.0 : (xmax-x0)/(Nx-1);
  double dy = Ny<=1 ? 0.0 : (ymax-y0)/(Ny-1);

  print(file,"%16d", Nx);
  for(int i=0; i<Nx; i++)
    print(file,"  %16e", x0+i*dx); 
  print(file,"\n");

  for(int j=0; j<(int)result.size(); j++) {
    print(file,"%16e", y0+j*dy);
    for(int i=0; i<(int)result[j].size(); i++)
      print(file,"  %16e", result[j][i]); 
    print(file,"\n");
  }

  fclose(file);
  delete file;

}

//------------------------------------------------------------------------------

