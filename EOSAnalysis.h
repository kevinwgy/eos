/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _EOS_ANALYSIS_
#define _EOS_ANALYSIS_

#include<VarFcnBase.h>
#include<IoData.h>
#include<vector>


/************************************************************************
 * class EOSAnalysis contains tools for analyzing EOS, including creating
 * 1D or 2D tables.
 ***********************************************************************/

class EOSAnalysis {

  std::vector<VarFcnBase*>& vf;

  ObjectMap<EOSTabulationData> &eos_tabulationMap; 

public:

  EOSAnalysis(ObjectMap<EOSTabulationData> &eos_tabulationMap_, std::vector<VarFcnBase*>& vf_);
  ~EOSAnalysis();

  void GenerateAllEOSTables();

private:

  void GenerateEOSTable(EOSTabulationData& tab);

};


#endif
