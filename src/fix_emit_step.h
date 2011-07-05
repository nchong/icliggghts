/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(emit_step,FixEmitStep)

#else

#ifndef LMP_FIX_EMIT_STEP_H
#define LMP_FIX_EMIT_STEP_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixEmitStep : public Fix {
 public:
  FixEmitStep(class LAMMPS *, int, char **);
  ~FixEmitStep();
  int setmask();
  void pre_force(int);
  void post_force(int);

 private:
  int step;   //number of simulation steps
  bool trigger(void);
  void emit(const char *, bool);
};

}

#endif
#endif
