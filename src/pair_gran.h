/* ----------------------------------------------------------------------
LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
Transfer Simulations

www.liggghts.com | www.cfdem.com
Christoph Kloss, christoph.kloss@cfdem.com

LIGGGHTS is based on LAMMPS
LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
http://lammps.sandia.gov, Sandia National Laboratories
Steve Plimpton, sjplimp@sandia.gov

Copyright (2003) Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software. This software is distributed under
the GNU General Public License.

See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors for original version: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

#else

#ifndef LMP_PAIR_GRAN_H
#define LMP_PAIR_GRAN_H

#include "pair.h"

namespace LAMMPS_NS {

class PairGran : public Pair {
 public:

  friend class FixWallGran;
  friend class FixCheckTimestepGran;
  friend class ComputePairGranLocal;

  PairGran(class LAMMPS *);
  ~PairGran();
  virtual void compute(int, int);
  virtual void compute(int, int,int)=0;
  virtual void settings(int, char **) =0;
  virtual void coeff(int, char **);
  virtual void init_style();
  virtual void init_substyle() {} 
  virtual void init_list(int, class NeighList *);
  virtual double init_one(int, int);
  virtual void write_restart(FILE *);
  virtual void read_restart(FILE *);
  virtual void write_restart_settings(FILE *){}
  virtual void read_restart_settings(FILE *){}
  virtual void reset_dt();
  int is_history(){return history;}

  class MechParamGran *mpg;

 protected:

  virtual void history_args(char**) =0;
  void register_compute_pair_local(class ComputePairGranLocal *,int&);
  void unregister_compute_pair_local(class ComputePairGranLocal *ptr);

  int cpl_enable;
  class FixRigid* fr;
  class ComputePairGranLocal *cpl;

  double dt;
  int freeze_group_bit;

  int history;
  int dnum;
  class FixContactHistory *fix_history;
  int shearupdate;
  int laststep;

  double *onerad_dynamic,*onerad_frozen;
  double *maxrad_dynamic,*maxrad_frozen;

  void allocate();
};

}

#endif
#endif
