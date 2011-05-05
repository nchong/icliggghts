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

#ifdef FIX_CLASS

FixStyle(wall/gran,FixWallGran)

#else

#ifndef LMP_FIX_WALL_GRAN_H
#define LMP_FIX_WALL_GRAN_H

#include "fix.h"

enum{MESHGRAN_FWG,XPLANE_FWG,YPLANE_FWG,ZPLANE_FWG,ZCYLINDER_FWG};

#define F_SHRINKAGE -0.000000001  

namespace LAMMPS_NS {

class FixWallGran : public Fix {

 friend class FixCheckTimestepGran;
 friend class FixTriNeighlist;
 friend class MechParamGran;

 public:
  FixWallGran(class LAMMPS *, int, char **);
  ~FixWallGran();
  virtual int setmask();
  void post_create();
  void pre_delete();
  void init();
  virtual void init_substyle(){}
  void setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);

  double memory_usage();
  void grow_arrays(int);
  virtual void grow_arrays_maxtritouch(int); 
  void copy_arrays(int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  void write_restart(FILE *);
  void restart(char *);      
  int size_restart(int);
  int maxsize_restart();
  void reset_dt();

  int n_contacts(); 

 protected:
  virtual void init_heattransfer();
  void registerTriNeighlist();

  void reset_wall_forces();
  int  add_to_contact_list(int, int, int);
  void history_transition(int,int,int);
  void reset_history(int,int);
  void copy_history(int,int);
  void remove_from_contact_list(int, int, int);
  void remove_from_contact_list_ext(int, int, int);

  virtual void compute_force(int i,double deltan,double rsq,double meff_wall,double dx,double dy,double dz,double *vwall,double *c_history,double area_ratio) {}
  virtual void addHeatFlux(int i,double rsq,double area_ratio){}

  int atom_type_wall;
  int wallstyle,wiggle,wshear,axis;
  double lo,hi,cylradius;
  double amplitude,period,omega,vshear;
  double dt;
  int nlevels_respa;
  int time_origin;

  int meshwall;
  int nFixMeshGran;
  class FixMeshGran **FixMeshGranList;
  class FixTriNeighlist *fix_tri_neighlist;

  char *pairstyle;
  class PairGran *pairgran;
  class FixRigid *fr;

  int history;
  int dnum;
  int maxpartners; 
  int *npartners;
  int ***partner;
  double ***contacthistory;
  int shearupdate;
  int laststep;

  double Temp_wall;
  class FixPropertyPerAtom *fppa_T;
  class FixPropertyPerAtom *fppa_hf;
  double *Temp_p;
  double *heatflux;
  const double *th_cond;
  double const* const* deltan_ratio;
};

}

#endif
#endif
