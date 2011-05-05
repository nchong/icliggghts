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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_heat_gran.h"
#include "atom.h"
#include "domain.h"
#include "group.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair_gran.h"
#include "math_extra.h"
#include "fix_propertyGlobal.h"
#include "fix_propertyPerAtom.h"
#include "fix_scalar_transport_equation.h"
#include "mech_param_gran.h"
#include "respa.h"

using namespace LAMMPS_NS;

#define SMALL 1e-8

/* ---------------------------------------------------------------------- */

FixHeatGran::FixHeatGran(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{

  if ((!atom->radius_flag)||(!atom->rmass_flag)) error->all("Fix heat/gran needs per particle radius and mass");

  //check if a fix of this style already exists
  for (int i=0;i<modify->nfix;i++)
      if (strcmp(modify->fix[i]->style,style) == 0) error->all("A fix of type heat/gran is already registered. Can not have more than one");

  if (narg < 4) error->all("Illegal fix heat/gran command, not enough arguments");
  T0 = atof(arg[3]);

  fix_temp = fix_heatFlux = fix_heatSource = NULL;
  fix_conductivity = NULL;

  conductivity = NULL;

  peratom_flag = 1;              
  size_peratom_cols = 0;         
  peratom_freq = 1;
  time_depend = 1;

  scalar_flag = 1; 
  global_freq = 1; 
}

/* ---------------------------------------------------------------------- */

FixHeatGran::~FixHeatGran()
{
    
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::post_create()
{
    fix_ste = modify->find_fix_scalar_transport_equation("heattransfer");

    if(!fix_ste)
    {
        char **newarg = new char*[15];
        newarg[0] = (char *) "ste_heattransfer";
        newarg[1] = group->names[igroup];
        newarg[2] = (char *) "transportequation/scalar";
        newarg[3] = (char *) "equation_id";
        newarg[4] = (char *) "heattransfer";
        newarg[5] = (char *) "quantity";
        newarg[6] = (char *) "Temp";
        newarg[7] = (char *) "default_value";
        newarg[8] = new char[30];
        sprintf(newarg[8],"%f",T0);
        newarg[9] = (char *) "flux_quantity";
        newarg[10] = (char *) "heatFlux";
        newarg[11] = (char *) "source_quantity";
        newarg[12] = (char *) "heatSource";
        newarg[13] = (char *) "capacity_quantity";
        newarg[14] = (char *) "thermalCapacity";
        modify->add_fix(15,newarg);
        delete [] newarg[8];
        delete [] newarg;
    }
}

/* ---------------------------------------------------------------------- */

int FixHeatGran::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}
/* ---------------------------------------------------------------------- */
void FixHeatGran::updatePtrs()
{
  Temp = fix_temp->vector_atom;
  vector_atom = Temp; 

  heatFlux = fix_heatFlux->vector_atom;
  heatSource = fix_heatSource->vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::init()
{
  if (!atom->radius_flag || !atom->rmass_flag) error->all("Please use a granular atom style for fix heat/gran");

  if(!force->pair_match("gran", 0)) error->all("Please use a granular pair style for fix heat/gran");

  pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
  history_flag = pair_gran->is_history();

  fix_ste = modify->find_fix_scalar_transport_equation("heattransfer");
  if(!fix_ste) error->all("Fix heat/gran needs a fix transportequation/scalar to work with");

  fix_temp = static_cast<FixPropertyPerAtom*>(modify->fix[modify->find_fix_property("Temp","property/peratom","scalar",0,0)]);
  fix_heatFlux = static_cast<FixPropertyPerAtom*>(modify->fix[modify->find_fix_property("heatFlux","property/peratom","scalar",0,0)]);
  fix_heatSource = static_cast<FixPropertyPerAtom*>(modify->fix[modify->find_fix_property("heatSource","property/peratom","scalar",0,0)]);

  int max_type = pair_gran->mpg->max_type();
  if(conductivity) delete []conductivity;
  conductivity = new double[max_type+1];
  fix_conductivity = static_cast<FixPropertyGlobal*>(modify->fix[modify->find_fix_property("thermalConductivity","property/global","peratomtype",atom->ntypes,0)]);

  //pre-calculate parameters for possible contact material combinations
  for(int i=1;i< max_type+1; i++)
      for(int j=1;j<max_type+1;j++)
          conductivity[i] = fix_conductivity->compute_vector(i-1);

  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::post_force(int vflag)
{
    //template function for using touchflag or not
    if(history_flag == 0) post_force_eval<0>(vflag);
    if(history_flag == 1) post_force_eval<1>(vflag);
}

/* ---------------------------------------------------------------------- */

template <int HISTFLAG>
void FixHeatGran::post_force_eval(int vflag)
{
  double hc,contactArea;
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv,tcoi,tcoj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;

  int newton_pair = force->newton_pair;

  if (strcmp(force->pair_style,"hybrid")==0) error->warning("Fix heat/gran implementation may not be valid for pair style hybrid");
  if (strcmp(force->pair_style,"hybrid/overlay")==0) error->warning("Fix heat/gran implementation may not be valid for pair style hybrid/overlay");

  inum = pair_gran->list->inum;
  ilist = pair_gran->list->ilist;
  numneigh = pair_gran->list->numneigh;
  firstneigh = pair_gran->list->firstneigh;
  if(HISTFLAG) firsttouch = pair_gran->listgranhistory->firstneigh;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  updatePtrs();

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    if(HISTFLAG) touch = firsttouch[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      if (!(mask[i] & groupbit && mask[j] & groupbit)) continue;

      if(!HISTFLAG)
      {
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          radj = radius[j];
          radsum = radi + radj;
      }

      if (HISTFLAG && touch[jj] || !HISTFLAG && (rsq < radsum*radsum)) {  //contact
         
         if(HISTFLAG)
         {
              delx = xtmp - x[j][0];
              dely = ytmp - x[j][1];
              delz = ztmp - x[j][2];
              rsq = delx*delx + dely*dely + delz*delz;
              radj = radius[j];
              radsum = radi + radj;
         }

         r = sqrt(rsq);
         contactArea = - M_PI/4 * ( (r-radi-radj)*(r+radi-radj)*(r-radi+radj)*(r+radi+radj) )/(r*r); //contact area of the two spheres
         tcoi = conductivity[type[i]];
         tcoj = conductivity[type[j]];

         if ((fabs(tcoi) < SMALL) || (fabs(tcoj) < SMALL)) hc = 0.;
         else hc=4.*tcoi*tcoj/(tcoi+tcoj)*sqrt(contactArea);

         heatFlux[i] += (Temp[j]-Temp[i])*hc;
         if (newton_pair||j<nlocal) heatFlux[j] += (Temp[i]-Temp[j])*hc;
      }
    }
  }

}

/* ---------------------------------------------------------------------- */
double FixHeatGran::compute_scalar()
{
    return fix_ste->compute_scalar();
}
