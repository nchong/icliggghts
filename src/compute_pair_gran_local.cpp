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
#include "string.h"
#include "compute_pair_gran_local.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "pair_gran.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "pair_gran.h"

using namespace LAMMPS_NS;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

ComputePairGranLocal::ComputePairGranLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all("Illegal compute pair/gran/local command");

  local_flag = 1;
  nmax = 0;
  array = NULL;

  //store everything by default
  posflag = idflag = fflag = tflag = hflag = 1;

  //if further args, store only the properties that are listed
  if(narg > 3)
     posflag = idflag = fflag = tflag = hflag = 0;

  for (int iarg = 3; iarg < narg; iarg++) {
    int i = iarg-3;
    if (strcmp(arg[iarg],"pos") == 0) posflag = 1;
    else if (strcmp(arg[iarg],"id") == 0) idflag = 1;
    else if (strcmp(arg[iarg],"force") == 0) fflag = 1;
    else if (strcmp(arg[iarg],"torque") == 0) tflag = 1;
    else if (strcmp(arg[iarg],"history") == 0) hflag = 1;
    else error->all("Invalid keyword in compute pair/local command");
  }
}

/* ---------------------------------------------------------------------- */

ComputePairGranLocal::~ComputePairGranLocal()
{
  memory->destroy_2d_double_array(array);
  pairgran->unregister_compute_pair_local(this);
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::post_create()
{
  if (force->pair == NULL)
    error->all("No pair style is defined for compute pair/gran/local");

  pairgran = (PairGran*)force->pair_match("gran/",0);

  if (pairgran == NULL)
    error->all("No valid granular pair style found for use with compute pair/gran/local");

  if (pairgran->cpl_enable == 0)
    error->all("Pair style does not support compute pair/gran/local");

  if (idflag && atom->tag_enable == 0)
      error->all("Compute pair/gran/local requested to compute IDs, this requires atoms have IDs.");

  pairgran->register_compute_pair_local(this,dnum);

  if(hflag && dnum ==0) error->warning("Compute pair/gran/local can not calculate pair history values since pair style does not compute them");

  //14 standard values: pos1,pos2,id1,id2,force,torque
  nvalues = posflag*6 + idflag*2 + fflag*3 + tflag*3 + dnum;
  size_local_cols = nvalues;

}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::init()
{
  
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->gran = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // count local entries and compute pair info

  ncount = count_pairs();
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;

  ipair = 0;
  pairgran->compute(0,0,0);
}

/* ----------------------------------------------------------------------
   count pairs on this proc
------------------------------------------------------------------------- */

int ComputePairGranLocal::count_pairs()
{
  int i,j,m,n,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double *fbuf,*tbuf,*hbuf;

  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  neighbor->build_one(list->index);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  // skip if I or J are not in group

  m = n = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      if (j >= nall) j %= nall;

      if (!(mask[j] & groupbit)) continue;
      if (newton_pair == 0 && j >= nlocal) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      if (rsq >= (radius[i]+radius[j])*(radius[i]+radius[j])) continue;
      m++;
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::add_pair(int i,int j,double fx,double fy,double fz,double tor1,double tor2,double tor3,double *hist)
{
    
    if (!(atom->mask[i] & groupbit)) return;
    if (!(atom->mask[j] & groupbit)) return;

    int n = 0;
    if(posflag)
    {
        array[ipair][n++] = atom->x[i][0];
        array[ipair][n++] = atom->x[i][1];
        array[ipair][n++] = atom->x[i][2];
        array[ipair][n++] = atom->x[j][0];
        array[ipair][n++] = atom->x[j][1];
        array[ipair][n++] = atom->x[j][2];
    }
    if(idflag)
    {
        array[ipair][n++] = static_cast<double>(atom->tag[i]);
        array[ipair][n++] = static_cast<double>(atom->tag[j]);
    }
    if(fflag)
    {
        array[ipair][n++] = fx;
        array[ipair][n++] = fy;
        array[ipair][n++] = fz;
    }
    if(tflag)
    {
        array[ipair][n++] = tor1;
        array[ipair][n++] = tor2;
        array[ipair][n++] = tor3;
    }
    if(hflag)
    {
        for(int d = 0; d < dnum; d++)
           array[ipair][n++] = hist[d];
    }
    ipair++;
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::reallocate(int n)
{
  // grow vector or array and indices array

  while (nmax < n) nmax += DELTA;

  memory->destroy_2d_double_array(array);
  array = memory->create_2d_double_array(nmax,nvalues,"pair/local:array");
  array_local = array;

}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputePairGranLocal::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  return bytes;
}
