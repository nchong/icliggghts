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
#include "fix_meshGran_analyze.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "fix_gravity.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "neighbor.h"
#include "triSpherePenetration.h"
#include "stl_tri.h"
#include "mpi.h"
#include "fix_propertyGlobal.h"

using namespace LAMMPS_NS;

#define EPSILON 0.001

#define FABS(a) ((a) > 0 ? (a) : (-a))

/* ---------------------------------------------------------------------- */

FixMeshGranAnalyze::FixMeshGranAnalyze(LAMMPS *lmp, int narg, char **arg) :
  FixMeshGran(lmp, narg, arg)
{
   analyseStress=true;

   vector_flag = 1;
   size_vector = 6;
   global_freq = 1;
   extvector = 1;

   finnie_flag = 0;
   k_finnie = NULL;

   bool hasargs = true;
   while(iarg < narg && hasargs)
   {
      hasargs = false;
      if(strcmp(arg[iarg],"finnie") == 0)
      {
          if (narg < iarg+2) error->all("Illegal fix mesh/gran/stressanalysis command, not enough arguments");
          iarg++;
          if(strcmp(arg[iarg],"yes") == 0) finnie_flag = 1;
          else if(strcmp(arg[iarg],"no") == 0) finnie_flag = 0;
          else error->all("Illegal fix mesh/gran/stressanalysis command, expecting 'yes' or 'no'");
          iarg++;
          hasargs = true;
      }
      else if(strcmp(style,"mesh/gran/stressanalysis") == 0) error->all("Illegal fix mesh/gran/stressanalysis command, unknown keyword");
   }

   if(finnie_flag) error->warning("You are using the wear model, which is currently in beta mode!");
}

FixMeshGranAnalyze::~FixMeshGranAnalyze()
{

}

/* ---------------------------------------------------------------------- */

int FixMeshGranAnalyze::setmask()
{
  int mask = 0;
  mask |=PRE_FORCE;
  mask |=FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMeshGranAnalyze::init()
{
    if(finnie_flag)
        k_finnie = static_cast<FixPropertyGlobal*>(modify->fix[modify->find_fix_property("k_finnie","property/global","peratomtypepair",atom->ntypes,atom->ntypes)])->get_array();
}

/* ---------------------------------------------------------------------- */

void FixMeshGranAnalyze::pre_force(int vflag)
{
    force_total[0]=0.;  force_total[1]=0.;  force_total[2]=0.;
    torque_total[0]=0.; torque_total[1]=0.; torque_total[2]=0.;

    for(int i=0;i<nTri;i++)
        for(int j=0;j<3;j++)
            STLdata->f_tri[i][j]=0.;

}

/* ---------------------------------------------------------------------- */
//   called during wall force calc
/* ---------------------------------------------------------------------- */

void FixMeshGranAnalyze::add_particle_contribution(double *frc,double* contactPoint,int iTri,int ip)
{
    double E,c[3],cmag,vmag,cos_gamma,sin_gamma,sin_2gamma,tan_gamma;

    //do not include if not in fix group
    if(!(atom->mask[ip] & groupbit)) return;

    vectorAdd3D(STLdata->f_tri[iTri],frc,STLdata->f_tri[iTri]);

    vectorAdd3D(force_total,frc,force_total);
    vectorSubtract3D(contactPoint,p_ref,tmp);
    vectorCross3D(tmp,frc,tmp2); //tmp2 is torque contrib
    vectorAdd3D(torque_total,tmp2,torque_total);

    if(finnie_flag)
    {
        
        vectorSubtract3D(contactPoint,atom->x[ip],c);
        cmag = vectorMag3D(c);

        if(vectorDot3D(c,atom->v[ip]) < 0.) return;

        vmag = vectorMag3D(atom->v[ip]);

        sin_gamma = FABS(vectorDot3D(atom->v[ip],STLdata->facenormal[iTri])) / (vmag);
        
        if(sin_gamma > 1.) sin_gamma = 1.;
        if(sin_gamma < -1.) sin_gamma = -1.;

        cos_gamma = sqrt(1. - sin_gamma * sin_gamma);
        if(cos_gamma > EPSILON || (sin_gamma < 3. * cos_gamma))
        {
            E = 0.33333 * cos_gamma * cos_gamma;
            
        }
        else
        {
            sin_2gamma = 2. * sin_gamma * cos_gamma;
            E = sin_2gamma - 3. * sin_gamma * sin_gamma;
            
        }
        E *= 2.*k_finnie[atom_type_wall-1][atom->type[ip]-1] * vmag * vectorMag3D(frc);
        
        STLdata->wear_step[iTri] += E;
    }
}

/* ---------------------------------------------------------------------- */

void FixMeshGranAnalyze::final_integrate()
{
    calc_total_force();
}

/* ---------------------------------------------------------------------- */

void FixMeshGranAnalyze::calc_total_force()
{
    int nTri = STLdata->nTri;

    //total force on tri
    double *force_total_all=new double[3];
    double *torque_total_all=new double[3];

    MPI_Allreduce(force_total,force_total_all,3, MPI_DOUBLE, MPI_SUM,world);
    MPI_Allreduce(torque_total,torque_total_all,3, MPI_DOUBLE, MPI_SUM,world);

    for(int j=0;j<3;j++)
    {
        force_total[j]=force_total_all[j];
        torque_total[j]=torque_total_all[j];
    }

    delete []force_total_all;
    delete []torque_total_all;

    double *wear = STLdata->wear;
    double *wear_step = STLdata->wear_step;
    double *wear_step_all = new double[nTri];
    MPI_Allreduce(wear_step,wear_step_all,nTri, MPI_DOUBLE, MPI_SUM,world);

    for(int i=0;i<nTri;i++)
    {
        wear[i] += wear_step_all[i];
        wear_step[i] = 0.;
    }
    delete []wear_step_all;

    //forces on tri

    MPI_Allreduce(&(STLdata->f_tri[0][0]),&(STLdata->fn_fshear[0][0]),3*nTri, MPI_DOUBLE, MPI_SUM,world);

    //switch fn_fshear and f_tri
    double **helper;
    helper=STLdata->f_tri;
    STLdata->f_tri=STLdata->fn_fshear;
    STLdata->fn_fshear=helper;

    double temp[3];
    double p,s;
    for(int i=0;i<nTri;i++)
    {
        //pressure
        STLdata->fn_fshear[i][0]=vectorDot3D(STLdata->f_tri[i],STLdata->facenormal[i]);
        vectorScalarMult3D(STLdata->facenormal[i],STLdata->fn_fshear[i][0],temp);
        vectorSubtract3D(STLdata->f_tri[i],temp,temp);
        //shear force
        STLdata->fn_fshear[i][1]=vectorMag3D(temp);
    }
}

/* ----------------------------------------------------------------------
   return force/torque on body
------------------------------------------------------------------------- */

double FixMeshGranAnalyze::compute_vector(int n)
{
  //return force/torque
  if(n<3) return force_total[n];
  else    return torque_total[n-3];
}
