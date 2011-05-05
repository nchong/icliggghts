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

#include "particleToInsert.h"
#include "math.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;

ParticleToInsert::ParticleToInsert(LAMMPS* lmp,int ns) : Pointers(lmp)
{
        nspheres=ns;
        x_ins=memory->create_2d_double_array(nspheres,3,"x_ins");
        radius_ins=new double[nspheres];
        x_bound=new double[3];
        xcm=new double[3];
        inertia=new double[3];
        ex_space=new double[3];
        ey_space=new double[3];
        ez_space=new double[3];
        displace=memory->create_2d_double_array(nspheres,3,"displace");
}

ParticleToInsert::~ParticleToInsert()
{
        memory->destroy_2d_double_array(x_ins);
        delete []radius_ins;
        delete []x_bound;
        delete []xcm;
        delete []inertia;
        delete []ex_space;
        delete []ey_space;
        delete []ez_space;
        memory->destroy_2d_double_array(displace);
}

void ParticleToInsert::random_rotate(double rn1,double rn2, double rn3)
{
    
    if(nspheres==1)return;

    double *vert_before_rot;
    double vert_after_rot[3];

    double phix=rn1*2.*M_PI;
    double phiy=rn2*2.*M_PI;
    double phiz=rn3*2.*M_PI;

    double cos_phix = cos(phix);
    double cos_phiy = cos(phiy);
    double cos_phiz = cos(phiz);
    double sin_phix = sin(phix);
    double sin_phiy = sin(phiy);
    double sin_phiz = sin(phiz);

    for(int i=0;i<3;i++)
    {
        if     (i==0) vert_before_rot=ex_space;
        else if(i==1) vert_before_rot=ey_space;
        else if(i==2) vert_before_rot=ez_space;

        vert_after_rot[0] = vert_before_rot[0]*cos_phiy*cos_phiz+vert_before_rot[1]*(cos_phiz*sin_phix*sin_phiy-cos_phix*sin_phiz)+vert_before_rot[2]*(cos_phix*cos_phiz*sin_phiy+sin_phix*sin_phiz);
        vert_after_rot[1] = vert_before_rot[0]*cos_phiy*sin_phiz+vert_before_rot[2]*(-cos_phiz*sin_phix+cos_phix*sin_phiy*sin_phiz)+vert_before_rot[1]*(cos_phix*cos_phiz+sin_phix*sin_phiy*sin_phiz);
        vert_after_rot[2] = vert_before_rot[2]*cos_phix*cos_phiy+vert_before_rot[1]*cos_phiy*sin_phix-vert_before_rot[0]*sin_phiy;

        if     (i==0) for(int j=0;j<3;j++) ex_space[j]=vert_after_rot[j];
        else if(i==1) for(int j=0;j<3;j++) ey_space[j]=vert_after_rot[j];
        else if(i==2) for(int j=0;j<3;j++) ez_space[j]=vert_after_rot[j];
    }

    for(int i=0;i<nspheres;i++)
    {
        x_ins[i][0] = xcm[0] + ex_space[0]*displace[i][0] +   ey_space[0]*displace[i][1] +   ez_space[0]*displace[i][2];
        x_ins[i][1] = xcm[1] + ex_space[1]*displace[i][0] +   ey_space[1]*displace[i][1] +   ez_space[1]*displace[i][2];
        x_ins[i][2] = xcm[2] + ex_space[2]*displace[i][0] +   ey_space[2]*displace[i][1] +   ez_space[2]*displace[i][2];
    }

}
