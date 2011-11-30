#include "assert.h"
#include "stdlib.h"
#include <sstream>

#include "atom.h"
#include "error.h"
#include "fix_emit_step.h"
#include "force.h"
#include "mech_param_gran.h"
#include "neigh_list.h"
#include "pair.h"
#include "pair_gran_hertz_history.h"

#include <fstream>

using namespace LAMMPS_NS;
using namespace std;

FixEmitStep::FixEmitStep(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {
  if (!(force->pair_match("gran/hertz/history",0))) {
    error->all("Fix emit_step can only be used with pair style gran/hertz/history");
  }
  nevery = atoi(arg[3]);
  step = 0;
}

FixEmitStep::~FixEmitStep() {
}

int FixEmitStep::setmask() {
  int mask = PRE_FORCE | POST_FORCE;
  return mask;
}

bool FixEmitStep::trigger() {
  return (step > 0 && step % nevery == 0);
}

void FixEmitStep::pre_force(int vflag) {
  if (trigger()) {
    std::stringstream fname;
    fname << step << ".step";
    emit(fname.str().c_str(), true);
  }
}

void FixEmitStep::post_force(int vflag) {
  if (trigger()) {
    std::stringstream fname;
    fname << step << ".step";
    emit(fname.str().c_str(), false);
  }
  step++;
}

void FixEmitStep::emit(const char *fname, bool is_pre_force=true) {
  PairGranHertzHistory *pair = static_cast<PairGranHertzHistory*>(force->pair);
  class NeighList *list = pair->list;
  class NeighList *listgranhistory = pair->listgranhistory;

  int inum = list->inum;
  int nall = atom->nlocal + atom->nghost;

  ofstream ofile;
  ofile.open(fname, 
      ofstream::binary | (is_pre_force ? ofstream::out : ofstream::app));

  if (is_pre_force) {
    //MAGIC VALUE
    unsigned int MAGIC = 0xDEADBEEF;
    ofile.write(reinterpret_cast<char *>(&MAGIC), sizeof(MAGIC));

    //CONSTANTS
    ofile.write(reinterpret_cast<char *>(&(pair->dt)), sizeof(pair->dt));
    ofile.write(reinterpret_cast<char *>(&(force->nktv2p)), sizeof(force->nktv2p));
    int ntype = pair->mpg->max_type() + 1;
    ofile.write(reinterpret_cast<char *>(&(ntype)), sizeof(ntype));
    for (int p=0; p<ntype; p++) {
      for (int q=0; q<ntype; q++) {
        ofile.write(reinterpret_cast<char *>(&(pair->Yeff[p][q])), sizeof(pair->Yeff[p][q]));
      }
    }
    for (int p=0; p<ntype; p++) {
      for (int q=0; q<ntype; q++) {
        ofile.write(reinterpret_cast<char *>(&(pair->Geff[p][q])), sizeof(pair->Geff[p][q]));
      }
    }
    for (int p=0; p<ntype; p++) {
      for (int q=0; q<ntype; q++) {
        ofile.write(reinterpret_cast<char *>(&(pair->betaeff[p][q])), sizeof(pair->betaeff[p][q]));
      }
    }
    for (int p=0; p<ntype; p++) {
      for (int q=0; q<ntype; q++) {
        ofile.write(reinterpret_cast<char *>(&(pair->coeffFrict[p][q])), sizeof(pair->coeffFrict[p][q]));
      }
    }
        
    //NODES
    ofile.write(reinterpret_cast<char *>(&(nall)), sizeof(nall));
    for (int i=0; i<nall; i++) {
      for (int j=0; j<3; j++) {
        ofile.write(reinterpret_cast<char *>(&(atom->x[i][j])), sizeof(atom->x[i][j]));
      }
    }
    for (int i=0; i<nall; i++) {
      for (int j=0; j<3; j++) {
        ofile.write(reinterpret_cast<char *>(&(atom->v[i][j])), sizeof(atom->v[i][j]));
      }
    }
    for (int i=0; i<nall; i++) {
      for (int j=0; j<3; j++) {
        ofile.write(reinterpret_cast<char *>(&(atom->omega[i][j])), sizeof(atom->omega[i][j]));
      }
    }
    for (int i=0; i<nall; i++) {
      ofile.write(reinterpret_cast<char *>(&(atom->radius[i])), sizeof(atom->radius[i]));
    }
    for (int i=0; i<nall; i++) {
      ofile.write(reinterpret_cast<char *>(&(atom->rmass[i])), sizeof(atom->rmass[i]));
    }
    for (int i=0; i<nall; i++) {
      ofile.write(reinterpret_cast<char *>(&(atom->type[i])), sizeof(atom->type[i]));
    }
  }
  for (int i=0; i<nall; i++) {
    for (int j=0; j<3; j++) {
      ofile.write(reinterpret_cast<char *>(&(atom->f[i][j])), sizeof(atom->f[i][j]));
    }
  }
  for (int i=0; i<nall; i++) {
    for (int j=0; j<3; j++) {
      ofile.write(reinterpret_cast<char *>(&(atom->torque[i][j])), sizeof(atom->torque[i][j]));
    }
  }

  //EDGES
  if (is_pre_force) {
    int nedge = 0;
    for (int ii=0; ii<inum; ii++) {
      int i = list->ilist[ii];
      int jnum = list->numneigh[i];
      nedge += jnum;
    }
    ofile.write(reinterpret_cast<char *>(&(nedge)), sizeof(nedge));

    for (int ii=0; ii<inum; ii++) {
      int i = list->ilist[ii];
      int jnum = list->numneigh[i];

      for (int jj = 0; jj<jnum; jj++) {
        int j = list->firstneigh[i][jj];
        ofile.write(reinterpret_cast<char *>(&(i)), sizeof(i));
        ofile.write(reinterpret_cast<char *>(&(j)), sizeof(j));
      }
    }
  }

  for (int ii=0; ii<inum; ii++) {
    int i = list->ilist[ii];
    int jnum = list->numneigh[i];

    for (int jj = 0; jj<jnum; jj++) {
      double *shear = &(listgranhistory->firstdouble[i][3*jj]);
      for (int k=0; k<3; k++) {
        ofile.write(reinterpret_cast<char *>(&(shear[k])), sizeof(shear[k]));
      }
    }
  }
  ofile.close();
}
