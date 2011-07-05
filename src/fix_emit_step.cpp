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

using namespace LAMMPS_NS;

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

  FILE *ofile = is_pre_force ? fopen(fname, "w") : fopen(fname, "a");

  if (is_pre_force) {
    //CONSTANTS
    fprintf(ofile, "%.16f\n", pair->dt);
    fprintf(ofile, "%.16f\n", force->nktv2p);
    int ntype = pair->mpg->max_type() + 1;
    fprintf(ofile, "%d\n", ntype);
    for (int p=0; p<ntype; p++) {
      for (int q=0; q<ntype; q++) {
        fprintf(ofile, "%.16f\n", pair->Yeff[p][q]);
      }
    }
    for (int p=0; p<ntype; p++) {
      for (int q=0; q<ntype; q++) {
        fprintf(ofile, "%.16f\n", pair->Geff[p][q]);
      }
    }
    for (int p=0; p<ntype; p++) {
      for (int q=0; q<ntype; q++) {
        fprintf(ofile, "%.16f\n", pair->betaeff[p][q]);
      }
    }
    for (int p=0; p<ntype; p++) {
      for (int q=0; q<ntype; q++) {
        fprintf(ofile, "%.16f\n", pair->coeffFrict[p][q]);
      }
    }
        
    //NODES
    fprintf(ofile, "%d\n", nall);
    for (int i=0; i<nall; i++) {
      for (int j=0; j<3; j++) {
        fprintf(ofile, "%.16f\n", atom->x[i][j]);
      }
    }
    for (int i=0; i<nall; i++) {
      for (int j=0; j<3; j++) {
        fprintf(ofile, "%.16f\n", atom->v[i][j]);
      }
    }
    for (int i=0; i<nall; i++) {
      for (int j=0; j<3; j++) {
        fprintf(ofile, "%.16f\n", atom->omega[i][j]);
      }
    }
    for (int i=0; i<nall; i++) {
      fprintf(ofile, "%.16f\n", atom->radius[i]);
    }
    for (int i=0; i<nall; i++) {
      fprintf(ofile, "%.16f\n", atom->rmass[i]);
    }
    for (int i=0; i<nall; i++) {
      fprintf(ofile, "%d\n", atom->type[i]);
    }
  }
  for (int i=0; i<nall; i++) {
    for (int j=0; j<3; j++) {
      fprintf(ofile, "%.16f\n", atom->f[i][j]);
    }
  }
  for (int i=0; i<nall; i++) {
    for (int j=0; j<3; j++) {
      fprintf(ofile, "%.16f\n", atom->torque[i][j]);
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
    fprintf(ofile, "%d\n", nedge);

    for (int ii=0; ii<inum; ii++) {
      int i = list->ilist[ii];
      int jnum = list->numneigh[i];

      for (int jj = 0; jj<jnum; jj++) {
        int j = list->firstneigh[i][jj];
        fprintf(ofile, "%d\n", i);
        fprintf(ofile, "%d\n", j);
      }
    }
  }

  for (int ii=0; ii<inum; ii++) {
    int i = list->ilist[ii];
    int jnum = list->numneigh[i];

    for (int jj = 0; jj<jnum; jj++) {
      double *shear = &(listgranhistory->firstdouble[i][3*jj]);
      for (int k=0; k<3; k++) {
        fprintf(ofile, "%.16f\n", shear[k]);
      }
    }
  }
  fflush(ofile);
  fclose(ofile);
}
