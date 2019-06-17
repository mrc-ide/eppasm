// Copyright (C) 2019  Kinh Nguyen

// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.

// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.

// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
#include "progression.hpp"

extern "C" SEXP eppasmOOpp(SEXP fp, SEXP MODEL, SEXP MIX) {
  int cMODEL = INTEGER_VALUE(MODEL); bool cMIX = LOGICAL_VALUE(MIX);
  popC pop(fp, cMODEL, cMIX);
  hivC hivpop(fp, cMODEL);
  artC artpop(fp, cMODEL);
  for (int i = 1; i < pop.p.SIM_YEARS; ++i) {
    pop.year = i; hivpop.year = i; artpop.year = i;
    epp_aging(pop, hivpop, artpop);
    epp_death(pop, hivpop, artpop);
    epp_migration(pop, hivpop, artpop);
    pop.update_fertile();
    if (cMODEL != 0) { // Disease model simulation: events at dt timestep
      epp_disease_model(pop, hivpop, artpop);
      if (pop.p.eppmod == 2) //// Direct incidence input model
        pop.epp_disease_model_direct(hivpop, artpop);
    }
    if (pop.p.popadjust) { // match target pop
      pop.adjust_pop();
      if (cMODEL != 0) {
        hivpop.adjust_pop(pop.adj_prob);
        if (i >= pop.p.tARTstart - 1)
          artpop.adjust_pop(pop.adj_prob);
      }
    }
    if (cMODEL != 0) {
      if (i + pop.AGE_START <= pop.PROJ_YEARS - 1)
        pop.cal_prev_pregant(hivpop, artpop); // prevalence among pregnant women
      pop.save_prev_n_inc(); // save prevalence and incidence 15 to 49
    }
  }
  pop.finalize(hivpop, artpop);
  return pop.pop_sexp;
}