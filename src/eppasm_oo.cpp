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
#include "Classes.h"
#include "Rdefines.h"
using namespace arma;

// Natural deaths
// -----------------------------------------------------------------------------
void epp_death (popC& pop, hivC& hivpop, artC& artpop) {
  pop.deaths();
  if (pop.MODEL!=0) {
    hivpop.deaths(pop.hiv_sx_prob);
    if (pop.year > pop.p.tARTstart - 1)
      artpop.deaths(pop.hiv_sx_prob);
  }
}

// Migration at year i
// -----------------------------------------------------------------------------
void epp_migration (popC& pop, hivC& hivpop, artC& artpop) {
  pop.migration();
  if (pop.MODEL!=0) {
    hivpop.migration(pop.hiv_mr_prob);
    if (pop.year > pop.p.tARTstart - 1)
      artpop.migration(pop.hiv_mr_prob);
  }
}

// EPP populations aging
// -----------------------------------------------------------------------------
void epp_aging (popC& pop, hivC& hivpop, artC& artpop) {
  pop.aging();
  pop.add_entrants();
  if (pop.MODEL==2)
    pop.sexual_debut();
  if (pop.MODEL!=0) {
    mat hiv_ag_prob = pop.hiv_aging_prob();
    vec artYesNo = pop.entrant_art();
    hivpop.aging(hiv_ag_prob);
    hivpop.add_entrants(artYesNo);
    if (pop.MODEL==2)
      hivpop.sexual_debut();
    if (pop.year > pop.p.tARTstart - 1) {
      artpop.aging(hiv_ag_prob);
      artpop.add_entrants(artYesNo);
      if (pop.MODEL==2)
        artpop.sexual_debut();
    }
  }
}

// calculate, distribute eligible for ART, update grad, gradART
// -----------------------------------------------------------------------------
void popC::epp_art_init (hivC& hivpop, artC& artpop, int time_step) {
  artpop.grad_progress();
  artpop.art_dropout(hivpop); // pass hivpop to receive the drop out
  vec eligible = hivpop.eligible_for_art();
  cube art_elig = hivpop.data(year);
  for (uword sex = 0; sex < 2; ++sex)
    art_elig.slice(sex).each_col() %= eligible;
  if ( (p.pw_artelig(year)==1) & (p.artcd4elig_idx(year) > 0) )
    art_elig = update_preg(art_elig, hivpop, artpop); // add pregnant?
  if (MODEL==2) { // add sexual inactive but eligible for treatment
    cube art_elig_db = hivpop.data_db(year);
    for (uword sex = 0; sex < 2; ++sex)
      art_elig_db.slice(sex).each_col() %= eligible;
    art_elig += art_elig_db;
  }
  // calculate number to initiate ART and distribute
  vec art_curr = artpop.current_on_art();
  vec artnum_ii = artInit(art_curr, art_elig, time_step);
  vec art15plus_inits = artnum_ii - art_curr;
  art15plus_inits.elem( find(art15plus_inits < 0) ).zeros();
  cube artinit = artDist(art_elig, art15plus_inits);
  if (MODEL==1) {
    cube artinit_2 = hivpop.data(year) + DT * hivpop.grad;
    uvec greater = find(artinit > artinit_2); 
    artinit(greater) = artinit_2(greater);
  }
  if (MODEL==2) // split the number proportionally for active and idle pop
    artinit = hivpop.distribute_artinit(artinit, artpop);
  hivpop.grad -= artinit / DT;
  artpop.grad_init(artinit);
}

// Disease model
// -----------------------------------------------------------------------------
void epp_disease_model (popC& pop, hivC& hivpop, artC& artpop) {
  for (int time_step = 0; time_step < pop.hiv_steps_per_year; ++time_step) {
    if (pop.p.eppmod == 0) { // != "directincid"
      pop.update_rvec(time_step);
      mat infect;
      if (pop.MIX)
        infect = pop.infect_mix(time_step);
      else
        infect = pop.infect_spec(hivpop, artpop, time_step);
      pop.update_infection(infect);
      hivpop.update_infection(infect);
    }
    cube cd4_mort = pop.scale_cd4_mort(hivpop, artpop);
    hivpop.grad_progress(cd4_mort); // cd4 disease progression and mortality
    pop.remove_hiv_death(cd4_mort, hivpop, artpop); // Remove hivdeaths from pop
    if (pop.year >= pop.p.tARTstart - 1) // ART initiation
      pop.epp_art_init(hivpop, artpop, time_step);
    hivpop.add_grad_to_pop();
  } // end time step
}

// [[Rcpp::export]]
SEXP eppasmOO(SEXP fp, int MODEL = 1, bool MIX = true) {
  popC pop(fp, MODEL, MIX);
  hivC hivpop(fp, MODEL);
  artC artpop(fp, MODEL);
  for (uword i = 1; i < pop.p.SIM_YEARS; ++i) {
    pop.year = i, hivpop.year = i, artpop.year = i;
    epp_aging(pop, hivpop, artpop);
    epp_death(pop, hivpop, artpop);
    epp_migration(pop, hivpop, artpop);
    pop.update_fertile();
    if (MODEL!=0) { // Disease model simulation: events at dt timestep
      epp_disease_model(pop, hivpop, artpop);
    //   if (fp$eppmod == "directincid") //// Direct incidence input model
    //     pop.epp_disease_model_direct(hivpop, artpop)
    }
    if (pop.p.popadjust) { // match target pop
      pop.adjust_pop();
      if (MODEL!=0) {
        hivpop.adjust_pop(pop.adj_prob);
        if (i >= pop.p.tARTstart - 1)
          artpop.adjust_pop(pop.adj_prob);
      }
    }
    if (MODEL!=0) {
      if (i + pop.AGE_START <= pop.PROJ_YEARS - 1)
        pop.cal_prev_pregant(hivpop, artpop); // prevalence among pregnant women
      pop.save_prev_n_inc(); // save prevalence and incidence 15 to 49
    }
  }
  pop.finalize(hivpop, artpop);
  return pop.data_sexp;
}