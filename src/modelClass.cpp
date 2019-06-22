#include "Classes.hpp"

void Model::initiate() {
  pop.initiate(p);
  pop.entrant_art_.resize(pop.NG * pop.pDS);
  if (!p.scale_cd4_mort)
    hivpop.cd4_mort_ = p.cd4_mort; // moved from scale_cd4_mort, do it once
  if (p.eppmod != 1) // mvoed from update_rvec, do this once
    for (int i = 0; i < pop.n_steps; ++i)
      pop.rvec[i] = p.rvec[i];
}

void Model::aging() {
  pop.aging();
  pop.add_entrants(p);
  if (pop.MODEL == 2)
    pop.sexual_debut(p);
  if (pop.MODEL != 0) {
    pop.update_hiv_aging_prob();
    hivpop.aging(pop.hiv_aging_prob_);
    hivpop.add_entrants(pop.entrant_art_, p);
    if (pop.MODEL == 2)
      hivpop.sexual_debut(p);
    if (pop.year > p.tARTstart - 1) {
      artpop.aging(pop.hiv_aging_prob_);
      artpop.add_entrants(pop.entrant_art_, p);
      if (pop.MODEL == 2)
        artpop.sexual_debut(p);
    }
  }
}
void Model::death() {
  pop.deaths(p);
  if (pop.MODEL != 0) {
    hivpop.deaths(pop.hiv_sx_prob);
    if (pop.year > p.tARTstart - 1)
      artpop.deaths(pop.hiv_sx_prob);
  }
}

void Model::migration() {
  pop.migration(p);
  if (pop.MODEL != 0) {
    hivpop.migration(pop.hiv_mr_prob);
    if (pop.year > p.tARTstart - 1)
      artpop.migration(pop.hiv_mr_prob);
  }
}

void Model::adjust_pop() {
  if (p.popadjust) { // match target pop
    pop.adjust_pop(p);
    if (pop.MODEL != 0) {
      hivpop.adjust_pop(pop.adj_prob);
      if (pop.year >= p.tARTstart - 1)
        artpop.adjust_pop(pop.adj_prob);
    }
  }
}

void Model::infection_process() {
  for (int time_step = 0; time_step < pop.hiv_steps_per_year; ++time_step) {
    if (p.eppmod != 2) { // != "directincid"
      pop.update_rvec(time_step, p);
      if (pop.MIX)
        pop.infect_mix(time_step, p);
      else
        pop.infect_spec(hivpop, artpop, time_step, p);
      pop.update_infection();
      hivpop.update_infection(pop.infections_, p);
    }
    if (p.scale_cd4_mort && pop.year >= p.tARTstart - 1)
      hivpop.scale_cd4_mort(artpop, p);
    hivpop.grad_progress(p); // cd4 disease progression and mortality
    if (pop.year >= p.tARTstart - 1)
      artpop.count_death(p);
    pop.remove_hiv_death(hivpop, artpop, p); // Remove hivdeaths from pop
    if (pop.year >= p.tARTstart - 1) // ART initiation
      pop.epp_art_init(hivpop, artpop, time_step, p);
    hivpop.add_grad_to_pop();
  } // end time step
}

void Model::save_outputs() {
  if (pop.MODEL != 0) {
    if (pop.year + pop.AGE_START <= pop.PROJ_YEARS - 1)
      pop.cal_prev_pregant(hivpop, artpop, p); // prevalence among pregnant women
    pop.save_prev_n_inc(); // save prevalence and incidence 15 to 49
  }
}

void Model::run(int t) {
  pop.year = t; hivpop.year = t; artpop.year = t;
  aging();
  death();
  migration();
  pop.update_fertile(p);
  if (pop.MODEL != 0)  {
    infection_process();
    if (p.eppmod == 2)
      pop.epp_disease_model_direct(hivpop, artpop, p);
  }
  adjust_pop();
  save_outputs();
}