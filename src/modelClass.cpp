#include "Classes.hpp"

void Model::initiate() {
  pop.initiate(p, s);
  pop.entrant_art_.resize(s.NG * s.pDS);
  if (!p.scale_cd4_mort)
    hivpop.cd4_mort_ = p.cd4_mort; // moved from scale_cd4_mort, do it once
  if (p.eppmod != 1) // mvoed from update_rvec, do this once
    for (int i = 0; i < s.n_steps; ++i)
      pop.rvec[i] = p.rvec[i];
}

void Model::aging() {
  pop.aging(s);
  pop.add_entrants(p, s);
  if (s.MODEL == 2)
    pop.sexual_debut(p, s);
  if (s.MODEL != 0) {
    pop.update_hiv_aging_prob(s);
    hivpop.aging(pop.hiv_aging_prob_, s);
    hivpop.add_entrants(pop.entrant_art_, p, s);
    if (s.MODEL == 2)
      hivpop.sexual_debut(p, s);
    if (s.year > p.tARTstart - 1) {
      artpop.aging(pop.hiv_aging_prob_, s);
      artpop.add_entrants(pop.entrant_art_, p, s);
      if (s.MODEL == 2)
        artpop.sexual_debut(p, s);
    }
  }
}
void Model::death() {
  pop.deaths(p, s);
  if (s.MODEL != 0) {
    hivpop.deaths(pop.hiv_sx_prob, s);
    if (s.year > p.tARTstart - 1)
      artpop.deaths(pop.hiv_sx_prob, s);
  }
}

void Model::migration() {
  pop.migration(p, s);
  if (s.MODEL != 0) {
    hivpop.migration(pop.hiv_mr_prob, s);
    if (s.year > p.tARTstart - 1)
      artpop.migration(pop.hiv_mr_prob, s);
  }
}

void Model::adjust_pop() {
  if (p.popadjust) { // match target pop
    pop.adjust_pop(p, s);
    if (s.MODEL != 0) {
      hivpop.adjust_pop(pop.adj_prob, s);
      if (s.year >= p.tARTstart - 1)
        artpop.adjust_pop(pop.adj_prob, s);
    }
  }
}

void Model::infection_process() {
  for (int time_step = 0; time_step < s.hiv_steps_per_year; ++time_step) {
    if (p.eppmod != 2) { // != "directincid"
      pop.update_rvec(time_step, p, s);
      if (s.MIX)
        pop.infect_mix(time_step, p, s);
      else
        pop.infect_spec(hivpop, artpop, time_step, p, s);
      pop.update_infection(s);
      hivpop.update_infection(pop.infections_, p, s);
    }
    if (p.scale_cd4_mort && s.year >= p.tARTstart - 1)
      hivpop.scale_cd4_mort(artpop, p, s);
    hivpop.grad_progress(p, s); // cd4 disease progression and mortality
    if (s.year >= p.tARTstart - 1)
      artpop.count_death(p, s);
    pop.remove_hiv_death(hivpop, artpop, p, s); // Remove hivdeaths from pop
    if (s.year >= p.tARTstart - 1) // ART initiation
      pop.epp_art_init(hivpop, artpop, time_step, p, s);
    hivpop.add_grad_to_pop(s);
  } // end time step
}

void Model::save_outputs() {
  if (s.MODEL != 0) {
    if (s.year + s.AGE_START <= s.PROJ_YEARS - 1)
      pop.cal_prev_pregant(hivpop, artpop, p, s); // prevalence among pregnant women
    pop.save_prev_n_inc(s); // save prevalence and incidence 15 to 49
  }
}

void Model::run(int t) {
  s.year = t;
  aging();
  death();
  migration();
  pop.update_fertile(p, s);
  if (s.MODEL != 0)  {
    infection_process();
    if (p.eppmod == 2)
      pop.epp_disease_model_direct(hivpop, artpop, p, s);
  }
  adjust_pop();
  save_outputs();
}