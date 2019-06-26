#include "Classes.hpp"

void Model::initiate() {
  pop.initiate(p, s);
  pop.entrant_art_.resize(s.NG * s.pDS);
  if (!p.ad.scale_cd4_mort)
    hivpop.cd4_mort_ = p.nh.cd4_mort; // moved from scale_cd4_mort, do it once
  if (p.ic.eppmod != 1) // mvoed from update_rvec, do this once
    for (int i = 0; i < s.n_steps; ++i)
      pop.rvec[i] = p.ic.rvec[i];
}

void Model::update_views() {
  artpop.at_prev     = artpop.at_this;
  artpop.at_this    += artpop.N;
  hivpop.at_prev     = hivpop.at_this;
  hivpop.at_this    += hivpop.N;
  pop.at_prev        = pop.at_this;
  pop.at_this       += pop.N;
  artpop.at_prev_db  = artpop.at_this_db;
  artpop.at_this_db += artpop.N;
}

void Model::aging(Views& v) {
  pop.aging(v, s);
  pop.add_entrants(v, p, s);
  if (s.MODEL == 2)
    pop.sexual_debut(p, s);
  if (s.MODEL != 0) {
    pop.update_hiv_aging_prob(v, s);
    hivpop.aging(pop.hiv_aging_prob_, v, s);
    hivpop.add_entrants(pop.entrant_art_, v, p, s);
    if (s.MODEL == 2)
      hivpop.sexual_debut(v, p, s);
    if (s.year > s.tARTstart - 1) {
      artpop.aging(pop.hiv_aging_prob_, v, s);
      artpop.add_entrants(pop.entrant_art_, v, p, s);
      if (s.MODEL == 2)
        artpop.sexual_debut(v, p, s);
    }
  }
}

void Model::death(Views& v) {
  pop.deaths(v, p, s);
  if (s.MODEL != 0) {
    hivpop.deaths(pop.hiv_sx_prob, v, s);
    if (s.year > s.tARTstart - 1)
      artpop.deaths(pop.hiv_sx_prob, v, s);
  }
}

void Model::migration(Views& v) {
  pop.migration(v, p, s);
  if (s.MODEL != 0) {
    hivpop.migration(pop.hiv_mr_prob, v, s);
    if (s.year > s.tARTstart - 1)
      artpop.migration(pop.hiv_mr_prob, v, s);
  }
}

void Model::adjust_pop(Views& v) {
  if (p.dm.flag_popadjust) { // match target pop
    pop.adjust_pop(v, p, s);
    if (s.MODEL != 0) {
      hivpop.adjust_pop(pop.adj_prob, v, s);
      if (s.year >= s.tARTstart - 1)
        artpop.adjust_pop(pop.adj_prob, v, s);
    }
  }
}

void Model::infection_process(Views& v) {
  for (int time_step = 0; time_step < s.steps_per_year; ++time_step) {
    if (p.ic.eppmod != 2) { // != "directincid"
      pop.update_rvec(time_step, p, s);
      if (s.MIX)
        pop.infect_mix(time_step, v, p, s);
      else
        pop.infect_spec(hivpop, artpop, time_step, v, p, s);
      pop.update_infection(v, s);
      hivpop.update_infection(pop.infections_, p, s);
    }
    if (p.ad.scale_cd4_mort && s.year >= s.tARTstart - 1)
      hivpop.scale_cd4_mort(artpop, v, p, s);
    hivpop.grad_progress(v, p, s); // cd4 disease progression and mortality
    if (s.year >= s.tARTstart - 1)
      artpop.count_death(v, p, s);
    pop.remove_hiv_death(hivpop, artpop, v, p, s); // Remove hivdeaths from pop
    if (s.year >= s.tARTstart - 1) // ART initiation
      pop.epp_art_init(hivpop, artpop, time_step, v, p, s);
    hivpop.add_grad_to_pop(v, s);
  } // end time step
}

void Model::save_outputs(Views& v) {
  if (s.MODEL != 0) {
    if (s.year + s.AGE_START <= s.PROJ_YEARS - 1)
      pop.cal_prev_pregant(hivpop, artpop, v, p, s); // prevalence among pregnant women
    pop.save_prev_n_inc(v, s); // save prevalence and incidence 15 to 49
  }
}

void Model::run(int t) {
  s.year = t;
  update_views();
  Views v(pop_start, hiv_start, art_start, s, t);
  aging(v);
  death(v);
  migration(v);
  pop.update_fertile(v, p, s);
  if (s.MODEL != 0)  {
    infection_process(v);
    if (p.ic.eppmod == 2)
      pop.epp_disease_model_direct(hivpop, artpop, v, p, s);
  }
  adjust_pop(v);
  save_outputs(v);
}

extern "C" SEXP eppasmOOpp(SEXP fp) {
  StateSpace s(fp); // read only state-space
  Parameters p(fp); // read only parameters
  outputSEXP O(s);  // create all SEXP for output based on state-space
  Model model(O, s, p); // pass O for model to map to output's address
  model.initiate();
  for (int i = 1; i < s.SIM_YEARS; ++i)
    model.run(i);
  O.finalize(s);
  return O.pop;
}