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

// Natural deaths
// -----------------------------------------------------------------------------
void epp_death (popC& pop, hivC& hivpop, artC& artpop) {
  pop.deaths();
  if (pop.MODEL != 0) {
    hivpop.deaths(pop.hiv_sx_prob);
    if (pop.year > pop.p.tARTstart - 1)
      artpop.deaths(pop.hiv_sx_prob);
  }
}

// Migration at year i
// -----------------------------------------------------------------------------
void epp_migration (popC& pop, hivC& hivpop, artC& artpop) {
  pop.migration();
  if (pop.MODEL != 0) {
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
  if (pop.MODEL == 2)
    pop.sexual_debut();
  if (pop.MODEL != 0) {
    boost2D hiv_ag_prob = pop.hiv_aging_prob();
    dvec artYesNo = pop.entrant_art();
    hivpop.aging(hiv_ag_prob);
    hivpop.add_entrants(artYesNo);
    if (pop.MODEL == 2)
      hivpop.sexual_debut();
    if (pop.year > pop.p.tARTstart - 1) {
      artpop.aging(hiv_ag_prob);
      artpop.add_entrants(artYesNo);
      if (pop.MODEL == 2)
        artpop.sexual_debut();
    }
  }
}

// Disease model
// -----------------------------------------------------------------------------
void epp_disease_model (popC& pop, hivC& hivpop, artC& artpop) {
  for (int time_step = 0; time_step < pop.hiv_steps_per_year; ++time_step) {
    if (pop.p.eppmod != 2) { // != "directincid"
      pop.update_rvec(time_step);
      if (pop.MIX)
        pop.infect_mix(time_step);
      else
        pop.infect_spec(hivpop, artpop, time_step);
      pop.update_infection();
      hivpop.update_infection(pop.infections_);
    }
    hivpop.scale_cd4_mort(artpop);
    hivpop.grad_progress(); // cd4 disease progression and mortality
    if (pop.year >= pop.p.tARTstart - 1)
      artpop.count_death();
    pop.remove_hiv_death(hivpop, artpop); // Remove hivdeaths from pop
    if (pop.year >= pop.p.tARTstart - 1) // ART initiation
      pop.epp_art_init(hivpop, artpop, time_step);
    hivpop.add_grad_to_pop();
  } // end time step
}

// calculate, distribute eligible for ART, update grad, gradART
// -----------------------------------------------------------------------------
void popC::epp_art_init (hivC& hivpop, artC& artpop, int time_step) {
  artpop.grad_progress();
  artpop.art_dropout(hivpop); // pass hivpop to receive the drop out
  dvec eligible = hivpop.eligible_for_art();
  zeroing(art_elig_);
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        art_elig_[sex][agr][cd4] =
          hivpop.data[year][sex][agr][cd4] * eligible[cd4];
  if ( (p.pw_artelig[year] == 1) & (p.artcd4elig_idx[year] > 1) )
    update_preg(hivpop, artpop); // add pregnant?
  if (MODEL == 2) // add sexual inactive but eligible for treatment
    for (int sex = 0; sex < NG; sex++)
      for (int agr = 0; agr < hDB; agr++)
        for (int cd4 = 0; cd4 < hDS; cd4++)
          art_elig_[sex][agr][cd4] +=
            hivpop.data_db[year][sex][agr][cd4] * eligible[cd4];
  // calculate number to initiate ART and distribute
  dvec art_curr = artpop.current_on_art();
  dvec artnum_ii = art_initiate(art_curr, time_step);
  dvec art15plus_inits(NG);
  for (int sex = 0; sex < NG; ++sex) {
    double n_afford = artnum_ii[sex] - art_curr[sex];
    art15plus_inits[sex] = (n_afford > 0) ? n_afford : 0;
  }
  art_distribute(art15plus_inits);
  if (MODEL == 1) {
    for (int sex = 0; sex < NG; sex++)
      for (int agr = 0; agr < hAG; agr++)
        for (int cd4 = 0; cd4 < hDS; cd4++) {
          double x =
            hivpop.data[year][sex][agr][cd4] + DT * hivpop.grad[sex][agr][cd4];
          if (art_init_[sex][agr][cd4] > x) art_init_[sex][agr][cd4] = x;
        }
  }
  if (MODEL == 2) // split the number proportionally for active and idle pop
    hivpop.distribute_artinit(art_init_, artpop);
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        hivpop.grad[sex][agr][cd4] -= art_init_[sex][agr][cd4] / DT;
  artpop.grad_init(art_init_);
}