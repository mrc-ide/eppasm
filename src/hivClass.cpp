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

void hivC::aging(mat ag_prob) {
  data(year) = data(year-1);
  cube nHup = sweepX23( data(year).cols(0, hAG-2), ag_prob.rows(0, hAG -2) );
  data(year).cols(0, hAG-2) -= nHup;
  data(year).cols(1, hAG-1) += nHup;
  if (MODEL==2) {
    data_db(year) = data_db(year-1);
    nHup = sweepX23( data_db(year).cols(0, hAG-2), ag_prob.rows(0, hAG -2) );
    data_db(year).cols(0, hAG-2) -= nHup;
    data_db(year).cols(1, hAG-1) += nHup;
  }
}

void hivC::add_entrants(vec artYesNo) { // see pop.entrant_art
  rowvec artNO = trans(artYesNo.tail(2));
  if (MODEL==1)
    data(year).col(0) += p.paedsurv_cd4dist.slice(year).each_row() % artNO;
  if (MODEL==2) // add to virgin then debut
    data_db(year).col(0) += p.paedsurv_cd4dist.slice(year).each_row() % artNO;
}

void hivC::sexual_debut() {
  cube debut_hiv = sweepX23(data_db(year).cols(0, pDB), p.db_pr);
     data(year).cols(0, pDB) += debut_hiv;
  data_db(year).cols(0, pDB) -= debut_hiv;
}

void hivC::deaths (mat survival_pr) {
  data(year) = sweepX23(data(year), survival_pr);
  if (MODEL==2)
    data_db(year) = sweepX23(data_db(year), survival_pr);
}

void hivC::migration (mat migration_pr) {
  data(year) = sweepX23(data(year), migration_pr);
  if (MODEL==2)
    data_db(year) = sweepX23(data_db(year), migration_pr);
}

void hivC::update_infection (mat new_infect) {
  grad.zeros(); // reset every time step
  grad += sweepX23(p.cd4_initdist, sumByAG(new_infect, ag_idx));
}

void hivC::grad_progress (cube mortality_rate) { // HIV gradient progress
  if (p.eppmod == 1)
    grad.zeros(); // reset every time step
  // remove cd4 stage progression (untreated)
  grad.rows(0, hDS-2) -= p.cd4_prog % data(year).rows(0, hDS-2);
  grad.rows(1, hDS-1) += p.cd4_prog % data(year).rows(0, hDS-2); // add 
  grad -= mortality_rate % data(year); // HIV mortality, untreated
  if (MODEL==2) {
    grad_db.zeros(); // reset, this's the 1st time grad_db is used
    grad_db.rows(0, hDS-2) -= p.cd4_prog % data_db(year).rows(0, hDS-2);
    grad_db.rows(1, hDS-1) += p.cd4_prog % data_db(year).rows(0, hDS-2);
    grad_db -= mortality_rate % data_db(year);
  }
}

vec hivC::eligible_for_art () {
  vec A = zeros(hDS); vec B = zeros(hDS); vec C = zeros(hDS);
  A.tail(hDS - (p.artcd4elig_idx(year) + 1) + 1).fill(1); // artcd4elig_idx
  B.tail(hDS - 2).fill(p.who34percelig);
  C.fill(p.specpop_percelig(year));
  return 1 - (1 - A) % (1 - B) % (1 - C);
}

cube hivC::distribute_artinit (cube artinit, artC& artpop) {
    cube debut_now = data_db(year) + DT * grad_db;
    cube all_hivpop = (data(year) + DT * grad) + debut_now;
    uvec greater = find(artinit > all_hivpop);
    artinit(greater) = all_hivpop(greater);
    cube pr_weight_db = debut_now / all_hivpop;
    cube artinit_db  = artinit % pr_weight_db;
    artinit -= artinit_db;
    grad_db -= artinit_db / DT;
    artpop.grad_db_init(artinit_db);
    return artinit;
}

void hivC::add_grad_to_pop () { // this is over-extraction
  data(year) += DT * grad;
  if (MODEL==2)
    data_db(year) += DT * grad_db;
}

void hivC::adjust_pop (mat adj_prob) {
  data(year) = sweepX23(data(year), adj_prob);
  if (MODEL==2)
    data_db(year) = sweepX23(data_db(year), adj_prob);
}

  // # avoid default set method conflict
  // set_data = function(FUN="+", x, DS=T, AG=T, NG=T, YEAR=NULL) {
  //     FUN <- match.fun(FUN)
  //     if (is.null(YEAR))
  //         data[DS, AG, NG] <<- FUN(data[DS, AG, NG], x)
  //     else 
  //         data[DS, AG, NG, YEAR] <<- FUN(data[DS, AG, NG, YEAR], x)
  // } 

  // get = function(YEAR=NULL, DS=T, AG=T, NG=T) {
  //     'get does not change anything, just return the data, use set to change'
  //     if (is.null(YEAR))
  //         data[DS, AG, NG, TRUE] # get all years
  //     else 
  //         data[DS, AG, NG, YEAR]
  // } 

  // sweep_sex = function(FUN="*", x, year) {
  //     data[,,,year] <<- sweep(data[,,,year], 2:3, x, FUN)
  // }