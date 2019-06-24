#include "Classes.hpp"

void artC::aging(const boost2D& ag_prob, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          data[s.year][sex][agr][cd4][dur] = data[s.year-1][sex][agr][cd4][dur];
          if (s.MODEL == 2 && agr < s.hDB)
            data_db[s.year][sex][agr][cd4][dur] =
              data_db[s.year-1][sex][agr][cd4][dur];
        }
  double nARTup;
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG-1; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          nARTup = data[s.year-1][sex][agr][cd4][dur] * ag_prob[sex][agr];
          data[s.year][sex][agr][cd4][dur]   -= nARTup;
          data[s.year][sex][agr+1][cd4][dur] += nARTup;
          if (s.MODEL == 2 && agr < s.hDB - 1) {
            nARTup = data_db[s.year-1][sex][agr][cd4][dur] * ag_prob[sex][agr];
            data_db[s.year][sex][agr][cd4][dur]   -= nARTup;
            data_db[s.year][sex][agr+1][cd4][dur] += nARTup;
          }
        }
}

void artC::add_entrants(const dvec& artYesNo, const Parameters& p, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int cd4 = 0; cd4 < s.hDS; cd4++)
      for (int dur = 0; dur < s.hTS; dur++) {
        double add = p.paedsurv_artcd4dist[s.year][sex][cd4][dur] * artYesNo[sex];
        if (s.MODEL == 1)
          data[s.year][sex][0][cd4][dur] += add;
        if (s.MODEL == 2) // add to virgin then debut
          data_db[s.year][sex][0][cd4][dur] += add;
      }
}

void artC::sexual_debut(const Parameters& p, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hDB; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          double n_db = data_db[s.year][sex][agr][cd4][dur] * p.db_pr[sex][agr];
          data[s.year][sex][agr][cd4][dur]    += n_db;
          data_db[s.year][sex][agr][cd4][dur] -= n_db;
        }
}

void artC::deaths(const boost2D& survival_pr, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          data[s.year][sex][agr][cd4][dur] *= survival_pr[sex][agr];
          if (s.MODEL == 2 && agr < s.hDB)
            data_db[s.year][sex][agr][cd4][dur] *= survival_pr[sex][agr];
        }
}

void artC::migration(const boost2D& migration_pr, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          data[s.year][sex][agr][cd4][dur] *= migration_pr[sex][agr];
          if (s.MODEL == 2 && agr < s.hDB)
            data_db[s.year][sex][agr][cd4][dur] *= migration_pr[sex][agr];
        }
}

void artC::grad_progress(const StateSpace& s) {
  zeroing(gradART); // reset gradient
  if (s.MODEL == 2)
    zeroing(gradART_db); // reset gradient
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
        for (int dur = 0; dur < s.hTS - 1; dur++) {
          double art_up = 2.0 * data[s.year][sex][agr][cd4][dur];
          gradART[sex][agr][cd4][dur] -= (art_up + death_[sex][agr][cd4][dur]);
          gradART[sex][agr][cd4][dur+1] += art_up;
          if (s.MODEL == 2 && agr < s.hDB) {
            art_up = 2.0 * data_db[s.year][sex][agr][cd4][dur];
            gradART_db[sex][agr][cd4][dur] -=
              (art_up + death_db_[sex][agr][cd4][dur]);
            gradART_db[sex][agr][cd4][dur+1] += art_up;
          }
        }
        gradART[sex][agr][cd4][s.hTS-1] -= death_[sex][agr][cd4][s.hTS-1];
        if (s.MODEL == 2 && agr < s.hDB)
          gradART_db[sex][agr][cd4][s.hTS-1] -= death_db_[sex][agr][cd4][s.hTS-1];
      }
}

void artC::art_dropout(hivC& hivpop, const Parameters& p, const StateSpace& s) {
  double n_dropout;
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          n_dropout = data[s.year][sex][agr][cd4][dur] * p.art_dropout[s.year];
          hivpop.grad[sex][agr][cd4]  += n_dropout;
          gradART[sex][agr][cd4][dur] -= n_dropout;
          if (s.MODEL == 2 && agr < s.hDB) {
            n_dropout = data_db[s.year][sex][agr][cd4][dur] * p.art_dropout[s.year];
            hivpop.grad_db[sex][agr][cd4]  += n_dropout;
            gradART_db[sex][agr][cd4][dur] -= n_dropout;
          }
        }
}

void artC::update_current_on_art(const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++) {
    art_by_sex_[sex] = .0; // reset when call
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          art_by_sex_[sex] += (data[s.year][sex][agr][cd4][dur] + 
                               gradART[sex][agr][cd4][dur] * s.DT);
          if (s.MODEL == 2 && agr < s.hDB)  // add art from virgin pop
            art_by_sex_[sex] += (data_db[s.year][sex][agr][cd4][dur] + 
                                 gradART_db[sex][agr][cd4][dur] * s.DT);          
        }
  }
}

void artC::grad_init(const boost3D& artinit, const StateSpace& s) { // 7x9x2
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
        gradART[sex][agr][cd4][0] += artinit[sex][agr][cd4] / s.DT;
        for (int dur = 0; dur < s.hTS; dur++)
          data[s.year][sex][agr][cd4][dur] += s.DT * gradART[sex][agr][cd4][dur];
      }
}

void artC::grad_db_init(const boost3D& artinit_db, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hDB; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
        gradART_db[sex][agr][cd4][0] += artinit_db[sex][agr][cd4] / s.DT;
        for (int dur = 0; dur < s.hTS; dur++)
          data_db[s.year][sex][agr][cd4][dur] += 
            s.DT * gradART_db[sex][agr][cd4][dur];
      }
}

void artC::adjust_pop(const boost2D& adj_prob, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          data[s.year][sex][agr][cd4][dur] *= adj_prob[sex][agr];
          if (s.MODEL == 2 && agr < s.hDB)
            data_db[s.year][sex][agr][cd4][dur] *= adj_prob[sex][agr];
        }
}

void artC::count_death(const Parameters& p, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          double x = p.art_mort[sex][agr][cd4][dur] * p.artmx_timerr[s.year][dur];
          death_[sex][agr][cd4][dur] = data[s.year][sex][agr][cd4][dur] * x;
          if (s.MODEL == 2 && agr < s.hDB)
            death_db_[sex][agr][cd4][dur] = data_db[s.year][sex][agr][cd4][dur] * x;
        }
}