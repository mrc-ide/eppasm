#include "Classes.hpp"

void artC::aging(const boost2D& ag_prob, Views& v, const StateSpace& s) {
  for (int i = 0; i < N; i++) {
    *(at_this + i) = *(at_prev + i);
    if (s.MODEL == 2)
      *(at_this_db + i) = *(at_prev_db + i);
  }
  double nARTup;
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG-1; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          nARTup = v.pre_art[sex][agr][cd4][dur] * ag_prob[sex][agr];
          v.now_art[sex][agr][cd4][dur]   -= nARTup;
          v.now_art[sex][agr+1][cd4][dur] += nARTup;
          if (s.MODEL == 2 && agr < s.hDB - 1) {
            nARTup = data_db[s.year-1][sex][agr][cd4][dur] * ag_prob[sex][agr];
            data_db[s.year][sex][agr][cd4][dur]   -= nARTup;
            data_db[s.year][sex][agr+1][cd4][dur] += nARTup;
          }
        }
}

void artC::add_entrants(const dvec& artYesNo, Views& v, const Parameters& p,
                        const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int cd4 = 0; cd4 < s.hDS; cd4++)
      for (int dur = 0; dur < s.hTS; dur++) {
        double add =
          p.ph.paedsurv_artcd4dist[s.year][sex][cd4][dur] * artYesNo[sex];
        if (s.MODEL == 1)
          v.now_art[sex][0][cd4][dur] += add;
        if (s.MODEL == 2) // add to virgin then debut
          data_db[s.year][sex][0][cd4][dur] += add;
      }
}

void artC::sexual_debut(Views& v, const Parameters& p, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hDB; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          double n_db = data_db[s.year][sex][agr][cd4][dur] * p.ic.db_rate[s.year][sex][agr];
          v.now_art[sex][agr][cd4][dur]       += n_db;
          data_db[s.year][sex][agr][cd4][dur] -= n_db;
        }
}

void artC::deaths(const boost2D& survival_pr, Views& v, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          v.now_art[sex][agr][cd4][dur] *= survival_pr[sex][agr];
          if (s.MODEL == 2 && agr < s.hDB)
            data_db[s.year][sex][agr][cd4][dur] *= survival_pr[sex][agr];
        }
}

void artC::migration(const boost2D& migration_pr, Views& v, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          v.now_art[sex][agr][cd4][dur] *= migration_pr[sex][agr];
          if (s.MODEL == 2 && agr < s.hDB)
            data_db[s.year][sex][agr][cd4][dur] *= migration_pr[sex][agr];
        }
}

void artC::grad_progress(Views& v, const StateSpace& s) {
  // int itemsPerCacheLine = cache_line_size() / sizeof(double);
  zeroing(gradART); // reset gradient
  if (s.MODEL == 2)
    zeroing(gradART_db); // reset gradient
  double art_up;
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
        for (int dur = 0; dur < s.hTS - 1; dur++) { // 2 x 9 x 7 x 2
          art_up = 2.0 * v.now_art[sex][agr][cd4][dur];
          gradART[sex][agr][cd4][dur]   -= (art_up + death_[sex][agr][cd4][dur]);
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

void artC::art_dropout(hivC& hivpop, Views& v,
                       const Parameters& p,
                       const StateSpace& s) {
  double n_dropout, p_dropout = p.ad.art_dropout[s.year];
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          n_dropout = v.now_art[sex][agr][cd4][dur] * p_dropout;
          hivpop.grad[sex][agr][cd4]  += n_dropout;
          gradART[sex][agr][cd4][dur] -= n_dropout;
          if (s.MODEL == 2 && agr < s.hDB) {
            n_dropout = data_db[s.year][sex][agr][cd4][dur] * p_dropout;
            hivpop.grad_db[sex][agr][cd4]  += n_dropout;
            gradART_db[sex][agr][cd4][dur] -= n_dropout;
          }
        }
}

void artC::update_current_on_art(Views& v, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++) {
    art_by_sex_[sex] = .0; // reset when call
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          art_by_sex_[sex] += (v.now_art[sex][agr][cd4][dur] + 
                               gradART[sex][agr][cd4][dur] * s.DT);
          if (s.MODEL == 2 && agr < s.hDB)  // add art from virgin pop
            art_by_sex_[sex] += (data_db[s.year][sex][agr][cd4][dur] + 
                                 gradART_db[sex][agr][cd4][dur] * s.DT);          
        }
  }
}

void artC::grad_init(const boost3D& artinit, Views& v, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
        gradART[sex][agr][cd4][0] += artinit[sex][agr][cd4] / s.DT;
        for (int dur = 0; dur < s.hTS; dur++)
          v.now_art[sex][agr][cd4][dur] += s.DT * gradART[sex][agr][cd4][dur];
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

void artC::adjust_pop(const boost2D& adj_prob, Views& v, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          v.now_art[sex][agr][cd4][dur] *= adj_prob[sex][agr];
          if (s.MODEL == 2 && agr < s.hDB)
            data_db[s.year][sex][agr][cd4][dur] *= adj_prob[sex][agr];
        }
}

void artC::count_death(Views& v, const Parameters& p, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          double x =
            p.nh.art_mort[sex][agr][cd4][dur] * p.nh.artmx_timerr[s.year][dur];
          death_[sex][agr][cd4][dur] = v.now_art[sex][agr][cd4][dur] * x;
          if (s.MODEL == 2 && agr < s.hDB)
            death_db_[sex][agr][cd4][dur] =
              data_db[s.year][sex][agr][cd4][dur] * x;
        }
}

boost2D artC::n_by_agr(Views& v, const Parameters& p, const StateSpace& s) {
  boost2D out(extents[s.NG][s.hAG]);
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          out[sex][agr] += v.now_art[sex][agr][cd4][dur];
          if (s.MODEL == 2)
            out[sex][agr] += data_db[s.year][sex][agr][cd4][dur];
        }
  return out;
}