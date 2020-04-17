#include "Classes.hpp"

void hivC::aging(const boost2D& ag_prob, Views& v, const StateSpace& s) {
  for (int i = 0; i < N; i++)
    *(at_this + i) = *(at_prev + i);
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        if (s.MODEL == 2 && agr < s.hDB)
          data_db[s.year][sex][agr][cd4] = data_db[s.year-1][sex][agr][cd4];
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG - 1; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
        double nHup = v.pre_hiv[sex][agr][cd4] * ag_prob[sex][agr];
        v.now_hiv[sex][agr][cd4]   -= nHup;
        v.now_hiv[sex][agr+1][cd4] += nHup;
        if (s.MODEL == 2 && agr < s.hDB - 1) {
          nHup = data_db[s.year-1][sex][agr][cd4] * ag_prob[sex][agr];
          data_db[s.year][sex][agr][cd4]   -= nHup;
          data_db[s.year][sex][agr+1][cd4] += nHup;
        }
      }
}

void hivC::add_entrants(const dvec& artYesNo, Views& v,
                        const Parameters& p, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int cd4 = 0; cd4 < s.hDS; cd4++) {
      double add = p.ph.paedsurv_cd4dist[s.year][sex][cd4] * artYesNo[sex+2];
      if (s.MODEL == 1)
        v.now_hiv[sex][0][cd4] += add;
      if (s.MODEL == 2) // add to virgin then debut
        data_db[s.year][sex][0][cd4] += add;
    }
}

void hivC::sexual_debut(Views& v, const Parameters& p, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int adb = 0; adb < s.hDB; adb++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
        double n_db = data_db[s.year][sex][adb][cd4] * p.ic.db_rate[s.year][sex][adb];
        v.now_hiv[sex][adb][cd4]       += n_db;
        data_db[s.year][sex][adb][cd4] -= n_db;
      }
}

void hivC::deaths (const boost2D& survival_pr, Views& v, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
        v.now_hiv[sex][agr][cd4] *= survival_pr[sex][agr];
        if (s.MODEL == 2 && agr < s.hDB)
          data_db[s.year][sex][agr][cd4] *= survival_pr[sex][agr];
      }
}

void hivC::migration (const boost2D& migration_pr, Views& v, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
        v.now_hiv[sex][agr][cd4] *= migration_pr[sex][agr];
        if (s.MODEL == 2 && agr < s.hDB)
          data_db[s.year][sex][agr][cd4] *= migration_pr[sex][agr];
      }
}

void hivC::update_infection (const boost2D& new_infect,
                             const Parameters& p, const StateSpace& s) {
  zeroing(grad); // reset every time step
  zeroing(infect_by_agrp_);
  for (int sex = 0; sex < s.NG; ++sex) {
    int current_age_group = s.ag_[0]; // first age group
    for (int age = 0; age < s.pAG; ++age) {
      if ( s.ag_[age] != current_age_group)
        ++current_age_group;
      infect_by_agrp_[sex][current_age_group-1] += new_infect[sex][age];
    } // end age-groups
  }
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        grad[sex][agr][cd4] += 
          p.nh.cd4_initdist[sex][agr][cd4] * infect_by_agrp_[sex][agr];
}

void hivC::scale_cd4_mort(artC& artpop, Views& v,
                           const Parameters& p, const StateSpace& s) {
  double num, den = 0;
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
        num = v.now_hiv[sex][agr][cd4];
        if (s.MODEL==2)
          num += data_db[s.year][sex][agr][cd4];
        for (int dur = 0; dur < s.hTS; dur++) {
          den += v.now_art[sex][agr][cd4][dur];
          if (s.MODEL==2) 
            den += artpop.data_db[s.year][sex][agr][cd4][dur];
        }
        num = ((num + den) == 0.0) ? 1.0 : num / (num + den);
        cd4_mort_[sex][agr][cd4] = num * p.nh.cd4_mort[sex][agr][cd4];
        den = 0.0;
      }
}

void hivC::grad_progress(Views& v, const Parameters& p, const StateSpace& s) {
  if (p.ic.eppmod == 2)
    zeroing(grad); // reset every time step
  if (s.MODEL == 2)
    zeroing(grad_db); // reset, this's the 1st time grad_db is used
  // remove cd4 stage progression (untreated)
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++) {
      for (int cd4 = 0; cd4 < s.hDS - 1; cd4++) {
        double nHup = v.now_hiv[sex][agr][cd4] * p.nh.cd4_prog[sex][agr][cd4];
        death_[sex][agr][cd4]  =
          v.now_hiv[sex][agr][cd4] * cd4_mort_[sex][agr][cd4];
        grad[sex][agr][cd4]   -= (nHup + death_[sex][agr][cd4]);
        grad[sex][agr][cd4+1] += nHup;
        if (s.MODEL == 2 && agr < s.hDB) {
          nHup = data_db[s.year][sex][agr][cd4] * p.nh.cd4_prog[sex][agr][cd4];
          death_db_[sex][agr][cd4] =
            data_db[s.year][sex][agr][cd4] * cd4_mort_[sex][agr][cd4];
          grad_db[sex][agr][cd4]   -= (nHup + death_db_[sex][agr][cd4]);
          grad_db[sex][agr][cd4+1] += nHup;
        }
      }
      death_[sex][agr][s.hDS-1] = 
        v.now_hiv[sex][agr][s.hDS-1] * cd4_mort_[sex][agr][s.hDS-1];
      grad[sex][agr][s.hDS-1] -= death_[sex][agr][s.hDS-1];
      if (s.MODEL == 2 && agr < s.hDB) {
        death_db_[sex][agr][s.hDS-1] =
          data_db[s.year][sex][agr][s.hDS-1] * cd4_mort_[sex][agr][s.hDS-1];
        grad_db[sex][agr][s.hDS-1] -= death_db_[sex][agr][s.hDS-1];
      }
    }
}

void hivC::distribute_artinit (boost3D& artinit, artC& artpop,
                               Views& v, const StateSpace& s) {
    double debut_now, all_hivpop, pr_weight_db, n_artinit_db;
    boost3D artinit_db(extents[s.NG][s.hAG][s.hDS]);
    for (int sex = 0; sex < s.NG; sex++)
      for (int agr = 0; agr < s.hAG; agr++)
        for (int cd4 = 0; cd4 < s.hDS; cd4++) {
          debut_now =
            data_db[s.year][sex][agr][cd4] + s.DT * grad_db[sex][agr][cd4];
          all_hivpop =
            (v.now_hiv[sex][agr][cd4] + s.DT * grad[sex][agr][cd4]) + debut_now;
          if (artinit[sex][agr][cd4] > all_hivpop)
            artinit[sex][agr][cd4]   = all_hivpop;
          pr_weight_db               = debut_now / all_hivpop;
          n_artinit_db               = artinit[sex][agr][cd4] * pr_weight_db;
          artinit_db[sex][agr][cd4]  = n_artinit_db;
          artinit[sex][agr][cd4]    -= n_artinit_db;
          grad_db[sex][agr][cd4]    -= n_artinit_db / s.DT;
        }
    artpop.grad_db_init(artinit_db, s);
}

void hivC::add_grad_to_pop (Views& v, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
        v.now_hiv[sex][agr][cd4] += s.DT * grad[sex][agr][cd4];
        if (s.MODEL == 2 && agr < s.hDB)
          data_db[s.year][sex][agr][cd4] += s.DT * grad_db[sex][agr][cd4];
      }
}

void hivC::adjust_pop(const boost2D& adj_prob, Views& v, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
        v.now_hiv[sex][agr][cd4] *= adj_prob[sex][agr];
        if (s.MODEL == 2 && agr < s.hDB)
          data_db[s.year][sex][agr][cd4] *= adj_prob[sex][agr];
      }
}

boost2D hivC::n_by_agr(Views& v, const Parameters& p, const StateSpace& s) {
  boost2D out(extents[s.NG][s.hAG]); 
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
        out[sex][agr] += v.now_hiv[sex][agr][cd4];
        if (s.MODEL == 2)
          out[sex][agr] += data_db[s.year][sex][agr][cd4];
      }
  return out;
}