#include "Classes.hpp"

// calculate, distribute eligible for ART, update grad, gradART
// -----------------------------------------------------------------------------
void popC::epp_art_init (hivC& hivpop, artC& artpop, int time_step,
                         const Parameters& p, const StateSpace& s) {
  artpop.grad_progress(s);
  artpop.art_dropout(hivpop, p, s); // pass hivpop to receive the drop out
  dvec eligible = hivpop.eligible_for_art(p, s);
  zeroing(art_elig_);
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        art_elig_[sex][agr][cd4] =
          hivpop.data[s.year][sex][agr][cd4] * eligible[cd4];
  if ( (p.pw_artelig[s.year] == 1) && (p.artcd4elig_idx[s.year] > 1) )
    update_preg(hivpop, artpop, p, s); // add pregnant?
  if (s.MODEL == 2) // add sexual inactive but eligible for treatment
    for (int sex = 0; sex < s.NG; sex++)
      for (int agr = 0; agr < s.hDB; agr++)
        for (int cd4 = 0; cd4 < s.hDS; cd4++)
          art_elig_[sex][agr][cd4] +=
            hivpop.data_db[s.year][sex][agr][cd4] * eligible[cd4];
  // calculate number to initiate ART and distribute
  artpop.update_current_on_art(s);
  dvec artnum_ii = art_initiate(artpop.art_by_sex_, time_step, p, s);
  dvec art15plus_inits(s.NG);
  for (int sex = 0; sex < s.NG; ++sex) {
    double n_afford = artnum_ii[sex] - artpop.art_by_sex_[sex];
    art15plus_inits[sex] = (n_afford > 0) ? n_afford : 0;
  }
  art_distribute(art15plus_inits, p, s);
  if (s.MODEL == 1) {
    for (int sex = 0; sex < s.NG; sex++)
      for (int agr = 0; agr < s.hAG; agr++)
        for (int cd4 = 0; cd4 < s.hDS; cd4++) {
          double x =
            hivpop.data[s.year][sex][agr][cd4] + s.DT * hivpop.grad[sex][agr][cd4];
          if (art_init_[sex][agr][cd4] > x) art_init_[sex][agr][cd4] = x;
        }
  }
  if (s.MODEL == 2) // split the number proportionally for active and idle pop
    hivpop.distribute_artinit(art_init_, artpop, s);
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        hivpop.grad[sex][agr][cd4] -= art_init_[sex][agr][cd4] / s.DT;
  artpop.grad_init(art_init_, s);
}