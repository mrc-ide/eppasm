#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>

#include "model.hpp"

#define AGE_START 15

#define NG 2
#define pAG 66
#define pDS 2

#define pIDX_FERT 0
#define pAG_FERT 35
#define pIDX_15TO49 0
#define pAG_15TO49  35
#define pIDX_15PLUS 0
#define pAG_15PLUS  66

#define hAG 9
#define hDS 7
#define hTS 3

#define hIDX_FERT 0
#define hAG_FERT 8
#define hIDX_15TO49 0
#define hAG_15TO49  8
#define hIDX_15PLUS 0
#define hAG_15PLUS  9

#define hIDX_CD4_350 2

#define MALE 0
#define FEMALE 1

#define HIVN 0
#define HIVP 1

#define ART0MOS 0
#define ART6MOS 1
#define ART1YR 2

#define ART_STAGE_PROG_RATE 2.0 // HARD CODED: ART stage progression rate


#define EPP_RSPLINE 0
#define EPP_RTREND 1
#define EPP_DIRECTINCID 2  // annual direct incidence inputs (as Spectrum)

#define INCIDMOD_EPPSPEC 0
#define INCIDMOD_TRANSM 1

#define INCIDPOP_15TO49 0 // age range corresponding to incidence input
#define INCIDPOP_15PLUS 1

using namespace boost;

double calc_rtrend_rt(const multi_array_ref<double, 4> &pop,
                      double rtrend_tstab, const double *rtrend_beta, double rtrend_r0,
                      double projstep, double tsEpidemicStart, double DT, int t, int hts, double rveclast,
                      double *prevlast, double *prevcurr);

void calc_infections_eppspectrum2(const multi_array_ref<double, 4> &pop,
                                  const multi_array_ref<double, 4> &hivpop,
                                  const multi_array_ref<double, 5> &artpop,
                                  const double r_ts,
                                  const double relinfectART,
                                  const double iota,
                                  const double *incrr_sex,
                                  const const_multi_array_ref<double, 3> &incrr_age,
                                  const int t_ART_start,
                                  const double DT,
                                  const int t,
                                  const int hts,
                                  const int *hAG_START,
                                  const int *hAG_SPAN,
                                  double *prevcurr,
                                  double *incid15to49_ts,
                                  double infections_ts[NG][pAG]);

void calc_infections_eppspectrum2(const multi_array_ref<double, 4> &pop,
                                  const multi_array_ref<double, 4> &hivpop,
                                  const multi_array_ref<double, 5> &artpop,
                                  const double r_ts,
                                  const double relinfectART,
                                  const double iota,
                                  const double *incrr_sex,
                                  const const_multi_array_ref<double, 3> &incrr_age,
                                  const int t_ART_start,
                                  const double DT,
                                  const int t,
                                  const int hts,
                                  const int *hAG_START,
                                  const int *hAG_SPAN,
                                  double *prevcurr,
                                  double *incid15to49_ts,
                                  double infections_ts[NG][pAG])
{

  // sum population sizes
  double Xhivn_g[NG], Xhivn_incagerr[NG], Xhivp_noart = 0.0, Xart = 0.0;
  for(int g = 0; g < NG; g++){
    Xhivn_g[g] = 0.0;
    Xhivn_incagerr[g] = 0.0;
    for(int a = pIDX_15TO49; a < pIDX_15TO49+pAG_15TO49; a++){
      Xhivn_g[g] += pop[t][HIVN][g][a];
      Xhivn_incagerr[g] += incrr_age[t][g][a] * pop[t][HIVN][g][a];
    }

    for(int ha = hIDX_15TO49; ha < hIDX_15TO49+hAG_15TO49+1; ha++){

      // adjustment to first and last age group for partial year time step
      // calculation proportion of HIV population to include / exclude based on hivpop in single-year ages.
      double prop_include;
      if(ha == hIDX_15TO49){
        double hivp_ha = 0.0;
        for(int a = hAG_START[ha]; a < hAG_START[ha]+hAG_SPAN[ha]; a++)
          hivp_ha += pop[t][HIVP][g][a];
        prop_include = (hivp_ha > 0) ? 1.0 - pop[t][HIVP][g][hAG_START[ha]] / hivp_ha * (1.0 - DT*hts) : 1.0;
      } else if(ha == hIDX_15TO49+hAG_15TO49) {
        double hivp_ha = 0.0;
        for(int a = hAG_START[ha]; a < hAG_START[ha]+hAG_SPAN[ha]; a++)
          hivp_ha += pop[t][HIVP][g][a];
        prop_include = (hivp_ha > 0) ? pop[t][HIVP][g][hAG_START[ha]] / hivp_ha * (1.0 - DT*hts) : 1.0;
      } else
        prop_include = 1.0;

      for(int hm = 0; hm < hDS; hm++){
        Xhivp_noart += hivpop[t][g][ha][hm] * prop_include;
        if(t >= t_ART_start)
          for(int hu = 0; hu < hTS; hu++)
            Xart += artpop[t][g][ha][hm][hu] * prop_include;
      }
    }
  } // end loop over g
  double Xhivn = Xhivn_g[MALE] + Xhivn_g[FEMALE];

  // adjust HIV negative population for partial year time step
  for(int g = 0; g < NG; g++){
    Xhivn -= pop[t][HIVN][g][pIDX_15TO49] * (1.0 - DT*hts);
    Xhivn += pop[t][HIVN][g][pIDX_15TO49+pAG_15TO49] * (1.0 - DT*hts);
  }

  double Xtot = Xhivn + Xhivp_noart + Xart;
  *prevcurr = (Xhivp_noart + Xart) / Xtot;

  *incid15to49_ts = r_ts * (Xhivp_noart + relinfectART * Xart) / Xtot + iota;

  // incidence by sex
  double incid15to49_g[NG];
  incid15to49_g[MALE] = *incid15to49_ts * (Xhivn_g[MALE]+Xhivn_g[FEMALE]) / (Xhivn_g[MALE] + incrr_sex[t]*Xhivn_g[FEMALE]);
  incid15to49_g[FEMALE] = *incid15to49_ts * incrr_sex[t]*(Xhivn_g[MALE]+Xhivn_g[FEMALE]) / (Xhivn_g[MALE] + incrr_sex[t]*Xhivn_g[FEMALE]);

  // annualized infections by age and sex
  for(int g = 0; g < NG; g++)
    for(int a = 0; a < pAG; a++){
      infections_ts[g][a] = pop[t][HIVN][g][a] * incid15to49_g[g] * incrr_age[t][g][a] * Xhivn_g[g] / Xhivn_incagerr[g];
    }

  return;
}


void simulate_eppasm2(Model &md,
		      const SsDim &ss,
		      const DemogParam &demp,
		      const PaediatricHivParam &paedhp,
		      const NaturalHistoryParam &nhp,
		      const ArtData &artp,
		      const IncidenceParam &incp) {

  // md.initialise();

  for(int g = 0; g < NG; g++)
    for(int a = 0; a < pAG; a++){
      md.pop[0][HIVN][g][a] = demp.basepop[g][a];
      if((a >= pIDX_15TO49) & (a < pIDX_15TO49 + pAG_15TO49))
        md.hivn15to49 += demp.basepop[g][a];
    }


  ////////////////////////////////////
  ////  do population projection  ////
  ////////////////////////////////////

  for(int t = 1; t < ss.sim_years; t++){

    // md.ageing(t);

    // age the population one year
    for(int m = 0; m < pDS; m++)
      for(int g = 0; g < NG; g++){
        for(int a = 1; a < pAG; a++)
          md.pop[t][m][g][a] = md.pop[t-1][m][g][a-1];
        md.pop[t][m][g][pAG-1] += md.pop[t-1][m][g][pAG-1]; // open age group
      }

    // calculate proportion to age within each HIV age group
    double hiv_ag_prob[NG][hAG];
    for(int g = 0; g < NG; g++){
      int a = 0;
      for(int ha = 0; ha < (hAG-1); ha++){
        hiv_ag_prob[g][ha] = 0;
        for(int i = 0; i < ss.h_ag_span[ha]; i++){
          hiv_ag_prob[g][ha] += md.pop[t-1][HIVP][g][a];
          a++;
        }
        hiv_ag_prob[g][ha] = (hiv_ag_prob[g][ha] > 0) ? md.pop[t-1][HIVP][g][a-1] / hiv_ag_prob[g][ha] : 0;
      }
      hiv_ag_prob[g][hAG-1] = 0.0; // no one ages out of the open-ended age group
    }

    for(int g = 0; g < NG; g++){

      // youngest age group: only ageing out; ageing in entered below
      for(int hm = 0; hm < hDS; hm++){
        md.hivpop[t][g][0][hm] = (1-hiv_ag_prob[g][0]) * md.hivpop[t-1][g][0][hm];
        if(t > ss.t_art_start){
          for(int hu = 0; hu < hTS; hu++)
            md.artpop[t][g][0][hm][hu] = (1-hiv_ag_prob[g][0]) * md.artpop[t-1][g][0][hm][hu];
        }
      }

      for(int ha = 1; ha < hAG; ha++)
        for(int hm = 0; hm < hDS; hm++){
          md.hivpop[t][g][ha][hm] = (1-hiv_ag_prob[g][ha]) * md.hivpop[t-1][g][ha][hm] + hiv_ag_prob[g][ha-1] * md.hivpop[t-1][g][ha-1][hm];
          if(t > ss.t_art_start)
            for(int hu = 0; hu < hTS; hu++)
              md.artpop[t][g][ha][hm][hu] = (1-hiv_ag_prob[g][ha]) * md.artpop[t-1][g][ha][hm][hu] + hiv_ag_prob[g][ha-1] * md.artpop[t-1][g][ha-1][hm][hu];
        }
    }


    // add lagged births to youngest age group
    // md.pop_entrants(t);

    for(int g = 0; g < NG; g++){

      double paedsurv_g;
      double entrant_prev;

      if(paedhp.use_entrantprev)
        entrant_prev = paedhp.entrantprev[t][g];
      else
        entrant_prev = md.pregprevlag[t-1] * paedhp.verttrans_lag[t-1] * paedhp.paedsurv_lag[t-1];

      if(demp.flag_popadjust){
        md.pop[t][HIVN][g][0] = demp.entrantpop[t-1][g] * (1.0-entrant_prev);
        paedsurv_g = demp.entrantpop[t-1][g] * entrant_prev;
      } else {
        md.pop[t][HIVN][g][0] = md.birthslag[t-1][g] * demp.cumsurv[t-1][g] * (1.0-entrant_prev / paedhp.paedsurv_lag[t-1]) + demp.cumnetmigr[t-1][g] * (1.0-md.pregprevlag[t-1] * paedhp.netmig_hivprob);
        paedsurv_g = md.birthslag[t-1][g] * demp.cumsurv[t-1][g] * entrant_prev + demp.cumnetmigr[t-1][g] * entrant_prev;
      }

      md.pop[t][HIVP][g][0] = paedsurv_g;

      md.entrantprev[t] = (md.pop[t][HIVP][MALE][0] + md.pop[t][HIVP][FEMALE][0]) / (md.pop[t][HIVN][MALE][0] + md.pop[t][HIVN][FEMALE][0] + md.pop[t][HIVP][MALE][0] + md.pop[t][HIVP][FEMALE][0]);

      for(int hm = 0; hm < hDS; hm++){
        md.hivpop[t][g][0][hm] += paedsurv_g * paedhp.paedsurv_cd4dist[t][g][hm] * (1.0 - paedhp.entrantartcov[t][g]);
        if(t > ss.t_art_start){
          for(int hu = 0; hu < hTS; hu++){
            md.artpop[t][g][0][hm][hu] += paedsurv_g * paedhp.paedsurv_artcd4dist[t][g][hm][hu] * paedhp.entrantartcov[t][g];
          }
        }
      }
    }

    // non-HIV mortality and netmigration
    // md.death_and_migration(t);

    for(int g = 0; g < NG; g++){
      int a = 0;
      for(int ha = 0; ha < hAG; ha++){
        double deathsmig_ha = 0, hivpop_ha = 0;
        for(int i = 0; i < ss.h_ag_span[ha]; i++){

          hivpop_ha += md.pop[t][HIVP][g][a];

          // non-HIV mortality
          double qx = 1.0 - demp.Sx[t][g][a];
          double ndeaths_a = md.pop[t][HIVN][g][a] * qx;
          md.pop[t][HIVN][g][a] -= ndeaths_a; // survival HIV- population
          double hdeaths_a = md.pop[t][HIVP][g][a] * qx;
          deathsmig_ha -= hdeaths_a;
          md.pop[t][HIVP][g][a] -= hdeaths_a;   // survival HIV+ population
          md.natdeaths[t][g][a] = ndeaths_a + hdeaths_a;

          // net migration
          double migrate_a = demp.netmigr[t][g][a] * (1+demp.Sx[t][g][a])/2.0 / (md.pop[t][HIVN][g][a] + md.pop[t][HIVP][g][a]);
          md.pop[t][HIVN][g][a] *= 1+migrate_a;
          double hmig_a = migrate_a * md.pop[t][HIVP][g][a];
          deathsmig_ha += hmig_a;
          md.pop[t][HIVP][g][a] += hmig_a;

          a++;
        }

        // migration and deaths for hivpop
        double deathmigrate_ha = hivpop_ha > 0 ? deathsmig_ha / hivpop_ha : 0.0;
        for(int hm = 0; hm < hDS; hm++){
          md.hivpop[t][g][ha][hm] *= 1+deathmigrate_ha;
          if(t > ss.t_art_start)
            for(int hu = 0; hu < hTS; hu++)
              md.artpop[t][g][ha][hm][hu] *= 1+deathmigrate_ha;
        } // loop over hm
      } // loop over ha
    } // loop over g


    // md.fertility(t);

    // fertility
    md.births = 0.0;
    std::fill(md.births_by_ha.begin(), md.births_by_ha.end(), 0);
    for(int m = 0; m < pDS; m++){
      int a = pIDX_FERT;
      for(int ha = hIDX_FERT; ha < hAG_FERT; ha++){
        for(int i = 0; i < ss.h_ag_span[ha]; i++){
          md.births_by_ha[ha-hIDX_FERT] += (md.pop[t-1][m][FEMALE][a] + md.pop[t][m][FEMALE][a])/2 * demp.asfr[t][a];
          a++;
        }
      }
    }
    for(int ha = hIDX_FERT; ha < hAG_FERT; ha++)
      md.births += md.births_by_ha[ha-hIDX_FERT];

    if(t + AGE_START < ss.proj_years)
      for(int g = 0; g < NG; g++)
        md.birthslag[t + AGE_START-1][g] = demp.srb[t][g] * md.births;

    
    ////////////////////////////////
    ////  HIV model simulation  ////
    ////////////////////////////////

    int cd4elig_idx = artp.artcd4elig_idx[t] - 1; // -1 for 0-based indexing vs. 1-based in R
    int anyelig_idx = (artp.special_pop_percelig[t] > 0 | artp.preg_women_artelig[t] > 0) ? 0 : (artp.who34percelig > 0) ? hIDX_CD4_350 : cd4elig_idx;
    md.everARTelig_idx = anyelig_idx < md.everARTelig_idx ? anyelig_idx : md.everARTelig_idx;

    for(int hts = 0; hts < ss.hiv_steps_per_year; hts++){

      int ts = (t-1)*ss.hiv_steps_per_year + hts;

      double hivdeaths_ha[NG][hAG];
      memset(hivdeaths_ha, 0, sizeof(double)*NG*hAG);

      // untreated population

      // disease progression and mortality
      double grad[NG][hAG][hDS];
      for(int g = 0; g < NG; g++)
        for(int ha = 0; ha < hAG; ha++){
          for(int hm = 0; hm < hDS; hm++){

            grad[g][ha][hm] = 0.0;

            double cd4mx_scale = 1.0;
            if(artp.scale_cd4_mort & (t >= ss.t_art_start) & (hm >= md.everARTelig_idx)){
              double artpop_hahm = 0.0;
              for(int hu = 0; hu < hTS; hu++)
                artpop_hahm += md.artpop[t][g][ha][hm][hu];
              cd4mx_scale = md.hivpop[t][g][ha][hm] / (md.hivpop[t][g][ha][hm] + artpop_hahm);
            }

            double deaths = cd4mx_scale * nhp.cd4_mort[g][ha][hm] * md.hivpop[t][g][ha][hm];
            hivdeaths_ha[g][ha] += ss.dt * deaths;
            grad[g][ha][hm] = -deaths;
          }
          for(int hm = 1; hm < hDS; hm++){
            grad[g][ha][hm-1] -= nhp.cd4_prog[g][ha][hm-1] * md.hivpop[t][g][ha][hm-1];
            grad[g][ha][hm] += nhp.cd4_prog[g][ha][hm-1] * md.hivpop[t][g][ha][hm-1];
          }
        }

      if(incp.eppmod != EPP_DIRECTINCID){
        // incidence

        // calculate r(t)
        if(incp.eppmod == EPP_RSPLINE)
          md.rvec_ts[ts] = incp.rspline_rvec[ts];
        else
          md.rvec_ts[ts] = calc_rtrend_rt(md.pop, incp.rtrend_tstab, incp.rtrend_beta, incp.rtrend_r0,
                                          ss.proj_start + 0.5 + ts * 1.0 / ss.hiv_steps_per_year,
                                          ss.proj_start + 0.5 + incp.ts_epidemic_start * 1.0 / ss.hiv_steps_per_year,
                                          ss.dt, t, hts,
                                          md.rvec_ts[ts-1], &md.prevlast, &md.prevcurr);

        // calculate new infections by sex and age
        double infections_ts[NG][pAG];
        calc_infections_eppspectrum2(md.pop, md.hivpop, md.artpop,
                                     md.rvec_ts[ts], incp.rel_infect_art, (ts == incp.ts_epidemic_start) ? incp.iota : 0.0,
                                     incp.incrr_sex, incp.incrr_age,
                                     ss.t_art_start, ss.dt, t, hts, ss.h_ag_start.data(), ss.h_ag_span,
                                     &md.prevcurr, &md.incid15to49_ts[ts], infections_ts);

        md.prev15to49_ts[ts] = md.prevcurr;

        // add new infections to HIV population
        for(int g = 0; g < NG; g++){
          int a = 0;
          for(int ha = 0; ha < hAG; ha++){
            double infections_a, infections_ha = 0.0;
            for(int i = 0; i < ss.h_ag_span[ha]; i++){
              infections_ha += infections_a = infections_ts[g][a];
              md.infections[t][g][a] += ss.dt*infections_a;
              md.pop[t][HIVN][g][a] -= ss.dt*infections_a;
              md.pop[t][HIVP][g][a] += ss.dt*infections_a;
              a++;
            }
            if(ha < hIDX_15TO49+hAG_15TO49 )
              md.incid15to49[t] += ss.dt*infections_ha;

            // add infections to grad hivpop
            for(int hm = 0; hm < hDS; hm++)
              grad[g][ha][hm] += infections_ha * nhp.cd4_initdist[g][ha][hm];
          }
        }
      }

      // ART progression, mortality, and initiation
      if(t >= ss.t_art_start){

        double gradART[NG][hAG][hDS][hTS];

        // progression and mortality
        for(int g = 0; g < NG; g++)
          for(int ha = 0; ha < hAG; ha++)
            for(int hm = md.everARTelig_idx; hm < hDS; hm++){

              for(int hu = 0; hu < hTS; hu++){
                double deaths = nhp.art_mort[g][ha][hm][hu] * nhp.artmx_timerr[t][hu] * md.artpop[t][g][ha][hm][hu];
                hivdeaths_ha[g][ha] += ss.dt*deaths;
                gradART[g][ha][hm][hu] = -deaths;
              }

              gradART[g][ha][hm][ART0MOS] += -ART_STAGE_PROG_RATE * md.artpop[t][g][ha][hm][ART0MOS];
              gradART[g][ha][hm][ART6MOS] += ART_STAGE_PROG_RATE * md.artpop[t][g][ha][hm][ART0MOS] - ART_STAGE_PROG_RATE * md.artpop[t][g][ha][hm][ART6MOS];
              gradART[g][ha][hm][ART1YR] += ART_STAGE_PROG_RATE * md.artpop[t][g][ha][hm][ART6MOS];

              // ART dropout
              if(artp.art_dropout[t] > 0)
                for(int hu = 0; hu < hTS; hu++){
                  grad[g][ha][hm] += artp.art_dropout[t] * md.artpop[t][g][ha][hm][hu];
                  gradART[g][ha][hm][hu] -= artp.art_dropout[t] * md.artpop[t][g][ha][hm][hu];
                }

            }


        // ART initiation
        for(int g = 0; g < NG; g++){

          double artelig_hahm[hAG_15PLUS][hDS], Xart_15plus = 0.0, Xartelig_15plus = 0.0, expect_mort_artelig15plus = 0.0;
          for(int ha = hIDX_15PLUS; ha < hAG; ha++){
            for(int hm = md.everARTelig_idx; hm < hDS; hm++){
              if(hm >= anyelig_idx){
                double prop_elig = (hm >= cd4elig_idx) ? 1.0 : (hm >= hIDX_CD4_350) ? 1.0 - (1.0-artp.special_pop_percelig[t])*(1.0-artp.who34percelig) : artp.special_pop_percelig[t];
                Xartelig_15plus += artelig_hahm[ha-hIDX_15PLUS][hm] = prop_elig * md.hivpop[t][g][ha][hm] ;
                expect_mort_artelig15plus += nhp.cd4_mort[g][ha][hm] * artelig_hahm[ha-hIDX_15PLUS][hm];
              }
              for(int hu = 0; hu < hTS; hu++)
                Xart_15plus += md.artpop[t][g][ha][hm][hu] + ss.dt * gradART[g][ha][hm][hu];
            }

            // if artp.preg_women_artelig, add pregnant women to artelig_hahm population
            if(g == FEMALE & artp.preg_women_artelig[t] > 0 & ha < hAG_FERT){
              double frr_pop_ha = 0;
              for(int a =  ss.h_ag_start[ha]; a < ss.h_ag_start[ha]+ss.h_ag_span[ha]; a++)
                frr_pop_ha += md.pop[t][HIVN][g][a]; // add HIV- population
              for(int hm = 0; hm < hDS; hm++){
                frr_pop_ha += nhp.frr_cd4[t][ha-hIDX_FERT][hm] * md.hivpop[t][g][ha][hm];
                for(int hu = 0; hu < hTS; hu++)
                  frr_pop_ha += nhp.frr_art[t][ha-hIDX_FERT][hm][hu] * md.artpop[t][g][ha][hm][hu];
              }
              for(int hm = anyelig_idx; hm < cd4elig_idx; hm++){
                double pw_elig_hahm = md.births_by_ha[ha-hIDX_FERT] * nhp.frr_cd4[t][ha-hIDX_FERT][hm] * md.hivpop[t][g][ha][hm] / frr_pop_ha;
                artelig_hahm[ha-hIDX_15PLUS][hm] += pw_elig_hahm;
                Xartelig_15plus += pw_elig_hahm;
                expect_mort_artelig15plus += nhp.cd4_mort[g][ha][hm] * pw_elig_hahm;
              }
            }
          } // loop over ha

            // calculate number on ART at end of ts, based on number or percent
          double artnum_hts = 0.0;
          if(ss.dt*(hts+1) < 0.5){
            if(!artp.art15plus_isperc[t-2][g] & !artp.art15plus_isperc[t-1][g]){ // both numbers
              artnum_hts = (0.5-ss.dt*(hts+1))*artp.artnum15plus[t-2][g] + (ss.dt*(hts+1)+0.5)*artp.artnum15plus[t-1][g];
            } else if(artp.art15plus_isperc[t-2][g] & artp.art15plus_isperc[t-1][g]){ // both percentages
              double artcov_hts = (0.5-ss.dt*(hts+1))*artp.artnum15plus[t-2][g] + (ss.dt*(hts+1)+0.5)*artp.artnum15plus[t-1][g];
              artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
            } else if((!artp.art15plus_isperc[t-2][g]) & artp.art15plus_isperc[t-1][g]){ // transition from number to percentage
              double curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus);
              double artcov_hts = curr_coverage + (artp.artnum15plus[t-1][g] - curr_coverage) * ss.dt / (0.5-ss.dt*hts);
              artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
            }
          } else {
            if(!artp.art15plus_isperc[t-1][g] & !artp.art15plus_isperc[t][g]){ // both numbers
              artnum_hts = (1.5-ss.dt*(hts+1))*artp.artnum15plus[t-1][g] + (ss.dt*(hts+1)-0.5)*artp.artnum15plus[t][g];
            } else if(artp.art15plus_isperc[t-1][g] & artp.art15plus_isperc[t][g]){ // both percentages
              double artcov_hts = (1.5-ss.dt*(hts+1))*artp.artnum15plus[t-1][g] + (ss.dt*(hts+1)-0.5)*artp.artnum15plus[t][g];
              artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
            } else if((!artp.art15plus_isperc[t-1][g]) & artp.art15plus_isperc[t][g]){ // transition from number to percentage
              double curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus);
              double artcov_hts = curr_coverage + (artp.artnum15plus[t][g] - curr_coverage) * ss.dt / (1.5-ss.dt*hts);
              artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
            }
          }

          double artinit_hts = artnum_hts > Xart_15plus ? artnum_hts - Xart_15plus : 0;

          // median CD4 at initiation inputs
          if(artp.med_cd4init_input[t]){

            const int CD4_LOW_LIM[hDS] = {500, 350, 250, 200, 100, 50, 0};
            const int CD4_UPP_LIM[hDS] = {1000, 500, 350, 250, 200, 100, 50};

            int medcd4_idx = artp.med_cd4init_cat[t] - 1; // -1 for 0-based indexing vs. 1-based in R
            double medcat_propbelow = (artp.median_cd4init[t] - CD4_LOW_LIM[medcd4_idx]) / (CD4_UPP_LIM[medcd4_idx] - CD4_LOW_LIM[medcd4_idx]);

            double elig_below = 0.0, elig_above = 0.0;
            for(int ha = hIDX_15PLUS; ha < hAG; ha++){
              for(int hm = anyelig_idx; hm < medcd4_idx; hm++)
                elig_above += artelig_hahm[ha-hIDX_15PLUS][hm];
              elig_above += (1.0 - medcat_propbelow) * artelig_hahm[ha-hIDX_15PLUS][medcd4_idx];
              elig_below += medcat_propbelow * artelig_hahm[ha-hIDX_15PLUS][medcd4_idx];
              for(int hm = medcd4_idx+1; hm < hDS; hm++)
                elig_below += artelig_hahm[ha-hIDX_15PLUS][hm];
            }

            double initprob_below = (elig_below > artinit_hts * 0.5) ? artinit_hts * 0.5 / elig_below : 1.0;
            double initprob_above = (elig_above > artinit_hts * 0.5) ? artinit_hts * 0.5 / elig_above : 1.0;
            double initprob_medcat = initprob_below * medcat_propbelow + initprob_above * (1.0-medcat_propbelow);

            for(int ha = hIDX_15PLUS; ha < hAG; ha++)
              for(int hm = anyelig_idx; hm < hDS; hm++){
                double artinit_hahm;
                if(hm < medcd4_idx)
                  artinit_hahm = artelig_hahm[ha-hIDX_15PLUS][hm] * initprob_above;
                else if(hm == medcd4_idx)
                  artinit_hahm = artelig_hahm[ha-hIDX_15PLUS][hm] * initprob_medcat;
                if(hm > medcd4_idx)
                  artinit_hahm = artelig_hahm[ha-hIDX_15PLUS][hm] * initprob_below;
                if(artinit_hahm > md.hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm])
                  artinit_hahm = md.hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm];
                grad[g][ha][hm] -= artinit_hahm / ss.dt;
                gradART[g][ha][hm][ART0MOS] += artinit_hahm / ss.dt;
              }

          } else if(artp.art_alloc_method == 4) {  // lowest CD4 first

            for(int hm = hDS-1; hm >= anyelig_idx; hm--){
              double artelig_hm = 0;
              for(int ha = hIDX_15PLUS; ha < hAG; ha++)
                artelig_hm += artelig_hahm[ha-hIDX_15PLUS][hm];
              double init_prop = (artelig_hm == 0 | artinit_hts > artelig_hm) ? 1.0 : artinit_hts / artelig_hm;

              for(int ha = hIDX_15PLUS; ha < hAG; ha++){
                double artinit_hahm = init_prop * artelig_hahm[ha-hIDX_15PLUS][hm];

                if(artinit_hahm > md.hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm])
                  artinit_hahm = md.hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm];

                grad[g][ha][hm] -= artinit_hahm / ss.dt;
                gradART[g][ha][hm][ART0MOS] += artinit_hahm / ss.dt;
              }
              if(init_prop < 1.0)
                break;
              artinit_hts -= init_prop * artelig_hm;
            }

          } else { // Use mixture of eligibility and expected mortality for initiation distribution

            for(int ha = hIDX_15PLUS; ha < hAG; ha++)
              for(int hm = anyelig_idx; hm < hDS; hm++){
                double artinit_hahm = artinit_hts * artelig_hahm[ha-hIDX_15PLUS][hm] * ((1.0 - artp.art_alloc_mxweight)/Xartelig_15plus + artp.art_alloc_mxweight * nhp.cd4_mort[g][ha][hm] / expect_mort_artelig15plus);
                if(artinit_hahm > artelig_hahm[ha-hIDX_15PLUS][hm])
                  artinit_hahm = artelig_hahm[ha-hIDX_15PLUS][hm];
                if(artinit_hahm > md.hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm])
                  artinit_hahm = md.hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm];
                grad[g][ha][hm] -= artinit_hahm / ss.dt;
                gradART[g][ha][hm][ART0MOS] += artinit_hahm / ss.dt;
              }
          }
        }

        for(int g = 0; g < NG; g++)
          for(int ha = 0; ha < hAG; ha++)
            for(int hm = md.everARTelig_idx; hm < hDS; hm++)
              for(int hu = 0; hu < hTS; hu++)
                md.artpop[t][g][ha][hm][hu] += ss.dt*gradART[g][ha][hm][hu];

      } // if(t >= ss.t_art_start)

      for(int g = 0; g < NG; g++)
        for(int ha = 0; ha < hAG; ha++)
          for(int hm = 0; hm < hDS; hm++)
            md.hivpop[t][g][ha][hm] += ss.dt*grad[g][ha][hm];


      // remove hivdeaths from pop
      for(int g = 0; g < NG; g++){

        // sum HIV+ population size in each hivpop age group
        double hivpop_ha[hAG];
        int a = 0;
        for(int ha = 0; ha < hAG; ha++){
          hivpop_ha[ha] = 0.0;
          for(int i = 0; i < ss.h_ag_span[ha]; i++){
            hivpop_ha[ha] += md.pop[t][HIVP][g][a];
            a++;
          }
        }

        // remove hivdeaths proportionally to age-distribution within each age group
        a = 0;
        for(int ha = 0; ha < hAG; ha++){
          if(hivpop_ha[ha] > 0){
            double hivqx_ha = hivdeaths_ha[g][ha] / hivpop_ha[ha];
            for(int i = 0; i < ss.h_ag_span[ha]; i++){
              md.hivdeaths[t][g][a] += md.pop[t][HIVP][g][a] * hivqx_ha;
              md.pop[t][HIVP][g][a] *= (1.0-hivqx_ha);
              a++;
            }
          } else {
            a += ss.h_ag_span[ha];
          }  // end if(pop_ha[ha] > 0)
        }
      }



    } // loop HIVSTEPS_PER_YEAR



    if(incp.eppmod == EPP_DIRECTINCID){
      // Calculating new infections once per year (like Spectrum)

      double Xhivp = 0.0, Xhivn[NG], Xhivn_incagerr[NG];

      for(int g = 0; g < NG; g++){
        Xhivn[g] = 0.0;
        Xhivn_incagerr[g] = 0.0;
        for(int a = incp.i_p_ag_incidpop; a < incp.i_p_ag_incidpop + incp.r_p_ag_incidpop; a++){
          // for(int a = pIDX_INCIDPOP; a < pIDX_INCIDPOP+pAG_INCIDPOP; a++){
          Xhivp += md.pop[t-1][HIVP][g][a];
          Xhivn[g] += md.pop[t-1][HIVN][g][a];
          Xhivn_incagerr[g] += incp.incrr_age[t][g][a] * md.pop[t-1][HIVN][g][a];
        }
      }
      // double prev_i = Xhivp / (Xhivn[MALE] + Xhivn[FEMALE] + Xhivp);
      // double incrate15to49_i = (prev15to49[t] - prev_i)/(1.0 - prev_i);
      double incrate_i = incp.incidinput[t];
      double incrate_g[NG];
      incrate_g[MALE] = incrate_i * (Xhivn[MALE]+Xhivn[FEMALE]) / (Xhivn[MALE] + incp.incrr_sex[t]*Xhivn[FEMALE]);
      incrate_g[FEMALE] = incrate_i * incp.incrr_sex[t]*(Xhivn[MALE]+Xhivn[FEMALE]) / (Xhivn[MALE] + incp.incrr_sex[t]*Xhivn[FEMALE]);

      for(int g = 0; g < NG; g++){
        int a = 0;
        for(int ha = 0; ha < hAG; ha++){
          double infections_a, infections_ha = 0.0;
          for(int i = 0; i < ss.h_ag_span[ha]; i++){
            infections_ha += infections_a = md.pop[t-1][HIVN][g][a] * incrate_g[g] * incp.incrr_age[t][g][a] * Xhivn[g] / Xhivn_incagerr[g];
            md.infections[t][g][a] += infections_a;
            md.pop[t][HIVN][g][a] -= infections_a;
            md.pop[t][HIVP][g][a] += infections_a;
            a++;
          }
          if(ha < hIDX_15TO49+hAG_15TO49)
            md.incid15to49[t] += infections_ha;

          // add infections to hivpop
          for(int hm = 0; hm < hDS; hm++)
            md.hivpop[t][g][ha][hm] += infections_ha * nhp.cd4_initdist[g][ha][hm];
        }
      }
    }


    // md.adjust_population(t);

    // adjust population to match target population
    if(demp.flag_popadjust){
      for(int g = 0; g < NG; g++){
        int a = 0;
        for(int ha = 0; ha < hAG; ha++){
          double popadj_ha = 0, hivpop_ha = 0;
          for(int i = 0; i < ss.h_ag_span[ha]; i++){

            hivpop_ha += md.pop[t][HIVP][g][a];

            double popadjrate_a = md.popadjust[t][g][a] = demp.targetpop[t][g][a] / (md.pop[t][HIVN][g][a] + md.pop[t][HIVP][g][a]);
            md.pop[t][HIVN][g][a] *= popadjrate_a;
            double hpopadj_a = (popadjrate_a-1.0) * md.pop[t][HIVP][g][a];
            popadj_ha += hpopadj_a;
            md.pop[t][HIVP][g][a] += hpopadj_a;
            a++;
          }

          // population adjustment for hivpop
          double popadjrate_ha = hivpop_ha > 0 ? popadj_ha / hivpop_ha : 0.0;
          for(int hm = 0; hm < hDS; hm++){
            md.hivpop[t][g][ha][hm] *= 1+popadjrate_ha;
            if(t >= ss.t_art_start)
              for(int hu = 0; hu < hTS; hu++)
                md.artpop[t][g][ha][hm][hu] *= 1+popadjrate_ha;
          } // loop over hm
        } // loop over ha
      } // loop over g
    } // if(flag_popadjust)

    // prevalence among pregnant women

    double hivbirths = 0;
    for(int ha = hIDX_FERT; ha < hIDX_FERT+hAG_FERT; ha++){
      double hivn_ha = 0, frr_hivpop_ha = 0;
      for(int a =  ss.h_ag_start[ha]; a < ss.h_ag_start[ha]+ss.h_ag_span[ha]; a++)
        hivn_ha += (md.pop[t-1][HIVN][FEMALE][a] + md.pop[t][HIVN][FEMALE][a])/2;
      for(int hm = 0; hm < hDS; hm++){
        frr_hivpop_ha += nhp.frr_cd4[t][ha-hIDX_FERT][hm] * (md.hivpop[t-1][FEMALE][ha][hm]+md.hivpop[t][FEMALE][ha][hm])/2;
        if(t == ss.t_art_start)
          for(int hu = 0; hu < hTS; hu++)
            frr_hivpop_ha += nhp.frr_art[t][ha-hIDX_FERT][hm][hu] * md.artpop[t][FEMALE][ha][hm][hu]/2;
        else if(t > ss.t_art_start)
          for(int hu = 0; hu < hTS; hu++)
            frr_hivpop_ha += nhp.frr_art[t][ha-hIDX_FERT][hm][hu] * (md.artpop[t-1][FEMALE][ha][hm][hu]+md.artpop[t][FEMALE][ha][hm][hu])/2;
      }
      hivbirths += md.births_by_ha[ha-hIDX_FERT] * frr_hivpop_ha / (hivn_ha + frr_hivpop_ha);
    }

    md.pregprev[t] = hivbirths/md.births;
    if(t + AGE_START < ss.proj_years)
      md.pregprevlag[t + AGE_START-1] = md.pregprev[t];

    md.incid15to49[t] /= md.hivn15to49;

    // prevalence 15 to 49
    md.hivn15to49 = 0;
    for(int g = 0; g < NG; g++)
      for(int a = pIDX_15TO49; a < pIDX_15TO49 + pAG_15TO49; a++){
        md.hivn15to49 += md.pop[t][HIVN][g][a];
        md.prev15to49[t] += md.pop[t][HIVP][g][a];
      }
    md.prev15to49[t] /= (md.hivn15to49 + md.prev15to49[t]);

  }
}



void simulate_eppasm3(Model &md,
		      const SsDim &ss,
		      const DemogParam &demp,
		      const PaediatricHivParam &paedhp,
		      const NaturalHistoryParam &nhp,
		      const ArtData &artp,
		      const IncidenceParam &incp) {

  // md.initialise();

  for(int g = 0; g < ss.ng; g++)
    for(int a = 0; a < ss.p_ag; a++){
      md.pop[0][ss.i_hivn][g][a] = demp.basepop[g][a];
      for(int a = ss.i_p_ag_15to49; a < ss.i_p_ag_15to49 + ss.r_p_ag_15to49; a++)
        md.hivn15to49 += demp.basepop[g][a];
    }


  ////////////////////////////////////
  ////  do population projection  ////
  ////////////////////////////////////

  for(int t = 1; t < ss.sim_years; t++){

    // md.ageing(t);

    // age the population one year
    for(int m = 0; m < ss.p_ds; m++)
      for(int g = 0; g < ss.ng; g++){
        for(int a = 1; a < ss.p_ag; a++)
          md.pop[t][m][g][a] = md.pop[t-1][m][g][a-1];
        md.pop[t][m][g][ss.p_ag-1] += md.pop[t-1][m][g][ss.p_ag-1]; // open age group
      }

    // calculate proportion to age within each HIV age group
    double hiv_ag_prob[NG][hAG];
    for(int g = 0; g < ss.ng; g++){
      int a = 0;
      for(int ha = 0; ha < (ss.h_ag-1); ha++){
        hiv_ag_prob[g][ha] = 0;
        for(int i = 0; i < ss.h_ag_span[ha]; i++){
          hiv_ag_prob[g][ha] += md.pop[t-1][ss.i_hivp][g][a];
          a++;
        }
        hiv_ag_prob[g][ha] = (hiv_ag_prob[g][ha] > 0) ? md.pop[t-1][ss.i_hivp][g][a-1] / hiv_ag_prob[g][ha] : 0;
      }
      hiv_ag_prob[g][ss.h_ag-1] = 0.0; // no one ages out of the open-ended age group
    }

    for(int g = 0; g < ss.ng; g++){

      // youngest age group: only ageing out; ageing in entered below
      for(int hm = 0; hm < ss.h_ds; hm++){
        md.hivpop[t][g][0][hm] = (1-hiv_ag_prob[g][0]) * md.hivpop[t-1][g][0][hm];
        if(t > ss.t_art_start){
          for(int hu = 0; hu < ss.h_ts; hu++)
            md.artpop[t][g][0][hm][hu] = (1-hiv_ag_prob[g][0]) * md.artpop[t-1][g][0][hm][hu];
        }
      }

      for(int ha = 1; ha < ss.h_ag; ha++)
        for(int hm = 0; hm < ss.h_ds; hm++){
          md.hivpop[t][g][ha][hm] = (1-hiv_ag_prob[g][ha]) * md.hivpop[t-1][g][ha][hm] + hiv_ag_prob[g][ha-1] * md.hivpop[t-1][g][ha-1][hm];
          if(t > ss.t_art_start)
            for(int hu = 0; hu < ss.h_ts; hu++)
              md.artpop[t][g][ha][hm][hu] = (1-hiv_ag_prob[g][ha]) * md.artpop[t-1][g][ha][hm][hu] + hiv_ag_prob[g][ha-1] * md.artpop[t-1][g][ha-1][hm][hu];
        }
    }


    // add lagged births to youngest age group
    // md.pop_entrants(t);

    for(int g = 0; g < ss.ng; g++){

      double paedsurv_g;
      double entrant_prev;

      if(paedhp.use_entrantprev)
        entrant_prev = paedhp.entrantprev[t][g];
      else
        entrant_prev = md.pregprevlag[t-1] * paedhp.verttrans_lag[t-1] * paedhp.paedsurv_lag[t-1];

      if(demp.flag_popadjust){
        md.pop[t][ss.i_hivn][g][0] = demp.entrantpop[t-1][g] * (1.0-entrant_prev);
        paedsurv_g = demp.entrantpop[t-1][g] * entrant_prev;
      } else {
        md.pop[t][ss.i_hivn][g][0] = md.birthslag[t-1][g] * demp.cumsurv[t-1][g] * (1.0-entrant_prev / paedhp.paedsurv_lag[t-1]) + demp.cumnetmigr[t-1][g] * (1.0-md.pregprevlag[t-1] * paedhp.netmig_hivprob);
        paedsurv_g = md.birthslag[t-1][g] * demp.cumsurv[t-1][g] * entrant_prev + demp.cumnetmigr[t-1][g] * entrant_prev;
      }

      md.pop[t][ss.i_hivp][g][0] = paedsurv_g;

      md.entrantprev[t] = (md.pop[t][ss.i_hivp][ss.i_male][0] + md.pop[t][ss.i_hivp][ss.i_female][0]) / (md.pop[t][ss.i_hivn][ss.i_male][0] + md.pop[t][ss.i_hivn][ss.i_female][0] + md.pop[t][ss.i_hivp][ss.i_male][0] + md.pop[t][ss.i_hivp][ss.i_female][0]);

      for(int hm = 0; hm < ss.h_ds; hm++){
        md.hivpop[t][g][0][hm] += paedsurv_g * paedhp.paedsurv_cd4dist[t][g][hm] * (1.0 - paedhp.entrantartcov[t][g]);
        if(t > ss.t_art_start){
          for(int hu = 0; hu < ss.h_ts; hu++){
            md.artpop[t][g][0][hm][hu] += paedsurv_g * paedhp.paedsurv_artcd4dist[t][g][hm][hu] * paedhp.entrantartcov[t][g];
          }
        }
      }
    }

    // non-HIV mortality and netmigration
    // md.death_and_migration(t);

    for(int g = 0; g < ss.ng; g++){
      int a = 0;
      for(int ha = 0; ha < ss.h_ag; ha++){
        double deathsmig_ha = 0, hivpop_ha = 0;
        for(int i = 0; i < ss.h_ag_span[ha]; i++){

          hivpop_ha += md.pop[t][ss.i_hivp][g][a];

          // non-HIV mortality
          double qx = 1.0 - demp.Sx[t][g][a];
          double ndeaths_a = md.pop[t][ss.i_hivn][g][a] * qx;
          md.pop[t][ss.i_hivn][g][a] -= ndeaths_a; // survival HIV- population
          double hdeaths_a = md.pop[t][ss.i_hivp][g][a] * qx;
          deathsmig_ha -= hdeaths_a;
          md.pop[t][ss.i_hivp][g][a] -= hdeaths_a;   // survival HIV+ population
          md.natdeaths[t][g][a] = ndeaths_a + hdeaths_a;

          // net migration
          double migrate_a = demp.netmigr[t][g][a] * (1+demp.Sx[t][g][a])/2.0 / (md.pop[t][ss.i_hivn][g][a] + md.pop[t][ss.i_hivp][g][a]);
          md.pop[t][ss.i_hivn][g][a] *= 1+migrate_a;
          double hmig_a = migrate_a * md.pop[t][ss.i_hivp][g][a];
          deathsmig_ha += hmig_a;
          md.pop[t][ss.i_hivp][g][a] += hmig_a;

          a++;
        }

        // migration and deaths for hivpop
        double deathmigrate_ha = hivpop_ha > 0 ? deathsmig_ha / hivpop_ha : 0.0;
        for(int hm = 0; hm < ss.h_ds; hm++){
          md.hivpop[t][g][ha][hm] *= 1+deathmigrate_ha;
          if(t > ss.t_art_start)
            for(int hu = 0; hu < ss.h_ts; hu++)
              md.artpop[t][g][ha][hm][hu] *= 1+deathmigrate_ha;
        } // loop over hm
      } // loop over ha
    } // loop over g


    // md.fertility(t);

    // fertility
    md.births = 0.0;
    std::fill(md.births_by_ha.begin(), md.births_by_ha.end(), 0);
    for(int m = 0; m < ss.p_ds; m++){
      int a = ss.i_p_ag_fert;
      for(int ha = ss.i_h_ag_fert; ha < ss.r_h_ag_fert; ha++){
        for(int i = 0; i < ss.h_ag_span[ha]; i++){
          md.births_by_ha[ha-ss.i_h_ag_fert] += (md.pop[t-1][m][ss.i_female][a] + md.pop[t][m][ss.i_female][a])/2 * demp.asfr[t][a];
          a++;
        }
      }
    }
    for(int ha = ss.i_h_ag_fert; ha < ss.r_h_ag_fert; ha++)
      md.births += md.births_by_ha[ha-ss.i_h_ag_fert];

    if(t + ss.pop_age_start < ss.proj_years)
      for(int g = 0; g < ss.ng; g++)
        md.birthslag[t + ss.pop_age_start-1][g] = demp.srb[t][g] * md.births;

    
    ////////////////////////////////
    ////  HIV model simulation  ////
    ////////////////////////////////

    int cd4elig_idx = artp.artcd4elig_idx[t] - 1; // -1 for 0-based indexing vs. 1-based in R
    int anyelig_idx = (artp.special_pop_percelig[t] > 0 | artp.preg_women_artelig[t] > 0) ? 0 : (artp.who34percelig > 0) ? nhp.i_h_cd4_350 : cd4elig_idx;
    md.everARTelig_idx = anyelig_idx < md.everARTelig_idx ? anyelig_idx : md.everARTelig_idx;

    for(int hts = 0; hts < ss.hiv_steps_per_year; hts++){

      int ts = (t-1)*ss.hiv_steps_per_year + hts;

      double hivdeaths_ha[NG][hAG];
      memset(hivdeaths_ha, 0, sizeof(double)*ss.ng*ss.h_ag);

      // untreated population

      // disease progression and mortality
      double grad[NG][hAG][hDS];
      for(int g = 0; g < ss.ng; g++)
        for(int ha = 0; ha < ss.h_ag; ha++){
          for(int hm = 0; hm < ss.h_ds; hm++){

            grad[g][ha][hm] = 0.0;

            double cd4mx_scale = 1.0;
            if(artp.scale_cd4_mort & (t >= ss.t_art_start) & (hm >= md.everARTelig_idx)){
              double artpop_hahm = 0.0;
              for(int hu = 0; hu < ss.h_ts; hu++)
                artpop_hahm += md.artpop[t][g][ha][hm][hu];
              cd4mx_scale = md.hivpop[t][g][ha][hm] / (md.hivpop[t][g][ha][hm] + artpop_hahm);
            }

            double deaths = cd4mx_scale * nhp.cd4_mort[g][ha][hm] * md.hivpop[t][g][ha][hm];
            hivdeaths_ha[g][ha] += ss.dt * deaths;
            grad[g][ha][hm] = -deaths;
          }
          for(int hm = 1; hm < ss.h_ds; hm++){
            grad[g][ha][hm-1] -= nhp.cd4_prog[g][ha][hm-1] * md.hivpop[t][g][ha][hm-1];
            grad[g][ha][hm] += nhp.cd4_prog[g][ha][hm-1] * md.hivpop[t][g][ha][hm-1];
          }
        }

      if(incp.eppmod != EPP_DIRECTINCID){
        // incidence

        // calculate r(t)
        if(incp.eppmod == EPP_RSPLINE)
          md.rvec_ts[ts] = incp.rspline_rvec[ts];
        else
          md.rvec_ts[ts] = calc_rtrend_rt(md.pop, incp.rtrend_tstab, incp.rtrend_beta, incp.rtrend_r0,
                                          ss.proj_start + 0.5 + ts * 1.0 / ss.hiv_steps_per_year,
                                          ss.proj_start + 0.5 + incp.ts_epidemic_start * 1.0 / ss.hiv_steps_per_year,
                                          ss.dt, t, hts,
                                          md.rvec_ts[ts-1], &md.prevlast, &md.prevcurr);

        // calculate new infections by sex and age
        double infections_ts[NG][pAG];
        calc_infections_eppspectrum2(md.pop, md.hivpop, md.artpop,
                                     md.rvec_ts[ts], incp.rel_infect_art, (ts == incp.ts_epidemic_start) ? incp.iota : 0.0,
                                     incp.incrr_sex, incp.incrr_age,
                                     ss.t_art_start, ss.dt, t, hts, ss.h_ag_start.data(), ss.h_ag_span,
                                     &md.prevcurr, &md.incid15to49_ts[ts], infections_ts);

        md.prev15to49_ts[ts] = md.prevcurr;

        // add new infections to HIV population
        for(int g = 0; g < ss.ng; g++){
          int a = 0;
          for(int ha = 0; ha < ss.h_ag; ha++){
            double infections_a, infections_ha = 0.0;
            for(int i = 0; i < ss.h_ag_span[ha]; i++){
              infections_ha += infections_a = infections_ts[g][a];
              md.infections[t][g][a] += ss.dt*infections_a;
              md.pop[t][ss.i_hivn][g][a] -= ss.dt*infections_a;
              md.pop[t][ss.i_hivp][g][a] += ss.dt*infections_a;
              a++;
            }
	    if(ha < ss.i_h_ag_15to49+ss.r_h_ag_15to49)
              md.incid15to49[t] += ss.dt*infections_ha;

            // add infections to grad hivpop
            for(int hm = 0; hm < ss.h_ds; hm++)
              grad[g][ha][hm] += infections_ha * nhp.cd4_initdist[g][ha][hm];
          }
        }
      }

      // ART progression, mortality, and initiation
      if(t >= ss.t_art_start){

        double gradART[NG][hAG][hDS][hTS];

        // progression and mortality
        for(int g = 0; g < ss.ng; g++)
          for(int ha = 0; ha < ss.h_ag; ha++)
            for(int hm = md.everARTelig_idx; hm < ss.h_ds; hm++){

              for(int hu = 0; hu < ss.h_ts; hu++){
                double deaths = nhp.art_mort[g][ha][hm][hu] * nhp.artmx_timerr[t][hu] * md.artpop[t][g][ha][hm][hu];
                hivdeaths_ha[g][ha] += ss.dt*deaths;
                gradART[g][ha][hm][hu] = -deaths;
              }

              gradART[g][ha][hm][nhp.art0mos] += -nhp.art_stage_prog_rate * md.artpop[t][g][ha][hm][nhp.art0mos];
              gradART[g][ha][hm][nhp.art6mos] += nhp.art_stage_prog_rate * md.artpop[t][g][ha][hm][nhp.art0mos] - nhp.art_stage_prog_rate * md.artpop[t][g][ha][hm][nhp.art6mos];
              gradART[g][ha][hm][nhp.art1yr] += nhp.art_stage_prog_rate * md.artpop[t][g][ha][hm][nhp.art6mos];

              // ART dropout
              if(artp.art_dropout[t] > 0)
                for(int hu = 0; hu < ss.h_ts; hu++){
                  grad[g][ha][hm] += artp.art_dropout[t] * md.artpop[t][g][ha][hm][hu];
                  gradART[g][ha][hm][hu] -= artp.art_dropout[t] * md.artpop[t][g][ha][hm][hu];
                }

            }


        // ART initiation
        for(int g = 0; g < ss.ng; g++){

          double artelig_hahm[hAG_15PLUS][hDS], Xart_15plus = 0.0, Xartelig_15plus = 0.0, expect_mort_artelig15plus = 0.0;
          for(int ha = ss.i_h_ag_15plus; ha < ss.h_ag; ha++){
            for(int hm = md.everARTelig_idx; hm < ss.h_ds; hm++){
              if(hm >= anyelig_idx){
                double prop_elig = (hm >= cd4elig_idx) ? 1.0 : (hm >= nhp.i_h_cd4_350) ? 1.0 - (1.0-artp.special_pop_percelig[t])*(1.0-artp.who34percelig) : artp.special_pop_percelig[t];
                Xartelig_15plus += artelig_hahm[ha-ss.i_h_ag_15plus][hm] = prop_elig * md.hivpop[t][g][ha][hm] ;
                expect_mort_artelig15plus += nhp.cd4_mort[g][ha][hm] * artelig_hahm[ha-ss.i_h_ag_15plus][hm];
              }
              for(int hu = 0; hu < ss.h_ts; hu++)
                Xart_15plus += md.artpop[t][g][ha][hm][hu] + ss.dt * gradART[g][ha][hm][hu];
            }

            // if artp.preg_women_artelig, add pregnant women to artelig_hahm population
            if(g == ss.i_female & artp.preg_women_artelig[t] > 0 & ha < ss.r_h_ag_fert){
              double frr_pop_ha = 0;
              for(int a =  ss.h_ag_start[ha]; a < ss.h_ag_start[ha]+ss.h_ag_span[ha]; a++)
                frr_pop_ha += md.pop[t][ss.i_hivn][g][a]; // add HIV- population
              for(int hm = 0; hm < ss.h_ds; hm++){
                frr_pop_ha += nhp.frr_cd4[t][ha-ss.i_h_ag_fert][hm] * md.hivpop[t][g][ha][hm];
                for(int hu = 0; hu < ss.h_ts; hu++)
                  frr_pop_ha += nhp.frr_art[t][ha-ss.i_h_ag_fert][hm][hu] * md.artpop[t][g][ha][hm][hu];
              }
              for(int hm = anyelig_idx; hm < cd4elig_idx; hm++){
                double pw_elig_hahm = md.births_by_ha[ha-ss.i_h_ag_fert] * nhp.frr_cd4[t][ha-ss.i_h_ag_fert][hm] * md.hivpop[t][g][ha][hm] / frr_pop_ha;
                artelig_hahm[ha-ss.i_h_ag_15plus][hm] += pw_elig_hahm;
                Xartelig_15plus += pw_elig_hahm;
                expect_mort_artelig15plus += nhp.cd4_mort[g][ha][hm] * pw_elig_hahm;
              }
            }
          } // loop over ha

            // calculate number on ART at end of ts, based on number or percent
          double artnum_hts = 0.0;
          if(ss.dt*(hts+1) < 0.5){
            if(!artp.art15plus_isperc[t-2][g] & !artp.art15plus_isperc[t-1][g]){ // both numbers
              artnum_hts = (0.5-ss.dt*(hts+1))*artp.artnum15plus[t-2][g] + (ss.dt*(hts+1)+0.5)*artp.artnum15plus[t-1][g];
            } else if(artp.art15plus_isperc[t-2][g] & artp.art15plus_isperc[t-1][g]){ // both percentages
              double artcov_hts = (0.5-ss.dt*(hts+1))*artp.artnum15plus[t-2][g] + (ss.dt*(hts+1)+0.5)*artp.artnum15plus[t-1][g];
              artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
            } else if((!artp.art15plus_isperc[t-2][g]) & artp.art15plus_isperc[t-1][g]){ // transition from number to percentage
              double curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus);
              double artcov_hts = curr_coverage + (artp.artnum15plus[t-1][g] - curr_coverage) * ss.dt / (0.5-ss.dt*hts);
              artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
            }
          } else {
            if(!artp.art15plus_isperc[t-1][g] & !artp.art15plus_isperc[t][g]){ // both numbers
              artnum_hts = (1.5-ss.dt*(hts+1))*artp.artnum15plus[t-1][g] + (ss.dt*(hts+1)-0.5)*artp.artnum15plus[t][g];
            } else if(artp.art15plus_isperc[t-1][g] & artp.art15plus_isperc[t][g]){ // both percentages
              double artcov_hts = (1.5-ss.dt*(hts+1))*artp.artnum15plus[t-1][g] + (ss.dt*(hts+1)-0.5)*artp.artnum15plus[t][g];
              artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
            } else if((!artp.art15plus_isperc[t-1][g]) & artp.art15plus_isperc[t][g]){ // transition from number to percentage
              double curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus);
              double artcov_hts = curr_coverage + (artp.artnum15plus[t][g] - curr_coverage) * ss.dt / (1.5-ss.dt*hts);
              artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
            }
          }

          double artinit_hts = artnum_hts > Xart_15plus ? artnum_hts - Xart_15plus : 0;

          // median CD4 at initiation inputs
          if(artp.med_cd4init_input[t]){

            const int CD4_LOW_LIM[hDS] = {500, 350, 250, 200, 100, 50, 0};
            const int CD4_UPP_LIM[hDS] = {1000, 500, 350, 250, 200, 100, 50};

            int medcd4_idx = artp.med_cd4init_cat[t] - 1; // -1 for 0-based indexing vs. 1-based in R
            double medcat_propbelow = (artp.median_cd4init[t] - CD4_LOW_LIM[medcd4_idx]) / (CD4_UPP_LIM[medcd4_idx] - CD4_LOW_LIM[medcd4_idx]);

            double elig_below = 0.0, elig_above = 0.0;
            for(int ha = ss.i_h_ag_15plus; ha < ss.h_ag; ha++){
              for(int hm = anyelig_idx; hm < medcd4_idx; hm++)
                elig_above += artelig_hahm[ha-ss.i_h_ag_15plus][hm];
              elig_above += (1.0 - medcat_propbelow) * artelig_hahm[ha-ss.i_h_ag_15plus][medcd4_idx];
              elig_below += medcat_propbelow * artelig_hahm[ha-ss.i_h_ag_15plus][medcd4_idx];
              for(int hm = medcd4_idx+1; hm < ss.h_ds; hm++)
                elig_below += artelig_hahm[ha-ss.i_h_ag_15plus][hm];
            }

            double initprob_below = (elig_below > artinit_hts * 0.5) ? artinit_hts * 0.5 / elig_below : 1.0;
            double initprob_above = (elig_above > artinit_hts * 0.5) ? artinit_hts * 0.5 / elig_above : 1.0;
            double initprob_medcat = initprob_below * medcat_propbelow + initprob_above * (1.0-medcat_propbelow);

            for(int ha = ss.i_h_ag_15plus; ha < ss.h_ag; ha++)
              for(int hm = anyelig_idx; hm < ss.h_ds; hm++){
                double artinit_hahm;
                if(hm < medcd4_idx)
                  artinit_hahm = artelig_hahm[ha-ss.i_h_ag_15plus][hm] * initprob_above;
                else if(hm == medcd4_idx)
                  artinit_hahm = artelig_hahm[ha-ss.i_h_ag_15plus][hm] * initprob_medcat;
                if(hm > medcd4_idx)
                  artinit_hahm = artelig_hahm[ha-ss.i_h_ag_15plus][hm] * initprob_below;
                if(artinit_hahm > md.hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm])
                  artinit_hahm = md.hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm];
                grad[g][ha][hm] -= artinit_hahm / ss.dt;
                gradART[g][ha][hm][nhp.art0mos] += artinit_hahm / ss.dt;
              }

          } else if(artp.art_alloc_method == 4) {  // lowest CD4 first

            for(int hm = ss.h_ds-1; hm >= anyelig_idx; hm--){
              double artelig_hm = 0;
              for(int ha = ss.i_h_ag_15plus; ha < ss.h_ag; ha++)
                artelig_hm += artelig_hahm[ha-ss.i_h_ag_15plus][hm];
              double init_prop = (artelig_hm == 0 | artinit_hts > artelig_hm) ? 1.0 : artinit_hts / artelig_hm;

              for(int ha = ss.i_h_ag_15plus; ha < ss.h_ag; ha++){
                double artinit_hahm = init_prop * artelig_hahm[ha-ss.i_h_ag_15plus][hm];

                if(artinit_hahm > md.hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm])
                  artinit_hahm = md.hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm];

                grad[g][ha][hm] -= artinit_hahm / ss.dt;
                gradART[g][ha][hm][nhp.art0mos] += artinit_hahm / ss.dt;
              }
              if(init_prop < 1.0)
                break;
              artinit_hts -= init_prop * artelig_hm;
            }

          } else { // Use mixture of eligibility and expected mortality for initiation distribution

            for(int ha = ss.i_h_ag_15plus; ha < ss.h_ag; ha++)
              for(int hm = anyelig_idx; hm < ss.h_ds; hm++){
                double artinit_hahm = artinit_hts * artelig_hahm[ha-ss.i_h_ag_15plus][hm] * ((1.0 - artp.art_alloc_mxweight)/Xartelig_15plus + artp.art_alloc_mxweight * nhp.cd4_mort[g][ha][hm] / expect_mort_artelig15plus);
                if(artinit_hahm > artelig_hahm[ha-ss.i_h_ag_15plus][hm])
                  artinit_hahm = artelig_hahm[ha-ss.i_h_ag_15plus][hm];
                if(artinit_hahm > md.hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm])
                  artinit_hahm = md.hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm];
                grad[g][ha][hm] -= artinit_hahm / ss.dt;
                gradART[g][ha][hm][nhp.art0mos] += artinit_hahm / ss.dt;
              }
          }
        }

        for(int g = 0; g < ss.ng; g++)
          for(int ha = 0; ha < ss.h_ag; ha++)
            for(int hm = md.everARTelig_idx; hm < ss.h_ds; hm++)
              for(int hu = 0; hu < ss.h_ts; hu++)
                md.artpop[t][g][ha][hm][hu] += ss.dt*gradART[g][ha][hm][hu];

      } // if(t >= ss.t_art_start)

      for(int g = 0; g < ss.ng; g++)
        for(int ha = 0; ha < ss.h_ag; ha++)
          for(int hm = 0; hm < ss.h_ds; hm++)
            md.hivpop[t][g][ha][hm] += ss.dt*grad[g][ha][hm];


      // remove hivdeaths from pop
      for(int g = 0; g < ss.ng; g++){

        // sum HIV+ population size in each hivpop age group
        double hivpop_ha[hAG];
        int a = 0;
        for(int ha = 0; ha < ss.h_ag; ha++){
          hivpop_ha[ha] = 0.0;
          for(int i = 0; i < ss.h_ag_span[ha]; i++){
            hivpop_ha[ha] += md.pop[t][ss.i_hivp][g][a];
            a++;
          }
        }

        // remove hivdeaths proportionally to age-distribution within each age group
        a = 0;
        for(int ha = 0; ha < ss.h_ag; ha++){
          if(hivpop_ha[ha] > 0){
            double hivqx_ha = hivdeaths_ha[g][ha] / hivpop_ha[ha];
            for(int i = 0; i < ss.h_ag_span[ha]; i++){
              md.hivdeaths[t][g][a] += md.pop[t][ss.i_hivp][g][a] * hivqx_ha;
              md.pop[t][ss.i_hivp][g][a] *= (1.0-hivqx_ha);
              a++;
            }
          } else {
            a += ss.h_ag_span[ha];
          }  // end if(pop_ha[ha] > 0)
        }
      }



    } // loop HIVSTEPS_PER_YEAR



    if(incp.eppmod == EPP_DIRECTINCID){
      // Calculating new infections once per year (like Spectrum)

      double Xhivp = 0.0, Xhivn[NG], Xhivn_incagerr[NG];

      for(int g = 0; g < ss.ng; g++){
        Xhivn[g] = 0.0;
        Xhivn_incagerr[g] = 0.0;
        for(int a = incp.i_p_ag_incidpop; a < incp.i_p_ag_incidpop + incp.r_p_ag_incidpop; a++){
          Xhivp += md.pop[t-1][ss.i_hivp][g][a];
          Xhivn[g] += md.pop[t-1][ss.i_hivn][g][a];
          Xhivn_incagerr[g] += incp.incrr_age[t][g][a] * md.pop[t-1][ss.i_hivn][g][a];
        }
      }
      // double prev_i = Xhivp / (Xhivn[ss.i_male] + Xhivn[ss.i_female] + Xhivp);
      // double incrate15to49_i = (prev15to49[t] - prev_i)/(1.0 - prev_i);
      double incrate_i = incp.incidinput[t];
      double incrate_g[NG];
      incrate_g[ss.i_male] = incrate_i * (Xhivn[ss.i_male]+Xhivn[ss.i_female]) / (Xhivn[ss.i_male] + incp.incrr_sex[t]*Xhivn[ss.i_female]);
      incrate_g[ss.i_female] = incrate_i * incp.incrr_sex[t]*(Xhivn[ss.i_male]+Xhivn[ss.i_female]) / (Xhivn[ss.i_male] + incp.incrr_sex[t]*Xhivn[ss.i_female]);

      for(int g = 0; g < ss.ng; g++){
        int a = 0;
        for(int ha = 0; ha < ss.h_ag; ha++){
          double infections_a, infections_ha = 0.0;
          for(int i = 0; i < ss.h_ag_span[ha]; i++){
            infections_ha += infections_a = md.pop[t-1][ss.i_hivn][g][a] * incrate_g[g] * incp.incrr_age[t][g][a] * Xhivn[g] / Xhivn_incagerr[g];
            md.infections[t][g][a] += infections_a;
            md.pop[t][ss.i_hivn][g][a] -= infections_a;
            md.pop[t][ss.i_hivp][g][a] += infections_a;
            a++;
          }
          if(ha < ss.i_h_ag_15to49+ss.r_h_ag_15to49)
            md.incid15to49[t] += infections_ha;

          // add infections to hivpop
          for(int hm = 0; hm < ss.h_ds; hm++)
            md.hivpop[t][g][ha][hm] += infections_ha * nhp.cd4_initdist[g][ha][hm];
        }
      }
    }


    // md.adjust_population(t);

    // adjust population to match target population
    if(demp.flag_popadjust){
      for(int g = 0; g < ss.ng; g++){
        int a = 0;
        for(int ha = 0; ha < ss.h_ag; ha++){
          double popadj_ha = 0, hivpop_ha = 0;
          for(int i = 0; i < ss.h_ag_span[ha]; i++){

            hivpop_ha += md.pop[t][ss.i_hivp][g][a];

            double popadjrate_a = md.popadjust[t][g][a] = demp.targetpop[t][g][a] / (md.pop[t][ss.i_hivn][g][a] + md.pop[t][ss.i_hivp][g][a]);
            md.pop[t][ss.i_hivn][g][a] *= popadjrate_a;
            double hpopadj_a = (popadjrate_a-1.0) * md.pop[t][ss.i_hivp][g][a];
            popadj_ha += hpopadj_a;
            md.pop[t][ss.i_hivp][g][a] += hpopadj_a;
            a++;
          }

          // population adjustment for hivpop
          double popadjrate_ha = hivpop_ha > 0 ? popadj_ha / hivpop_ha : 0.0;
          for(int hm = 0; hm < ss.h_ds; hm++){
            md.hivpop[t][g][ha][hm] *= 1+popadjrate_ha;
            if(t >= ss.t_art_start)
              for(int hu = 0; hu < ss.h_ts; hu++)
                md.artpop[t][g][ha][hm][hu] *= 1+popadjrate_ha;
          } // loop over hm
        } // loop over ha
      } // loop over g
    } // if(flag_popadjust)

    // prevalence among pregnant women

    double hivbirths = 0;
    for(int ha = ss.i_h_ag_fert; ha < ss.i_h_ag_fert+ss.r_h_ag_fert; ha++){
      double hivn_ha = 0, frr_hivpop_ha = 0;
      for(int a =  ss.h_ag_start[ha]; a < ss.h_ag_start[ha]+ss.h_ag_span[ha]; a++)
        hivn_ha += (md.pop[t-1][ss.i_hivn][ss.i_female][a] + md.pop[t][ss.i_hivn][ss.i_female][a])/2;
      for(int hm = 0; hm < ss.h_ds; hm++){
        frr_hivpop_ha += nhp.frr_cd4[t][ha-ss.i_h_ag_fert][hm] * (md.hivpop[t-1][ss.i_female][ha][hm]+md.hivpop[t][ss.i_female][ha][hm])/2;
        if(t == ss.t_art_start)
          for(int hu = 0; hu < ss.h_ts; hu++)
            frr_hivpop_ha += nhp.frr_art[t][ha-ss.i_h_ag_fert][hm][hu] * md.artpop[t][ss.i_female][ha][hm][hu]/2;
        else if(t > ss.t_art_start)
          for(int hu = 0; hu < ss.h_ts; hu++)
            frr_hivpop_ha += nhp.frr_art[t][ha-ss.i_h_ag_fert][hm][hu] * (md.artpop[t-1][ss.i_female][ha][hm][hu]+md.artpop[t][ss.i_female][ha][hm][hu])/2;
      }
      hivbirths += md.births_by_ha[ha-ss.i_h_ag_fert] * frr_hivpop_ha / (hivn_ha + frr_hivpop_ha);
    }

    md.pregprev[t] = hivbirths/md.births;
    if(t + ss.pop_age_start < ss.proj_years)
      md.pregprevlag[t + ss.pop_age_start-1] = md.pregprev[t];

    md.incid15to49[t] /= md.hivn15to49;

    // prevalence 15 to 49
    md.hivn15to49 = 0;
    for(int g = 0; g < ss.ng; g++)
      for(int a = ss.i_p_ag_15to49; a < ss.i_p_ag_15to49 + ss.r_p_ag_15to49; a++){
        md.hivn15to49 += md.pop[t][ss.i_hivn][g][a];
        md.prev15to49[t] += md.pop[t][ss.i_hivp][g][a];
      }
    md.prev15to49[t] /= (md.hivn15to49 + md.prev15to49[t]);

  }
}



void simulate_eppasm4(Model &md,
		      const SsDim &ss,
		      const DemogParam &demp,
		      const PaediatricHivParam &paedhp,
		      const NaturalHistoryParam &nhp,
		      const ArtData &artp,
		      const IncidenceParam &incp) {

  // md.initialise();

  for(int g = 0; g < ss.ng; g++)
    for(int a = 0; a < ss.p_ag; a++){
      md.pop[0][ss.i_hivn][g][a] = demp.basepop[g][a];
      for(int a = ss.i_p_ag_15to49; a < ss.i_p_ag_15to49 + ss.r_p_ag_15to49; a++)
        md.hivn15to49 += demp.basepop[g][a];
    }


  ////////////////////////////////////
  ////  do population projection  ////
  ////////////////////////////////////

  for(int t = 1; t < ss.sim_years; t++){

    // md.ageing(t);

    // age the population one year
    for(int m = 0; m < ss.p_ds; m++)
      for(int g = 0; g < ss.ng; g++){
        for(int a = 1; a < ss.p_ag; a++)
          md.pop[t][m][g][a] = md.pop[t-1][m][g][a-1];
        md.pop[t][m][g][ss.p_ag-1] += md.pop[t-1][m][g][ss.p_ag-1]; // open age group
      }

    // calculate proportion to age within each HIV age group
    multi_array<double, 2> hiv_ag_prob(extents[ss.ng][ss.h_ag]);
    for(int g = 0; g < ss.ng; g++){
      int a = 0;
      for(int ha = 0; ha < (ss.h_ag-1); ha++){
        hiv_ag_prob[g][ha] = 0;
        for(int i = 0; i < ss.h_ag_span[ha]; i++){
          hiv_ag_prob[g][ha] += md.pop[t-1][ss.i_hivp][g][a];
          a++;
        }
        hiv_ag_prob[g][ha] = (hiv_ag_prob[g][ha] > 0) ? md.pop[t-1][ss.i_hivp][g][a-1] / hiv_ag_prob[g][ha] : 0;
      }
      hiv_ag_prob[g][ss.h_ag-1] = 0.0; // no one ages out of the open-ended age group
    }

    for(int g = 0; g < ss.ng; g++){

      // youngest age group: only ageing out; ageing in entered below
      for(int hm = 0; hm < ss.h_ds; hm++){
        md.hivpop[t][g][0][hm] = (1-hiv_ag_prob[g][0]) * md.hivpop[t-1][g][0][hm];
        if(t > ss.t_art_start){
          for(int hu = 0; hu < ss.h_ts; hu++)
            md.artpop[t][g][0][hm][hu] = (1-hiv_ag_prob[g][0]) * md.artpop[t-1][g][0][hm][hu];
        }
      }

      for(int ha = 1; ha < ss.h_ag; ha++)
        for(int hm = 0; hm < ss.h_ds; hm++){
          md.hivpop[t][g][ha][hm] = (1-hiv_ag_prob[g][ha]) * md.hivpop[t-1][g][ha][hm] + hiv_ag_prob[g][ha-1] * md.hivpop[t-1][g][ha-1][hm];
          if(t > ss.t_art_start)
            for(int hu = 0; hu < ss.h_ts; hu++)
              md.artpop[t][g][ha][hm][hu] = (1-hiv_ag_prob[g][ha]) * md.artpop[t-1][g][ha][hm][hu] + hiv_ag_prob[g][ha-1] * md.artpop[t-1][g][ha-1][hm][hu];
        }
    }


    // add lagged births to youngest age group
    // md.pop_entrants(t);

    for(int g = 0; g < ss.ng; g++){

      double paedsurv_g;
      double entrant_prev;

      if(paedhp.use_entrantprev)
        entrant_prev = paedhp.entrantprev[t][g];
      else
        entrant_prev = md.pregprevlag[t-1] * paedhp.verttrans_lag[t-1] * paedhp.paedsurv_lag[t-1];

      if(demp.flag_popadjust){
        md.pop[t][ss.i_hivn][g][0] = demp.entrantpop[t-1][g] * (1.0-entrant_prev);
        paedsurv_g = demp.entrantpop[t-1][g] * entrant_prev;
      } else {
        md.pop[t][ss.i_hivn][g][0] = md.birthslag[t-1][g] * demp.cumsurv[t-1][g] * (1.0-entrant_prev / paedhp.paedsurv_lag[t-1]) + demp.cumnetmigr[t-1][g] * (1.0-md.pregprevlag[t-1] * paedhp.netmig_hivprob);
        paedsurv_g = md.birthslag[t-1][g] * demp.cumsurv[t-1][g] * entrant_prev + demp.cumnetmigr[t-1][g] * entrant_prev;
      }

      md.pop[t][ss.i_hivp][g][0] = paedsurv_g;

      md.entrantprev[t] = (md.pop[t][ss.i_hivp][ss.i_male][0] + md.pop[t][ss.i_hivp][ss.i_female][0]) / (md.pop[t][ss.i_hivn][ss.i_male][0] + md.pop[t][ss.i_hivn][ss.i_female][0] + md.pop[t][ss.i_hivp][ss.i_male][0] + md.pop[t][ss.i_hivp][ss.i_female][0]);

      for(int hm = 0; hm < ss.h_ds; hm++){
        md.hivpop[t][g][0][hm] += paedsurv_g * paedhp.paedsurv_cd4dist[t][g][hm] * (1.0 - paedhp.entrantartcov[t][g]);
        if(t > ss.t_art_start){
          for(int hu = 0; hu < ss.h_ts; hu++){
            md.artpop[t][g][0][hm][hu] += paedsurv_g * paedhp.paedsurv_artcd4dist[t][g][hm][hu] * paedhp.entrantartcov[t][g];
          }
        }
      }
    }

    // non-HIV mortality and netmigration
    // md.death_and_migration(t);

    for(int g = 0; g < ss.ng; g++){
      int a = 0;
      for(int ha = 0; ha < ss.h_ag; ha++){
        double deathsmig_ha = 0, hivpop_ha = 0;
        for(int i = 0; i < ss.h_ag_span[ha]; i++){

          hivpop_ha += md.pop[t][ss.i_hivp][g][a];

          // non-HIV mortality
          double qx = 1.0 - demp.Sx[t][g][a];
          double ndeaths_a = md.pop[t][ss.i_hivn][g][a] * qx;
          md.pop[t][ss.i_hivn][g][a] -= ndeaths_a; // survival HIV- population
          double hdeaths_a = md.pop[t][ss.i_hivp][g][a] * qx;
          deathsmig_ha -= hdeaths_a;
          md.pop[t][ss.i_hivp][g][a] -= hdeaths_a;   // survival HIV+ population
          md.natdeaths[t][g][a] = ndeaths_a + hdeaths_a;

          // net migration
          double migrate_a = demp.netmigr[t][g][a] * (1+demp.Sx[t][g][a])/2.0 / (md.pop[t][ss.i_hivn][g][a] + md.pop[t][ss.i_hivp][g][a]);
          md.pop[t][ss.i_hivn][g][a] *= 1+migrate_a;
          double hmig_a = migrate_a * md.pop[t][ss.i_hivp][g][a];
          deathsmig_ha += hmig_a;
          md.pop[t][ss.i_hivp][g][a] += hmig_a;

          a++;
        }

        // migration and deaths for hivpop
        double deathmigrate_ha = hivpop_ha > 0 ? deathsmig_ha / hivpop_ha : 0.0;
        for(int hm = 0; hm < ss.h_ds; hm++){
          md.hivpop[t][g][ha][hm] *= 1+deathmigrate_ha;
          if(t > ss.t_art_start)
            for(int hu = 0; hu < ss.h_ts; hu++)
              md.artpop[t][g][ha][hm][hu] *= 1+deathmigrate_ha;
        } // loop over hm
      } // loop over ha
    } // loop over g


    // md.fertility(t);

    // fertility
    md.births = 0.0;
    std::fill(md.births_by_ha.begin(), md.births_by_ha.end(), 0);
    for(int m = 0; m < ss.p_ds; m++){
      int a = ss.i_p_ag_fert;
      for(int ha = ss.i_h_ag_fert; ha < ss.r_h_ag_fert; ha++){
        for(int i = 0; i < ss.h_ag_span[ha]; i++){
          md.births_by_ha[ha-ss.i_h_ag_fert] += (md.pop[t-1][m][ss.i_female][a] + md.pop[t][m][ss.i_female][a])/2 * demp.asfr[t][a];
          a++;
        }
      }
    }
    for(int ha = ss.i_h_ag_fert; ha < ss.r_h_ag_fert; ha++)
      md.births += md.births_by_ha[ha-ss.i_h_ag_fert];

    if(t + ss.pop_age_start < ss.proj_years)
      for(int g = 0; g < ss.ng; g++)
        md.birthslag[t + ss.pop_age_start-1][g] = demp.srb[t][g] * md.births;

    
    ////////////////////////////////
    ////  HIV model simulation  ////
    ////////////////////////////////

    int cd4elig_idx = artp.artcd4elig_idx[t] - 1; // -1 for 0-based indexing vs. 1-based in R
    int anyelig_idx = (artp.special_pop_percelig[t] > 0 | artp.preg_women_artelig[t] > 0) ? 0 : (artp.who34percelig > 0) ? nhp.i_h_cd4_350 : cd4elig_idx;
    md.everARTelig_idx = anyelig_idx < md.everARTelig_idx ? anyelig_idx : md.everARTelig_idx;

    for(int hts = 0; hts < ss.hiv_steps_per_year; hts++){

      int ts = (t-1)*ss.hiv_steps_per_year + hts;

      multi_array<double, 2> hivdeaths_ha(extents[ss.ng][ss.h_ag]);
      memset(hivdeaths_ha.data(), 0, sizeof(double)*ss.ng*ss.h_ag);

      // untreated population

      // disease progression and mortality
      multi_array<double, 3> grad(extents[ss.ng][ss.h_ag][ss.h_ds]);
      for(int g = 0; g < ss.ng; g++)
        for(int ha = 0; ha < ss.h_ag; ha++){
          for(int hm = 0; hm < ss.h_ds; hm++){

            grad[g][ha][hm] = 0.0;

            double cd4mx_scale = 1.0;
            if(artp.scale_cd4_mort & (t >= ss.t_art_start) & (hm >= md.everARTelig_idx)){
              double artpop_hahm = 0.0;
              for(int hu = 0; hu < ss.h_ts; hu++)
                artpop_hahm += md.artpop[t][g][ha][hm][hu];
              cd4mx_scale = md.hivpop[t][g][ha][hm] / (md.hivpop[t][g][ha][hm] + artpop_hahm);
            }

            double deaths = cd4mx_scale * nhp.cd4_mort[g][ha][hm] * md.hivpop[t][g][ha][hm];
            hivdeaths_ha[g][ha] += ss.dt * deaths;
            grad[g][ha][hm] = -deaths;
          }
          for(int hm = 1; hm < ss.h_ds; hm++){
            grad[g][ha][hm-1] -= nhp.cd4_prog[g][ha][hm-1] * md.hivpop[t][g][ha][hm-1];
            grad[g][ha][hm] += nhp.cd4_prog[g][ha][hm-1] * md.hivpop[t][g][ha][hm-1];
          }
        }

      if(incp.eppmod != EPP_DIRECTINCID){
        // incidence

        // calculate r(t)
        if(incp.eppmod == EPP_RSPLINE)
          md.rvec_ts[ts] = incp.rspline_rvec[ts];
        else
          md.rvec_ts[ts] = calc_rtrend_rt(md.pop, incp.rtrend_tstab, incp.rtrend_beta, incp.rtrend_r0,
                                          ss.proj_start + 0.5 + ts * 1.0 / ss.hiv_steps_per_year,
                                          ss.proj_start + 0.5 + incp.ts_epidemic_start * 1.0 / ss.hiv_steps_per_year,
                                          ss.dt, t, hts,
                                          md.rvec_ts[ts-1], &md.prevlast, &md.prevcurr);

        // calculate new infections by sex and age
        double infections_ts[NG][pAG];
        calc_infections_eppspectrum2(md.pop, md.hivpop, md.artpop,
                                     md.rvec_ts[ts], incp.rel_infect_art, (ts == incp.ts_epidemic_start) ? incp.iota : 0.0,
                                     incp.incrr_sex, incp.incrr_age,
                                     ss.t_art_start, ss.dt, t, hts, ss.h_ag_start.data(), ss.h_ag_span,
                                     &md.prevcurr, &md.incid15to49_ts[ts], infections_ts);

        md.prev15to49_ts[ts] = md.prevcurr;

        // add new infections to HIV population
        for(int g = 0; g < ss.ng; g++){
          int a = 0;
          for(int ha = 0; ha < ss.h_ag; ha++){
            double infections_a, infections_ha = 0.0;
            for(int i = 0; i < ss.h_ag_span[ha]; i++){
              infections_ha += infections_a = infections_ts[g][a];
              md.infections[t][g][a] += ss.dt*infections_a;
              md.pop[t][ss.i_hivn][g][a] -= ss.dt*infections_a;
              md.pop[t][ss.i_hivp][g][a] += ss.dt*infections_a;
              a++;
            }
	    if(ha < ss.i_h_ag_15to49+ss.r_h_ag_15to49)
              md.incid15to49[t] += ss.dt*infections_ha;

            // add infections to grad hivpop
            for(int hm = 0; hm < ss.h_ds; hm++)
              grad[g][ha][hm] += infections_ha * nhp.cd4_initdist[g][ha][hm];
          }
        }
      }

      // ART progression, mortality, and initiation
      if(t >= ss.t_art_start){

	multi_array<double, 4> gradART(extents[ss.ng][ss.h_ag][ss.h_ds][ss.h_ts]);

        // progression and mortality
        for(int g = 0; g < ss.ng; g++)
          for(int ha = 0; ha < ss.h_ag; ha++)
            for(int hm = md.everARTelig_idx; hm < ss.h_ds; hm++){

              for(int hu = 0; hu < ss.h_ts; hu++){
                double deaths = nhp.art_mort[g][ha][hm][hu] * nhp.artmx_timerr[t][hu] * md.artpop[t][g][ha][hm][hu];
                hivdeaths_ha[g][ha] += ss.dt*deaths;
                gradART[g][ha][hm][hu] = -deaths;
              }

              gradART[g][ha][hm][nhp.art0mos] += -nhp.art_stage_prog_rate * md.artpop[t][g][ha][hm][nhp.art0mos];
              gradART[g][ha][hm][nhp.art6mos] += nhp.art_stage_prog_rate * md.artpop[t][g][ha][hm][nhp.art0mos] - nhp.art_stage_prog_rate * md.artpop[t][g][ha][hm][nhp.art6mos];
              gradART[g][ha][hm][nhp.art1yr] += nhp.art_stage_prog_rate * md.artpop[t][g][ha][hm][nhp.art6mos];

              // ART dropout
              if(artp.art_dropout[t] > 0)
                for(int hu = 0; hu < ss.h_ts; hu++){
                  grad[g][ha][hm] += artp.art_dropout[t] * md.artpop[t][g][ha][hm][hu];
                  gradART[g][ha][hm][hu] -= artp.art_dropout[t] * md.artpop[t][g][ha][hm][hu];
                }

            }


        // ART initiation
        for(int g = 0; g < ss.ng; g++){

	  multi_array<double, 2> artelig_hahm(extents[ss.r_h_ag_15plus][ss.h_ds]);
	  double Xart_15plus = 0.0, Xartelig_15plus = 0.0, expect_mort_artelig15plus = 0.0;
          for(int ha = ss.i_h_ag_15plus; ha < ss.h_ag; ha++){
            for(int hm = md.everARTelig_idx; hm < ss.h_ds; hm++){
              if(hm >= anyelig_idx){
                double prop_elig = (hm >= cd4elig_idx) ? 1.0 : (hm >= nhp.i_h_cd4_350) ? 1.0 - (1.0-artp.special_pop_percelig[t])*(1.0-artp.who34percelig) : artp.special_pop_percelig[t];
                Xartelig_15plus += artelig_hahm[ha-ss.i_h_ag_15plus][hm] = prop_elig * md.hivpop[t][g][ha][hm] ;
                expect_mort_artelig15plus += nhp.cd4_mort[g][ha][hm] * artelig_hahm[ha-ss.i_h_ag_15plus][hm];
              }
              for(int hu = 0; hu < ss.h_ts; hu++)
                Xart_15plus += md.artpop[t][g][ha][hm][hu] + ss.dt * gradART[g][ha][hm][hu];
            }

            // if artp.preg_women_artelig, add pregnant women to artelig_hahm population
            if(g == ss.i_female & artp.preg_women_artelig[t] > 0 & ha < ss.r_h_ag_fert){
              double frr_pop_ha = 0;
              for(int a =  ss.h_ag_start[ha]; a < ss.h_ag_start[ha]+ss.h_ag_span[ha]; a++)
                frr_pop_ha += md.pop[t][ss.i_hivn][g][a]; // add HIV- population
              for(int hm = 0; hm < ss.h_ds; hm++){
                frr_pop_ha += nhp.frr_cd4[t][ha-ss.i_h_ag_fert][hm] * md.hivpop[t][g][ha][hm];
                for(int hu = 0; hu < ss.h_ts; hu++)
                  frr_pop_ha += nhp.frr_art[t][ha-ss.i_h_ag_fert][hm][hu] * md.artpop[t][g][ha][hm][hu];
              }
              for(int hm = anyelig_idx; hm < cd4elig_idx; hm++){
                double pw_elig_hahm = md.births_by_ha[ha-ss.i_h_ag_fert] * nhp.frr_cd4[t][ha-ss.i_h_ag_fert][hm] * md.hivpop[t][g][ha][hm] / frr_pop_ha;
                artelig_hahm[ha-ss.i_h_ag_15plus][hm] += pw_elig_hahm;
                Xartelig_15plus += pw_elig_hahm;
                expect_mort_artelig15plus += nhp.cd4_mort[g][ha][hm] * pw_elig_hahm;
              }
            }
          } // loop over ha

            // calculate number on ART at end of ts, based on number or percent
          double artnum_hts = 0.0;
          if(ss.dt*(hts+1) < 0.5){
            if(!artp.art15plus_isperc[t-2][g] & !artp.art15plus_isperc[t-1][g]){ // both numbers
              artnum_hts = (0.5-ss.dt*(hts+1))*artp.artnum15plus[t-2][g] + (ss.dt*(hts+1)+0.5)*artp.artnum15plus[t-1][g];
            } else if(artp.art15plus_isperc[t-2][g] & artp.art15plus_isperc[t-1][g]){ // both percentages
              double artcov_hts = (0.5-ss.dt*(hts+1))*artp.artnum15plus[t-2][g] + (ss.dt*(hts+1)+0.5)*artp.artnum15plus[t-1][g];
              artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
            } else if((!artp.art15plus_isperc[t-2][g]) & artp.art15plus_isperc[t-1][g]){ // transition from number to percentage
              double curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus);
              double artcov_hts = curr_coverage + (artp.artnum15plus[t-1][g] - curr_coverage) * ss.dt / (0.5-ss.dt*hts);
              artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
            }
          } else {
            if(!artp.art15plus_isperc[t-1][g] & !artp.art15plus_isperc[t][g]){ // both numbers
              artnum_hts = (1.5-ss.dt*(hts+1))*artp.artnum15plus[t-1][g] + (ss.dt*(hts+1)-0.5)*artp.artnum15plus[t][g];
            } else if(artp.art15plus_isperc[t-1][g] & artp.art15plus_isperc[t][g]){ // both percentages
              double artcov_hts = (1.5-ss.dt*(hts+1))*artp.artnum15plus[t-1][g] + (ss.dt*(hts+1)-0.5)*artp.artnum15plus[t][g];
              artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
            } else if((!artp.art15plus_isperc[t-1][g]) & artp.art15plus_isperc[t][g]){ // transition from number to percentage
              double curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus);
              double artcov_hts = curr_coverage + (artp.artnum15plus[t][g] - curr_coverage) * ss.dt / (1.5-ss.dt*hts);
              artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
            }
          }

          double artinit_hts = artnum_hts > Xart_15plus ? artnum_hts - Xart_15plus : 0;

          // median CD4 at initiation inputs
          if(artp.med_cd4init_input[t]){

            const int CD4_LOW_LIM[hDS] = {500, 350, 250, 200, 100, 50, 0};
            const int CD4_UPP_LIM[hDS] = {1000, 500, 350, 250, 200, 100, 50};

            int medcd4_idx = artp.med_cd4init_cat[t] - 1; // -1 for 0-based indexing vs. 1-based in R
            double medcat_propbelow = (artp.median_cd4init[t] - CD4_LOW_LIM[medcd4_idx]) / (CD4_UPP_LIM[medcd4_idx] - CD4_LOW_LIM[medcd4_idx]);

            double elig_below = 0.0, elig_above = 0.0;
            for(int ha = ss.i_h_ag_15plus; ha < ss.h_ag; ha++){
              for(int hm = anyelig_idx; hm < medcd4_idx; hm++)
                elig_above += artelig_hahm[ha-ss.i_h_ag_15plus][hm];
              elig_above += (1.0 - medcat_propbelow) * artelig_hahm[ha-ss.i_h_ag_15plus][medcd4_idx];
              elig_below += medcat_propbelow * artelig_hahm[ha-ss.i_h_ag_15plus][medcd4_idx];
              for(int hm = medcd4_idx+1; hm < ss.h_ds; hm++)
                elig_below += artelig_hahm[ha-ss.i_h_ag_15plus][hm];
            }

            double initprob_below = (elig_below > artinit_hts * 0.5) ? artinit_hts * 0.5 / elig_below : 1.0;
            double initprob_above = (elig_above > artinit_hts * 0.5) ? artinit_hts * 0.5 / elig_above : 1.0;
            double initprob_medcat = initprob_below * medcat_propbelow + initprob_above * (1.0-medcat_propbelow);

            for(int ha = ss.i_h_ag_15plus; ha < ss.h_ag; ha++)
              for(int hm = anyelig_idx; hm < ss.h_ds; hm++){
                double artinit_hahm;
                if(hm < medcd4_idx)
                  artinit_hahm = artelig_hahm[ha-ss.i_h_ag_15plus][hm] * initprob_above;
                else if(hm == medcd4_idx)
                  artinit_hahm = artelig_hahm[ha-ss.i_h_ag_15plus][hm] * initprob_medcat;
                if(hm > medcd4_idx)
                  artinit_hahm = artelig_hahm[ha-ss.i_h_ag_15plus][hm] * initprob_below;
                if(artinit_hahm > md.hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm])
                  artinit_hahm = md.hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm];
                grad[g][ha][hm] -= artinit_hahm / ss.dt;
                gradART[g][ha][hm][nhp.art0mos] += artinit_hahm / ss.dt;
              }

          } else if(artp.art_alloc_method == 4) {  // lowest CD4 first

            for(int hm = ss.h_ds-1; hm >= anyelig_idx; hm--){
              double artelig_hm = 0;
              for(int ha = ss.i_h_ag_15plus; ha < ss.h_ag; ha++)
                artelig_hm += artelig_hahm[ha-ss.i_h_ag_15plus][hm];
              double init_prop = (artelig_hm == 0 | artinit_hts > artelig_hm) ? 1.0 : artinit_hts / artelig_hm;

              for(int ha = ss.i_h_ag_15plus; ha < ss.h_ag; ha++){
                double artinit_hahm = init_prop * artelig_hahm[ha-ss.i_h_ag_15plus][hm];

                if(artinit_hahm > md.hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm])
                  artinit_hahm = md.hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm];

                grad[g][ha][hm] -= artinit_hahm / ss.dt;
                gradART[g][ha][hm][nhp.art0mos] += artinit_hahm / ss.dt;
              }
              if(init_prop < 1.0)
                break;
              artinit_hts -= init_prop * artelig_hm;
            }

          } else { // Use mixture of eligibility and expected mortality for initiation distribution

            for(int ha = ss.i_h_ag_15plus; ha < ss.h_ag; ha++)
              for(int hm = anyelig_idx; hm < ss.h_ds; hm++){
                double artinit_hahm = artinit_hts * artelig_hahm[ha-ss.i_h_ag_15plus][hm] * ((1.0 - artp.art_alloc_mxweight)/Xartelig_15plus + artp.art_alloc_mxweight * nhp.cd4_mort[g][ha][hm] / expect_mort_artelig15plus);
                if(artinit_hahm > artelig_hahm[ha-ss.i_h_ag_15plus][hm])
                  artinit_hahm = artelig_hahm[ha-ss.i_h_ag_15plus][hm];
                if(artinit_hahm > md.hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm])
                  artinit_hahm = md.hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm];
                grad[g][ha][hm] -= artinit_hahm / ss.dt;
                gradART[g][ha][hm][nhp.art0mos] += artinit_hahm / ss.dt;
              }
          }
        }

        for(int g = 0; g < ss.ng; g++)
          for(int ha = 0; ha < ss.h_ag; ha++)
            for(int hm = md.everARTelig_idx; hm < ss.h_ds; hm++)
              for(int hu = 0; hu < ss.h_ts; hu++)
                md.artpop[t][g][ha][hm][hu] += ss.dt*gradART[g][ha][hm][hu];

      } // if(t >= ss.t_art_start)

      for(int g = 0; g < ss.ng; g++)
        for(int ha = 0; ha < ss.h_ag; ha++)
          for(int hm = 0; hm < ss.h_ds; hm++)
            md.hivpop[t][g][ha][hm] += ss.dt*grad[g][ha][hm];


      // remove hivdeaths from pop
      for(int g = 0; g < ss.ng; g++){

        // sum HIV+ population size in each hivpop age group
	std::vector<double> hivpop_ha(ss.h_ag);
        int a = 0;
        for(int ha = 0; ha < ss.h_ag; ha++){
          hivpop_ha[ha] = 0.0;
          for(int i = 0; i < ss.h_ag_span[ha]; i++){
            hivpop_ha[ha] += md.pop[t][ss.i_hivp][g][a];
            a++;
          }
        }

        // remove hivdeaths proportionally to age-distribution within each age group
        a = 0;
        for(int ha = 0; ha < ss.h_ag; ha++){
          if(hivpop_ha[ha] > 0){
            double hivqx_ha = hivdeaths_ha[g][ha] / hivpop_ha[ha];
            for(int i = 0; i < ss.h_ag_span[ha]; i++){
              md.hivdeaths[t][g][a] += md.pop[t][ss.i_hivp][g][a] * hivqx_ha;
              md.pop[t][ss.i_hivp][g][a] *= (1.0-hivqx_ha);
              a++;
            }
          } else {
            a += ss.h_ag_span[ha];
          }  // end if(pop_ha[ha] > 0)
        }
      }



    } // loop HIVSTEPS_PER_YEAR



    if(incp.eppmod == EPP_DIRECTINCID){
      // Calculating new infections once per year (like Spectrum)

      double Xhivp = 0.0;
      std::vector<double> Xhivn(ss.ng), Xhivn_incagerr(ss.ng);

      for(int g = 0; g < ss.ng; g++){
        Xhivn[g] = 0.0;
        Xhivn_incagerr[g] = 0.0;
        for(int a = incp.i_p_ag_incidpop; a < incp.i_p_ag_incidpop + incp.r_p_ag_incidpop; a++){
          Xhivp += md.pop[t-1][ss.i_hivp][g][a];
          Xhivn[g] += md.pop[t-1][ss.i_hivn][g][a];
          Xhivn_incagerr[g] += incp.incrr_age[t][g][a] * md.pop[t-1][ss.i_hivn][g][a];
        }
      }
      // double prev_i = Xhivp / (Xhivn[ss.i_male] + Xhivn[ss.i_female] + Xhivp);
      // double incrate15to49_i = (prev15to49[t] - prev_i)/(1.0 - prev_i);
      double incrate_i = incp.incidinput[t];
      std::vector<double> incrate_g(ss.ng);
      incrate_g[ss.i_male] = incrate_i * (Xhivn[ss.i_male]+Xhivn[ss.i_female]) / (Xhivn[ss.i_male] + incp.incrr_sex[t]*Xhivn[ss.i_female]);
      incrate_g[ss.i_female] = incrate_i * incp.incrr_sex[t]*(Xhivn[ss.i_male]+Xhivn[ss.i_female]) / (Xhivn[ss.i_male] + incp.incrr_sex[t]*Xhivn[ss.i_female]);

      for(int g = 0; g < ss.ng; g++){
        int a = 0;
        for(int ha = 0; ha < ss.h_ag; ha++){
          double infections_a, infections_ha = 0.0;
          for(int i = 0; i < ss.h_ag_span[ha]; i++){
            infections_ha += infections_a = md.pop[t-1][ss.i_hivn][g][a] * incrate_g[g] * incp.incrr_age[t][g][a] * Xhivn[g] / Xhivn_incagerr[g];
            md.infections[t][g][a] += infections_a;
            md.pop[t][ss.i_hivn][g][a] -= infections_a;
            md.pop[t][ss.i_hivp][g][a] += infections_a;
            a++;
          }
          if(ha < ss.i_h_ag_15to49+ss.r_h_ag_15to49)
            md.incid15to49[t] += infections_ha;

          // add infections to hivpop
          for(int hm = 0; hm < ss.h_ds; hm++)
            md.hivpop[t][g][ha][hm] += infections_ha * nhp.cd4_initdist[g][ha][hm];
        }
      }
    }


    // md.adjust_population(t);

    // adjust population to match target population
    if(demp.flag_popadjust){
      for(int g = 0; g < ss.ng; g++){
        int a = 0;
        for(int ha = 0; ha < ss.h_ag; ha++){
          double popadj_ha = 0, hivpop_ha = 0;
          for(int i = 0; i < ss.h_ag_span[ha]; i++){

            hivpop_ha += md.pop[t][ss.i_hivp][g][a];

            double popadjrate_a = md.popadjust[t][g][a] = demp.targetpop[t][g][a] / (md.pop[t][ss.i_hivn][g][a] + md.pop[t][ss.i_hivp][g][a]);
            md.pop[t][ss.i_hivn][g][a] *= popadjrate_a;
            double hpopadj_a = (popadjrate_a-1.0) * md.pop[t][ss.i_hivp][g][a];
            popadj_ha += hpopadj_a;
            md.pop[t][ss.i_hivp][g][a] += hpopadj_a;
            a++;
          }

          // population adjustment for hivpop
          double popadjrate_ha = hivpop_ha > 0 ? popadj_ha / hivpop_ha : 0.0;
          for(int hm = 0; hm < ss.h_ds; hm++){
            md.hivpop[t][g][ha][hm] *= 1+popadjrate_ha;
            if(t >= ss.t_art_start)
              for(int hu = 0; hu < ss.h_ts; hu++)
                md.artpop[t][g][ha][hm][hu] *= 1+popadjrate_ha;
          } // loop over hm
        } // loop over ha
      } // loop over g
    } // if(flag_popadjust)

    // prevalence among pregnant women

    double hivbirths = 0;
    for(int ha = ss.i_h_ag_fert; ha < ss.i_h_ag_fert+ss.r_h_ag_fert; ha++){
      double hivn_ha = 0, frr_hivpop_ha = 0;
      for(int a =  ss.h_ag_start[ha]; a < ss.h_ag_start[ha]+ss.h_ag_span[ha]; a++)
        hivn_ha += (md.pop[t-1][ss.i_hivn][ss.i_female][a] + md.pop[t][ss.i_hivn][ss.i_female][a])/2;
      for(int hm = 0; hm < ss.h_ds; hm++){
        frr_hivpop_ha += nhp.frr_cd4[t][ha-ss.i_h_ag_fert][hm] * (md.hivpop[t-1][ss.i_female][ha][hm]+md.hivpop[t][ss.i_female][ha][hm])/2;
        if(t == ss.t_art_start)
          for(int hu = 0; hu < ss.h_ts; hu++)
            frr_hivpop_ha += nhp.frr_art[t][ha-ss.i_h_ag_fert][hm][hu] * md.artpop[t][ss.i_female][ha][hm][hu]/2;
        else if(t > ss.t_art_start)
          for(int hu = 0; hu < ss.h_ts; hu++)
            frr_hivpop_ha += nhp.frr_art[t][ha-ss.i_h_ag_fert][hm][hu] * (md.artpop[t-1][ss.i_female][ha][hm][hu]+md.artpop[t][ss.i_female][ha][hm][hu])/2;
      }
      hivbirths += md.births_by_ha[ha-ss.i_h_ag_fert] * frr_hivpop_ha / (hivn_ha + frr_hivpop_ha);
    }

    md.pregprev[t] = hivbirths/md.births;
    if(t + ss.pop_age_start < ss.proj_years)
      md.pregprevlag[t + ss.pop_age_start-1] = md.pregprev[t];

    md.incid15to49[t] /= md.hivn15to49;

    // prevalence 15 to 49
    md.hivn15to49 = 0;
    for(int g = 0; g < ss.ng; g++)
      for(int a = ss.i_p_ag_15to49; a < ss.i_p_ag_15to49 + ss.r_p_ag_15to49; a++){
        md.hivn15to49 += md.pop[t][ss.i_hivn][g][a];
        md.prev15to49[t] += md.pop[t][ss.i_hivp][g][a];
      }
    md.prev15to49[t] /= (md.hivn15to49 + md.prev15to49[t]);

  }
}


void simulate_eppasm5(Model &md,
		      const SsDim &ss,
		      const DemogParam &demp,
		      const PaediatricHivParam &paedhp,
		      const NaturalHistoryParam &nhp,
		      const ArtData &artp,
		      const IncidenceParam &incp) {

  md.initialise();


  ////////////////////////////////////
  ////  do population projection  ////
  ////////////////////////////////////

  for(int t = 1; t < ss.sim_years; t++){

    md.ageing(t);
    md.pop_entrants(t);
    md.death_and_migration(t);
    md.fertility(t);
    
    ////////////////////////////////
    ////  HIV model simulation  ////
    ////////////////////////////////


    for(int hts = 0; hts < ss.hiv_steps_per_year; hts++){

      md.initialise_hiv_ts(t); 
      md.hiv_progression_mortality(t);
      
      if(incp.eppmod != EPP_DIRECTINCID){

	int ts = (t-1)*ss.hiv_steps_per_year + hts;

        // calculate r(t)
        if(incp.eppmod == EPP_RSPLINE)
          md.rvec_ts[ts] = incp.rspline_rvec[ts];
        else  // r-trend
          md.rvec_ts[ts] = calc_rtrend_rt(md.pop, incp.rtrend_tstab, incp.rtrend_beta, incp.rtrend_r0,
                                          ss.proj_start + 0.5 + ts * 1.0 / ss.hiv_steps_per_year,
                                          ss.proj_start + 0.5 + incp.ts_epidemic_start * 1.0 / ss.hiv_steps_per_year,
                                          ss.dt, t, hts,
                                          md.rvec_ts[ts-1], &md.prevlast, &md.prevcurr);

        // calculate new infections by sex and age
        double infections_ts[NG][pAG];
        calc_infections_eppspectrum2(md.pop, md.hivpop, md.artpop,
                                     md.rvec_ts[ts], incp.rel_infect_art, (ts == incp.ts_epidemic_start) ? incp.iota : 0.0,
                                     incp.incrr_sex, incp.incrr_age,
                                     ss.t_art_start, ss.dt, t, hts, ss.h_ag_start.data(), ss.h_ag_span,
                                     &md.prevcurr, &md.incid15to49_ts[ts], infections_ts);

        md.prev15to49_ts[ts] = md.prevcurr;

        // add new infections to HIV population
        for(int g = 0; g < ss.ng; g++){
          int a = 0;
          for(int ha = 0; ha < ss.h_ag; ha++){
            double infections_a, infections_ha = 0.0;
            for(int i = 0; i < ss.h_ag_span[ha]; i++){
              infections_ha += infections_a = infections_ts[g][a];
              md.infections[t][g][a] += ss.dt*infections_a;
              md.pop[t][ss.i_hivn][g][a] -= ss.dt*infections_a;
              md.pop[t][ss.i_hivp][g][a] += ss.dt*infections_a;
              a++;
            }
	    if(ha < ss.i_h_ag_15to49+ss.r_h_ag_15to49)
              md.incid15to49[t] += ss.dt*infections_ha;

            // add infections to grad hivpop
            for(int hm = 0; hm < ss.h_ds; hm++)
              md.grad[g][ha][hm] += infections_ha * nhp.cd4_initdist[g][ha][hm];
          }
        }
      }

      md.art_initiation(t, hts);

      md.do_hts_step(t);
      md.remove_hiv_deaths(t);

    } // loop HIVSTEPS_PER_YEAR



    if(incp.eppmod == EPP_DIRECTINCID){
      // Calculating new infections once per year (like Spectrum)

      double Xhivp = 0.0;
      std::vector<double> Xhivn(ss.ng), Xhivn_incagerr(ss.ng);

      for(int g = 0; g < ss.ng; g++){
        Xhivn[g] = 0.0;
        Xhivn_incagerr[g] = 0.0;
        for(int a = incp.i_p_ag_incidpop; a < incp.i_p_ag_incidpop + incp.r_p_ag_incidpop; a++){
          Xhivp += md.pop[t-1][ss.i_hivp][g][a];
          Xhivn[g] += md.pop[t-1][ss.i_hivn][g][a];
          Xhivn_incagerr[g] += incp.incrr_age[t][g][a] * md.pop[t-1][ss.i_hivn][g][a];
        }
      }
      // double prev_i = Xhivp / (Xhivn[ss.i_male] + Xhivn[ss.i_female] + Xhivp);
      // double incrate15to49_i = (prev15to49[t] - prev_i)/(1.0 - prev_i);
      double incrate_i = incp.incidinput[t];
      std::vector<double> incrate_g(ss.ng);
      incrate_g[ss.i_male] = incrate_i * (Xhivn[ss.i_male]+Xhivn[ss.i_female]) / (Xhivn[ss.i_male] + incp.incrr_sex[t]*Xhivn[ss.i_female]);
      incrate_g[ss.i_female] = incrate_i * incp.incrr_sex[t]*(Xhivn[ss.i_male]+Xhivn[ss.i_female]) / (Xhivn[ss.i_male] + incp.incrr_sex[t]*Xhivn[ss.i_female]);

      for(int g = 0; g < ss.ng; g++){
        int a = 0;
        for(int ha = 0; ha < ss.h_ag; ha++){
          double infections_a, infections_ha = 0.0;
          for(int i = 0; i < ss.h_ag_span[ha]; i++){
            infections_ha += infections_a = md.pop[t-1][ss.i_hivn][g][a] * incrate_g[g] * incp.incrr_age[t][g][a] * Xhivn[g] / Xhivn_incagerr[g];
            md.infections[t][g][a] += infections_a;
            md.pop[t][ss.i_hivn][g][a] -= infections_a;
            md.pop[t][ss.i_hivp][g][a] += infections_a;
            a++;
          }
          if(ha < ss.i_h_ag_15to49+ss.r_h_ag_15to49)
            md.incid15to49[t] += infections_ha;

          // add infections to hivpop
          for(int hm = 0; hm < ss.h_ds; hm++)
            md.hivpop[t][g][ha][hm] += infections_ha * nhp.cd4_initdist[g][ha][hm];
        }
      }
    }

    md.adjust_population(t);

    // prevalence among pregnant women

    double hivbirths = 0;
    for(int ha = ss.i_h_ag_fert; ha < ss.i_h_ag_fert+ss.r_h_ag_fert; ha++){
      double hivn_ha = 0, frr_hivpop_ha = 0;
      for(int a =  ss.h_ag_start[ha]; a < ss.h_ag_start[ha]+ss.h_ag_span[ha]; a++)
        hivn_ha += (md.pop[t-1][ss.i_hivn][ss.i_female][a] + md.pop[t][ss.i_hivn][ss.i_female][a])/2;
      for(int hm = 0; hm < ss.h_ds; hm++){
        frr_hivpop_ha += nhp.frr_cd4[t][ha-ss.i_h_ag_fert][hm] * (md.hivpop[t-1][ss.i_female][ha][hm]+md.hivpop[t][ss.i_female][ha][hm])/2;
        if(t == ss.t_art_start)
          for(int hu = 0; hu < ss.h_ts; hu++)
            frr_hivpop_ha += nhp.frr_art[t][ha-ss.i_h_ag_fert][hm][hu] * md.artpop[t][ss.i_female][ha][hm][hu]/2;
        else if(t > ss.t_art_start)
          for(int hu = 0; hu < ss.h_ts; hu++)
            frr_hivpop_ha += nhp.frr_art[t][ha-ss.i_h_ag_fert][hm][hu] * (md.artpop[t-1][ss.i_female][ha][hm][hu]+md.artpop[t][ss.i_female][ha][hm][hu])/2;
      }
      hivbirths += md.births_by_ha[ha-ss.i_h_ag_fert] * frr_hivpop_ha / (hivn_ha + frr_hivpop_ha);
    }

    md.pregprev[t] = hivbirths/md.births;
    if(t + ss.pop_age_start < ss.proj_years)
      md.pregprevlag[t + ss.pop_age_start-1] = md.pregprev[t];

    md.incid15to49[t] /= md.hivn15to49;

    // prevalence 15 to 49
    md.hivn15to49 = 0;
    for(int g = 0; g < ss.ng; g++)
      for(int a = ss.i_p_ag_15to49; a < ss.i_p_ag_15to49 + ss.r_p_ag_15to49; a++){
        md.hivn15to49 += md.pop[t][ss.i_hivn][g][a];
        md.prev15to49[t] += md.pop[t][ss.i_hivp][g][a];
      }
    md.prev15to49[t] /= (md.hivn15to49 + md.prev15to49[t]);

  }
}
