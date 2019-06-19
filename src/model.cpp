#include <boost/multi_array.hpp>
#include "model.hpp"


void Model::initialise() {

  for(int g = 0; g < ss.ng; g++)
    for(int a = 0; a < ss.p_ag; a++){
      pop[0][ss.i_hivn][g][a] = demp.basepop[g][a];
      if((a >= ss.i_p_ag_15to49) & (a < ss.i_p_ag_15to49 + ss.r_p_ag_15to49))
        hivn15to49 += demp.basepop[g][a];
    }
}

void Model::ageing(int t) {

  // age the population one year
  for(int m = 0; m < ss.p_ds; m++)
    for(int g = 0; g < ss.ng; g++){
      for(int a = 1; a < ss.p_ag; a++)
        pop[t][m][g][a] = pop[t-1][m][g][a-1];
      pop[t][m][g][ss.p_ag-1] += pop[t-1][m][g][ss.p_ag-1]; // open age group
    }

  // calculate proportion to age within each HIV age group
  boost::multi_array<double, 2> hiv_ag_prob(boost::extents[ss.ng][ss.h_ag]);
  for(int g = 0; g < ss.ng; g++){
    int a = 0;
    for(int ha = 0; ha < (ss.h_ag-1); ha++){
      hiv_ag_prob[g][ha] = 0;
      for(int i = 0; i < ss.h_ag_span[ha]; i++){
        hiv_ag_prob[g][ha] += pop[t-1][ss.i_hivp][g][a];
        a++;
      }
      hiv_ag_prob[g][ha] = (hiv_ag_prob[g][ha] > 0) ? pop[t-1][ss.i_hivp][g][a-1] / hiv_ag_prob[g][ha] : 0;
    }
    hiv_ag_prob[g][ss.h_ag-1] = 0.0; // no one ages out of the open-ended age group
  }

  for(int g = 0; g < ss.ng; g++){

    // youngest age group: only ageing out; ageing in entered below
    for(int hm = 0; hm < ss.h_ds; hm++){
      hivpop[t][g][0][hm] = (1-hiv_ag_prob[g][0]) * hivpop[t-1][g][0][hm];
      if(t > ss.t_art_start){
        for(int hu = 0; hu < ss.h_ts; hu++)
          artpop[t][g][0][hm][hu] = (1-hiv_ag_prob[g][0]) * artpop[t-1][g][0][hm][hu];
      }
    }

    for(int ha = 1; ha < ss.h_ag; ha++)
      for(int hm = 0; hm < ss.h_ds; hm++){
        hivpop[t][g][ha][hm] = (1-hiv_ag_prob[g][ha]) * hivpop[t-1][g][ha][hm] + hiv_ag_prob[g][ha-1] * hivpop[t-1][g][ha-1][hm];
        if(t > ss.t_art_start)
          for(int hu = 0; hu < ss.h_ts; hu++)
            artpop[t][g][ha][hm][hu] = (1-hiv_ag_prob[g][ha]) * artpop[t-1][g][ha][hm][hu] + hiv_ag_prob[g][ha-1] * artpop[t-1][g][ha-1][hm][hu];
      }
  }

}

void Model::death_and_migration(int t) {

  for(int g = 0; g < ss.ng; g++){
    int a = 0;
    for(int ha = 0; ha < ss.h_ag; ha++){
      double deathsmig_ha = 0, hivpop_ha = 0;
      for(int i = 0; i < ss.h_ag_span[ha]; i++){

        hivpop_ha += pop[t][ss.i_hivp][g][a];

        // non-HIV mortality
        double qx = 1.0 - demp.Sx[t][g][a];
        double ndeaths_a = pop[t][ss.i_hivn][g][a] * qx;
        pop[t][ss.i_hivn][g][a] -= ndeaths_a; // survival HIV- population
        double hdeaths_a = pop[t][ss.i_hivp][g][a] * qx;
        deathsmig_ha -= hdeaths_a;
        pop[t][ss.i_hivp][g][a] -= hdeaths_a;   // survival HIV+ population
        natdeaths[t][g][a] = ndeaths_a + hdeaths_a;

        // net migration
        double migrate_a = demp.netmigr[t][g][a] * (1+demp.Sx[t][g][a])/2.0 / (pop[t][ss.i_hivn][g][a] + pop[t][ss.i_hivp][g][a]);
        pop[t][ss.i_hivn][g][a] *= 1+migrate_a;
        double hmig_a = migrate_a * pop[t][ss.i_hivp][g][a];
        deathsmig_ha += hmig_a;
        pop[t][ss.i_hivp][g][a] += hmig_a;

        a++;
      }

      // migration and deaths for hivpop
      double deathmigrate_ha = hivpop_ha > 0 ? deathsmig_ha / hivpop_ha : 0.0;
      for(int hm = 0; hm < ss.h_ds; hm++){
        hivpop[t][g][ha][hm] *= 1+deathmigrate_ha;
        if(t > ss.t_art_start)
          for(int hu = 0; hu < ss.h_ts; hu++)
            artpop[t][g][ha][hm][hu] *= 1+deathmigrate_ha;
      } // loop over hm
    } // loop over ha
  } // loop over g

}


void Model::pop_entrants(int t) {

  // add lagged births to youngest age group
  for(int g = 0; g < ss.ng; g++){

    double paedsurv_g;
    double entrant_prev;

    if(paedhp.use_entrantprev)
      entrant_prev = paedhp.entrantprev[t][g];
    else
      entrant_prev = pregprevlag[t-1] * paedhp.verttrans_lag[t-1] * paedhp.paedsurv_lag[t-1];

    if(demp.flag_popadjust){
      pop[t][ss.i_hivn][g][0] = demp.entrantpop[t-1][g] * (1.0-entrant_prev);
      paedsurv_g = demp.entrantpop[t-1][g] * entrant_prev;
    } else {
      pop[t][ss.i_hivn][g][0] = birthslag[t-1][g] * demp.cumsurv[t-1][g] * (1.0-entrant_prev / paedhp.paedsurv_lag[t-1]) + demp.cumnetmigr[t-1][g] * (1.0-pregprevlag[t-1] * paedhp.netmig_hivprob);
      paedsurv_g = birthslag[t-1][g] * demp.cumsurv[t-1][g] * entrant_prev + demp.cumnetmigr[t-1][g] * entrant_prev;
    }

    pop[t][ss.i_hivp][g][0] = paedsurv_g;

    entrantprev[t] = (pop[t][ss.i_hivp][ss.i_male][0] + pop[t][ss.i_hivp][ss.i_female][0]) / (pop[t][ss.i_hivn][ss.i_male][0] + pop[t][ss.i_hivn][ss.i_female][0] + pop[t][ss.i_hivp][ss.i_male][0] + pop[t][ss.i_hivp][ss.i_female][0]);

    for(int hm = 0; hm < ss.h_ds; hm++){
      hivpop[t][g][0][hm] += paedsurv_g * paedhp.paedsurv_cd4dist[t][g][hm] * (1.0 - paedhp.entrantartcov[t][g]);
      if(t > ss.t_art_start){
        for(int hu = 0; hu < ss.h_ts; hu++){
          artpop[t][g][0][hm][hu] += paedsurv_g * paedhp.paedsurv_artcd4dist[t][g][hm][hu] * paedhp.entrantartcov[t][g];
        }
      }
    }
  }
}

void Model::fertility(int t) {

  births = 0.0;
  std::fill(births_by_ha.begin(), births_by_ha.end(), 0);
  for(int m = 0; m < ss.p_ds; m++){
    int a = ss.i_p_ag_fert;
    for(int ha = ss.i_h_ag_fert; ha < ss.i_h_ag_fert+ss.r_h_ag_fert; ha++){
      for(int i = 0; i < ss.h_ag_span[ha]; i++){
        births_by_ha[ha-ss.i_h_ag_fert] += (pop[t-1][m][ss.i_female][a] + pop[t][m][ss.i_female][a])/2 * demp.asfr[t][a];
        a++;
      }
    }
  }
  for(int ha = ss.i_h_ag_fert; ha < ss.r_h_ag_fert; ha++)
    births += births_by_ha[ha-ss.i_h_ag_fert];

  if(t + ss.pop_age_start < ss.proj_years)
    for(int g = 0; g < ss.ng; g++)
      birthslag[t + ss.pop_age_start-1][g] = demp.srb[t][g] * births;

}

void Model::adjust_population(int t) {

  if(demp.flag_popadjust){
    for(int g = 0; g < ss.ng; g++){
      int a = 0;
      for(int ha = 0; ha < ss.h_ag; ha++){
        double popadj_ha = 0, hivpop_ha = 0;
        for(int i = 0; i < ss.h_ag_span[ha]; i++){

          hivpop_ha += pop[t][ss.i_hivp][g][a];

          double popadjrate_a = popadjust[t][g][a] = demp.targetpop[t][g][a] / (pop[t][ss.i_hivn][g][a] + pop[t][ss.i_hivp][g][a]);
          pop[t][ss.i_hivn][g][a] *= popadjrate_a;
          double hpopadj_a = (popadjrate_a-1.0) * pop[t][ss.i_hivp][g][a];
          popadj_ha += hpopadj_a;
          pop[t][ss.i_hivp][g][a] += hpopadj_a;
          a++;
        }

        // population adjustment for hivpop
        double popadjrate_ha = hivpop_ha > 0 ? popadj_ha / hivpop_ha : 0.0;
        for(int hm = 0; hm < ss.h_ds; hm++){
          hivpop[t][g][ha][hm] *= 1+popadjrate_ha;
          if(t >= ss.t_art_start)
            for(int hu = 0; hu < ss.h_ts; hu++)
              artpop[t][g][ha][hm][hu] *= 1+popadjrate_ha;
        } // loop over hm
      } // loop over ha
    } // loop over g
  } // if(flag_popadjust)
}


void Model::initialise_hiv_ts(int t) {

  cd4elig_idx = artp.artcd4elig_idx[t] - 1; // -1 for 0-based indexing vs. 1-based in R
  anyelig_idx = (artp.special_pop_percelig[t] > 0 | artp.preg_women_artelig[t] > 0) ? 0 : (artp.who34percelig > 0) ? nhp.i_h_cd4_350 : cd4elig_idx;
  everARTelig_idx = anyelig_idx < everARTelig_idx ? anyelig_idx : everARTelig_idx;

  memset(hivdeaths_ha.data(), 0, hivdeaths_ha.num_elements() * sizeof(double));
  memset(grad.data(), 0, grad.num_elements() * sizeof(double));

  if(t >= ss.t_art_start)
    memset(gradART.data(), 0, gradART.num_elements() * sizeof(double));
}

void Model::hiv_progression_mortality(int t) {

  // untreated population
  for(int g = 0; g < ss.ng; g++)
    for(int ha = 0; ha < ss.h_ag; ha++){
      for(int hm = 0; hm < ss.h_ds; hm++){

        double cd4mx_scale = 1.0;
        if(artp.scale_cd4_mort & (t >= ss.t_art_start) & (hm >= everARTelig_idx)){
          double artpop_hahm = 0.0;
          for(int hu = 0; hu < ss.h_ts; hu++)
            artpop_hahm += artpop[t][g][ha][hm][hu];
          cd4mx_scale = hivpop[t][g][ha][hm] / (hivpop[t][g][ha][hm] + artpop_hahm);
        }

        double deaths = cd4mx_scale * nhp.cd4_mort[g][ha][hm] * hivpop[t][g][ha][hm];
        hivdeaths_ha[g][ha] += ss.dt * deaths;
        grad[g][ha][hm] = -deaths;
      }
      for(int hm = 1; hm < ss.h_ds; hm++){
        grad[g][ha][hm-1] -= nhp.cd4_prog[g][ha][hm-1] * hivpop[t][g][ha][hm-1];
        grad[g][ha][hm] += nhp.cd4_prog[g][ha][hm-1] * hivpop[t][g][ha][hm-1];
      }
    }

  // ART population

  if(t >= ss.t_art_start){

    for(int g = 0; g < ss.ng; g++)
      for(int ha = 0; ha < ss.h_ag; ha++)
        for(int hm = everARTelig_idx; hm < ss.h_ds; hm++){

          for(int hu = 0; hu < ss.h_ts; hu++){
            double deaths = nhp.art_mort[g][ha][hm][hu] * nhp.artmx_timerr[t][hu] * artpop[t][g][ha][hm][hu];
            hivdeaths_ha[g][ha] += ss.dt*deaths;
            gradART[g][ha][hm][hu] = -deaths;
          }

          gradART[g][ha][hm][nhp.art0mos] += -nhp.art_stage_prog_rate * artpop[t][g][ha][hm][nhp.art0mos];
          gradART[g][ha][hm][nhp.art6mos] += nhp.art_stage_prog_rate * artpop[t][g][ha][hm][nhp.art0mos] - nhp.art_stage_prog_rate * artpop[t][g][ha][hm][nhp.art6mos];
          gradART[g][ha][hm][nhp.art1yr] += nhp.art_stage_prog_rate * artpop[t][g][ha][hm][nhp.art6mos];

          // ART dropout
          if(artp.art_dropout[t] > 0)
            for(int hu = 0; hu < ss.h_ts; hu++){
              grad[g][ha][hm] += artp.art_dropout[t] * artpop[t][g][ha][hm][hu];
              gradART[g][ha][hm][hu] -= artp.art_dropout[t] * artpop[t][g][ha][hm][hu];
            }
        }
  }
}

void Model::art_initiation(int t, int hts){

  if(t >= ss.t_art_start){

    for(int g = 0; g < ss.ng; g++){

      boost::multi_array<double, 2> artelig_hahm(boost::extents[ss.r_h_ag_15plus][ss.h_ds]);
      double Xart_15plus = 0.0, Xartelig_15plus = 0.0, expect_mort_artelig15plus = 0.0;
      for(int ha = ss.i_h_ag_15plus; ha < ss.h_ag; ha++){
        for(int hm = everARTelig_idx; hm < ss.h_ds; hm++){
          if(hm >= anyelig_idx){
            double prop_elig = (hm >= cd4elig_idx) ? 1.0 : (hm >= nhp.i_h_cd4_350) ? 1.0 - (1.0-artp.special_pop_percelig[t])*(1.0-artp.who34percelig) : artp.special_pop_percelig[t];
            Xartelig_15plus += artelig_hahm[ha-ss.i_h_ag_15plus][hm] = prop_elig * hivpop[t][g][ha][hm] ;
            expect_mort_artelig15plus += nhp.cd4_mort[g][ha][hm] * artelig_hahm[ha-ss.i_h_ag_15plus][hm];
          }
          for(int hu = 0; hu < ss.h_ts; hu++)
            Xart_15plus += artpop[t][g][ha][hm][hu] + ss.dt * gradART[g][ha][hm][hu];
        }

        // if artp.preg_women_artelig, add pregnant women to artelig_hahm population
        if(g == ss.i_female & artp.preg_women_artelig[t] > 0 & ha < ss.r_h_ag_fert){
          double frr_pop_ha = 0;
          for(int a =  ss.h_ag_start[ha]; a < ss.h_ag_start[ha]+ss.h_ag_span[ha]; a++)
            frr_pop_ha += pop[t][ss.i_hivn][g][a]; // add HIV- population
          for(int hm = 0; hm < ss.h_ds; hm++){
            frr_pop_ha += nhp.frr_cd4[t][ha-ss.i_h_ag_fert][hm] * hivpop[t][g][ha][hm];
            for(int hu = 0; hu < ss.h_ts; hu++)
              frr_pop_ha += nhp.frr_art[t][ha-ss.i_h_ag_fert][hm][hu] * artpop[t][g][ha][hm][hu];
          }
          for(int hm = anyelig_idx; hm < cd4elig_idx; hm++){
            double pw_elig_hahm = births_by_ha[ha-ss.i_h_ag_fert] * nhp.frr_cd4[t][ha-ss.i_h_ag_fert][hm] * hivpop[t][g][ha][hm] / frr_pop_ha;
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

	const std::vector<double> CD4_LOW_LIM({500, 350, 250, 200, 100, 50, 0});
	const std::vector<double> CD4_UPP_LIM({1000, 500, 350, 250, 200, 100, 50});

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
            if(artinit_hahm > hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm])
              artinit_hahm = hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm];
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

            if(artinit_hahm > hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm])
              artinit_hahm = hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm];

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
            if(artinit_hahm > hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm])
              artinit_hahm = hivpop[t][g][ha][hm] + ss.dt * grad[g][ha][hm];
            grad[g][ha][hm] -= artinit_hahm / ss.dt;
            gradART[g][ha][hm][nhp.art0mos] += artinit_hahm / ss.dt;
          }
      }
    }
  }
}

void Model::do_hts_step(int t) {
  
      for(int g = 0; g < ss.ng; g++)
        for(int ha = 0; ha < ss.h_ag; ha++)
          for(int hm = 0; hm < ss.h_ds; hm++)
            hivpop[t][g][ha][hm] += ss.dt * grad[g][ha][hm];

      if(t >= ss.t_art_start)	            
        for(int g = 0; g < ss.ng; g++)
          for(int ha = 0; ha < ss.h_ag; ha++)
            for(int hm = everARTelig_idx; hm < ss.h_ds; hm++)
              for(int hu = 0; hu < ss.h_ts; hu++)
                artpop[t][g][ha][hm][hu] += ss.dt * gradART[g][ha][hm][hu];
}


void Model::remove_hiv_deaths(int t) {

  // remove hivdeaths from pop
  for(int g = 0; g < ss.ng; g++){

    // sum HIV+ population size in each hivpop age group
    std::vector<double> hivpop_ha(ss.h_ag);
    int a = 0;
    for(int ha = 0; ha < ss.h_ag; ha++){
      hivpop_ha[ha] = 0.0;
      for(int i = 0; i < ss.h_ag_span[ha]; i++){
        hivpop_ha[ha] += pop[t][ss.i_hivp][g][a];
        a++;
      }
    }

    // remove hivdeaths proportionally to age-distribution within each age group
    a = 0;
    for(int ha = 0; ha < ss.h_ag; ha++){
      if(hivpop_ha[ha] > 0){
        double hivqx_ha = hivdeaths_ha[g][ha] / hivpop_ha[ha];
        for(int i = 0; i < ss.h_ag_span[ha]; i++){
          hivdeaths[t][g][a] += pop[t][ss.i_hivp][g][a] * hivqx_ha;
          pop[t][ss.i_hivp][g][a] *= (1.0-hivqx_ha);
          a++;
        }
      } else {
        a += ss.h_ag_span[ha];
      }  // end if(pop_ha[ha] > 0)
    }
  }
}
