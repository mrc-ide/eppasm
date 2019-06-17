# Natural deaths
# -----------------------------------------------------------------------------
epp_death <- function(pop, hivpop, artpop) {
  pop$deaths()
  if (pop$MODEL!=0) {
    hivpop$deaths(pop$hiv_sx_prob)
    if (pop$year > pop$p$tARTstart)
      artpop$deaths(pop$hiv_sx_prob)
  }
}

# Migration at year i
# -----------------------------------------------------------------------------
epp_migration <- function(pop, hivpop, artpop) {
  pop$migration()
  if (pop$MODEL!=0) {
    hivpop$migration(pop$hiv_mr_prob)
    if (pop$year > pop$p$tARTstart)
      artpop$migration(pop$hiv_mr_prob)
  }
}

# EPP populations aging
# -----------------------------------------------------------------------------
epp_aging <- function(pop, hivpop, artpop) {
  pop$aging()
  pop$add_entrants()
  if (pop$MODEL==2) 
    pop$sexual_debut()
  if (pop$MODEL!=0) {
    hiv.ag.prob <- pop$hiv_aging_prob()
    artYesNo    <- pop$entrant_art()
    hivpop$aging(hiv.ag.prob)
    hivpop$add_entrants(artYesNo)
    if (pop$MODEL==2) 
      hivpop$sexual_debut()
    if (pop$year > pop$p$tARTstart) {
      artpop$aging(hiv.ag.prob)
      artpop$add_entrants(artYesNo)
      if (pop$MODEL == 2) 
        artpop$sexual_debut()
    }
  }
}

# Disease model
# -----------------------------------------------------------------------------
epp_disease_model <- function(pop, hivpop, artpop) {
  for (time_step in seq_len(pop$hiv_steps_per_year)) {
    if (pop$p$eppmod != "directincid") {
      pop$update_rvec(time_step)
      if (pop$MIX)
        infect <- pop$infect_mix(hivpop, artpop, time_step)
      else
        infect <- pop$infect_spec(hivpop, artpop, time_step)
      pop$update_infection(infect)
      hivpop$update_infection(infect)
    }
    cd4_mort <- scale_cd4_mort(hivpop, artpop)
    hivpop$grad_progress(cd4_mort) # cd4 disease progression and mortality
    artpop$count_death()
    pop$remove_hiv_death(cd4_mort, hivpop, artpop) # Remove hivdeaths from pop
    if (pop$year >= pop$p$tARTstart) # ART initiation
      pop$epp_art_init(hivpop, artpop, time_step)
    hivpop$add_grad_to_pop()
  } # end time step
}

# calculate, distribute eligible for ART, update grad, gradART
# -----------------------------------------------------------------------------
.epp_art_init <- function(hivpop, artpop, time_step) {
  artpop$grad_progress()
  artpop$art_dropout(hivpop) # pass hivpop to receive the drop out
  eligible <- hivpop$eligible_for_art()
  art_elig <- sweep(hivpop$get(year), 1, eligible, "*")
  if (p$pw_artelig[year] & p$artcd4elig_idx[year] > 1)
    art_elig <- update_preg(art_elig, hivpop, artpop) ## add pregnant?
  if (MODEL == 2) # add sexual inactive but eligible for treatment
    art_elig <- art_elig + sweep(hivpop$data_db[,,,year], 1, eligible, "*")
  # calculate number to initiate ART and distribute
  art_curr        <- artpop$current_on_art()
  artnum_ii       <- art_initiate(art_curr, art_elig, time_step)
  art15plus.inits <- pmax(artnum_ii - art_curr, 0)
  artinit         <- art_distribute(art_elig, art15plus.inits)
  if (MODEL == 1) 
    artinit <- pmin(artinit, hivpop$get(year) + DT * hivpop$grad)
  if (MODEL == 2) # split the number proportionally for active and idle pop
    artinit <- hivpop$distribute_artinit(artinit, artpop)
  hivpop$grad <- hivpop$grad - artinit / DT
  artpop$grad_init(artinit)
}
setMembers(popEPP, "public", "epp_art_init", .epp_art_init)