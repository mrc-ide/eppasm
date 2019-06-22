#include "Classes.hpp"

extern "C" SEXP eppasmOOpp(SEXP fp) {
  oSEXP O(fp);
  Model model(O, fp);
  model.initiate();
  for (int i = 1; i < model.p.SIM_YEARS; ++i)
    model.run(i);
  O.finalize();
  return O.pop;
}