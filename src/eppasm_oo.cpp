#include "Classes.hpp"

extern "C" SEXP eppasmOOpp(SEXP fp) {
  StateSpace s(fp); // read only state-space
  Parameters p(fp); // read only parameters
  outputSEXP O(s); // create all SEXP for output based on state-space
  Model model(O, s, p); // pass O for model to map to output's address
  model.initiate();
  for (int i = 1; i < p.SIM_YEARS; ++i)
    model.run(i);
  O.finalize(s);
  return O.pop;
}