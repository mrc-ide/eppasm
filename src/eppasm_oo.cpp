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
#include "progression.hpp"

extern "C" SEXP eppasmOOpp(SEXP fp) {
  oSEXP O(fp);
  Model model(O, fp);
  model.initiate();
  for (int i = 1; i < model.p.SIM_YEARS; ++i)
    model.run(i);
  O.finalize();
  return O.pop;
}