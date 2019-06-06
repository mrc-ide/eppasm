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
#include "Classes.hpp"

// Natural deaths
// -----------------------------------------------------------------------------
void epp_death (popC& pop, hivC& hivpop, artC& artpop) ;

// Migration at year i
// -----------------------------------------------------------------------------
void epp_migration (popC& pop, hivC& hivpop, artC& artpop) ;

// EPP populations aging
// -----------------------------------------------------------------------------
void epp_aging (popC& pop, hivC& hivpop, artC& artpop) ;

// Disease model
// -----------------------------------------------------------------------------
void epp_disease_model (popC& pop, hivC& hivpop, artC& artpop) ;