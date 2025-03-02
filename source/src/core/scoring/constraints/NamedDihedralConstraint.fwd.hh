// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/constraints/NamedDihedralConstraint.fwd.hh
/// @brief Provides a dihedral constraint based on atom names, rather than numbers. Useful for when atom numbers change
/// @author Tom Linsky ( tlinsky at uw dot edu )


#ifndef INCLUDED_core_scoring_constraints_NamedDihedralConstraint_fwd_hh
#define INCLUDED_core_scoring_constraints_NamedDihedralConstraint_fwd_hh

// Utility header
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace constraints {

class NamedDihedralConstraint;

typedef utility::pointer::shared_ptr< NamedDihedralConstraint > NamedDihedralConstraintOP;
typedef utility::pointer::shared_ptr< NamedDihedralConstraint const > NamedDihedralConstraintCOP;

}
}
}

#endif
