// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/AddMembraneSpanTermZConstraintCreator.hh
/// @brief      add membrant span constraint
/// @author     Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)

#ifndef INCLUDED_protocols_membrane_AddMembraneSpanTermZConstraint_hh
#define INCLUDED_protocols_membrane_AddMembraneSpanTermZConstraint_hh

// Unit Headers
#include <protocols/membrane/AddMembraneSpanTermZConstraint.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers

// Package Headers
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace membrane {

class AddMembraneSpanTermZConstraint : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	AddMembraneSpanTermZConstraint();

	/// @brief Copy Constructor
	AddMembraneSpanTermZConstraint( AddMembraneSpanTermZConstraint const & src );

	/// @brief Assignment Operator
	AddMembraneSpanTermZConstraint & operator = ( AddMembraneSpanTermZConstraint const & src );

	/// @brief Destructor
	~AddMembraneSpanTermZConstraint() override;

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	protocols::moves::MoverOP clone() const override;

	/// @brief Create a Fresh Instance of this Mover
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Pase Rosetta Scripts Options for this Mover
	void parse_my_tag(
		utility::tag::TagCOP,
		basic::datacache::DataMap &
	) override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (AddMembraneSpanTermZConstraint)

	/// @brief Flip the downstream partner in the membrane
	void apply( core::pose::Pose & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private: // methods

	/////////////////////
	/// Setup Methods ///
	/////////////////////
};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_AddMembraneSpanTermZConstraint_hh
