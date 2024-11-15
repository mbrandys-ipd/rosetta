// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/enzymatic_movers/NTerminalAcetyltransferaseMover.cc
/// @brief   Method definitions for NTerminalAcetyltransferaseMover.
/// @author  Labonte  <JWLabonte@jhu.edu>


// Unit headers
#include <protocols/enzymatic_movers/NTerminalAcetyltransferaseMover.hh>
#include <protocols/enzymatic_movers/NTerminalAcetyltransferaseMoverCreator.hh>
#include <protocols/enzymatic_movers/EnzymaticMover.hh>

// Project headers
#include <protocols/moves/mover_schemas.hh>

#include <core/chemical/VariantType.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// Utility header
#include <utility/tag/XMLSchemaGeneration.hh>

// Basic header
#include <basic/Tracer.hh>


// Construct tracers.
static basic::Tracer TR( "protocols.enzymatic_movers.NTerminalAcetyltransferaseMover" );


namespace protocols {
namespace enzymatic_movers {

// Public methods /////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////
// Default constructor
NTerminalAcetyltransferaseMover::NTerminalAcetyltransferaseMover(): EnzymaticMover( "N-terminal_acetyltransferases" )
{
	type( "NTerminalAcetyltransferaseMover" );
}


// Standard Rosetta methods ///////////////////////////////////////////////////
// General methods
void
NTerminalAcetyltransferaseMover::register_options()
{
	EnzymaticMover::register_options();
}


// Mover methods
std::string
NTerminalAcetyltransferaseMover::get_name() const {
	return mover_name();
}

moves::MoverOP
NTerminalAcetyltransferaseMover::clone() const
{
	return utility::pointer::make_shared< NTerminalAcetyltransferaseMover >( *this );
}

moves::MoverOP
NTerminalAcetyltransferaseMover::fresh_instance() const
{
	return utility::pointer::make_shared< NTerminalAcetyltransferaseMover >();
}


void
NTerminalAcetyltransferaseMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	EnzymaticMover::xml_schema_complex_type_generator()->element_name( mover_name() )
		.complex_type_naming_func( & moves::complex_type_name_for_mover )
		.description( "Enzymatic mover to acetylate the N-termini of a peptide-containing pose" )
		.write_complex_type_to_schema( xsd );
}

// Protected methods //////////////////////////////////////////////////////////
void
NTerminalAcetyltransferaseMover::perform_reaction(
	core::pose::Pose & input_pose,
	core::uint const site,
	std::string const & /*cosubstrate*/ )
{
	using namespace core::chemical;
	using namespace core::pose;

	core::uint const seqpos( get_reactive_site_sequence_position( site ) );
	if ( input_pose.residue( seqpos ).has_variant_type( LOWER_TERMINUS_VARIANT ) ) {
		// We need to check for the LOWER_TERMINUS_VARIANT, not the LOWER_TERMINUS property,
		// because the property can also be used for other PTMs of the N-terminus, e.g., METHYLATED_NTERM_VARIANT
		remove_lower_terminus_type_from_pose_residue( input_pose, seqpos );
		add_variant_type_to_pose_residue( input_pose, ACETYLATED_NTERMINUS_VARIANT, seqpos );
	}
}


// Creator methods ////////////////////////////////////////////////////////////
std::string
NTerminalAcetyltransferaseMoverCreator::keyname() const {
	return NTerminalAcetyltransferaseMover::mover_name();
}

// Return an up-casted owning pointer (MoverOP) to the mover.
protocols::moves::MoverOP
NTerminalAcetyltransferaseMoverCreator::create_mover() const {
	return utility::pointer::make_shared< NTerminalAcetyltransferaseMover >();
}

void
NTerminalAcetyltransferaseMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	NTerminalAcetyltransferaseMover::provide_xml_schema( xsd );
}


// Helper methods /////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that NTerminalAcetyltransferaseMover can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, NTerminalAcetyltransferaseMover const & object_to_output )
{
	object_to_output.show( output );
	return output;
}

}  // namespace enzymatic_movers
}  // namespace protocols
