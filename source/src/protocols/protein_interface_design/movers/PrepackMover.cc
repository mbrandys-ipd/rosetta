// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/PrepackMover.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/PrepackMover.hh>
#include <protocols/protein_interface_design/movers/PrepackMoverCreator.hh>
#include <protocols/rosetta_scripts/util.hh>

// Project headers
#include <core/scoring/ScoreFunction.hh>
#include <utility>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <protocols/simple_task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <protocols/minimization_packing/RotamerTrialsMinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>

#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>

#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <core/conformation/Conformation.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <basic/options/option.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.PrepackMover" );




PrepackMover::PrepackMover() :
	protocols::minimization_packing::PackRotamersMover( PrepackMover::mover_name() ),
	scorefxn_( /* NULL */ ),
	jump_num_( 0 ),
	min_bb_( false ),
	mm_( /* NULL */ )
{}

PrepackMover::PrepackMover(
	core::scoring::ScoreFunctionCOP scorefxn,
	core::Size jump_num
) :
	protocols::minimization_packing::PackRotamersMover( PrepackMover::mover_name() ),
	scorefxn_(std::move( scorefxn )),
	jump_num_( jump_num )
{}


protocols::moves::MoverOP
PrepackMover::clone() const {
	return( utility::pointer::make_shared< PrepackMover >( *this ) );
}

protocols::moves::MoverOP
PrepackMover::fresh_instance() const {
	return utility::pointer::make_shared< PrepackMover >();
}

PrepackMover::~PrepackMover() = default;

/// @details separate bound partners (if any), minimize, do rotamer trials, and re-minimize.
void PrepackMover::apply( pose::Pose & pose )
{
	// make a local packertask, reading resfiles and including current rotamers, excluding disulfides
	TR << "Performing repack..." << std::endl;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	TaskFactoryOP tf;
	if ( task_factory() ) tf = utility::pointer::make_shared< TaskFactory >( *task_factory() );
	else tf = utility::pointer::make_shared< TaskFactory >();
	tf->push_back( utility::pointer::make_shared< operation::InitializeFromCommandline >() );
	tf->push_back( utility::pointer::make_shared< operation::IncludeCurrent >() );
	tf->push_back( utility::pointer::make_shared< operation::RestrictToRepacking >() );
	tf->push_back( utility::pointer::make_shared< operation::NoRepackDisulfides >() );

	// incorporating Ian's UnboundRotamer operation.
	// note that nothing happens if unboundrot option is inactive!
	core::pack::rotamer_set::UnboundRotamersOperationOP unboundrot( new core::pack::rotamer_set::UnboundRotamersOperation() );
	unboundrot->initialize_from_command_line();
	core::pack::task::operation::AppendRotamerSetOP unboundrot_operation( new core::pack::task::operation::AppendRotamerSet( unboundrot ) );
	tf->push_back( unboundrot_operation );
	core::pack::dunbrack::load_unboundrot(pose); // adds scoring bonuses for the "unbound" rotamers, if any

	using namespace protocols::simple_task_operations;
	if ( basic::options::option[ basic::options::OptionKeys::docking::norepack1 ]() ) tf->push_back( utility::pointer::make_shared< DockingNoRepack1 >( jump_num_) );
	if ( basic::options::option[ basic::options::OptionKeys::docking::norepack2 ]() ) tf->push_back( utility::pointer::make_shared< DockingNoRepack2 >( jump_num_) );

	//in case there is a resfile, information in this resfile overrides the computed task
	if ( basic::options::option[basic::options::OptionKeys::packing::resfile].user() ) {
		tf->push_back( utility::pointer::make_shared< operation::ReadResfile >() );
	}
	PackerTaskOP task = tf->create_task_and_apply_taskoperations( pose );

	TR << "Pre-minimizing structure..." << std::endl;
	core::kinematics::MoveMapOP mm_general;
	if ( min_bb() ) {
		mm_general = mm(pose)->clone();
		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			if ( !pose.residue(i).is_protein() ) {
				mm_general->set_chi( i, false );
				continue;
			}
			// Check for disulfide bonded cysteines
			if ( pose.residue(i).type().is_disulfide_bonded() ) mm_general->set_chi( i, false );
		}

		protocols::minimization_packing::MinMover min_bb_mover( mm_general, scorefxn_, "lbfgs_armijo_nonmonotone", 1e-5, true/*nblist*/, false/*deriv_check*/  );
		min_bb_mover.apply( pose );
	} else {
		mm_general = utility::pointer::make_shared< core::kinematics::MoveMap >();
		mm_general->clear();
	}

	// separate any bound partners
	protocols::rigid::RigidBodyTransMoverOP translate;
	if ( (jump_num_ > 0) && (pose.conformation().num_chains() > 1) ) {
		TR<<"Translating along jump #"<<jump_num_<<std::endl;
		translate = utility::pointer::make_shared< protocols::rigid::RigidBodyTransMover >( pose, jump_num_ ) ;
		translate->step_size( 1000.0 );
		translate->apply( pose );
	}
	mm_general->set_bb( false );
	mm_general->set_jump( false );
	protocols::minimization_packing::MinMover min_mover( mm_general, scorefxn_, "lbfgs_armijo_nonmonotone", 1e-5, true/*nblist*/, false/*deriv_check*/  );
	// pre-minimize sidechains
	min_mover.apply( pose );

	if ( basic::options::option[basic::options::OptionKeys::docking::dock_rtmin].value() ) {
		protocols::minimization_packing::RotamerTrialsMinMover rtmin( scorefxn_, tf );
		rtmin.apply( pose );
	} else {
		protocols::minimization_packing::PackRotamersMover pack( scorefxn_, task );
		pack.apply( pose );
	}

	// post-minimize
	// using packer include_current() will make sure that these minimized rotamers are used.
	TR << "Post-minimizing structure..." << std::endl;
	min_mover.apply( pose );
	TR << "Done!\n";

	// move back together
	if ( (jump_num_ > 0) && (pose.conformation().num_chains() > 1) ) {
		translate->trans_axis().negate();
		translate->apply( pose );
	}

	//final rescore to get everyone on the same page
	(*scorefxn_)(pose);
	TR.flush();
}


void
PrepackMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data )
{
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	jump_num_ = tag->getOption<core::Size>("jump_number", 1 );
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	min_bb( tag->getOption< bool >( "min_bb", false ));
	if ( min_bb() ) {
		mmf_ = protocols::rosetta_scripts::parse_movemap_factory_legacy( tag, data );
	}
	TR << "Prepack mover with scorefxn " << rosetta_scripts::get_score_function_name(tag) << " over jump number " << jump_num_ << "with min_bb "<<min_bb()<<std::endl;
}

void
PrepackMover::min_bb( bool const m ){
	min_bb_ = m;
}

bool
PrepackMover::min_bb() const{
	return min_bb_;
}

core::kinematics::MoveMapOP
PrepackMover::mm(core::pose::Pose const & pose) const{
	if ( !min_bb() ) TR.Warning << "movemap requested but min_bb is set to false. This is probably wrong!"<<std::endl;
	if ( mm_ != nullptr ) {
		return mm_;
	} else if ( mmf_  != nullptr ) {
		return mmf_->create_movemap_from_pose( pose );
	} else {
		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
		mm->clear();
		mm->set_bb( true );
		return mm;
	}
}

void
PrepackMover::mm( core::kinematics::MoveMapOP mm ){
	mm_ = mm;
}

void
PrepackMover::mmf( core::select::movemap::MoveMapFactoryCOP mmf ){
	mmf_ = mmf;
}

std::string PrepackMover::get_name() const {
	return mover_name();
}

std::string PrepackMover::mover_name() {
	return "Prepack";
}

void PrepackMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	rosetta_scripts::attributes_for_parse_score_function( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "jump_number", xsct_non_negative_integer, "Number of the jump interface to prepack", "1" );
	rosetta_scripts::attributes_for_parse_task_operations( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "min_bb", xsct_rosetta_bool, "Minimize the backbone", "0" );

	XMLSchemaSimpleSubelementList ssl;
	rosetta_scripts::append_subelement_for_parse_movemap_factory_legacy( xsd, ssl );


	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW TO DO", attlist, ssl );
}

std::string PrepackMoverCreator::keyname() const {
	return PrepackMover::mover_name();
}

protocols::moves::MoverOP
PrepackMoverCreator::create_mover() const {
	return utility::pointer::make_shared< PrepackMover >();
}

void PrepackMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PrepackMover::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols
