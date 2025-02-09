// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/md/CartesianMD.cc
/// @brief   Cartesian Mover based MD-sim; capable of running FEP simulations and calculations with flag -fep_on true
/// @details
/// @author  Hahnbeom Park, Marisa Brandys

#include <protocols/md/CartesianMDCreator.hh>
#include <protocols/md/CartesianMD.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/pose/symmetry/util.hh>

// Constraints
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>


//Optimization
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMultifunc.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/CartesianMinimizerMap.hh>
#include <core/optimization/types.hh>

//Mass setup
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/Element.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>

// Rosetta Scripts
#include <utility/tag/Tag.hh>
#include <protocols/rosetta_scripts/util.hh>

//Constraints
#include <core/scoring/constraints/ConstraintSet.fwd.hh>

//For reading native
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>

#include <utility/vector1.hh>
#include <numeric/random/random.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <fstream>
// #include <iostream>

//Temporary
#include <core/scoring/rms_util.hh>
#ifdef WIN32
#include <time.h>
#else
#include <sys/time.h>
#endif
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
//mbedit:
#include <core/scoring/Energies.hh>
#include <basic/options/keys/md.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

#include <core/scoring/DerivVectorPair.hh> // MANUAL IWYU
#include <protocols/md/Rattle.hh> // AUTO IWYU For Rattle
#include <protocols/md/thermostat.hh> // AUTO IWYU For Thermostat
#include <basic/datacache/DataMap.hh> // AUTO IWYU For DataMap

namespace protocols {
namespace md {

static basic::Tracer TR( "protocols.md.Cartesian" );

using namespace devel::md;
using namespace core::optimization;
using namespace core;
using namespace ObjexxFCL::format;

// creator




// mover
CartesianMD::CartesianMD()
{
	init();
}

CartesianMD::CartesianMD( core::pose::Pose const & pose,
	core::scoring::ScoreFunctionCOP sfxn,
	core::kinematics::MoveMapCOP movemap )
{
	if ( movemap == nullptr ) {
		core::kinematics::MoveMapOP mmloc( new core::kinematics::MoveMap );
		mmloc->set_jump( true ); mmloc->set_bb( true ); mmloc->set_chi( true );
		mmloc->set( core::id::THETA, true ); mmloc->set( core::id::D, true);
		set_movemap( pose, mmloc );
	} else {
		set_movemap( pose, movemap );
	}

	set_scorefxn( sfxn );
	set_scorefxn_obj( sfxn );
	init();
}

CartesianMD::CartesianMD( core::pose::Pose const & pose,
	core::scoring::ScoreFunction const &sfxn )
{
	core::kinematics::MoveMapOP mmloc( new core::kinematics::MoveMap );
	mmloc->set_jump( true ); mmloc->set_bb( true ); mmloc->set_chi( true );
	mmloc->set( core::id::THETA, true ); mmloc->set( core::id::D, true );
	set_movemap( pose, mmloc );

	set_scorefxn( sfxn );
	set_scorefxn_obj( sfxn );
	init();
}

CartesianMD::CartesianMD( core::pose::Pose const & pose,
	core::scoring::ScoreFunction const &sfxn,
	core::kinematics::MoveMap const &movemap )
{
	set_movemap( pose, movemap.clone() );

	set_scorefxn( sfxn );
	set_scorefxn_obj( sfxn );
	init(); // MDBase
}

CartesianMD::~CartesianMD() = default;

protocols::moves::MoverOP
CartesianMD::clone() const {
	return utility::pointer::make_shared< CartesianMD >(*this);
}


void CartesianMD::get_native_info( core::pose::Pose const &pose )
{

	native_given_ = false;
	if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		native_given_ = true;

		std::string nativepdb = basic::options::option[ basic::options::OptionKeys::in::file::native ]();
		core::chemical::ResidueTypeSetCOP rsd_set( pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t ) );
		core::import_pose::pose_from_file( native_, *rsd_set, nativepdb, core::import_pose::PDB_file );

		// Set resmap
		std::map< core::Size, core::Size > resmap;
		for ( core::Size ires = 1; ires <= pose.size(); ++ires ) {
			if ( !pose.residue( ires ).is_protein() ) continue;
			core::Size ii_pdb( pose.pdb_info()->number( ires ) );

			for ( core::Size jres = 1; jres <= native_.size(); ++jres ) {
				if ( !native_.residue( jres ).is_protein() ) continue;
				core::Size jj_pdb( native_.pdb_info()->number( jres ) );
				if ( ii_pdb == jj_pdb ) {
					resmap[ires] = jres;
					break;
				}
			}
		}
		native_resmap_ = resmap;
	}
}

void CartesianMD::do_initialize( core::pose::Pose &pose )
{

	// Steal utilities in CartesianMinimizer.cc
	CartesianMinimizerMap min_map;

	min_map.setup( pose, *movemap() );

	set_n_dof( min_map.ndofs() );
	set_n_dof_temp( n_dof() - 6 );
	// TR << "mb debug: ndof/ndof_temp:" << n_dof() << "/" << n_dof_temp() << std::endl;
	set_cummulative_time( 0.0 );
	set_pose0( pose );

	// Check initial time
#ifndef WIN32
	gettimeofday(&inittime_, nullptr );
#endif

	// Reallocate
	resize_natm_variables();

	// for adaptive rsr
	Multivec xyz_loc( n_dof() );
	min_map.copy_dofs_from_pose( pose, xyz_loc );
	set_ref_xyz( xyz_loc );
	set_prv_eqxyz( xyz_loc );

	// Mass setup
	core::chemical::ResidueTypeSetCOP rsdtype_set( pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t ) );
	core::chemical::ElementSetCOP element_set( rsdtype_set->element_set() );

	for ( core::Size iatm = 1; iatm <= (core::Size)(min_map.natoms()); ++iatm ) {
		core::id::AtomID AtomID = min_map.get_atom( iatm );
		core::Size resno = AtomID.rsd();
		core::Size atmno = AtomID.atomno();

		set_mass( iatm, pose.residue_type(resno).element_type(atmno)->weight() );
	}

}

// void CartesianMD::use_rattle( bool const value )
// {
// 	//mb TR:
// 	// TR << "use_rattle function is being called." << std::endl;
// 	use_rattle_ = value;
// 	if ( use_rattle_ ) {
// 		set_dt( 0.002 );
// 		TR << "Set Rattle on, changing dt as " << dt() << std::endl;
// 	} else {
// 		set_dt( 0.001 );
// 		TR << "Set Rattle off, changing dt as " << dt() << std::endl;
// 	}
// }

void
CartesianMD::update_restraint( core::pose::Pose & pose,
	CartesianMinimizerMap const &min_map )
{

	if ( Gamma() == 0.0 || !uniform_coord_constrained() ) return;
	std::cout << "mbedit: inside update_restraint!" << std::endl;
	TR.Debug << "Update restraints with Kappa/Gamma = " << Kappa() << "/" << Gamma() << std::endl;
	Multivec curr_eqxyz = get_current_eqxyz();
	Multivec prv_eqxyz_loc( prv_eqxyz() );
	cst_on_pose_dynamic( pose, ref_xyz(), curr_eqxyz, prv_eqxyz_loc, min_map );
	set_prv_eqxyz( prv_eqxyz_loc );

	// clear temporary trj
	renew_trj_scratch();
}

Multivec
CartesianMD::get_current_eqxyz() const
{
	Multivec curr_eqxyz;
	curr_eqxyz.resize( n_dof() );
	utility::vector1< Multivec > const trj_tmp = trj_scratch();
	core::Size const ntrj( trj_tmp.size() );

	for ( core::Size i_dof = 1; i_dof <= n_dof(); ++i_dof ) {
		curr_eqxyz[i_dof] = 0.0;
		for ( core::Size i_trj = 1; i_trj <= ntrj; ++i_trj ) {
			curr_eqxyz[i_dof] += trj_tmp[i_trj][i_dof];
		}

		if ( ntrj > 0 ) curr_eqxyz[i_dof] /= (core::Real)(ntrj);
	}
	return curr_eqxyz;
}

void
CartesianMD::cst_on_pose_simple ( 
	core::pose::Pose &pose ) const
{
	using namespace core::scoring::constraints;
	using namespace core;

	//mb TR line:
	// TR << "cst on pose fxn called; if next TR is Reporting initial... then this fxn did nothing." << std::endl;
	//mb note 12/18/23: this fxn did not get used

	// First, add cst_file info into pose
	if ( basic::options::option[ basic::options::OptionKeys::constraints::cst_fa_file ].user() ) {
		TR << "mb debug Set constraints from input file..." << std::endl;
		pose.remove_constraints();
		scoring::constraints::add_fa_constraints_from_cmdline_to_pose( pose );
	}

	// Next, set coordinate constraint if specified
	if ( uniform_coord_constrained() ) {
		TR << "mb debug cst_on_pose_simple called, with std dev: " << cst_sdev() << std::endl;
		// pose.remove_constraints();

		for ( core::Size i_res = 1; i_res <= pose.size(); ++i_res ) {
			std::string resname = pose.residue(i_res).name();

			if ( pose.residue(i_res).is_ligand() ) continue; //mb added, skip ligand res

			core::Size iatm;
			if ( pose.residue(i_res).has(" CA " ) ) {
				iatm = pose.residue( i_res ).atom_index(" CA ");
				id::AtomID atomID( iatm, i_res );
				core::Vector xyz = pose.xyz( atomID );
				scoring::func::FuncOP fx( new scoring::func::HarmonicFunc( 0.0, cst_sdev() ) );
				// std::out os; //mb added
				// fx.show_definition(os); // mb added
				pose.add_constraint(  ConstraintCOP( ConstraintOP(
					new CoordinateConstraint( atomID, atomID, xyz, fx )
					)));
				// TR << "mb debug alpha carbon constrained." << std::endl; //comment out when running, this is just for testing purposes
			
			} 
			
			if ( cst_beta_ && (pose.residue(i_res).has(" CB " )) ) {
				iatm = pose.residue( i_res ).atom_index(" CB ");
				id::AtomID atomID_cb( iatm, i_res );
				core::Vector xyz_cb = pose.xyz( atomID_cb );
				scoring::func::FuncOP fx_cb( new scoring::func::HarmonicFunc( 0.0, cst_sdev() ) );
				// scoring::func::HarmonicFunc::show_definition(*fx); // mb added
				pose.add_constraint(  ConstraintCOP( ConstraintOP(
					new CoordinateConstraint( atomID_cb, atomID_cb, xyz_cb, fx_cb )
					)));
				// TR << "mb debug beta carbon constrained." << std::endl; //comment out when running, this is just for testing purposes
			}

			
		}
	}
}

void
CartesianMD::cst_on_pose_dynamic( core::pose::Pose &pose,
	Multivec const &ref_xyz,
	Multivec const &curr_eqxyz,
	Multivec &prv_eqxyz,
	CartesianMinimizerMap const &min_map ) const
{
	using namespace core::scoring::constraints;

	// not supporting cst_fa_file yet
	if ( basic::options::option[ basic::options::OptionKeys::constraints::cst_fa_file ].user() ) return;

	// Remove all the constraints first
	//mb TR line:
	// TR << "inside cst on pose dynamic, constraints are about to be removed!" << std::endl;
	//mb note 12-18 23-this was not called.

	pose.remove_constraints();

	//CartesianMinimizerMap min_map;
	//min_map.setup( pose, *movemap() );

	// Next, set coordinate constraint if specified
	Multivec prv_eqxyz0( prv_eqxyz );

	std::ofstream rsrfile( rsrfilename().c_str(), std::ios_base::app );
	core::Size modality = (core::Size)(cummulative_time()*100+1)%100;
	bool write_rsr = write_dynamic_rsr() && ( modality <= 2);

	TR.Debug << cummulative_time() << " " << ((core::Size)(cummulative_time()*100))%100 << " " << write_rsr << std::endl;
	TR.Debug << "ref/curr/prv_eqxyz0? " << ref_xyz.size() << " " << curr_eqxyz.size() << " " << prv_eqxyz0.size() << std::endl;

	if ( uniform_coord_constrained() ) {
		if ( write_rsr ) {
			rsrfile << "dynamic rsr: " << cummulative_time() << std::endl;
			rsrfile << "Res | X Y Z | prv_eq X Y Z | delta X Y Z" << std::endl;
		}
		core::Real rmsd( 0.0 );
		core::Real rmsd_rsr2ref( 0.0 );
		core::Real rmsd_rsr2crd( 0.0 );
		core::Size nrsr_dof( 0 );
		//for ( core::Size i_res = 1; i_res <= pose.size(); ++i_res ) {
		for ( core::Size i_atm = 1; i_atm <= n_dof()/3; ++i_atm ) {
			core::id::AtomID atomID = min_map.get_atom( i_atm );
			core::Size resno = atomID.rsd();
			core::Size atmno = atomID.atomno();

			std::string resname = pose.residue(resno).name();
			std::string atmname = pose.residue(resno).atom_name(atmno);

			bool is_protein_ca = ( atmname == " CA " );
			bool is_water = ( resname == "TP3" && atmname == " O  " );
			if ( !(is_protein_ca || is_water ) ) continue;

			nrsr_dof++;
			TR.Debug << "on " << resno << " " << atmname << std::endl;

			core::Vector xyzmix( 0.0 );
			core::Real dist_res( 0.0 );
			core::Real dist_togo( 0.0 );
			core::Real dist_rsr( 0.0 );

			for ( core::Size i = 1; i <= 3; ++i ) {
				core::Size i_dof = (i_atm-1)*3 + i;
				xyzmix[i-1] = (1.0 - Kappa())*prv_eqxyz[ i_dof ];
				xyzmix[i-1] += Kappa()*(Gamma()*curr_eqxyz[ i_dof ] +
					(1.0-Gamma())*ref_xyz[ i_dof ] );

				// update prv_eq crd
				prv_eqxyz[ i_dof ] = xyzmix[i-1];
				dist_res  += (curr_eqxyz[ i_dof ] - ref_xyz[ i_dof ])  *(curr_eqxyz[ i_dof ] - ref_xyz[ i_dof ]);
				dist_togo += (curr_eqxyz[ i_dof ] - prv_eqxyz[ i_dof ])*(curr_eqxyz[ i_dof ] - prv_eqxyz[ i_dof ]);
				dist_rsr  += (prv_eqxyz[ i_dof ]  - ref_xyz[ i_dof ])  *(prv_eqxyz[ i_dof ] - ref_xyz[ i_dof ]);
			}
			rmsd += dist_res;
			rmsd_rsr2crd += dist_togo;
			rmsd_rsr2ref += dist_rsr;
			dist_res = std::sqrt(dist_res);

			core::Real sdev( cst_sdev() );

			core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0,  sdev ) );
			pose.add_constraint( ConstraintCOP( ConstraintOP(
				new CoordinateConstraint( atomID, atomID, xyzmix, fx )
				)));

			if ( write_rsr ) {
				rsrfile << I(4,resno);
				rsrfile << " | " << F(8,3,xyzmix[0])  << " " << F(8,3,xyzmix[1]) << " " << F(8,3,xyzmix[2]);
				rsrfile << " | " << F(8,3,prv_eqxyz0[3*(i_atm-1)+1]);
				rsrfile << " "   << F(8,3,prv_eqxyz0[3*(i_atm-1)+2]);
				rsrfile << " "   << F(8,3,prv_eqxyz0[3*(i_atm-1)+3]);
				rsrfile << " | " << F(8,3,dist_res);
				rsrfile << " "   << F(8,3,prv_eqxyz[3*(i_atm-1)+1] - ref_xyz[3*(i_atm-1)+1]);
				rsrfile << " "   << F(8,3,prv_eqxyz[3*(i_atm-1)+2] - ref_xyz[3*(i_atm-1)+2]);
				rsrfile << " "   << F(8,3,prv_eqxyz[3*(i_atm-1)+3] - ref_xyz[3*(i_atm-1)+3]);
				rsrfile << std::endl;
			}

		} //iatm

		if ( write_rsr ) {
			rmsd         /= (core::Real)(nrsr_dof);        rmsd = std::sqrt(rmsd);
			rmsd_rsr2ref /= (core::Real)(nrsr_dof); rmsd_rsr2ref = std::sqrt(rmsd_rsr2ref);
			rmsd_rsr2crd /= (core::Real)(nrsr_dof); rmsd_rsr2crd = std::sqrt(rmsd_rsr2crd);

			rsrfile << "At time " << cummulative_time() << ", RMSD of rsr to ref: " << F(8,3,rmsd_rsr2ref);
			rsrfile << " , crd to rsr: " << F(8,3,rmsd_rsr2crd);
			rsrfile << " , crd to ref:" << F(8,3,rmsd) << std::endl;
		}
	} //if uniform_coordinate_constrained

	rsrfile.close();

}

void CartesianMD::apply( core::pose::Pose & pose ) {
	using namespace core::optimization;
	// TR << "cartmd apply fxn start." << std::endl;
	get_native_info( pose );

	if ( movemap_factory_ ) {
		// reset the movemap if we have a valid MoveMapFactory
		set_movemap( pose, movemap_factory_->create_movemap_from_pose( pose ) );
	}

	//fpd we have to do this here since this the first time "seeing" the symm pose
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		core::pose::symmetry::make_symmetric_movemap( pose, *movemap() );
	}

	// setup the map of the degrees of freedom
	const std::string minopt( "lbfgs_armijo_nonmonotone" );
	MinimizerOptions options_init( minopt, 0.0001, true, false, false );
	options_init.max_iter( ncyc_premin() );
	MinimizerOptions options_final( minopt, 0.000001, true, false, false );
	options_final.max_iter( ncyc_postmin() );

	// Get initial crd info
	const core::pose::Pose pose0( pose );


	// Set constraint
	// TR << "Set constraint on initial." << std::endl;
	if ( (*scorefxn())[ core::scoring::coordinate_constraint ] != 0.0 ) {
		set_constraint( 1.0 ); //sets sdev to 1.0, and uniform_coord_constraint to true
		TR << "mb debug: coord cst wt is non 0, cst on pose simple fxn should be called next. " << std::endl;
		cst_on_pose_simple( pose );
	}
	

	//check if boolean calc_intE_ is true, and if so, confirm a ligand is present
	//if no ligand when calc_intE_ is true, change to false and print message
	if (calc_intE_) {
		bool has_ligand (false);
		//loop thru resis; if any resi is ligand, has_ligand becomes true
		for ( core::Size i_res = 1; i_res <= pose.size(); ++i_res ) {
			if ( pose.residue(i_res).is_ligand() ) {
				has_ligand = true;
			} //mb added, skip ligand res
		}
		if ( !has_ligand ) { //if has_ligand true, then !true becomes false and calc_intE_ stays true
			TR << "[WARNING] Even though calc_intE was set to True, pose does not have ligands. Switching calc_intE to false." << std::endl;
			calc_intE_ = false;
		}
	}
	
	TR << "Reporting initial..." << std::endl;
	
	CartesianMinimizerMap min_map;
	min_map.setup( pose, *movemap() );
	// TR << "apply fxn 1st report md call start:" << std::endl;
	core::Size istep = 0; //so that reportMD for initial stuff has nstep 0
	report_MD( pose, min_map, true, istep);
	header_on_ = false;
	// TR << "apply fxn 1st report md call end." << std::endl;

	//if premin steps are not 0, do premin:
	if (ncyc_premin() != 0) {
		do_minimize( pose, options_init, true );
		//dump the pose after premin:
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		std::stringstream ss_mar1;
		std::string mar_name1 = option[out::prefix];
		ss_mar1 << mar_name1 << "_preMD_minimized.pdb";
		pose.dump_pdb( ss_mar1.str() );


		TR << "Reporting after pre-MD minimization..." << std::endl;
		report_MD( pose, min_map, true, istep);
	}
		

	
	// TR << "apply fxn 2nd report md call start:" << std::endl;
	
	// TR << "apply fxn 2nd report md call end." << std::endl;
	// TR << "scorecomp of total complex:" << std::endl;
	if ( report_scorecomp() ) scorefxn()->show( TR, pose );

	// Set initial minobj_ after minimization
	// These are used to pick the best pose along trajectory (if selectmode == minobj )
	set_pose_minobj( pose );
	set_Emin_obj( 1.0e6 );
	set_time_minobj( 0.0 );


	// TR << "dT is set to:" << dt() << std::endl; //mbedit TR statement to check dT when needed.
	
	// Main
	if ( !scheduled() ) { // typical run - mbedit note: this is my run
		TR << "dt() is: " << dt() << std::endl;
		do_MD( pose, nstep(), temp0(), true );
	} else { // scheduled run

		for ( core::Size i_step = 1; i_step <= mdsch().size(); ++i_step ) {
			std::string const runtype( mdsch(i_step).type );
			if ( runtype == "sch" ) {
				TR << "Changing schedule, Nstep/Temp:";
				TR << mdsch(i_step).nstep << " " << mdsch(i_step).temp0 << std::endl;
				if ( i_step == 1 ) {
					do_MD( pose, mdsch(i_step).nstep, mdsch(i_step).temp0, true );
				} else {
					do_MD( pose, mdsch(i_step).nstep, mdsch(i_step).temp0, false );
				}
			} else if ( runtype == "repack" ) {
				//
			}

		}
	}

	/// Selection for returning structure
	if ( selectmode() == "final" ) {
		// just return final pose
		TR << "Returning final structure from MD sim..." << std::endl;
	} else if ( selectmode() == "minobj" ) {
		pose = pose_minobj();
		TR << "Returning minimum objective function structure at ";
		TR << time_minobj() << " in MD trajectory..." << std::endl;
	}

	if (ncyc_postmin() != 0) {
		TR << "Running minimize on final MD structure..." << std::endl;
		do_minimize( pose, options_final, true );
	}
	
	if ( report_scorecomp() ) scorefxn()->show( TR, pose );
	TR << "MD Done! " << std::endl;
}


void CartesianMD::do_minimize( core::pose::Pose &pose,
	core::optimization::MinimizerOptions const &options,
	bool const &show_energy )
{
	CartesianMinimizer minimizer;

	if ( show_energy ) {
		core::Real score_before = scorefxn()->score( pose );

		minimizer.run( pose, *movemap(), *scorefxn(), options );

		core::Real score_after = scorefxn()->score( pose );
		TR << "Energy before/after Min: " << score_before << " " << score_after << std::endl;

	} else {
		minimizer.run( pose, *movemap(), *scorefxn(), options );
	}
}

Real
CartesianMD::calculate_free_receptor_score(
	core::pose::Pose pose //this needs to be the pose with protein-ligA-ligB
)
{
	// delete ligand (modifying fxn from GALigandDock.cc version)
	core::Size lastres = pose.total_residue(); //because my ligand is the last residue
	pose.delete_residue_range_slow( (lastres - 1), lastres );

	//mbedit - using this to confirm i deleted ligand
	// std::stringstream ss_mar2;
	// ss_mar2 << "noLigand.pdb";
	// pose.dump_pdb( ss_mar2.str() );
	//mbedit end

	// get sum of just lambda-scaled e-terms
	core::Real total_E ( scorefxn()->score( pose ) );

	// std::cout << "protein score breakdown:" << std::endl; // debugging
	// scorefxn()->show( pose ); //for debugging
	
	return total_E; 

}

Real
CartesianMD::calculate_free_ligandA_score(
	core::pose::Pose pose_ref //this needs to be the pose with protein-ligA-ligB
)
{
	// just get ligand (modifying fxn from GALigandDock.cc version)
	core::Size ligA_resnum = ( pose_ref.total_residue() - 1 );

	core::pose::PoseOP pose( new core::pose::Pose );
    pose->append_residue_by_jump( pose_ref.residue(ligA_resnum), 0);
	core::pose::initialize_disulfide_bonds( *pose );

	// //mbedit confirming i have just ligand
	// std::stringstream ss_mar2;
	// ss_mar2 << "justLigandA.pdb";
	// pose->dump_pdb( ss_mar2.str() );
	// //mbedit end

	// get sum of just lambda-scaled e-terms
	core::Real total_E ( scorefxn()->score( *pose ) );
	
	return total_E;

}

Real
CartesianMD::calculate_free_ligandB_score(
	core::pose::Pose pose_ref //this needs to be the pose with protein-ligA-ligB
)
{
	// just get ligand (modifying fxn from GALigandDock.cc version)
	core::Size ligB_resnum = pose_ref.total_residue();

	core::pose::PoseOP pose( new core::pose::Pose );
    pose->append_residue_by_jump( pose_ref.residue(ligB_resnum), 0);
	core::pose::initialize_disulfide_bonds( *pose );

	// //mbedit confirming i have just ligand
	// std::stringstream ss_mar2;
	// ss_mar2 << "justLigandB.pdb";
	// pose->dump_pdb( ss_mar2.str() );
	// //mbedit end

	// get sum of just lambda-scaled e-terms
	core::Real total_E ( scorefxn()->score( *pose ) );
	
	return total_E;

}

utility::vector1<double>
CartesianMD::calculate_complexligA_score(
	core::pose::Pose pose //this needs to be the pose with protein-ligA-ligB
)
{
	
	core::Size ligB_resnum = pose.total_residue(); //removes ligB
	pose.delete_residue_slow( ligB_resnum ); 

	//get just the protein - ligA interaction energy
	core::Real complex_energy ( scorefxn()->score( pose ) );
	core::Real apc = 0.0;
	if ( (*scorefxn())[ core::scoring::atom_pair_constraint ] != 0.0 ) {
		apc = pose.energies().total_energies_weighted()[ core::scoring::atom_pair_constraint ];
	}
	// std::cout << "within complex ligA calc fxn, complex_energy without apc: " << complex_energy << ", apc: " << apc << std::endl;
	complex_energy = complex_energy - apc;
	utility::vector1<double> energies = {complex_energy,apc};

	return energies; 

}

utility::vector1<double>
CartesianMD::calculate_complexligB_score(
	core::pose::Pose pose //this needs to be the pose with protein-ligA-ligB
)
{
	
	core::Size ligA_resnum = (pose.total_residue() - 1); //removes ligB
	pose.delete_residue_slow( ligA_resnum ); 

	// //mbedit - using this to confirm i deleted ligandA only
	// std::stringstream ss_mar2;
	// ss_mar2 << "protein-ligB_complex.pdb";
	// pose.dump_pdb( ss_mar2.str() );
	// //mbedit end

	//get just the protein - ligA interaction energy
	core::Real complex_energy ( scorefxn()->score( pose ) );
	core::Real apc = 0.0;
	if ( (*scorefxn())[ core::scoring::atom_pair_constraint ] != 0.0 ) {
		apc = pose.energies().total_energies_weighted()[ core::scoring::atom_pair_constraint ];
	}	
	// std::cout << "within complex ligB calc fxn, complex_energy without apc: " << complex_energy << ", apc: " << apc << std::endl;
	complex_energy = complex_energy - apc;
	utility::vector1<double> energies = {complex_energy,apc};

	return energies; 

}

void CartesianMD::do_MD( core::pose::Pose & pose,
	core::Size const &nstep,
	core::Real const &temp0,
	bool const &initialize )
{
	if ( initialize ) {
		TR << "do_MD start: Running MD simulations for " << nstep*dt() << " ps at " << temp0 << " K... " << std::endl;
		do_initialize( pose );
	}

	// // Reporting about adaptive rsr
	// if ( uniform_coord_constrained() ) {
	// 	TR << "Run uniform restrained simulation ";
	// 	if ( Gamma() == 0.0 ) {
	// 		TR << " with static restraints on starting pose." << std::endl;
	// 	} else {
	// 		TR << " with dynamic restraints with Kappa/Gamma = ";
	// 		TR << Kappa() << " " << Gamma() << std::endl;
	// 	}
	// }

	// Set dof variables
	CartesianMinimizerMap min_map;
	Multivec xyz_loc( n_dof() );
	min_map.setup( pose, *movemap() );
	min_map.copy_dofs_from_pose( pose, xyz_loc );
	set_xyz( xyz_loc );

	// Setup RATTLE using min_map
	md::Rattle rattle( pose, min_map );
	
	if ( use_rattle_ ) set_n_dof_temp( n_dof() - 6 - rattle.ncst() );

	// This should come after Rattle setup to get n_dof_temp_
	core::pose::PoseOP pose_temp ( new core::pose::Pose ( pose )); //create pointer to pass into init_vel mb edit
	if ( initialize ) initialize_velocity( temp0, pose_temp );

	// Set thermostat
	Thermostat thermostat( temp0, n_dof_temp() );

	// Start MD integrator
	scorefxn()->setup_for_minimizing( pose, min_map );
	for ( core::Size istep = 1; istep <= nstep; istep++ ) {
		set_cummulative_time( cummulative_time() + dt() );
		// Dump pdb
		if ( istep%dumpstep_ == 0 ) {

			using namespace basic::options;
			using namespace basic::options::OptionKeys;

			float marTime = static_cast<float>(cummulative_time());
			std::stringstream ss_mar;
			std::string mar_name = option[out::prefix];
			ss_mar << mar_name << "_" << marTime << "ps.pdb";
			pose.dump_pdb( ss_mar.str() );
		}
		
		// Report
		if ( istep%report_step() == 0 ) {
			report_MD( pose, min_map, false, istep); // if 3rd item is true, then including trajectory
		}

		bool update_score( false );
		// if ( istep%context_update_step() == 0 ) update_score = true;

		// // For adaptive restraint
		// if ( istep%md_rsr_update_stepsize() == 0 ) {
		// 	// std::cout << "mbedit: Adaptive Restraint is updating." << std::endl;
		// 	update_restraint( pose, min_map );
		// } else if ( trj_scratch().size() < 100 ) {
		// 	// make sure scratch space doesn't use too much memory
		// 	// std::cout << "mbedit: trj going to scratch?." << std::endl;
		// 	add_trj_scratch( xyz() );
		// }

		// Integrate
		VelocityVerlet_Integrator( pose, min_map, rattle, update_score, istep); //mbedit, added nstep as a final passed in variable for my reporting

		/*
		// let's get avrg velocities
		core::Real v2sum( 0.0 ), temperature( 0.0 );
		for (core::Size i_dof = 1; i_dof<=vel_.size(); ++i_dof) {
		int i_atm = (i_dof+2)/3;
		// pass Virtual atoms
		if ( mass_[i_atm] < 1e-3 ) continue;
		TR << "vel: " << i_atm << " " << vel_[i_dof] << std::endl;
		temperature += mass_[i_atm]*vel_[i_dof]*vel_[i_dof];
		v2sum += vel_[i_dof]*vel_[i_dof];
		}
		temperature /= n_dof_temp_*Boltzmann;
		*/

		// Calculate/re-eval temperature
		set_temperature( thermostat.get_temperature( vel(), mass() ) );

		if ( istep%thermostat.nstep_per_update() == 0 ) {
			// TR << "mb debug, thermostat is rescaling at end of step: " << istep << std::endl;
			Multivec vel_loc( vel() );
			// TR << "mb debug, vel_loc[1] before rescale: " << vel_loc[1] << std::endl;
			thermostat.rescale( vel_loc, dt(), mass(), nve_mode_ );
			set_vel( vel_loc );
			// TR << ", vel_loc[1] after rescale: " << vel_loc[1] << std::endl;
			set_temperature( thermostat.get_temperature( vel(), mass() ) );
		}
		//TR << "v2avrg/Temp/Temp2: " << std::sqrt(v2sum/vel_.size()) << " " << temperature << " " << temperature_ << std::endl;
		set_kinetic_energy( 0.5*temperature()*n_dof()*GasConst );

	}
	// TR << "do MD: minmap copy dofs starts (do_MD ends statement is just after.)" << std::endl;
	min_map.copy_dofs_to_pose( pose, xyz() );
	// TR << "do_MD ends." << std::endl;
}

void CartesianMD::VelocityVerlet_Integrator( core::pose::Pose &pose,
	CartesianMinimizerMap &min_map,
	md::Rattle &rattle,
	bool const update_score,
	core::Size istep ) //mbedit added istep
{
	core::Real dt2_2 = dt()*dt()*0.5;
	// TR << "start of VVI fxn.//////////" << std::endl;
	// Use previous acceleration here
	// and integrate first half of the velocity
	Multivec xyz_loc( xyz() ), vel_loc( vel() ), acc_loc( acc() );
	
	for ( core::Size i_dof = 1; i_dof <= n_dof(); ++i_dof ) {

		xyz_loc[i_dof] += vel_loc[i_dof]*dt() + acc_loc[i_dof]*dt2_2;
		
		vel_loc[i_dof] += 0.5*acc_loc[i_dof]*dt();

	}

	//rattle is off by default
	if ( use_rattle_ ) {
		rattle.run_rattle1( dt(), xyz_loc, vel_loc, mass() );
	}

	// Reflect change in coordinates into pose
	min_map.copy_dofs_to_pose( pose, xyz_loc );
	
	// Don't need this unless context needs to be updated
	if ( update_score ) scorefxn()->score( pose );
	Multivec force;
	// TR << "in VVI, before f_ros is called (and takes in *scorefxn) /////////////////////" << std::endl;
	CartesianMultifunc f_ros( pose, min_map, *scorefxn(), false, false ); //this line yielded none of my fa_elec print statments.
	// TR << "in VVI, after f_ros is called but before .dfunc fxn is called /////////////" << std::endl;
	f_ros.dfunc( xyz_loc, force ); // only saw "elec deriv calc: applying scale_factor" print statement here.

	// //mb debug statement block - not production code start
	// std::cout << "mb debug, ndofs: " << n_dof() << ":force-before-changes:" << force << std::endl;
	// for ( core::Size i_res = 1; i_res <= pose.size(); ++i_res ) {
	// 	for (core::Size i_atm = 1; i_atm <= pose.residue(i_res).natoms(); ++i_atm) {
	// 		// std::cout << "mb debug, ligand? " << pose.residue(i_res).is_ligand() << ", res " << i_res << " atm  " << i_atm << " name:" << pose.residue(i_res).atom_name(i_atm) << std::endl;
	// 		id::AtomID atomID( i_atm, i_res );
	// 		core::Vector xyz = pose.xyz( atomID );
	// 		std::cout << "mb debug, ligand? " << pose.residue(i_res).is_ligand() << ",resnum: " << i_res << " atomnum: " << i_atm << " xyz coord: " << xyz[0] << "," << xyz[1] << "," << xyz[2] << std::endl;
	// 	}
	// } //mb debug statement block - not production code end


	//mb edit: if fep on, do gradient matching
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if (option[ basic::options::OptionKeys::md::fep_on ]) {
		
		// mapping / idof force update code, hardcoded for now! mb edit
		// utility::vector1<int> atm_map = {1, 68, 2, 67, 3, 60, 4, 59, 5, 57, 6, 56, 7, 63, 8, 62, 9, 51, 10, 53, 11, 49, 12, 61, 13, 65, 14, 64, 15, 58, 16, 66, 17, 50, 18, 54, 19, 55, 20, 52, 21, 47, 22, 48, 23, 46, 24, 44, 26, 84, 27, 85, 28, 86, 29, 81, 30, 82, 31, 83, 32, 78, 33, 77, 34, 76, 35, 75, 36, 80, 37, 79, 38, 72, 39, 71, 40, 74, 41, 73, 42, 70, 43, 69};
		utility::vector1<int> atm_map = option[ basic::options::OptionKeys::md::md_atm_map ];
		// std::vector<int> atm_map = {1, 68, 2, 67};
		// std::cout << "mb debug, atm_map: " << atm_map << ", size of atm_map: " << atm_map.size() << std::endl;

		//calculate shift so that my hard-coded atom map makes sense regardless of how many non-ligand atoms are present.
		int total_atoms = n_dof()/3;
		int total_lig_atoms = pose.residue( pose.size() - 1 ).natoms() + pose.residue( pose.size() ).natoms(); //adds natom number of ligA and ligB
		int idx_shift = total_atoms - total_lig_atoms;
		// std::cout << "mb debug, total_atoms:" << total_atoms << ":total_lig_atoms:" << total_lig_atoms << ":idx_shift:" << idx_shift << std::endl;

		for ( core::Size i = 1; i <= atm_map.size() - 1; i+=2 ) {
			// std::cout << "i is: " << i << std::endl;
			// std::cout << "BEFORE SHIFT:i_atm A: " << atm_map[i] << ",i_atm B: " << atm_map[i+1] << std::endl;
			int atm_a = atm_map[i] + idx_shift; //stores atom number for ligA atom
			int atm_b = atm_map[i + 1] + idx_shift; //stores atom number for ligB atom
			// std::cout << "AFTER SHIFT:i_atm A: " << atm_a << ",i_atm B: " << atm_b << std::endl;
			float x_a = force[(atm_a * 3) - 2]; //gets gradient for x-coord of atm_a
			float y_a = force[(atm_a * 3) - 1]; //gets gradient for y-coord of atm_a
			float z_a = force[(atm_a * 3)]; //gets gradient for z-coord of atm_a

			float x_b = force[(atm_b * 3) - 2]; //gets gradient for x-coord of atm_b
			float y_b = force[(atm_b * 3) - 1]; //gets gradient for y-coord of atm_b
			float z_b = force[(atm_b * 3)]; //gets gradient for z-coord of atm_b

			force[(atm_a * 3) - 2] = x_a + x_b; //sets force xcoord for ligA atom to sum of x_a and x_b
			force[(atm_a * 3) - 1] = y_a + y_b; //sets force ycoord for ligA atom to sum of y_a and y_b
			force[(atm_a * 3)] = z_a + z_b; //sets force zcoord for ligA atom to sum of z_a and z_b

			force[(atm_b * 3) - 2] = x_a + x_b; //sets force xcoord for ligB atom to sum of x_a and x_b
			force[(atm_b * 3) - 1] = y_a + y_b; //sets force ycoord for ligB atom to sum of y_a and y_b
			force[(atm_b * 3)] = z_a + z_b; //sets force zcoord for ligB atom to sum of z_a and z_b

			// std::cout << "xa/ya/za " << x_a << "/" << y_a << "/" << z_a << " xb/yb/zb " << x_b << "/" << y_b << "/" << z_b << std::endl;
		}
		// std::cout << "mb debug,force-afterFEP-changes:" << force << std::endl;
	}

	// TR << "in VVI, after f_ros.dfunc is called. End of debug checking block in cartmd./////////" << std::endl;
	/*mb note - about f_ros.dfunc line:
	Here, core.optimization.CartesianMultifunc.cc::dfunc() calls core.optimization.cartesian_minimize.cc::cartesian_dfunc()
	which in turn calls core.optimization.cartesian_minimize.cc::cartesian_collect_atompairE_deriv
	and cartesian_collect_torsional_deriv to fill the variable here called force with values
	from those functions' latter two functions. So this is the line that fills empty force with values!
	I think this means there's no point doing a dforce check because every spin of this Verlet fxn reinitializes
	force before filling it up again.
	*/
	
	// //mb block start (variables being initialized for debug_mode_ comments)
	// //store values before force gets applied in next for loop
	// Multivec vel_loc_pre( vel_loc ), acc_loc_pre( acc_loc );
	// double avg_dvel, avg_dacc, avg_vel_pre, avg_vel_post, avg_acc_pre, avg_acc_post;
	// avg_dvel = avg_dacc = avg_vel_pre = avg_vel_post = avg_acc_pre = avg_acc_post = 0.0;
	// double max_dvel, min_dvel;
	// max_dvel = min_dvel = 0.0;
	// double max_vel_pre, min_vel_pre, max_vel_post, min_vel_post;
	// max_vel_pre = min_vel_pre = max_vel_post = min_vel_post = vel_loc_pre[1];
	// double max_dacc, min_dacc;
	// max_dacc = min_dacc = 0.0; 
	// double max_acc_pre, min_acc_pre, max_acc_post, min_acc_post;
	// max_acc_pre = min_acc_pre = max_acc_post = min_acc_post = acc_loc_pre[1];
	// //mb block end

	// Here, convert force into acceleration
	// and integrate remaining half of velocity
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	bool cap_on = option[basic::options::OptionKeys::md::cap_movements]; //default false
	Real vel_cap = option[basic::options::OptionKeys::md::max_vel]; //default 10.0
	Real acc_cap = option[basic::options::OptionKeys::md::max_acc]; //default 3e3

	for ( core::Size i_dof = 1; i_dof <= n_dof(); ++i_dof ) {
		core::Size i_atm = (i_dof+2)/3;
		
		// pass Virtual atoms
		if ( mass(i_atm) < 1e-3 ) continue;

		//acc(i_dof) = -MDForceFactor*force[i_dof]/mass(i_atm);
		//vel(i_dof) += 0.5*acc(i_dof)*dt();

		// if (debug_mode_) {
		// 	// //accumulating to get pre-calculation average vel and acc
		// 	// avg_vel_pre += vel_loc[i_dof];
		// 	// avg_acc_pre += acc_loc[i_dof];
		// 	//get max and min of pre-calculation average vel and acc
		// 	if ( vel_loc[i_dof] > max_vel_pre ) {
		// 		max_vel_pre = vel_loc[i_dof];
		// 	} 
		// 	if ( vel_loc[i_dof] < min_vel_pre ) {
		// 		min_vel_pre = vel_loc[i_dof];
		// 	}
		// 	if ( acc_loc[i_dof] > max_acc_pre ) {
		// 		max_acc_pre = acc_loc[i_dof];
		// 	}
		// 	if ( acc_loc[i_dof] < min_acc_pre ) {
		// 		min_acc_pre = acc_loc[i_dof];
		// 	}

		// }
		// std::cout << "mb debug loop(acc/vel from force loop), iatm: " << i_atm << " idof: " << i_dof << " force[idof]: " << force[i_dof] << std::endl;

		acc_loc[i_dof] = -MDForceFactor*force[i_dof]/mass(i_atm); 
		vel_loc[i_dof] += 0.5*acc_loc[i_dof]*dt();

		if (cap_on) {
			if (vel_loc[i_dof] > vel_cap) {
				vel_loc[i_dof] = vel_cap;
			}

			if (vel_loc[i_dof] < (vel_cap * -1.0)) {
				vel_loc[i_dof] = (vel_cap * -1.0);
			}

			if (acc_loc[i_dof] > acc_cap) {
				acc_loc[i_dof] = acc_cap;
			}

			if (acc_loc[i_dof] < (acc_cap * -1.0)) {
				acc_loc[i_dof] = (acc_cap * -1.0);
			}
		}

		// if (debug_mode_) {
		// 	// //get values related to delta variables now that it's post-calculation/post-cap code
		// 	// double dvel = vel_loc[i_dof] - vel_loc_pre[i_dof];
		// 	// double dacc = acc_loc[i_dof] - acc_loc_pre[i_dof];
		// 	// //accumulating for averages:
		// 	// avg_dvel += dvel;
		// 	// avg_dacc += dacc;
		// 	// avg_vel_post += vel_loc[i_dof];
		// 	// avg_acc_post += acc_loc[i_dof];
		// 	// //getting max and min values:
		// 	// if (dvel > max_dvel) {
		// 	// 	max_dvel = dvel;
		// 	// }
		// 	// if (dvel < min_dvel) {
		// 	// 	min_dvel = dvel;
		// 	// }
		// 	// if (dacc > max_dacc) {
		// 	// 	max_dacc = dacc;
		// 	// }
		// 	// if (dacc < min_dacc) {
		// 	// 	min_dacc = dacc;
		// 	// }

		// 	if (vel_loc[i_dof] > max_vel_post) {
		// 		max_vel_post = vel_loc[i_dof];
		// 	}
		// 	if (vel_loc[i_dof] < min_vel_post) {
		// 		min_vel_post = vel_loc[i_dof];
		// 	}
		// 	if (acc_loc[i_dof] > max_acc_post) {
		// 		max_acc_post = acc_loc[i_dof];
		// 	}
		// 	if (acc_loc[i_dof] < min_acc_post) {
		// 		min_acc_post = acc_loc[i_dof];
		// 	}
		// }

	}

	// if ( (debug_mode_) && (istep%report_step() == 0 ) )  {
	// if (debug_mode_)  {
	// 	// //getting averages
	// 	// avg_dvel = avg_dvel/(n_dof());
	// 	// avg_dacc = avg_dacc/(n_dof());
	// 	// avg_vel_pre = avg_vel_pre/(n_dof());
	// 	// avg_vel_post = avg_vel_post/(n_dof());
	// 	// avg_acc_pre = avg_acc_pre/(n_dof());
	// 	// avg_acc_post = avg_acc_post/(n_dof());

	// 	// //mbedit - reporting values
	// 	// TR << "mbedit_vvi nstep:" << istep
	// 	// << ":avg_dvel:" << avg_dvel << ":avg_dacc:" << avg_dacc << ":avg_vel_pre:" << avg_vel_pre 
	// 	// << ":avg_vel_post:" << avg_vel_post << ":avg_acc_pre:" << avg_acc_pre << ":avg_acc_post:" << avg_acc_post
	// 	// << ":max_dvel:" << max_dvel << ":min_dvel:" << min_dvel
	// 	// << ":max_vel_pre:" << max_vel_pre << ":min_vel_pre:" << min_vel_pre << ":max_vel_post:" << max_vel_post << ":min_vel_post:" << min_vel_post
	// 	// << ":max_dacc:" << max_dacc << ":min_dacc:" << min_dacc
	// 	// << ":max_acc_pre:" << max_acc_pre << ":min_acc_pre:" << min_acc_pre << ":max_acc_post:" << max_acc_post << ":min_acc_post:" << min_acc_post
	// 	// << std::endl;
	// 	// //mbedit

	// 	//mbedit - reporting values for debug_mode_ output
	// 	TR << "mbedit_vvi nstep:" << istep
	// 	<< ":max_vel_pre:" << max_vel_pre << ":min_vel_pre:" << min_vel_pre << ":max_vel_post:" << max_vel_post << ":min_vel_post:" << min_vel_post
	// 	<< ":max_acc_pre:" << max_acc_pre << ":min_acc_pre:" << min_acc_pre << ":max_acc_post:" << max_acc_post << ":min_acc_post:" << min_acc_post
	// 	<< std::endl;
	// 	//mbedit
	// }

	// Multivec xyz_loc_temp4( xyz_loc ), vel_loc_temp4( vel_loc ), acc_loc_temp4( acc_loc ); //mbedit line - storing pre-rattle2 values

	if ( use_rattle_ ) {
	
		rattle.run_rattle2( dt(), xyz_loc, vel_loc, mass() );

	}
	
	set_xyz( xyz_loc );
	set_vel( vel_loc );
	set_acc( acc_loc );

	//Stop rotation and translation
	//if ((step%nrottrans)==0) {
	//  stop_rot_trans(xyz, vel, size);
	//  thermostat_.get_temperature();
	//}
	// TR << "end of VVI fxn." << std::endl;
}

void CartesianMD::initialize_velocity( core::Real const &temperature, core::pose::PoseOP const pose_temp )
{
	//mbedit:
    using namespace basic::options;
    using namespace basic::options::OptionKeys;

	bool cap_on = option[basic::options::OptionKeys::md::cap_movements]; //default false
	Real vel_cap = option[basic::options::OptionKeys::md::max_vel]; //default 10.0
	Real acc_cap = option[basic::options::OptionKeys::md::max_acc]; //default 3e3

	if (cap_on) {
		TR << "capping on: max/min vel = " << vel_cap << "/" << (vel_cap * -1.0) << ", max/min acc = " << acc_cap << "/" << (acc_cap * -1.0) << std::endl;
	}
	//mbedit end


	TR << "Setting initial velocity with temp = " << temperature << std::endl;
	//mbedit_T -> set temperature to 300 here if xml flag on, check what init temp is after doing so, potentially change initial temp to 300 as well with flag

	// Supposed to make Maxwell-Boltzmann distribution... is this really working?
	// To make sure we should use error-function, but too lazy to do that...
	Multivec vel_loc( n_dof() ), acc_loc( n_dof(), 0.0 );
	for ( core::Size i_dof=1; i_dof<=n_dof(); i_dof++ ) {
		core::Size i_atm = (i_dof+2)/3;
		
		// pass Virtual atoms
		if ( mass(i_atm) < 1e-3 ) continue;
	
		Real scalar = sqrt(2.0*temperature*Boltzmann/mass(i_atm))*numeric::random::rg().gaussian(); //mbedit_T this should be 300 based on start of function
		if ( numeric::random::rg().uniform() > 0.5 ) scalar *= -1.0;

		vel_loc[i_dof] = scalar; //default behaviour is this line
		// vel_loc[i_dof] = 1.0; //mb debug, this is for debugging purposes

		//mbedit cap_on:
		if ( cap_on ) {
			if (vel_loc[i_dof] > vel_cap) {
				vel_loc[i_dof] = vel_cap;
			}

			if (vel_loc[i_dof] < (vel_cap * -1.0)) {
				vel_loc[i_dof] = (vel_cap * -1.0);
			}
		}
		//mbedit end

		//acc(i_dof) = 0.0;
		//printf("%4d %4d %8.3f %8.3f %8.3f\n",int(i_dof),int(i_atm),scalar,mass(i_atm),vel_loc[i_dof]);
	}

	//mb edit: place velocity matching code here that will equate ligA/ligB substructure atom velocities
	// std::cout << "init vel BEFORE idof matching:" << vel_loc << std::endl; //for debugging, comment out otherwise
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if (option[ basic::options::OptionKeys::md::fep_on ]) {
		// mapping / idof force update code, hardcoded for now! mb edit
		// utility::vector1<int> atm_map = {1, 68, 2, 67, 3, 60, 4, 59, 5, 57, 6, 56, 7, 63, 8, 62, 9, 51, 10, 53, 11, 49, 12, 61, 13, 65, 14, 64, 15, 58, 16, 66, 17, 50, 18, 54, 19, 55, 20, 52, 21, 47, 22, 48, 23, 46, 24, 44, 26, 84, 27, 85, 28, 86, 29, 81, 30, 82, 31, 83, 32, 78, 33, 77, 34, 76, 35, 75, 36, 80, 37, 79, 38, 72, 39, 71, 40, 74, 41, 73, 42, 70, 43, 69};
		// std::vector<int>
		// std::cout << "mb debug, atm_map[1]: " << atm_map[1] << ", atm_map[2]: " << atm_map[2] << ", atm_map.size(): " << atm_map.size() << ", atm_map[atm_map.size()]: " << atm_map[atm_map.size()] << ", atm_map[atm_map.size()-1]: " << atm_map[atm_map.size()-1] << std::endl;
		//calculate shift so that my hard-coded atom map makes sense regardless of how many non-ligand atoms are present.
		utility::vector1<int> atm_map = option[ basic::options::OptionKeys::md::md_atm_map ];
		int total_atoms = n_dof()/3;
		int total_lig_atoms = pose_temp->residue( pose_temp->size() - 1 ).natoms() + pose_temp->residue( pose_temp->size() ).natoms(); //adds natom number of ligA and ligB
		int idx_shift = total_atoms - total_lig_atoms;
		// std::cout << "mb debug, total_atoms:" << total_atoms << ":total_lig_atoms:" << total_lig_atoms << ":idx_shift:" << idx_shift << std::endl;
		for ( core::Size i = 1; i <= atm_map.size() - 1; i+=2 ) {
			int atm_a = atm_map[i] + idx_shift; //stores atom number for ligA atom
			int atm_b = atm_map[i + 1] + idx_shift; //stores atom number for ligB atom

			float x_a = vel_loc[(atm_a * 3) - 2]; //gets velocity for x-coord of atm_a
			float y_a = vel_loc[(atm_a * 3) - 1]; //gets velocity for y-coord of atm_a
			float z_a = vel_loc[(atm_a * 3)]; //gets velocity for z-coord of atm_a

			vel_loc[(atm_b * 3) - 2] = x_a; //sets velocity xcoord for ligB atom to vel of ligA atoms
			vel_loc[(atm_b * 3) - 1] = y_a; //sets velocity ycoord for ligB atom to vel of ligA atoms
			vel_loc[(atm_b * 3)] = z_a; //sets velocity zcoord for ligB atom to vel of ligA atoms
		}
	}
	// std::cout << "init vel AFTER idof matching:" << vel_loc << std::endl; //for debugging, comment out otherwise
	// Uniformly scale down to make sure init temperature assigned correctly
	Thermostat thermostat( temperature, n_dof_temp() );
	Real init_temp = thermostat.get_temperature( vel_loc, mass() );
	Real const scale( temperature/init_temp );
	// double max_vel, min_vel, avg_vel;
	// max_vel = min_vel = avg_vel = 0;

	for ( core::Size i_dof=1; i_dof<=n_dof(); i_dof++ ) {
		vel_loc[i_dof] *= scale;
		// if (debug_mode_) {
		// 	avg_vel += vel_loc[i_dof];
		// 	if (vel_loc[i_dof] > max_vel) {
		// 		max_vel = vel_loc[i_dof];
		// 	}
		// 	if (vel_loc[i_dof] < min_vel) {
		// 		min_vel = vel_loc[i_dof];
		// 	}
		// }
	}
	
	// if (debug_mode_) {
	// 	avg_vel = avg_vel/(n_dof());
	// 	TR << "mb_temp:max_scaled_vel:" << max_vel << ":min_scaled_vel:" << min_vel << ":avg_scaled_vel:" << avg_vel << std::endl;
	// }

	set_vel( vel_loc );
	Real init_temp2 = thermostat.get_temperature( vel_loc, mass() );

	TR << "mb_temp: initial temp after initial velocities chosen " << init_temp;
	TR << ", scaling down by factor " << scale << ". ";
	TR << "After scaling, init temp of system is now: " << init_temp2 << std::endl;
	// std::cout << "init vel:" << vel_loc << std::endl; //for debugging, comment out otherwise
}

void CartesianMD::report_MD( core::pose::Pose &pose,
	CartesianMinimizerMap const &min_map,
	bool const report_trj,
	core::Size istep )
{
	core::Real const rmsd( core::scoring::CA_rmsd( pose0(), pose ));

	// core::scoring::constraints::ConstraintSetCOP cstset( pose.constraint_set() );
	//TR << "Is there CST? " << std::endl;
	//cstset->show_numbers( TR );



// 	timeval currtime;
// #ifndef WIN32
// 	gettimeofday(&currtime, nullptr );
// #endif
// 	Real elapsedTime = (currtime.tv_sec - inittime_.tv_sec) * 1000.0;
// 	elapsedTime += (currtime.tv_usec - inittime_.tv_usec) / 1000.0;
// 	elapsedTime /= 60000.0; // in minute

	core::Real Epot( scorefxn()->score( pose ) );
	// std::cout << "total complex E breakdown:" << std::endl;
	// scorefxn()->show( pose ); //for debugging

	//mbedit - put code here to print specific energies; format is pose.energies().total_energies()[ sc_type ]:
	core::Real cb_angle ( pose.energies().total_energies_weighted()[ core::scoring::cart_bonded_angle ] );
	core::Real cb_length ( pose.energies().total_energies_weighted()[ core::scoring::cart_bonded_length ] );
	core::Real cb_torsion ( pose.energies().total_energies_weighted()[ core::scoring::cart_bonded_torsion ] );
	core::Real fa_rep ( pose.energies().total_energies_weighted()[ core::scoring::fa_rep ] );
	core::Real apc = 0;
	core::Real coordcst = 0;
	if ( (*scorefxn())[ core::scoring::coordinate_constraint ] != 0.0 ) {
		coordcst = pose.energies().total_energies_weighted()[ core::scoring::coordinate_constraint ];
	}

	if ( (*scorefxn())[ core::scoring::atom_pair_constraint ] != 0.0 ) {
		 apc = pose.energies().total_energies_weighted()[ core::scoring::atom_pair_constraint ];
	}

	core::Real fa_atr ( pose.energies().total_energies_weighted()[ core::scoring::fa_atr ] );
	core::Real fa_sol ( pose.energies().total_energies_weighted()[ core::scoring::fa_sol ] );
	
	core::Real fa_dun_dev ( pose.energies().total_energies_weighted()[ core::scoring::fa_dun_dev ] );
	core::Real fa_dun_rot ( pose.energies().total_energies_weighted()[ core::scoring::fa_dun_rot ] );
	core::Real fa_dun_semi ( pose.energies().total_energies_weighted()[ core::scoring::fa_dun_semi ] );

	core::Real lk_ball ( pose.energies().total_energies_weighted()[ core::scoring::lk_ball ] );
	core::Real lk_ball_iso ( pose.energies().total_energies_weighted()[ core::scoring::lk_ball_iso ] );
	core::Real lk_ball_bridge ( pose.energies().total_energies_weighted()[ core::scoring::lk_ball_bridge ] );
	core::Real lk_ball_bridge_uncpl ( pose.energies().total_energies_weighted()[ core::scoring::lk_ball_bridge_uncpl ] );
	
	core::Real hbond_sr_bb ( pose.energies().total_energies_weighted()[ core::scoring::hbond_sr_bb ] );
	core::Real hbond_lr_bb ( pose.energies().total_energies_weighted()[ core::scoring::hbond_lr_bb ] );
	core::Real hbond_bb_sc ( pose.energies().total_energies_weighted()[ core::scoring::hbond_bb_sc ] );
	core::Real hbond_sc ( pose.energies().total_energies_weighted()[ core::scoring::hbond_sc ] );

	
	// Getting interaction energies:
	core::Real proteinscore = 0.12345; //made this a weirdly suspicious number so that if somethings off, odds are I'll see it while debugging
	core::Real ligAscore = 0.12345;
	core::Real ligBscore = 0.12345;
	core::Real complexligA_score = 0.12345;
	core::Real ligA_apc = 0.12345;
	core::Real ligB_apc = 0.12345;
	core::Real complexligB_score = 0.12345;
	core::Real ligAinteraction_E = 0.12345; 
	core::Real ligBinteraction_E = 0.12345;

	//Getting scale factor controlling the sim to get un-scaled interaction energy
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	double scaleFactor_A = option[ corrections::score::elec_md_scale_factor ]; //scale factor to be used for ligand A unscaled intE calc
	double scaleFactor_B = (1.0 - scaleFactor_A); //scale factor to be used for ligand B unscaled intE calc
	
	if ( scaleFactor_A == 0.0 ) {
		scaleFactor_A = 1e-9;
	}

	if ( scaleFactor_B == 0.0 ) {
		scaleFactor_B = 1e-9;
	}

	double unscaled_ligA_intE = 0.12345;
	double unscaled_ligB_intE = 0.12345;

	if (calc_intE_) {
		proteinscore = calculate_free_receptor_score(pose);

		ligAscore = calculate_free_ligandA_score(pose);
		utility::vector1<double> temp_e1 = calculate_complexligA_score(pose);
		// std::cout << "calc intE temp_e1[1]:" << temp_e1[1] << ", temp_e1[2]:" << temp_e1[2] << std::endl;
		complexligA_score = temp_e1[1];
		ligA_apc = temp_e1[2];
		ligAinteraction_E = complexligA_score - proteinscore - ligAscore;
		unscaled_ligA_intE = ligAinteraction_E/scaleFactor_A;

		ligBscore = calculate_free_ligandB_score(pose);
		utility::vector1<double> temp_e2 = calculate_complexligB_score(pose);

		complexligB_score = temp_e2[1];
		ligB_apc = temp_e2[2];
		ligBinteraction_E = complexligB_score - proteinscore - ligBscore;
		unscaled_ligB_intE = ligBinteraction_E/scaleFactor_B;

	}

	//store user entered logfile path/file name
	std::string logfile_path = option[ basic::options::OptionKeys::md::md_outfile_name ];
	// std::cout << "logfile path is: " << logfile_path << std::endl;
	//check if the logfile exists already; if not, create it and write to it & if so, append to it
	std::ofstream ofs; //create outfile writer
	std::ifstream file(logfile_path); //returns true if file exists and is accessible
	if ( file.good() ) {
		// std::cout << "logfile exists!" << std::endl;
		ofs.open(logfile_path, std::ios::app);
		// ofs << "logfile exists so this line gets append. \n";
	} else {
		// std::cout << "logfile doesnt exist, writing to." << std::endl;
		ofs.open(logfile_path, std::ofstream::out);
		// ofs << "logfile doesnt exist so this gets written at the top. \n";
	}
	

	//TR if header_on_ true, print this one time header; else just do the values.
	if ( header_on_ ) {
		// TR << "md_labels nstep cb_ang cb_len cb_tor fa_rep fa_atr fa_sol lkb lkbi lkbb lkbbu fa_dun_dev fa_dun_rot fa_dun_semi hb_sr_bb hb_lr_bb hb_bb_sc hb_sc";
		ofs << "md_labels nstep cb_ang cb_len cb_tor fa_rep fa_atr fa_sol lkb lkbi lkbb lkbbu fa_dun_dev fa_dun_rot fa_dun_semi hb_sr_bb hb_lr_bb hb_bb_sc hb_sc";

		if ( (*scorefxn())[ core::scoring::atom_pair_constraint ] != 0.0 ) {
			// TR << " apc";
			ofs << " apc";
		}

		if ( (*scorefxn())[ core::scoring::coordinate_constraint ] != 0.0 ) {
			// TR << " coordcst";
			ofs << " coordcst";
		}

		if (calc_intE_) {
			// TR << " unscaled_ligA_intE unscaled_ligB_intE ligA_apc ligB_apc";
			ofs << " unscaled_ligA_intE unscaled_ligB_intE ligA_apc ligB_apc";

		}
		
		// TR << " time total_e temp rmsd";
		ofs << " time total_e temp rmsd\n";
	}

	//TR block to output the information; first pose then resi specific
	// TR << istep << " " << cb_angle << " " << cb_length << " " << cb_torsion << " " << fa_rep << " " << fa_atr << " " << fa_sol << " " << lk_ball << " " << lk_ball_iso << " " << lk_ball_bridge << " " << lk_ball_bridge_uncpl << " " << fa_dun_dev << " " << fa_dun_rot << " " << fa_dun_semi << " " << hbond_sr_bb << " " << hbond_lr_bb << " " << hbond_bb_sc << " " << hbond_sc;
	ofs << istep << " " << cb_angle << " " << cb_length << " " << cb_torsion << " " << fa_rep << " " << fa_atr << " " << fa_sol << " " << lk_ball << " " << lk_ball_iso << " " << lk_ball_bridge << " " << lk_ball_bridge_uncpl << " " << fa_dun_dev << " " << fa_dun_rot << " " << fa_dun_semi << " " << hbond_sr_bb << " " << hbond_lr_bb << " " << hbond_bb_sc << " " << hbond_sc;
	
	if ( (*scorefxn())[ core::scoring::atom_pair_constraint ] != 0.0 ) {
		// TR << " " << apc;
		ofs << " " << apc;
	}

	if ( (*scorefxn())[ core::scoring::coordinate_constraint ] != 0.0 ) {
		// TR << " " << coordcst;
		ofs << " " << coordcst;
	}

	if (calc_intE_) {
		// TR << " " << unscaled_ligA_intE << " " << unscaled_ligB_intE << " " << ligA_apc << " " << ligB_apc;
		ofs << " " << unscaled_ligA_intE << " " << unscaled_ligB_intE << " " << ligA_apc << " " << ligB_apc;
	}
	
	// TR << " " << cummulative_time() << " " << Epot << " " << temperature() << " " << rmsd;
	ofs << " " << cummulative_time() << " " << Epot << " " << temperature() << " " << rmsd;
	
	// core::Real rmsd_native( 0.0 ), gdttm_native( 0.0 ), gdtha_native( 0.0 );
	// if ( native_given_ ) {
	// 	rmsd_native = core::scoring::CA_rmsd( pose, native_, native_resmap_ );
	// 	core::scoring::CA_gdttm( pose, native_, gdttm_native, gdtha_native, native_resmap_ );

	// 	TR << ", RMSD/GDTtoNative ";
	// 	TR << std::setw(8) << std::setprecision(4) << rmsd_native;
	// 	TR << std::setw(8) << std::setprecision(4) << gdttm_native;
	// 	TR << std::setw(8) << std::setprecision(4) << gdtha_native;
	// 	ofs << ", RMSD/GDTtoNative ";
	// 	ofs << std::setw(8) << std::setprecision(4) << rmsd_native;
	// 	ofs << std::setw(8) << std::setprecision(4) << gdttm_native;
	// 	ofs << std::setw(8) << std::setprecision(4) << gdtha_native;
	// }

	// core::Real Eobj( scorefxn_obj()->score( pose ) );
	// TR << "Eobj: " << Eobj << " Emin_obj: " << Emin_obj() ;
	// ofs << "Eobj: " << Eobj << " Emin_obj: " << Emin_obj() ;
	// TR << std::endl; //final TR line ends here
	ofs << "\n";
	ofs.close(); //logfile gets closed here.

	// if ( cummulative_time() > 0.1 && // Truncate initial 1ps to remove minimization memory
	// 		selectmode() == "minobj" && Eobj < Emin_obj() ) {
	// 	set_pose_minobj( pose );
	// 	set_Emin_obj( scorefxn_obj()->score( pose ) );
	// 	set_time_minobj(  cummulative_time() );
	// 	TR << "Updating minimum objective score value / pose at time " << cummulative_time() << std::endl;
	// }

	// // Store trj
	// if ( report_trj && store_trj() ) {
	// 	//CartesianMinimizerMap min_map;
	// 	//min_map.setup( pose, *movemap() );

	// 	Multivec xyz_loc( min_map.ndofs() );
	// 	min_map.copy_dofs_from_pose( pose, xyz_loc );
	// 	add_trj( xyz_loc );
	// 	if ( report_as_silent() ) {
	// 		if ( native_given_ ) {
	// 			report_silent( pose, rmsd, gdttm_native, gdtha_native );
	// 		} else {
	// 			report_silent( pose );
	// 		}
	// 	}
	// }

	/*
	std::stringstream ss;
	core::Size timesize( (int)(cummulative_time()*1000.0) );
	ss << "md" << timesize << ".pdb";
	pose.dump_pdb( ss.str() );
	*/

	/*
	//Check
	for ( core::Size i_atm = 1; i_atm <= n_dof()/3; ++i_atm ) {
	Real vv = std::sqrt(vel_[3*i_atm-2]*vel_[3*i_atm-2] + vel_[3*i_atm-1]*vel_[3*i_atm-1] + vel_[3*i_atm]*vel_[3*i_atm]);
	Real av = std::sqrt(acc_[3*i_atm-2]*vel_[3*i_atm-2] + acc_[3*i_atm-1]*acc_[3*i_atm-1] + acc_[3*i_atm]*acc_[3*i_atm]);
	id::AtomID AtomID = min_map.get_atom( i_atm );
	core::Size resno = AtomID.rsd();
	core::Size atmno = AtomID.atomno();

	std::cout << "I/x/y/z/v/a: " << std::setw(6) << i_atm;
	std::cout << " " << std::setw(5) << pose.residue(resno).name();
	std::cout << " " << std::setw(4) << pose.residue(resno).atom_name(atmno);
	std::cout << " " << std::setw(10) << xyz_[3*i_atm-2];
	std::cout << " " << std::setw(10) << xyz_[3*i_atm-1];
	std::cout << " " << std::setw(10) << xyz_[3*i_atm];
	std::cout << " " << std::setw(10) << vv;
	std::cout << " " << std::setw(10) << av;
	std::cout << std::endl;
	}
	*/
} //end of report_MD()

// pose_ref should have exactly same molecule as stored in trj
utility::vector1< pose::Pose >
CartesianMD::dump_poses( pose::Pose const &pose_ref ) const {

	utility::vector1< pose::Pose > poses;

	// Don't know why min_map doesn't support const pose,
	// so let's just copy for now
	pose::Pose pose_tmp( pose_ref );

	// note that min_map is not const function
	CartesianMinimizerMap min_map;
	min_map.setup( pose_tmp, *movemap() );

	for ( core::Size i_trj = 1; i_trj <= trj().size(); ++i_trj ) {
		min_map.copy_dofs_to_pose( pose_tmp, trj(i_trj) );
		poses.push_back( pose_tmp );
	}

	return poses;
}

void  CartesianMD::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data
)
{
	// Other options
	parse_opts( tag, data );

	// Movemap
	parse_movemap( tag, data );
}

void CartesianMD::parse_opts(
	TagCOP const tag,
	basic::datacache::DataMap & data
) {
	using namespace core;
	using namespace scoring;


	std::string const scorefxn_name( tag->getOption< std::string >( "scorefxn" ) );
	set_scorefxn( data.get_ptr< ScoreFunction >( "scorefxns", scorefxn_name ) );

	std::string const scoreobj_name( tag->getOption< std::string >( "scorefxn_obj","" ) );

	if ( scoreobj_name == "" ) {
		set_scorefxn_obj( scorefxn() );
	} else {
		set_scorefxn_obj( data.get_ptr< ScoreFunction >( "scorefxns", scoreobj_name ) );
	}
	cst_beta_ = tag->getOption< bool >("cst_beta", false);
	use_rattle_ = tag->getOption< bool >("rattle", false );
	debug_mode_ = tag->getOption< bool >("debugMD", false);
	nve_mode_ = tag->getOption< bool >("nve_mode", false);
	calc_intE_ = tag->getOption< bool >("calc_intE", false);
	dumpstep_ =tag->getOption< core::Size >( "dumpstep", 100000 );
	// if ( use_rattle_ ) set_dt( 0.002 );
	// TR << "mb debug: use_rattle_ /debug_mode_/nve_mode_ values after tag get option:" << use_rattle_ << "/" << debug_mode_ << "/" << nve_mode_ << std::endl;

	//flags i added - mb:
	set_dt( tag->getOption< core::Real >("timestep", 0.002) );
	set_report_step( tag->getOption< core::Size >( "reportstep", 50 ) );

	set_nstep( tag->getOption< core::Size >( "nstep", 100 ) );
	set_temp0( tag->getOption<core::Real>("temp", 300.0) );
	set_scheduled( false );

	set_ncyc_premin( tag->getOption< core::Size >( "premin", 50 ) );
	set_ncyc_postmin( tag->getOption< core::Size >( "postmin", 200 ) );

	// set_md_report_stepsize( tag->getOption< core::Size >( "report", 100 ) );
	set_report_scorecomp( tag->getOption< bool >( "report_scorecomp", false ) ); //mbedit note: this just gives pre and post sim scoreterm breakdown; not during-run scoreterm breakdown
	set_selectmode( tag->getOption< std::string >( "selectmode", "final" ) );

	// Use parsed schedule file - this will overload nstep, temperature, etc. defined above
	std::string const schfile( tag->getOption< std::string >( "schfile","" ) );
	if ( schfile != "" ) {
		parse_schfile( schfile );
		set_scheduled( true );
	}
}

void CartesianMD::parse_movemap(
	TagCOP const tag,
	basic::datacache::DataMap & data
) {
	// set initial guess
	core::select::movemap::MoveMapFactoryOP mmf( new core::select::movemap::MoveMapFactory );
	bool const chi( tag->getOption< bool >( "chi", true ) ), bb( tag->getOption< bool >( "bb", true ) );
	mmf->all_chi( chi );
	mmf->all_bb( bb );

	movemap_factory_ = protocols::rosetta_scripts::parse_movemap_factory_legacy( tag, data, false, mmf );
	// std::cout << "in parse_movemap: movemap_factor_:" << movemap_factory_ << std::endl;
}

std::string CartesianMD::get_name() const {
	return mover_name();
}

std::string CartesianMD::mover_name() {
	return "CartesianMD";
}

void CartesianMD::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "timestep", xsct_real,
		"Sets the value of the timestep. Timestep is in ps. Default value is 0.002 ps (or 2 fs)." );
	attlist + XMLSchemaAttribute( "report_scorecomp", xsct_rosetta_bool,
		"shows scorefxn breakdown, useful for debugging. Default is false." );
	attlist + XMLSchemaAttribute( "cst_beta", xsct_rosetta_bool,
		"Turns on coord csts for all beta carbons in the pose. Default is false."
		"Note that the cst is harmonic, and the std dev is 1.0, therefore strength is tuned via CoordCst weight." );
	attlist + XMLSchemaAttribute( "dumpstep", xs_integer,
		"Sets the nstep interval at which a pdb will be dumped; so if dumpstep is 10, every 10 nsteps a pdb will dump."
		"Default is 100,000 because at timestep of 0.002, this dumps a pdb every 200 ps. This is an integer." );
	attlist + XMLSchemaAttribute( "reportstep", xs_integer,
		"When to report energies; so report happens every [reportstep] nsteps. "
		"If reportstep = 50, then report happens every 50 nsteps. This is an integer. 50 is the default." );
	attlist + XMLSchemaAttribute( "rattle", xsct_rosetta_bool,
		"Use Rattle algorithm to constraint hydrogen locations. "
		"Default is false." );
	attlist + XMLSchemaAttribute( "calc_intE", xsct_rosetta_bool,
		"boolean with default to false. When this is true, fxns are run inside reportMD()"
		"(which reports energies) in order to calculate just the interaction energy of the ligand to the protein." );
	attlist + XMLSchemaAttribute( "nve_mode", xsct_rosetta_bool,
		"Turns on NVE mode. Default is false. "
		"When this is true, the thermostat is allowed to run for initial velocity rescaling (aka only within the initalize_velocity fxn of CartesianMD.cc)"
		"but is afterwards 'frozen'; temperature is still measured but no rescaling of velocity occurs i.e. temperature is never 'set', just monitored. " );
	attlist + XMLSchemaAttribute( "debugMD", xsct_rosetta_bool,
		"Turns on debug mode. Default is false. "
		"When this is true, additional information is printed every reportstep, such as average changes in acc, vel, or position. "
		"Also gets additional score term values and prints THOSE energies such as fa_sol, lkball, hbond and other 2b terms. " );
	attlist + XMLSchemaAttribute( "scorefxn", xs_string,
		"Specify a scorefunction to run MD simulation with" );
	attlist + XMLSchemaAttribute( "scorefxn_obj", xs_string,
		"Optional, identical to scorefxn unless specified. "
		"Specify a scorefunction to use as objective function "
		"for selecting a pose from trajectory. "
		"This will be used only when selectmode=\"minobj\"" );
	attlist + XMLSchemaAttribute( "nstep", xs_integer,
		"Total number of steps to simulate. "
		"Since 1 nstep corresponds to a value of simulation time assigned to timestep, to estimate how many nsteps to run"
		"just do X ps / X timestep. Example: 'I want to run a 500 ps simulation with timestep 0.002."
		"500 ps / 0.002 = 250,000, therefore I should set nstep to 250,000!'. nstep is an integer." );
	attlist + XMLSchemaAttribute( "temp", xsct_real,
		"Reference temperature for constant temperature simulation. "
		"Recommended values: "
		"150~200K for talaris2014_cart and ~250 for beta_nov15_cart" );
	attlist + XMLSchemaAttribute( "premin", xs_integer,
		"Steps of Cartesian minimization before MD simulation" );
	attlist + XMLSchemaAttribute( "postmin", xs_integer,
		"Steps of Cartesian minimization after MD simulation" );
	attlist + XMLSchemaAttribute( "selectmode", xs_string,
		"How to select single pose from the trajectory. "
		"\"final\" to take the final pose, \"minobj\" to take the "
		"lowest objective function (by scorefxn_obj) pose" );
	attlist + XMLSchemaAttribute( "schfile", xs_string,
		"Use user-defined schedule file. "
		"This overrides any other flags or options. "
		"Syntax: \"sch [temperature] [nsteps]\" to run simulation, "
		"or \"repack\" to repack side-chains" );
	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"CartesianMD calls Molecular Dynamics simulation in Rosetta "
		"with user-defined energy function. Runs NVT simulation "
		"(constant volume and temperature) with Berendsen thermostat. "
		"Integrator uses Velocity Verlet algorithm", attlist );
}

std::string CartesianMDCreator::keyname() const {
	return CartesianMD::mover_name();
}

protocols::moves::MoverOP
CartesianMDCreator::create_mover() const {
	return utility::pointer::make_shared< CartesianMD >();
}

void CartesianMDCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CartesianMD::provide_xml_schema( xsd );
}


}
}
