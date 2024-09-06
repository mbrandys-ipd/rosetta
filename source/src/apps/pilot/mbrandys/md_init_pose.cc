// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW CoMotion, email: license@uw.edu.

/// @file   md_init_pose.cc
/// @brief  piloting app to prep double ligand pose for FEP md sims + atom map key gen
/// @author Marisa Brandys (mbrandys@uw.edu)

//devel
#include <devel/init.hh>

// core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/Element.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>


#include <basic/options/option_macros.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/Energies.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/import_pose/import_pose.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/pose/PDBInfo.hh>

#include <basic/Tracer.hh>

#include <core/conformation/Residue.hh>

// utility headers
#include <utility/vector1.hh>
#include <utility/file/FileName.fwd.hh>

// c++ headers
#include <string>
#include <protocols/jd2/JobDistributor.hh>

// protocol headers
#include <protocols/ligand_docking/GALigandDock/MCSAligner.hh>
#include <protocols/ligand_docking/GALigandDock/LigandConformer.hh>
#include <protocols/ligand_docking/GALigandDock/GALigandDock.hh>

//option stuff
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/md.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

static basic::Tracer TR("apps.md_init_pose");

//fxn to create mcs_aligner
protocols::ligand_docking::ga_ligand_dock::MCSAligner
get_mcs_aligner(
	core::pose::PoseOP const pose_template,
	int const ligid 
) {
	std::cout << "mb debug, pilot get_mcs_aligner start." << std::endl; //mb edit
	
	std::cout << "mb debug, get_mcs_aligner, dumping pose_template: mcs_pose_template_pilotapp.pdb" << std::endl; //mb edit
	pose_template->dump_pdb("getmcsfxn_pose_template_pilotapp.pdb");
	protocols::ligand_docking::ga_ligand_dock::MCSAlignerOptions aligner_options = protocols::ligand_docking::ga_ligand_dock::MCSAlignerOptions();

    aligner_options.perturb_rb = false;
    aligner_options.restype_to_rdmol_options.noImplicitHs = true;
    aligner_options.restype_to_rdmol_options.skipHs = false;

	std::cout << "mb debug pose_template size: " << pose_template->size() << ", ligid: " << ligid << std::endl;
	protocols::ligand_docking::ga_ligand_dock::MCSAligner aligner( *pose_template, ligid, aligner_options );
	return aligner;
}

int
main( int argc, char * argv [] ) {
	
	try {
		devel::init(argc, argv);
		using namespace core;

        std::cout << "mb debug, pilot app starting." << std::endl;

        core::pose::Pose input_pose;
        core::import_pose::pose_from_file(input_pose, basic::options::option[ basic::options::OptionKeys::in::file::s ]()[1], core::import_pose::PDB_file);

		//set up pose to be aligned aka protein + ligB only
		core::pose::Pose pose_to_align;
		pose_to_align.detached_copy(input_pose);

		//remove the ligands appropriately
		core::Size ligA_ID = (input_pose.total_residue() - 1); //store resnum of ligA which is second to last
		core::Size ligB_ID = (input_pose.total_residue()); //store resnum of ligB which is last

		pose_to_align.delete_residue_slow(ligA_ID);
		input_pose.delete_residue_slow(ligB_ID);

		//debug: dumping poses to make sure they're correct
		pose_to_align.dump_pdb("pose_to_align_cmpd67_pilotapp.pdb");
		input_pose.dump_pdb("input_pose_cmpd60_pilotapp.pdb");


		core::pose::PoseOP template_pose( new core::pose::Pose( input_pose ) ); //points to input_pose which at this point is just protein-ligA
		protocols::ligand_docking::ga_ligand_dock::MCSAligner mcs_aligner = get_mcs_aligner(template_pose, template_pose->total_residue());

		//here for ligandconformer to work i need to pass in a core::pose::PoseCOP pose, not a regular pose.
		core::pose::PoseOP aligned_pose( new core::pose::Pose( pose_to_align ) ); //make pointer to protein-ligB complex aka pose_to_align
		aligned_pose->energies().clear();
		aligned_pose->data().clear();
		
		utility::vector1 < core::Size > lig_resids; //just for gene creation later
		lig_resids.push_back( aligned_pose->total_residue() );

		utility::vector1< core::Size > movable_scs; //creating empty object so gene can be initialized
		std::cout << "mb debug, pilot app, gonna dump/print the items in: gene_test( aligned_pose, lig_resids, movable_scs, false, false ):" << std::endl; //mb edit

		aligned_pose->dump_pdb("aligned_pose_preMCSapply_pilotapp.pdb"); //mb edit
		std::cout << "mb debug:lig_resids: " << lig_resids << ":movable_scs:" << movable_scs << std::endl; //mb edit
		
		protocols::ligand_docking::ga_ligand_dock::LigandConformer gene_test( aligned_pose, lig_resids, movable_scs, false, false );

		// LigandConformer gene_initial( pose_working, lig_resids, movable_scs, freeze_ligand_backbone_, freeze_ligand_ );

		mcs_aligner.apply(gene_test);

		std::cout << "mb debug, pilot app, post apply step, dumping gene_test as pose_FINAL_mbdebug.pdb to see how it compares..." << std::endl; //mb edit
		core::pose::PoseOP pose_mbdebug( new core::pose::Pose() ); //mb edit
		gene_test.to_pose( pose_mbdebug ); //mb edit -> note, pose_mbdebug here has the final aligned ligB onto protein pose; aligned_pose does not.

		aligned_pose->dump_pdb( "pose_FINAL_aligned_pose.pdb" ); //mb edit
		pose_mbdebug->dump_pdb( "pose_FINAL_mbdebug.pdb" ); //mb edit
		pose_to_align.dump_pdb("pose_FINAL_pose_to_align.pdb");

		std::cout << "mb debug, pilot app end." << std::endl;

		//next to do: put new ligB back onto ligA to get the proper input
		//next to do: save map object so that we can use it in MD
        
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}