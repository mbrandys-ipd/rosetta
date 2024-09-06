// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/md/thermostat.hh
/// @brief  Thermostat for MD simulation
/// @author Hahnbeom Park

#ifndef INCLUDED_protocols_md_thermostat_hh
#define INCLUDED_protocols_md_thermostat_hh

#include <protocols/md/MDConstants.hh>
#include <core/types.hh>
#include <core/optimization/types.hh>

#include <cmath>

namespace protocols {
namespace md {

class Thermostat {

public:

	Thermostat( core::Real const &temperature0, core::Size const ndof ){
		temp0_ = temperature0*Boltzmann;
		tau_ = 0.015;
		nstep_per_update_ = 10;
		ndof_ = ndof;
	}

	Thermostat( core::Real const &temperature0, core::Real const &tau, core::Size const ndof ){
		temp0_ = temperature0*Boltzmann;
		tau_ = tau;
		nstep_per_update_ = 10;
		ndof_ = ndof;
	}

	~Thermostat()= default;

	void
	rescale( 
		core::optimization::Multivec &vel, 
		core::Real const& dt, 
		core::optimization::Multivec const &mass,
		bool const nve_mode )
	{
		core::Real lambda = 1.0;
		
		//if nve_mode = false, aka temperature SHOULD be rescaled, calculate lambda (rescale value for velocity in this context)
		if (!nve_mode) {
			core::Real curr_temperature = get_temperature( vel, mass );
			core::Real delta  = (temp0_/Boltzmann)/curr_temperature;
			lambda = std::sqrt(1.0 + dt/tau_*(delta-1.0));
		}
		
		for ( core::Size i_dof = 1; i_dof <= vel.size(); ++i_dof ) {
			vel[i_dof] *= lambda;
		}
		// std::cout << "mb debug, inside thermostat rescale(); nve_mode is: " << nve_mode << "and rescale factor (lambda) is: " << lambda << std::endl;
	}

	core::Real
	get_temperature( core::optimization::Multivec const &vel, core::optimization::Multivec const &mass)
	{
		core::Real temperature = 0.0;

		for ( core::Size i_dof = 1; i_dof<=vel.size(); ++i_dof ) {
			int i_atm = (i_dof+2)/3;
			// pass Virtual atoms
			if ( mass[i_atm] < 1e-3 ) continue;

			temperature += mass[i_atm]*vel[i_dof]*vel[i_dof];
		}
		temperature /= ndof_*Boltzmann; //mbedit question: is this supposed to be temperature /= (ndof_*Boltzmann)?
		return temperature;
	}

	core::Size
	nstep_per_update(){ return nstep_per_update_; }

private:
	core::Real temp0_;
	core::Real tau_;
	core::Size nstep_per_update_;
	core::Size ndof_;
};
}
}

#endif
