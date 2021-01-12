/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

// declare cell definitions here


void create_cell_types( void )
{

	// housekeeping

	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

	// Name the default cell type

	cell_defaults.type = 0;
	cell_defaults.name = "motile tumor cell";

	// set default cell cycle model

	cell_defaults.functions.cycle_model = Ki67_basic;

	// set default_cell_functions;
	cell_defaults.functions.update_phenotype = update_cell_and_death_parameters_O2_based;

	// only needed for a 2-D simulation:

	/*
	cell_defaults.functions.set_orientation = up_orientation;
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true;
	*/

	// make sure the defaults are self-consistent.

	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.molecular.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions );

	// set the rate terms in the default phenotype

	// first find index for a few key variables.
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" );

	// set oxygen uptake / secretion parameters for the default cell type
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 10;
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 38;

	// cell velocity is influenced by drag forces
	cell_defaults.functions.update_velocity = drag_update_cell_velocity;

	// enable random motility
	cell_defaults.phenotype.motility.is_motile = true;
	cell_defaults.phenotype.motility.migration_bias = 0.0; // completely random

	cell_defaults.phenotype.motility.persistence_time =
		parameters.doubles( "cell_persistence_time" );

	double initial_velocity = locomotive_forces_generator();
	cell_defaults.phenotype.motility.migration_speed = initial_velocity;

	// mechanics parameters
	cell_defaults.phenotype.mechanics.cell_cell_repulsion_strength =
	   parameters.doubles( "cell_cell_cell_repulsion" );

	cell_defaults.phenotype.mechanics.cell_cell_adhesion_strength =
	   parameters.doubles( "cell_cell_cell_adhesion" );

	cell_defaults.phenotype.mechanics.relative_maximum_adhesion_distance =
	   parameters.doubles("cell_rel_max_adhesion_distance");


	build_cell_definitions_maps();
	display_cell_definitions( std::cout );

	return;
}


double locomotive_forces_generator( )
{
	// random number generator to define cell locomotive forces
	// based on anempirically obtained velocity distribution

	double random_value, force_value;

	random_value = UniformRandom();

  force_value =
		0.1157 * pow(random_value, 3) + 0.2422 * pow(random_value, 2) +
		0.0053 * random_value + 0.0048;

	force_value = 13.5*force_value;

	return force_value;
}


void drag_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt )
{

	// find location of variables and base parameter values
	int ECM_density_index = microenvironment.find_density_index( "ECM" );

	// sample ECM
	double ECM_density = pCell->nearest_density_vector()[ECM_density_index];
	double dyn_viscosity;

	// get viscosity based on concentration
	if(ECM_density == 2.5)
	{
		dyn_viscosity = 7.96;
	}
	else if(ECM_density == 4.0)
	{
		dyn_viscosity = 18.42;
	}
	else if(ECM_density == 6.0)
	{
		dyn_viscosity = 39.15;
	}

	// update the speed value
	pCell->phenotype.motility.migration_speed = locomotive_forces_generator();

  // update velocity
	standard_update_cell_velocity(pCell, phenotype, dt);

	// include the 1/vu (1/ECM density) term to consider friction
	pCell->velocity /= dyn_viscosity;

	return;

}

void setup_microenvironment( void )
{
	// set domain parameters

/*
	default_microenvironment_options.X_range = {-500, 500};
	default_microenvironment_options.Y_range = {-500, 500};
	default_microenvironment_options.Z_range = {-500, 500};
*/
	// make sure to override and go back to 2D
	if( default_microenvironment_options.simulate_2D == true )
	{
		std::cout << "Warning: overriding XML config option and setting to 3D!" << std::endl;
		default_microenvironment_options.simulate_2D = false;
	}


/*
	// all this is in XML as of August 2019 (1.6.0)
	// no gradients need for this example

	default_microenvironment_options.calculate_gradients = false;

	// set Dirichlet conditions

	default_microenvironment_options.outer_Dirichlet_conditions = true;

	// if there are more substrates, resize accordingly
	std::vector<double> bc_vector( 1 , 38.0 ); // 5% o2
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;

	// set initial conditions
	default_microenvironment_options.initial_condition_vector = { 38.0 };
*/

	// initialize BioFVM

	initialize_microenvironment();

	return;
}

void setup_tissue( void )
{

	// create spaced cells
	double cell_radius = cell_defaults.phenotype.geometry.radius;
	double cell_spacing = 32.0 * cell_radius;

	Cell* pC;

	for(double x_value = -cell_spacing; x_value <= cell_spacing; x_value += cell_spacing)
	{
		for(double y_value = -cell_spacing; y_value <= cell_spacing; y_value += cell_spacing)
		{
			pC = create_cell(  );
			pC->assign_position( x_value, y_value, 0.0 );
		}
	}

	return;
}
