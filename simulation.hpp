/** 
 *  @file
 *  @author Fabian Bösch
 *  @brief simulation
 *
 *  Last edit: Feb 2. 2016. ful. Working
 */

#ifndef LB_SIMULATION_HPP_INCLUDED
#define LB_SIMULATION_HPP_INCLUDED

#include "H_root.hpp"
#include "lattice.hpp"
#include <sstream>
#include "quadrilateral_2D.hpp"
#include "quadrilateral_cantilever_2D.h"


namespace lb {
	
/**
 *  @brief Simulation class implementing LB
 * 
 *  This class holds a lattice as member (see @ref simulation::l) and 
 *  carries out the simulation steps on top of it. The main methods of 
 *  this class are @ref simulation::advect() and 
 *  @ref simulation::collide().
 */
class simulation {
public: // ctor

	/**
	 *  @brief Construct from domain size and flow parameters
	 *  @param[in] nx    extent in x direction
	 *  @param[in] ny    extent in y direction
	 *  @param[in] _Re   Reynolds number
	 *  @param[in] _Vmax mean flow velocity
	 */
	simulation(unsigned int nx, unsigned int ny, float_type _Re, float_type _Vmax, const double _D)

			: l(nx, ny),
			  shift(velocity_set().size),
			  Re(_Re),
			  Vmax(_Vmax),
			  D(_D),
			  visc(_Vmax * _D / _Re),
			  beta(1 / (2 * visc / (velocity_set().cs * velocity_set().cs) + 1)),
			  time(0),
			  file_output(true), // set to true if you want to write files
			  output_freq(10),
			  output_index(0) {
		// define amount to shift populations for advection
		for (unsigned int i = 0; i < velocity_set().size; ++i) {

			shift[i] = velocity_set().c[0][i] + velocity_set().c[1][i] * l.real_nx;


		}
	}

	/**
	 *  @brief Initialize the flow field
	 *  
	 *  Initialization includes defining initial density, velocity and
	 *  populations. You can use Taylor-Green vortex flow conditions.
	 */
	void initialize() {
		// Variables
		const float_type pi(std::acos(-1.0));

		int lambdax = 1;
		int lambday = 1;
		const float_type Kx = (2 * pi) / (lambdax * l.nx); // Set L to nx
		const float_type Ky = (2 * pi) / (lambday * l.nx); // Set L to nx
		const float_type K = std::sqrt(Kx * Kx + Ky * Ky);

		//Add time dependence to u,v,rho0 or not since it is the initial condition?
		// Initial velocities and density

		//  initialize nodes [channel flow]
		for (int j = -1; j <= static_cast<int>(l.ny); ++j) {
			for (int i = -1; i <= static_cast<int>(l.nx); ++i) {

				// double y = j - 0.5; //for zou-hhe
				// double L_y = l.ny-2; // for Zou-he
				l.get_node(i, j).u() = Vmax;
				l.get_node(i, j).v() = 0;
				l.get_node(i, j).rho() = 1;

				lb::velocity_set().equilibrate(l.get_node(i, j));
			}


		}

		/* *************
		* Add shapes 
		******************* */
		  l.add_to_shapes(new cylinder_2D(lb::coordinate<double>(10*D,10*D),0,D/2));
//		l.add_to_shapes(new quadrilateral_2D(lb::coordinate<double>(10*D-10,10*D-10), 0, D, D));
//		l.add_to_shapes(
//				new quadrilateral_cantilever_2D(lb::coordinate<double>(10 * D, 10 * D), 0, 10 * D, D / 10, 10, 1));


		fix_missing_populations();

	}

	/** 
	 *  @brief advect the populations
	 *  
	 *  Include periodic boundary conditions here also
	 */
	void advect() {

		// Store populations for quench
		std::vector<float_type> outlet3(l.ny);
		std::vector<float_type> outlet6(l.ny);
		std::vector<float_type> outlet7(l.ny);
		for (int j = 0; j <= l.ny - 1; j++) {
			outlet3[j] = l.f[3][l.index(l.nx - 1, j)];
			outlet6[j] = l.f[6][l.index(l.nx - 1, j)];
			outlet7[j] = l.f[7][l.index(l.nx - 1, j)];
		}


		// Advection with shift				
		for (int k = 0; k < velocity_set().size; k++) {
			if (shift[k] < 0) {
				for (int i = 0; i < l.real_size + shift[k]; i++) {
					l.f[k][i] = l.f[k][i - shift[k]];
				}
			}
			else {
				for (int i = l.real_size - 1; i >= 0 + shift[k]; i--) {
					l.f[k][i] = l.f[k][i - shift[k]];
				}
			}
		}
		// South Wall --> North moving particles
		for (int i = 0; i <= l.nx - 1; ++i) {
			l.f[4][l.index(i, l.ny - 1) - shift[4]] = l.f[2][l.index(i, l.ny)]; // moving north
			l.f[8][l.index(i, l.ny - 1) - shift[8]] = l.f[5][l.index(i, l.ny)]; // moving NE
			l.f[7][l.index(i, l.ny - 1) - shift[7]] = l.f[6][l.index(i, l.ny)]; // moving NW
		}
		// North Wall --> South moving particles
		for (int i = 0; i <= l.nx - 1; ++i) {
			l.f[2][l.index(i, 0) - shift[2]] = l.f[4][l.index(i, -1)]; // moving south
			l.f[6][l.index(i, 0) - shift[6]] = l.f[7][l.index(i, -1)]; // moving SW
			l.f[5][l.index(i, 0) - shift[5]] = l.f[8][l.index(i, -1)]; // moving SE
		}
		// Buffer CORNERS
		l.f[5][l.index(0, 0) - shift[5]] = l.f[5][l.index(l.nx, l.ny)]; // sw corner
		l.f[6][l.index(l.nx - 1, 0) - shift[6]] = l.f[6][l.index(-1, l.ny)]; // se corner
		l.f[7][l.index(l.nx - 1, l.ny - 1) - shift[7]] = l.f[7][l.index(-1, -1)]; // ne corner
		l.f[8][l.index(0, l.ny - 1) - shift[8]] = l.f[8][l.index(l.nx, -1)]; // NW corner
		// INLET & OUT BC
		for (int j = 0; j <= l.ny - 1; j++) {
			// double y = j - 0.5;
			// double L_y = l.ny-2;

			double u_inlet = Vmax; //Uniform flow
			double v_inlet = 0;
			double rho_inlet = 1;

			//inlet
			l.get_node(-1, j).u() = u_inlet;
			l.get_node(-1, j).v() = v_inlet;
			l.get_node(-1, j).rho() = rho_inlet;
			lb::velocity_set().equilibrate(l.get_node(-1, j));

			//outlet (quench)
			// 'refill ghost boundary nodes (right wall) with previous time step'
			l.f[3][l.index(l.nx - 1, j)] = outlet3[j];
			l.f[6][l.index(l.nx - 1, j)] = outlet6[j];
			l.f[7][l.index(l.nx - 1, j)] = outlet7[j];

			// // Zou-he Inlet (Microscopic BC)
			// 	l.f[1][l.index(-1,j)] = l.f[3][l.index(-1,j)]+2/3*rho_inlet*u_inlet;
			// 	l.f[5][l.index(-1,j)] = l.f[7][l.index(-1,j)] +	0.5*(l.f[4][l.index(-1,j)]-l.f[2][l.index(-1,j)]) +\
				// 							0.5*rho_inlet*v_inlet +1/6*rho_inlet*u_inlet ;
			// 	l.f[8][l.index(-1,j)] = l.f[6][l.index(-1,j)] - 0.5*(l.f[4][l.index(-1,j)] - l.f[2][l.index(-1,j)]) -\
				// 							0.5*rho_inlet*v_inlet + 1/6 * rho_inlet*u_inlet ;
			// //outlet
			// 	double v_outlet = 0;
			// 	double rho_outlet = 1;
			// 	double u_outlet =  -1.0 + (l.f[0][l.index(l.nx,j)] + l.f[2][l.index(l.nx,j)]+l.f[4][l.index(l.nx,j)] +\
				// 	 2.0*(l.f[1][l.index(l.nx,j)]+l.f[5][l.index(l.nx,j)]+l.f[8][l.index(l.nx,j)]))/rho_outlet ;

			// // Zou-he Outlet (Microscopic Oulet BC)
			// 	l.f[3][l.index(l.nx,j)] = l.f[1][l.index(l.nx,j)]-2/3*rho_outlet*u_outlet ;
			// 	l.f[7][l.index(l.nx,j)] = l.f[5][l.index(l.nx,j)]-0.5*(l.f[4][l.index(l.nx,j)]-l.f[2][l.index(l.nx,j)])-\
				// 							0.5*rho_outlet*v_outlet-rho_outlet*u_outlet/6 ;		
			// 	l.f[6][l.index(l.nx,j)] = l.f[8][l.index(l.nx,j)]+0.5*(l.f[4][l.index(l.nx,j)]-l.f[2][l.index(l.nx,j)])+\
				// 							0.5*rho_outlet*v_outlet-rho_outlet*u_outlet/6 ;	

		}

		// (no slip) top and bottom
		// include corners??
		// South wall
		// 	for (int i = 0; i <= l.nx-1; ++i)
		// 	{
		// 		l.f[2][l.index(i,0)-shift[2]] = l.f[4][l.index(i,-1)]; // moving north
		// 		l.f[5][l.index(i,0)-shift[5]] = l.f[8][l.index(i,-1)]; // moving NE
		// 		l.f[6][l.index(i,0)-shift[6]] = l.f[7][l.index(i,-1)]; // moving NW
		// 	}

		// // North Wall
		// 	for (int i = l.nx-1; i >=0 ; --i)
		// 	{
		// 		l.f[4][l.index(i,l.ny-1)-shift[4]] = l.f[2][l.index(i,l.ny)]; // moving south
		// 		l.f[8][l.index(i,l.ny-1)-shift[8]] = l.f[5][l.index(i,l.ny)]; // moving SW
		// 		l.f[7][l.index(i,l.ny-1)-shift[7]] = l.f[6][l.index(i,l.ny)]; // moving SE
		// 	}
		fix_missing_populations();

	}

	/**  @brief apply wall boundary conditions */
	void wall_bc() {
		// fix_missing_populations()
		// Set the type of simulation of objects (dynamic, static)
		std::string sim_type = "dynamic";

		for (auto i = l.shapes.begin(); i != l.shapes.end(); i++) {

			std::vector<lb::coordinate<int>> currentBoundaryNodes = (*i)->get_boundary_nodes();
			std::vector<lb::coordinate<int>> currentSolidNodes = (*i)->get_internal_nodes();

			if (sim_type == "dynamic") {
				// Rotate (all shapes)	
				// (*i)->set_orientation((*i)->get_orientation()+0.1);
				auto COM = (*i)->get_center_of_mass();
//				(*i)->create_sinusoidal_motion( sin(time * 2 * 3.141596265 / 10) * 0.005);
//				(*i)->set_center_of_mass(lb::coordinate<double>(COM.i+0.1,COM.j));

				double F = 0.6;
				double f_s = 0.165*Vmax/D;
				double y_max = 0.25*D;
				(*i)->set_linear_velocity(lb::coordinate<double>(0.0,y_max*2*3.141592*f_s*F*sin(time*2*3.141592*f_s*F)));
				(*i)->step_forward_one_step();



				for (auto j = currentSolidNodes.begin(); j != currentSolidNodes.end(); j++) {
					//Delete shapes in visualization (to refresh)
					lb::coordinate<int> a_coordinate(j->i, j->j);
					l.unset_is_wall_node(a_coordinate);

				}

//				for (auto j = currentBoundaryNodes.begin(); j != currentBoundaryNodes.end(); j++) {
//					lb::coordinate<int> a_coordinate(j->i, j->j);
//					l.unset_is_boundary_node(a_coordinate);
//					l.get_node(j->i, j->j).rho() = 1; //test visualization
//				}
//
//				for (auto j = l.refill_nodes.begin(); j != l.refill_nodes.end(); j++) {
//
//					l.get_node(j->first.first, j->first.second).rho() = 1;
//				}

				l.refill_nodes.clear();
				//check_node_status();
				// UPDATE -Find the new wall & boundary nodes
				(*i)->update_shape();
				currentBoundaryNodes = (*i)->get_boundary_nodes();
				std::vector<lb::coordinate<int>> newSolidNodes = (*i)->get_internal_nodes();

				for (auto j = newSolidNodes.begin(); j != newSolidNodes.end(); j++) {
					//redraw shapes
					lb::coordinate<int> a_coordinate(j->i, j->j);
					l.set_is_wall_node(a_coordinate);
				}

//				for (auto j = currentBoundaryNodes.begin(); j != currentBoundaryNodes.end(); j++) {
//					lb::coordinate<int> a_coordinate(j->i, j->j);
//					l.set_is_boundary_node(a_coordinate);
//					l.get_node(j->i, j->j).rho() = 100; //test visualization
//				}

				std::vector<bool> test_refill_nodes = l.find_refill_nodes(currentSolidNodes, newSolidNodes);

				for (auto j = 0; j < test_refill_nodes.size(); j++) {
					if (!test_refill_nodes[j]) {
						std::cout << "Has " << test_refill_nodes.size() << " refill nodes" << std::endl;
						int fluid_direction_index = ((*i)->find_missing_populations(currentSolidNodes[j])[0]);
						lb::coordinate<int> fluid_direction(lb::velocity_set().c[0][fluid_direction_index],lb::velocity_set().c[1][fluid_direction_index]);
						l.get_node(currentSolidNodes[j].i, currentSolidNodes[j].j).rho() = 1;
						l.get_node(currentSolidNodes[j].i, currentSolidNodes[j].j).u() = l.get_node(currentSolidNodes[j].i+fluid_direction.i,currentSolidNodes[j].j+fluid_direction.j).u();
						l.get_node(currentSolidNodes[j].i, currentSolidNodes[j].j).v() = l.get_node(currentSolidNodes[j].i+fluid_direction.i,currentSolidNodes[j].j+fluid_direction.j).v();
						lb::velocity_set().equilibrate(l.get_node(currentSolidNodes[j].i, currentSolidNodes[j].j));
						auto a_refill_node = lb::coordinate<int>(currentSolidNodes[j].i, currentSolidNodes[j].j);
						l.set_is_refill_node(a_refill_node);
					}

				}

			}

		}


#pragma omp parallel for
		for (unsigned int i = 0; i < l.wall_nodes.size(); ++i) {
			// **************************

			// **************************
		}
	}

	/** @brief collide the populations */
	void collide() {
		// **************************
		for (int i = 0; i <= l.nx - 1; ++i) {
			for (int j = 0; j <= l.ny - 1; ++j) {

				// Calculate new density
				l.get_node(i, j).rho() = 0;
				for (int k = 0; k < velocity_set().size; ++k) {
					l.get_node(i, j).rho() += l.get_node(i, j).f(k); // l.f[k][l.index(i,j)]

				}
				// Calculate new velocities
				l.get_node(i, j).u() = 0;
				l.get_node(i, j).v() = 0;
				for (int k = 0; k < velocity_set().size; ++k) {
					l.get_node(i, j).u() += l.get_node(i, j).f(k) * velocity_set().c[0][k];
					l.get_node(i, j).v() += l.get_node(i, j).f(k) * velocity_set().c[1][k];

				}

				l.get_node(i, j).u() /= l.get_node(i, j).rho();
				l.get_node(i, j).v() /= l.get_node(i, j).rho();

				// Collision step 
				lb::float_type feqLocal[9];
				velocity_set().f_eq(feqLocal, l.get_node(i, j).rho(), l.get_node(i, j).u(), l.get_node(i, j).v());

				for (int k = 0; k < velocity_set().size; ++k) {

					l.get_node(i, j).f(k) = l.get_node(i, j).f(k) +
											H_root(l.get_node(i, j)) * beta * (feqLocal[k] - l.get_node(i, j).f(k));
				}

			}
		}

		// **************************

	}

	/** @brief LB step */
	void step() {

		advect();
		wall_bc();
		force_evaluation(FxSingleBody,FySingleBody);
		CdSingleBody = get_aerodynamic_coeffcient_from_force(FxSingleBody);
		ClSingleBody = get_aerodynamic_coeffcient_from_force(FySingleBody);
		print_aerodynamic_info();

		collide();


		// file io
		if (file_output && (((time + 1) % output_freq) == 0 || time == 0)) {
			//write_fields();
			++output_index;
		}

		++time;
	}

public: // write to file

	/** write macroscopic variables to ascii file */
	void write_fields() {
		std::stringstream fns;
		fns << "output/data_" << std::setfill('0') << std::setw(4) << output_index << "Re_" << Re << ".txt";
		l.write_fields(fns.str());
	}

public: // print

	/** print to output stream */
	friend std::ostream &operator<<(std::ostream &os, const simulation &sim) {
		os << "simulation parameters\n"
		<< "---------------------\n";
		os << "domain: " << sim.l.nx << " x " << sim.l.ny << "\n";
		os << "Re:     " << sim.Re << "\n";
		os << "Vmax:   " << sim.Vmax << "\n";
		os << "visc:   " << sim.visc << "\n";
		os << "beta:   " << sim.beta << "\n";
		return os;
	}

public:    // for interaction with immersed shape

	void set_shape(geometry_2D *a_shape) {
		mSingleImmersedBody = a_shape;
	}

	void check_node_status() {
		// boundary nodes iterator	
		auto mSingleImmersedBody = *l.shapes.begin(); // More general
		auto WallNodes = mSingleImmersedBody->get_internal_nodes();

		for (auto wall_iter = WallNodes.begin(); wall_iter <= WallNodes.end(); wall_iter++) {
			int iIndex = wall_iter->i;
			int jIndex = wall_iter->j;
			lb::coordinate<int> currentCoordinate(iIndex, jIndex);

			if (l.get_node(iIndex, jIndex).has_flag_property("fluid") ||
				l.get_node(iIndex, jIndex).has_flag_property("boundary")) {
				// l.set_is_refill_node(currentCoordinate);
				//l.get_node(wall_iter->i, wall_iter->j).rho() = 100;
			}

		}
	}
		void force_evaluation(double &Fx, double &Fy) {

			Fx = 0;
			Fy = 0;
			auto mSingleImmersedBody = *l.shapes.begin(); // More general
			auto currentBoundary = mSingleImmersedBody->get_boundary_nodes();

			for (auto it = currentBoundary.begin(); it != currentBoundary.end(); it++) {
				int iIndex = it->i;
				int jIndex = it->j;
				lb::coordinate<int> currentCoordinate(iIndex, jIndex);
				std::vector<int> currentMissingPopulations = mSingleImmersedBody->find_missing_populations(currentCoordinate);

				for (auto mIt = currentMissingPopulations.begin(); mIt != currentMissingPopulations.end(); mIt++) {
					lb::coordinate<int> currentVelocity(lb::velocity_set().c[0][*mIt], lb::velocity_set().c[1][*mIt]);
					lb::coordinate<int> adjacent_solidnode(iIndex - currentVelocity.i, jIndex - currentVelocity.j);

					Fx += -lb::velocity_set().c[0][*mIt] *
						  (l.get_node(adjacent_solidnode.i, adjacent_solidnode.j).f(lb::velocity_set().incoming_velocity_to_outgoing_velocity(*mIt))
						   + l.get_node(iIndex, jIndex).f(*mIt));
					Fy += -lb::velocity_set().c[1][*mIt] *
						  (l.get_node(adjacent_solidnode.i, adjacent_solidnode.j).f(lb::velocity_set().incoming_velocity_to_outgoing_velocity(*mIt))
						   + l.get_node(iIndex, jIndex).f(*mIt));
				}
			}


			// STORE Fx and Fy in data file for post-processing
		}

		double get_aerodynamic_coeffcient_from_force(double &Force) {
			double result = Force;
			result /= (0.5 * 1 * Vmax * Vmax * D);

			return result;
		}

		void print_aerodynamic_info() {

			std::cout << "At timestep: " << time << std::endl;
			std::cout << '\t' << "Drag coefficient of: " << CdSingleBody << '\n';
			std::cout << '\t' << "Lift coefficient of: " << ClSingleBody << '\n';
			std::cout << '\t' << "Change in magnitude from previous step: " <<
			std::abs(CdSingleBody - CdSingleBodyPreviousTimeStep) << std::endl;

			std::vector<int> recirc_length(1, ceil(10 * D + D / 2));
			for (int i = ceil(10 * D + D / 2); i <= l.nx - 1; i++) {
				if (l.get_node(i, 10 * D).u() <= 0) {
					recirc_length.push_back(i);
				}
			}

			double non_dim_length = 2 * (recirc_length.back() - ceil(10 * D + D / 2)) / D;
			CdSingleBodyPreviousTimeStep = CdSingleBody;
			ClSingleBodyPreviousTimeStep = ClSingleBody;
			std::cout << "Recirculation length of : " << non_dim_length << std::endl;

		}


		void fix_missing_populations() {
			string approx = "grads";

			// boundary nodes iterator
			auto mSingleImmersedBody = *l.shapes.begin(); // More general
			auto currentBoundary = mSingleImmersedBody->get_boundary_nodes();

			// mSingleImmersedBody->print_boundary_nodes(); // Test boundary nodes

			// get_boundary_nodes();
			// // clear target rho and u
			l.clear_rho_target_for_all_boundary_nodes();
			l.clear_u_target_for_all_boundary_nodes();

			newRhoMap.clear();
			newUMap.clear();
			newVMap.clear();
			// set new densities and velocities for boundary nodes and their fluid neighbours
			for (auto it = currentBoundary.begin(); it != currentBoundary.end(); it++) {
				//lattice representation indices for current node
				int iIndex = it->i;
				int jIndex = it->j;
				lb::coordinate<int> currentCoordinate(iIndex, jIndex);

				// find missing populations at current node
				std::vector<int> currentMissingPopulations = mSingleImmersedBody->find_missing_populations(
						currentCoordinate);

				for (auto mIt = currentMissingPopulations.begin(); mIt != currentMissingPopulations.end(); mIt++)
				{

					lb::coordinate<int> currentVelocity(lb::velocity_set().c[0][*mIt], lb::velocity_set().c[1][*mIt]);        // get the velocity represented by the current missing population index

					// Bounce-back (simple implementation)
					l.get_node(iIndex, jIndex).f(*mIt) = l.get_node(iIndex, jIndex).f(lb::velocity_set().incoming_velocity_to_outgoing_velocity(*mIt));

					if (newRhoMap.count(std::make_pair(iIndex + currentVelocity.i, jIndex + currentVelocity.j))) continue;

					double newRho = 0; double newU = 0; double newV = 0;

					for (int veloUpdateIndex = 0; veloUpdateIndex < 9; veloUpdateIndex++)
					{
						newRho += l.get_node(iIndex + currentVelocity.i, jIndex + currentVelocity.j).f(veloUpdateIndex);
						newU += l.get_node(iIndex + currentVelocity.i, jIndex + currentVelocity.j).f(veloUpdateIndex)*lb::velocity_set().c[0][veloUpdateIndex];
						newV += l.get_node(iIndex + currentVelocity.i, jIndex + currentVelocity.j).f(veloUpdateIndex)*lb::velocity_set().c[1][veloUpdateIndex];
					}

					newU /= newRho;
					newV /= newRho;

					newRhoMap[std::make_pair(iIndex + currentVelocity.i, jIndex + currentVelocity.j)] = newRho;
					newUMap[std::make_pair(iIndex + currentVelocity.i, jIndex + currentVelocity.j)] = newU;
					newVMap[std::make_pair(iIndex + currentVelocity.i, jIndex + currentVelocity.j)] = newV;
				}

			}


			for (auto it = currentBoundary.begin(); it != currentBoundary.end(); it++) {
				//lattice representation indices for current node
				int iIndex = it->i;
				int jIndex = it->j;
				lb::coordinate<int> currentCoordinate(iIndex, jIndex);

				// find missing populations at current node
				std::vector<int> currentMissingPopulations = mSingleImmersedBody->find_missing_populations(
						currentCoordinate);


				//-- The following naming is in accordance with equation 14 from reference paper of Dorschner et al. -------//
				if (approx == "grads") {
					lb::coordinate<double> u_tgt(0.0, 0.0);
					std::vector<bool> rho_bb_has_added_index_i(lb::velocity_set().size, false);
					double rho_s = 0;
					double rho_bb = 0;

					lb::node current_node = l.get_node(iIndex, jIndex);
					double dvdy = 0;
					double dvdx = 0;
					double dudx = 0;
					double dudy = 0;

					// iterate through all missing populations at current node
					for (auto mIt = currentMissingPopulations.begin(); mIt != currentMissingPopulations.end(); mIt++) 
					{

						lb::coordinate<int> currentVelocity(lb::velocity_set().c[0][*mIt], lb::velocity_set().c[1][*mIt]);        // get the velocity represented by the current missing population index
						lb::node fluidNeighborNode = l.get_node(iIndex + currentVelocity.i, jIndex + currentVelocity.j);        // get the fluid node to interpolate velocity from
						lb::coordinate<double> adjacentFluidVelocity(fluidNeighborNode.u(), fluidNeighborNode.v());                // get its velocity
						lb::coordinate<double> newAdjacentFluidVelocity(newUMap[make_pair(iIndex + currentVelocity.i, jIndex + currentVelocity.j)],
																		newVMap[make_pair(iIndex + currentVelocity.i, jIndex + currentVelocity.j)]);

						double q_i = (mSingleImmersedBody->get_ray_length_at_intersection(currentCoordinate, *mIt)) / lb::velocity_set().magnitude_c[*mIt];
						lb::coordinate<double> adjacentWallVelocity = mSingleImmersedBody->get_velocity_at_intersection(currentCoordinate, *mIt);

						u_tgt.i += (q_i * newAdjacentFluidVelocity.i + adjacentWallVelocity.i) / (1 + q_i);
						u_tgt.j += (q_i * newAdjacentFluidVelocity.j + adjacentWallVelocity.j) / (1 + q_i);

						// work on rho bb

						// add the bounce back using the incoming to outgoing velocity swapper
						rho_bb += l.get_node(iIndex, jIndex).f(lb::velocity_set().incoming_velocity_to_outgoing_velocity(*mIt));
						rho_bb_has_added_index_i[*mIt] = true;            // set velocity has been added to rho_bb

						// work on rho s
						double ci_dot_u_wi = currentVelocity.i * adjacentWallVelocity.i + currentVelocity.j * adjacentWallVelocity.j;

						rho_s += lb::velocity_set().W[*mIt] * ci_dot_u_wi;

						//velocity gradient
						if (currentVelocity.i == 0) 
						{
							dvdy = currentVelocity.j * (adjacentFluidVelocity.j - current_node.v());
							dudy = currentVelocity.j * (adjacentFluidVelocity.i - current_node.u());
						}
						else if (currentVelocity.j == 0) 
						{
							dudx = currentVelocity.i * (adjacentFluidVelocity.i - current_node.u());// velocity comp. 'i', derived in 'i' direction
							dvdx = currentVelocity.i * (adjacentFluidVelocity.j - current_node.v());
						}
					}

					double rho_0 = 0;

					for (int rho_0_index = 0; rho_0_index < 9; rho_0_index++)
					{
						rho_0 += l.get_node(iIndex, jIndex).f(rho_0_index);
					}

					rho_s *= 6 * rho_0;

					// add remaining densities to rho_bb
					for (int remaining_indices = 0;
						 remaining_indices < rho_bb_has_added_index_i.size(); remaining_indices++)
					{
						if (!rho_bb_has_added_index_i[remaining_indices])
						{
							rho_bb += l.get_node(iIndex, jIndex).f(remaining_indices);
						}
					}
					// complete rho_tgt update
					double rho_tgt = rho_bb + rho_s;
					l.set_rho_target_at_node(currentCoordinate, rho_tgt);
					//complete u_tgt update
					u_tgt.i /= currentMissingPopulations.size();
					u_tgt.j /= currentMissingPopulations.size();
					// l.set_u_target_at_node(currentCoordinate, u_tgt);

					// Pressure tensor
					double Pxxeq = rho_tgt * lb::velocity_set().cs * lb::velocity_set().cs + rho_tgt * u_tgt.i * u_tgt.i;
					double Pyyeq = rho_tgt * lb::velocity_set().cs * lb::velocity_set().cs + rho_tgt * u_tgt.j * u_tgt.j;
					double Pxyeq = rho_tgt * u_tgt.i * u_tgt.j; // also --> Pyxeq

					double Pxxneq = -rho_tgt * lb::velocity_set().cs * lb::velocity_set().cs / (2 * beta) * (dudx + dudx);
					double Pyyneq = -rho_tgt * lb::velocity_set().cs * lb::velocity_set().cs / (2 * beta) * (dvdy + dvdy);
					double Pxyneq = -rho_tgt * lb::velocity_set().cs * lb::velocity_set().cs / (2 * beta) * (dudy + dvdx);

					double Pxx = Pxxeq + Pxxneq;
					double Pyy = Pyyeq + Pyyneq;
					double Pxy = Pxyeq + Pxyneq;

					// Grad's Approximation
					for (auto mIt = currentMissingPopulations.begin(); mIt != currentMissingPopulations.end(); mIt++) 
					{
						lb::coordinate<int> currentVelocity(lb::velocity_set().c[0][*mIt],lb::velocity_set().c[1][*mIt]);

						l.f[*mIt][l.index(iIndex, jIndex)] = lb::velocity_set().W[*mIt] *
						(rho_tgt + 
						rho_tgt * u_tgt.i * currentVelocity.i / (lb::velocity_set().cs * lb::velocity_set().cs) + 
						rho_tgt * u_tgt.j * currentVelocity.j /(lb::velocity_set().cs * lb::velocity_set().cs) +
					   	(1 /(2 * lb::velocity_set().cs *lb::velocity_set().cs * lb::velocity_set().cs *lb::velocity_set().cs)) *
					   	((Pxx - rho_tgt *lb::velocity_set().cs *lb::velocity_set().cs) *(currentVelocity.i *currentVelocity.i -lb::velocity_set().cs *lb::velocity_set().cs) +
						(Pyy - rho_tgt *lb::velocity_set().cs * lb::velocity_set().cs) *(currentVelocity.j *currentVelocity.j -lb::velocity_set().cs *lb::velocity_set().cs) +
						2 * (Pxy) *(currentVelocity.i *currentVelocity.j)));
					}


				}
			}


		}



		public: // members

		lattice l;                 ///< lattice
		std::vector<int> shift;    ///< amount of nodes to shift each population in data structure during advection
		const float_type Re;       ///< Reynolds number
		const float_type Vmax;     ///< mean flow velocity
		const float_type visc;     ///< viscosity
		const float_type beta;     ///< LB parameter beta
		const double D;
		unsigned int time;         ///< simulation time
		bool file_output;          ///< flag whether to write files
		unsigned int output_freq;  ///< file output frequency
		unsigned int output_index; ///< index for file naming
		double FxSingleBody;        // reaction force on a single submerged body
		double FySingleBody;
		double CdSingleBody;
		double ClSingleBody;
		double CdSingleBodyPreviousTimeStep = 0;
		double ClSingleBodyPreviousTimeStep = 0;
		std::map<std::pair<int,int>,double> newRhoMap;
		std::map<std::pair<int,int>,double> newUMap;
		std::map<std::pair<int,int>,double> newVMap;



		geometry_2D *mSingleImmersedBody = nullptr; //An immersed object interacting with the fluid


	};


} // lb

#endif // LB_SIMULATION_HPP_INCLUDED

