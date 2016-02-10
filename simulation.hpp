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

namespace lb {
	
/**
 *  @brief Simulation class implementing LB
 * 
 *  This class holds a lattice as member (see @ref simulation::l) and 
 *  carries out the simulation steps on top of it. The main methods of 
 *  this class are @ref simulation::advect() and 
 *  @ref simulation::collide().
 */
class simulation
{
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
			// visc(_Vmax*nx/_Re ), // visc = Vmax*D/Re
			  visc(_Vmax*D/_Re),
			  beta( 1/(2*visc/(velocity_set().cs*velocity_set().cs) + 1) ),
			  time(0),
			  file_output(true), // set to true if you want to write files
			  output_freq(100),
			  output_index(0),
			  mDensityRho(1.0),
			  mHorizontalChannelVelocity(0.5)

	{ 
		// define amount to shift populations for advection
		for (unsigned int i=0; i<velocity_set().size; ++i)
		{

			shift[i] = velocity_set().c[0][i] + velocity_set().c[1][i]*l.real_nx;

			
		}
	}
	
	/**
	 *  @brief Initialize the flow field
	 *  
	 *  Initialization includes defining initial density, velocity and
	 *  populations. You can use Taylor-Green vortex flow conditions.
	 */
	void initialize()
	{
		// Variables
		const float_type pi(std::acos(-1.0));
		
		int lambdax = 1; 
		int lambday = 1;
		const float_type Kx = (2*pi)/(lambdax*l.nx); // Set L to nx
		const float_type Ky = (2*pi)/(lambday*l.nx); // Set L to nx
		const float_type K = std::sqrt(Kx*Kx + Ky*Ky);
		
		//Add time dependence to u,v,rho0 or not since it is the initial condition?
		// Initial velocities and density

		//  initialize nodes
		for (int j=-1; j<=static_cast<int>(l.ny); ++j)
		{
			for (int i=-1; i<=static_cast<int>(l.nx); ++i)
			{

				l.get_node(i,j).u()   = mHorizontalChannelVelocity; // -((Vmax*Ky/std::sqrt(Kx*Kx+Ky*Ky))*std::sin(Ky*j)*std::cos(Kx*i));
				l.get_node(i,j).v()   = 0;//((Vmax*Ky/std::sqrt(Kx*Kx+Ky*Ky))*std::sin(Kx*i)*std::cos(Ky*j));
				l.get_node(i,j).rho() = mDensityRho;// - (Vmax/velocity_set().cs)*(Vmax/velocity_set().cs)/(2*K*K)*(Ky*Ky*std::cos(2*Kx*i)+Kx*Kx*std::cos(2*Ky*j));
				lb::velocity_set().equilibrate(l.get_node(i,j));
								
			}
		}
		
	
	}
	
	/** 
	 *  @brief advect the populations
	 *  
	 *  Include periodic boundary conditions here also
	 */
	void advect()
	{		
		// **************************
		
		// Advection with shift				
			for (int k=0; k<velocity_set().size; k++)
			{							
				if (shift[k]<0)
				{
					for (int i = 0 ; i < l.real_size+shift[k]; i ++)
					{
						l.f[k][i] = l.f[k][i-shift[k]];
					}
				}
				else
				{
					for (int i = l.real_size-1; i >= 0+shift[k]; i--)
					{
						l.f[k][i] = l.f[k][i-shift[k]];
					}
				}
			}


			// Variables
			const float_type pi(std::acos(-1.0));

			int lambdax = 1;
			int lambday = 1;
			const float_type Kx = (2*pi)/(lambdax*l.nx); // Set L to nx
			const float_type Ky = (2*pi)/(lambday*l.nx); // Set L to nx
			const float_type K = std::sqrt(Kx*Kx + Ky*Ky);

			// Periodic boundaries
			// East moving particles
			for (int j = 0; j <= l.ny-1; j++)
			{
				l.get_node(0,j).u()   = mHorizontalChannelVelocity;// -((Vmax*Ky/std::sqrt(Kx*Kx+Ky*Ky))*std::sin(Ky*j)*std::cos(Kx*0));
				l.get_node(0,j).v()   = 0;//((Vmax*Ky/std::sqrt(Kx*Kx+Ky*Ky))*std::sin(Kx*i)*std::cos(Ky*j));
				l.get_node(0,j).rho() = mDensityRho;// - (Vmax/velocity_set().cs)*(Vmax/velocity_set().cs)/(2*K*K)*(Ky*Ky*std::cos(2*Kx*i)+Kx*Kx*std::cos(2*Ky*j));
				lb::velocity_set().equilibrate(l.get_node(0,j));

			}
//
			// West moving particles
			for (int j = 0; j <= l.ny-1; j++)
			{
				lb::float_type feqLocal[9];
				velocity_set().f_eq(feqLocal,l.get_node(l.nx-1,j).rho(),l.get_node(l.nx-1,j).u(),l.get_node(l.nx-1,j).v());

				for(int k=0; k<velocity_set().size; ++k)
				{

					l.get_node(l.nx-1,j).f(k) = feqLocal[k];
				}

			}

			// North moving particles
			for (int i = 0; i <= l.nx-1; ++i)
			{
			 l.f[4][l.index(i,l.ny-1)-shift[4]] = l.f[2][l.index(i,l.ny)]; // moving north
			 l.f[8][l.index(i,l.ny-1)-shift[8]] = l.f[5][l.index(i,l.ny)]; // moving NE
			 l.f[7][l.index(i,l.ny-1)-shift[7]] = l.f[6][l.index(i,l.ny)]; // moving NW
			}
			
			// South moving particles
			for (int i = 0; i <= l.nx-1; ++i)
			{
			 l.f[2][l.index(i,0)-shift[2]] = l.f[4][l.index(i,-1)]; // moving south
			 l.f[6][l.index(i,0)-shift[6]] = l.f[7][l.index(i,-1)]; // moving SW
			 l.f[5][l.index(i,0)-shift[5]] = l.f[8][l.index(i,-1)]; // moving SE
			}



			// Buffer CORNERS
			// ignored on purpose

		// **************************
		
	}
	
	/**  @brief apply wall boundary conditions */
	void wall_bc()
	{
		// fix_missing_populations()


		#pragma omp parallel for
		for (unsigned int i=0; i<l.wall_nodes.size(); ++i)
		{
			// **************************

			// **************************
		}
	}
	
	/** @brief collide the populations */
	void collide()
	{
		// **************************
		for(int i=0; i<=l.nx-1; ++i)
		{
			for(int j=0 ; j<=l.ny-1; ++j)
			{

				// Calculate new density
				l.get_node(i,j).rho() = 0;
				for(int k=0; k<velocity_set().size; ++k)
				{
					l.get_node(i,j).rho() +=  l.get_node(i,j).f(k); // l.f[k][l.index(i,j)]
					
				}
				// Calculate new velocities
				l.get_node(i,j).u() = 0;
				l.get_node(i,j).v() = 0;
				for(int k=0; k<velocity_set().size; ++k)
				{
					l.get_node(i,j).u() += l.get_node(i,j).f(k) * velocity_set().c[0][k];
					l.get_node(i,j).v() += l.get_node(i,j).f(k) * velocity_set().c[1][k];
					
				}
				
				l.get_node(i,j).u() /= l.get_node(i,j).rho() ;
				l.get_node(i,j).v() /= l.get_node(i,j).rho() ;
				
				// Collision step 
				lb::float_type feqLocal[9];
				velocity_set().f_eq(feqLocal,l.get_node(i,j).rho(),l.get_node(i,j).u(),l.get_node(i,j).v());			
				
				for(int k=0; k<velocity_set().size; ++k)
				{
							
					l.get_node(i,j).f(k) = l.get_node(i,j).f(k) + H_root(l.get_node(i,j))*beta*(feqLocal[k]-l.get_node(i,j).f(k));
				}
				
			}
		}
		// **************************
		
	}
	
	/** @brief LB step */
	void step()
	{

		advect();
		wall_bc();
		collide();

		// eval_forces();
		// move_shape();
		// find_and_fix_refill_nodes();


		// file io
		if ( file_output && ( ((time+1) % output_freq) == 0 || time == 0 ) )
		{
			write_fields();
			++output_index;
		}
		
		++time;
	}
	
public: // write to file

	/** write macroscopic variables to ascii file */
	void write_fields()
	{
 /*std::stringstream fns;
		fns << "output/data_" << std::setfill('0') << std::setw(4) << output_index << ".txt";
		l.write_fields(fns.str());*/
	}

public: // print

	/** print to output stream */
	friend std::ostream& operator<<(std::ostream& os, const simulation& sim)
	{
		os << "simulation parameters\n" 
		   << "---------------------\n";
		os << "domain: " << sim.l.nx << " x " << sim.l.ny << "\n";
		os << "Re:     " << sim.Re << "\n";
		os << "Vmax:   " << sim.Vmax << "\n";
		os << "visc:   " << sim.visc << "\n";
		os << "beta:   " << sim.beta << "\n";
		return os;
	}

public:	// for interaction with immersed shape

	void set_shape(geometry_2D* a_shape)
	{
		mSingleImmersedBody = a_shape;
	}

	void fix_missing_populations()
	{
		// boundary nodes iterator
		std::vector<lb::coordinate<int>> currentBoundary = mSingleImmersedBody->get_boundary_nodes();

		// clear target rho and u
		l.clear_rho_target_for_all_boundary_nodes();
		l.clear_rho_target_for_all_boundary_nodes();

		for (auto it = currentBoundary.begin(); it != currentBoundary.end(); it++)
		{
			// lattice representation indices for current node
			int iIndex = it->i; int jIndex = it->j;
			lb::coordinate<int> currentCoordinate(iIndex,jIndex);

			// find missing populations at current node
			std::vector<int> currentMissingPopulations = mSingleImmersedBody->find_missing_populations(currentCoordinate);

			//-- The following naming is in accordance with equation 14 from reference paper of Dorschner et al. -------//

			lb::coordinate<double> u_tgt(0.0,0.0);
			double rho_bb = 0;
			double rho_s = 0;
			std::vector<bool> rho_bb_has_added_index_i(lb::velocity_set().size,false);

			// iterate through all missing populations at current node
			for (auto mIt = currentMissingPopulations.begin(); mIt != currentMissingPopulations.end(); mIt++)
			{
				lb::coordinate<int> currentVelocity(lb::velocity_set().c[0][*mIt] , lb::velocity_set().c[1][*mIt]);		// get the velocity represented by the current missing population index
				lb::node fluidNeighborNode = l.get_node(iIndex + currentVelocity.i , jIndex + currentVelocity.j); 		// get the fluid node to interpolate velocity from
				lb::coordinate<double> adjacentFluidVelocity(fluidNeighborNode.u(),fluidNeighborNode.v());				// get its velocity

				double q_i =  mSingleImmersedBody->get_ray_length_at_intersection(currentCoordinate,*mIt)/lb::velocity_set().magnitude_c[*mIt];
				lb::coordinate<double> adjacentWallVelocity = mSingleImmersedBody->get_velocity_at_intersection(currentCoordinate,*mIt);

				u_tgt.i += (q_i * adjacentFluidVelocity.i + adjacentWallVelocity.i) / (1 + q_i);
				u_tgt.j += (q_i * adjacentFluidVelocity.j + adjacentWallVelocity.j) / (1 + q_i);
//				currentUTarget += mSingleImmersedBody->

				// work on rho bb
				// add the bounce back using the incoming to outgoing velocity swapper
				rho_bb += l.get_node(iIndex,jIndex).f(lb::velocity_set().incoming_velocity_to_outgoing_velocity(*mIt));
				rho_bb_has_added_index_i[*mIt] = true;			// set velocity has been added to rho_bb

				// work on rho s
				double ci_dot_u_wi = currentVelocity.i * adjacentWallVelocity.i + currentVelocity.j + adjacentWallVelocity.j;
				rho_s += lb::velocity_set().W[*mIt]*ci_dot_u_wi;
			}

			rho_s *= 6*mDensityRho;
			// add remaining densities to rho_bb
			for (int remaining_indices = 0; remaining_indices < rho_bb_has_added_index_i.size(); remaining_indices++)
			{
				if (!rho_bb_has_added_index_i[remaining_indices])
				{
					rho_bb += l.get_node(iIndex,jIndex).f(remaining_indices);
				}
			}
			// complete rho_tgt update
			l.set_rho_target_at_node(currentCoordinate,rho_bb+rho_s);

			u_tgt.i /= currentMissingPopulations.size(); u_tgt.j /= currentMissingPopulations.size();
			l.set_u_target_at_node(currentCoordinate,u_tgt);
		}


		//l.f[missing_pop_index][node_index] = ...
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

	geometry_2D* mSingleImmersedBody = nullptr; //An immersed object interacting with the fluid

	double mDensityRho;
	double mHorizontalChannelVelocity;

};

} // lb

#endif // LB_SIMULATION_HPP_INCLUDED
