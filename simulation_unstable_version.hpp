/** 
 *  @file
 *  @author Fabian Bösch
 *  @brief simulation
 */

#ifndef LB_SIMULATION_HPP_INCLUDED
#define LB_SIMULATION_HPP_INCLUDED

#include "H_root.hpp"
#include "lattice.hpp"
#include <sstream>

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
		simulation(unsigned int nx, unsigned int ny, float_type _Re, float_type _Vmax)
				: l(nx, ny),
				  shift(velocity_set().size),
				  Re(_Re),
				  Vmax(_Vmax),
				  visc(_Vmax*nx/_Re ), //
				  beta( 1/(2*visc/(velocity_set().cs*velocity_set().cs) + 1) ),
				  time(0),
				  file_output(false), // set to true if you want to write files
				  output_freq(100),
				  output_index(0)
		{
			// define amount to shift populations for advection
			for (unsigned int i=0; i<velocity_set().size; ++i)
			{
				// **************************
				// * fill in your code here *
				// Can be positive or negative --> determines the loop order
				// **************************
				shift[i] = 0;
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
			// **************************
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
			for (int j=0; j<static_cast<int>(l.ny); ++j)
			{
				for (int i=0; i<static_cast<int>(l.nx); ++i)
				{
					l.get_node(i,j).u()   =  -((Vmax*Ky/std::sqrt(Kx*Kx+Ky*Ky))*std::sin(Ky*j)*std::cos(Kx*i));
					l.get_node(i,j).v()   = ((Vmax*Ky/std::sqrt(Kx*Kx+Ky*Ky))*std::sin(Kx*i)*std::cos(Ky*j));
					l.get_node(i,j).rho() = 1 - Vmax/velocity_set().cs*Vmax/velocity_set().cs/(2*K*K)*(Ky*Ky*std::cos(2*Kx*i)+Kx*Kx*std::cos(2*Ky*j));
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
			// Periodic boundaries
			// East moving particles
			for (int j = 0; j <= l.ny-1; j++)
			{
				l.f[1][l.index(-1,j)] = l.f[1][l.index(l.nx-1,j)];	//moving east
				l.f[5][l.index(-1,j)] = l.f[5][l.index(l.nx-1,j)];	// moving NE
				l.f[8][l.index(-1,j)] = l.f[8][l.index(l.nx-1,j)];	// moving SE
			}

			// West moving particles
			for (int j = 0; j <= l.ny-1; j++)
			{
				l.f[3][l.index(l.nx,j)] = l.f[3][l.index(0,j)];	//moving West
				l.f[6][l.index(l.nx,j)] = l.f[6][l.index(0,j)];	// moving NW
				l.f[7][l.index(l.nx,j)] = l.f[7][l.index(0,j)];	// moving SW
			}

			// North moving particles
			for (int i = 0; i <= l.nx; ++i)
			{
				l.f[2][l.index(i,-1)] = l.f[2][l.index(i,l.ny-1)]; // moving north
				l.f[5][l.index(i,-1)] = l.f[5][l.index(i,l.ny-1)]; // moving NE
				l.f[6][l.index(i,-1)] = l.f[6][l.index(i,l.ny-1)]; // moving NW
			}

			// South moving particles
			for (int i = 0; i <= l.nx; ++i)
			{
				l.f[4][l.index(i,0)] = l.f[4][l.index(i,l.ny)]; // moving south
				l.f[7][l.index(i,0)] = l.f[7][l.index(i,l.ny)]; // moving SW
				l.f[8][l.index(i,0)] = l.f[8][l.index(i,l.ny)]; // moving SE
			}

			// Buffer CORNERS// FLOW OR Advection whatever you want to call it
			// // Start at NW corner, calculate North and NW streaming particles
			// for (int i=0; i <= l.nx-1; ++i)
			// {
			// for (int j=l.ny-1; j>=0; --j)
			// {
			// l.f[2][l.index(i,j)] = l.f[2][l.index(i,j-1)];		// N movement
			// l.f[6][l.index(i,j)] = l.f[6][l.index(i+1,j-1)];	//NW movement
			// }
			// }

			// // Start at NE corner, calculate East and NE streaming particles
			// for (int i=l.nx-1; i >= 0; --i)
			// {
			// for (int j=l.ny-1; j>=0; --j)
			// {
			// l.f[1][l.index(i,j)] = l.f[1][l.index(i-1,j)]; 		// East movement
			// l.f[5][l.index(i,j)] = l.f[5][l.index(i-1,j-1)];	// NE movement
			// }
			// }

			// // Start at SE corner, calculate South and SE streaming particles
			// for (int i=l.nx-1; i >= 0; --i)
			// {
			// for (int j=0; j<=l.ny-1; ++j)
			// {
			// l.f[4][l.index(i,j)] = l.f[4][l.index(i,j+1)]; 		// South movement
			// l.f[8][l.index(i,j)] = l.f[8][l.index(i-1,j+1)];	// SE movement
			// }
			// }

			// // Start at SW corner, calculate West and SW streaming particles
			// for (int i=0; i <= l.nx-1; ++i)
			// {
			// for (int j=0; j<=l.ny-1; ++j)
			// {
			// l.f[3][l.index(i,j)] = l.f[3][l.index(i+1,j)]; 		// West movement
			// l.f[7][l.index(i,j)] = l.f[7][l.index(i+1,j+1)];	// SW movement
			// }
			// }
			//corners
			l.f[5][l.index(-1,-1)] = l.f[5][l.index(l.nx-1,l.ny-1)] ; // sw corner
			l.f[6][l.index(l.nx,-1)] = l.f[6][l.index(0,l.ny-1)] ; // se corner
			l.f[7][l.index(l.nx,l.ny)] = l.f[7][l.index(0,0)] ; // ne corner
			l.f[8][l.index(-1,l.ny)] = l.f[8][l.index(l.nx-1,0)] ; // NW corner


			//

			// *** Advect right
			for (int i=l.nx-1; i >= 0; --i)
			{
				for (int j = -1; j <= l.ny; j++)
				{
					l.f[1][l.index(i+1,j)] = l.f[1][l.index(i,j)];	//moving east
					l.f[5][l.index(i+1,j)] = l.f[5][l.index(i,j)];	// moving NE
					l.f[8][l.index(i+1,j)] = l.f[8][l.index(i,j)];	// moving SE
				}
			}
			// *** Advect left
			for (int i=0; i<l.nx-1; ++i)
			{
				for (int j = -1; j <= l.ny; j++)
				{
					l.f[3][l.index(l.nx-1,j)] = l.f[3][l.index(i,j)];	//moving West
					l.f[6][l.index(l.nx-1,j)] = l.f[6][l.index(i,j)];	// moving NW
					l.f[7][l.index(l.nx-1,j)] = l.f[7][l.index(i,j)];	// moving SW
				}
			}
			// *** Advect Up
			for (int j = 0; j <= l.ny-1; j++)
			{
				for (int i = -1; i <= l.nx; ++i)
				{
					l.f[2][l.index(i,j+1)] = l.f[2][l.index(i,j)]; // moving north
					l.f[5][l.index(i,j+1)] = l.f[5][l.index(i,j)]; // moving NE
					l.f[6][l.index(i,j+1)] = l.f[6][l.index(i,j)]; // moving NW
				}

			}
			// *** Advect Down
			for (int j = l.ny-1; j >= 0; j--)
			{
				for (int i = -1; i <= l.nx; ++i)
				{
					l.f[4][l.index(i,j-1)] = l.f[4][l.index(i,j)]; // moving south
					l.f[7][l.index(i,j-1)] = l.f[7][l.index(i,j)]; // moving SW
					l.f[8][l.index(i,j-1)] = l.f[8][l.index(i,j)]; // moving SE
				}

			}
			// **************************

		}

		/**  @brief apply wall boundary conditions */
		void wall_bc()
		{
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
			for(int i= -1; i<=l.nx; ++i)
			{
				for(int j= -1 ; j<=l.ny; ++j)
				{

					// Calculate new density
//					l.rho[l.index(i,j)] = 0;
					l.get_node(i,j).rho() = 0;
					for(int k=0; k<velocity_set().size; ++k)
					{
						l.get_node(i,j).rho() += l.f[k][l.index(i,j)];
					}
					// Calculate new velocities
//					l.get_node(i,j).u() = 0;
//					l.get_node(i,j).v() = 0;
					l.u[l.index(i,j)] = 0;
					l.v[l.index(i,j)] = 0;
					for(int k=0; k<velocity_set().size; ++k)
					{
						l.get_node(i,j).u() += l.get_node(i,j).f(k) * velocity_set().c[0][k];
//						l.u[l.index(i,j)] += l.f[k][l.index(i,j)] * velocity_set().c[0][k];
						l.get_node(i,j).v() += l.get_node(i,j).f(k) * velocity_set().c[1][k];
//						l.v[l.index(i,j)] += l.f[k][l.index(i,j)] * velocity_set().c[1][k];

					}

					l.get_node(i,j).u() /= l.get_node(i,j).rho() ;
					l.get_node(i,j).v() /= l.get_node(i,j).rho() ;
//
//					l.u[l.index(i,j)] /= l.rho[l.index(i,j)];
//					l.v[l.index(i,j)] /= l.rho[l.index(i,j)];

					// Collision step
					lb::float_type feqLocal[9];
					velocity_set().f_eq(feqLocal,l.rho[l.index(i,j)],l.u[l.index(i,j)],l.v[l.index(i,j)]);

					for(int k=0; k<velocity_set().size; ++k)
					{

						l.get_node(i,j).f(k) +=  H_root(l.get_node(i,j))*beta*(feqLocal[k]-l.get_node(i,j).f(k));
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
			std::stringstream fns;
			fns << "output/data_" << std::setfill('0') << std::setw(4) << output_index << ".txt";
			l.write_fields(fns.str());
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

	public: // members

		lattice l;                 ///< lattice
		std::vector<int> shift;    ///< amount of nodes to shift each population in data structure during advection
		const float_type Re;       ///< Reynolds number
		const float_type Vmax;     ///< mean flow velocity
		const float_type visc;     ///< viscosity
		const float_type beta;     ///< LB parameter beta
		unsigned int time;         ///< simulation time
		bool file_output;          ///< flag whether to write files
		unsigned int output_freq;  ///< file output frequency
		unsigned int output_index; ///< index for file naming
	};

} // lb

#endif // LB_SIMULATION_HPP_INCLUDED