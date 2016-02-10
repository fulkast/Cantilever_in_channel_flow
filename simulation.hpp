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
	  visc(_Vmax*D/_Re),
	  beta( 1/(2*visc/(velocity_set().cs*velocity_set().cs) + 1) ), 
	  time(0),
	  file_output(true), // set to true if you want to write files
	  output_freq(100),
	  output_index(0)

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

		//  initialize nodes [channel flow]
		for (int j=-1; j<=static_cast<int>(l.ny); ++j)
		{
			for (int i=-1; i<=static_cast<int>(l.nx); ++i)
			{
				// double y = j - 0.5;
				// double L_y = l.ny-2;
				l.get_node(i,j).u()   =Vmax ;// 4* Vmax / (L_y*L_y) * (L_y*y - y*y); // poiseulle flow
				l.get_node(i,j).v()   = 0;
				l.get_node(i,j).rho() = 1;
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
		// ADDRESS THIS PROBLEM --> same type? 
		// Store populations for quench
		std::vector<float_type> outlet3(l.ny);
		std::vector<float_type> outlet6(l.ny);
		std::vector<float_type> outlet7(l.ny);
		for (int j=0; j<=l.ny-1; j++)
		{
		 	outlet3[j] = l.f[3][l.index(l.nx-1,j)];
		 	outlet6[j]=  l.f[6][l.index(l.nx-1,j)];
		 	outlet7[j]=  l.f[7][l.index(l.nx-1,j)];
		 }


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
			
			// Periodic boundaries 
			// West wall --> East moving particles 
			for (int j = 0; j <= l.ny-1; j++)
			{
			 l.f[1][l.index(0,j)-shift[1]] = l.f[1][l.index(l.nx,j)];	//moving east
			 l.f[5][l.index(0,j)-shift[5]] = l.f[5][l.index(l.nx,j)];	// moving NE
			 l.f[8][l.index(0,j)-shift[8]] = l.f[8][l.index(l.nx,j)];	// moving SE
			}
			
			// East wall --> West moving particles 
			for (int j = 0; j <= l.ny-1; j++)
			{
			 l.f[3][l.index(l.nx-1,j)-shift[3]] = l.f[3][l.index(-1,j)];	//moving West
			 l.f[6][l.index(l.nx-1,j)-shift[6]] = l.f[6][l.index(-1,j)];	// moving NW
			 l.f[7][l.index(l.nx-1,j)-shift[7]] = l.f[7][l.index(-1,j)];	// moving SW
			}
			
			// South Wall --> North moving particles
			for (int i = 0; i <= l.nx-1; ++i)
			{
			 l.f[2][l.index(i,0)-shift[2]] = l.f[2][l.index(i,l.ny)]; // moving north
			 l.f[5][l.index(i,0)-shift[5]] = l.f[5][l.index(i,l.ny)]; // moving NE
			 l.f[6][l.index(i,0)-shift[6]] = l.f[6][l.index(i,l.ny)]; // moving NW
			}
			
			// North Wall --> South moving particles
			for (int i = 0; i <= l.nx-1; ++i)
			{
			 l.f[4][l.index(i,l.ny-1)-shift[4]] = l.f[4][l.index(i,-1)]; // moving south
			 l.f[7][l.index(i,l.ny-1)-shift[7]] = l.f[7][l.index(i,-1)]; // moving SW
			 l.f[8][l.index(i,l.ny-1)-shift[8]] = l.f[8][l.index(i,-1)]; // moving SE
			}
		
			// Buffer CORNERS
			l.f[5][l.index(0,0)-shift[5]] = l.f[5][l.index(l.nx,l.ny)] ; // sw corner
			l.f[6][l.index(l.nx-1,0)-shift[6]] = l.f[6][l.index(-1,l.ny)] ; // se corner
			l.f[7][l.index(l.nx-1,l.ny-1)-shift[7]] = l.f[7][l.index(-1,-1)] ; // ne corner
			l.f[8][l.index(0,l.ny-1)-shift[8]] = l.f[8][l.index(l.nx,-1)] ; // NW corner

		
			// INLET & OUT BC
			 for (int j = 0; j <= l.ny-1; j++)
			 {
				// double y = j - 0.5;
				// double L_y = l.ny-2;

				double u_inlet = Vmax ; //Uniform flow //4* Vmax / (L_y*L_y) * (L_y*y - y*y);
				double v_inlet = 0;
				double rho_inlet =  1;
								
				//inlet		
					l.get_node(-1,j).u()   = u_inlet; 
					l.get_node(-1,j).v()   = v_inlet;
				 	l.get_node(-1,j).rho() = rho_inlet;
				 	lb::velocity_set().equilibrate(l.get_node(-1,j));

				//outlet (quench)
				// 'refill ghost boundary nodes (right wall) with previous time step' 	
				 	// l.f[3][l.index(l.nx-1,j)] = outlet3;
				 	// l.f[6][l.index(l.nx-1,j)] = outlet6;
				 	// l.f[7][l.index(l.nx-1,j)] = outlet7;

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

		
	}
	
	/**  @brief apply wall boundary conditions */
	void wall_bc()
	{

		for (auto i = l.shapes.begin(); i != l.shapes.end(); i++)
		{

			std::vector<lb::coordinate<int>> currentBoundaryNodes = (*i)->get_boundary_nodes();
			std::vector<lb::coordinate<int>> currentSolidNodes = (*i)->get_internal_nodes();
			(*i)->set_orientation((*i)->get_orientation()+0.1);
			(*i)->update_shape();

			for(auto j = currentSolidNodes.begin(); j != currentSolidNodes.end(); j++)
			{
				l.set_is_wall_node(lb::coordinate<int>(j->i,j->j));
			}

			for(auto j = currentBoundaryNodes.begin(); j != currentBoundaryNodes.end(); j++)
			{
				//
				l.get_node(j->i,j->j).rho() = 100;
			}


		}

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
};

} // lb

#endif // LB_SIMULATION_HPP_INCLUDED
