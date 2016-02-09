
#include "simulation.hpp"
#ifdef USE_OPENGL_VISUALIZATION
#include "visualization.hpp"
#endif
#include <omp.h>
#include "global.hpp"
#include <iostream>



int main(int argc, char *argv[])
{
	omp_set_num_threads(std::max(omp_get_max_threads(),omp_get_num_procs()));
	const double D = 10; //Diameter of cylinder in cross flow

	lb::simulation* sim = new lb::simulation(40*D,20*D,20,0.05); // nx, ny, Re, v_max
	sim->initialize();


	// Some stuff to test functionalities
	// sim->l.print_bounding_nodes();
	// sim->l.print_out_going_velocities(lb::coordinate<int>(90,89));
	// std::cout << "The distance to true boundary at " << lb::coordinate<int>(89,100) <<
	// 		" is: " << sim->l.get_shortest_distance_to_true_boundary(lb::coordinate<int>(89,100)) << std::endl;



#ifdef USE_OPENGL_VISUALIZATION
	
		lb::visualization::initialize(sim,argc,argv);
		lb::visualization::get_instance().run();



	#else
	
		// Here are some hints for getting aquainted with the lattice class
		// ================================================================
		
//		// how to print the lattice:
//		// -------------------------
//
//		std::cout << sim->l << std::endl;
//
//		// how to access the lattice:
//		// --------------------------
//
//		// 1) access via node proxy
//		sim->l.get_node(1,0).f(0) = 2;
//
//		// 2) access data directly (make sure you know what you're doing)
//		sim->l.f[0][sim->l.index(2,0)] = 3;
//
//		// 3) using iterators to nodes
//		(sim->l.begin() + sim->l.index(0,0))->f(0) = 1;
//
//
//
//		std::cout << sim->l << std::endl;
		
		
		
		// use a loop like this to run the simulation
		
		for (unsigned int i=0; i<500; ++i)
		{
			std::cout << "Starting step: " << i << std::endl;
			sim->step();
			std::cout << "Made a step: " << i << std::endl;
		}

//	std::cout<<sim->l << std::endl;
	#endif
	
	return 0;
}
