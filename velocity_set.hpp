/** 
 *  @file
 *  @author Fabian Bösch
 *  @brief velocity set
 */

#ifndef LB_VELOCITY_SET_HPP_INCLUDED
#define LB_VELOCITY_SET_HPP_INCLUDED

#include "global.hpp"
#include <array>
#include <cmath> 

namespace lb {

struct v9;                // forward declaration
const v9& velocity_set(); // forward declaration
	
/**
 *  @brief Lattice parameters for 9 velocity model.
 *  
 *  This class models a the singleton design pattern. That means there 
 *  exists only one single instance throughout the lifetime of the 
 *  program. To instantiate and access this object use the free function
 *  @ref velocity_set.
 * 
 *  This class holds parameters like lattice weights, molecular 
 *  velocities and speed of sound. It also exposes member functions to 
 *  compute the equilibrium populations.
 */
struct v9 // singleton
{
private:
	
	/** @brief Default constructor */
	v9(){};
	/** @brief Function for instantiating the singleton is a friend */
	friend const v9& lb::velocity_set(); 

public:
	
	v9(const v9&) = delete;
	v9& operator=(const v9&) = delete;

	
	//                                                     0,       1,       2,       3,       4,       5,       6,       7,       8
	const std::array<float_type, 9>         W =   {{ 16.0f/36.0f,  4.0f/36.0f,  4.0f/36.0f,  4.0f/36.0f,  4.0f/36.0f,  1.0f/36.0f,  1.0f/36.0f,  1.0f/36.0f,  1.0f/36.0f}};   ///< Lattice weights
	
	const std::array<std::array<int, 9>, 2> c = {{{{       0,       1,       0,      -1,       0,       1,      -1,      -1,       1}}, 
	                                              {{       0,       0,       1,       0,      -1,       1,       1,      -1,      -1}}}}; ///< Molecular velocities

	const std::array<float_type, 9> magnitude_c = {{0, 1, 1, 1, 1,std::sqrt(2.0f), std::sqrt(2.0f), std::sqrt(2.0f), std::sqrt(2.0f)}};	                                              
	// returns the reflected index of a given velocity index. Used for bounce back calculation against walls
	int incoming_velocity_to_outgoing_velocity(int incoming_velocity_index) const
	{
		switch (incoming_velocity_index)
			{
				case 1:
					return 3;
				case 2:
					return 4;
				case 3:
					return 1;
				case 4:
					return 2;
				case 5:
					return 8;
				case 6:
					return 7;
				case 7:
					return 6;
				case 8:
					return 5;
				default:
				std::cerr << "Invalid velocity switching index called... " << std::endl;
			}
	}
	const float_type cs = 1.0/std::sqrt(3.0);   ///< Speed of sound
	
	const unsigned int size = 9;                ///< Number of velocities

	/** 
	 *  @brief Compute equilibrium.
	 * 
	 *  Compute f_eq from the locally conserved quantities rho, u and v 
	 *  (see also @ref v9::equilibrate).
	 *  @param[in,out] f_eq Pointer to an array of size 9 to store the computed values
	 *  @param[in]     rho  Local density
	 *  @param[in]     u    Local flow velocity in x-direction
	 *  @param[in]     v    Local flow velocity in y-direction
	 */
	inline void f_eq(float_type* f_eq, float_type rho, float_type u, float_type v) const
	{
		// **************************
		// * fill in your code here *
		// **************************
		float a = rho*(2-std::sqrt(1+3*u*u))*(2-std::sqrt(1+3*v*v));
		float ueqx = (2*u+std::sqrt(1+3*u*u))/(1-u);
		float ueqxm = 1/((2*u+std::sqrt(1+3*u*u))/(1-u));
		float ueqy = (2*v + std::sqrt(1+3*v*v))/(1-v);
		float ueqym = 1/((2*v + std::sqrt(1+3*v*v))/(1-v));


		f_eq[0] = W[0]*a;
		f_eq[1] = W[1]*a*ueqx;
		f_eq[2] = W[2]*a*ueqy;
		f_eq[3] = W[3]*a*ueqxm;
		f_eq[4] = W[4]*a*ueqym;
		f_eq[5] = W[5]*a*ueqx*ueqy;
		f_eq[6] = W[6]*a*ueqxm*ueqy;
		f_eq[7] = W[7]*a*ueqxm*ueqym;
		f_eq[8] = W[8]*a*ueqx*ueqym;
	}

	/** 
	 *  @brief Equilibrate a node.
	 * 
	 *  Compute f_eq from the locally conserved quantities rho, u and v
	 *  and set the node's population to that equilibrium ( see also 
	 *  @ref v9::f_eq).
	 *  @tparam        Node A node type
	 *  @param[in,out] n    Reference to a Node object
	 *  @param[in]     rho  Local density
	 *  @param[in]     u    Local flow velocity in x-direction
	 *  @param[in]     v    Local flow velocity in y-direction
	 */
	template <typename Node>
	inline void equilibrate(Node& n, float_type rho, float_type u, float_type v) const
	{
		// **************************
		// * fill in your code here *
		// **************************
		float a = rho*(2-std::sqrt(1+3*u*u))*(2-std::sqrt(1+3*v*v));
		float ueqx = (2*u+std::sqrt(1+3*u*u))/(1-u);
		float ueqxm = 1/((2*u+std::sqrt(1+3*u*u))/(1-u));
		float ueqy = (2*v + std::sqrt(1+3*v*v))/(1-v);
		float ueqym = 1/((2*v + std::sqrt(1+3*v*v))/(1-v));


		n.f(0) = W[0]*a;
		n.f(1)  = W[1]*a*ueqx;
		n.f(2)  = W[2]*a*ueqy;
		n.f(3)   = W[3]*a*ueqxm;
		n.f(4)  = W[4]*a*ueqym;
		n.f(5)  = W[5]*a*ueqx*ueqy;
		n.f(6)  = W[6]*a*ueqxm*ueqy;
		n.f(7)  = W[7]*a*ueqxm*ueqym;
		n.f(8) = W[8]*a*ueqx*ueqym;
	}

	/** 
	 *  @brief Equilibrate a node.
	 * 
	 *  Compute f_eq from the locally conserved quantities rho, u and v
	 *  and set the node's population to that equilibrium ( see also 
	 *  @ref v9::f_eq and v9::equilibrate). The locally conserved 
	 *  quantities are taken form the node object itself.
	 *  @tparam        Node A node type
	 *  @param[in,out] n    Reference to a Node object
	 */
	template <typename Node>
	inline void equilibrate(Node& n) const
	{
		return equilibrate(n, n.rho(), n.u(), n.v());
	}
};

/**
 *  @brief Get a reference single instance of the velocity set.
 *  @return 9-velocity set
 */
inline const v9& velocity_set()
{
	static v9 v_set;
	return v_set;
}

} // lb

#endif // LB_VELOCITY_SET_HPP_INCLUDED
