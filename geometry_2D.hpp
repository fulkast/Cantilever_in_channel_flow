#include "global.hpp"
#include <vector>
#include "velocity_set.hpp"
#include <map>
#include <cassert>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Circle_2.h>
#include <algorithm>

#pragma once

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Aff_transformation_2<K> Transformation;
typedef CGAL::Vector_2<K> Vector;
typedef CGAL::Segment_2<K> Segment;
typedef CGAL::Circle_2<K> Circle;


class geometry_2D {
    // Define a 2D geometry by its center of mass and counter clock wise orientation relative to the
    // x-axis
public:
    // constructors
    geometry_2D();
    geometry_2D(lb::coordinate<double> position, double orientation);

    // setters for the geometry's position and orientation
    void set_center_of_mass(lb::coordinate<double> position);
    void set_orientation(double orientation);

    // getters for the position and orientation of the shape
    lb::coordinate<double> get_center_of_mass();
    double get_orientation();

    // print position orientation and if available, dimensions of the shape
    virtual void print();

    // clears current vector of boundary nodes and refills their entities
    // given the current state of the shape. Also refreshes the internal nodes!
    virtual void update_boundary_and_internal_nodes() = 0;

    // getter for the boundary nodes
    std::vector<lb::coordinate<int>> get_boundary_nodes();

    // getter for the internal nodes
    std::vector<lb::coordinate<int>> get_internal_nodes();

    // print position of current boundary nodes
    void print_boundary_nodes();

    // print out going velocity indices at position
    void print_out_going_velocities(lb::coordinate<int> position);

    // for a given velocity index emanating from a boundary node,
    // find the length at intersection point with the true shape boundary
    virtual double get_ray_length_at_intersection(lb::coordinate<int> boundary_node, int lb_velocity_index) = 0;

    // for a given velocity index emanating from a boundary node,
    // find the velocity at the intersection point with the true shape boundary
    virtual lb::coordinate<double> get_velocity_at_intersection(lb::coordinate<int> boundary_node, int lb_velocity_index ) = 0;

    // get out-going velocity set indices
    virtual std::vector<int> find_missing_populations(lb::coordinate<int> position) = 0;

    virtual void update_shape();

    bool isCurrentlyAnInternalNode(lb::coordinate<int>& tuple);

    double get_angular_velocity();

    void set_omega(double a_value);

    lb::coordinate<double> get_linear_velocity();

    void set_linear_velocity(lb::coordinate<double> a_tuple);

    double get_angular_acceleration();

    void set_angular_acceleration(double a_value);

    lb::coordinate<double> get_linear_acceleration();

    void set_linear_acceleration(lb::coordinate<double> a_tuple);

    void step_forward_one_step();

    virtual void create_sinusoidal_motion(double segment_deflection) = 0;


protected:
    lb::coordinate<double> mCenterOfMass;                           // current object's center of mass
    double mOrientation;                                            // current object's orientation
    std::vector<lb::coordinate<int>> mBoundaryNodes;                // current object's boundaries
    std::vector<lb::coordinate<int>> mInternalNodes;                // current object's wall internal nodes
    std::map<std::pair<int,int>,std::vector<int>>
            mMissingPopulationIndexMap;                             // current boundary's out going velocity
    std::vector<std::pair<int,int>> mIsCurrentlyAnInternalNode;     //
    // indices

    lb::coordinate<double> mLinearVelocity = lb::coordinate<double>(0.0,0.0);
    lb::coordinate<double> mLinearAcceleration = lb::coordinate<double>(0.0,0.0);
    double mAngularVelocity = 0;
    double mAngularAcceleration = 0;

};
