#include "global.hpp"
#include <vector>
#include "velocity_set.hpp"
#include <map>
#pragma once

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

    // calculate proximity of boundary node to the shape's true boundary
    virtual double get_shortest_distance_to_true_boundary(lb::coordinate<int> position) = 0;

    // get out-going velocity set indices
    virtual std::vector<int> find_missing_populations(lb::coordinate<int> position) = 0;

    virtual void update_shape();


protected:
    lb::coordinate<double> mCenterOfMass;                             // current object's center of mass
    double mOrientation;                                           // current object's orientation
    std::vector<lb::coordinate<int>> mBoundaryNodes;               // current object's boundaries
    std::vector<lb::coordinate<int>> mInternalNodes;               // current object's wall internal nodes
    std::map<std::pair<int,int>,std::vector<int>>
            mMissingPopulationIndexMap;                             // curre boundary's out going velocity
                                                                    // indices

};


