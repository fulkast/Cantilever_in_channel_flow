#include "geometry_2D.hpp"




geometry_2D::geometry_2D(){}

geometry_2D::geometry_2D(lb::coordinate<double> centerOfMass, double orientation) : mCenterOfMass(centerOfMass.i,centerOfMass.j)
        , mOrientation(orientation)  {}

void geometry_2D::set_center_of_mass(lb::coordinate<double> position)
{
    mCenterOfMass = position;
}

void geometry_2D::set_orientation(double orientation)
{
    mOrientation = orientation;
}

lb::coordinate<double> geometry_2D::get_center_of_mass(){return mCenterOfMass;}

double geometry_2D::get_orientation(){return mOrientation;}

void geometry_2D::print()
{
    std::cout << "centered at: " << mCenterOfMass << " with orientation " << mOrientation;
}

// clears current internal nodes and generates new ones from the current state of the object

std::vector<lb::coordinate<int>> geometry_2D::get_boundary_nodes()
{
    return mBoundaryNodes;
}

std::vector<lb::coordinate<int>> geometry_2D::get_internal_nodes()
{
    return mInternalNodes;
}

void geometry_2D::print_boundary_nodes()
{
    std::cout << "Current boundary nodes are at: " << std::endl;
    for (std::vector<lb::coordinate<int>>::iterator i = mBoundaryNodes.begin(); i != mBoundaryNodes.end(); i++)
    {
        std::cout << *i << std::endl;
    }
}

void geometry_2D::print_out_going_velocities(lb::coordinate<int> position)
{
    if (mMissingPopulationIndexMap.count(std::make_pair(position.i,position.j)) == 0)
    {
        std::cerr << "Wrong boundary situation: No outgoing velocities at queried boundary point." << std::endl;
        return;
    }
    std::cout << "Boundary node at " << position << " has velocities: " << std::endl;
    for (auto i = mMissingPopulationIndexMap[std::make_pair(position.i,position.j)].begin();
         i != mMissingPopulationIndexMap[std::make_pair(position.i,position.j)].end(); i++)
    {
        std::cout << lb::velocity_set().c[0][*i] << " " << lb::velocity_set().c[1][*i] << std::endl;
    }
}

bool geometry_2D::isCurrentlyAnInternalNode(lb::coordinate<int>& tuple)
{
    auto isFoundAt = find(mIsCurrentlyAnInternalNode.begin(),mIsCurrentlyAnInternalNode.end(),std::make_pair(tuple.i,tuple.j));
    return !(isFoundAt == mIsCurrentlyAnInternalNode.end());
}


void geometry_2D::update_shape()
{

}

double geometry_2D::get_angular_velocity()
{
    return mAngularVelocity;
}

void geometry_2D::set_omega(double a_value)
{
    mAngularVelocity = a_value;
}

lb::coordinate<double> geometry_2D::get_linear_velocity()
{
    return mLinearVelocity;
}

void geometry_2D::set_linear_velocity(lb::coordinate<double> a_tuple)
{
    mLinearVelocity = a_tuple;
}

double geometry_2D::get_angular_acceleration()
{
    return mAngularAcceleration;
}

void geometry_2D::set_angular_acceleration(double a_value)
{
    mAngularAcceleration = a_value;
}

lb::coordinate<double> geometry_2D::get_linear_acceleration()
{
    return mLinearAcceleration;
}

void geometry_2D::set_linear_acceleration(lb::coordinate<double> a_tuple)
{
    mLinearAcceleration = a_tuple;
}