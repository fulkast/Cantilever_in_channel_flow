#include "geometry_2D.hpp"


class cylinder_2D : public geometry_2D {
    // Define a 2D cylinder derived from the base geometry_2D

public:

    cylinder_2D(lb::coordinate<double> centerOfMass, double orientation, double radius) :
            geometry_2D(centerOfMass,orientation), mRadius(radius) {
        update_boundary_and_internal_nodes();
    }

    void print()
    {
        std::cout << "A cylinder ";
        geometry_2D::print();
        std::cout << " and radius of: " << mRadius << std::endl;
    }

    // clears current internal nodes and generates new ones from the current state of the object
    void update_boundary_and_internal_nodes()
    {
        mBoundaryNodes.clear();
        mInternalNodes.clear();

        // push in points that are between the current radius and 1 integer radial distance further away exclusive of the boundaries
        // we already entered the 4 cardinal points earlier so every other boundary node should be strictly between
        // the current radius and current radius plus one ## note the integer nature of the radius
        for (int i = std::floor(-mRadius-1); i < std::ceil(mRadius+1); i++) {
            for (int j = std::floor(-mRadius - 1); j < std::ceil(mRadius + 1); j++) {
                if (i*i + j*j <= mRadius*mRadius)
                {
                    // set as wall internal node
                    mInternalNodes.push_back(lb::coordinate<int>(mCenterOfMass.i+i,mCenterOfMass.j+j));
                    continue;
                }
                if (i*i + j*j > (mRadius+1)*(mRadius+1)) continue;

                // set as boundary node
                mBoundaryNodes.push_back(lb::coordinate<int>(mCenterOfMass.i+i,mCenterOfMass.j+j));
            }
        }

    }

    double get_ray_length_at_intersection(lb::coordinate<int> boundary_node, int lb_velocity_index)
    {
        int i = boundary_node.i - mCenterOfMass.i;
        int j = boundary_node.j - mCenterOfMass.j;

        return sqrt(i*i + j*j) - mRadius;
    }

    std::vector<int> find_missing_populations(lb::coordinate<int> position)
    {
        int i = position.i - mCenterOfMass.i;
        int j = position.j - mCenterOfMass.j;

        std::vector<int> out_going_velocity_indices;

        for (int index = 0; index < lb::velocity_set().size; i++)
        {
            // get indices of the velocity set members pointing
            // outwards to the query point "position"
            if((lb::velocity_set().c[0][index] * i) + (lb::velocity_set().c[1][index] * j) > 0)
            {
                out_going_velocity_indices.push_back(index);
            }
        }

        return out_going_velocity_indices;
    }

    virtual lb::coordinate<double> get_velocity_at_intersection(lb::coordinate<int> boundary_node, int lb_velocity_index )
    {


    }


private:
    double mRadius = 0;
};
