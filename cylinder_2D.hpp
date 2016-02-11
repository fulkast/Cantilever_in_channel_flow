#include "geometry_2D.hpp"
#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/Line_arc_2.h>
#include <iterator>

typedef CGAL::Exact_circular_kernel_2 Circle_K;
typedef CGAL::Circle_2<Circle_K> Circle_K_Circle;
typedef CGAL::Line_arc_2<Circle_K> Circle_K_Line_Arc;
typedef CGAL::Point_2<Circle_K> Circle_K_Point;


using namespace std;
class cylinder_2D : public geometry_2D {
    // Define a 2D cylinder derived from the base geometry_2D

public:

    cylinder_2D(lb::coordinate<double> centerOfMass, double orientation, double radius) :
            geometry_2D(centerOfMass,orientation), mRadius(radius) {

        mCircle = new Circle(Point(mCenterOfMass.i,mCenterOfMass.j), mRadius);
        update_shape();
    }

    void update_shape()
    {
        *mCircle = Circle(Point(mCenterOfMass.i,mCenterOfMass.j), mRadius);
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
        mMissingPopulationIndexMap.clear();

        CGAL::Bbox_2 boundingBox = mCircle->bbox();

        // iterate through bounding box nodes with 1 node lee way on all 4 sides
        for (int i = std::floor(CGAL::to_double(boundingBox.xmin())) - 2; i <= std::ceil(CGAL::to_double(boundingBox.xmax())) + 2; i++) {
            for (int j = std::floor(CGAL::to_double(boundingBox.ymin())) - 2; j <= std::ceil(CGAL::to_double(boundingBox.ymax())) + 2; j++) {

                // current point being checked
                Point query_point(i,j);

                // check if point is within the shape or right on the boundary
                if(mCircle->bounded_side(query_point) == CGAL::ON_BOUNDED_SIDE ||
                   mCircle->bounded_side(query_point) == CGAL::ON_BOUNDARY)
                {
                    mInternalNodes.push_back(lb::coordinate<int>(i,j));
                    continue;
                }

                bool is_boundary = false;
                vector<int> outGoingVelocityIndices;

                // check if the shape is on a boundary
                for (int velocity_index = 0; velocity_index < lb::velocity_set().size; velocity_index++)
                {
                    query_point = Point(i+lb::velocity_set().c[0][velocity_index],j+lb::velocity_set().c[1][velocity_index]);

                    if(mCircle->bounded_side(query_point) == CGAL::ON_BOUNDED_SIDE ||
                       mCircle->bounded_side(query_point) == CGAL::ON_BOUNDARY)
                    {
                        // reflect the velocity index to get the corresponding outgoing velocity at the node
                        // and add it to the outgoing velocity indices of the current (boundary) node
                        outGoingVelocityIndices.push_back(lb::velocity_set().incoming_velocity_to_outgoing_velocity(velocity_index));
                        is_boundary = true;
                    }

                }

                // update boundary data
                if (is_boundary)
                {
                    mBoundaryNodes.push_back(lb::coordinate<int>(i,j));
                    mMissingPopulationIndexMap.insert(make_pair(make_pair(i,j),outGoingVelocityIndices));
                }

            }
        }


    }


    double get_ray_length_at_intersection(lb::coordinate<int> boundary_node, int lb_velocity_index)
    {
        double squaredDistance = 0;

        // flip directions to get the ray opposite of the missing population index's ray
        lb_velocity_index = lb::velocity_set().incoming_velocity_to_outgoing_velocity(lb_velocity_index);

        Circle_K_Point boundaryNodePoint(boundary_node.i,boundary_node.j);
        Circle_K_Line_Arc projectingRay(boundaryNodePoint,
                                Circle_K_Point(boundary_node.i+lb::velocity_set().c[0][lb_velocity_index],
                                    boundary_node.j+lb::velocity_set().c[1][lb_velocity_index]));

        std::ostream_iterator<CGAL::Object> out_it (std::cout,", ");

//        out_it = CGAL::intersection(projectingRay,mCircle);





        int i = boundary_node.i - mCenterOfMass.i;
        int j = boundary_node.j - mCenterOfMass.j;

        return sqrt(i*i + j*j) - mRadius;
    }

    std::vector<int> find_missing_populations(lb::coordinate<int> position)
    {
        return mMissingPopulationIndexMap[make_pair(position.i,position.j)];
    }

    lb::coordinate<double> get_velocity_at_intersection(lb::coordinate<int> boundary_node, int lb_velocity_index )
    {
        return lb::coordinate<double>(0.0,0.0);
    }


private:
    double mRadius = 0;
    Circle *mCircle;
};
