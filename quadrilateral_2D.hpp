#include "geometry_2D.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Bbox_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;

using namespace std;

class quadrilateral_2D : public geometry_2D {
    // Define a 2D cylinder derived from the base geometry_2D

public:

    quadrilateral_2D(lb::coordinate<int> centerOfMass, double orientation, double width, double height) :
            geometry_2D(centerOfMass,orientation), mWidth(width), mHeight(height)
    {
        mPoints[0] = Point(centerOfMass.i-width/2,centerOfMass.j+height/2);
        mPoints[1] = Point(centerOfMass.i+width/2,centerOfMass.j+height/2);
        mPoints[2] = Point(centerOfMass.i+width/2,centerOfMass.j-height/2);
        mPoints[3] = Point(centerOfMass.i-width/2,centerOfMass.j-height/2);
        mPolygon =  new Polygon_2(mPoints, mPoints+4);
        update_boundary_and_internal_nodes();

    }

    void print()
    {
        std::cout << "A quadrilateral ";
        geometry_2D::print();
        std::cout << " and width: " << mWidth << " and height: " << mHeight << std::endl;
    }

    // clears current internal nodes and generates new ones from the current state of the object
    void update_boundary_and_internal_nodes()
    {
        mBoundaryNodes.clear();
        mInternalNodes.clear();
        mOutGoingVelocityIndexMap.clear();

        // get bounding box of current shape
        CGAL::Bbox_2 boundingBox = mPolygon->bbox();

        // iterate through bounding box nodes with 1 node lee way on all 4 sides
        for (int i = std::floor(CGAL::to_double(boundingBox.xmin())) - 1; i <= std::ceil(CGAL::to_double(boundingBox.xmax())) + 1; i++) {
            for (int j = std::floor(CGAL::to_double(boundingBox.ymin())) - 1; j <= std::ceil(CGAL::to_double(boundingBox.ymax())) + 1; j++) {

                // current point being checked
                Point query_point(i,j);

                // check if point is within the shape or right on the boundary
                if(CGAL::bounded_side_2(mPoints,mPoints+4,query_point, K()) == CGAL::ON_BOUNDED_SIDE ||
                        CGAL::bounded_side_2(mPoints,mPoints+4,query_point, K()) == CGAL::ON_BOUNDARY)
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

                    if (CGAL::bounded_side_2(mPoints,mPoints+4,query_point, K()) == CGAL::ON_BOUNDARY)
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
                    mOutGoingVelocityIndexMap.insert(make_pair(make_pair(i,j),outGoingVelocityIndices));
                }

            }
        }

    }

    double get_shortest_distance_to_true_boundary(lb::coordinate<int> position)
    {
       /* int i = position.i - mCenterOfMass.i;
        int j = position.j - mCenterOfMass.j;

        return sqrt(i*i + j*j) - mRadius;
        */
    }

    std::vector<int> get_out_going_velocity_indices(lb::coordinate<int> position)
    {
         return mOutGoingVelocityIndexMap[make_pair(position.i,position.j)];
    }


private:
    double mWidth = 0;
    double mHeight = 0;
    Polygon_2* mPolygon;
    Point mPoints[4];
};

// CGAL trial
