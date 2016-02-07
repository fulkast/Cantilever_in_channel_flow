#include "geometry_2D.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Bbox_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Aff_transformation_2<K> Transformation;
typedef CGAL::Vector_2<K> Vector;

using namespace std;

class quadrilateral_2D : public geometry_2D {
    // Define a 2D cylinder derived from the base geometry_2D

public:

    quadrilateral_2D(lb::coordinate<double> centerOfMass, double orientation, double width, double height) :
            geometry_2D(centerOfMass,orientation), mWidth(width), mHeight(height)
    {
        mPolygon =  new Polygon_2(mPoints, mPoints+4);

        update_shape();
    }

    void update_shape()
    {
        mPoints[0] = Point(-mWidth/2,+mHeight/2);
        mPoints[1] = Point(+mWidth/2,+mHeight/2);
        mPoints[2] = Point(+mWidth/2,-mHeight/2);
        mPoints[3] = Point(-mWidth/2,-mHeight/2);

        mRotate = Transformation(CGAL::ROTATION, sin(mOrientation), cos(mOrientation));
        Transformation translate(CGAL::TRANSLATION, Vector(mCenterOfMass.i, mCenterOfMass.j));
        for (int i = 0; i < 4; i++)
        {
            mPoints[i] = mRotate(mPoints[i]);
            mPoints[i] = translate(mPoints[i]);
        }

        *mPolygon = Polygon_2(mPoints, mPoints+4);
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
        for (int i = std::floor(CGAL::to_double(boundingBox.xmin())) - 2; i <= std::ceil(CGAL::to_double(boundingBox.xmax())) + 2; i++) {
            for (int j = std::floor(CGAL::to_double(boundingBox.ymin())) - 2; j <= std::ceil(CGAL::to_double(boundingBox.ymax())) + 2; j++) {

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

                    if (CGAL::bounded_side_2(mPoints,mPoints+4,query_point, K()) == CGAL::ON_BOUNDED_SIDE ||
                            CGAL::bounded_side_2(mPoints,mPoints+4,query_point, K()) == CGAL::ON_BOUNDARY)
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
        double result = std::numeric_limits<double>::max();
        Point currentPoint(position.i,position.j);
        for (auto anEdge = mPolygon->edges_begin(); anEdge != mPolygon->edges_end(); anEdge++)
        {
            result = std::min(result,
                              CGAL::to_double(CGAL::squared_distance(currentPoint,*anEdge))
            );
        }
        return std::sqrt(result);
    }

    std::vector<int> find_missing_populations(lb::coordinate<int> position)
    {
         return mMissingPopulationIndexMap[make_pair(position.i,position.j)];
    }


private:
    double mWidth = 0;
    double mHeight = 0;
    Polygon_2* mPolygon;
    Point mPoints[4];
    Transformation mRotate;
};



