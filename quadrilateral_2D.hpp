#include "geometry_2D.hpp"

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
        mMissingPopulationIndexMap.clear();
        mIsCurrentlyAnInternalNode.clear();

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
                    mIsCurrentlyAnInternalNode.push_back(make_pair(i,j));
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

        Point boundaryNodePoint(boundary_node.i,boundary_node.j);
        Segment projectingRay(boundaryNodePoint,
                              Point(boundary_node.i+lb::velocity_set().c[0][lb_velocity_index],
                                    boundary_node.j+lb::velocity_set().c[1][lb_velocity_index]));

        for (auto edge = mPolygon->edges_begin(); edge != mPolygon->edges_end(); edge++)
        {

            if (CGAL::do_intersect(projectingRay,*edge))
            {
                CGAL::Object o = CGAL::intersection(projectingRay,*edge);
                if(const Point* op = CGAL::object_cast<Point>(&o))
                {
                    //std::cout << "point type " << std::endl;
                    squaredDistance = CGAL::to_double(CGAL::squared_distance(boundaryNodePoint, *op));
                } else if (const Segment* os = CGAL::object_cast<Segment>(&o))
                {
                    //std::cout << "segment type " << std::endl;
                    squaredDistance =
                            std::min(CGAL::to_double(CGAL::squared_distance(boundaryNodePoint, os->target())),
                                     CGAL::to_double(CGAL::squared_distance(boundaryNodePoint, os->source())));
                } else
                {
                    std::runtime_error(" no intersection at __function__");
                }

            }


        }

        return std::sqrt(squaredDistance);
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
    double mWidth = 0;
    double mHeight = 0;
    Polygon_2* mPolygon;
    Point mPoints[4];
    Transformation mRotate;
};
