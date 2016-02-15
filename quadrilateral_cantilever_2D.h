//
// Created by frank on 15.02.16.
//

#ifndef LB2D_FRAMEWORK_QUADRILATERAL_CANTILEVER_2D_H
#define LB2D_FRAMEWORK_QUADRILATERAL_CANTILEVER_2D_H

#include "geometry_2D.hpp"
#include <map>


class quadrilateral_cantilever_2D : public geometry_2D {

public:

//    quadrilateral_cantilever_2D(lb::coordinate<double> centerOfMass, double orientation, double width, double height) :
    quadrilateral_cantilever_2D(lb::coordinate<double> leftMostSegmentCenterOfMass, double orientation, double width, double height, int nSegments, double springWidth) :
     mWidth(width), mHeight(height), mNSegments(nSegments), mSegmentWidth(width/nSegments), mSpringWidth(springWidth), mLeftMostSegmentOrientation(orientation)
    {
        mLeftMostSegment = new quadrilateral_2D(leftMostSegmentCenterOfMass,mLeftMostSegmentOrientation,width/nSegments,height);
        mQuadrilateralSegments.push_back(mLeftMostSegment);

        mUnionPolygon = new Polygon_2(mPoints,mPoints+4);

        for (int i = 1; i < mNSegments;i++)
        {
            lb::coordinate<double> previousSegmentCenterOfMass = mQuadrilateralSegments[i-1]->get_center_of_mass();
            double previousSegmentOrientation = mQuadrilateralSegments[i-1]->get_orientation();
            double segmentCenterX = cos(previousSegmentOrientation)*(mSpringWidth+mSegmentWidth);
            double segmentCenterY = sin(previousSegmentOrientation)*(mSpringWidth+mSegmentWidth);

            double currentSegmentOrientation = 0.05;

            segmentCenterX = cos(currentSegmentOrientation)*segmentCenterX + -sin(currentSegmentOrientation)*segmentCenterY;
            segmentCenterY = sin(currentSegmentOrientation)*segmentCenterX + cos(currentSegmentOrientation)*segmentCenterY;

            lb::coordinate<double> segmentCenter(previousSegmentCenterOfMass.i + segmentCenterX,
                                                  (previousSegmentCenterOfMass.j + segmentCenterY)
            );

            mQuadrilateralSegments.push_back(new quadrilateral_2D(segmentCenter,previousSegmentOrientation+currentSegmentOrientation,mSegmentWidth,mHeight));
        }

        update_shape();
    }

    void update_shape() {

        for (int i = 0; i < mNSegments; i++)
        {
            mQuadrilateralSegments[i]->update_shape();

        }
        update_boundary_and_internal_nodes();



    }

    void print()
    {
        std::cout << "A cantilever ";
        geometry_2D::print();
        std::cout << " and segment width: " << mSegmentWidth << " and height: " << mHeight << std::endl;
    }

    // clears current internal nodes and generates new ones from the current state of the object
    void update_boundary_and_internal_nodes()
    {
        mUnionPoints.clear();

        for (int i = 0; i < mNSegments; i++)
        {
            mUnionPoints.push_back(mQuadrilateralSegments[i]->get_top_left_point());
            mUnionPoints.push_back(mQuadrilateralSegments[i]->get_top_right_point());
//            std::cout << mQuadrilateralSegments[i]->get_top_left_point()
//            << " " << mQuadrilateralSegments[i]->get_top_right_point() << std::endl;
        }

        for (int i = mNSegments-1; i >= 0 ; i--)
        {
            mUnionPoints.push_back(mQuadrilateralSegments[i]->get_bottom_right_point());
            mUnionPoints.push_back(mQuadrilateralSegments[i]->get_bottom_left_point());

//            std::cout << mQuadrilateralSegments[i]->get_bottom_right_point()
//            << " " << mQuadrilateralSegments[i]->get_bottom_left_point() << std::endl;
        }



        *mUnionPolygon = Polygon_2(mUnionPoints.begin(),mUnionPoints.end());

        for (int i = 0; i < mNSegments; i++)
        {
            mQuadrilateralSegments[i]->update_boundary_and_internal_nodes();
        }

        mBoundaryNodes.clear();
        mInternalNodes.clear();
        mMissingPopulationIndexMap.clear();

        // get bounding box of current shape
        CGAL::Bbox_2 boundingBox = mUnionPolygon->bbox();

//        std::cout << boundingBox.xmin() << " " << boundingBox.xmax() << std::endl;

        // iterate through bounding box nodes with 1 node lee way on all 4 sides
        for (int i = std::floor(CGAL::to_double(boundingBox.xmin())) - 2; i <= std::ceil(CGAL::to_double(boundingBox.xmax())) + 2; i++) {
            for (int j = std::floor(CGAL::to_double(boundingBox.ymin())) - 2; j <= std::ceil(CGAL::to_double(boundingBox.ymax())) + 2; j++) {

                // current point being checked
                Point query_point(i,j);

                // check if point is within the shape or right on the boundary
                if(CGAL::bounded_side_2(mUnionPoints.begin(),mUnionPoints.end(),query_point, K()) == CGAL::ON_BOUNDED_SIDE ||
                   CGAL::bounded_side_2(mUnionPoints.begin(),mUnionPoints.end(),query_point, K()) == CGAL::ON_BOUNDARY)
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

                    if (CGAL::bounded_side_2(mUnionPoints.begin(),mUnionPoints.end(),query_point, K()) == CGAL::ON_BOUNDED_SIDE ||
                        CGAL::bounded_side_2(mUnionPoints.begin(),mUnionPoints.end(),query_point, K()) == CGAL::ON_BOUNDARY)
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

        for (auto edge = mUnionPolygon->edges_begin(); edge != mUnionPolygon->edges_end(); edge++)
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

    void update_external_moments()
    {
        for (auto i = mQuadrilateralSegments.begin()+1; i != mQuadrilateralSegments.end()-1; i++)
        {
            double currentSegmentOrientation = (*i)->get_orientation();
            double previousSegmentOrientation = (*i-1)->get_orientation();
            double nextSegmentOrientation = (*i+1)->get_orientation();
            (*i)->set_angular_acceleration(-9.8*cos((*i)->get_orientation())
                                           -1*(currentSegmentOrientation-previousSegmentOrientation)
                                           +1*(nextSegmentOrientation-currentSegmentOrientation));

        }
        auto i = mQuadrilateralSegments.end()-1;
        double currentSegmentOrientation = (*i)->get_orientation();
        double previousSegmentOrientation = (*i-1)->get_orientation();
        (*i)->set_angular_acceleration(-9.8*cos((*i)->get_orientation())
                                       -1*(currentSegmentOrientation-previousSegmentOrientation));

    }



    std::vector<int> find_missing_populations(lb::coordinate<int> position)
    {
        return mMissingPopulationIndexMap[make_pair(position.i,position.j)];
    }

    lb::coordinate<double> get_velocity_at_intersection(lb::coordinate<int> boundary_node, int lb_velocity_index )
    {

        return lb::coordinate<double>(0.0,0.0);

    }


    void create_sinusoidal_motion(double segment_deflection)
    {

        for (auto i = 1; i < mQuadrilateralSegments.size();i ++)
        {
            mQuadrilateralSegments[i]->set_orientation(segment_deflection);
            lb::coordinate<double> previousSegmentCenterOfMass = mQuadrilateralSegments[i-1]->get_center_of_mass();
            double previousSegmentOrientation = mQuadrilateralSegments[i-1]->get_orientation();
            double segmentCenterX = cos(previousSegmentOrientation)*(mSpringWidth+mSegmentWidth);
            double segmentCenterY = sin(previousSegmentOrientation)*(mSpringWidth+mSegmentWidth);
            double currentSegmentOrientation = segment_deflection;

            segmentCenterX = cos(currentSegmentOrientation)*segmentCenterX + -sin(currentSegmentOrientation)*segmentCenterY;
            segmentCenterY = sin(currentSegmentOrientation)*segmentCenterX + cos(currentSegmentOrientation)*segmentCenterY;

            mQuadrilateralSegments[i]->set_center_of_mass(lb::coordinate<double>(previousSegmentCenterOfMass.i + segmentCenterX,
                                                 (previousSegmentCenterOfMass.j + segmentCenterY)
            ));

            mQuadrilateralSegments[i]->set_orientation(segment_deflection+previousSegmentOrientation);

        }

        update_shape();

    }

private:
    double mWidth = 0;
    double mHeight = 0;
    int mNSegments = 0;
    double mSegmentWidth = 0;
    double mSpringWidth;


    Polygon_2* mUnionPolygon;
    std::vector<Point> mUnionPoints;
    quadrilateral_2D* mLeftMostSegment;
    double mLeftMostSegmentOrientation;
    std::vector<quadrilateral_2D*> mQuadrilateralSegments;
    Point mPoints[4];
    Transformation mRotate;


};


#endif //LB2D_FRAMEWORK_CANTILEVER_2D_H
