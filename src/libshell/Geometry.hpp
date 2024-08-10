//
//  Geometry.hpp
//  Elasticity
//
//  Created by Wim van Rees on 2/24/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef Geometry_h
#define Geometry_h

#include "common.hpp"

#include <functional>

#include <igl/boundary_facets.h>
#include <igl/remove_unreferenced.h>
#include <igl/upsample.h>

#include "GeometryBase.hpp"


class Geometry_Dummy : public PlateGeometry
{
    Eigen::MatrixXd my_vertices;
    Eigen::MatrixXi my_face2vertices;

    void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        vertices = my_vertices;
        face2vertices = my_face2vertices;
        initVertexBoundaries(vertices.rows(), vertices_bc);
    }

public:
    Geometry_Dummy(const Eigen::Ref<const Eigen::MatrixXd> my_vertices_in, const Eigen::Ref<const Eigen::MatrixXi> my_face2vertices_in)
    {
        my_vertices = my_vertices_in;
        my_face2vertices = my_face2vertices_in;
    }
};


/*! \class TwoTriangles
 * \brief Geometry consisting of two equilateral triangles  (L=1) sharing one edge
 *
 */
class TwoTriangles : public PlateGeometry
{
    /* the picture is this

         2
     0   |   3
         1

     triangle 0 -> 0,1,2
     triangle 1 -> 1,3,2
     edge 0 -> v0, v1
     edge 1 -> v0, v2
     edge 2 -> v1, v2
     edge 3 -> v1, v3
     edge 4 -> v2, v3
     */
    void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        vertices.resize(4,3);
        face2vertices.resize(2,3);

        vertices.row(0) << 1 - 0.5*std::sqrt(3),0.5,0;
        vertices.row(1) << 1,0,0;
        vertices.row(2) << 1,1,0;
        vertices.row(3) << 1 + 0.5*std::sqrt(3),0.5,0;

        face2vertices.row(0) << 0,1,2;
        face2vertices.row(1) << 1,3,2;

        // set all free
        initVertexBoundaries(4, vertices_bc);
    }

public:

};


class RectangularPlate : public PlateGeometry
{
protected:
    const Real halfX, halfY, res;
    const std::array<bool,2> fixedEdgeX, fixedEdgeY;

    virtual void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        Eigen::MatrixXd polygon_vertices(4,2);
        Eigen::MatrixXi polygon_edges(4,2);
        polygon_vertices(0,0) = -halfX;
        polygon_vertices(0,1) = -halfY;
        polygon_vertices(1,0) = +halfX;
        polygon_vertices(1,1) = -halfY;
        polygon_vertices(2,0) = +halfX;
        polygon_vertices(2,1) = +halfY;
        polygon_vertices(3,0) = -halfX;
        polygon_vertices(3,1) = +halfY;
        polygon_edges(0,0) = 0;
        polygon_edges(0,1) = 1;
        polygon_edges(1,0) = 1;
        polygon_edges(1,1) = 2;
        polygon_edges(2,0) = 2;
        polygon_edges(2,1) = 3;
        polygon_edges(3,0) = 3;
        polygon_edges(3,1) = 0;

        const Real totalArea = halfX*halfY;
        const std::string flags = "q20a"+helpers::ToStringPrecision(std::sqrt(3.0) * 0.25 * res * res, 20)+(verbose ? "" : "Q");
        Triangulate triangulate(polygon_vertices, polygon_edges, flags);
        triangulate.get(vertices, face2vertices, vertices_bc);

        const int nV = vertices.rows();
        for(int i=0;i<nV;++i)
        {
            // check if we are on the boundary
            if(std::abs(vertices(i,0) + halfX) < 1e-6*halfX and fixedEdgeX[0])
            {
                vertices_bc(i,0) = true;
                vertices_bc(i,1) = true;
                vertices_bc(i,2) = true;
            }

            if(std::abs(vertices(i,0) - halfX) < 1e-6*halfX and fixedEdgeX[1])
            {
                vertices_bc(i,0) = true;
                vertices_bc(i,1) = true;
                vertices_bc(i,2) = true;
            }

            if(std::abs(vertices(i,1) + halfY) < 1e-6*halfY and fixedEdgeY[0])
            {
                vertices_bc(i,0) = true;
                vertices_bc(i,1) = true;
                vertices_bc(i,2) = true;
            }

            if(std::abs(vertices(i,1) - halfY) < 1e-6*halfY and fixedEdgeY[1])
            {
                vertices_bc(i,0) = true;
                vertices_bc(i,1) = true;
                vertices_bc(i,2) = true;
            }
        }
    }

public:
    RectangularPlate(const Real halfX, const Real halfY, const Real res, const std::array<bool,2> fixedEdgeX_, const std::array<bool,2> fixedEdgeY_):
    halfX(halfX),
    halfY(halfY),
    res(res),
    fixedEdgeX(fixedEdgeX_),
    fixedEdgeY(fixedEdgeY_)
    {}

    RectangularPlate(ArgumentParser & parser) :
    halfX(parser.parse<Real>("-halfX", 0.5)),
    halfY(parser.parse<Real>("-halfY", 0.5)),
    res(parser.parse<Real>("-res", 0.1)),
    fixedEdgeX({parser.parse<bool>("-fixedEdgeXi", false), parser.parse<bool>("-fixedEdgeXf", false)}),
    fixedEdgeY({parser.parse<bool>("-fixedEdgeYi", false), parser.parse<bool>("-fixedEdgeYf", false)})
    {}
};


class Cylinder : public ShellGeometry
{
protected:
    const Real radius;
    const Real length;
    const Real res;
    const bool fixBottomEdge, fixTopEdge;

    virtual void getShellGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        const int nx = (int)(2 * M_PI * radius / res) + 1;
        const int ny = (int)(length / res) + 1;
        const Real dx = 2 * M_PI * radius / (nx * 1.0);
        const Real dy = length / (ny * 1.0);

        Eigen::MatrixXd polygon_vertices(2 * nx + 2 * ny, 2);
        Eigen::MatrixXi polygon_edges(2 * nx + 2 * ny, 2);

        // rectangular boundary
        for (int i = 0; i < 2 * nx + 2 * ny; i++) {
            if (i < nx) {
                polygon_vertices(i, 0) = i * dx;
                polygon_vertices(i, 1) = 0.0;
            }
            else if (i < nx + ny) {
                const int j = i - nx;
                polygon_vertices(i, 0) = 2 * M_PI * radius;
                polygon_vertices(i, 1) = j * dy;
            }
            else if (i < 2 * nx + ny) {
                const int j = i - nx - ny;
                polygon_vertices(i, 0) = 2 * M_PI * radius - j * dx;
                polygon_vertices(i, 1) = length;
            }
            else {
                const int j = i - 2 * nx - ny;
                polygon_vertices(i, 0) = 0.0;
                polygon_vertices(i, 1) = length - j * dy;
            }

            polygon_edges(i, 0) = i;
            polygon_edges(i, 1) = (i < 2 * nx + 2 * ny - 1) ? i + 1 : 0;
        }

        // triangulate in plane
        const std::string flags = "q20a"+std::to_string(std::sqrt(3.0)/4.0*res*res)+(verbose ? "" : "Q");

        Eigen::MatrixXd vertices_plane;
        Eigen::MatrixXi face2vertices_plane;
        Eigen::MatrixXb vertices_bc_plane;
        Triangulate triangulate(polygon_vertices, polygon_edges, flags);
        triangulate.get(vertices_plane, face2vertices_plane, vertices_bc_plane);

        const int nFaces = face2vertices_plane.rows();
        const int nVertices = vertices_plane.rows();

        // Get faces belonging to the "patch" (left quarter + right quarter)
        // as well as the "inverse patch"
        std::vector<int> patch_faces_list;
        std::vector<int> patchInverse_faces_list;
        std::vector<int> right_patch_faces_list;

        for (int i = 0; i < nFaces; i++) {
            const Real xc = (vertices_plane(face2vertices_plane(i, 0), 0) + vertices_plane(face2vertices_plane(i, 1), 0) + vertices_plane(face2vertices_plane(i, 2), 0)) / 3.0;
            if (xc < 0.5 * M_PI * radius) {
                patch_faces_list.push_back(i);
            }
            else if (xc > 1.5 * M_PI * radius) {
                patch_faces_list.push_back(i);
                right_patch_faces_list.push_back(i);
            }
            else {
                patchInverse_faces_list.push_back(i);
            }
        }

        // move right patch to left side of mesh
        Eigen::VectorXb moved = Eigen::VectorXb::Constant(vertices_plane.rows(), false);
        for (int i : right_patch_faces_list) {
            for (int j = 0; j < 3; j++) {
                if (moved(face2vertices_plane(i, j))) {
                    continue;
                }
                else {
                    vertices_plane(face2vertices_plane(i, j), 0) -= 2.0 * M_PI * radius;
                    moved(face2vertices_plane(i, j)) = true;
                }
            }
        }

        // remove (unreference) vertices at x = 2pi * radius and identify them
        // with vertices at x = -pi * radius
        for (int i = 0; i < nFaces; i++) {
            for (int j = 0; j < 3; j++) {
                if (face2vertices_plane(i, j) == nx) {
                    face2vertices_plane(i, j) = 0;
                }
                if (nx < face2vertices_plane(i, j) && face2vertices_plane(i, j) <= nx + ny) {
                    face2vertices_plane(i, j) = 2 * nx + 2 * ny - (face2vertices_plane(i, j) - nx);
                }
            }
        }

        // remesh this patch
        Eigen::VectorXi indexList;
        Eigen::MatrixXi patch_edges_;
        Eigen::MatrixXi patch_edges;
        Eigen::MatrixXd patch_vertices_;
        igl::boundary_facets(extractRows(patch_faces_list, face2vertices_plane), patch_edges_);
        igl::remove_unreferenced(vertices_plane, patch_edges_, patch_vertices_, patch_edges, indexList);

        Eigen::MatrixXd patch_vertices;
        Eigen::MatrixXi patch_faces;
        patch_vertices_ = patch_vertices_.block(0, 0, patch_vertices_.rows(), 2);
        Triangulate triangulate_patch(patch_vertices_, patch_edges, flags);
        triangulate_patch.get(patch_vertices, patch_faces, vertices_bc);

        // other (non-patch) part of the mesh
        Eigen::MatrixXd patchInverse_vertices;
        Eigen::MatrixXi patchInverse_faces;
        igl::remove_unreferenced(vertices_plane, extractRows(patchInverse_faces_list, face2vertices_plane), patchInverse_vertices, patchInverse_faces, indexList);

        for (int i = 0; i < patchInverse_vertices.rows(); i++) {
            if (patchInverse_vertices(i, 0) < 0.0) {
                patchInverse_vertices(i, 0) += 2 * M_PI * radius;
            }
        }

        // roll up cylinder
        for (int i = 0; i < patch_vertices.rows(); i++) {
            const Real x = patch_vertices(i, 0);
            const Real y = patch_vertices(i, 1);
            patch_vertices(i, 0) = radius * std::cos(x / radius);
            patch_vertices(i, 1) = radius * std::sin(x / radius);
            patch_vertices(i, 2) = -0.5 * length + y;
        }
        for (int i = 0; i < patchInverse_vertices.rows(); i++) {
            const Real x = patchInverse_vertices(i, 0);
            const Real y = patchInverse_vertices(i, 1);
            patchInverse_vertices(i, 0) = radius * std::cos(x / radius);
            patchInverse_vertices(i, 1) = radius * std::sin(x / radius);
            patchInverse_vertices(i, 2) = -0.5 * length + y;
        }

        // recombine patch with main part of cylinder
        combineMeshes(patch_vertices, patch_faces, patchInverse_vertices, patchInverse_faces, vertices, face2vertices, res * 0.001);

        // boundary conditions
        const Real eps = 0.001 * res;
        vertices_bc.resize(vertices.rows(), vertices.cols());
        for (int i = 0; i < vertices_bc.rows(); i++) {
            if ((fixTopEdge && vertices(i, 2) > length - eps) ||
                (fixBottomEdge && vertices(i, 2) < eps)) {
                for (int j = 0; j < 3; j++) {
                    vertices_bc(i, j) = true;
                }
            } else {
                for (int j = 0; j < 3; j++) {
                    vertices_bc(i, j) = false;
                }
            }
        }
    }

public:
    Cylinder(const Real R, const Real L, const Real res, const bool fixBottomEdge = false, const bool fixTopEdge = false, const bool cyl_verbose = false):
    ShellGeometry(),
    radius(R),
    length(L),
    res(res),
    fixBottomEdge(fixBottomEdge),
    fixTopEdge(fixTopEdge)
    {}

    Cylinder(ArgumentParser & parser):
    radius(parser.parse<Real>("-radius", 1.0)),
    length(parser.parse<Real>("-length", 3.0)),
    res(parser.parse<Real>("-res", 0.1)),
    fixBottomEdge(parser.parse<bool>("-fixBottomEdge", false)),
    fixTopEdge(parser.parse<bool>("-fixTopEdge", false))
    {}

    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d pos) const override
    {
        const Real x = pos(0);
        const Real y = pos(1);
        const Real rad = std::sqrt(x*x + y*y);
        const Real nX = -x/rad; // minus since normal is actually pointed inwards
        const Real nY = -y/rad;

        Eigen::Vector3d n;
        n << nX, nY, 0.0;
        return n.normalized();
    }
};


class ConicalFrustrum : public ShellGeometry
{
protected:
    const Real height;
    const Real bottomRadius;
    const Real topRadius;
    const Real edgeLength;
    const Real fixedTop;
    const Real fixedBottom;

    virtual void getShellGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc ) const override
    {
        // Conformal map from conical frustrum to annulus
        const Real halfConeAngle = std::atan((bottomRadius - topRadius) / height);
        const Real sinHalfConeAngle = std::sin(halfConeAngle);

        const Real innerR = std::pow(topRadius, 1.0 / sinHalfConeAngle);
        const Real outerR = std::pow(bottomRadius, 1.0 / sinHalfConeAngle);
        const int nPointsInner = std::max((int)(2.0 * M_PI * topRadius / edgeLength) + 1, 3);
        const int nPointsOuter = std::max((int)(2.0 * M_PI * bottomRadius / edgeLength) + 1, 3);
        const Real dPhiInner = 2.0 * M_PI / (1.0 * nPointsInner);
        const Real dPhiOuter = 2.0 * M_PI / (1.0 * nPointsOuter);

        Eigen::MatrixXd polygon_vertices(nPointsInner + nPointsOuter, 2);
        Eigen::MatrixXi polygon_edges(nPointsInner + nPointsOuter, 2);
        Eigen::MatrixXd polygon_holes(1, 2);

        for (int i = 0; i < nPointsOuter; i++) {
            polygon_vertices(i, 0) = outerR * std::cos(dPhiOuter * i);
            polygon_vertices(i, 1) = outerR * std::sin(dPhiOuter * i);
            polygon_edges(i, 0) = i;
            polygon_edges(i, 1) = (i != nPointsOuter - 1) ? i + 1 : 0;
        }

        for (int i = 0; i < nPointsInner; i++) {
            int j = nPointsOuter + i;
            polygon_vertices(j, 0) = innerR * std::cos(dPhiInner * i);
            polygon_vertices(j, 1) = innerR * std::sin(dPhiInner * i);
            polygon_edges(j, 0) = j;
            polygon_edges(j, 1) = (i != nPointsInner - 1) ? j + 1 : nPointsOuter;
        }

        polygon_holes(0, 0) = 0.0;
        polygon_holes(0, 1) = 0.0;

        // Triangulation flags
        const std::string flags = "q20";

        // Set triangulation resolution function
        Real args[2];
        args[0] = edgeLength;
        args[1] = sinHalfConeAngle;
        ::setTriangleResolutionFunction([](Real x, Real y, Real* arg) {
            return arg[0] * std::pow(x * x + y * y, 0.5 * (1.0 - arg[1]));
        }, args);

        Triangulate triangulate(polygon_vertices, polygon_edges, polygon_holes, flags);
        triangulate.get(vertices, face2vertices, vertices_bc);

        // Project into 3D
        Real r, p;
        for (int i = 0; i < vertices.rows(); i++) {
            r = std::pow(vertices(i, 0) * vertices(i, 0) + vertices(i, 1) * vertices(i, 1), 0.5 * sinHalfConeAngle);
            p = std::atan2(vertices(i, 1), vertices(i, 0));
            vertices(i, 0) = r * std::cos(p);
            vertices(i, 1) = r * std::sin(p);
            vertices(i, 2) = height * (bottomRadius - r) / (bottomRadius - topRadius);
        }

        // Set boundary conditions
        for (int i = 0; i < vertices.rows(); i++) {
            if (fixedTop && vertices(i, 2) > height - 0.1 * edgeLength) {
                vertices_bc(i, 0) = true;
                vertices_bc(i, 1) = true;
                vertices_bc(i, 2) = true;
            }
            if (fixedBottom && vertices(i, 0) < 0.1 * edgeLength) {
                vertices_bc(i, 0) = true;
                vertices_bc(i, 1) = true;
                vertices_bc(i, 2) = true;
            }
        }
    }

public:
    ConicalFrustrum(const Real height, const Real bottomRadius, const Real topRadius, const Real edgeLength, const bool fixedTop = false, const bool fixedBottom = false):
    height(height),
    bottomRadius(bottomRadius),
    topRadius(topRadius),
    edgeLength(edgeLength),
    fixedTop(fixedTop),
    fixedBottom(fixedBottom)
    {
        if (topRadius >= bottomRadius) {
            std::cout << "ConicalFrustrum geometry: topRadius must be less than bottomRadius!" << std::endl;
            exit(1);
        }
    }

    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d pos) const override
    {
        const Real halfConeAngle = std::atan((bottomRadius - topRadius) / height);
        const Real ct = std::cos(halfConeAngle);
        const Real st = std::sin(halfConeAngle);
        const Real phi = std::atan2(pos(1), pos(0));
        return Eigen::Vector3d(ct * std::cos(phi), ct * std::sin(phi), st);
    }
};


class CircularPlate : public PlateGeometry
{
protected:
    const Real radius;
    const Real res;
    const bool fixedBoundary;

    virtual void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        const int nPointsAlongPerimeter = (int)(2.0 * M_PI * radius / res) + 1;
        const Real dTheta = 2.0 * M_PI / (nPointsAlongPerimeter * 1.0);

        Eigen::MatrixXd polygon_vertices(nPointsAlongPerimeter,2);
        Eigen::MatrixXi polygon_edges(nPointsAlongPerimeter,2);

        polygon_edges(0,0) = 0;
        polygon_edges(nPointsAlongPerimeter-1,1) = 0;

        for(int i=0;i<nPointsAlongPerimeter;++i)
        {
            const Real theta = i*dTheta;
            const Real myX = std::cos(theta)*radius;
            const Real myY = std::sin(theta)*radius;
            polygon_vertices(i,0) = myX;
            polygon_vertices(i,1) = myY;
            if(i<nPointsAlongPerimeter-1)
            {
                polygon_edges(i,1) = i+1;
                polygon_edges(i+1,0) = i+1;
            }
        }

        const std::string flags = "q20a"+std::to_string(std::sqrt(3.0)/4.0*res*res)+(verbose ? "" : "Q");;

        Triangulate triangulate(polygon_vertices, polygon_edges, flags);
        triangulate.get(vertices, face2vertices, vertices_bc);

        if(fixedBoundary)
        {
            const int nV = vertices.rows();
            for(int i=0;i<nV;++i)
            {
                const Real radSq = vertices(i,0)*vertices(i,0) + vertices(i,1)*vertices(i,1);
                const Real diff = std::abs(radSq - radius*radius);
                if(diff < 1e-6*radius)
                {
                    vertices_bc(i,0) = true;
                    vertices_bc(i,1) = true;
                    vertices_bc(i,2) = true;
                }
            }
        }
    }

public:
    CircularPlate(const Real radius, const Real res, const bool fixedBoundary):
    radius(radius),
    res(res),
    fixedBoundary(fixedBoundary)
    {}

    CircularPlate(ArgumentParser & parser):
    radius(parser.parse<Real>("-radius", 1.0)),
    res(parser.parse<Real>("-res", 0.1)),
    fixedBoundary(parser.parse<bool>("-fixedBoundary", false))
    {}
};


class AnnulusPlate : public PlateGeometry
{
protected:
    const Real outerRadius;
    const Real innerRadius;
    const Real edgeLength;
    const bool fixedOuterBoundary;
    const bool fixedInnerBoundary;

    virtual void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        const Real perimeter_inner = 2.0*M_PI*innerRadius;
        const Real perimeter_outer = 2.0*M_PI*outerRadius;

        const int Npoints_inner = (int)(perimeter_inner/edgeLength);
        const int Npoints_outer = (int)(perimeter_outer/edgeLength);

        const Real dTheta_inner = 2.0*M_PI/Npoints_inner;
        const Real dTheta_outer = 2.0*M_PI/Npoints_outer;

        Eigen::MatrixXd polygon_vertices(Npoints_inner + Npoints_outer,2);
        Eigen::MatrixXi polygon_edges(Npoints_inner + Npoints_outer ,2);

        // fill the inner radius
        for(int i=0;i<Npoints_inner;++i)
        {
            const Real theta = i*dTheta_inner;
            const Real myX = std::cos(theta)*innerRadius;
            const Real myY = std::sin(theta)*innerRadius;
            polygon_vertices(i,0) = myX;
            polygon_vertices(i,1) = myY;
        }

        // fill the outer radius
        for(int i=0;i<Npoints_outer;++i)
        {
            const int idx = i + Npoints_inner;
            const Real theta = i*dTheta_outer;
            const Real myX = std::cos(theta)*outerRadius;
            const Real myY = std::sin(theta)*outerRadius;
            polygon_vertices(idx,0) = myX;
            polygon_vertices(idx,1) = myY;
        }

        // edges for inner radius
        for(int i=0;i<Npoints_inner-1;++i)
        {
            polygon_edges(i,0) = i;
            polygon_edges(i,1) = i+1;
        }
        // last edge
        polygon_edges(Npoints_inner-1,0) = Npoints_inner-1;
        polygon_edges(Npoints_inner-1,1) = 0;

        // edges for outer radius : go the other way around (not sure if I need to)
        for(int i=0;i<Npoints_outer;++i)
        {
            const int edge_idx = i + Npoints_inner;
            polygon_edges(edge_idx,0) = edge_idx;
            polygon_edges(edge_idx,1) = edge_idx + 1;
        }
        // last connection
        polygon_edges(Npoints_inner + Npoints_outer - 1, 0) = Npoints_inner + Npoints_outer - 1;
        polygon_edges(Npoints_inner + Npoints_outer - 1, 1) = Npoints_inner;

        // hole
        Eigen::MatrixXd polygon_holes(1,3);
        polygon_holes.row(0) << 0,0,0; // center
        const std::string flags = "q20a"+std::to_string(std::sqrt(3.0)/4.0*edgeLength*edgeLength)+(verbose ? "" : "Q");;

        Triangulate triangulate(polygon_vertices, polygon_edges, polygon_holes, flags);
        triangulate.get(vertices, face2vertices, vertices_bc);

        // boundary conditions
        if(fixedInnerBoundary or fixedOuterBoundary)
        {
            const int nVertices = vertices.rows();
            for(int i=0;i<nVertices;++i)
            {
                const Real rSq = std::pow(vertices(i,0),2) + std::pow(vertices(i,1),2);

                const bool innerR = std::abs(rSq - innerRadius*innerRadius) < 1e-6;
                const bool outerR = std::abs(rSq - outerRadius*outerRadius) < 1e-6;

                if((innerR and fixedInnerBoundary) or (outerR and fixedOuterBoundary))
                {
                    vertices_bc(i,0) = true;
                    vertices_bc(i,1) = true;
                    vertices_bc(i,2) = true;
                }
            }
        }

    }

public:
    AnnulusPlate(const Real outerRadius, const Real innerRadius, const Real edgeLength, const bool fixedOuterBoundary = false, const bool fixedInnerBoundary = false):
    PlateGeometry(),
    outerRadius(outerRadius),
    innerRadius(innerRadius),
    edgeLength(edgeLength),
    fixedOuterBoundary(fixedOuterBoundary),
    fixedInnerBoundary(fixedInnerBoundary)
    {}

    AnnulusPlate(ArgumentParser & parser):
    PlateGeometry(),
    outerRadius(parser.parse<Real>("-outerRadius", 1.0)),
    innerRadius(parser.parse<Real>("-innerRadius", 0.7)),
    edgeLength(parser.parse<Real>("-res", 0.01)),
    fixedOuterBoundary(parser.parse<bool>("-fixedOuterBoundary", false)),
    fixedInnerBoundary(parser.parse<bool>("-fixedInnerBoundary", false))
    {}
};


class SphericalCap : public ShellGeometry
{
protected:
    const Real radius;
    const Real angularRadius;
    int nPointsAroundPerimeter;
    Real res;

    void getShellGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        // Outer edge
        // const int nPointsAroundPerimeter = (int)(2.0 * M_PI * radius * std::sin(angularRadius) / res) + 1;
        const Real dAngle = 2.0 * M_PI / (nPointsAroundPerimeter * 1.0);
        const Real outerRadius = radius * 2.0 * std::sin(angularRadius) / (1.0 + std::cos(angularRadius));

        Eigen::MatrixXd polygon_vertices(nPointsAroundPerimeter, 2);
        Eigen::MatrixXi polygon_edges(nPointsAroundPerimeter, 2);

        for (int i = 0; i < nPointsAroundPerimeter; i++) {
            polygon_vertices(i, 0) = outerRadius * std::cos(i * dAngle);
            polygon_vertices(i, 1) = outerRadius * std::sin(i * dAngle);
            polygon_edges(i, 0) = i;
            polygon_edges(i, 1) = (i < nPointsAroundPerimeter - 1) ? i + 1 : 0;
        }

        // Triangulation flags. Include "YY" option to prevent Steiner points on specified segments.
        const std::string flags = "q20YYa"+std::to_string(0.25 * outerRadius)+(verbose ? "" : "Q");

        // Set triangulation resolution function (from conformal mapping)
        Real args[1];
        args[0] = res;
        ::setTriangleResolutionFunction([](Real u, Real v, Real* arg) {
            Real jacobian_sqrt = 1.0 + 0.25 * (u * u + v * v);
            return arg[0] * jacobian_sqrt;
        }, args);

        // Store projected vertices
        Eigen::MatrixXd projectedVertices;
        Triangulate triangulate(polygon_vertices, polygon_edges, flags);
        triangulate.get(projectedVertices, face2vertices, vertices_bc);

        // Map 2D triangulated surface to points on 3D sphere
        inverseConformalMap(radius, projectedVertices, vertices);
    }

public:
    SphericalCap(const Real radius, const Real angularRadius, const int nPointsAroundPerimeter):
    radius(radius),
    angularRadius(angularRadius),
    nPointsAroundPerimeter(nPointsAroundPerimeter)
    {
        assert(angularRadius > 0.0);
        assert(angularRadius < M_PI);
        res = 2.0 * M_PI * radius * std::sin(angularRadius) / (nPointsAroundPerimeter * 1.0);
    }

    SphericalCap(const Real radius, const Real angularRadius, const Real res):
    radius(radius),
    angularRadius(angularRadius),
    res(res)
    {
        assert(angularRadius > 0.0);
        assert(angularRadius < M_PI);
        nPointsAroundPerimeter = (int)(2.0 * M_PI * radius * std::sin(angularRadius) / res) + 1;
    }

    SphericalCap(ArgumentParser & parser):
    radius(parser.parse<Real>("-radius", 1.0)),
    angularRadius(parser.parse<Real>("-angularRadius", 1.0)),
    res(parser.parse<Real>("-res", 0.01))
    {
        assert(angularRadius > 0.0);
        assert(angularRadius < M_PI);
        nPointsAroundPerimeter = (int)(2.0 * M_PI * radius * std::sin(angularRadius) / res) + 1;
    }

    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d pos) const override
    {
        return pos.normalized();
    }
};


class Sphere : public ShellGeometry
{
protected:
    const Real radius;
    const Real res;

    void getShellGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        const int nPointsAroundPerimeter = (int)(2.0 * M_PI * radius * std::sin(0.25 * M_PI) / res) + 1;

        // Create two spherical caps on either side of the +45˚ latitude line
        Eigen::MatrixXd north_vertices;
        Eigen::MatrixXi north_faces;
        Eigen::MatrixXb north_bc;
        SphericalCap north(radius, 0.25 * M_PI, nPointsAroundPerimeter);
        north.get(north_vertices, north_faces, north_bc);

        Eigen::MatrixXd south_vertices;
        Eigen::MatrixXi south_faces;
        Eigen::MatrixXb south_bc;
        SphericalCap south(radius, 0.75 * M_PI, nPointsAroundPerimeter);
        south.get(south_vertices, south_faces, south_bc);
        south_vertices.col(2) = -south_vertices.col(2);

        // Put them together into one mesh.
        int nnv = north_vertices.rows();
        int nsv = south_vertices.rows();
        int nnf = north_faces.rows();
        int nsf = south_faces.rows();
        north_vertices.conservativeResize(nnv + nsv - nPointsAroundPerimeter, 3);
        north_bc.conservativeResize(nnv + nsv - nPointsAroundPerimeter, 3);
        north_faces.conservativeResize(nnf + nsf, 3);

        for (int i = 0; i < nPointsAroundPerimeter; i++) {
            assert(std::abs(north_vertices(i, 0) - south_vertices(i, 0)) < 1e-6);
            assert(std::abs(north_vertices(i, 1) - south_vertices(i, 1)) < 1e-6);
            assert(std::abs(north_vertices(i, 2) - south_vertices(i, 2)) < 1e-6);
        }

        for (int i = nPointsAroundPerimeter; i < nsv; i++) {
            north_vertices(nnv + i - nPointsAroundPerimeter, 0) = south_vertices(i, 0);
            north_vertices(nnv + i - nPointsAroundPerimeter, 1) = south_vertices(i, 1);
            north_vertices(nnv + i - nPointsAroundPerimeter, 2) = south_vertices(i, 2);
            north_bc(nnv + i, 0) = south_bc(i, 0);
            north_bc(nnv + i, 1) = south_bc(i, 1);
            north_bc(nnv + i, 2) = south_bc(i, 2);
        }

        for (int i = 0; i < nsf; i++) {
            north_faces(nnf + i, 0) = (south_faces(i, 0) >= nPointsAroundPerimeter) ? nnv + south_faces(i, 0) - nPointsAroundPerimeter : south_faces(i, 0);
            north_faces(nnf + i, 1) = (south_faces(i, 2) >= nPointsAroundPerimeter) ? nnv + south_faces(i, 2) - nPointsAroundPerimeter : south_faces(i, 2); // flip triangle
            north_faces(nnf + i, 2) = (south_faces(i, 1) >= nPointsAroundPerimeter) ? nnv + south_faces(i, 1) - nPointsAroundPerimeter : south_faces(i, 1);
        }

        vertices = north_vertices;
        face2vertices = north_faces;
        vertices_bc = north_bc;

        // Find the northern hemisphere.
        std::vector<int> keepFaces;
        std::vector<int> removeFaces;
        for (int i = 0; i < face2vertices.rows(); i++) {
            Eigen::Vector3d faceCenter = (vertices.row(face2vertices(i, 0)) + vertices.row(face2vertices(i, 1)) + vertices.row(face2vertices(i, 2))) / 3.0;
            if (faceCenter(2) > radius * std::cos(0.5 * M_PI)) {
                keepFaces.push_back(i);
            }
            else {
                removeFaces.push_back(i);
            }
        }

        // Extract northern/southern hemisphere patches and the bounding edge vertices (for the northern hemisphere)
        Eigen::VectorXi indexList;
        igl::remove_unreferenced(vertices, extractRows(keepFaces, face2vertices), north_vertices, north_faces, indexList);
        igl::remove_unreferenced(vertices, extractRows(removeFaces, face2vertices), south_vertices, south_faces, indexList);

        // Map the northern hemisphere to the plane
        Eigen::MatrixXd north_vertices_uv;
        conformalMap(radius, north_vertices, north_vertices_uv);

        // Remesh this patch using Triangle. Area constraints remain the same.
        Eigen::MatrixXi north_boundary_edges_;
        igl::boundary_facets(north_faces, north_boundary_edges_);

        Eigen::MatrixXi north_boundary_edges;
        Eigen::MatrixXd north_boundary_vertices;
        igl::remove_unreferenced(north_vertices_uv, north_boundary_edges_, north_boundary_vertices, north_boundary_edges, indexList);

        // Reset resolution function (due to lifetime of args)
        ::clearTriangleResolutionFunction();
        Real args[1];
        args[0] = res;
        ::setTriangleResolutionFunction([](Real u, Real v, Real* arg) {
            Real jacobian_sqrt = 1.0 + 0.25 * (u * u + v * v);
            return arg[0] * jacobian_sqrt;
        }, args);

        const std::string patch_flags = "q20a"+std::to_string(0.25 * radius)+(verbose ? "" : "Q");
        Triangulate triangulate(north_boundary_vertices, north_boundary_edges, patch_flags);
        triangulate.get(north_vertices_uv, north_faces, north_bc);

        inverseConformalMap(radius, north_vertices_uv, north_vertices);

        // Now we need to combine the northern mesh: north_vertices, north_faces
        // and the southern mesh: south_vertices, south_faces (removing excess
        // points from the southern mesh).
        nnv = north_vertices.rows();
        nnf = north_faces.rows();
        nsv = south_vertices.rows();
        nsf = south_faces.rows();

        Eigen::VectorXi south2north(nsv);

        for (int j = 0; j < nsv; j++) {
            south2north(j) = -1;
            for (int i = 0; i < nnv; i++) {
                if ((north_vertices.row(i) - south_vertices.row(j)).norm() < 0.01 * res) {
                    south2north(j) = i;
                    break;
                }
            }
        }

        // Replace all southern hemisphere face2vertex references.
        for (int f = 0; f < nsf; f++) {
            if (south2north(south_faces(f, 0)) >= 0) {
                south_faces(f, 0) = south2north(south_faces(f, 0));
            }
            else {
                south_faces(f, 0) = nnv + south_faces(f, 0);
            }

            if (south2north(south_faces(f, 1)) >= 0) {
                south_faces(f, 1) = south2north(south_faces(f, 1));
            }
            else {
                south_faces(f, 1) = nnv + south_faces(f, 1);
            }

            if (south2north(south_faces(f, 2)) >= 0) {
                south_faces(f, 2) = south2north(south_faces(f, 2));
            }
            else {
                south_faces(f, 2) = nnv + south_faces(f, 2);
            }
        }

        // Combine all southern hemisphere vertices and faces into the northern ones.
        north_vertices.conservativeResize(nnv + nsv, 3);
        north_faces.conservativeResize(nnf + nsf, 3);
        north_vertices.block(nnv, 0, nsv, 3) = south_vertices;
        north_faces.block(nnf, 0, nsf, 3) = south_faces;

        igl::remove_unreferenced(north_vertices, north_faces, vertices, face2vertices, indexList);

        vertices_bc = Eigen::MatrixXb::Constant(vertices.rows(), 3, false);
    }

public:
    Sphere(const Real radius, const Real res):
    radius(radius),
    res(res)
    {}

    Sphere(ArgumentParser & parser):
    radius(parser.parse<Real>("-radius", 1.0)),
    res(parser.parse<Real>("-res", 0.01))
    {}

    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d pos) const override
    {
        return pos.normalized();
    }
};


class Sphere_Icosahedron : public ShellGeometry
{
protected:
    const Real radius;
    const Real res;

    void getShellGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        // icosahedron
        Eigen::MatrixXd polygon_vertices(12, 3);
        const Real cancer = 0.5 * M_PI - std::atan(0.5);
        const Real capricorn = 0.5 * M_PI + std::atan(0.5);
        const Real spacing = 0.1 * 2.0 * M_PI;
        polygon_vertices << 0.0, 0.0, 1.0,
                            std::cos(0 * spacing) * std::sin(cancer), std::sin(0 * spacing) * std::sin(cancer), std::cos(cancer),
                            std::cos(1 * spacing) * std::sin(capricorn), std::sin(1 * spacing) * std::sin(capricorn), std::cos(capricorn),
                            std::cos(2 * spacing) * std::sin(cancer), std::sin(2 * spacing) * std::sin(cancer), std::cos(cancer),
                            std::cos(3 * spacing) * std::sin(capricorn), std::sin(3 * spacing) * std::sin(capricorn), std::cos(capricorn),
                            std::cos(4 * spacing) * std::sin(cancer), std::sin(4 * spacing) * std::sin(cancer), std::cos(cancer),
                            std::cos(5 * spacing) * std::sin(capricorn), std::sin(5 * spacing) * std::sin(capricorn), std::cos(capricorn),
                            std::cos(6 * spacing) * std::sin(cancer), std::sin(6 * spacing) * std::sin(cancer), std::cos(cancer),
                            std::cos(7 * spacing) * std::sin(capricorn), std::sin(7 * spacing) * std::sin(capricorn), std::cos(capricorn),
                            std::cos(8 * spacing) * std::sin(cancer), std::sin(8 * spacing) * std::sin(cancer), std::cos(cancer),
                            std::cos(9 * spacing) * std::sin(capricorn), std::sin(9 * spacing) * std::sin(capricorn), std::cos(capricorn),
                            0.0, 0.0, -1.0;
        Eigen::MatrixXi polygon_faces(20, 2);
        polygon_faces << 1, 3, 0, // top cap
                         3, 5, 0,
                         5, 7, 0,
                         7, 9, 0,
                         9, 1, 0,
                         1, 2, 3, // top teeth
                         3, 4, 5,
                         5, 6, 7,
                         7, 8, 9,
                         9, 10, 1,
                         2, 4, 3, // bottom teeth
                         4, 6, 5,
                         6, 8, 7,
                         8, 10, 9,
                         10, 2, 1,
                         11, 4, 2, // bottom cap
                         11, 6, 4,
                         11, 8, 6,
                         11, 10, 8,
                         11, 2, 10;

        const int nSubdivisions = (int)(std::log2(radius * cancer / res)) + 1;

        vertices = radius * polygon_vertices.rowwise().normalized();

        face2vertices = polygon_faces;
        vertices_bc = Eigen::MatrixXb::Constant(vertices.rows(), 3, false);
    }

public:
    Sphere_Icosahedron(const Real radius, const Real res):
    radius(radius),
    res(res)
    {}

    Sphere_Icosahedron(ArgumentParser & parser):
    radius(parser.parse<Real>("-radius", 1.0)),
    res(parser.parse<Real>("-res", 0.01))
    {}

    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d pos) const override
    {
        return pos.normalized();
    }
};


class Ellipsoid : public ShellGeometry
{
protected:
    const Real rx;
    const Real ry;
    const Real rz;
    const Real res;

    void getShellGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        // icosahedron
        Eigen::MatrixXd polygon_vertices(12, 3);
        const Real cancer = 0.5 * M_PI - std::atan(0.5);
        const Real capricorn = 0.5 * M_PI + std::atan(0.5);
        const Real spacing = 0.1 * 2.0 * M_PI;
        polygon_vertices << 0.0, 0.0, 1.0,
                            std::cos(0 * spacing) * std::sin(cancer), std::sin(0 * spacing) * std::sin(cancer), std::cos(cancer),
                            std::cos(1 * spacing) * std::sin(capricorn), std::sin(1 * spacing) * std::sin(capricorn), std::cos(capricorn),
                            std::cos(2 * spacing) * std::sin(cancer), std::sin(2 * spacing) * std::sin(cancer), std::cos(cancer),
                            std::cos(3 * spacing) * std::sin(capricorn), std::sin(3 * spacing) * std::sin(capricorn), std::cos(capricorn),
                            std::cos(4 * spacing) * std::sin(cancer), std::sin(4 * spacing) * std::sin(cancer), std::cos(cancer),
                            std::cos(5 * spacing) * std::sin(capricorn), std::sin(5 * spacing) * std::sin(capricorn), std::cos(capricorn),
                            std::cos(6 * spacing) * std::sin(cancer), std::sin(6 * spacing) * std::sin(cancer), std::cos(cancer),
                            std::cos(7 * spacing) * std::sin(capricorn), std::sin(7 * spacing) * std::sin(capricorn), std::cos(capricorn),
                            std::cos(8 * spacing) * std::sin(cancer), std::sin(8 * spacing) * std::sin(cancer), std::cos(cancer),
                            std::cos(9 * spacing) * std::sin(capricorn), std::sin(9 * spacing) * std::sin(capricorn), std::cos(capricorn),
                            0.0, 0.0, -1.0;
        Eigen::MatrixXi polygon_faces(20, 3);
        polygon_faces << 1, 3, 0, // top cap
                         3, 5, 0,
                         5, 7, 0,
                         7, 9, 0,
                         9, 1, 0,
                         1, 2, 3, // top teeth
                         3, 4, 5,
                         5, 6, 7,
                         7, 8, 9,
                         9, 10, 1,
                         2, 4, 3, // bottom teeth
                         4, 6, 5,
                         6, 8, 7,
                         8, 10, 9,
                         10, 2, 1,
                         11, 4, 2, // bottom cap
                         11, 6, 4,
                         11, 8, 6,
                         11, 10, 8,
                         11, 2, 10;

        const Real radius = std::pow(rx * ry * rz, 0.33333333333333);
        const int nSubdivisions = (int)(std::log2(radius * cancer / res)) + 1;
        igl::upsample(polygon_vertices, polygon_faces, nSubdivisions);

        vertices = polygon_vertices.rowwise().normalized();

        face2vertices = polygon_faces;
        vertices_bc = Eigen::MatrixXb::Constant(vertices.rows(), 3, false);

        // stretch to ellipsoid
        for (int i = 0; i < vertices.rows(); ++i) {
            vertices(i, 0) = vertices(i, 0) * rx;
            vertices(i, 1) = vertices(i, 1) * ry;
            vertices(i, 2) = vertices(i, 2) * rz;
        }
    }

public:
    Ellipsoid(const Real rx, const Real ry, const Real rz, const Real res):
    rx(rx),
    ry(ry),
    rz(rz),
    res(res)
    {}

    Ellipsoid(ArgumentParser & parser):
    rx(parser.parse<Real>("-rx", 1.0)),
    ry(parser.parse<Real>("-ry", 1.0)),
    rz(parser.parse<Real>("-rz", 1.0)),
    res(parser.parse<Real>("-res", 0.01))
    {}

    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d pos) const override
    {
        Eigen::Vector3d normal(pos(0) / rx / rx, pos(1) / ry / ry, pos(2) / rz / rz);
        return normal.normalized();
    }
};


class MonkeySaddle : public ShellGeometry
{
protected:
    const Real radius;
    const Real amplitude;
    const Real res;

    void getShellGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        // Create a circular disk
        CircularPlate plate(radius, res, false);
        plate.get(vertices, face2vertices, vertices_bc);

        const int nVertices = vertices.rows();
        for (int i = 0; i < nVertices; ++i) {
            vertices(i, 2) = amplitude * (vertices(i, 0) * vertices(i, 0) * vertices(i, 0) - 3.0 * vertices(i, 0) * vertices(i, 1) * vertices(i, 1));
        }
    }

public:
    MonkeySaddle(const Real radius, const Real amplitude, const Real res) :
    radius(radius),
    amplitude(amplitude),
    res(res)
    {}

    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d pos) const override
    {
        Eigen::Vector3d grad(3.0 * amplitude * pos(0) * pos(0) - 3.0 * pos(1) * pos(1),
                             -6.0 * pos(0) * pos(1),
                             -1.0);
        return grad.normalized();
    }
};


class Torus : public ShellGeometry
{
protected:
    const Real a;
    const Real b;
    const Real res;
    Real s;
    Real t;

    virtual void getShellGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        // rectangular boundary
        Eigen::MatrixXd polygon_vertices(4, 2);
        Eigen::MatrixXi polygon_edges(4, 2);
        {
            polygon_vertices(0, 0) = 0.0;
            polygon_vertices(0, 1) = 0.0;
            polygon_vertices(1, 0) = 2 * M_PI * s;
            polygon_vertices(1, 1) = 0.0;
            polygon_vertices(2, 0) = 2 * M_PI * s;
            polygon_vertices(2, 1) = 2 * M_PI;
            polygon_vertices(3, 0) = 0.0;
            polygon_vertices(3, 1) = 2 * M_PI;

            for (int i = 0; i < 4; i++) {
                polygon_edges(i, 0) = i;
                polygon_edges(i, 1) = (i < 3) ? i + 1 : 0;
            }
        }

        // triangulate in plane
        Real args[3];
        args[0] = res;
        args[1] = s;
        args[2] = t;
        ::setTriangleResolutionFunction([](Real u, Real v, Real* arg) {
            Real res = arg[0];
            Real s = arg[1];
            Real t = arg[2];
            Real sqrt_g = t / (std::sqrt(s * s + 1.0) - std::cos(v));
            return res / sqrt_g;
        }, args);

        const std::string flags = "q20";

        Eigen::MatrixXd vertices_plane;
        Eigen::MatrixXi face2vertices_plane;
        Eigen::MatrixXb vertices_bc_plane;
        Triangulate triangulate(polygon_vertices, polygon_edges, flags);
        triangulate.get(vertices_plane, face2vertices_plane, vertices_bc_plane);

        // Move quadrants around
        std::vector<int> NE_faces_list;
        std::vector<int> NW_faces_list;
        std::vector<int> SW_faces_list;
        std::vector<int> SE_faces_list;
        for (int i = 0; i < face2vertices_plane.rows(); i++) {
            const Real xc = (vertices_plane(face2vertices_plane(i, 0), 0) + vertices_plane(face2vertices_plane(i, 1), 0) + vertices_plane(face2vertices_plane(i, 2), 0)) / 3.0;
            const Real yc = (vertices_plane(face2vertices_plane(i, 0), 1) + vertices_plane(face2vertices_plane(i, 1), 1) + vertices_plane(face2vertices_plane(i, 2), 1)) / 3.0;
            if (yc > M_PI) {
                if (xc > M_PI * s) {
                    NE_faces_list.push_back(i);
                } else {
                    NW_faces_list.push_back(i);
                }
            } else {
                if (xc > M_PI * s) {
                    SE_faces_list.push_back(i);
                } else {
                    SW_faces_list.push_back(i);
                }
            }
        }

        Eigen::VectorXi indexList;
        Eigen::MatrixXd NE_vertices;
        Eigen::MatrixXi NE_faces;
        Eigen::MatrixXd NW_vertices;
        Eigen::MatrixXi NW_faces;
        Eigen::MatrixXd SW_vertices;
        Eigen::MatrixXi SW_faces;
        Eigen::MatrixXd SE_vertices;
        Eigen::MatrixXi SE_faces;
        igl::remove_unreferenced(vertices_plane, extractRows(NE_faces_list, face2vertices_plane), NE_vertices, NE_faces, indexList);
        igl::remove_unreferenced(vertices_plane, extractRows(NW_faces_list, face2vertices_plane), NW_vertices, NW_faces, indexList);
        igl::remove_unreferenced(vertices_plane, extractRows(SW_faces_list, face2vertices_plane), SW_vertices, SW_faces, indexList);
        igl::remove_unreferenced(vertices_plane, extractRows(SE_faces_list, face2vertices_plane), SE_vertices, SE_faces, indexList);

        for (int i = 0; i < NE_vertices.rows(); i++) {
            NE_vertices(i, 0) -= 2 * M_PI * s;
            NE_vertices(i, 1) -= 2 * M_PI;
        }

        for (int i = 0; i < NW_vertices.rows(); i++) {
            NW_vertices(i, 1) -= 2 * M_PI;
        }

        for (int i = 0; i < SE_vertices.rows(); i++) {
            SE_vertices(i, 0) -= 2 * M_PI * s;
        }

        Eigen::MatrixXd E_vertices;
        Eigen::MatrixXi E_faces;
        Eigen::MatrixXd W_vertices;
        Eigen::MatrixXi W_faces;
        Eigen::MatrixXd vertices_;
        Eigen::MatrixXi faces_;
        combineMeshes(NE_vertices, NE_faces, SE_vertices, SE_faces, E_vertices, E_faces, 0.001 * t / (std::sqrt(s * s + 1.0) - 1.0));
        combineMeshes(NW_vertices, NW_faces, SW_vertices, SW_faces, W_vertices, W_faces, 0.001 * t / (std::sqrt(s * s + 1.0) - 1.0));
        combineMeshes(E_vertices, E_faces, W_vertices, W_faces, vertices_, faces_, 0.001 * t / (std::sqrt(s * s + 1.0) - 1.0));

        // Remesh
        Eigen::MatrixXi edges_;
        igl::boundary_facets(faces_, edges_);

        Eigen::MatrixXi edges;
        igl::remove_unreferenced(vertices_, edges_, vertices, edges, indexList);

        const std::string flags2 = "q20Q";
        Triangulate triangulate2(vertices, edges, flags2);
        triangulate2.get(vertices, face2vertices, vertices_bc);

        // Roll up the torus
        for (int i = 0; i < vertices.rows(); i++) {
            const Real u = vertices(i, 0);
            const Real v = vertices(i, 1);
            const Real d = std::sqrt(s * s + 1.0) - std::cos(v);
            vertices(i, 0) = s * t * std::cos(u / s) / d;
            vertices(i, 1) = s * t * std::sin(u / s) / d;
            vertices(i, 2) = t * std::sin(v) / d;
        }

        // Sew up the torus
        combineMesh(vertices, face2vertices, 0.001 * res);

        // boundary conditions
        vertices_bc.resize(vertices.rows(), vertices.cols());
        for (int i = 0; i < vertices.rows(); i++) {
            for (int j = 0; j < vertices.cols(); j++) {
                vertices_bc(i, j) = false;
            }
        }
    }

public:
    Torus(const Real a, const Real b, const Real res):
    a(a),
    b(b),
    res(res)
    {
        s = std::sqrt((a / b) * (a / b) - 1.0);
        t = std::sqrt(a * a - b * b);
    }

    Torus(ArgumentParser & parser):
    a(parser.parse<Real>("-a", 1.0)),
    b(parser.parse<Real>("-b", 0.3)),
    res(parser.parse<Real>("-res", 0.01))
    {
        s = std::sqrt((a / b) * (a / b) - 1.0);
        t = std::sqrt(a * a - b * b);
    }

    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d pos) const override
    {
        const Real r_sq = pos(0) * pos(0) + pos(1) * pos(1);
        const Real p = std::sqrt(r_sq) - a;
        Eigen::Vector3d vec(pos(0) * p, pos(1) * p, std::sqrt(r_sq) * pos(2));
        return vec.normalized();
    }
};


class OpenTorus : public ShellGeometry
{
protected:
    const Real a;
    const Real b;
    const Real res;
    Real s;
    Real t;

    virtual void getShellGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        // rectangular boundary
        Eigen::MatrixXd polygon_vertices(4, 2);
        Eigen::MatrixXi polygon_edges(4, 2);
        {
            polygon_vertices(0, 0) = 0.0;
            polygon_vertices(0, 1) = 0.0;
            polygon_vertices(1, 0) = 2 * M_PI * s;
            polygon_vertices(1, 1) = 0.0;
            polygon_vertices(2, 0) = 2 * M_PI * s;
            polygon_vertices(2, 1) = 2 * M_PI;
            polygon_vertices(3, 0) = 0.0;
            polygon_vertices(3, 1) = 2 * M_PI;

            for (int i = 0; i < 4; i++) {
                polygon_edges(i, 0) = i;
                polygon_edges(i, 1) = (i < 3) ? i + 1 : 0;
            }
        }

        // triangulate in plane
        Real args[3];
        args[0] = res;
        args[1] = s;
        args[2] = t;
        ::setTriangleResolutionFunction([](Real u, Real v, Real* arg) {
            Real res = arg[0];
            Real s = arg[1];
            Real t = arg[2];
            Real sqrt_g = t / (std::sqrt(s * s + 1.0) - std::cos(v));
            return res / sqrt_g;
        }, args);

        const std::string flags = "q20";

        // Eigen::MatrixXd vertices_plane;
        // Eigen::MatrixXi face2vertices_plane;
        // Eigen::MatrixXb vertices_bc_plane;
        Triangulate triangulate(polygon_vertices, polygon_edges, flags);
        triangulate.get(vertices, face2vertices, vertices_bc);

        // // Move quadrants around
        // std::vector<int> NE_faces_list;
        // std::vector<int> NW_faces_list;
        // std::vector<int> SW_faces_list;
        // std::vector<int> SE_faces_list;
        // for (int i = 0; i < face2vertices_plane.rows(); i++) {
        //     const Real xc = (vertices_plane(face2vertices_plane(i, 0), 0) + vertices_plane(face2vertices_plane(i, 1), 0) + vertices_plane(face2vertices_plane(i, 2), 0)) / 3.0;
        //     const Real yc = (vertices_plane(face2vertices_plane(i, 0), 1) + vertices_plane(face2vertices_plane(i, 1), 1) + vertices_plane(face2vertices_plane(i, 2), 1)) / 3.0;
        //     if (yc > M_PI) {
        //         if (xc > M_PI * s) {
        //             NE_faces_list.push_back(i);
        //         } else {
        //             NW_faces_list.push_back(i);
        //         }
        //     } else {
        //         if (xc > M_PI * s) {
        //             SE_faces_list.push_back(i);
        //         } else {
        //             SW_faces_list.push_back(i);
        //         }
        //     }
        // }
        //
        // Eigen::VectorXi indexList;
        // Eigen::MatrixXd NE_vertices;
        // Eigen::MatrixXi NE_faces;
        // Eigen::MatrixXd NW_vertices;
        // Eigen::MatrixXi NW_faces;
        // Eigen::MatrixXd SW_vertices;
        // Eigen::MatrixXi SW_faces;
        // Eigen::MatrixXd SE_vertices;
        // Eigen::MatrixXi SE_faces;
        // igl::remove_unreferenced(vertices_plane, extractRows(NE_faces_list, face2vertices_plane), NE_vertices, NE_faces, indexList);
        // igl::remove_unreferenced(vertices_plane, extractRows(NW_faces_list, face2vertices_plane), NW_vertices, NW_faces, indexList);
        // igl::remove_unreferenced(vertices_plane, extractRows(SW_faces_list, face2vertices_plane), SW_vertices, SW_faces, indexList);
        // igl::remove_unreferenced(vertices_plane, extractRows(SE_faces_list, face2vertices_plane), SE_vertices, SE_faces, indexList);
        //
        // for (int i = 0; i < NE_vertices.rows(); i++) {
        //     NE_vertices(i, 0) -= 2 * M_PI * s;
        //     NE_vertices(i, 1) -= 2 * M_PI;
        // }
        //
        // for (int i = 0; i < NW_vertices.rows(); i++) {
        //     NW_vertices(i, 1) -= 2 * M_PI;
        // }
        //
        // for (int i = 0; i < SE_vertices.rows(); i++) {
        //     SE_vertices(i, 0) -= 2 * M_PI * s;
        // }
        //
        // Eigen::MatrixXd E_vertices;
        // Eigen::MatrixXi E_faces;
        // Eigen::MatrixXd W_vertices;
        // Eigen::MatrixXi W_faces;
        // Eigen::MatrixXd vertices_;
        // Eigen::MatrixXi faces_;
        // combineMeshes(NE_vertices, NE_faces, SE_vertices, SE_faces, E_vertices, E_faces, 0.001 * t / (std::sqrt(s * s + 1.0) - 1.0));
        // combineMeshes(NW_vertices, NW_faces, SW_vertices, SW_faces, W_vertices, W_faces, 0.001 * t / (std::sqrt(s * s + 1.0) - 1.0));
        // combineMeshes(E_vertices, E_faces, W_vertices, W_faces, vertices_, faces_, 0.001 * t / (std::sqrt(s * s + 1.0) - 1.0));
        //
        // // Remesh
        // Eigen::MatrixXi edges_;
        // igl::boundary_facets(faces_, edges_);
        //
        // Eigen::MatrixXi edges;
        // igl::remove_unreferenced(vertices_, edges_, vertices, edges, indexList);
        //
        // const std::string flags2 = "q20Q";
        // Triangulate triangulate2(vertices, edges, flags2);
        // triangulate2.get(vertices, face2vertices, vertices_bc);

        // Roll up the torus
        for (int i = 0; i < vertices.rows(); i++) {
            const Real u = vertices(i, 0);
            const Real v = vertices(i, 1);
            const Real d = std::sqrt(s * s + 1.0) - std::cos(v);
            vertices(i, 0) = s * t * std::cos(u / s) / d;
            vertices(i, 1) = s * t * std::sin(u / s) / d;
            vertices(i, 2) = t * std::sin(v) / d;
        }

        // // Sew up the torus
        // combineMesh(vertices, face2vertices, 0.001 * res);

        // // boundary conditions
        // vertices_bc.resize(vertices.rows(), vertices.cols());
        // for (int i = 0; i < vertices.rows(); i++) {
        //     for (int j = 0; j < vertices.cols(); j++) {
        //         vertices_bc(i, j) = false;
        //     }
        // }
    }

public:
    OpenTorus(const Real a, const Real b, const Real res):
    a(a),
    b(b),
    res(res)
    {
        s = std::sqrt((a / b) * (a / b) - 1.0);
        t = std::sqrt(a * a - b * b);
    }

    OpenTorus(ArgumentParser & parser):
    a(parser.parse<Real>("-a", 1.0)),
    b(parser.parse<Real>("-b", 0.3)),
    res(parser.parse<Real>("-res", 0.01))
    {
        s = std::sqrt((a / b) * (a / b) - 1.0);
        t = std::sqrt(a * a - b * b);
    }

    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d pos) const override
    {
        const Real r_sq = pos(0) * pos(0) + pos(1) * pos(1);
        const Real p = std::sqrt(r_sq) - a;
        Eigen::Vector3d vec(pos(0) * p, pos(1) * p, std::sqrt(r_sq) * pos(2));
        return vec.normalized();
    }
};


// For testing
class SphericalShell_Lambert : public ShellGeometry
{
protected:
    const Real sphereRadius;
    const Real height; // should be smaller than radius
    const Real holeHeight; // should be smaller than height
    const Real edgeLength;
    const bool fixedPerimeter;

    void getShellGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        // sphere radius is always one (rescaling later)

        // compute the radius of the disk that will be mapped onto the sphere (should have the same area)
        const Real diskRadius = std::sqrt(2.0 * 1.0 * height / sphereRadius); // surface area of cap is 2 pi R h

        // create disk
        const bool withHole = holeHeight > std::numeric_limits<Real>::epsilon();
        Geometry * geometry = nullptr;
        if(withHole)
        {
            const Real innerRadius = std::sqrt(2.0 * 1.0 * holeHeight / sphereRadius);
            geometry = new AnnulusPlate (diskRadius, innerRadius, edgeLength / sphereRadius, fixedPerimeter, false);
        }
        else
        {
            // compute the points along the perimeter of the final cap
            const Real cap_radius = sphereRadius*std::sqrt(1.0 - std::min(1.0, std::pow((sphereRadius-height)/sphereRadius,2)));
            geometry = new CircularPlate(diskRadius, edgeLength, fixedPerimeter);
        }
        if(not verbose) geometry->setQuiet();

        geometry->get(vertices, face2vertices, vertices_bc);
        delete geometry;

        // apply the projection and rescale
        const int nVertices = vertices.rows();
        for(int i=0;i<nVertices;++i)
        {
            Eigen::Vector3d v_old = vertices.row(i);
            const Real rSq = v_old(0)*v_old(0) + v_old(1)*v_old(1);

            const Real newX = std::sqrt(1 - rSq/4)*v_old(0);
            const Real newY = std::sqrt(1 - rSq/4)*v_old(1);
            const Real newZ = 0.5*rSq - 1;

            vertices(i,0) = newX*sphereRadius;
            vertices(i,1) = newY*sphereRadius;
            vertices(i,2) = -newZ*sphereRadius;
        }
    }

public:
    SphericalShell_Lambert(const Real sphereRadius, const Real height, const Real holeHeight, const Real edgeLength, const bool fixedPerimeter):
    sphereRadius(sphereRadius),
    height(height),
    holeHeight(holeHeight),
    edgeLength(edgeLength),
    fixedPerimeter(fixedPerimeter)
    {
    }

    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d pos) const override
    {
        Eigen::Vector3d n;

        // move the position vector so that z=0 coincides with the sphere origin
        n(0) = pos(0);
        n(1) = pos(1);
        n(2) = pos(2);

        return n.normalized();
    }
};


#endif /* Geometry_h */
