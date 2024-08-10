//
//  GeometryBase.hpp
//  Elasticity
//
//  Created by Wim van Rees on 2/24/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef Geometry_Base_h
#define Geometry_Base_h

#include <functional>

#include "common.hpp"
#include "ArgumentParser.hpp"
#include "VTKData.hpp"

#ifdef USETRILIBRARY
#include "igl/triangle/triangulate.h"
#include "Triangulate_TRIANGLE.hpp"
#endif

#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/edge_topology.h>


/*! \class Geometry
 * \brief Base class for different basic triangle mesh geometries.
 *
 * @details
 * This class provides one (pure virtual) method that fills the lists of vertices and faces. Every geometry that derives from here represents a different triangle geometry. Boundary conditions for now are not yet taken into account. Specifics for each geometry are initialized in the constructor.
 */
class Geometry
{
protected:
    /**
     * Whether or not to print stuff during running
     */
    mutable bool verbose;

    /**
     * A helper method for initializing the vertices_bc array to false
     */
    void initVertexBoundaries(const int nVertices, Eigen::MatrixXb & vertices_bc) const
    {
        vertices_bc.resize(nVertices,3);
        for(int i=0;i<nVertices;++i) vertices_bc.row(i) << false,false,false;
    }

public:
    Geometry():
    verbose(true)
    {}

    /**
     * The main method of this class: it fills/overwrites the input arrays vertices (Reals, nVertices x 3 entries -- vertex locations in 3D space); face2vertices(integers, nFaces x 3 entries -- each row contains the indices of the 3 vertices that make up that face .. ordering has to be consistent with the face normal);  and vertices_bc(booleans, nVertices x 3 entries -- each row says whether the vertex of that row is fixed or not in each of the 3 cartesian directions)
     */
    virtual void get(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const = 0;

    /**
     * whether or not we can provide an analytic expression of the normal over the entire mid-surface
     */
    virtual bool analyticNormals() const = 0;

    /**
     * returns analytic expression of the normal vector (in world-coordinates) at location pos on the mesh - only used in case analyticNormals returns true
     */
    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d pos) const = 0;

    virtual ~Geometry()
    {}

    void setVerbose() const{verbose=true;}
    void setQuiet() const{verbose=false;}
};


/*! \class IOGeometry
 * \brief Geometry class that defines the topology and vertex locations from an input file
 *
 * @details
 * This class reads an external file from disk that contains the vertex location in three dimensions, as well as the vertx indices for each face. The external file can be an OBJ or OFF file (read using the libIGL methods), or a VTP file (read using the VTK library). Whichever method is used to read the file is chosen based on the file extension. The analyticNormals method returns false : this way the initial edge directors are aligned with the dihedral angles.
 */
class IOGeometry : public Geometry
{
protected:
    enum class MESHTYPE {OBJ, OFF, VTP};

    const std::string filename;
    MESHTYPE meshtype;

    std::string getFileExtension(const std::string& filename) const
    {
        if(filename.find_last_of(".") != std::string::npos)
            return filename.substr(filename.find_last_of(".")+1);
        return "";
    }

public:
    IOGeometry(const std::string filename):
    filename(filename)
    {
        std::string extension = getFileExtension(filename);
        std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower); // make everything lower case
        if(extension == "obj")
            meshtype = MESHTYPE::OBJ;
        else if(extension == "off")
            meshtype = MESHTYPE::OFF;
        else if(extension == "vtp")
            meshtype = MESHTYPE::VTP;
        else
        {
            std::cout << "invalid extension ( " << extension << " ) : no mesh reader available" << std::endl;
        }
    }

    virtual void get(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        // read the vertices and faces into the input arrays
        if (meshtype == MESHTYPE::OBJ) {
            igl::readOBJ(filename, vertices, face2vertices);
        }
        else if (meshtype == MESHTYPE::OFF) {
            igl::readOFF(filename, vertices, face2vertices);
        }
        else if (meshtype == MESHTYPE::VTP) {
            VTKData data(filename);
            vertices = data.getVertices();
            face2vertices = data.getFaces();
        }
        else {
            std::cout << "Invalid meshtype" << std::endl;
        }

        // set the boundary conditions
        initVertexBoundaries(vertices.rows(), vertices_bc);
    }

    /**
     * we have no analytic expression for the normals in this discrete mesh
     */
    bool analyticNormals() const override
    {
        return false;
    }

    /**
     * getNormal is not used - returns zero vector everywhere
     */
    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d ) const override
    {
        return (Eigen::Vector3d() << 0,0,0).finished();
    }
};


/*! \class PlateGeometry
 * \brief Specialization of Geometry class for geometries that are flat and aligned with the xy-plane, so that the normal vector of each face points in the positive z direction
 */
class PlateGeometry : public Geometry
{
    virtual void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const = 0;

public:
    virtual void get(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        getPlateGeometry(vertices, face2vertices, vertices_bc);
    }

    /**
     * for a plate, the normals are oriented vertically
     */
    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d ) const override
    {
        return (Eigen::Vector3d() << 0,0,1).finished();
    }

    bool analyticNormals() const override
    {
        return true;
    }
};


class Triangulate : public PlateGeometry
{
protected:
    const Eigen::MatrixXd & polygon_vertices; // 2D vertex position array
    const Eigen::MatrixXi & polygon_boundary_edges; // list of (2) vertex ids forming unoriented edges
    const Eigen::MatrixXd polygon_holes;
    const Eigen::MatrixXd polygon_areas;
    const std::string & option_flags;
    std::function<Real(Real, Real)> resFunction;

#ifdef USETRILIBRARY
    virtual void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        Eigen::MatrixXd vertices_2D;
        TRIANGLE::triangulate(polygon_vertices, polygon_boundary_edges, polygon_holes, polygon_areas, option_flags, vertices_2D, face2vertices);
        const int nVertices = vertices_2D.rows();
        vertices.resize(nVertices,3);
        for(int i=0;i<nVertices;++i)
        {
            vertices(i,0) = vertices_2D(i,0);
            vertices(i,1) = vertices_2D(i,1);
            vertices(i,2) = 0.0;
        }
        initVertexBoundaries(vertices.rows(), vertices_bc);
    }
#else
    virtual void getPlateGeometry(Eigen::MatrixXd & , Eigen::MatrixXi & , Eigen::MatrixXb & ) const override
    {
        std::cout << "Need USETRILIBRARY compilation flag to use the Triangulate geometry class" << std::endl;
        std::exit(1);
    }
#endif

public:

    Triangulate(const Eigen::MatrixXd & p_vertices, const Eigen::MatrixXi & p_boundaries, const std::string & flags):
    PlateGeometry(),
    polygon_vertices(p_vertices),
    polygon_boundary_edges(p_boundaries),
    option_flags(flags)
    {
    }

    Triangulate(const Eigen::MatrixXd & p_vertices, const Eigen::MatrixXi & p_boundaries, const Eigen::MatrixXd p_holes, const std::string & flags):
    PlateGeometry(),
    polygon_vertices(p_vertices),
    polygon_boundary_edges(p_boundaries),
    polygon_holes(p_holes),
    option_flags(flags)
    {
    }

    Triangulate(const Eigen::MatrixXd & p_vertices, const Eigen::MatrixXi & p_boundaries, const Eigen::MatrixXd p_holes, const Eigen::MatrixXd p_areas, const std::string & flags):
    PlateGeometry(),
    polygon_vertices(p_vertices),
    polygon_boundary_edges(p_boundaries),
    polygon_holes(p_holes),
    polygon_areas(p_areas),
    option_flags(flags)
    {
    }

    Triangulate(const Eigen::MatrixXd & p_vertices, const Eigen::MatrixXi & p_boundaries, const std::string & flags, std::function<Real(Real, Real)> resFunction):
    PlateGeometry(),
    polygon_vertices(p_vertices),
    polygon_boundary_edges(p_boundaries),
    option_flags(flags),
    resFunction(resFunction)
    {
    }

    Triangulate(const Eigen::MatrixXd & p_vertices, const Eigen::MatrixXi & p_boundaries, const Eigen::MatrixXd p_holes, const std::string & flags, std::function<Real(Real, Real)> resFunction):
    PlateGeometry(),
    polygon_vertices(p_vertices),
    polygon_boundary_edges(p_boundaries),
    polygon_holes(p_holes),
    option_flags(flags),
    resFunction(resFunction)
    {
    }

    Triangulate(const Eigen::MatrixXd & p_vertices, const Eigen::MatrixXi & p_boundaries, const Eigen::MatrixXd p_holes, const Eigen::MatrixXd p_areas, const std::string & flags, std::function<Real(Real, Real)> resFunction):
    PlateGeometry(),
    polygon_vertices(p_vertices),
    polygon_boundary_edges(p_boundaries),
    polygon_holes(p_holes),
    polygon_areas(p_areas),
    option_flags(flags),
    resFunction(resFunction)
    {
    }
};


/*! \class ShellGeometry
 * \brief Specialization of Geometry class for geometries that are curved. Sets analyticNormals true so that getNormal needs to be implemented - also provides a helper method
 */
class ShellGeometry : public Geometry
{
protected:
    virtual void getShellGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const  = 0;

    int orientNormalsOutwards(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices) const
    {
        const int nFaces = face2vertices.rows();
        int nFlippedNormals = 0;

        for(int i=0;i<nFaces;++i)
        {
            const int vIdx_0 = face2vertices(i,0);
            const int vIdx_1 = face2vertices(i,1);
            const int vIdx_2 = face2vertices(i,2);

            const Eigen::Vector3d v0 = vertices.row(vIdx_0);
            const Eigen::Vector3d v1 = vertices.row(vIdx_1);
            const Eigen::Vector3d v2 = vertices.row(vIdx_2);

            // get target normal
            const Eigen::Vector3d facepos = (v0 + v1 + v2)/3.0;
            const Eigen::Vector3d target_normal = getNormal(facepos);

            // compute actual normal
            const Eigen::Vector3d e1 = v1 - v0;
            const Eigen::Vector3d e2 = v2 - v0;
            const Eigen::Vector3d normal = e1.cross(e2);//.normalized();
            const Real r = normal.norm();
            if(r < 1e-9)
            {
                std::cout << "DEGENERATE FACE DETECTED : " << r << "\t" << i << "\t" << vIdx_0 << "\t" << vIdx_1 << "\t" << vIdx_2 << std::endl;
            }
            const Real dotProd = (normal.normalized()).dot(target_normal);
            if(dotProd < 0)
            {
                face2vertices(i,0) = vIdx_2;
                face2vertices(i,2) = vIdx_0;
                nFlippedNormals++;
            }
        }

        return nFlippedNormals;
    }

public:
    void get(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        getShellGeometry(vertices, face2vertices, vertices_bc);
    }

    bool analyticNormals() const override
    {
        return true;
    }
};


#endif /* Geometry_Base_h */
