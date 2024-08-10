//
//  common.hpp
//  Elasticity
//
//  Created by Wim van Rees on 2/15/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef common_hpp
#define common_hpp

#include <iostream>
#include <cassert>
#include <random>
#include <limits>
#include <utility>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <map>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <array>
#include <set>
#include <memory>
#include <ctime>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/StdVector>

#include <igl/barycentric_coordinates.h>
#include <igl/facet_components.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/remove_unreferenced.h>

//#include "ArgumentParser.h"

typedef double Real;
typedef unsigned long long tUint;
typedef std::vector<Eigen::Matrix2d, Eigen::aligned_allocator<Eigen::Matrix2d>> tVecMat2d;
typedef std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>> tVecMat3d;
enum MeshLayer {single, bottom, top};

//#define USETANTHETA

// little helper to avoid unused warnings when something is called in release mode and variable only checked in debug
#define _unused(x) ((void)(x))


namespace rnd
{
    /*! global random number generator */
    static std::mt19937 gen;
}


namespace helpers
{
    /*! convert an integer to a string with a specified number of zeros */
    inline std::string ToString(int value,int digitsCount)
    {
        std::ostringstream os;
        os<<std::setfill('0')<<std::setw(digitsCount)<<value;
        return os.str();
    }

    /*! convert an real to a string with a specified number of decimals */
    inline std::string ToStringPrecision(Real value, int digitsCount)
    {
        std::ostringstream out;
        out.precision(digitsCount);
        out << std::fixed << value;
        return out.str();
    }

    /*! check if two numbers are close (within some accuracy) */
    inline bool isclose(const Real val1, const Real val2, const bool printOnly=false, const Real rtol=1e-5, const Real atol=1e-8)
    {
        const Real diff = std::abs(val1-val2) ;
        bool retval = true;
        if(diff > (atol + rtol * std::abs(val2)))
            retval = false;
        if(diff > (atol + rtol * std::abs(val1)))
            retval = false;

        if(printOnly && not retval)
        {
            std::cout << "PROBLEM : " << val1 << "\t" << val2 << "\t" << diff << std::endl;
            retval = true;
        }

        return retval;
    }


    inline bool isPointInProjectedTriangle(const Real x, const Real y, const Eigen::Vector3d & v0, const Eigen::Vector3d & v1, const Eigen::Vector3d & v2)
    {
        // dont use the third dimension
        const Real xp[3] = {v0(0),v1(0),v2(0)};
        const Real yp[3] = {v0(1),v1(1),v2(1)};

        const int npoly = 3; // triangle

        int i,j,c=0;
        for (i = 0, j = npoly-1; i < npoly; j = i++) {
            if ((((yp[i] <= y) && (y < yp[j])) ||
                 ((yp[j] <= y) && (y < yp[i]))) &&
                (x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))
                c = !c;
        }
        return c;
    }


    inline Eigen::Vector3d getBaryCentricWeights2D(const Eigen::Vector3d & v0, const Eigen::Vector3d & v1, const Eigen::Vector3d & v2, const Real x, const Real y)
    {
        assert(isPointInProjectedTriangle(x, y, v0, v1, v2));

        Eigen::Vector3d me;
        me << x, y, 0;
        // compute my point in barycentric coordinates of undeformed configuration
        // Compute barycentric coordinates (u, v, w) for
        // point p with respect to triangle (a, b, c)
        // http://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
        const Eigen::Vector3d tmp0 = v1 - v0;
        const Eigen::Vector3d tmp1 = v2 - v0;
        const Eigen::Vector3d tmp2 = me - v0;
        const Real d00 = tmp0.dot(tmp0);
        const Real d01 = tmp0.dot(tmp1);
        const Real d11 = tmp1.dot(tmp1);
        const Real d20 = tmp2.dot(tmp0);
        const Real d21 = tmp2.dot(tmp1);
        const Real denom = d00 * d11 - d01 * d01;
        const Real lambda1 = (d11 * d20 - d01 * d21) / denom;
        const Real lambda2 = (d00 * d21 - d01 * d20) / denom;
        const Real lambda0 = 1.0 - lambda1 - lambda2;
        // if me == v0 --> lambda0 = 1
        // if me == v1 --> lambda1 = 1
        // if me == v2 --> lambda2 = 1

        Eigen::Vector3d retval;
        retval << lambda0, lambda1, lambda2;
        return retval;
    }

    template<int component>
    inline Real getDisplacementOfPointInRestTriangle(const Eigen::Vector3d & rv0, const Eigen::Vector3d & rv1, const Eigen::Vector3d & rv2, const Eigen::Vector3d & v0, const Eigen::Vector3d & v1, const Eigen::Vector3d & v2, const Real x, const Real y)
    {

        const Eigen::Vector3d lambdas = getBaryCentricWeights2D(rv0, rv1, rv2, x, y);

        const Real displ_0 = (v0 - rv0)(component);
        const Real displ_1 = (v1 - rv1)(component);
        const Real displ_2 = (v2 - rv2)(component);

        // get my displacement
        return -(lambdas(0)*displ_0 + lambdas(1)*displ_1 + lambdas(2)*displ_2);
    }


    inline Real interpolateOverTriangle(const Eigen::Vector3d & v0, const Eigen::Vector3d & v1, const Eigen::Vector3d & v2, const Real x, const Real y)
    {
        // linear interpolation over the triangle
        const Real A = v0(1)*(v1(2)-v2(2)) + v1(1)*(v2(2)-v0(2)) + v2(1)*(v0(2)-v1(2));
        const Real B = v0(2)*(v1(0)-v2(0)) + v1(2)*(v2(0)-v0(0)) + v2(2)*(v0(0)-v1(0));
        const Real C = v0(0)*(v1(1)-v2(1)) + v1(0)*(v2(1)-v0(1)) + v2(0)*(v0(1)-v1(1));
        const Real D = v0(0)*(v1(1)*v2(2)-v1(2)*v2(1)) + v1(0)*(v2(1)*v0(2)-v2(2)*v0(1)) + v2(0)*(v0(1)*v1(2)-v0(2)*v1(1));

        const Real ptZ = (D-A*x-B*y)/C;

        return ptZ;
    }

    inline Real computeDistanceToLineSegment(const Real x, const Real y, const Real xi, const Real yi, const Real xf, const Real yf)
    {
        //http://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment

        // case: line segment is degenerate
        const Real Lsq = (xf - xi) * (xf - xi) + (yf - yi) * (yf - yi);
        if(Lsq < std::numeric_limits<Real>::epsilon())
            return std::sqrt((xi - x) * (xi - x) + (yi - y) * (yi - y));

        Real a = x - xi;
        Real b = y - yi;
        Real c = xf - xi;
        Real d = yf - yi;

        Real dot = a * c + b * d;
        Real len_sq = c * c + d * d;
        Real temp = -1.0;
        if (len_sq != 0.0) { temp = dot / len_sq; }

        // Find the location on the segment that the point is closest to.
        Real xx;
        Real yy;

        if (temp < 0.0) {
            xx = xi;
            yy = yi;
        }
        else if (temp > 1.0) {
            xx = xf;
            yy = yf;
        }
        else {
            xx = xi + temp * c;
            yy = yi + temp * d;
        }

        Real dx = x - xx;
        Real dy = y - yy;
        Real dist = sqrt(dx * dx + dy * dy);

        return dist;
    }

    inline int getIntersectionPointBetweenVectors(const Eigen::Vector2d & pt1, const Eigen::Vector2d & n1, const Eigen::Vector2d & pt2, const Eigen::Vector2d & n2, Eigen::Vector2d & result)
    {
        // return 0 if intersection is found
        // return 1 if an intersection is found, but beyond the extent of pt1+n1 and pt2+n2
        // return 2 if no intersection is found (parallel vectors)

        const Eigen::Vector2d n1_hat = n1.normalized();
        const Eigen::Vector2d n2_hat = n2.normalized();

        const Real pt1_vec1 = pt1.dot(n1_hat);
        const Real pt1_vec2 = pt1.dot(n2_hat);
        const Real pt2_vec1 = pt2.dot(n1_hat);
        const Real pt2_vec2 = pt2.dot(n2_hat);
        const Real vec1_vec2 = n1_hat.dot(n2_hat);

        // deal with orthogonal
        if(std::abs(vec1_vec2 - 1) < std::numeric_limits<Real>::epsilon())
        {
            if((pt1-pt2).norm()  < std::numeric_limits<Real>::epsilon())
            {
                result = pt1;
                return 0;
            }
            else
            {
                return 2;
            }
        }

        const Real beta1 = (pt2_vec1 + (pt1_vec2 - pt2_vec2) * vec1_vec2 - pt1_vec1) / (1.0 - std::pow(vec1_vec2,2));
        const Real beta2 = (pt1_vec2 + (pt2_vec1 - pt1_vec1) * vec1_vec2 - pt2_vec2) / (1.0 - std::pow(vec1_vec2,2));

        //const Eigen::Vector2d retval1 = pt1 + beta1*n1_hat;
        const Eigen::Vector2d retval2 = pt2 + beta2*n2_hat;

        // set the intersection point
        result = retval2;

        const bool inRange_1 = (beta1 >= 0 && beta1 <= n1.norm());
        const bool inRange_2 = (beta2 >= 0 && beta2 <= n2.norm());
        return (inRange_1 && inRange_2) ? 0 : 1;
    }

    inline std::string removeExtension(const std::string input, const std::string extension)
    {
        // remove extension from input. Assume extension is given as (eg) .vtp (including the period)
        // if not found, return input unaltered
        const std::size_t found_ext = input.rfind(extension);
        const std::string output = (found_ext != std::string::npos ? input.substr(0, found_ext) : input);
        return output;
    }

    inline void catastrophe(std::string s, const char* file, const int line, const bool isFatal=true)    // write ``error: s and exit program
    {
        std::cerr << "Something went wrong, error: " << s << std::endl;
        std::cerr << "Called from file : " << file << " at line number " << line << std::endl;
        if(isFatal)
        {
            std::cout << "Exiting..." << std::endl;
            std::exit(1);
        }
    }

    template<class Matrix>
    inline void write_matrix_binary(const std::string & filename, const Matrix& matrix)
    {
        // from https://stackoverflow.com/a/25389481
        std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
        typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
        out.write((char*) (&rows), sizeof(typename Matrix::Index));
        out.write((char*) (&cols), sizeof(typename Matrix::Index));
        out.write((char*) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar) );
        out.close();
    }

    template<class Matrix>
    inline void read_matrix_binary(const std::string & filename, Matrix& matrix)
    {
        // from https://stackoverflow.com/a/25389481
        std::ifstream in(filename, std::ios::in | std::ios::binary);
        typename Matrix::Index rows=0, cols=0;
        in.read((char*) (&rows),sizeof(typename Matrix::Index));
        in.read((char*) (&cols),sizeof(typename Matrix::Index));
        matrix.resize(rows, cols);
        in.read( (char *) matrix.data() , rows*cols*sizeof(typename Matrix::Scalar) );
        in.close();
    }

    template<class Matrix>
    inline void write_vecmat_binary(const std::string & filename, const std::vector<Matrix, Eigen::aligned_allocator<Matrix>> & vecmat)
    {
        // from https://stackoverflow.com/a/25389481
        std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
        const size_t vecsize = vecmat.size();
        if(vecsize==0) return;
        typename Matrix::Index rows=vecmat[0].rows(), cols=vecmat[0].cols();
        out.write((char*) &vecsize, sizeof(vecsize));
        out.write((char*) (&rows), sizeof(typename Matrix::Index));
        out.write((char*) (&cols), sizeof(typename Matrix::Index));
        out.write((char*)&vecmat[0], vecsize * rows*cols*sizeof(typename Matrix::Scalar));
        out.close();
    }

    template<class Matrix>
    inline void read_vecmat_binary(const std::string & filename, std::vector<Matrix, Eigen::aligned_allocator<Matrix>> & vecmat)
    {
        std::ifstream in(filename, std::ios::in | std::ios::binary);
        size_t vecsize = 0;
        typename Matrix::Index rows=0, cols=0;

        in.read((char*) (&vecsize),sizeof(vecsize));
        assert(vecsize > 0);

        in.read((char*) (&rows),sizeof(typename Matrix::Index));
        in.read((char*) (&cols),sizeof(typename Matrix::Index));

        assert(rows == Matrix::RowsAtCompileTime);
        assert(cols == Matrix::ColsAtCompileTime);

        vecmat.resize(vecsize);
        in.read( (char *) &vecmat[0] , vecsize * rows*cols*sizeof(typename Matrix::Scalar) );
        in.close();
    }
}


namespace Eigen
{
    /*! typedef vector of booleans */
    typedef Matrix<bool, 3, 1> Vector3b;

    /*! typedef vector of booleans */
    typedef Matrix<bool, Dynamic, 1> VectorXb;
    typedef Matrix<bool, Dynamic, Dynamic> MatrixXb;
}


inline std::string getDateTime()
{
   std::time_t now = std::time(0);
   std::tm *ltm = std::localtime(&now);
   std::ostringstream os;
   os << std::setfill('0') << std::setw(2) << 1 + ltm->tm_mon
      << std::setfill('0') << std::setw(2) << ltm->tm_mday
      << std::setfill('0') << std::setw(2) << 1 + ltm->tm_hour
      << std::setfill('0') << std::setw(2) << 1 + ltm->tm_min
      << std::setfill('0') << std::setw(2) << 1 + ltm->tm_sec;
   return os.str();
}


inline Eigen::MatrixXd getDataFromFile(const std::string filename, const int nColumns)
{
    std::ifstream file(filename);
    double temp;
    std::vector<Real> temp_values;
    if (file.is_open()) {
        while (true) {
            file >> temp;
            temp_values.push_back(temp);
            if (file.eof()) { break; }
        }
        file.close();
    }
    else {
        std::cout << "Error reading file '" << filename << "'." << std::endl;
        exit(1);
    }

    Eigen::MatrixXd data(temp_values.size() / nColumns, nColumns);
    for (int i = 0; i < data.rows(); i++) {
        for (int j = 0; j < data.cols(); j++) {
            data(i, j) = temp_values[i * data.cols() + j];
        }
    }

    return data;
}


inline int which(const int element, const Eigen::VectorXi list)
{
    for (int i = 0; i < list.size(); ++i) {
        if (list(i) == element) { return i; }
    }
    return -1;
}


inline int which(const int element, std::vector<int> list)
{
    for (int i = 0; i < list.size(); ++i) {
        if (list[i] == element) { return i; }
    }
    return -1;
}


inline void printMathematicaArray(const Eigen::MatrixXd & matrix)
{
    int nCols = 2;
    std::cout << std::fixed << "{";
    for (int i = 0; i < matrix.rows(); i++) {
        std::cout << "{";
        for (int j = 0; j < nCols; j++) {
            std::cout << matrix(i, j) << (j < nCols - 1 ? "," : "");
        }
        std::cout << "}" << (i < matrix.rows() - 1 ? "," : "}");
    }
    std::cout << std::endl;
}


inline void dumpOBJ(const std::string tag, const Eigen::MatrixXd & vertices, const Eigen::MatrixXi & faces)
{
    FILE * f = fopen((tag + ".obj").c_str(), "w");

    for (int i = 0; i < vertices.rows(); i++) {
        fprintf(f, "v ");
        for (int j = 0; j < vertices.cols(); j++) {
            fprintf(f, "%f ", vertices(i, j));
        }
        if (vertices.cols() < 3) {
            fprintf(f, "0.0");
        }
        fprintf(f, "\n");
    }

    for (int i = 0; i < faces.rows(); i++) {
        fprintf(f, "f ");
        for (int j = 0; j < faces.cols(); j++) {
            fprintf(f, "%d ", faces(i, j) + 1); // faces reference zero-indexed vertices, so we need to increment them by one
        }
        fprintf(f, "\n");
    }

    fclose(f);
}


// sphere centered on origin with radius = R projected onto plane z = R
// from the point (0, 0, -R)
inline void conformalMap(const Real radius, const Eigen::MatrixXd & xyz, Eigen::MatrixXd & uv)
{
    // (x, y, z) --> (R, theta, phi) --> (u, v)
    uv.resize(xyz.rows(), 2);

    Eigen::ArrayXd rho = (xyz.array().col(0).square() + xyz.array().col(1).square()).sqrt();
    Eigen::ArrayXd ratio = ((radius - xyz.array().col(2)) / (radius + xyz.array().col(2))).sqrt();
    Eigen::ArrayXd temp = 2.0 * radius * ratio / rho;

    uv.col(0) = xyz.array().col(0) * temp;
    uv.col(1) = xyz.array().col(1) * temp;
}


inline void inverseConformalMap(const Real radius, const Eigen::MatrixXd & uv, Eigen::MatrixXd & xyz)
{
    xyz.resize(uv.rows(), 3);

    Real t = 4.0 * radius * radius;
    Eigen::ArrayXd temp = t + uv.array().col(0).square() + uv.array().col(1).square();
    xyz.col(0) = t * uv.array().col(0) / temp;
    xyz.col(1) = t * uv.array().col(1) / temp;
    xyz.col(2) = radius * (2.0 * t / temp - 1.0);
}


inline bool pointInsidePolygon(const Eigen::MatrixXd & polygonVertices, const Eigen::VectorXd & point)
{
    // Polygon must be defined in an oriented fashion. See https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html
    int j = polygonVertices.rows() - 1;
    bool inside = false;
    bool flippedYet = false;
    for (int i = 0; i < polygonVertices.rows(); i++) {
        if ((polygonVertices(i, 1) > point(1)) != (polygonVertices(j, 1) > point(1))) {
            Real x = (polygonVertices(j, 0) - polygonVertices(i, 0)) * (point(1) - polygonVertices(i, 1)) / (polygonVertices(j, 1) - polygonVertices(i, 1)) + polygonVertices(i, 0);
            if (point(0) < x) {
                inside = !inside;
            }
        }
        j = i;
    }
    return inside;
}


inline bool getSegmentIntersection(const Real p1x, const Real p1y, const Real p2x, const Real p2y, const Real q1x, const Real q1y, const Real q2x, const Real q2y, Real &x, Real &y)
{
    // Check for an intersection between the segments P and Q. If one is found,
    // store the intersection point in (x, y). From https://stackoverflow.com/a/1968345
    Real s1_x = p2x - p1x;
    Real s1_y = p2y - p1y;
    Real s2_x = q2x - q1x;
    Real s2_y = q2y - q1y;

    Real s = (-s1_y * (p1x - q1x) + s1_x * (p1y - q1y)) / (-s2_x * s1_y + s1_x * s2_y);
    Real t = ( s2_x * (p1y - q1y) - s2_y * (p1x - q1x)) / (-s2_x * s1_y + s1_x * s2_y);

    if (s >= 0 && s <= 1 && t >= 0 && t <= 1) { // Collision detected
        x = p1x + (t * s1_x);
        y = p1y + (t * s1_y);
        return true;
    }

    return false; // No collision
}


inline bool getSegmentIntersection(const Real p1x, const Real p1y, const Real p2x, const Real p2y, const Real q1x, const Real q1y, const Real q2x, const Real q2y)
{
    Real x, y;
    return getSegmentIntersection(p1x, p1y, p2x, p2y, q1x, q1y, q2x, q2y, x, y);
}


inline bool getSegmentIntersection(const Eigen::VectorXd & p1, const Eigen::VectorXd & p2, const Eigen::VectorXd & q1, const Eigen::VectorXd & q2, Real &x, Real &y)
{
    return getSegmentIntersection(p1(0), p1(1), p2(0), p2(1), q1(0), q1(1), q2(0), q2(1), x, y);
}


inline bool getSegmentIntersection(const Eigen::VectorXd & p1, const Eigen::VectorXd & p2, const Eigen::VectorXd & q1, const Eigen::VectorXd & q2)
{
    Real x, y;
    return getSegmentIntersection(p1(0), p1(1), p2(0), p2(1), q1(0), q1(1), q2(0), q2(1), x, y);
}


inline void removePolygonFromMesh(Eigen::MatrixXd & vertices, Eigen::MatrixXi & faces, Eigen::MatrixXb & vertices_bc, const Eigen::MatrixXd & polygonVertices)
{
    // Remove the part of a mesh that lies within some polygon
    int nKeepFaces = 0;
    Eigen::VectorXi trianglesPerVertex = Eigen::VectorXi::Constant(vertices.rows(), 0);
    int nKeepVertices = 0;
    Eigen::VectorXi vertexMap(vertices.rows());

    // Only keep faces that don't lie within the polygon
    for (int i = 0; i < faces.rows(); i++) {
        Eigen::Vector2d faceCenter((vertices(faces(i, 0), 0) + vertices(faces(i, 1), 0) + vertices(faces(i, 2), 0)) / 3.0,
                                   (vertices(faces(i, 0), 1) + vertices(faces(i, 1), 1) + vertices(faces(i, 2), 1)) / 3.0);
        if (!pointInsidePolygon(polygonVertices, faceCenter)) {
            for (int j = 0; j < 3; j++) { trianglesPerVertex(faces(i, j)) += 1; }
            faces.row(nKeepFaces) = faces.row(i);
            nKeepFaces++;
        }
    }
    faces.conservativeResize(nKeepFaces, faces.cols());

    // Only keep vertices (and their boundary conditions) that are still referenced by a triangle
    for (int i = 0; i < vertices.rows(); i++) {
        if (trianglesPerVertex(i) > 0) {
            vertices.row(nKeepVertices) = vertices.row(i);
            vertices_bc.row(nKeepVertices) = vertices_bc.row(i);
            vertexMap(i) = nKeepVertices; // new index
            nKeepVertices++;
        }
        else {
            vertexMap(i) = -1;
        }
    }
    vertices.conservativeResize(nKeepVertices, vertices.cols());
    vertices_bc.conservativeResize(nKeepVertices, vertices_bc.cols());

    // Update vertex indices referenced by the faces
    for (int i = 0; i < faces.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            if (vertexMap(faces(i, j)) < 0) {
                std::cout << "Error: old vertex " << faces(i, j) << " was indexed by an updated face!" << std::endl;
                exit(1);
            }
            faces(i, j) = vertexMap(faces(i, j));
        }
    }
}


inline void keepPolygonFromMesh(Eigen::MatrixXd & vertices, Eigen::MatrixXi & faces, Eigen::MatrixXb & vertices_bc, const Eigen::MatrixXd & polygonVertices)
{
    // Keep only the part of a mesh that lies within some polygon
    int nKeepFaces = 0;
    Eigen::VectorXi trianglesPerVertex = Eigen::VectorXi::Constant(vertices.rows(), 0);
    int nKeepVertices = 0;
    Eigen::VectorXi vertexMap(vertices.rows());

    // Only keep faces that lie within the polygon
    for (int i = 0; i < faces.rows(); i++) {
        Eigen::Vector2d faceCenter((vertices(faces(i, 0), 0) + vertices(faces(i, 1), 0) + vertices(faces(i, 2), 0)) / 3.0,
                                   (vertices(faces(i, 0), 1) + vertices(faces(i, 1), 1) + vertices(faces(i, 2), 1)) / 3.0);
        if (pointInsidePolygon(polygonVertices, faceCenter)) {
            for (int j = 0; j < 3; j++) { trianglesPerVertex(faces(i, j)) += 1; }
            faces.row(nKeepFaces) = faces.row(i);
            nKeepFaces++;
        }
    }
    faces.conservativeResize(nKeepFaces, faces.cols());

    // Only keep vertices (and their boundary conditions) that are still referenced by a triangle
    for (int i = 0; i < vertices.rows(); i++) {
        if (trianglesPerVertex(i) > 0) {
            vertices.row(nKeepVertices) = vertices.row(i);
            vertices_bc.row(nKeepVertices) = vertices_bc.row(i);
            vertexMap(i) = nKeepVertices; // new index
            nKeepVertices++;
        }
        else {
            vertexMap(i) = -1;
        }
    }
    vertices.conservativeResize(nKeepVertices, vertices.cols());
    vertices_bc.conservativeResize(nKeepVertices, vertices_bc.cols());

    // Update vertex indices referenced by the faces
    for (int i = 0; i < faces.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            if (vertexMap(faces(i, j)) < 0) {
                std::cout << "Error: old vertex " << faces(i, j) << " was indexed by an updated face!" << std::endl;
                exit(1);
            }
            faces(i, j) = vertexMap(faces(i, j));
        }
    }
}


inline void removeUnconnectedComponents(Eigen::MatrixXd & vertices, Eigen::MatrixXi & faces, Eigen::MatrixXb & vertices_bc, int vIndex = 0)
{
    Eigen::VectorXi F1;
    igl::facet_components(faces, F1);
    int fComponent = 0;
    for (int i = 0; i < faces.rows(); i++) {
        if (faces(i, 0) == vIndex || faces(i, 1) == vIndex || faces(i, 2) == vIndex) {
            fComponent = F1(i);
            break;
        }
    }

    // Only keep faces that correspond to this component
    int nKeepFaces = 0;
    Eigen::VectorXi trianglesPerVertex = Eigen::VectorXi::Constant(vertices.rows(), 0);
    for (int i = 0; i < faces.rows(); i++) {
        if (F1(i) == fComponent) {
            for (int j = 0; j < 3; j++) { trianglesPerVertex(faces(i, j)) += 1; }
            faces.row(nKeepFaces) = faces.row(i);
            nKeepFaces++;
        }
    }
    faces.conservativeResize(nKeepFaces, faces.cols());

    // Only keep vertices (and their boundary conditions) that are still referenced by a triangle
    int nKeepVertices = 0;
    Eigen::VectorXi vertexMap(vertices.rows());
    for (int i = 0; i < vertices.rows(); i++) {
        if (trianglesPerVertex(i) > 0) {
            vertices.row(nKeepVertices) = vertices.row(i);
            vertices_bc.row(nKeepVertices) = vertices_bc.row(i);
            vertexMap(i) = nKeepVertices; // new index
            nKeepVertices++;
        }
        else {
            vertexMap(i) = -1;
        }
    }
    vertices.conservativeResize(nKeepVertices, vertices.cols());
    vertices_bc.conservativeResize(nKeepVertices, vertices_bc.cols());

    // Update vertex indices referenced by the faces
    for (int i = 0; i < faces.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            if (vertexMap(faces(i, j)) < 0) {
                std::cout << "Error: old vertex " << faces(i, j) << " was indexed by an updated face!" << std::endl;
                exit(1);
            }
            faces(i, j) = vertexMap(faces(i, j));
        }
    }
}


// Open up a cut in the mesh along a string of vertices.
inline void cutMesh(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc, std::vector< std::vector<int> > & vertex2faces, const Eigen::VectorXi & cutVertices, Eigen::VectorXi & duplicateMeshIndices, bool duplicateRightSide = true)
{
    const int nOriginalVertices = vertices.rows();
    const int nCutVertices = cutVertices.size();
    vertices.conservativeResize(nOriginalVertices + nCutVertices - 2, vertices.cols());
    vertices_bc.conservativeResize(nOriginalVertices + nCutVertices - 2, vertices_bc.cols());
    duplicateMeshIndices.conservativeResize(nCutVertices - 2);

    int prev, curr, next, dupl, ip, ic, io, prev_checkpoint;
    prev = cutVertices(0);

    bool run_backward = false;

    for (int i = 1; i < cutVertices.size() - 1; i++) {
        curr = cutVertices(i);
        // std::cout << "\ncurrent vertex: " << curr << " at " << vertices.row(curr) << std::endl;
        next = cutVertices(i + 1);

        // duplicate current vertex
        dupl = nOriginalVertices + i - 1;
        vertices.row(dupl) = vertices.row(curr);
        vertices_bc.row(dupl) = vertices_bc.row(curr);
        vertex2faces.push_back(std::vector<int>());
        if (vertices.cols() == 3) {
            vertices(curr, 2) = 0.001;
            vertices(dupl, 2) = -0.001;
        }

        duplicateMeshIndices(i - 1) = dupl;
        // std::cout << "previous vertex: " << prev << " at " << vertices.row(prev) << std::endl;

        while (prev != next) {
            prev_checkpoint = prev;

            for (auto f : vertex2faces[curr]) {
                // std::cout << "checking face with vertices " << face2vertices.row(f) << std::endl;

                ip = -1;
                ic = -1;

                for (int j = 0; j < 3; j++) {
                    if (face2vertices(f, j) == prev) { ip = j; }
                    else if (face2vertices(f, j) == curr) { ic = j; }
                }

                // face is not adjacent to this segment (prev --> curr), skip
                if (ip < 0) {
                    // std::cout << "\tface not adjacent to prev --> curr" << std::endl;
                    continue;
                }

                // third corner
                io = 3 - ip - ic;

                // only want faces on right/left side of cut; parity of triangle (prev --> io --> curr) should be even/odd
                Real parity = (vertices(face2vertices(f, io), 0) - vertices(face2vertices(f, ip), 0)) * (vertices(face2vertices(f, ic), 1) - vertices(face2vertices(f, io), 1))
                              - (vertices(face2vertices(f, io), 1) - vertices(face2vertices(f, ip), 1)) * (vertices(face2vertices(f, ic), 0) - vertices(face2vertices(f, io), 0));
                if (duplicateRightSide && parity < 0.0) { continue; }
                if (!duplicateRightSide && parity > 0.0) { continue; }

                // update so this face points to duplicate vertex instead of curr
                face2vertices(f, ic) = dupl;
                // std::cout << "\t-- updated face to point to vertices " << face2vertices.row(f) << std::endl;
                vertex2faces[curr].erase(std::remove(vertex2faces[curr].begin(), vertex2faces[curr].end(), f), vertex2faces[curr].end());
                vertex2faces[dupl].push_back(f);
                prev = face2vertices(f, io);
                // std::cout << "\t-- new 'prev' vertex is " << prev << " at " << vertices.row(prev) << std::endl;
                break;
            }

            // if we are stuck, set prev = next, flip the parity requirement, and go backward
            if (prev != next && prev == prev_checkpoint) {
                // std::cout << std::endl << "prev is not updating from the forward pass. trying to go backward." << std::endl << std::endl;
                run_backward = true;
                break;
            }
        }

        // we need to do a second pass, going backward
        if (run_backward) {
            prev = next; // starting point
            // std::cout << "previous vertex: " << prev << " now at " << vertices.row(prev) << std::endl;

            while (true) {
                prev_checkpoint = prev;

                for (auto f : vertex2faces[curr]) {
                    // std::cout << "checking face with vertices " << face2vertices.row(f) << std::endl;

                    ip = -1;
                    ic = -1;

                    for (int j = 0; j < 3; j++) {
                        if (face2vertices(f, j) == prev) { ip = j; }
                        else if (face2vertices(f, j) == curr) { ic = j; }
                    }

                    // face is not adjacent to this segment (prev --> curr), skip
                    if (ip < 0) {
                        // std::cout << "\tface not adjacent to prev --> curr" << std::endl;
                        continue;
                    }

                    // third corner
                    io = 3 - ip - ic;

                    // parity conditions flipped wrt. above
                    Real parity = (vertices(face2vertices(f, io), 0) - vertices(face2vertices(f, ip), 0)) * (vertices(face2vertices(f, ic), 1) - vertices(face2vertices(f, io), 1))
                                  - (vertices(face2vertices(f, io), 1) - vertices(face2vertices(f, ip), 1)) * (vertices(face2vertices(f, ic), 0) - vertices(face2vertices(f, io), 0));
                    if (duplicateRightSide && parity > 0.0) { continue; }
                    if (!duplicateRightSide && parity < 0.0) { continue; }

                    // update so this face points to duplicate vertex instead of curr
                    face2vertices(f, ic) = dupl;
                    // std::cout << "\t-- updated face to point to vertices " << face2vertices.row(f) << std::endl;
                    vertex2faces[curr].erase(std::remove(vertex2faces[curr].begin(), vertex2faces[curr].end(), f), vertex2faces[curr].end());
                    vertex2faces[dupl].push_back(f);
                    prev = face2vertices(f, io);
                    // std::cout << "\t-- new 'prev' vertex is " << prev << " at " << vertices.row(prev) << std::endl;
                    break;
                }

                // we looked through all the faces but prev didn't update ==> we're done
                if (prev == prev_checkpoint) {
                    break;
                }
            }
        }

        prev = dupl;
    }
}


inline Eigen::MatrixXd displaceCurve(const Eigen::MatrixXd & points, const Real displacement)
{
    // Displacement is by default toward the right
    int nPoints = points.rows();
    Eigen::MatrixXd displacedPoints(nPoints, points.cols());

    // Propagate each point in the curve toward the right
    for (int i = 0; i < nPoints; i++) {
        Real dx, dy;

        if (i == 0) {
            Real next_length = (points.row(i + 1) - points.row(i)).norm();
            dx = displacement * (points(i + 1, 1) - points(i, 1)) / next_length;
            dy = -displacement * (points(i + 1, 0) - points(i, 0)) / next_length;
        }
        else if (i < nPoints - 1) {
            Real prev_angle = std::atan2(points(i - 1, 1) - points(i, 1), points(i - 1, 0) - points(i, 0)); // angle to previous point
            Real next_angle = std::atan2(points(i + 1, 1) - points(i, 1), points(i + 1, 0) - points(i, 0)); // angle to next point
            Real diff_angle = std::fmod(next_angle - prev_angle + 2.0 * M_PI, 2.0 * M_PI);
            Real theta = prev_angle + 0.5 * diff_angle;

            dx = displacement * std::cos(theta);
            dy = displacement * std::sin(theta);
        }
        else {
            Real prev_length = (points.row(i) - points.row(i - 1)).norm();
            dx = displacement * (points(i, 1) - points(i - 1, 1)) / prev_length;
            dy = -displacement * (points(i, 0) - points(i - 1, 0)) / prev_length;
        }

        displacedPoints(i, 0) = points(i, 0) + dx;
        displacedPoints(i, 1) = points(i, 1) + dy;
    }

    // Check new curve for loops
    int i = 0;
    while (i < displacedPoints.rows() - 1) {
        // Check for intersections, starting with end segment
        for (int j = displacedPoints.rows() - 2; j > i + 1; j--) {
            Real x, y;
            bool intersection = getSegmentIntersection(displacedPoints.row(i), displacedPoints.row(i + 1), displacedPoints.row(j), displacedPoints.row(j + 1), x, y);

            if (intersection) { // replace segments in between with intersection point
                displacedPoints(j, 0) = x;
                displacedPoints(j, 1) = y;

                int n = j - i - 1;
                for (int k = i + 1; k < displacedPoints.rows() - n; k++) {
                    displacedPoints.row(k) = displacedPoints.row(k + n);
                }

                displacedPoints.conservativeResize(displacedPoints.rows() - n, displacedPoints.cols());
                break;
            }
        }

        i++;
    }

    return displacedPoints;
}


inline Eigen::MatrixXd rectifyCurve(const Eigen::MatrixXd & points)
{
    // Rectify y >= 0
    Eigen::MatrixXd newPoints = points;

    Real xmin = newPoints.col(0).minCoeff() - 1.0;
    Real xmax = newPoints.col(0).maxCoeff() + 1.0;

    // Insert points where segments intersect the y-axis
    int i = 0;
    while (i < newPoints.rows() - 1) {
        Real x, y;
        bool intersection = getSegmentIntersection(newPoints.row(i), newPoints.row(i + 1), Eigen::Vector2d(xmin, 0.0), Eigen::Vector2d(xmax, 0.0), x, y);
        if (intersection) {
            int k = newPoints.rows() - i - 1;
            newPoints.conservativeResize(newPoints.rows() + 1, newPoints.cols());
            newPoints.block(i + 2, 0, k, newPoints.cols()) = newPoints.block(i + 1, 0, k, newPoints.cols());
            newPoints(i + 1, 0) = x;
            newPoints(i + 1, 1) = y;
            i++;
        }
        i++;
    }

    // Remove all points that sit below the y-axis
    i = 0;
    while (i < newPoints.rows()) {
        if (newPoints(i, 1) < -1e-6) {
            int k = newPoints.rows() - i - 1;
            newPoints.block(i, 0, k, newPoints.cols()) = newPoints.block(i + 1, 0, k, newPoints.cols());
            newPoints.conservativeResize(newPoints.rows() - 1, newPoints.cols());
            continue;
        }
        i++;
    }

    return newPoints;
}


inline void deleteSmallSegments(Eigen::MatrixXd & points, const Real eps = 1e-6)
{
    int i = 0;
    while (i < points.rows() - 1) {
        Real length = (points.row(i + 1) - points.row(i)).norm();
        if (length < eps) {
            points.row(i) = 0.5 * (points.row(i) + points.row(i + 1));

            // Delete (i + 1)th row
            int k = points.rows() - i - 2;
            points.block(i + 1, 0, k, points.cols()) = points.block(i + 2, 0, k, points.cols());
            points.conservativeResize(points.rows() - 1, points.cols());
        }
        else {
            i++;
        }
    }
}


template <typename tContainer>
inline Eigen::VectorXi getIndexMap(tContainer & selectedIndices, int nTotalElements)
{
    Eigen::VectorXi map = Eigen::VectorXi::Constant(nTotalElements, -1);
    int ind = 0;
    for (auto selectedInd : selectedIndices) {
        map(selectedInd) = ind;
        ind++;
    }
    return map;
}


// Given an index map, returns the collapsed indices of the resulting vector
// once the map has been applied.
inline Eigen::VectorXi collapseMap(const Eigen::VectorXi & indexMap)
{
    Eigen::VectorXi map = Eigen::VectorXi::Constant(indexMap.size(), -1);
    int ci = 0;
    for (int mi = 0; mi < indexMap.size(); mi++) {
        if (indexMap(mi) < 0) { // this element will NOT be collapsed
            map(mi) = ci;
            ci++;
        }
    }
    return map;
}


// Go through all elements of a matrix, re-mapping their values.
inline void mapValues(const Eigen::VectorXi & indexMap, Eigen::MatrixXi & elements)
{
    for (int i = 0; i < elements.rows(); i++) {
        for (int j = 0; j < elements.cols(); j++) {
            if (indexMap(elements(i, j)) >= 0) {
                elements(i, j) = indexMap(elements(i, j));
            }
        }
    }
}


template <typename DerivedType>
inline void mapRows(const Eigen::VectorXi & indexMap, Eigen::PlainObjectBase<DerivedType> & elements, const bool removeDuplicates = false)
{
    assert(indexMap.size() == elements.rows());

    if (!removeDuplicates) {
        for (int i = 0; i < indexMap.size(); i++) {
            if (indexMap(i) >= 0) {
                elements.row(i) = elements.row(indexMap(i));
            }
        }
    }
    else {
        int nDuplicates = 0;
        for (int i = 0; i < indexMap.size(); i++) {
            if (indexMap(i) >= 0) {
                elements.row(indexMap(i)) = elements.row(i); // this is reversed from the previous case
            }
            else {
                nDuplicates++;
            }
        }
        elements.conservativeResize(elements.rows() - nDuplicates, elements.cols());
    }
}


inline Eigen::MatrixXi extractRows(const std::vector<int> & selectedRows, const Eigen::MatrixXi & elements)
{
    Eigen::MatrixXi extractedElements(selectedRows.size(), elements.cols());
    int ind = 0;
    for (auto mappedInd : selectedRows) {
        extractedElements.row(ind) = elements.row(mappedInd);
        ind++;
    }
    return extractedElements;
}


template <typename tMatrix>
inline tMatrix extractRows(const Eigen::VectorXi & selectedRows, const tMatrix & elements)
{
    tMatrix extractedElements(selectedRows.size(), elements.cols());
    int ind = 0;
    for (int i = 0; i < selectedRows.size(); i++) {
        int mappedInd = selectedRows(i);
        extractedElements.row(ind) = elements.row(mappedInd);
        ind++;
    }
    return extractedElements;
}


template <typename tMatrix>
inline tMatrix keepRows(const Eigen::VectorXi & keepRows_list, const tMatrix & elements)
{
    int nRows = elements.rows();
    tMatrix extractedElements(nRows, elements.cols());
    int ind = 0;
    for (int i = 0; i < nRows; i++) {
        if (keepRows_list(i) >= 0) {
            extractedElements.row(ind) = elements.row(i);
            ind++;
        }
    }
    extractedElements.conservativeResize(ind, elements.cols());
    return extractedElements;
}


// combine two meshes into one, zipping up vertices that are very close to each other.
inline void combineMeshes(const Eigen::MatrixXd & vertices1, const Eigen::MatrixXi & faces1,
                          const Eigen::MatrixXd & vertices2, const Eigen::MatrixXi & faces2,
                          Eigen::MatrixXd & vertices, Eigen::MatrixXi & faces,
                          const Real eps = 1e-6)
{
    assert(vertices1.cols() == vertices2.cols());
    assert(faces1.cols() == faces2.cols());
    Eigen::MatrixXd vertices_(vertices1.rows() + vertices2.rows(), vertices1.cols());
    Eigen::MatrixXi faces_(faces1.rows() + faces2.rows(), faces1.cols());

    vertices_.block(0, 0, vertices1.rows(), vertices1.cols()) = vertices1;
    vertices_.block(vertices1.rows(), 0, vertices2.rows(), vertices1.cols()) = vertices2;

    faces_.block(0, 0, faces1.rows(), faces1.cols()) = faces1;
    faces_.block(faces1.rows(), 0, faces2.rows(), faces1.cols()) = faces2;

    for (int i = 0; i < faces2.rows(); i++) {
        for (int j = 0; j < faces2.cols(); j++) {
            faces_(faces1.rows() + i, j) += vertices1.rows();
        }
    }

    // remove "duplicate" vertices
    Eigen::VectorXi indexMap = Eigen::VectorXi::Constant(vertices_.rows(), -1);
    for (int i = 0; i < vertices_.rows(); i++) {
        for (int j = 0; j < i; j++) {
            if (vertices1.cols() == 2) {
                if (std::abs(vertices_(i, 0) - vertices_(j, 0)) < eps
                    && std::abs(vertices_(i, 1) - vertices_(j, 1)) < eps) {
                    indexMap(j) = i;
                }
            }
            else if (vertices1.cols() == 3) {
                if (std::abs(vertices_(i, 0) - vertices_(j, 0)) < eps
                    && std::abs(vertices_(i, 1) - vertices_(j, 1)) < eps
                    && std::abs(vertices_(i, 2) - vertices_(j, 2)) < eps) {
                    indexMap(j) = i;
                }
            }
        }
    }

    // Re-map each face's referenced vertices
    mapValues(indexMap, faces_);

    // Remove unreferenced vertices (indexMap used as a dummy variable)
    igl::remove_unreferenced(vertices_, faces_, vertices, faces, indexMap);
}


// zip up vertices that are very close to each other.
inline void combineMesh(Eigen::MatrixXd & vertices, Eigen::MatrixXi & faces, const Real eps = 1e-6)
{
    Eigen::MatrixXd vertices_ = vertices;
    Eigen::MatrixXi faces_ = faces;

    Eigen::VectorXi indexMap = Eigen::VectorXi::Constant(vertices_.rows(), -1);
    for (int i = 0; i < vertices_.rows(); i++) {
        for (int j = i + 1; j < vertices_.rows(); j++) {
            if (indexMap(i) < 0) {
                if (std::abs(vertices_(i, 0) - vertices_(j, 0)) < eps
                    && std::abs(vertices_(i, 1) - vertices_(j, 1)) < eps
                    && std::abs(vertices_(i, 2) - vertices_(j, 2)) < eps) {
                    indexMap(j) = i;
                }
            }
        }
    }

    // Re-map each face's referenced vertices
    mapValues(indexMap, faces_);

    // Remove unreferenced vertices (indexMap used as a dummy variable)
    igl::remove_unreferenced(vertices_, faces_, vertices, faces, indexMap);
}


inline Eigen::MatrixXd mapPointsByMesh(const Eigen::MatrixXd & rpoints, const Eigen::MatrixXd & rvertices_old, const Eigen::MatrixXd & cvertices_old, const Eigen::MatrixXi & faces_old)
{
    Eigen::MatrixXd newPoints;
    const int n = rpoints.rows();
    const int m = rpoints.cols();

    // Get the closest positions of these points to the old mesh
    Eigen::VectorXd sq_distance(n);
    Eigen::VectorXi faceIndices(n);
    Eigen::MatrixXd closestPositions(n, 3);
    if (m == 2) {
        Eigen::MatrixXd points_ = Eigen::MatrixXd::Constant(n, 3, 0.0);
        points_.leftCols(2) = rpoints;
        igl::point_mesh_squared_distance(points_, rvertices_old, faces_old, sq_distance, faceIndices, closestPositions);
    }
    else {
        igl::point_mesh_squared_distance(rpoints, rvertices_old, faces_old, sq_distance, faceIndices, closestPositions);
    }

    // Project from the old mesh to the new mesh
    Eigen::MatrixXi faceVertices = extractRows(faceIndices, faces_old);

    Eigen::MatrixXd barycentricPositions;
    Eigen::MatrixXd a = extractRows(faceVertices.col(0), rvertices_old);
    Eigen::MatrixXd b = extractRows(faceVertices.col(1), rvertices_old);
    Eigen::MatrixXd c = extractRows(faceVertices.col(2), rvertices_old);
    igl::barycentric_coordinates(closestPositions, a, b, c, barycentricPositions);

    Eigen::MatrixXd aa = extractRows(faceVertices.col(0), cvertices_old);
    Eigen::MatrixXd bb = extractRows(faceVertices.col(1), cvertices_old);
    Eigen::MatrixXd cc = extractRows(faceVertices.col(2), cvertices_old);

    newPoints.resize(n, 3);
    for (int i = 0; i < n; i++) {
        newPoints.row(i) =   barycentricPositions(i, 0) * aa.row(i)
                           + barycentricPositions(i, 1) * bb.row(i)
                           + barycentricPositions(i, 2) * cc.row(i);
    }

    return newPoints;
}


inline Real enclosedVolume(const Eigen::MatrixXd & vertices, const Eigen::MatrixXi & faces)
{
    // https://rosenzweig.io/blog/hilariously-fast-volume-computation-with-the-divergence-theorem.html
    Real volume = 0.0;
    for (int i = 0; i < faces.rows(); i++) {
        int i0 = faces(i, 0);
        int i1 = faces(i, 1);
        int i2 = faces(i, 2);

        Real cross = (vertices(i1, 1) - vertices(i0, 1)) * (vertices(i2, 2) - vertices(i0, 2))
                     - (vertices(i1, 2) - vertices(i0, 2)) * (vertices(i2, 1) - vertices(i0, 1));
        Real integral = vertices(i0, 0) + vertices(i1, 0) + vertices(i2, 0);

        volume = volume + cross * integral;
    }
    return volume / 6.0;
}

inline Eigen::MatrixXi toInt(const Eigen::MatrixXb & mat)
{
    Eigen::MatrixXi res = Eigen::MatrixXi::Constant(mat.rows(), mat.cols(), 0);
    for (int i = 0; i < mat.rows(); i++) {
        for (int j = 0; j < mat.cols(); j++) {
            if (mat(i, j)) { res(i, j) = 1; }
        }
    }
    return res;
}

inline Eigen::MatrixXb toBool(const Eigen::MatrixXi & mat)
{
    Eigen::MatrixXb res = Eigen::MatrixXb::Constant(mat.rows(), mat.cols(), false);
    for (int i = 0; i < mat.rows(); i++) {
        for (int j = 0; j < mat.cols(); j++) {
            if (mat(i, j) != 0) { res(i, j) = true; }
        }
    }
    return res;
}

inline int argmax(const Eigen::VectorXd & v)
{
    const int n = v.size();

    int i = 0;
    for (int j = 1; j < n; j++) {
        if (v(j) > v(i)) {
            i = j;
        }
    }
    return i;
}

inline int argmin(const Eigen::VectorXd & v)
{
    const int n = v.size();

    int i = 0;
    for (int j = 1; j < n; j++) {
        if (v(j) < v(i)) {
            i = j;
        }
    }
    return i;
}

#endif /* common_hpp */
