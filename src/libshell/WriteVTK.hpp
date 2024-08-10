//
//  WriteVTK.hpp
//  Elasticity
//
//  Created by Wim van Rees on 2/17/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef WriteVTK_hpp
#define WriteVTK_hpp

#include "common.hpp"

#ifdef USEVTK
#include <vtkVersion.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkStringArray.h>
#include <vtkBitArray.h>
#include <vtkFieldData.h>
#include <vtkTriangle.h>
#include <vtkLine.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#endif

/*! \class WriteVTK
 * \brief Write a triangle mesh in VTK output to be read by Paraview or Visit
 *
 * If USEVTK is defined, this method will use VTK libraries to write output in binary format using VTP (polygonal) format. If USEVTK is not defined, it will write plain ASCII using inefficient legacy VTK format.
 */

#ifdef USEVTK
class WriteVTK
{
protected:
    int nVertices;
    int nFaces;
    vtkSmartPointer<vtkPolyData> polydata;

public:
    WriteVTK(const Eigen::Ref<const Eigen::MatrixXd> verts, const Eigen::Ref<const Eigen::MatrixXi> faces):
    nVertices(verts.rows()),
    nFaces(faces.rows())
    {
        assert(nVertices > 0);
        assert(nFaces > 0);
        assert(verts.cols() == 3);
        assert(faces.cols() == 3);

        // create points
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        points->SetNumberOfPoints(nVertices);
        for(int i=0;i<nVertices;++i)
        {
            points->SetPoint(i, verts(i,0), verts(i,1), verts(i,2) );
        }

        // Create our polydata object and add the points to it.
        polydata = vtkSmartPointer<vtkPolyData>::New();
        polydata->SetPoints(points);

        // create faces
        vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
        for(int i=0;i<nFaces;++i)
        {
            // Create a triangle
            vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();

            triangle->GetPointIds()->SetId(0,faces(i,0));
            triangle->GetPointIds()->SetId(1,faces(i,1));
            triangle->GetPointIds()->SetId(2,faces(i,2));

            triangles->InsertNextCell(triangle);
        }

        polydata->SetPolys(triangles);

        // Container for extra attributes
        polydata->SetFieldData(vtkSmartPointer<vtkFieldData>::New());
    }

    WriteVTK(const Eigen::Ref<const Eigen::MatrixXd> verts):
    nVertices(verts.rows()),
    nFaces(-1)
    {
        assert(nVertices > 0);
        assert(verts.cols() == 3);

        // create points
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        points->SetNumberOfPoints(nVertices);
        for(int i=0;i<nVertices;++i)
        {
            points->SetPoint(i, verts(i,0), verts(i,1), verts(i,2) );
        }

        // Create our polydata object and add the points to it.
        polydata = vtkSmartPointer<vtkPolyData>::New();
        polydata->SetPoints(points);

        // No faces

        // Container for extra attributes
        polydata->SetFieldData(vtkSmartPointer<vtkFieldData>::New());
    }

    WriteVTK():
    nVertices(-1),
    nFaces(-1)
    {
        // Create our polydata object.
        polydata = vtkSmartPointer<vtkPolyData>::New();

        // Container for extra attributes
        polydata->SetFieldData(vtkSmartPointer<vtkFieldData>::New());
    }

    bool initialized()
    {
        return nVertices >= 0;
    }

    void reset()
    {
        nVertices = -1;
        nFaces = -1;
        polydata = nullptr; // calls destructor

        // Create our polydata object.
        polydata = vtkSmartPointer<vtkPolyData>::New();

        // Container for extra attributes
        polydata->SetFieldData(vtkSmartPointer<vtkFieldData>::New());
    }

    void updateMesh(const Eigen::Ref<const Eigen::MatrixXd> verts, const Eigen::Ref<const Eigen::MatrixXi> faces)
    {
        assert(verts.cols() == 3);
        assert(faces.cols() == 3);

        nVertices = verts.rows();
        nFaces = faces.rows();

        // create points
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        points->SetNumberOfPoints(nVertices);
        for(int i=0;i<nVertices;++i)
        {
            points->SetPoint(i, verts(i,0), verts(i,1), verts(i,2) );
        }

        polydata->SetPoints(points);

        // create faces
        vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
        for(int i=0;i<nFaces;++i)
        {
            // Create a triangle
            vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();

            triangle->GetPointIds()->SetId(0,faces(i,0));
            triangle->GetPointIds()->SetId(1,faces(i,1));
            triangle->GetPointIds()->SetId(2,faces(i,2));

            triangles->InsertNextCell(triangle);
        }

        polydata->SetPolys(triangles);
    }

    vtkSmartPointer<vtkDoubleArray> convertVectorFieldToVTK(const Eigen::Ref<const Eigen::MatrixXd> field, const std::string fieldname)
    {
        const int N = field.rows();
        assert(field.cols() == 3);
        vtkSmartPointer<vtkDoubleArray> vector_array = vtkSmartPointer<vtkDoubleArray>::New();
        vector_array->SetNumberOfComponents(3);
        vector_array->SetNumberOfTuples(N);
        vector_array->SetName(fieldname.c_str());
        for(int i=0;i<N;++i)
            vector_array->SetTuple3(i, field(i,0),field(i,1),field(i,2));

        return vector_array;
    }

    vtkSmartPointer<vtkDoubleArray> convertScalarFieldToVTK(const Eigen::Ref<const Eigen::VectorXd> field, const std::string fieldname)
    {
        const int N = field.rows();
        vtkSmartPointer<vtkDoubleArray> scalar_array = vtkSmartPointer<vtkDoubleArray>::New();
        scalar_array->SetNumberOfComponents(1);
        scalar_array->SetNumberOfTuples(N);
        scalar_array->SetName(fieldname.c_str());
        for(int i=0;i<N;++i)
            scalar_array->SetTuple1(i, field(i));

        return scalar_array;
    }

    void addTensorFieldToFaces(const tVecMat2d & field, const std::string fieldname)
    {
        const int N = (int)field.size();
        assert(N == nFaces);
        // assuming the tensor is symmetric

        vtkSmartPointer<vtkDoubleArray> tensor_array = vtkSmartPointer<vtkDoubleArray>::New();
        tensor_array->SetNumberOfComponents(4);
        tensor_array->SetNumberOfTuples(N);
        tensor_array->SetName(fieldname.c_str());
        for(int i=0;i<N;++i)
        {
            const Eigen::Matrix2d & tensor = field[i];
            tensor_array->SetTuple4(i, tensor(0,0), tensor(0,1), tensor(1,0), tensor(1,1));
        }
        polydata->GetCellData()->AddArray(tensor_array);
    }

    void addTensorFieldToFaces(const tVecMat3d & field, const std::string fieldname)
    {
        const int N = (int)field.size();
        assert(N == nFaces);
        // assuming the tensor is symmetric

        vtkSmartPointer<vtkDoubleArray> tensor_array = vtkSmartPointer<vtkDoubleArray>::New();
        tensor_array->SetNumberOfComponents(9);
        tensor_array->SetNumberOfTuples(N);
        tensor_array->SetName(fieldname.c_str());
        for(int i=0;i<N;++i)
        {
            const Eigen::Matrix3d & tensor = field[i];
            tensor_array->SetTuple9(i, tensor(0,0), tensor(0,1), tensor(0,2), tensor(1,0), tensor(1,1), tensor(1,2), tensor(2,0), tensor(2,1), tensor(2,2));
        }
        polydata->GetCellData()->AddArray(tensor_array);
    }

    void addVectorFieldToFaces(const Eigen::Ref<const Eigen::MatrixXd> field, const std::string fieldname)
    {
        assert(field.rows() == nFaces);
        assert(field.cols() == 3);

        vtkSmartPointer<vtkDoubleArray> vector_array = convertVectorFieldToVTK(field, fieldname);
        polydata->GetCellData()->AddArray(vector_array);
    }

    void addScalarFieldToFaces(const Eigen::Ref<const Eigen::VectorXd> field, const std::string fieldname)
    {
        assert(field.rows() == nFaces);

        vtkSmartPointer<vtkDoubleArray> scalar_array = convertScalarFieldToVTK(field, fieldname);
        polydata->GetCellData()->AddArray(scalar_array);
    }

    void addVectorFieldToVertices(const Eigen::Ref<const Eigen::MatrixXd> field, const std::string fieldname)
    {
        assert(field.rows() == nVertices);
        assert(field.cols() == 3);

        vtkSmartPointer<vtkDoubleArray> vector_array = convertVectorFieldToVTK(field, fieldname);
        polydata->GetPointData()->AddArray(vector_array);
    }

    void addScalarFieldToVertices(const Eigen::Ref<const Eigen::VectorXd> field, const std::string fieldname)
    {
        assert(field.rows() == nVertices);

        vtkSmartPointer<vtkDoubleArray> scalar_array = convertScalarFieldToVTK(field, fieldname);
        polydata->GetPointData()->AddArray(scalar_array);
    }

    void addScalarAttribute(const Real & scalar, const std::string fieldname)
    {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfComponents(1);
        data->SetNumberOfTuples(1);
        data->SetName(fieldname.c_str());
        data->SetTuple1(0, scalar);

        vtkSmartPointer<vtkFieldData> field = polydata->GetFieldData();
        field->AddArray(data);
    }

    void addScalarAttribute(const int & scalar, const std::string fieldname)
    {
        vtkSmartPointer<vtkIntArray> data = vtkSmartPointer<vtkIntArray>::New();
        data->SetNumberOfComponents(1);
        data->SetNumberOfTuples(1);
        data->SetName(fieldname.c_str());
        data->SetTuple1(0, scalar);

        vtkSmartPointer<vtkFieldData> field = polydata->GetFieldData();
        field->AddArray(data);
    }

    void addStringAttribute(const std::string str, const std::string fieldname)
    {
        vtkSmartPointer<vtkStringArray> data = vtkSmartPointer<vtkStringArray>::New();
        data->SetNumberOfValues(1);
        data->SetName(fieldname.c_str());
        data->SetValue(0, str.c_str());

        vtkSmartPointer<vtkFieldData> field = polydata->GetFieldData();
        field->AddArray(data);
    }

    int write(const std::string filename)
    {
        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
#if VTK_MAJOR_VERSION <= 5
        writer->SetInput(polydata);
#else
        writer->SetInputData(polydata);
#endif
        writer->SetFileName( (filename+".vtp").c_str());
        const int retval = writer->Write(); // return value 1 for success, 0 for failure

        return retval;
    }
};


// class WriteVTKWithEdges
// {
// protected:
//     bool rest = true;
//     int nVertices = -1;
//     int nFaces = -1;
//     int nEdges = -1;
//
//     vtkSmartPointer<vtkPolyData> polydata;
//
//     // Helper structs for dealing with VTK types
//     template<typename T>
//     struct VTKArrayType {};
//
//     template<>
//     struct VTKArrayType<int> { typedef vtkIntArray type; };
//
//     template<>
//     struct VTKArrayType<Real> { typedef vtkDoubleArray type; };
//
//     template<>
//     struct VTKArrayType<bool> { typedef vtkBitArray type; };
//
// public:
//     enum DataType { VERTEX, FACE, EDGE };
//
//     WriteVTKWithEdges(const Eigen::Ref<const Eigen::MatrixXd> verts, const Eigen::Ref<const Eigen::MatrixXi> edges, const Eigen::Ref<const Eigen::MatrixXi> faces):
//     nVertices(verts.rows()),
//     nEdges(edges.rows()),
//     nFaces(faces.rows())
//     {
//         assert(nVertices > 0);
//         assert(nEdges > 0);
//         assert(nFaces > 0);
//         assert(verts.cols() == 3);
//         assert(edges.cols() == 2);
//         assert(faces.cols() == 3);
//
//         // Create our polydata object and add the points to it.
//         polydata = vtkSmartPointer<vtkPolyData>::New();
//
//         // create points
//         vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
//         points->SetNumberOfPoints(nVertices);
//         for(int i=0;i<nVertices;++i)
//         {
//             points->SetPoint(i, verts(i,0), verts(i,1), verts(i,2) );
//         }
//         polydata->SetPoints(points);
//
//         // create edges/lines
//         vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
//         for(int i=0;i<nEdges;++i)
//         {
//             // Create a line
//             vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
//
//             line->GetPointIds()->SetId(0,edges(i,0));
//             line->GetPointIds()->SetId(1,edges(i,1));
//
//             lines->InsertNextCell(line);
//         }
//
//         polydata->SetLines(lines);
//
//         // create faces
//         vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
//         for(int i=0;i<nFaces;++i)
//         {
//             // Create a triangle
//             vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
//
//             triangle->GetPointIds()->SetId(0,faces(i,0));
//             triangle->GetPointIds()->SetId(1,faces(i,1));
//             triangle->GetPointIds()->SetId(2,faces(i,2));
//
//             triangles->InsertNextCell(triangle);
//         }
//
//         polydata->SetPolys(triangles);
//
//         // Container for extra attributes
//         polydata->SetFieldData(vtkSmartPointer<vtkFieldData>::New());
//     }
//
//     bool initialized()
//     {
//         return nVertices >= 0;
//     }
//
//     void reset()
//     {
//         nVertices = -1;
//         nEdges = -1;
//         nFaces = -1;
//         polydata = nullptr; // calls destructor
//
//         // Create our polydata object.
//         polydata = vtkSmartPointer<vtkPolyData>::New();
//
//         // Container for extra attributes
//         polydata->SetFieldData(vtkSmartPointer<vtkFieldData>::New());
//     }
//
//     void updateMesh(const Eigen::Ref<const Eigen::MatrixXd> verts, const Eigen::Ref<const Eigen::MatrixXi> edges, const Eigen::Ref<const Eigen::MatrixXi> faces)
//     {
//         assert(verts.cols() == 3);
//         assert(edges.cols() == 2);
//         assert(faces.cols() == 3);
//
//         nVertices = verts.rows();
//         nEdges = edges.rows();
//         nFaces = faces.rows();
//
//         // create points
//         vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
//         points->SetNumberOfPoints(nVertices);
//         for(int i=0;i<nVertices;++i)
//         {
//             points->SetPoint(i, verts(i,0), verts(i,1), verts(i,2) );
//         }
//
//         polydata->SetPoints(points);
//
//         // create edges/lines
//         vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
//         for(int i=0;i<nEdges;++i)
//         {
//             // Create a line
//             vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
//
//             line->GetPointIds()->SetId(0,edges(i,0));
//             line->GetPointIds()->SetId(1,edges(i,1));
//
//             lines->InsertNextCell(line);
//         }
//
//         polydata->SetLines(lines);
//
//         // create faces
//         vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
//         for(int i=0;i<nFaces;++i)
//         {
//             // Create a triangle
//             vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
//
//             triangle->GetPointIds()->SetId(0,faces(i,0));
//             triangle->GetPointIds()->SetId(1,faces(i,1));
//             triangle->GetPointIds()->SetId(2,faces(i,2));
//
//             triangles->InsertNextCell(triangle);
//         }
//
//         polydata->SetPolys(triangles);
//     }
//
//     vtkSmartPointer<vtkDoubleArray> convertTensorFieldToVTK(const tVecMat2d field, const std::string fieldname, const DataType datatype)
//     {
//         const int N = (int)field.size();
//         vtkSmartPointer<vtkDoubleArray> vector_array = vtkSmartPointer<vtkDoubleArray>::New();
//         vector_array->SetNumberOfComponents(3);
//         vector_array->SetNumberOfTuples(datatype == VERTEX ? nVertices : nEdges + nFaces);
//         vector_array->SetName(fieldname.c_str());
//
//         if (datatype == FACE) {
//             assert(N == nFaces);
//             for (int i=0;i<nEdges;++i) {
//                 vector_array->SetTuple3(i, 0, 0, 0);
//             }
//             for (int i=0;i<N;++i) {
//                 vector_array->SetTuple3(nEdges + i, field[i](0, 0), field[i](0, 1), field[i](1, 1));
//             }
//         }
//         else if (datatype == EDGE) {
//             assert(N == nEdges);
//             for (int i=0;i<N;++i) {
//                 vector_array->SetTuple3(i, field[i](0, 0), field[i](0, 1), field[i](1, 1));
//             }
//             for (int i=0;i<nFaces;++i) {
//                 vector_array->SetTuple3(nEdges + i, 0, 0, 0);
//             }
//         }
//         else if (datatype == VERTEX) {
//             assert(N == nVertices);
//             for (int i=0;i<N;++i) {
//                 vector_array->SetTuple3(i, field[i](0, 0), field[i](0, 1), field[i](1, 1));
//             }
//         }
//         else {
//             helpers::catastrophe("Datatype of mesh data must be one of { VERTEX, FACE, EDGE }.", __FILE__, __LINE__);
//         }
//
//         return vector_array;
//     }
//
//     template <typename Derived>
//     vtkSmartPointer<typename VTKArrayType<typename Derived::Scalar>::type> convertFieldToVTK(const Eigen::MatrixBase<Derived> & field, const std::string fieldname, const DataType datatype=VERTEX)
//     {
//         const int N = field.rows();
//
//         assert(N == (datatype == VERTEX ? nVertices : (datatype == EDGE ? nEdges : nFaces)));
//
//         vtkSmartPointer<typename VTKArrayType<typename Derived::Scalar>::type> vector_array = vtkSmartPointer<typename VTKArrayType<typename Derived::Scalar>::type>::New();
//         vector_array->SetNumberOfComponents(field.cols());
//         vector_array->SetNumberOfTuples(datatype == VERTEX ? nVertices : nEdges + nFaces);
//         vector_array->SetName(fieldname.c_str());
//
//         if (field.cols() == 1) { // Scalar data
//             if (datatype == FACE) {
//                 for (int i=0;i<nEdges;++i) {
//                     vector_array->SetTuple1(i, 0);
//                 }
//                 for (int i=0;i<N;++i) {
//                     vector_array->SetTuple1(nEdges + i, field(i, 0));
//                 }
//             }
//             else if (datatype == EDGE) {
//                 for (int i=0;i<N;++i) {
//                     vector_array->SetTuple1(i, field(i, 0));
//                 }
//                 for (int i=0;i<nFaces;++i) {
//                     vector_array->SetTuple1(nEdges + i, 0);
//                 }
//             }
//             else if (datatype == VERTEX) {
//                 for (int i=0;i<N;++i) {
//                     vector_array->SetTuple1(i, field(i, 0));
//                 }
//             }
//             else {
//                 helpers::catastrophe("Datatype of mesh data must be one of { VERTEX, FACE, EDGE }.", __FILE__, __LINE__);
//             }
//         }
//         else if (field.cols() == 3) { // Vector data
//             if (datatype == FACE) {
//                 for (int i=0;i<nEdges;++i) {
//                     vector_array->SetTuple3(i, 0, 0, 0);
//                 }
//                 for (int i=0;i<N;++i) {
//                     vector_array->SetTuple3(nEdges + i, field(i, 0), field(i, 1), field(i, 2));
//                 }
//             }
//             else if (datatype == EDGE) {
//                 for (int i=0;i<N;++i) {
//                     vector_array->SetTuple3(i, field(i, 0), field(i, 1), field(i, 2));
//                 }
//                 for (int i=0;i<nFaces;++i) {
//                     vector_array->SetTuple3(nEdges + i, 0, 0, 0);
//                 }
//             }
//             else if (datatype == VERTEX) {
//                 for (int i=0;i<N;++i) {
//                     vector_array->SetTuple3(i, field(i, 0), field(i, 1), field(i, 2));
//                 }
//             }
//             else {
//                 helpers::catastrophe("Datatype of mesh data must be one of { VERTEX, FACE, EDGE }.", __FILE__, __LINE__);
//             }
//         }
//         else {
//             helpers::catastrophe("Only fields with 1 and 3 components are supported.", __FILE__, __LINE__);
//         }
//
//         return vector_array;
//     }
//
//     void addTensorField(const tVecMat2d & field, const std::string fieldname, const DataType datatype)
//     {
//         vtkSmartPointer<vtkDoubleArray> vector_array = convertTensorFieldToVTK(field, fieldname, datatype);
//         polydata->GetCellData()->AddArray(vector_array);
//     }
//
//     template <typename Derived>
//     void addField(const Eigen::MatrixBase<Derived> & field, const std::string fieldname, const DataType datatype)
//     {
//         vtkSmartPointer<typename VTKArrayType<typename Derived::Scalar>::type> vector_array = convertFieldToVTK(field, fieldname, datatype);
//         polydata->GetCellData()->AddArray(vector_array);
//     }
//
//     void addScalarAttribute(const Real & scalar, const std::string fieldname)
//     {
//         vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
//         data->SetNumberOfComponents(1);
//         data->SetNumberOfTuples(1);
//         data->SetName(fieldname.c_str());
//         data->SetTuple1(0, scalar);
//
//         vtkSmartPointer<vtkFieldData> field = polydata->GetFieldData();
//         field->AddArray(data);
//     }
//
//     void addScalarAttribute(const int & scalar, const std::string fieldname)
//     {
//         vtkSmartPointer<vtkIntArray> data = vtkSmartPointer<vtkIntArray>::New();
//         data->SetNumberOfComponents(1);
//         data->SetNumberOfTuples(1);
//         data->SetName(fieldname.c_str());
//         data->SetTuple1(0, scalar);
//
//         vtkSmartPointer<vtkFieldData> field = polydata->GetFieldData();
//         field->AddArray(data);
//     }
//
//     void addStringAttribute(const std::string str, const std::string fieldname)
//     {
//         vtkSmartPointer<vtkStringArray> data = vtkSmartPointer<vtkStringArray>::New();
//         data->SetNumberOfValues(1);
//         data->SetName(fieldname.c_str());
//         data->SetValue(0, str.c_str());
//
//         vtkSmartPointer<vtkFieldData> field = polydata->GetFieldData();
//         field->AddArray(data);
//     }
//
//     int write(const std::string filename)
//     {
//         vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
// #if VTK_MAJOR_VERSION <= 5
//         writer->SetInput(polydata);
// #else
//         writer->SetInputData(polydata);
// #endif
//         writer->SetFileName( (filename+".vtp").c_str());
//         const int retval = writer->Write(); // return value 1 for success, 0 for failure
//
//         return retval;
//     }
// };

#else /* USEVTK */

    // ASCII! not for larger domains

class WriteVTK
{
protected:
    const int nVertices;
    const int nFaces;
    const Eigen::Ref<const Eigen::MatrixXd> verts;
    const Eigen::Ref<const Eigen::MatrixXi> faces;

public:
    WriteVTK(const Eigen::Ref<const Eigen::MatrixXd> verts, const Eigen::Ref<const Eigen::MatrixXi> faces):
    nVertices(verts.rows()),
    nFaces(faces.rows()),
    verts(verts),
    faces(faces)
    {
        assert(nVertices > 0);
        assert(nFaces > 0);
        assert(verts.cols() == 3);
        assert(faces.cols() == 3);

    }

    void addVectorFieldToFaces(const Eigen::Ref<const Eigen::MatrixXd> , const std::string )
    {
        std::cout << "addVectorFieldToFaces not implemented for WriteVTK if not compiled with VTK support\n";
    }

    void addScalarFieldToFaces(const Eigen::Ref<const Eigen::VectorXd> , const std::string )
    {
        std::cout << "addScalarFieldToFaces not implemented for WriteVTK if not compiled with VTK support\n";
    }

    void addVectorFieldToVertices(const Eigen::Ref<const Eigen::MatrixXd> , const std::string )
    {
        std::cout << "addVectorFieldToVertices not implemented for WriteVTK if not compiled with VTK support\n";
    }

    void addScalarFieldToVertices(const Eigen::Ref<const Eigen::VectorXd> , const std::string )
    {
        std::cout << "addScalarFieldToVertices not implemented for WriteVTK if not compiled with VTK support\n";
    }

    void addScalarAttribute(const Real &, const std::string)
    {
        std::cout << "addScalarAttribute not implemented for WriteVTK if not compiled with VTK support\n";
    }

    int write(const std::string filename)
    {
        std::ofstream ofs(filename+".vtk");
        if (not ofs.is_open()) return 0;

        ofs << "# vtk DataFile Version 1.0" << std::endl;
        ofs << "2D Unstructured Grid of Linear Triangles" << std::endl;
        ofs << "ASCII" << std::endl;
        ofs << std::endl;
        ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;
        ofs << "POINTS "+std::to_string(nVertices)+" float" << std::endl;

        for(int i=0;i<nVertices;++i)
            ofs << verts(i,0) << " " << verts(i,1) << " " << verts(i,2) << std::endl;

        ofs << std::endl;
        ofs << "CELLS "+std::to_string(nFaces)+" "+std::to_string(nFaces*4) << std::endl;
        for(int i=0;i<nFaces;++i)
            ofs << 3 << " " << faces(i,0) << " " << faces(i,1) << " " << faces(i,2) << std::endl;

        ofs<<std::endl;
        ofs << "CELL_TYPES "+std::to_string(nFaces) << std::endl;
        for(int i=0;i<nFaces;++i)
            ofs << "5" << std::endl;

        ofs.close();

        return 1;
    }
};
#endif /* USEVTK */

#endif /* WriteVTK_hpp */
