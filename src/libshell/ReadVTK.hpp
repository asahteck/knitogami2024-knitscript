//
//  ReadVTK.hpp
//  Elasticity
//
//  Created by Wim van Rees on 7/11/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef ReadVTK_h
#define ReadVTK_h

#include "common.hpp"

#ifdef USEVTK
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkXMLReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkDataSetReader.h>
#include <vtkDataSet.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkStringArray.h>
#include <vtkBitArray.h>
#include <vtkFieldData.h>
#include <vtkTriangle.h>
#include <vtkLine.h>
#include <vtkDataArrayAccessor.h>
#endif

class ReadVTK
{

public:

#ifdef USEVTK

    static void read(const std::string filename, Eigen::MatrixXd & verts, Eigen::MatrixXi & faces)
    {
        vtkDataSet *dataSet;

        // open file and assign pointer
        {
            vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();
            reader->GetOutput()->Register(reader);
            dataSet = vtkDataSet::SafeDownCast(reader->GetOutput());
        }

        const int nVertices = dataSet->GetNumberOfPoints();
        const int nFaces = dataSet->GetNumberOfCells();

        std::cout << "input file \t" << filename << "\t contains " << nVertices << " vertices and " << nFaces << " faces" << std::endl;
        verts.resize(nVertices,3);
        faces.resize(nFaces,3);

        // get the points
        for(int i=0;i<nVertices;++i)
        {
            double x[3];
            dataSet->GetPoint(i, x);
            for(int d=0;d<3;++d)
                verts(i,d) = x[d];
        }

        // get the face indices
        for(int i=0;i<nFaces;++i)
        {
            vtkSmartPointer<vtkIdList> ids_list = vtkSmartPointer<vtkIdList>::New();
            dataSet->GetCellPoints(i, ids_list);
            for(int d=0;d<3;++d)
                faces(i,d) = ids_list->GetId(d);
        }

        dataSet->Delete();
    }

    static void read(const std::string filename_in, Eigen::MatrixXd & verts, Eigen::MatrixXi & faces, Eigen::MatrixXd & attributes, std::vector<std::pair<std::string, int>> & attribute_entries)
    {
        vtkDataSet *dataSet;

        // check whether or not extension is already there
        const std::size_t found_vtp_ext = filename_in.rfind(".vtp");
        const std::string filename = (found_vtp_ext != std::string::npos ? filename_in.substr(0, found_vtp_ext) : filename_in);

        // open file and assign pointer
        {
            vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
            reader->SetFileName((filename+".vtp").c_str());
            reader->Update();
            reader->GetOutput()->Register(reader);
            dataSet = vtkDataSet::SafeDownCast(reader->GetOutput());
        }

        const int nVertices = dataSet->GetNumberOfPoints();
        const int nFaces = dataSet->GetNumberOfCells();

        std::cout << "input file \t" << filename << "\t contains " << nVertices << " vertices and " << nFaces << " faces" << std::endl;
        verts.resize(nVertices,3);
        faces.resize(nFaces,3);

        // get the points
        for(int i=0;i<nVertices;++i)
        {
            double x[3];
            dataSet->GetPoint(i, x);
            for(int d=0;d<3;++d)
                verts(i,d) = x[d];
        }

        // get the face indices
        for(int i=0;i<nFaces;++i)
        {
            vtkSmartPointer<vtkIdList> ids_list = vtkSmartPointer<vtkIdList>::New();
            dataSet->GetCellPoints(i, ids_list);
            for(int d=0;d<3;++d)
                faces(i,d) = ids_list->GetId(d);
        }

        // get the attributes
        vtkCellData * cd = dataSet->GetCellData();
        const int nArrays = cd->GetNumberOfArrays();
        if(nArrays>0)
        {
            attribute_entries.resize(nArrays);
            // count the number of columns we need in our matrix
            int nCols = 0;
            for(int i=0;i<nArrays;++i)
            {
                const std::string arrayName = cd->GetArrayName(i);
                vtkSmartPointer<vtkDataArray> data_array = cd->GetArray(arrayName.c_str());
                const int nTuples = data_array->GetNumberOfTuples();
                const int nComponents = data_array->GetNumberOfComponents();
                std::cout << "found array " << arrayName << " with " << nComponents << " components" << std::endl;
                if(nTuples!=nFaces)
                {
                    std::cout << "ReadVTK::read : problem with input file, number of faces = " << nFaces << " and number of tuples for array " << arrayName << " is " << nTuples << std::endl;
                    continue;
                }
                nCols += nComponents;
                attribute_entries[i] = std::make_pair(arrayName, nComponents);
            }

            // allocate memory
            attributes.resize(nFaces, nCols);

            // fill the array
            int colidx = 0;
            for(int i=0;i<nArrays;++i)
            {
                const std::string arrayName = cd->GetArrayName(i);
                vtkSmartPointer<vtkDataArray> data_array = cd->GetArray(arrayName.c_str());
                const int nTuples = data_array->GetNumberOfTuples();
                const int nComponents = data_array->GetNumberOfComponents();
                if(nTuples!=nFaces) continue;

                for(int t=0;t<nTuples;++t)
                {
                    double x[nComponents];
                    data_array->GetTuple(t, x);
                    for(int c=0;c<nComponents;++c)
                        attributes(t,colidx+c) = x[c];
                }
                colidx += nComponents;
            }
        }

        dataSet->Delete();
    }

#else
    static void read(const std::string , Eigen::MatrixXd & , Eigen::MatrixXi & )
    {
        std::cout << "ReadVTK::read() only implemented for compilation with VTK (flag USEVTK)\n";
    }

    static void read(const std::string , Eigen::MatrixXd & , Eigen::MatrixXi & , Eigen::MatrixXd & , std::vector<std::pair<std::string, int>> & )
    {
        std::cout << "ReadVTK::read() only implemented for compilation with VTK (flag USEVTK)\n";
    }
#endif
};


// class ReadVTKWithEdges
// {
// public:
//     static void read(const std::string filename_in,
//                      Eigen::MatrixXd & vertices,
//                      Eigen::MatrixXd & vertices_rest,
//                      Eigen::MatrixXb & vertices_bc,
//                      Eigen::MatrixXi & edges,
//                      Eigen::VectorXd & edgeDirectors,
//                      Eigen::VectorXd & edgeDirectors_rest,
//                      Eigen::VectorXb & edges_bc,
//                      Eigen::MatrixXi & faces,
//                      bool & bilayer)
//     {
//         vtkPolyData *dataSet;
//
//         // check whether or not extension is already there
//         const std::size_t found_vtp_ext = filename_in.rfind(".vtp");
//         const std::string filename = (found_vtp_ext != std::string::npos ? filename_in.substr(0, found_vtp_ext) : filename_in);
//
//         // open file and assign pointer
//         {
//             vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
//             reader->SetFileName((filename+".vtp").c_str());
//             reader->Update();
//             reader->GetOutput()->Register(reader);
//             dataSet = vtkPolyData::SafeDownCast(reader->GetOutput());
//         }
//
//         const int nVertices = dataSet->GetNumberOfPoints();
//         const int nEdges = dataSet->GetNumberOfLines();
//         const int nFaces = dataSet->GetNumberOfCells() - nEdges;
//
//         // get the edges and faces
//         std::cout << "input file \t" << filename << "\t contains "
//                   << nVertices << " vertices and "
//                   << nEdges << " edges and "
//                   << nFaces << " faces." << std::endl;
//
//         vertices.resize(nVertices, 3);
//         edges.resize(nEdges, 2);
//         faces.resize(nFaces, 3);
//
//         bilayer = false; // TODO
//
//         // get the current vertices
//         for(int i=0;i<nVertices;++i)
//         {
//             double p[3];
//             dataSet->GetPoint(i, p);
//             for(int d=0;d<3;++d) {
//                 vertices(i, d) = p[d];
//             }
//         }
//
//         // get the edge indices
//         for(int i=0;i<nEdges;++i)
//         {
//             vtkSmartPointer<vtkIdList> ids_list = vtkSmartPointer<vtkIdList>::New();
//             dataSet->GetCellPoints(i, ids_list);
//             for(int d=0;d<2;++d)
//                 edges(i,d) = ids_list->GetId(d);
//         }
//
//         // get the face indices
//         for(int i=0;i<nFaces;++i)
//         {
//             vtkSmartPointer<vtkIdList> ids_list = vtkSmartPointer<vtkIdList>::New();
//             dataSet->GetCellPoints(nEdges + i, ids_list);
//             for(int d=0;d<3;++d)
//                 faces(i,d) = ids_list->GetId(d);
//         }
//
//         // get the cell data (directors, rest directors, edges_bc)
//         vtkCellData * cd = dataSet->GetCellData();
//         vtkSmartPointer<vtkAbstractArray> array;
//
//         array = cd->GetArray("currentDirectors");
//         if (array) {
//             vtkSmartPointer<vtkDoubleArray> data_array = vtkDoubleArray::SafeDownCast(array);
//             if (data_array) {
//                 assert(data_array->GetNumberOfTuples() == nEdges);
//                 assert(data_array->GetNumberOfComponents() == 1);
//
//                 Real x;
//                 for (int i = 0; i < nEdges; i++) {
//                     data_array->GetTuple(i, &x);
//                     edgeDirectors(i) = x;
//                 }
//                 std::cout << "read currentDirectors" << std::endl;
//             }
//         }
//
//         array = cd->GetArray("restDirectors");
//         if (array) {
//             vtkSmartPointer<vtkDoubleArray> data_array = vtkDoubleArray::SafeDownCast(array);
//             if (data_array) {
//                 assert(data_array->GetNumberOfTuples() == nEdges);
//                 assert(data_array->GetNumberOfComponents() == 1);
//
//                 Real x;
//                 for (int i = 0; i < nEdges; i++) {
//                     data_array->GetTuple(i, &x);
//                     edgeDirectors_rest(i) = x;
//                 }
//                 std::cout << "read restDirectors" << std::endl;
//             }
//         }
//
//         array = cd->GetArray("edges_bc");
//         if (array) {
//             vtkSmartPointer<vtkIntArray> data_array = vtkIntArray::SafeDownCast(array);
//             if (data_array) {
//                 assert(data_array->GetNumberOfTuples() == nEdges);
//                 assert(data_array->GetNumberOfComponents() == 1);
//
//                 int x;
//                 for (int i = 0; i < nEdges; i++) {
//                     data_array->GetTypedTuple(i, &x);
//                     edges_bc(i) = x;
//                 }
//                 std::cout << "read edges_bc" << std::endl;
//             }
//         }
//
//         // get the point data
//         vtkPointData * pd = dataSet->GetPointData();
//
//         array = pd->GetArray("vertices_rest");
//         if (array) {
//             vtkSmartPointer<vtkDoubleArray> data_array = vtkDoubleArray::SafeDownCast(array);
//             if (data_array) {
//                 assert(data_array->GetNumberOfTuples() == nVertices);
//                 assert(data_array->GetNumberOfComponents() == 3);
//
//                 Real x[3];
//                 for (int i = 0; i < nEdges; i++) {
//                     data_array->GetTuple(i, x);
//                     vertices_rest(i, 0) = x[0];
//                     vertices_rest(i, 1) = x[1];
//                     vertices_rest(i, 2) = x[2];
//                 }
//                 std::cout << "read vertices_rest" << std::endl;
//             }
//         }
//
//         array = pd->GetArray("vertices_bc");
//         if (array) {
//             vtkSmartPointer<vtkIntArray> data_array = vtkIntArray::SafeDownCast(array);
//             if (data_array) {
//                 assert(data_array->GetNumberOfTuples() == nVertices);
//                 assert(data_array->GetNumberOfComponents() == 3);
//
//                 int x[3];
//                 for (int i = 0; i < nEdges; i++) {
//                     data_array->GetTypedTuple(i, x);
//                     vertices_bc(i, 0) = x[0];
//                     vertices_bc(i, 1) = x[1];
//                     vertices_bc(i, 2) = x[2];
//                 }
//                 std::cout << "read vertices_bc" << std::endl;
//             }
//         }
//
//         dataSet->Delete();
//     }
// };


#endif /* ReadVTK_h */
