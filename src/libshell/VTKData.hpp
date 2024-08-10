#ifndef VTKData_hpp
#define VTKData_hpp

#include "common.hpp"

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
#include <vtkXMLPolyDataReader.h>


// Helper structs for dealing with VTK types
#ifndef VTK_ARRAY_TYPES
#define VTK_ARRAY_TYPES

template<typename T>
struct VTKArrayType {};

template<>
struct VTKArrayType<int> { typedef vtkIntArray type; };

template<>
struct VTKArrayType<Real> { typedef vtkDoubleArray type; };

template<>
struct VTKArrayType<bool> { typedef vtkBitArray type; };

template<>
struct VTKArrayType<std::string> { typedef vtkStringArray type; };

#endif


class VTKData
{
protected:
    int nVertices = 0;
    int nEdges = 0;
    int nFaces = 0;

    vtkSmartPointer<vtkPolyData> polydata;

    template <typename Derived>
    vtkSmartPointer<typename VTKArrayType<typename Derived::Scalar>::type> convertFieldToVTK(const Eigen::MatrixBase<Derived> & field, const std::string fieldname, const std::string datatype)
    {
        if (nEdges == 0 && datatype == "edge") {
            std::cout << "Edge data not supported for this instance of VTKData because there are no edges." << std::endl;
            exit(1);
        }

        const int N = field.rows();
        assert(N == (datatype == "vertex" ? nVertices : (datatype == "edge" ? nEdges : nFaces)));

        vtkSmartPointer<typename VTKArrayType<typename Derived::Scalar>::type> vector_array = vtkSmartPointer<typename VTKArrayType<typename Derived::Scalar>::type>::New();
        vector_array->SetNumberOfComponents(field.cols());
        vector_array->SetNumberOfTuples(datatype == "vertex" ? nVertices : (nEdges > 0 ? nEdges + nFaces : (datatype == "edge" ? nEdges : nFaces)));
        vector_array->SetName(fieldname.c_str());

        if (field.cols() == 1) { // Scalar data
            if (datatype == "face") {
                for (int i=0;i<nEdges;++i) {
                    vector_array->SetTuple1(i, 0);
                }
                for (int i=0;i<N;++i) {
                    vector_array->SetTuple1(nEdges + i, field(i, 0));
                }
            }
            else if (datatype == "edge") { //
                for (int i=0;i<N;++i) {
                    vector_array->SetTuple1(i, field(i, 0));
                }
                for (int i=0;i<nFaces;++i) {
                    vector_array->SetTuple1(nEdges + i, 0);
                }
            }
            else if (datatype == "vertex") {
                for (int i=0;i<N;++i) {
                    vector_array->SetTuple1(i, field(i, 0));
                }
            }
            else {
                helpers::catastrophe("Datatype of mesh data must be one of { 'vertex', 'face', 'edge' }.", __FILE__, __LINE__);
            }
        }
        else if (field.cols() == 3) { // Vector data
            if (datatype == "face") {
                for (int i=0;i<nEdges;++i) {
                    vector_array->SetTuple3(i, 0, 0, 0);
                }
                for (int i=0;i<N;++i) {
                    vector_array->SetTuple3(nEdges + i, field(i, 0), field(i, 1), field(i, 2));
                }
            }
            else if (datatype == "edge") {
                for (int i=0;i<N;++i) {
                    vector_array->SetTuple3(i, field(i, 0), field(i, 1), field(i, 2));
                }
                for (int i=0;i<nFaces;++i) {
                    vector_array->SetTuple3(nEdges + i, 0, 0, 0);
                }
            }
            else if (datatype == "vertex") {
                for (int i=0;i<N;++i) {
                    vector_array->SetTuple3(i, field(i, 0), field(i, 1), field(i, 2));
                }
            }
            else {
                helpers::catastrophe("Datatype of mesh data must be one of { 'vertex', 'face', 'edge' }.", __FILE__, __LINE__);
            }
        }
        else {
            helpers::catastrophe("Only fields with 1 and 3 components are supported.", __FILE__, __LINE__);
        }

        return vector_array;
    }

public:
    VTKData()
    {
        clear();
    }

    VTKData(const std::string filename_in)
    {
        load(filename_in);
    }

    void load(const std::string filename_in)
    {
        clear();

        // check whether or not extension is already there
        const std::size_t found_vtp_ext = filename_in.rfind(".vtp");
        const std::string filename = (found_vtp_ext != std::string::npos ? filename_in.substr(0, found_vtp_ext) : filename_in);

        // open file and assign pointer
        vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        reader->SetFileName((filename+".vtp").c_str());
        reader->Update();
        reader->GetOutput()->Register(reader);
        polydata = vtkPolyData::SafeDownCast(reader->GetOutput());

        nVertices = polydata->GetNumberOfPoints();
        nEdges = polydata->GetNumberOfLines();
        nFaces = polydata->GetNumberOfCells() - nEdges;
    }

    void clear()
    {
        nVertices = 0;
        nEdges = 0;
        nFaces = 0;

        polydata = nullptr; // calls destructor

        // Create our polydata object.
        polydata = vtkSmartPointer<vtkPolyData>::New();

        // Container for extra attributes
        polydata->SetFieldData(vtkSmartPointer<vtkFieldData>::New());
    }

    bool any(const std::list<std::string> & fieldnames)
    {
        return std::any_of(fieldnames.begin(), fieldnames.end(), [this](const std::string f){ return has(f); });
    }

    bool all(const std::list<std::string> & fieldnames)
    {
        return std::all_of(fieldnames.begin(), fieldnames.end(), [this](const std::string f){ return has(f); });
    }

    bool has(const std::string fieldname)
    {
        return isAttribute(fieldname) || isField(fieldname);
    }

    bool isAttribute(const std::string fieldname)
    {
        vtkSmartPointer<vtkAbstractArray> arr = polydata->GetFieldData()->GetAbstractArray(fieldname.c_str());
        if (!arr) {
            return false;
        }
        return true;
    }

    bool isField(const std::string fieldname)
    {
        vtkSmartPointer<vtkAbstractArray> arr = polydata->GetCellData()->GetAbstractArray(fieldname.c_str());

        if (!arr) {
            arr = polydata->GetPointData()->GetAbstractArray(fieldname.c_str());
            if (!arr) {
                return false;
            }
        }

        return true;
    }

    template <typename ScalarT>
    ScalarT getAttribute(const std::string fieldname, const ScalarT & defval)
    {
        vtkSmartPointer<vtkAbstractArray> arr = polydata->GetFieldData()->GetAbstractArray(fieldname.c_str());
        if (!arr) {
            return defval;
        }

        assert(arr->GetNumberOfTuples() == 1);
        assert(arr->GetNumberOfComponents() == 1);

        ScalarT field;

        if (arr->GetDataType() == 1) { // bit data
            vtkSmartPointer<vtkBitArray> arr2 = vtkBitArray::SafeDownCast(arr);
            field = arr2->GetValue(0);
        }
        else {
            vtkSmartPointer<typename VTKArrayType<ScalarT>::type> arr2 = VTKArrayType<ScalarT>::type::SafeDownCast(arr);
            field = arr2->GetValue(0);
        }

        return field;
    }

    template <typename ScalarT>
    void setScalarAttribute(const ScalarT scalar, const std::string fieldname)
    {
        vtkSmartPointer<typename VTKArrayType<ScalarT>::type> data = vtkSmartPointer<typename VTKArrayType<ScalarT>::type>::New();
        data->SetNumberOfComponents(1);
        data->SetNumberOfTuples(1);
        data->SetName(fieldname.c_str());
        data->SetTuple1(0, scalar);

        vtkSmartPointer<vtkFieldData> field = polydata->GetFieldData();
        field->AddArray(data);
    }

    void setStringAttribute(const std::string str, const std::string fieldname)
    {
        vtkSmartPointer<vtkStringArray> data = vtkSmartPointer<vtkStringArray>::New();
        data->SetNumberOfValues(1);
        data->SetName(fieldname.c_str());
        data->SetValue(0, str.c_str());

        vtkSmartPointer<vtkFieldData> field = polydata->GetFieldData();
        field->AddArray(data);
    }

    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getField(const std::string fieldname, const std::string datatype)
    {
        if (nEdges == 0 && datatype == "edge") {
            std::cout << "Edge data '" << fieldname << "' not supported for this instance of VTKData because there are no edges." << std::endl;
            exit(1);
        }

        vtkSmartPointer<vtkAbstractArray> arr = polydata->GetCellData()->GetAbstractArray(fieldname.c_str());
        if (!arr) {
            arr = polydata->GetPointData()->GetAbstractArray(fieldname.c_str());
            if (!arr) {
                std::cout << "fieldname " << fieldname << " not found!" << std::endl;
                exit(1);
            }
        }

        const int N = datatype == "vertex" ? nVertices : (nEdges > 0 ? nEdges + nFaces : (datatype == "edge" ? nEdges : nFaces));
        const int D = arr->GetNumberOfComponents();

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> field = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Constant(N, D, 0);

        if (arr->GetDataType() == 1) { // bit data
            vtkSmartPointer<vtkBitArray> arr2 = vtkBitArray::SafeDownCast(arr);

            bool value;
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < D; j++) {
                    field(i, j) = arr2->GetValue(D * (i + (datatype == "face" ? nEdges : 0)) + j);
                }
            }
        }
        else {
            vtkSmartPointer<typename VTKArrayType<T>::type> arr2 = VTKArrayType<T>::type::SafeDownCast(arr);
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < D; j++) {
                    field(i, j) = arr2->GetValue(D * (i + (datatype == "face" ? nEdges : 0)) + j);
                }
            }
        }

        return field;
    }

    template <typename Derived>
    void setField(const Eigen::MatrixBase<Derived> & field, const std::string fieldname, const std::string datatype)
    {
        if (nEdges == 0 && datatype == "edge") {
            std::cout << "Edge data '" << fieldname << "' not supported for this instance of VTKData because there are no edges." << std::endl;
            exit(1);
        }

        vtkSmartPointer<typename VTKArrayType<typename Derived::Scalar>::type> vector_array = convertFieldToVTK(field, fieldname, datatype);

        if (datatype == "vertex") {
            polydata->GetPointData()->AddArray(vector_array);
        }
        else {
            polydata->GetCellData()->AddArray(vector_array);
        }
    }

    Eigen::MatrixXd getVertices()
    {
        Eigen::MatrixXd vertices(nVertices, 3);
        for (int i = 0; i < nVertices; ++i) {
            double x[3];
            polydata->GetPoint(i, x);
            for (int d = 0; d < 3; ++d) {
                vertices(i, d) = x[d];
            }
        }
        return vertices;
    }

    Eigen::MatrixXi getEdges()
    {
        if (nEdges == 0) {
            std::cout << "Cannot get edges from VTKData because there are none." << std::endl;
            exit(1);
        }

        Eigen::MatrixXi edges(nEdges, 2);
        for (int i = 0; i < nEdges; ++i) {
            vtkSmartPointer<vtkIdList> ids_list = vtkSmartPointer<vtkIdList>::New();
            polydata->GetCellPoints(i, ids_list);
            for (int d = 0; d < 2; ++d) {
                edges(i, d) = ids_list->GetId(d);
            }
        }
        return edges;
    }

    Eigen::MatrixXi getFaces()
    {
        Eigen::MatrixXi faces(nFaces, 3);
        for (int i = 0; i < nFaces; ++i) {
            vtkSmartPointer<vtkIdList> ids_list = vtkSmartPointer<vtkIdList>::New();
            polydata->GetCellPoints(nEdges + i, ids_list);
            for (int d = 0; d < 3; ++d) {
                faces(i, d) = ids_list->GetId(d);
            }
        }
        return faces;
    }

    void setVertices(const Eigen::MatrixXd & vertices)
    {
        nVertices = vertices.rows();
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        points->SetNumberOfPoints(nVertices);
        for(int i=0;i<nVertices;++i) {
            points->SetPoint(i, vertices(i,0), vertices(i,1), vertices(i,2) );
        }
        polydata->SetPoints(points);
    }

    void setEdges(const Eigen::MatrixXi & edges)
    {
        nEdges = edges.rows();

        // create edges/lines
        vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
        for(int i=0;i<nEdges;++i) {
            // Create a line
            vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();

            line->GetPointIds()->SetId(0, edges(i,0));
            line->GetPointIds()->SetId(1, edges(i,1));

            lines->InsertNextCell(line);
        }

        polydata->SetLines(lines);
    }

    void setFaces(const Eigen::MatrixXi & faces)
    {
        nFaces = faces.rows();
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

    int write(const std::string filename_out)
    {
        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
#if VTK_MAJOR_VERSION <= 5
        writer->SetInput(polydata);
#else
        writer->SetInputData(polydata);
#endif

        // check whether or not extension is already there
        const std::size_t found_vtp_ext = filename_out.rfind(".vtp");
        const std::string filename = (found_vtp_ext != std::string::npos ? filename_out.substr(0, found_vtp_ext) : filename_out);
        writer->SetFileName( (filename+".vtp").c_str());
        const int retval = writer->Write(); // return value 1 for success, 0 for failure

        return retval;
    }
};


#endif /* VTKData_hpp */
