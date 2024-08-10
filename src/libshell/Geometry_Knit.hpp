//
//  Geometry_Knit.hpp
//  Elasticity
//
//  Created by Wim van Rees on 2/24/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef Geometry_Knit_h
#define Geometry_Knit_h

#include <functional>

#include "common.hpp"
#include "GeometryBase.hpp"


Eigen::MatrixXi getKnitData(const std::string filename)
{
    Eigen::MatrixXi data;
    std::vector<std::vector<int> > jagged_data;
    int cols = 0;
    int nx = 1;
    int ny = 1;

    std::ifstream file(filename);
    if (file.is_open()) {
        std::string line;

        while (std::getline(file, line)) {
            std::vector<int> row;
            std::istringstream ss(line);
            std::string temp;
            bool macro_line = false;

            ss >> temp;
            while (ss) {
                if (!temp.empty()) {
                    if (temp == "K") {
                        row.push_back(1);
                    }
                    else if (temp == "P") {
                        row.push_back(0);
                    }
                    else if (temp == "#") {
                        macro_line = true;
                        ss >> ny >> nx;
                    }
                }

                ss >> temp;
            }

            if (!macro_line) {
                jagged_data.push_back(row);

                if (cols < row.size()) {
                    cols = row.size();
                }
            }
        }
        file.close();
    }
    else {
        std::cout << "Error reading file '" << filename << "'." << std::endl;
        exit(1);
    }

    data.resize(jagged_data.size(), cols);
    for (int i = 0; i < jagged_data.size(); i++) {
        for (int j = 0; j < jagged_data[i].size(); j++) {
            data(i, j) = jagged_data[i][j];
        }
    }

    std::cout << nx << ", " << ny << std::endl;

    Eigen::MatrixXi data2 = data.replicate(ny, nx);

    std::cout << "getKnitData:" << std::endl;
    std::cout << data2 << std::endl;

    return data2;
}


class KnitPlate : public PlateGeometry
{
protected:
    const Real dx;
    const Real dy;
    const Real res;
    const std::string filename;

    Eigen::MatrixXi knitData;

    virtual void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        const int rows = knitData.rows();
        const int cols = knitData.cols();

        const Real halfX = 0.5 * dx * cols;
        const Real halfY = 0.5 * dy * rows;

        Eigen::MatrixXd polygon_vertices((rows + 1) * (cols + 1), 2);
        Eigen::MatrixXi polygon_edges(rows * (cols + 1) + (rows + 1) * cols, 2);

        {
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    polygon_vertices(i * (cols + 1) + j, 0) = -halfX + j * dx;
                    polygon_vertices(i * (cols + 1) + j, 1) = -halfY + i * dy;
                    polygon_edges(i * cols + j, 0) = i * (cols + 1) + j;
                    polygon_edges(i * cols + j, 1) = i * (cols + 1) + j + 1;
                    polygon_edges((rows + 1) * cols + i * (cols + 1) + j, 0) = i * (cols + 1) + j;
                    polygon_edges((rows + 1) * cols + i * (cols + 1) + j, 1) = (i + 1) * (cols + 1) + j;
                }

                polygon_vertices((i + 1) * (cols + 1) - 1, 0) = halfX;
                polygon_vertices((i + 1) * (cols + 1) - 1, 1) = -halfY + i * dy;
                polygon_edges((rows + 1) * cols + i * (cols + 1) + cols, 0) = (i + 1) * (cols + 1) - 1;
                polygon_edges((rows + 1) * cols + i * (cols + 1) + cols, 1) = (i + 2) * (cols + 1) - 1;
            }

            for (int j = 0; j < cols; j++) {
                polygon_vertices(rows * (cols + 1) + j, 0) = -halfX + j * dx;
                polygon_vertices(rows * (cols + 1) + j, 1) = halfY;
                polygon_edges(rows * cols + j, 0) = rows * (cols + 1) + j;
                polygon_edges(rows * cols + j, 1) = rows * (cols + 1) + j + 1;
            }

            polygon_vertices((rows + 1) * (cols + 1) - 1, 0) = halfX;
            polygon_vertices((rows + 1) * (cols + 1) - 1, 1) = halfY;
        }

        // Triangulation flags
        const std::string flags = "q20a"+std::to_string(std::sqrt(3.0)/4.0*res*res)+(verbose ? "" : "Q");

        Triangulate triangulate(polygon_vertices, polygon_edges, flags);
        triangulate.get(vertices, face2vertices, vertices_bc);
    }

public:
    KnitPlate(ArgumentParser & parser):
    PlateGeometry(),
    dx(parser.parse<Real>("-dx", 1.4)),
    dy(parser.parse<Real>("-dy", 1.2)),
    res(parser.parse<Real>("-res", 0.2)),
    filename(parser.parse<std::string>("-filename", "knit.txt"))
    {
        knitData = getKnitData(filename);
    }

    Eigen::VectorXd getCurvatureData(const Eigen::MatrixXd & rvertices, const Eigen::MatrixXi & faces)
    {
        Eigen::VectorXd curvatures(faces.rows());
        curvatures.setZero();

        const int rows = knitData.rows();
        const int cols = knitData.cols();

        const Real halfX = 0.5 * dx * cols;
        const Real halfY = 0.5 * dy * rows;

        Real xc, yc;
        int row, col;
        for (int i = 0; i < faces.rows(); i++) {
            xc = (rvertices(faces(i, 0), 0) + rvertices(faces(i, 1), 0) + rvertices(faces(i, 2), 0)) / 3.0;
            yc = (rvertices(faces(i, 0), 1) + rvertices(faces(i, 1), 1) + rvertices(faces(i, 2), 1)) / 3.0;
            row = std::floor((yc + halfY) / dy);
            col = std::floor((xc + halfX) / dx);
            curvatures(i) = (knitData(row, col) == 1) ? 1.0 : -1.0;
        }

        return curvatures;
    }
};


class MiuraOriPlate : public PlateGeometry
{
protected:
    const int nx;
    const int ny;
    const Real dx;
    const Real dy;
    const Real res;

    virtual void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        Eigen::MatrixXd polygon_vertices((nx + 1) * (ny + 1), 2);
        Eigen::MatrixXi polygon_edges(nx * ny + nx * (ny + 1) + (nx + 1) * ny, 2);

        {
            for (int i = 0; i <= ny; i++) {
                for (int j = 0; j <= nx; j++) {
                    const int idx = i * (nx + 1) + j;

                    polygon_vertices(idx, 0) = j * dx;
                    polygon_vertices(idx, 1) = i * dy;

                    if (j < nx && i < ny) {
                        polygon_edges(i * nx + j, 0) = idx + (i % 2 == 0 ? 0 : 1);
                        polygon_edges(i * nx + j, 1) = idx + (nx + 1) + (i % 2 == 0 ? 1 : 0);
                    }

                    if (j < nx) {
                        polygon_edges(nx * ny + i * nx + j, 0) = idx;
                        polygon_edges(nx * ny + i * nx + j, 1) = idx + 1;
                    }

                    if (i < ny) {
                        polygon_edges(nx * ny + nx * (ny + 1) + i * (nx + 1) + j, 0) = idx;
                        polygon_edges(nx * ny + nx * (ny + 1) + i * (nx + 1) + j, 1) = idx + (nx + 1);
                    }
                }
            }
        }

        // Triangulation flags
        const std::string flags = "q20a"+std::to_string(std::sqrt(3.0)/4.0*res*res)+(verbose ? "" : "Q");

        Triangulate triangulate(polygon_vertices, polygon_edges, flags);
        triangulate.get(vertices, face2vertices, vertices_bc);
    }

public:
    MiuraOriPlate(ArgumentParser & parser):
    PlateGeometry(),
    nx(parser.parse<int>("-nx", 8)),
    ny(parser.parse<int>("-ny", 8)),
    dx(parser.parse<Real>("-dx", 1.4)),
    dy(parser.parse<Real>("-dy", 1.2)),
    res(parser.parse<Real>("-res", 0.2))
    {}

    Eigen::VectorXd getCurvatureData(const Eigen::MatrixXd & rvertices, const Eigen::MatrixXi & faces)
    {
        Eigen::VectorXd curvatures(faces.rows());
        curvatures.setZero();

        Real xc, yc;
        int row, col;
        for (int i = 0; i < faces.rows(); i++) {
            xc = (rvertices(faces(i, 0), 0) + rvertices(faces(i, 1), 0) + rvertices(faces(i, 2), 0)) / 3.0;
            yc = (rvertices(faces(i, 0), 1) + rvertices(faces(i, 1), 1) + rvertices(faces(i, 2), 1)) / 3.0;
            row = std::floor(yc / dy);
            col = std::floor(xc / dx);

            if (row % 2 == 0 && col % 2 == 0) {
                curvatures(i) = ((yc - row * dy) <= (dy / dx) * (xc - col * dx)) ? 1.0 : -1.0;
            }
            else if (row % 2 == 0 && col % 2 == 1) {
                curvatures(i) = ((yc - row * dy) > (dy / dx) * (xc - col * dx)) ? 1.0 : -1.0;
            }
            else if (row % 2 == 1 && col % 2 == 0) {
                curvatures(i) = ((yc - row * dy) > -(dy / dx) * (xc - (col + 1) * dx)) ? 1.0 : -1.0;
            }
            else if (row % 2 == 1 && col % 2 == 1) {
                curvatures(i) = ((yc - row * dy) <= -(dy / dx) * (xc - (col + 1) * dx)) ? 1.0 : -1.0;
            }
        }

        return curvatures;
    }
};


#endif /* Geometry_Knit_h */
