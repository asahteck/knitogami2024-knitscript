//
//  Sim_Knit.cpp
//  Elasticity
//
//  Created by Wim van Rees on 9/14/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#include "Sim_Knit.hpp"
#include "Geometry_Knit.hpp"
#include "MaterialProperties.hpp"
#include "CombinedOperator_Parametric.hpp"
#include "GrowthHelper.hpp"
#include "ArgumentParser.hpp"
#include <igl/boundary_loop.h>

#ifdef KNIT_BILAYER_MESH
#include "EnergyOperatorList.hpp"
#endif


#include "Parametrizer.hpp"
#include "HLBFGS_Wrapper_Parametrized.hpp"


template<typename tMesh>
class Parametrizer_Knit : public Parametrizer<tMesh>
{
public:
    using Parametrizer<tMesh>::mesh; // to avoid having to type this-> all the time
    using Parametrizer<tMesh>::data; // to avoid having to type this-> all the time

protected:
    const Real springConstant;
    const bool zOnly;

    // Containers for precomputed values that don't change after initialization
    Eigen::ArrayXd vertexAreas;

public:
    Parametrizer_Knit(tMesh & mesh_in, const Real springConstant, const bool zOnly):
    Parametrizer<tMesh>(mesh_in),
    springConstant(springConstant),
    zOnly(zOnly)
    {
        assert(springConstant >= 0.0);

        // precompute
        vertexAreas = mesh.getVertexAreas().array();
    }

    void initSolution(const Eigen::Ref<const Eigen::VectorXd>, const bool) override
    {
        // do nothing : we are use the mesh data
    }

    Real * getDataPointer() const override
    {
        return mesh.getDataPointer();
    }

    int getNumberOfVariables() const override
    {
        return 3 * mesh.getNumberOfVertices() + mesh.getNumberOfEdges(); // DCS energy minimization
    }

    void updateSolution() override
    {
        mesh.updateDeformedConfiguration();
    }

    Real computeEnergyContribution() override
    {
        // here we have a function that penalizes
        const int nVertices = mesh.getNumberOfVertices();
        const auto vertices = mesh.getCurrentConfiguration().getVertices();
        const auto rvertices = mesh.getRestConfiguration().getVertices();

        Eigen::MatrixXd displacements = vertices - rvertices;

        if (!zOnly) {
            Eigen::ArrayXd springLengths = displacements.rowwise().norm().array();
            return 0.5 * springConstant * (springLengths * springLengths * vertexAreas).sum();
        }
        else {
            Eigen::ArrayXd zDisplacements = displacements.col(2).array();
            return 0.5 * springConstant * (zDisplacements * zDisplacements * vertexAreas).sum();
        }
    }

    void updateGradient(const int nVars, const Eigen::Ref<const Eigen::VectorXd> energyGradient, Real * const grad_ptr) override
    {
        assert(nVars == getNumberOfVariables());

        // set grad_ptr equal to the actual gradient
        for (int i=0; i<nVars; ++i) {
            grad_ptr[i] = energyGradient(i);
        }

        const int nVertices = mesh.getNumberOfVertices();
        Eigen::Map<Eigen::VectorXd> gradVertices_x(grad_ptr, nVertices);
        Eigen::Map<Eigen::VectorXd> gradVertices_y(grad_ptr + nVertices, nVertices);
        Eigen::Map<Eigen::VectorXd> gradVertices_z(grad_ptr + 2*nVertices, nVertices);

        const auto vertices = mesh.getCurrentConfiguration().getVertices();
        const auto rvertices = mesh.getRestConfiguration().getVertices();

        // add the contribution from the mollifier
        Eigen::MatrixXd displacements = vertices - rvertices;

        if (!zOnly) {
            Eigen::VectorXd springLengths = displacements.rowwise().norm();
            Eigen::MatrixXd normalizedDisplacements = displacements.rowwise().normalized();
            Eigen::ArrayXd springMagnitudes = springConstant * springLengths.array() * vertexAreas;

            gradVertices_x += (springMagnitudes * normalizedDisplacements.col(0).array()).matrix();
            gradVertices_y += (springMagnitudes * normalizedDisplacements.col(1).array()).matrix();
            gradVertices_z += (springMagnitudes * normalizedDisplacements.col(2).array()).matrix();
        }
        else {
            gradVertices_z += (springConstant * displacements.col(2).array() * vertexAreas).matrix();
        }
    }
};


template<typename tMesh>
class Parametrizer_Knit_Boundary : public Parametrizer<tMesh>
{
public:
    using Parametrizer<tMesh>::mesh; // to avoid having to type this-> all the time
    using Parametrizer<tMesh>::data; // to avoid having to type this-> all the time

protected:
    const Real springConstant;

    // Containers for precomputed values that don't change after initialization
    Eigen::ArrayXd boundaryMask = Eigen::ArrayXd::Zero(mesh.getNumberOfVertices());

public:
    Parametrizer_Knit_Boundary(tMesh & mesh_in, const Real springConstant, const Real res):
    Parametrizer<tMesh>(mesh_in),
    springConstant(springConstant * res)
    {
        assert(springConstant >= 0.0);

        Eigen::ArrayXi boundaryVertices;
        igl::boundary_loop(mesh.getTopology().getFace2Vertices(), boundaryVertices);

        for (const auto i : boundaryVertices) {
            boundaryMask(i) = 1;
        }
    }

    void initSolution(const Eigen::Ref<const Eigen::VectorXd>, const bool) override
    {
        // do nothing : we are use the mesh data
    }

    Real * getDataPointer() const override
    {
        return mesh.getDataPointer();
    }

    int getNumberOfVariables() const override
    {
        return 3 * mesh.getNumberOfVertices() + mesh.getNumberOfEdges(); // DCS energy minimization
    }

    void updateSolution() override
    {
        mesh.updateDeformedConfiguration();
    }

    Real computeEnergyContribution() override
    {
        // here we have a function that penalizes
        const int nVertices = mesh.getNumberOfVertices();
        const auto vertices = mesh.getCurrentConfiguration().getVertices();
        const auto rvertices = mesh.getRestConfiguration().getVertices();

        Eigen::MatrixXd displacements = vertices - rvertices;
        Eigen::ArrayXd springLengths = displacements.rowwise().norm().array();
        return 0.5 * springConstant * (springLengths * springLengths * boundaryMask).sum();
    }

    void updateGradient(const int nVars, const Eigen::Ref<const Eigen::VectorXd> energyGradient, Real * const grad_ptr) override
    {
        assert(nVars == getNumberOfVariables());

        // set grad_ptr equal to the actual gradient
        for (int i=0; i<nVars; ++i) {
            grad_ptr[i] = energyGradient(i);
        }

        const int nVertices = mesh.getNumberOfVertices();
        Eigen::Map<Eigen::VectorXd> gradVertices_x(grad_ptr, nVertices);
        Eigen::Map<Eigen::VectorXd> gradVertices_y(grad_ptr + nVertices, nVertices);
        Eigen::Map<Eigen::VectorXd> gradVertices_z(grad_ptr + 2*nVertices, nVertices);

        const auto vertices = mesh.getCurrentConfiguration().getVertices();
        const auto rvertices = mesh.getRestConfiguration().getVertices();

        // add the contribution from the mollifier
        Eigen::MatrixXd displacements = vertices - rvertices;
        Eigen::VectorXd springLengths = displacements.rowwise().norm();
        Eigen::MatrixXd normalizedDisplacements = displacements.rowwise().normalized();
        Eigen::ArrayXd springMagnitudes = springConstant * springLengths.array() * boundaryMask;

        gradVertices_x += (springMagnitudes * normalizedDisplacements.col(0).array()).matrix();
        gradVertices_y += (springMagnitudes * normalizedDisplacements.col(1).array()).matrix();
        gradVertices_z += (springMagnitudes * normalizedDisplacements.col(2).array()).matrix();
    }
};


void Sim_Knit::init()
{}


void Sim_Knit::run()
{
    runKnit();
}


void Sim_Knit::runKnit()
{
    const std::string filename = parser.parse<std::string>("-filename", "knit.txt");
    tag = filename.substr(0, filename.length() - 4);

    // Initialize geometry/mesh
    const Real res = parser.parse<Real>("-res", 0.1);
    KnitPlate geometry(parser);

    if (!load()) {
        mesh.init(geometry);
    }

    // Define material parameters
    const Real E = parser.parse<Real>("-E", 1.0);
    const Real nu = parser.parse<Real>("-nu", 0.4);
    const Real h = parser.parse<Real>("-h", 0.5);

    const Real boundarySpringConstant = parser.parse<Real>("-boundarySpringConstant", 0.0);
    const Real springConstant = parser.parse<Real>("-springConstant", 0.0);
    const Real zOnly = parser.parse<bool>("-zOnly", false);

#ifdef KNIT_BILAYER_MESH
    MaterialProperties_Iso_Constant matprop_bottom(E, nu, h);
    MaterialProperties_Iso_Constant matprop_top(E, nu, h);
    CombinedOperator_Parametric<tMesh, Material_Isotropic, bottom> engOp_bot(matprop_bottom);
    CombinedOperator_Parametric<tMesh, Material_Isotropic, top> engOp_top(matprop_top);
    EnergyOperatorList<tMesh> engOp({&engOp_bot, &engOp_top});
#else
    MaterialProperties_Iso_Constant matprop(E, nu, h);
    CombinedOperator_Parametric<tMesh, Material_Isotropic> engOp(matprop);
#endif

    // Mesh info
    auto rvertices = mesh.getRestConfiguration().getVertices();
    auto faces = mesh.getTopology().getFace2Vertices();
    const int nVertices = rvertices.rows();
    const int nFaces = faces.rows();
#ifdef KNIT_BILAYER_MESH
    tVecMat2d & aforms_top = mesh.getRestConfiguration().getFirstFundamentalForms<top>();
    tVecMat2d & aforms_bottom = mesh.getRestConfiguration().getFirstFundamentalForms<bottom>();
#else
    tVecMat2d & bforms = mesh.getRestConfiguration().getSecondFundamentalForms();
#endif

    // Get geometric data
    Eigen::VectorXd curvatureData = geometry.getCurvatureData(rvertices, faces);

    // Define constant natural curvatures
    const Real xCurv = parser.parse<Real>("-xCurv", -2.0);
    const Real yCurv = parser.parse<Real>("-yCurv", 2.0);
    Eigen::VectorXd growthAngles = Eigen::VectorXd::Constant(nFaces, 0.0);
    Eigen::VectorXd growthCurvatures_p = xCurv * curvatureData;
    Eigen::VectorXd growthCurvatures_o = yCurv * curvatureData;

    // Curve incrementally
    const int nSteps = parser.parse<int>("-nSteps", 20);
    for (int step = 0; step <= nSteps; step++) {
        const Real s = (1.0 * step) / (1.0 * nSteps);

        if ((step > 0 && !cached()) || (step == nSteps)) {
            // Update "bform" on folded faces
#ifdef KNIT_BILAYER_MESH
            GrowthHelper<tMesh>::computeBbarsOrthoGrowth(mesh, h, growthAngles, s * growthCurvatures_p, s * growthCurvatures_o, aforms_top, aforms_bottom);
#else
            GrowthHelper<tMesh>::computeBbarsOrthoGrowthViaBbar(mesh, growthAngles, s * growthCurvatures_p, s * growthCurvatures_o, bforms);
#endif

            // Load lowres points (if requested)
            {
                const std::string lowresDir = parser.parse<std::string>("-lowresDir", "");
                if (lowresDir.size() > 1 && step < 2) {
                const std::string initFilename = lowresDir + "/" + tag + "_f" + std::to_string(s) + ".vtp";
                std::cout << "initFilename: " << initFilename << std::endl;

                std::ifstream f(initFilename.c_str());
                if (f.good()) {
                    VTKData data(initFilename);
                    const Eigen::MatrixXd initRestVertices = data.getField<Real>("rvertices", "vertex");
                    const Eigen::MatrixXd initVertices = data.getVertices();
                    const Eigen::MatrixXi initFaces = data.getFaces();

                    // load data
                    auto vertices = mesh.getCurrentConfiguration().getVertices();
                    vertices = mapPointsByMesh(rvertices, initRestVertices, initVertices, initFaces);
                }
            }
            }

            addNoiseToVertices_c<2>(0.01 * h);

            // Minimize energy (potentially using parametrizer)
            if (springConstant == 0.0 && boundarySpringConstant == 0.0) {
                Real eps = 1e-2;
                minimizeEnergy(engOp, eps);
            }
            else if (boundarySpringConstant != 0.0) {
                Real eps = 1e-3;

                Parametrizer_Knit_Boundary<tMesh> parametrizer(mesh, boundarySpringConstant, res);
                HLBFGS_Methods::HLBFGS_EnergyOp_Parametrized<tMesh, Parametrizer_Knit_Boundary, true> hlbfgs_wrapper(mesh, engOp, parametrizer);

                const Real epsMin = std::numeric_limits<Real>::epsilon();
                hlbfgs_wrapper.minimize(tag + "_diagnostics.dat", epsMin);
                eps = hlbfgs_wrapper.get_lastnorm();

                { // store energies
                    std::vector<std::pair<std::string, Real>> energies;
                    engOp.addEnergy(energies);
                    currentEnergies = energies;
                    FILE * f = fopen((tag+"_energies.dat").c_str(), "a");
                    for(const auto & eng : energies)
                    {
                        fprintf(f, "%s \t\t %10.10e\n", eng.first.c_str(), eng.second);
                    }
                    fclose(f);
                }
            }
            else if (springConstant != 0.0) {
                Real eps = 1e-3;

                Parametrizer_Knit<tMesh> parametrizer(mesh, springConstant, zOnly);
                HLBFGS_Methods::HLBFGS_EnergyOp_Parametrized<tMesh, Parametrizer_Knit, true> hlbfgs_wrapper(mesh, engOp, parametrizer);

                const Real epsMin = std::numeric_limits<Real>::epsilon();
                hlbfgs_wrapper.minimize(tag + "_diagnostics.dat", epsMin);
                eps = hlbfgs_wrapper.get_lastnorm();

                { // store energies
                    std::vector<std::pair<std::string, Real>> energies;
                    engOp.addEnergy(energies);
                    currentEnergies = energies;
                    FILE * f = fopen((tag+"_energies.dat").c_str(), "a");
                    for(const auto & eng : energies)
                    {
                        fprintf(f, "%s \t\t %10.10e\n", eng.first.c_str(), eng.second);
                    }
                    fclose(f);
                }
            }
        }

        dump(tag + "_f" + std::to_string(s), "curvatureData", curvatureData);
    }
}
