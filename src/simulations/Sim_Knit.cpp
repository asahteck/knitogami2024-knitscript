#include "Sim_Knit.hpp"
#include "Geometry_Knit.hpp"
#include "MaterialProperties.hpp"
#include "CombinedOperator_Parametric.hpp"
#include "GrowthHelper.hpp"
#include "ArgumentParser.hpp"
#include <igl/boundary_loop.h>

#include "Parametrizer.hpp"
#include "HLBFGS_Wrapper_Parametrized.hpp"


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
    mesh.init(geometry);

    // Define material parameters
    const Real E = parser.parse<Real>("-E", 1.0);
    const Real nu = parser.parse<Real>("-nu", 0.4);
    const Real h = parser.parse<Real>("-h", 0.5);

    const Real boundarySpringConstant = parser.parse<Real>("-boundarySpringConstant", 0.0);

    MaterialProperties_Iso_Constant matprop(E, nu, h);
    CombinedOperator_Parametric<tMesh, Material_Isotropic> engOp(matprop);

    // Mesh info
    auto rvertices = mesh.getRestConfiguration().getVertices();
    auto faces = mesh.getTopology().getFace2Vertices();
    const int nVertices = rvertices.rows();
    const int nFaces = faces.rows();
    tVecMat2d & bforms = mesh.getRestConfiguration().getSecondFundamentalForms();

    // Get geometric data
    Eigen::VectorXd stitchData = geometry.getCurvatureData(rvertices, faces);

    // Define constant natural curvatures
    const Real xCurv = parser.parse<Real>("-xCurv", -1.0);
    const Real yCurv = parser.parse<Real>("-yCurv", 1.0);
    Eigen::VectorXd curvatureAngles = Eigen::VectorXd::Constant(nFaces, 0.0);
    Eigen::VectorXd curvatures_p = xCurv * stitchData;
    Eigen::VectorXd curvatures_o = yCurv * stitchData;

    // Curve incrementally
    const int nSteps = parser.parse<int>("-nSteps", 10);
    for (int step = 0; step <= nSteps; step++) {
        const Real s = (1.0 * step) / (1.0 * nSteps);

        if (step > 0) {
            // Update bform on faces
            GrowthHelper<tMesh>::computeBbarsOrthoGrowthViaBbar(mesh, curvatureAngles, s * curvatures_p, s * curvatures_o, bforms);

            addNoiseToVertices_c<2>(0.01 * h);

            // Minimize energy
            if (boundarySpringConstant == 0.0) { // Free boundary condition
                Real eps = 1e-2;
                minimizeEnergy(engOp, eps);
            } else { // Boundary springs
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
        }

        dump(tag + "_f" + std::to_string(s), "stitchData", stitchData);
    }
}
