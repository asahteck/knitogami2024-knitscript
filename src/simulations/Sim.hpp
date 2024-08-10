//
//  Sim.hpp
//  Elasticity
//
//  Created by Wim van Rees on 21/02/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef Sim_hpp
#define Sim_hpp

#include <iostream>

#include "common.hpp"
#include "ArgumentParser.hpp"

#ifdef USEVTK
#include "ReadVTK.hpp"
#include "WriteVTK.hpp"
#include "VTKData.hpp"
#endif

#include "WriteSTL.hpp"
#include "ComputeCurvatures.hpp"
#include "Umbilics.hpp"

#include "Mesh.hpp"
#include "MaterialProperties.hpp"
#include "EnergyOperator.hpp"
#include "CombinedOperator_Parametric.hpp"

#ifdef USELIBLBFGS
#include "LBFGS_Wrapper.hpp"
#endif

#ifdef USEHLBFGS
#include "HLBFGS_Wrapper.hpp"
#endif

#include <igl/writeOBJ.h>


/*! \class Sim
 * \brief Base class for simulations.
 *
 * This class is called from main, and performs the main simulation. Every class that derives from here can implement a simulation case.
 */
class BaseSim
{
protected:
    ArgumentParser & parser;

public:
    BaseSim(ArgumentParser & parser_in):
    parser(parser_in)
    {
        parser.save_defaults(); // make sure all default values are printed as well from now on
        parser.save_options();
    }

    virtual void init() = 0;
    virtual void run() = 0;
    virtual int optimize()
    {
        std::cout << "Optimize not (yet) implemented for this class " << std::endl;
        return -1;
    };

    virtual ~BaseSim()
    {}
};


template<typename tMesh>
class Sim : public BaseSim
{
public:
    typedef typename tMesh::tCurrentConfigData tCurrentConfigData;
    typedef typename tMesh::tReferenceConfigData tReferenceConfigData;

private: // helper functions
    std::string withVTPExtension(const std::string filename)
    {
        const std::size_t found_vtp_ext = filename.rfind(".vtp");
        return (found_vtp_ext != std::string::npos ? filename.substr(0, found_vtp_ext) : filename) + ".vtp";
    }


    std::string withoutVTPExtension(const std::string filename)
    {
        const std::size_t found_vtp_ext = filename.rfind(".vtp");
        return (found_vtp_ext != std::string::npos ? filename.substr(0, found_vtp_ext) : filename);
    }

    // Returns false if no reference file exists
    bool loadFromReferenceFile(const std::string referenceFilename)
    {
        // Check if file exists
        const std::string fname = withVTPExtension(referenceFilename);
        if (FILE *file = fopen(fname.c_str(), "r")) {
            fclose(file);
        }
        else {
            std::cout << " -- Loading from reference file: no reference file found." << std::endl;
            return false;
        }

        try {
            auto cvertices = mesh.getCurrentConfiguration().getVertices();
            const int nVertices = mesh.getNumberOfVertices();

            vtkData.load(fname);

            Eigen::MatrixXd vertices = vtkData.getVertices();
            assert(vertices.rows() == nVertices && vertices.cols() == 3);

            Eigen::MatrixXi faces = vtkData.getFaces();
            assert(faces.array() == mesh.getTopology().getFace2Vertices().array());

            for (int i = 0; i < nVertices; i++) {
                for (int j = 0; j < 3; j++) {
                    cvertices(i, j) = vertices(i, j);
                }
            }
        }
        catch (...) {
            std::cout << " -- Loading from reference file: something went wrong! Aborting load." << std::endl;
            return false;
        }

        return true;
    }

protected:
    std::string tag;
    tMesh mesh;
    std::vector< std::pair<std::string, Real> > currentEnergies;
    const bool verbose = parser.parse<bool>("-verbose", true);

    VTKData vtkData;
    std::string filename;

    // `simStep` is the step in the current simulation run. `loadedStep` can be
    // set by loading a cached mesh with the proper attribute (aka "simStep").
    // Only `dump_` is responsible for keeping track of `simStep` and its
    // relation to `loadedStep`, since the cache file always mirrors a mesh
    // that `dump_` puts out.
    int simStep = 0;
    int loadedStep = -1;

    // Is a cache successfully loaded? (may not need)
    bool loaded() {
        return loadedStep > -1;
    }

    // Is the current sim step already computed?
    bool cached() {
        return simStep <= loadedStep;
    }

    virtual void writeSTL(const std::string filename)
    {
        WriteSTL::write(mesh.getTopology(), mesh.getCurrentConfiguration(), filename);
    }

    // save everything necessary to recreate the mesh. TODO: want ability to add extra fields
    void save(const std::string filename, const std::string filename_mirror)
    {
        VTKData data;

        const Eigen::MatrixXd vertices = mesh.getCurrentConfiguration().getVertices();
        const Eigen::MatrixXi edges = mesh.getTopology().getEdge2Vertices();
        const Eigen::MatrixXi faces = mesh.getTopology().getFace2Vertices();

        data.setVertices(vertices);
        data.setEdges(edges);
        data.setFaces(faces);

        // add arguments
        data.setStringAttribute(parser.options_string(), "arguments");

        // add mesh data
        const Eigen::MatrixXd rvertices = mesh.getRestConfiguration().getVertices();
        const Eigen::VectorXd directors = mesh.getCurrentConfiguration().getEdgeDirectors();
        const Eigen::VectorXd rdirectors = mesh.getRestConfiguration().getEdgeDirectors();
        data.setField(rvertices, "rvertices", "vertex");
        data.setField(directors, "directors", "edge");
        data.setField(rdirectors, "rdirectors", "edge");

        data.setField(toInt(mesh.getBoundaryConditions().getVertexBoundaryConditions()), "vertices_bc", "vertex");

        const bool clamped = mesh.getBoundaryConditions().getClampedEdgeNormals().size() != 0;
        data.setScalarAttribute(clamped, "clamped");
        if (clamped) {
            const Eigen::MatrixXd clampedEdgeNormals = mesh.getBoundaryConditions().getClampedEdgeNormals();
            const std::map<int, int> clampedEdgeNormalsIndices = mesh.getBoundaryConditions().getClampedEdgeNormalsIndices();
            Eigen::MatrixXd reduced_edge_normals = Eigen::MatrixXd::Constant(mesh.getNumberOfEdges(), 3, 0);

            for (auto const & map : clampedEdgeNormalsIndices) {
                reduced_edge_normals.row(map.first) = clampedEdgeNormals.row(map.second);
            }

            data.setField(reduced_edge_normals, "clampedEdgeNormals", "edge");
        }

        { // log
            const std::string filepath = "argumentparser.log";
            FILE * f = fopen(filepath.c_str(), "a");
            if (f == nullptr) {
                printf("Can not open file %s.\n", filepath.data());
            }
            else {
                fprintf(f, "\n<cache:%d:%s:%s>", simStep, withVTPExtension(filename).c_str(), withVTPExtension(filename_mirror).c_str());
                fclose(f);
            }
        }

        // update cache file to reflect the mesh it's intended to mirror
        data.setScalarAttribute(simStep, "simStep");
        data.setStringAttribute(withVTPExtension(filename_mirror), "filename_mirror");

        // write
        data.write(filename);
    }

    // Only works on files with edges and proper fields/attributes set!
    bool load(const std::string fname = "")
    {
        std::string filename = fname;

        // By default, load the cache file
        if (fname.size() == 0) {
            filename = tag + "_cache";
        }

        // Check if file exists
        if (FILE *file = fopen(withVTPExtension(filename).c_str(), "r")) {
            fclose(file);
        }
        else {
            std::cout << " -- Loading mesh: file " << withVTPExtension(filename) << " does not exist! Aborting load." << std::endl;
            return false;
        }

        try {
            vtkData.load(filename);

            Eigen::MatrixXd vertices = vtkData.getVertices();
            Eigen::MatrixXi edges = vtkData.getEdges();
            Eigen::MatrixXi faces = vtkData.getFaces();

            Eigen::MatrixXd rvertices;
            if (!vtkData.isField("rvertices")) {
                std::cout << "\t -- Loading mesh: rvertices not found! Defaulting to current vertices." << std::endl;
                rvertices = vertices;
            }
            else {
                rvertices = vtkData.getField<Real>("rvertices", "vertex");
            }

            Eigen::VectorXd rdirectors;
            if (!vtkData.isField("rdirectors")) {
                std::cout << "\t -- Loading mesh: rdirectors not found! Defaulting to all zeros." << std::endl;
                rdirectors = Eigen::VectorXd::Constant(edges.rows(), 0.0);
            }
            else {
                rdirectors = vtkData.getField<Real>("rdirectors", "edge");
            }

            Eigen::VectorXd directors;
            if (!vtkData.isField("directors")) {
                std::cout << "\t -- Loading mesh: directors not found! Defaulting to rest directors." << std::endl;
                directors = rdirectors;
            }
            else {
                directors = vtkData.getField<Real>("directors", "edge");
            }

            Eigen::MatrixXb vertices_bc;
            if (!vtkData.isField("vertices_bc")) {
                std::cout << "\t -- Loading mesh: vertices_bc not found! Defaulting to all false." << std::endl;
                vertices_bc = Eigen::MatrixXb::Constant(vertices.rows(), 3, false);
            }
            else {
                vertices_bc = toBool(vtkData.getField<int>("vertices_bc", "vertex"));
            }

            bool clamped;
            if (!vtkData.isAttribute("clamped")) {
                std::cout << "\t -- Loading mesh: 'clamped' attribute not found! Defaulting to false." << std::endl;
                clamped = false;
            }
            else {
                clamped = vtkData.getAttribute<bool>("clamped", false); // technically defval is redundant here
            }

            if (clamped && vtkData.isField("clampedEdgeNormals")) {
                Eigen::MatrixXd clampedEdgeNormals = vtkData.getField<Real>("clampedEdgeNormals", "edge");
                mesh.init(vertices, edges, faces, rvertices, directors, rdirectors, vertices_bc, clampedEdgeNormals);
            }
            else {
                if (clamped) { std::cout << "\t -- Loading mesh: 'clampedEdgeNormals' not found despite clamped attribute. Defaulting to none." << std::endl; }
                mesh.init(vertices, edges, faces, rvertices, directors, rdirectors, vertices_bc);
            }

            if (vtkData.isAttribute("simStep")) {
                loadedStep = vtkData.getAttribute<int>("simStep", -1);
                std::cout << " -- Loading mesh: iteration " << std::to_string(loadedStep) << std::endl;
            }

            return true;
        }
        catch (...) {
            std::cout << " -- Loading mesh: something went wrong! Aborting load." << std::endl;
            return false;
        }
    }

    void dump(const std::string filename, const bool dumpVerbose = false)
    {
        vtkData.clear();

        vtkData.setVertices(mesh.getCurrentConfiguration().getVertices());
        vtkData.setFaces(mesh.getTopology().getFace2Vertices());

        dump_(filename, dumpVerbose);
    }

    void dump(const std::string filename, const std::string fieldName, const Eigen::VectorXd & field, const bool dumpVerbose = false)
    {
        vtkData.clear();

        vtkData.setVertices(mesh.getCurrentConfiguration().getVertices());
        vtkData.setFaces(mesh.getTopology().getFace2Vertices());

        if (field.size() == mesh.getNumberOfFaces()) {
            vtkData.setField(field, fieldName, "face");
        }
        else if (field.size() == mesh.getNumberOfVertices()) {
            vtkData.setField(field, fieldName, "vertex");
        }
        else {
            std::cout << "Cannot dump with field `" << fieldName << "`; continuing without." << std::endl;
        }

        dump_(filename, dumpVerbose);
    }

    // add energy to dump for monolayer only
    template <typename T = tMesh, typename tMeshOperator>
    typename std::enable_if<std::is_same<T, Mesh>::value, void>::type
    dump(const std::string filename, const tMeshOperator & engOp, const bool dumpVerbose = false)
    {
        vtkData.clear();

        vtkData.setVertices(mesh.getCurrentConfiguration().getVertices());
        vtkData.setFaces(mesh.getTopology().getFace2Vertices());

        const int nFaces = mesh.getNumberOfFaces();

        // Compute energies for each face
        Eigen::MatrixXd perFaceEnergies(nFaces, 3);
        engOp.computePerFaceEnergies(mesh, perFaceEnergies);
        Eigen::VectorXd stretchingEnergy = perFaceEnergies.col(0);
        Eigen::VectorXd bendingEnergy = perFaceEnergies.col(1);

        vtkData.setField(stretchingEnergy, "stretchingEnergy", "face");
        vtkData.setField(bendingEnergy, "bendingEnergy", "face");

        dump_(filename, dumpVerbose);
    }

    // add energy to dump for bilayer only
    template <typename T = tMesh, typename tEngOp_bottom, typename tEngOp_top>
    typename std::enable_if<std::is_same<T, BilayerMesh>::value, void>::type
    dump(const std::string filename, const tEngOp_bottom & engOp_bottom, const tEngOp_top & engOp_top, const Real E, const Real nu, const Real h, const bool dumpVerbose = false)
    {
        vtkData.clear();

        vtkData.setVertices(mesh.getCurrentConfiguration().getVertices());
        vtkData.setFaces(mesh.getTopology().getFace2Vertices());

        const int nFaces = mesh.getNumberOfFaces();

        // Compute energies for each face
        Eigen::MatrixXd energy_bi_top(nFaces, 3);
        Eigen::MatrixXd energy_bi_bottom(nFaces, 3);
        energy_bi_top.setZero();
        energy_bi_bottom.setZero();

        const Real eng_bi_top = engOp_top.computePerFaceEnergies(mesh, energy_bi_top);
        const Real eng_bi_bottom = engOp_bottom.computePerFaceEnergies(mesh, energy_bi_bottom);

        vtkData.setField(energy_bi_bottom.col(0), "bi_bot_aa", "face");
        vtkData.setField(energy_bi_bottom.col(1), "bi_bot_bb", "face");
        vtkData.setField(energy_bi_bottom.col(2), "bi_bot_ab", "face");
        vtkData.setField(energy_bi_top.col(0), "bi_top_aa", "face");
        vtkData.setField(energy_bi_top.col(1), "bi_top_bb", "face");
        vtkData.setField(energy_bi_top.col(2), "bi_top_ab", "face");

        // copy bilayer mesh to a single layer mesh
        Mesh mesh_dummy;
        Geometry_Dummy geometry_dummy(mesh.getCurrentConfiguration().getVertices(), mesh.getTopology().getFace2Vertices());
        mesh_dummy.init(geometry_dummy);
        mesh_dummy.getCurrentConfiguration().getEdgeDirectors() = mesh.getCurrentConfiguration().getEdgeDirectors();
        mesh_dummy.updateDeformedConfiguration();

        // create new material properties and energy operator
        MaterialProperties_Iso_Constant matprop_mono(E, nu, h);
        Eigen::MatrixXd energy_mono(nFaces, 3);
        energy_mono.setZero();
        CombinedOperator_Parametric<Mesh, Material_Isotropic> engOp_mono(matprop_mono);

        // set metrics
        tVecMat2d & aform_tmp = mesh_dummy.getRestConfiguration().getFirstFundamentalForms();
        tVecMat2d & bform_tmp = mesh_dummy.getRestConfiguration().getSecondFundamentalForms();
        const tVecMat2d & aform_bot = mesh.getRestConfiguration().template getFirstFundamentalForms<bottom>();
        const tVecMat2d & aform_top = mesh.getRestConfiguration().template getFirstFundamentalForms<top>();
        for (int i = 0; i < nFaces; ++i) {
            const Real h1 = h;
            const Real h2 = h;
            const Real abot_prefac = h1*(4.0*h2*h2 - h1*h2 + h1*h1)/std::pow(h1 + h2,3);
            const Real atop_prefac = h2*(4.0*h1*h1 - h1*h2 + h2*h2)/std::pow(h1 + h2,3);
            const Real b_prefac = 6.0*h1*h2/std::pow(h1 + h2,3);

            aform_tmp[i] = abot_prefac*aform_bot[i] + atop_prefac*aform_top[i];
            bform_tmp[i] = b_prefac * (aform_bot[i] - aform_top[i]);
        }

        const Real eng_mono = engOp_mono.computePerFaceEnergies(mesh_dummy, energy_mono);
        vtkData.setField(energy_mono.col(0), "mono_aa", "face");
        vtkData.setField(energy_mono.col(1), "mono_bb", "face");

        dump_(filename, dumpVerbose);
    }

    virtual void dump_(const std::string filename, const bool dumpVerbose)
    {
        if (!cached()) { // this is a new mesh that we need to save
            save(tag + "_cache", filename);
            std::cout << " -- dump: iteration " << simStep << " saved in " << tag << "_cache to reflect " << withVTPExtension(filename) << std::endl;
            simStep++;
        }
        else { // this is an old mesh; do nothing and increment our internal sim step
            std::cout << " -- dump: iteration " << simStep << " has already been run; incrementing simStep" << std::endl;
            simStep++;
            return;
        }

        // add arguments
        vtkData.setStringAttribute(parser.options_string(), "arguments");

        // add total energies
        for (const auto & eng : currentEnergies) {
            vtkData.setScalarAttribute(eng.second, eng.first);
        }

        // add curvatures
        const int nFaces = mesh.getNumberOfFaces();
        Eigen::VectorXd gauss(nFaces);
        Eigen::VectorXd mean(nFaces);
        ComputeCurvatures<tMesh> computeCurvatures;
        if (dumpVerbose) {
            Eigen::MatrixXd pCurvature1(nFaces, 3);
            Eigen::MatrixXd pCurvature2(nFaces, 3);
            computeCurvatures.compute(mesh, gauss, mean, pCurvature1, pCurvature2, true);
            vtkData.setField(pCurvature1, "pCurvature1", "face");
            vtkData.setField(pCurvature2, "pCurvature2", "face");
        }
        else {
            computeCurvatures.compute(mesh, gauss, mean);
        }
        vtkData.setField(mean, "mean", "face");
        vtkData.setField(gauss, "gauss", "face");

        // rest vertices (always include)
        vtkData.setField(mesh.getRestConfiguration().getVertices(), "rvertices", "vertex");

        // write
        vtkData.write(filename);
    }

    void computeMomenta(const Eigen::VectorXd & massVertices, const Eigen::MatrixXd & vVertices, std::array<Real,3> & lin, std::array<Real,3> & ang)
    {
        lin[0] = lin[1] = lin[2] = 0.0;
        ang[0] = ang[1] = ang[2] = 0.0;

        const int nVertices = mesh.getNumberOfVertices();
        const auto vertices = mesh.getCurrentConfiguration().getVertices();

        const Eigen::Vector3d org_vec = (Eigen::Vector3d() << 0,0,0).finished();

        for(int i=0;i<nVertices;++i)
        {
            const Real mass = massVertices(i);
            const Eigen::Vector3d vertex_pos = vertices.row(i);
            const Eigen::Vector3d position = vertex_pos - org_vec;
            const Eigen::Vector3d velocity = vVertices.row(i);
            const Eigen::Vector3d linMom = mass*velocity;
            const Eigen::Vector3d angMom = position.cross(linMom);
            for(int d=0;d<3;++d)
            {
                lin[d] += linMom(d);
                ang[d] += angMom(d);
            }
        }
    }

    template<int component>
    void addNoiseToVertices_c(const Real ampl)
    {
        std::mt19937 gen;
        gen.seed(42);
        std::uniform_real_distribution<Real> distV(-ampl, ampl);
        auto perturb_v = [&](Eigen::Vector3d in)
        {
            Eigen::Vector3d retval;
            for(int d=0;d<3;++d)
                retval(d) = in(d) + (d == component ? distV(gen) : 0.0);
            return retval;
        };
        mesh.changeVertices(perturb_v);
    }

    void addNoiseToVertices(const Real ampl)
    {
        std::mt19937 gen;
        gen.seed(42);
        std::uniform_real_distribution<Real> distV(-ampl, ampl);
        auto perturb_v = [&](Eigen::Vector3d in)
        {
            Eigen::Vector3d retval;
            retval << in(0)+distV(gen), in(1)+distV(gen), in(2)+distV(gen);
            return retval;
        };
        mesh.changeVertices(perturb_v);
    }

    void addNoiseToEdgeDirectors(const Real ampl)
    {
        std::mt19937 gen;
        gen.seed(42);
        std::uniform_real_distribution<Real> distE(-ampl, ampl);
        auto perturb_e = [&](Real in)
        {
            return in + distE(gen);
        };

        mesh.changeEdgeDirectors(perturb_e);
    }

    template<typename tMeshOperator, bool verbose = true>
    int minimizeEnergy(const tMeshOperator & op, Real & eps, const std::string referenceFilename = "", const Real epsMin=std::numeric_limits<Real>::epsilon(), const bool stepWise = false)
    {
        if (cached()) {
            std::cout << " -- minimizeEnergy: iteration " << simStep << " is already finished; continuing without minimization" << std::endl;
            return 0;
        }
        else if (referenceFilename.size() > 0) {
            std::cout << " -- minimizeEnergy: loading from reference filename " << referenceFilename << std::endl;
            loadFromReferenceFile(referenceFilename);
        }

#ifdef USELIBLBFGS
        // use the LBFGS_Energy class to directly minimize the energy on the mesh with these operators
        // LBFGS does not use the hessian
        LBFGS::LBFGS_Energy<Real, tMesh, tMeshOperator, true> lbfgs_energy(mesh, op);
        int retval = 0;
        while(retval == 0 && eps > epsMin)
        {
            eps *= 0.1;
            retval = lbfgs_energy.minimize(eps);
        }
#else
#ifdef USEHLBFGS
        HLBFGS_Methods::HLBFGS_Energy<tMesh, tMeshOperator, true> hlbfgs_wrapper(mesh, op);
        int retval = 0;
        if(stepWise)
        {
            while(retval == 0 && eps > epsMin)
            {
                eps *= 0.1;
                retval = hlbfgs_wrapper.minimize(tag+"_diagnostics.dat", eps);
            }
        }
        else
        {
            // retval = hlbfgs_wrapper.minimize(tag+"_diagnostics.dat", epsMin);
            retval = hlbfgs_wrapper.minimize(tag+"_diagnostics.dat", epsMin); // 0 - pure gradient descent
            eps = hlbfgs_wrapper.get_lastnorm();
        }
#else
        std::cout << "should use liblbfgs or hlbfgs\n";
#endif
#endif

        // store energies
        {
            std::vector<std::pair<std::string, Real>> energies;
            op.addEnergy(energies);
            currentEnergies = energies;
            FILE * f = fopen((tag+"_energies.dat").c_str(), "a");
            for(const auto & eng : energies)
            {
                fprintf(f, "%s \t\t %10.10e\n", eng.first.c_str(), eng.second);
            }
            fclose(f);
        }
        return retval;
    }

public:
    Sim(ArgumentParser & parser_in):
    BaseSim(parser_in),
    tag("Sim")
    {}

    virtual ~Sim()
    {}
};

#endif /* Sim_hpp */
