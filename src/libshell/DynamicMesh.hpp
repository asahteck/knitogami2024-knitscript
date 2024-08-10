//
//  DynamicMesh.hpp
//  Elasticity
//
//  Created by Wim van Rees on 9/6/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef DynamicMesh_hpp
#define DynamicMesh_hpp

#include "Mesh.hpp"
#include "MaterialProperties.hpp"

template<typename tCCD, typename tRCD>
class DynamicBaseMesh : public BaseMesh<tCCD, tRCD>
{
protected:
    Eigen::VectorXd massVector;

    virtual void initMassVectorFromMaterialVectors(const Eigen::Ref<const Eigen::VectorXd> per_face_thickness, const Eigen::Ref<const Eigen::VectorXd> per_face_density) = 0;
public:

    virtual void updateMesh(
        const Eigen::MatrixXd & newRestVertices,
        const Eigen::MatrixXd & newCurrentVertices,
        const Eigen::MatrixXd & vNewCurrentVertices,
        const Eigen::MatrixXb & newVertexBoundaryConditions,
        const Eigen::MatrixXi & newFaces,
        const Eigen::VectorXd & newEdgeDirectors,
        const Eigen::VectorXd & vNewEdgeDirectors)
    {
        // keep a copy of old edges in a vector of triplets
        // ...

        // update topology
        this->topology.clear();
        this->topology.init(newRestVertices, newFaces);

        // update boundary conditions
        // clampFixedEdges ...
        this->boundaryConditions.clear();
        this->boundaryConditions.init(newVertexBoundaryConditions, this->topology.getNumberOfEdges());

        // update rest state
        this->restState.clear();
        this->restState.init(newRestVertices, newEdgeDirectors.size()); // rest edge directors automatically set to zero
        this->restState.update(this->topology, this->boundaryConditions);

        // update current state (with velocities)
        this->currentState.clear();
        this->currentState.init(newCurrentVertices, newEdgeDirectors, vNewCurrentVertices, vNewEdgeDirectors);
        this->currentState.update(this->topology, this->boundaryConditions);
    }

    // little ugly but we only have two materials for now, so it should be fine
    void initMassVector(const MaterialProperties<Material_Isotropic> & matprop)
    {
        const int nFaces = this->getNumberOfFaces();
        Eigen::VectorXd per_face_thickness(nFaces), per_face_density(nFaces);

        for(int i=0;i<nFaces;++i)
        {
            per_face_thickness(i) = matprop.getFaceMaterial(i).getThickness();
            per_face_density(i) = matprop.getFaceMaterial(i).getDensity();
        }

        initMassVectorFromMaterialVectors(per_face_thickness, per_face_density);
    }

    void initMassVector(const MaterialProperties<Material_Orthotropic> & matprop)
    {
        const int nFaces = this->getNumberOfFaces();
        Eigen::VectorXd per_face_thickness(nFaces), per_face_density(nFaces);

        for(int i=0;i<nFaces;++i)
        {
            per_face_thickness(i) = matprop.getFaceMaterial(i).getThickness();
            per_face_density(i) = matprop.getFaceMaterial(i).getDensity();
        }

        initMassVectorFromMaterialVectors(per_face_thickness, per_face_density);
    }

    const Eigen::Ref<const Eigen::VectorXd> getMassVector() const
    {
        assert(massVector.rows() > 0);
        return massVector;
    }

    virtual int getNumberOfVariables() const = 0;

    virtual Real computeKineticEnergy() const
    {
        const int nVariables = getNumberOfVariables();
        const auto velocities = this->currentState.getVelocities();
        Real kinEng = 0.0;
        for(int i=0;i<nVariables;++i)
            kinEng += 0.5 * velocities(i) * velocities(i) * massVector(i);

        return kinEng;
    }
};

class DynamicMesh : public DynamicBaseMesh<DCSDynamicCurrentConfiguration, DCSRestConfiguration>
{
protected:
    void initMassVectorFromMaterialVectors(const Eigen::Ref<const Eigen::VectorXd> per_face_thickness, const Eigen::Ref<const Eigen::VectorXd> per_face_density) override
    {
        const int nVertices = this->getNumberOfVertices();
        const int nEdges = this->getNumberOfEdges();

        this->massVector.resize(3*nVertices + nEdges);

        // first compute the areas
        Eigen::Map<Eigen::VectorXd> mass_vertices(this->massVector.data(), nVertices);
        Eigen::Map<Eigen::VectorXd> mass_edges(this->massVector.data() + 3*nVertices, nEdges);
        this->restState.computeAreaVectors(this->topology, mass_vertices, mass_edges);

        // repeat the vertices masses
        for(int i=0;i<nVertices;++i)
        {
            massVector(i + nVertices) = massVector(i);
            massVector(i + 2*nVertices) = massVector(i);
        }
        // compute the per-vertex and per-edge density and thickness
        Eigen::VectorXd per_vertex_density(nVertices), per_vertex_thickness(nVertices);
        Eigen::VectorXd per_edge_density(nEdges), per_edge_thickness(nEdges);
        {
            const auto face2vertices = this->getTopology().getFace2Vertices();
            const auto face2edges = this->getTopology().getFace2Edges();
            const int nFaces = this->getNumberOfFaces();

            const Eigen::VectorXd doubleFaceAreas = this->getRestConfiguration().computeDoubleFaceAreas(this->topology);
            per_vertex_density.setZero();
            per_vertex_thickness.setZero();
            per_edge_density.setZero();
            per_edge_thickness.setZero();

            for(int i=0;i<nFaces;++i)
            {
                const Real areaPart = doubleFaceAreas(i)/6.0;

                const Real density = per_face_density(i);
                const Real thickness = per_face_thickness(i);

                const int idx_v0 = face2vertices(i,0);
                const int idx_v1 = face2vertices(i,1);
                const int idx_v2 = face2vertices(i,2);

                per_vertex_density(idx_v0) += density*areaPart;
                per_vertex_density(idx_v1) += density*areaPart;
                per_vertex_density(idx_v2) += density*areaPart;

                per_vertex_thickness(idx_v0) += thickness*areaPart;
                per_vertex_thickness(idx_v1) += thickness*areaPart;
                per_vertex_thickness(idx_v2) += thickness*areaPart;

                const int idx_e0 = face2edges(i,0);
                const int idx_e1 = face2edges(i,1);
                const int idx_e2 = face2edges(i,2);

                per_edge_density(idx_e0) += density*areaPart;
                per_edge_density(idx_e1) += density*areaPart;
                per_edge_density(idx_e2) += density*areaPart;

                per_edge_thickness(idx_e0) += thickness*areaPart;
                per_edge_thickness(idx_e1) += thickness*areaPart;
                per_edge_thickness(idx_e2) += thickness*areaPart;
            }

            // normalize
            for(int i=0;i<nVertices;++i)
            {
                per_vertex_density(i) /= mass_vertices(i);
                per_vertex_thickness(i) /= mass_vertices(i);
            }
            for(int i=0;i<nEdges;++i)
            {
                per_edge_density(i) /= mass_edges(i);
                per_edge_thickness(i) /= mass_edges(i);
            }
        }

        // now we use the material data
        for(int i=0;i<nVertices;++i)
        {
            const Real dens = per_vertex_density(i);
            const Real thickness = per_vertex_thickness(i);
            for(int d=0;d<3;++d)
                massVector(i + d*nVertices) *= dens * thickness;
        }

        const auto edge2vertices = this->getTopology().getEdge2Vertices();
        for(int i=0;i<nEdges;++i)
        {
            const Real dens = per_edge_density(i);
            const Real thickness = per_edge_thickness(i);
            massVector(i + 3*nVertices) *= dens * std::pow(thickness,3) / 12.;
        }
    }

public:
    int getNumberOfVariables() const override
    {
        const int nVertices = this->getNumberOfVertices();
        const int nEdges = this->getNumberOfEdges();
        return 3*nVertices + nEdges;
    }

    std::array<Real, 3> computeEdgeLengths() const
    {
        const int nEdges = this->getNumberOfEdges();
        const auto edge2vertices = this->topology.getEdge2Vertices();
        const auto vertices = this->getCurrentConfiguration().getVertices();

        std::array<Real, 3> retval = {std::numeric_limits<Real>::max(), 0, -1}; // min/avg/max
        for(int i=0;i<nEdges;++i)
        {
            const int idx_v0 = edge2vertices(i,0);
            const int idx_v1 = edge2vertices(i,1);
            const Real el = (vertices.row(idx_v0) - vertices.row(idx_v1)).norm();

            retval[0] = std::min(retval[0], el);
            retval[1] += el;
            retval[2] = std::max(retval[2], el);
        }
        retval[1] /= nEdges;

        return retval;
    }
};

#endif /* DynamicMesh_hpp */
