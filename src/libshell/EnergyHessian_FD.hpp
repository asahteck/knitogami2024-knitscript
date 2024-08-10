//
//  EnergyHessian_FD.hpp
//  Elasticity
//
//  Created by Wim van Rees on 12/31/17.
//  Copyright © 2017 Wim van Rees. All rights reserved.
//

#ifndef EnergyHessian_FD_hpp
#define EnergyHessian_FD_hpp

#include "EnergyOperator.hpp"
#include "Profiler.hpp"
#include "MaterialProperties.hpp"
#include <set>


template<typename tMesh>
struct ComputeHessian
{
    tMesh & mesh;
    const EnergyOperator<tMesh> & energyoperator;
    const int order;
    const Real eps;
    const std::vector< std::vector <Real>> FD1_prefac;
    const std::vector< std::vector <Real>> FD1_offset;
    const std::vector<Real> FD1_denum;
    const int nPoints;
    
    ComputeHessian(tMesh & mesh, const EnergyOperator<tMesh> &  engOp, const int order, const Real eps):
    mesh(mesh),
    energyoperator(engOp),
    order(order),
    eps(eps),
    FD1_prefac({ {1, -1}, {1, -8, 8, -1}, {-1, 9, -45, 45, -9, 1}, {3, -32, 168, -672, 672, -168, 32, -3} }),
    FD1_offset({ {-1, 1}, {-2, -1, 1, 2}, {-3, -2, -1, 1, 2, 3}, {-4, -3, -2, -1, 1, 2, 3, 4} }),
    FD1_denum({2, 12, 60, 840}), // denumerator for first derivative
    nPoints(FD1_prefac[order].size())
    {}
    
    inline Real compute(Real & val1, Real & val2)
    {
        Real retval = 0;
        
        for(int i0=0;i0<nPoints;++i0)
        {
            // perturb the first vertex
            val1 += FD1_offset[order][i0]*eps;
            
            for(int i1=0;i1<nPoints;++i1)
            {
                // perturb
                val2 += FD1_offset[order][i1]*eps;
                mesh.updateDeformedConfiguration();
                
                retval += FD1_prefac[order][i0] * FD1_prefac[order][i1] * energyoperator.compute(mesh);
                
                // unperturb
                val2 -= FD1_offset[order][i1]*eps;
            }
            
            // unperturb
            val1 -= FD1_offset[order][i0]*eps;
        }
        
        // normalize one
        retval /= (eps*FD1_denum[order]);
        // normalize two
        retval /= (eps*FD1_denum[order]);
        
        return retval;
    }
};


template<typename tMesh>
struct ComputeHessian_Subset
{
    tMesh & mesh;
    const EnergyOperator<tMesh> & energyoperator;
    const int order;
    const Real eps;
    const std::vector< std::vector <Real>> FD1_prefac;
    const std::vector< std::vector <Real>> FD1_offset;
    const std::vector<Real> FD1_denum;
    const int nPoints;
    
    ComputeHessian_Subset(tMesh & mesh, const EnergyOperator<tMesh> & engOp, const int order, const Real eps):
    mesh(mesh),
    energyoperator(engOp),
    order(order),
    eps(eps),
    FD1_prefac({ {1, -1}, {1, -8, 8, -1}, {-1, 9, -45, 45, -9, 1}, {3, -32, 168, -672, 672, -168, 32, -3} }),
    FD1_offset({ {-1, 1}, {-2, -1, 1, 2}, {-3, -2, -1, 1, 2, 3}, {-4, -3, -2, -1, 1, 2, 3, 4} }),
    FD1_denum({2, 12, 60, 840}), // denumerator for first derivative
    nPoints(FD1_prefac[order].size())
    {
    }
    
    inline Real compute(const std::vector<int> & face_indices, Real & val1, Real & val2)
    {
        Real retval = 0;
        
        for(int i0=0;i0<nPoints;++i0)
        {
            // perturb the first vertex
            val1 += FD1_offset[order][i0]*eps;
            
            for(int i1=0;i1<nPoints;++i1)
            {
                // perturb
                val2 += FD1_offset[order][i1]*eps;
                mesh.updateDeformedConfiguration(face_indices);
                
                retval += FD1_prefac[order][i0] * FD1_prefac[order][i1] * energyoperator.computeSubsetEnergies(mesh, face_indices);
                
                // unperturb
                val2 -= FD1_offset[order][i1]*eps;
            }
            
            // unperturb
            val1 -= FD1_offset[order][i0]*eps;
        }
        
        // normalize one
        retval /= (eps*FD1_denum[order]);
        // normalize two
        retval /= (eps*FD1_denum[order]);
        
        // reset for the next call
        mesh.updateDeformedConfiguration(face_indices);
        
        return retval;
    }
};

template<typename tMesh>
class EnergyHessian_FD
{
protected:
    
    const Real FD_eps;
    const bool useSubset;
    
    mutable Profiler profiler;
    
    std::vector<int> getArrayIndicesWithoutFixedBoundaries(const tMesh & mesh) const;
    
    std::vector<std::set<int>> computeVertexVertexList(const tMesh & mesh) const;
    std::vector<std::set<int>> computeEdgeVertexList(const tMesh & mesh) const;
    std::vector<std::set<int>> computeEdgeEdgeList(const tMesh & mesh) const;
    
    inline size_t countListEntries(const std::vector<std::set<int>> & list) const
    {
        size_t retval=0;
        for(const auto & sublist : list)
            retval += sublist.size();
        return retval;
    }
    
    // basic interaction methods: compute energy for all faces
    int computeVertexVertexInteractions(tMesh & mesh, const EnergyOperator<tMesh> & energyoperator, const std::vector<std::set<int>> & vertex_vertex_list, const std::vector<int> & indicesWithoutFixedBoundaries, std::vector<Eigen::Triplet<Real>> & hessian_triplets) const;
    int computeEdgeVertexInteractions(tMesh & mesh, const EnergyOperator<tMesh> & energyoperator, const std::vector<std::set<int>> & edge_vertex_list, const std::vector<int> & indicesWithoutFixedBoundaries, std::vector<Eigen::Triplet<Real>> & hessian_triplets) const;
    int computeEdgeEdgeInteractions(tMesh & mesh, const EnergyOperator<tMesh> & energyoperator, const std::vector<std::set<int>> & edge_edge_list, const std::vector<int> & indicesWithoutFixedBoundaries, std::vector<Eigen::Triplet<Real>> & hessian_triplets) const;
    
    // faster interaction methods: compute energy only for faces that change
    int computeVertexVertexInteractionsSubset(tMesh & mesh, const EnergyOperator<tMesh> & energyoperator, const std::vector<std::set<int>> & vertex_vertex_list, const std::vector<int> & indicesWithoutFixedBoundaries, std::vector<Eigen::Triplet<Real>> & hessian_triplets) const;
    int computeEdgeVertexInteractionsSubset(tMesh & mesh, const EnergyOperator<tMesh> & energyoperator, const std::vector<std::set<int>> & edge_vertex_list, const std::vector<int> & indicesWithoutFixedBoundaries, std::vector<Eigen::Triplet<Real>> & hessian_triplets) const;
    int computeEdgeEdgeInteractionsSubset(tMesh & mesh, const EnergyOperator<tMesh> & energyoperator, const std::vector<std::set<int>> & edge_edge_list, const std::vector<int> & indicesWithoutFixedBoundaries, std::vector<Eigen::Triplet<Real>> & hessian_triplets) const;
    
    void compute(tMesh & mesh, const EnergyOperator<tMesh> & energyoperator);
public:
    
    EnergyHessian_FD(tMesh & mesh, const EnergyOperator<tMesh> & energyoperator, const Real FD_eps = 1e-3, const bool useSubset = true):
    FD_eps(FD_eps),
    useSubset(useSubset)
    {
        compute(mesh, energyoperator);
    }
    
    // public stuff
    std::vector<int> indicesWithoutFixedBoundaries;
    int numberOfNonFixedUnknowns;
    Eigen::SparseMatrix<Real> hessianMatrix;
    //

    
    void printProfilerSummary() const
    {
        profiler.printSummary();
    }    
};

#endif /* EnergyHessian_FD_hpp */
