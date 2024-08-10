//
//  EnergyHessian_FD.cpp
//  Elasticity
//
//  Created by Wim van Rees on 12/31/17.
//  Copyright © 2017 Wim van Rees. All rights reserved.
//

#include "EnergyHessian_FD.hpp"

template<typename tMesh>
void EnergyHessian_FD<tMesh>::compute(tMesh & mesh, const EnergyOperator<tMesh> & energyoperator)
{    
    // create the interaction lists
    profiler.push_start("list vertex-vertex");
    const std::vector<std::set<int>> vertex_vertex_list = computeVertexVertexList(mesh);
    profiler.pop_stop();
    
    // interaction list between vertices and edges : for each edge: all vertices for the attached faces plus the opposite vertices of the three edges of each of the two faces
    profiler.push_start("list edge-vertex");
    const std::vector<std::set<int>> edge_vertex_list = computeEdgeVertexList(mesh);
    profiler.pop_stop();
    
    // interaction list between edges and edges : for each edge: all three edges of both attached faces
    profiler.push_start("list edge-edge");
    const std::vector<std::set<int>> edge_edge_list = computeEdgeEdgeList(mesh);
    profiler.pop_stop();
    
    // count number of interactions
    profiler.push_start("count lists");
    const size_t vv_interactions = countListEntries(vertex_vertex_list);
    const size_t ve_interactions = countListEntries(edge_vertex_list);
    const size_t ee_interactions = countListEntries(edge_edge_list);
    profiler.pop_stop();
    
    // get the list of array indices
    profiler.push_start("compute non-fixed unknowns");
    indicesWithoutFixedBoundaries = getArrayIndicesWithoutFixedBoundaries(mesh);
    numberOfNonFixedUnknowns = std::count_if(indicesWithoutFixedBoundaries.begin(), indicesWithoutFixedBoundaries.end(), [](int i){return i>=0;});
    profiler.pop_stop();
    
    //    std::cout << "number of vv/ve/ee interactions = \t " << vv_interactions << "\t" << ve_interactions << "\t" << ee_interactions << std::endl;
    
    // allocate triplets : twice as much since we only need lower diagonal but ok
    std::vector<Eigen::Triplet<Real>> hessian_triplets;
    const int nTripletsApprox = vv_interactions + ve_interactions + ee_interactions;
    hessian_triplets.reserve(nTripletsApprox);
    
    // go through the interactions
    // can not do in parallel since mesh is shared
    // but energy computation should be in parallel
    
    // check if we will use the subset approach
    const bool shouldUseSubset = energyoperator.isSubsetImplemented() ? (useSubset ? true : false) : false;
    if(shouldUseSubset)
        std::cout << "EnergyHessian_FD : Using subset hessian energy computation" << std::endl;
    
    profiler.push_start("vertex-vertex");
    const int cnt_VV = shouldUseSubset ?
    computeVertexVertexInteractionsSubset(mesh, energyoperator, vertex_vertex_list, indicesWithoutFixedBoundaries, hessian_triplets) :
    computeVertexVertexInteractions(mesh, energyoperator, vertex_vertex_list, indicesWithoutFixedBoundaries, hessian_triplets);
    profiler.pop_stop();
    
    // edge-vertex
    profiler.push_start("edge-vertex");
    const int cnt_EV = shouldUseSubset ?
    computeEdgeVertexInteractionsSubset(mesh, energyoperator, edge_vertex_list, indicesWithoutFixedBoundaries, hessian_triplets) :
    computeEdgeVertexInteractions(mesh, energyoperator, edge_vertex_list, indicesWithoutFixedBoundaries, hessian_triplets);
    profiler.pop_stop();
    
    // edge-edge
    profiler.push_start("edge-edge");
    const int cnt_EE = shouldUseSubset ?
    computeEdgeEdgeInteractionsSubset(mesh, energyoperator, edge_edge_list, indicesWithoutFixedBoundaries, hessian_triplets) :
    computeEdgeEdgeInteractions(mesh, energyoperator, edge_edge_list, indicesWithoutFixedBoundaries, hessian_triplets);
    profiler.pop_stop();
    
    std::cout << "number of computed interactions = \t " << cnt_VV << "\t" << cnt_EV << "\t" << cnt_EE << std::endl;
    hessianMatrix.resize(numberOfNonFixedUnknowns, numberOfNonFixedUnknowns);
    hessianMatrix.setFromTriplets(hessian_triplets.begin(), hessian_triplets.end());
    
    hessianMatrix.prune(1e-12);
    

}

template<typename tMesh>
std::vector<std::set<int>> EnergyHessian_FD<tMesh>::computeVertexVertexList(const tMesh & mesh) const
{
    const int nVertices = mesh.getNumberOfVertices();
    const auto face2edges = mesh.getTopology().getFace2Edges();
    const auto face2vertices = mesh.getTopology().getFace2Vertices();
    const auto edge2vertices = mesh.getTopology().getEdge2Vertices();
    const auto edge2faces = mesh.getTopology().getEdge2Faces();
    const auto edge2oppositevertices = mesh.getTopology().getEdge2OppositeVertices();
    const auto & vertex2faces = mesh.getTopology().getVertex2Faces();
    
    std::vector<std::set<int>> vertex_vertex_list(nVertices);
    
    // first tier = all the vertices that share a face with my vertex :
    // basically : need vertex2edges and the edge2vertices of the other vertex on the edge
    // (or vertex2faces and then face2vertices)
    
    for(int v0_idx=0;v0_idx<nVertices;++v0_idx)
    {
        const std::vector<int> & v0_faces = vertex2faces[v0_idx];
        
        for(size_t i=0;i<v0_faces.size();++i)
        {
            const int f_idx = v0_faces[i];
            
            // insert all three indices of this face (also my own : set will take care of duplicates)
            vertex_vertex_list[v0_idx].insert(face2vertices(f_idx,0));
            vertex_vertex_list[v0_idx].insert(face2vertices(f_idx,1));
            vertex_vertex_list[v0_idx].insert(face2vertices(f_idx,2));
            
            
            // the faces for which we are the 'opposite' vertex will matter : anything that changes the energy of those faces will affect us.
            // what will change the energy of those faces is :
            //  - their vertices (we already have them above)
            //  - their edge directors
            //  - their own opposite vertices
            // here we deal with the last one : all opposite vertices of all faces for which we are the opposite vertex
            
            
            // find the edge which is opposite to us : the one that we are not attached to
            for(int e_idx_rel=0; e_idx_rel<3; ++e_idx_rel)
            {
                const int e_idx = face2edges(f_idx, e_idx_rel);
                if((edge2vertices(e_idx,0) == v0_idx) or (edge2vertices(e_idx,1) == v0_idx)) continue;
                
                // now we know this e_idx is opposite to us : we are one of its opposite vertices
                assert((edge2oppositevertices(e_idx,0) == v0_idx) or (edge2oppositevertices(e_idx,1) == v0_idx));
                assert((edge2faces(e_idx,0) == f_idx) or (edge2faces(e_idx,1) == f_idx));
                
                // find the other face
                const int f_idx_other = (edge2faces(e_idx,0) == f_idx) ? edge2faces(e_idx,1) : edge2faces(e_idx,0);
                if(f_idx_other < 0) continue; // in case the edge is on the boundary
                
                // ok finally we are in business: we add all 6 opposite vertices to this edge (3 of them we already have : those two that share the edge with f_idx and myself (v0_idx)
                for(int d=0;d<3;++d)
                    for(int j=0;j<2;++j)
                    {
                        const int v_idx_opp = edge2oppositevertices(face2edges(f_idx_other,d),j);
                        if(v_idx_opp >= 0) vertex_vertex_list[v0_idx].insert(v_idx_opp);
                    }
            }
        }
    }
    
    return vertex_vertex_list;
}

template<typename tMesh>
std::vector<std::set<int>> EnergyHessian_FD<tMesh>::computeEdgeVertexList(const tMesh & mesh) const
{
    const int nEdges = mesh.getNumberOfEdges();
    
    const auto face2edges = mesh.getTopology().getFace2Edges();
    const auto face2vertices = mesh.getTopology().getFace2Vertices();
    const auto edge2faces = mesh.getTopology().getEdge2Faces();
    const auto edge2oppositevertices = mesh.getTopology().getEdge2OppositeVertices();
    
    std::vector<std::set<int>> edge_vertex_list(nEdges);
    for(int idx_e=0;idx_e<nEdges;++idx_e)
    {
        // loop over the faces adjacent to this edge
        for(int idx_f_rel=0;idx_f_rel<2;++idx_f_rel)
        {
            const int idx_f = edge2faces(idx_e,idx_f_rel);
            if(idx_f >= 0) // only consider if there is actually a face
            {
                // loop over the vertices of this face
                for(int d=0;d<3;++d)
                {
                    // add the direct vertices
                    edge_vertex_list[idx_e].insert(face2vertices(idx_f,d));
                    
                    // add the opposite vertices : add both and let the set figure out which one is duplicated
                    // (now d is used as edge index)
                    const int idx_v_opp_0 = edge2oppositevertices(face2edges(idx_f, d),0);
                    const int idx_v_opp_1 = edge2oppositevertices(face2edges(idx_f, d),1);
                    if(idx_v_opp_0 >= 0) edge_vertex_list[idx_e].insert(idx_v_opp_0);
                    if(idx_v_opp_1 >= 0) edge_vertex_list[idx_e].insert(idx_v_opp_1);
                }
            }
        }
    }
    
    return edge_vertex_list;
    
}

template<typename tMesh>
std::vector<std::set<int>> EnergyHessian_FD<tMesh>::computeEdgeEdgeList(const tMesh & mesh) const
{
    const int nEdges = mesh.getNumberOfEdges();
    
    const auto face2edges = mesh.getTopology().getFace2Edges();
    const auto edge2faces = mesh.getTopology().getEdge2Faces();
    
    std::vector<std::set<int>> edge_edge_list(nEdges);
    for(int idx_e=0;idx_e<nEdges;++idx_e)
        // loop over the faces adjacent to this edge
        for(int idx_f_rel=0;idx_f_rel<2;++idx_f_rel)
        {
            const int idx_f = edge2faces(idx_e,idx_f_rel);
            if(idx_f >= 0) // only consider if there is actually a face
                for(int d=0;d<3;++d)
                {
                    const int idx_e_other = face2edges(idx_f, d);
                    edge_edge_list[idx_e].insert(idx_e_other);
                }
        }
    
    return edge_edge_list;
    
}


template<typename tMesh>
int EnergyHessian_FD<tMesh>::computeVertexVertexInteractions(tMesh & mesh, const EnergyOperator<tMesh> & energyoperator, const std::vector<std::set<int>> & vertex_vertex_list, const std::vector<int> & indicesWithoutFixedBoundaries, std::vector<Eigen::Triplet<Real>> & hessian_triplets) const
{
    ComputeHessian<tMesh> computeHessian(mesh, energyoperator, 1, FD_eps);
    
    const int nVertices = mesh.getNumberOfVertices();
    const auto vertices_bc = mesh.getBoundaryConditions().getVertexBoundaryConditions();
    
    auto vertices = mesh.getCurrentConfiguration().getVertices();
    
    int cnt_VV = 0;
    // vertex-vertex
    for(int idx_v0=0;idx_v0<nVertices;++idx_v0)
        for(int d0=0;d0<3;++d0)
        {
            const int hessIdx_0 = indicesWithoutFixedBoundaries[idx_v0 + nVertices*d0];
            if(vertices_bc(idx_v0,d0)) continue; // fixed: zero derivative
            
            for(const int idx_v1 : vertex_vertex_list[idx_v0])
                for(int d1=0;d1<3;++d1)
                {
                    const int hessIdx_1 = indicesWithoutFixedBoundaries[idx_v1 + nVertices*d1]; // col-major
                    if(hessIdx_1 > hessIdx_0) continue; // only need lower triangular
                    if(vertices_bc(idx_v1,d1)) continue; // fixed: zero derivative
                    
                    // compute
                    const Real hessianFD = computeHessian.compute(vertices(idx_v0,d0), vertices(idx_v1,d1));
                    
                    // store
                    hessian_triplets.push_back(Eigen::Triplet<Real>(hessIdx_0, hessIdx_1, hessianFD));
                    
                    cnt_VV++;
                }
        }
    return cnt_VV;
}

template<typename tMesh>
int EnergyHessian_FD<tMesh>::computeEdgeVertexInteractions(tMesh & mesh, const EnergyOperator<tMesh> & energyoperator, const std::vector<std::set<int>> & edge_vertex_list, const std::vector<int> & indicesWithoutFixedBoundaries, std::vector<Eigen::Triplet<Real>> & hessian_triplets) const
{
    ComputeHessian<tMesh> computeHessian(mesh, energyoperator, 1, FD_eps);
    
    const int nVertices = mesh.getNumberOfVertices();
    const int nEdges = mesh.getNumberOfEdges();
    
    const auto vertices_bc = mesh.getBoundaryConditions().getVertexBoundaryConditions();
    const auto edges_bc = mesh.getBoundaryConditions().getEdgeBoundaryConditions();
    
    auto vertices = mesh.getCurrentConfiguration().getVertices();
    auto edgedirectors = mesh.getCurrentConfiguration().getEdgeDirectors();
    
    int cnt_EV = 0;
    for(int idx_e0=0;idx_e0<nEdges;++idx_e0)
    {
        const int hessIdx_0 = indicesWithoutFixedBoundaries[idx_e0 + 3*nVertices]; // col-major
        if(edges_bc(idx_e0)) continue; // fixed: zero derivative
        
        for(const int idx_v1 : edge_vertex_list[idx_e0])
            for(int d1=0;d1<3;++d1)
            {
                const int hessIdx_1 = indicesWithoutFixedBoundaries[idx_v1 + nVertices*d1]; // col-major
                // edge is always in the lower triangular part wrt vertex index
                if(vertices_bc(idx_v1,d1)) continue; // fixed: zero derivative
                
                // compute
                const Real hessianFD = computeHessian.compute(edgedirectors(idx_e0), vertices(idx_v1,d1));
                
                // store
                hessian_triplets.push_back(Eigen::Triplet<Real>(hessIdx_0, hessIdx_1, hessianFD));
                
                cnt_EV++;
            }
    }
    
    return cnt_EV;
}

template<typename tMesh>
int EnergyHessian_FD<tMesh>::computeEdgeEdgeInteractions(tMesh & mesh, const EnergyOperator<tMesh> & energyoperator, const std::vector<std::set<int>> & edge_edge_list, const std::vector<int> & indicesWithoutFixedBoundaries, std::vector<Eigen::Triplet<Real>> & hessian_triplets) const
{
    ComputeHessian<tMesh> computeHessian(mesh, energyoperator, 1, FD_eps);
    
    const int nVertices = mesh.getNumberOfVertices();
    const int nEdges = mesh.getNumberOfEdges();
    
    const auto edges_bc = mesh.getBoundaryConditions().getEdgeBoundaryConditions();
    
    auto edgedirectors = mesh.getCurrentConfiguration().getEdgeDirectors();
    
    int cnt_EE = 0;
    for(int idx_e0=0;idx_e0<nEdges;++idx_e0)
    {
        const int hessIdx_0 = indicesWithoutFixedBoundaries[idx_e0 + 3*nVertices]; // col-major
        if(edges_bc(idx_e0)) continue; // fixed: zero derivative
        
        for(const int idx_e1 : edge_edge_list[idx_e0])
        {
            const int hessIdx_1 = indicesWithoutFixedBoundaries[idx_e1 + nVertices*3]; // col-major
            if(hessIdx_1 > hessIdx_0) continue; // only need lower triangular
            if(edges_bc(idx_e1)) continue; // fixed: zero derivative
            
            // compute
            const Real hessianFD = computeHessian.compute(edgedirectors(idx_e0), edgedirectors(idx_e1));
            
            // store
            hessian_triplets.push_back(Eigen::Triplet<Real>(hessIdx_0, hessIdx_1, hessianFD));
            
            cnt_EE++;
        }
    }
    
    return cnt_EE;
}


template<typename tMesh>
int EnergyHessian_FD<tMesh>::computeVertexVertexInteractionsSubset(tMesh & mesh, const EnergyOperator<tMesh> & energyoperator, const std::vector<std::set<int>> & vertex_vertex_list, const std::vector<int> & indicesWithoutFixedBoundaries, std::vector<Eigen::Triplet<Real>> & hessian_triplets) const
{
    ComputeHessian_Subset<tMesh> computeHessian(mesh, energyoperator, 1, FD_eps);
    
    const int nVertices = mesh.getNumberOfVertices();
    const auto vertices_bc = mesh.getBoundaryConditions().getVertexBoundaryConditions();
    
    const auto & vertex2faces = mesh.getTopology().getVertex2Faces();
    const auto face2edges = mesh.getTopology().getFace2Edges();
    const auto edge2faces = mesh.getTopology().getEdge2Faces();
    
    auto vertices = mesh.getCurrentConfiguration().getVertices();
    
    int cnt_VV = 0;
    // vertex-vertex
    for(int idx_v0=0;idx_v0<nVertices;++idx_v0)
    {
        // create list of faces we care about : use a set to avoid duplicates
        std::set<int> v0_faces;
        
        {
            const std::vector<int> & v0_direct_faces = vertex2faces[idx_v0];
            
            for(const auto & faceidx : v0_direct_faces)
            {
                if(faceidx < 0) continue;
                
                // first tier: all faces adjacent to the vertex
                v0_faces.insert(faceidx);
                
                for(int e0=0;e0<3;++e0)
                    for(int f0=0;f0<2;++f0)
                    {
                        // second tier : all faces adjacent to the faces adjacent to the vertex
                        const int adj_faceidx = edge2faces(face2edges(faceidx, e0), f0);
                        if(adj_faceidx < 0) continue;
                        
                        v0_faces.insert(adj_faceidx);
                        
                        // third tier : all faces adjacent to the face adjacent to the face adjacent to the vertex (to get two 'opposite' vertices)
                        for(int e1=0;e1<3;++e1)
                            for(int f1=0;f1<2;++f1)
                            {
                                const int adj_adj_faceidx = edge2faces(face2edges(adj_faceidx, e1), f1);
                                if(adj_adj_faceidx >=0) v0_faces.insert(adj_adj_faceidx);
                            }
                    }
            }
        }
        
        // convert back to vector
        const std::vector<int> v0_subset( v0_faces.begin(), v0_faces.end() );
        
        for(int d0=0;d0<3;++d0)
        {
            const int hessIdx_0 = indicesWithoutFixedBoundaries[idx_v0 + nVertices*d0]; // col-major
            if(vertices_bc(idx_v0,d0)) continue; // fixed: zero derivative
            
            for(const int idx_v1 : vertex_vertex_list[idx_v0])
                for(int d1=0;d1<3;++d1)
                {
                    const int hessIdx_1 = indicesWithoutFixedBoundaries[idx_v1 + nVertices*d1]; // col-major
                    if(hessIdx_1 > hessIdx_0) continue; // only need lower triangular
                    if(vertices_bc(idx_v1,d1)) continue; // fixed: zero derivative
                    
                    // compute
                    const Real hessianFD = computeHessian.compute(v0_subset, vertices(idx_v0,d0), vertices(idx_v1,d1));
                    
                    // store
                    hessian_triplets.push_back(Eigen::Triplet<Real>(hessIdx_0, hessIdx_1, hessianFD));
                    
                    cnt_VV++;
                }
        }
    }
    return cnt_VV;
}


template<typename tMesh>
int EnergyHessian_FD<tMesh>::computeEdgeVertexInteractionsSubset(tMesh & mesh, const EnergyOperator<tMesh> & energyoperator, const std::vector<std::set<int>> & edge_vertex_list, const std::vector<int> & indicesWithoutFixedBoundaries, std::vector<Eigen::Triplet<Real>> & hessian_triplets) const
{
    ComputeHessian_Subset<tMesh> computeHessian(mesh, energyoperator, 1, FD_eps);
    
    const int nVertices = mesh.getNumberOfVertices();
    const int nEdges = mesh.getNumberOfEdges();
    
    const auto face2edges = mesh.getTopology().getFace2Edges();
    const auto edge2faces = mesh.getTopology().getEdge2Faces();
    const auto vertices_bc = mesh.getBoundaryConditions().getVertexBoundaryConditions();
    const auto edges_bc = mesh.getBoundaryConditions().getEdgeBoundaryConditions();
    
    auto vertices = mesh.getCurrentConfiguration().getVertices();
    auto edgedirectors = mesh.getCurrentConfiguration().getEdgeDirectors();
    
    int cnt_EV = 0;
    for(int idx_e0=0;idx_e0<nEdges;++idx_e0)
    {
        const int hessIdx_0 = indicesWithoutFixedBoundaries[idx_e0 + 3*nVertices]; // col-major
        if(edges_bc(idx_e0)) continue; // fixed: zero derivative
        
        // create list of faces we care about : use a set to avoid duplicates
        std::set<int> e0_faces;
        
        {
            const std::vector<int> e0_direct_faces = {edge2faces(idx_e0,0), edge2faces(idx_e0,1)};
            
            for(const auto & faceidx : e0_direct_faces)
            {
                if(faceidx < 0) continue;
                
                // first tier: all faces adjacent to the edge
                e0_faces.insert(faceidx);
                
                for(int e0=0;e0<3;++e0)
                    for(int f0=0;f0<2;++f0)
                    {
                        // second tier : all faces adjacent to the faces adjacent to the edge
                        const int adj_faceidx = edge2faces(face2edges(faceidx, e0), f0);
                        if(adj_faceidx >= 0) e0_faces.insert(adj_faceidx);
                    }
            }
        }
        
        // convert to vector
        const std::vector<int> e0_subset( e0_faces.begin(), e0_faces.end() );
        
        for(const int idx_v1 : edge_vertex_list[idx_e0])
            for(int d1=0;d1<3;++d1)
            {
                const int hessIdx_1 = indicesWithoutFixedBoundaries[idx_v1 + nVertices*d1]; // col-major
                // edge is always in the lower triangular part wrt vertex index
                if(vertices_bc(idx_v1,d1)) continue; // fixed: zero derivative
                
                // compute
                const Real hessianFD = computeHessian.compute(e0_subset, edgedirectors(idx_e0), vertices(idx_v1,d1));
                
                // store
                hessian_triplets.push_back(Eigen::Triplet<Real>(hessIdx_0, hessIdx_1, hessianFD));
                
                cnt_EV++;
            }
    }
    
    return cnt_EV;
}

template<typename tMesh>
int EnergyHessian_FD<tMesh>::computeEdgeEdgeInteractionsSubset(tMesh & mesh, const EnergyOperator<tMesh> & energyoperator, const std::vector<std::set<int>> & edge_edge_list, const std::vector<int> & indicesWithoutFixedBoundaries, std::vector<Eigen::Triplet<Real>> & hessian_triplets) const
{
    ComputeHessian_Subset<tMesh> computeHessian(mesh, energyoperator, 1, FD_eps);
    
    const int nVertices = mesh.getNumberOfVertices();
    const int nEdges = mesh.getNumberOfEdges();
    
    const auto face2edges = mesh.getTopology().getFace2Edges();
    const auto edge2faces = mesh.getTopology().getEdge2Faces();
    
    const auto edges_bc = mesh.getBoundaryConditions().getEdgeBoundaryConditions();
    
    auto edgedirectors = mesh.getCurrentConfiguration().getEdgeDirectors();
    
    int cnt_EE = 0;
    for(int idx_e0=0;idx_e0<nEdges;++idx_e0)
    {
        const int hessIdx_0 = indicesWithoutFixedBoundaries[idx_e0 + 3*nVertices]; // col-major
        if(edges_bc(idx_e0)) continue; // fixed: zero derivative
        
        // create list of faces we care about : use a set to avoid duplicates
        std::set<int> e0_faces;
        {
            const std::vector<int> e0_direct_faces = {edge2faces(idx_e0,0), edge2faces(idx_e0,1)};
            
            for(const auto & faceidx : e0_direct_faces)
            {
                // first tier: all faces adjacent to the edge : that's it
                if(faceidx >= 0) e0_faces.insert(faceidx);
            }
        }
        
        // convert to vector
        const std::vector<int> e0_subset( e0_faces.begin(), e0_faces.end() );
        
        for(const int idx_e1 : edge_edge_list[idx_e0])
        {
            const int hessIdx_1 = indicesWithoutFixedBoundaries[idx_e1 + nVertices*3]; // col-major
            if(hessIdx_1 > hessIdx_0) continue; // only need lower triangular
            if(edges_bc(idx_e1)) continue; // fixed: zero derivative
            
            // compute
            const Real hessianFD = computeHessian.compute(e0_subset, edgedirectors(idx_e0), edgedirectors(idx_e1));
            
            // store
            hessian_triplets.push_back(Eigen::Triplet<Real>(hessIdx_0, hessIdx_1, hessianFD));
            
            cnt_EE++;
        }
    }
    
    return cnt_EE;
}

template<typename tMesh>
std::vector<int> EnergyHessian_FD<tMesh>::getArrayIndicesWithoutFixedBoundaries(const tMesh & mesh) const
{
    // what we want is an array of consecutive vertex and edge indices that are free to move
    // if we have the following (vertices and BC)
    // v0_x, v0_y, v0_z, v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, e0  , e1  , e2
    // true, fals, fals, fals, fals, true, true, fals, fals, fals, true, fals
    
    // we want the array to be
    // unk:  (v0_x, v0_y, v0_z, v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, e0  , e1  , e2)
    
    // idx:   0  , 1   , 2   , 3   , 4   ,  5  ,  6  , 7   , 8   , 9   , 10  , 11
    // arr : -1  , 0   , 1   , 2   , 3   , -1  , -1  , 4   , 5   , 6   , -1  , 7
    
    
    const int nVertices = mesh.getNumberOfVertices();
    const int nEdges = mesh.getNumberOfEdges();
    
    const auto vertices_bc = mesh.getBoundaryConditions().getVertexBoundaryConditions();
    const auto edges_bc = mesh.getBoundaryConditions().getEdgeBoundaryConditions();
    
    std::vector<int> retval(3*nVertices + nEdges, -1);
    
    int cntr = 0;
    for(int d=0;d<3;++d) // loop over d first to make it consistent with column major form
        for(int i=0;i<nVertices;++i)
        {
            if(vertices_bc(i,d)) continue;
            
            retval[i + d*nVertices] = cntr;
            cntr++;
        }
    
    for(int i=0;i<nEdges;++i)
    {
        if(edges_bc(i)) continue;
        
        retval[3*nVertices + i] = cntr;
        cntr++;
    }
    
    return retval;
}

#include "DynamicMesh.hpp"
template class EnergyHessian_FD<DynamicMesh>;
