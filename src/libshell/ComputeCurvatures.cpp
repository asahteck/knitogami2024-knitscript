//
//  ComputeCurvatures.cpp
//  Elasticity
//
//  Created by Wim van Rees on 5/20/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#include "ComputeCurvatures.hpp"
#include "ExtendedTriangleInfo.hpp"

template<typename tMesh>
void ComputeCurvatures<tMesh>::compute(const tMesh & mesh, Eigen::VectorXd & gauss, Eigen::VectorXd & mean) const
{
    // get current state quantities
    const int nFaces = mesh.getNumberOfFaces();
    const auto & currentState = mesh.getCurrentConfiguration();

    for(int i=0;i<nFaces;++i)
    {
        const ExtendedTriangleInfo & info = currentState.getTriangleInfo(i);

        const Eigen::Matrix2d aform = info.computeFirstFundamentalForm();
        const Eigen::Matrix2d bform = info.computeSecondFundamentalForm();

        const Eigen::Matrix2d shapeOp = aform.inverse() * bform;

        const Real gauss_curv = shapeOp.determinant();
        const Real mean_curv = 0.5*shapeOp.trace();

        gauss(i) = gauss_curv;
        mean(i) = mean_curv;
    }
}


template<typename tMesh>
void ComputeCurvatures<tMesh>::compute(const tMesh & mesh, Eigen::VectorXd & gauss, Eigen::VectorXd & mean,
                                       Eigen::MatrixXd & pCurvature1, Eigen::MatrixXd & pCurvature2, const bool orderByAbs) const
{
    // get current state quantities
    const int nFaces = mesh.getNumberOfFaces();
    const auto & currentState = mesh.getCurrentConfiguration();

    mean.resize(nFaces);
    gauss.resize(nFaces);
    pCurvature1.resize(nFaces, 3);
    pCurvature2.resize(nFaces, 3);

    for(int i=0;i<nFaces;++i)
    {
        const ExtendedTriangleInfo & info = currentState.getTriangleInfo(i);

        const Eigen::Matrix2d aform = info.computeFirstFundamentalForm();
        const Eigen::Matrix2d bform = info.computeSecondFundamentalForm();

        const Eigen::Matrix2d shapeOp = aform.inverse() * bform;

        const Real gauss_curv = shapeOp.determinant();
        const Real mean_curv = 0.5*shapeOp.trace();

        Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(shapeOp);
        Eigen::VectorXd evalues = eigenSolver.eigenvalues().real();
        Eigen::MatrixXd evectors = eigenSolver.eigenvectors().real();

        Eigen::MatrixXd basisVectors(3, 2);
        basisVectors << info.e1(0), info.e2(0),
                        info.e1(1), info.e2(1),
                        info.e1(2), info.e2(2);

        if (orderByAbs ? (std::abs(evalues(0)) > std::abs(evalues(1))) : (evalues(0) > evalues(1))) {
            pCurvature1.row(i) = (orderByAbs ? std::abs(evalues(0)) : evalues(0)) * (basisVectors * evectors.col(0)).normalized().eval();
            pCurvature2.row(i) = (orderByAbs ? std::abs(evalues(1)) : evalues(1)) * (basisVectors * evectors.col(1)).normalized().eval();
        }
        else {
            pCurvature1.row(i) = (orderByAbs ? std::abs(evalues(1)) : evalues(1)) * (basisVectors * evectors.col(1)).normalized().eval();
            pCurvature2.row(i) = (orderByAbs ? std::abs(evalues(0)) : evalues(0)) * (basisVectors * evectors.col(0)).normalized().eval();
        }

        gauss(i) = gauss_curv;
        mean(i) = mean_curv;
    }
}


template<typename tMesh>
void ComputeCurvatures<tMesh>::computePrincipalVectors(const tMesh & mesh, Eigen::MatrixXd & pCurvature1, Eigen::MatrixXd & pCurvature2, const bool orderByAbs) const
{
    // get current state quantities
    const int nFaces = mesh.getNumberOfFaces();
    const auto & currentState = mesh.getCurrentConfiguration();

    pCurvature1.resize(nFaces, 3);
    pCurvature2.resize(nFaces, 3);

    for(int i=0;i<nFaces;++i)
    {
        const ExtendedTriangleInfo & info = currentState.getTriangleInfo(i);

        const Eigen::Matrix2d aform = info.computeFirstFundamentalForm();
        const Eigen::Matrix2d bform = info.computeSecondFundamentalForm();

        const Eigen::Matrix2d shapeOp = aform.inverse() * bform;

        Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(shapeOp);
        Eigen::VectorXd evalues = eigenSolver.eigenvalues().real();
        Eigen::MatrixXd evectors = eigenSolver.eigenvectors().real();

        Eigen::MatrixXd basisVectors(3, 2);
        basisVectors << info.e1(0), info.e2(0),
                        info.e1(1), info.e2(1),
                        info.e1(2), info.e2(2);

        if (orderByAbs ? (std::abs(evalues(0)) > std::abs(evalues(1))) : (evalues(0) > evalues(1))) {
            pCurvature1.row(i) = (orderByAbs ? std::abs(evalues(0)) : evalues(0)) * (basisVectors * evectors.col(0)).normalized().eval();
            pCurvature2.row(i) = (orderByAbs ? std::abs(evalues(1)) : evalues(1)) * (basisVectors * evectors.col(1)).normalized().eval();
        }
        else {
            pCurvature1.row(i) = (orderByAbs ? std::abs(evalues(1)) : evalues(1)) * (basisVectors * evectors.col(1)).normalized().eval();
            pCurvature2.row(i) = (orderByAbs ? std::abs(evalues(0)) : evalues(0)) * (basisVectors * evectors.col(0)).normalized().eval();
        }
    }
}


template<typename tMesh>
void ComputeCurvatures<tMesh>::computePrincipalCurvatures(const tMesh & mesh, Eigen::VectorXd & pCurvature1, Eigen::VectorXd & pCurvature2, const bool orderByAbs) const
{
    // get current state quantities
    const int nFaces = mesh.getNumberOfFaces();
    const auto & currentState = mesh.getCurrentConfiguration();

    pCurvature1.resize(nFaces);
    pCurvature2.resize(nFaces);

    for(int i=0;i<nFaces;++i)
    {
        const ExtendedTriangleInfo & info = currentState.getTriangleInfo(i);

        const Eigen::Matrix2d aform = info.computeFirstFundamentalForm();
        const Eigen::Matrix2d bform = info.computeSecondFundamentalForm();

        const Eigen::Matrix2d shapeOp = aform.inverse() * bform;

        Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(shapeOp);
        Eigen::VectorXd evalues = eigenSolver.eigenvalues().real();

        if (orderByAbs ? (std::abs(evalues(0)) > std::abs(evalues(1))) : (evalues(0) > evalues(1))) {
            pCurvature1(i) = orderByAbs ? std::abs(evalues(0)) : evalues(0);
            pCurvature2(i) = orderByAbs ? std::abs(evalues(1)) : evalues(1);
        }
        else {
            pCurvature1(i) = orderByAbs ? std::abs(evalues(1)) : evalues(1);
            pCurvature2(i) = orderByAbs ? std::abs(evalues(0)) : evalues(0);
        }
    }
}


#include "Mesh.hpp"
#include "DynamicMesh.hpp"
template class ComputeCurvatures<Mesh>;
template class ComputeCurvatures<BilayerMesh>;
template class ComputeCurvatures<DynamicMesh>;
