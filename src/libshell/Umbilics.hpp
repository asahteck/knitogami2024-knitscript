#ifndef Umbilics_hpp
#define Umbilics_hpp

#include <iostream>

#include "common.hpp"

#ifdef USEVTK
#include "ReadVTK.hpp"
#include "WriteVTK.hpp"
#endif

// #include "WriteSTL.hpp"
#include "ComputeCurvatures.hpp"

#include <igl/triangle_triangle_adjacency.h>
#include <igl/boundary_loop.h>
#include <igl/ears.h>


// Struct to store information about a patch; all methods provided in class
// `UmbilicApproximation`.
struct UmbilicPatch {
    int centerTriangle;
    std::vector<int> triangles;
    int doubleIndex;

    // temporary calculated quantities
    std::set<int> edges;
};


template <typename tMesh>
struct centroidDistanceFunction {
    tMesh & mesh;
    const int centerTriangle;
    Eigen::Vector3d tCentroid;

    centroidDistanceFunction(tMesh & mesh, const int centerTriangle) :
    mesh(mesh),
    centerTriangle(centerTriangle)
    {
        tCentroid = mesh.getCurrentConfiguration().getTriangleInfoLite(mesh.getTopology(), centerTriangle).computeFaceCenter();
    }

    Real operator()(const int t) const {
        const Eigen::Vector3d tC1 = mesh.getCurrentConfiguration().getTriangleInfoLite(mesh.getTopology(), t).computeFaceCenter();
        return (tC1 - tCentroid).norm();
    }
};


template <typename tMesh>
class UmbilicApproximation
{
protected:
    tMesh & mesh;
    const Real patchSize;
    const bool verbose;

    std::vector<UmbilicPatch> umbilicPatches;

    // Pre-computed values
    Eigen::MatrixXi f2f;
    Eigen::MatrixXd pv1;
    Eigen::MatrixXd pv2;
    Eigen::VectorXd p1;
    Eigen::VectorXd p2;
    Eigen::VectorXd umbilicRank;

    // Compute all curvatures needed: principal curvatures k1 >= k2,
    // interpolated on vertices.
    void computeCurvatures() {
        const int nFaces = mesh.getNumberOfFaces();

        pv1.resize(nFaces, 3);
        pv2.resize(nFaces, 3);
        p1.resize(nFaces);
        p2.resize(nFaces);
        ComputeCurvatures<tMesh> computeCurvatures;
        computeCurvatures.computePrincipalVectors(mesh, pv1, pv2, false);
        computeCurvatures.computePrincipalCurvatures(mesh, p1, p2, false);

        umbilicRank = p1 - p2;
        assert((umbilicRank.array() >= 0.0).all()); // all values non-negative
    }

    // Generate the surrounding triangles of a patch, returning `true` if the
    // patch is umbilic.
    bool generatePatch(UmbilicPatch & patch, const int centerTriangle) {

        std::vector<int> triangleQueue;
        Eigen::VectorXd triangleQueueDistances(mesh.getNumberOfFaces());
        centroidDistanceFunction<tMesh> distance(mesh, centerTriangle);

        const auto f2e = mesh.getTopology().getFace2Edges();
        const auto f2v = mesh.getTopology().getFace2Vertices();

        patch.centerTriangle = centerTriangle;
        patch.triangles.push_back(centerTriangle);
        patch.edges.insert(f2e(centerTriangle, 0));
        patch.edges.insert(f2e(centerTriangle, 1));
        patch.edges.insert(f2e(centerTriangle, 2));

        // Initiate queue with neighbors of the center triangle.
        for (int j = 0; j < 3; j++) {
            int adjacentFace = f2f(centerTriangle, j);
            if (adjacentFace >= 0) {
                triangleQueue.push_back(adjacentFace);
                triangleQueueDistances(adjacentFace) = distance(adjacentFace);
            }
        }

        while (triangleQueue.size() > 0) {
            // Sort queue by increasing centroid distance to `centerTriangle`.
            std::sort(triangleQueue.begin(), triangleQueue.end(), [triangleQueueDistances](const int t1, const int t2) {
                return triangleQueueDistances(t1) < triangleQueueDistances(t2);
            });

            // Remove duplicate faces.
            triangleQueue.erase(std::unique(triangleQueue.begin(), triangleQueue.end()), triangleQueue.end());

            // Remove any faces already in the patch.
            for (int patchTriangle : patch.triangles) {
                triangleQueue.erase(std::remove(triangleQueue.begin(), triangleQueue.end(), patchTriangle), triangleQueue.end());
            }

            // Queue is emptied
            if (triangleQueue.size() == 0) {
                return false;
            }

            // Nearest triangle to `centerTriangle` (by centroid)
            int tAdd = triangleQueue[0];

            // No triangles can be added; we are done.
            if (triangleQueueDistances(tAdd) >= patchSize) {
                break;
            }

            // Adding this triangle would violate the patch's topological disk condition.
            if (!topologicalDisk(patch, tAdd)) {
                triangleQueue.erase(std::remove(triangleQueue.begin(), triangleQueue.end(), tAdd), triangleQueue.end());
                continue;
            }

            // Adding this triangle makes the patch non-umbilic.
            if (umbilicRank(tAdd) < umbilicRank(centerTriangle)) {
                return false;
            }

            // All conditions met; add this triangle and remove it from the queue.
            patch.triangles.push_back(tAdd);
            patch.edges.insert(f2e(tAdd, 0));
            patch.edges.insert(f2e(tAdd, 1));
            patch.edges.insert(f2e(tAdd, 2));
            triangleQueue.erase(std::remove(triangleQueue.begin(), triangleQueue.end(), tAdd), triangleQueue.end());

            // Add neighboring triangles to queue and continue.
            for (int j = 0; j < 3; j++) {
                int adjacentFace = f2f(tAdd, j);
                if (adjacentFace >= 0) {
                    triangleQueue.push_back(adjacentFace);
                    triangleQueueDistances(adjacentFace) = distance(adjacentFace);
                }
            }
        }

        // Remove all ears.
        Eigen::VectorXi earFaces;
        Eigen::VectorXi dummy;

        Eigen::MatrixXi patch_triangles = extractRows(patch.triangles, f2v);
        igl::ears(patch_triangles, earFaces, dummy);

        std::vector<int> removeFaces;
        for (int i = 0; i < earFaces.size(); i++) {
            removeFaces.push_back(patch.triangles[earFaces(i)]);
        }

        for (int i = 0; i < removeFaces.size(); i++) {
            patch.triangles.erase(std::remove(patch.triangles.begin(), patch.triangles.end(), removeFaces[i]), patch.triangles.end());
        }

        // Reset the edges
        patch.edges.clear();
        for (int i = 0; i < patch.triangles.size(); i++) {
            patch.edges.insert(f2e(patch.triangles[i], 0));
            patch.edges.insert(f2e(patch.triangles[i], 1));
            patch.edges.insert(f2e(patch.triangles[i], 2));
        }

        return true;
    }

    // Check whether adding triangle `tAdd` preserves the patch as a topological disk.
    bool topologicalDisk(const UmbilicPatch & patch, const int tAdd) {
        const auto f2e = mesh.getTopology().getFace2Edges();
        const auto f2v = mesh.getTopology().getFace2Vertices();
        const auto e2v = mesh.getTopology().getEdge2Vertices();

        int edge = -1;

        for (int i = 0; i < 3; i++) {
            std::set<int>::iterator it = patch.edges.find(f2e(tAdd, i));
            if (it != patch.edges.end()) {
                if (edge >= 0) { // two shared edges
                    return true;
                }

                edge = *it;
            }
        }

        if (edge < 0) { // no shared edges
            return false;
        }

        // Find the vertex opposite to the single shared edge.
        int oppVertex = -1;
        for (int i = 0; i < 3; i++) {
            if (f2v(tAdd, i) != e2v(edge, 0) && f2v(tAdd, i) != e2v(edge, 1)) {
                oppVertex = f2v(tAdd, i);
            }
        }

        // Ensure that this vertex is not part of the extant patch.
        for (int f : patch.triangles) {
            for (int i = 0; i < 3; i++) {
                if (f2v(f, i) == oppVertex) {
                    return false;
                }
            }
        }

        return true;
    }

    // Compute the index of an umbilic patch, using the angle rule and midpoint-to-midpoint segments on the outer triangles.
    void computeIndex(UmbilicPatch & patch) {
        if (patch.triangles.size() <= 2) {
            patch.doubleIndex = 0;
            return;
        }

        const auto e2v = mesh.getTopology().getEdge2Vertices();
        const auto e2f = mesh.getTopology().getEdge2Faces();
        const auto f2e = mesh.getTopology().getFace2Edges();
        const auto cvertices = mesh.getCurrentConfiguration().getVertices();

        std::vector<int> boundaryVertices;
        igl::boundary_loop(extractRows(patch.triangles, mesh.getTopology().getFace2Vertices()), boundaryVertices);

        auto isInSet = [](const int element, std::set<int> container) {
            return container.find(element) != container.end();
        };

        auto isInVector = [](const int element, std::vector<int> container) {
            return std::find(container.begin(), container.end(), element) != container.end();
        };

        // Projection operator onto `centerTriangle`
        Eigen::Vector3d centerTriangleNormal = (mesh.getCurrentConfiguration().getTriangleInfo(patch.centerTriangle).face_normal).normalized();
        Eigen::Matrix3d projectionOperator = Eigen::Matrix3d::Identity() - (centerTriangleNormal * centerTriangleNormal.transpose());

        // Use an edge (and its adjacent faces) iff it has one boundary and one
        // non-boundary vertex.
        std::set<int> quasiBoundaryEdges;
        int first = -1; // start with some edge
        for (auto e : patch.edges) {
            if (isInVector(e2v(e, 0), boundaryVertices) != isInVector(e2v(e, 1), boundaryVertices)) {
                quasiBoundaryEdges.insert(e);
                first = e;
            }
        }

        if (first < 0) {
            patch.doubleIndex = 0;
            return;
        }

        // Use orientation to figure out which face of `first` edge to use.
        bool b0 = isInVector(e2v(first, 0), boundaryVertices);
        Eigen::Vector3d vo = cvertices.row(e2v(first, b0 ? 0 : 1));
        Eigen::Vector3d vi = cvertices.row(e2v(first, b0 ? 1 : 0));
        Eigen::Vector3d vm = 0.5 * (vi + vo);
        Eigen::Vector3d vc = mesh.getCurrentConfiguration().getTriangleInfoLite(mesh.getTopology(), e2f(first, 0)).computeFaceCenter();
        Eigen::Vector3d faceNormal = mesh.getCurrentConfiguration().getTriangleInfo(e2f(first, 0)).face_normal;
        int face = e2f(first, ((vo - vi).cross(vc - vm).eval().dot(faceNormal) > 0.0) ? 0 : 1);

        Eigen::Vector3d unprojectedVector = pv1.row(face).eval();
        Eigen::Vector3d prevVector(0.0, 0.0, 0.0);
        Eigen::Vector3d currVector = projectionOperator * unprojectedVector;

        // Iterate through the quasi-boundary edges.
        Real integratedAngle = 0.0;

        int prev = first;
        int edgeIndex = which(prev, f2e.row(face));
        int curr = isInSet(f2e(face, (edgeIndex + 1) % 3), quasiBoundaryEdges) ? f2e(face, (edgeIndex + 1) % 3) : f2e(face, (edgeIndex + 2) % 3);

        do {
            // Align the principal curvature vectors using the acute angle rule,
            // projecting them onto the plane of `centerTriangle`. Skip the first
            // triangle.
            if (prevVector.norm() > 1e-6) {
                unprojectedVector = pv1.row(face).eval();
                Eigen::Vector3d unprojectedVector2 = pv2.row(face).eval();
                currVector = projectionOperator * unprojectedVector;
                Eigen::Vector3d currVector2 = projectionOperator * unprojectedVector2;

                if (currVector.dot(prevVector) < 0.0) { currVector = -currVector; }

                // Calculate the two angles and clamp to [-1, 1]
                Real cosAngle11 = currVector.dot(prevVector) / currVector.norm() / prevVector.norm();
                cosAngle11 = std::max(-1.0, std::min(1.0, cosAngle11));
                Real cosAngle12 = currVector2.dot(prevVector) / currVector2.norm() / prevVector.norm();
                cosAngle12 = std::max(-1.0, std::min(1.0, cosAngle12));

                Real sign;
                bool switch12 = std::abs(cosAngle12) > std::abs(cosAngle11);

                sign = ((prevVector.cross(currVector)).dot(faceNormal) > 0.0) ? 1.0 : -1.0;

                // Compute the angular difference between this curvature vector and the previous curvature vector angle.
                integratedAngle = integratedAngle + sign * std::acos(cosAngle11);
            }

            // Update iteration values
            prev = curr;
            face = e2f(prev, (which(face, e2f.row(prev)) + 1) % 2);
            edgeIndex = which(prev, f2e.row(face));
            curr = isInSet(f2e(face, (edgeIndex + 1) % 3), quasiBoundaryEdges) ? f2e(face, (edgeIndex + 1) % 3) : f2e(face, (edgeIndex + 2) % 3);

            prevVector = currVector;
        } while (prev != first);

        if (verbose) {
            std::cout << "integrated angle: " << integratedAngle / 2.0 / M_PI << std::endl;
        }

        patch.doubleIndex = (int)(std::round(2.0 * integratedAngle / 2.0 / M_PI));
    }

    // Compute umbilic patches of the mesh.
    void computeUmbilicPatches(const Real patchSize) {
        const int nFaces = mesh.getNumberOfFaces();
        for (int i = 0; i < nFaces; ++i) {
            UmbilicPatch patch;
            if (generatePatch(patch, i)) {
                computeIndex(patch);
                umbilicPatches.push_back(patch);

                if (verbose) {
                    std::cout << "Patch for triangle " << i << " has index " << 0.5 * (Real)(patch.doubleIndex) << std::endl;
                }
            }
        }
    }

public:
    UmbilicApproximation(tMesh & mesh, const Real patchSize, const bool verbose = false) :
    mesh(mesh),
    patchSize(patchSize),
    verbose(verbose)
    {
        // Pre-compute values
        igl::triangle_triangle_adjacency(mesh.getTopology().getFace2Vertices(), f2f);
        computeCurvatures();
    }

    Eigen::VectorXd umbilics()
    {
        computeUmbilicPatches(patchSize);

        Eigen::VectorXd umbilics_ = Eigen::VectorXd::Constant(mesh.getNumberOfFaces(), 0.0);
        for (int i = 0; i < umbilicPatches.size(); i++) {
            umbilics_(umbilicPatches[i].centerTriangle) = 0.5 * umbilicPatches[i].doubleIndex;
        }

        return umbilics_;
    }

    // Check: sum of umbilic indices should equal the Euler characteristic.
    bool checkEulerCharacteristic(const Eigen::VectorXd & umbilics_, const int EulerCharacteristic = 2) {
        return (int)(std::round(2 * umbilics_.sum())) == 2 * EulerCharacteristic;
    }
};

#endif /* Umbilics_hpp */
