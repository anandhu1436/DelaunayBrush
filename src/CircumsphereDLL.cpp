#include "CircumsphereDLL.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <igl/copyleft/marching_cubes.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/barycentric_coordinates.h>
#include <igl/cotmatrix.h>
#include <igl/harmonic.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<int, K> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Tetrahedron_3 Tetrahedron;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef CGAL::Surface_mesh<Point> SurfaceMesh;

double compute_angle(const Vector& v1, const Vector& v2) {
    double dot_product = v1 * v2;
    double magnitude_v1 = std::sqrt(v1.squared_length());
    double magnitude_v2 = std::sqrt(v2.squared_length());
    return std::acos(dot_product / (magnitude_v1 * magnitude_v2));
}

Eigen::MatrixXd unorderedSetToEigenMatrix(const std::unordered_set<Point>& unique_vertices) {
    int n = static_cast<int>(unique_vertices.size());
    Eigen::MatrixXd V(n, 3);
    
    int i = 0;
    for (const auto& p : unique_vertices) {
        V(i, 0) = CGAL::to_double(p.x());
        V(i, 1) = CGAL::to_double(p.y());
        V(i, 2) = CGAL::to_double(p.z());
        ++i;
    }

    return V;
}


void uniform_remesh(
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    double target_edge_length,
    int iterations = 3,
    int max_vertices = 5000)
{
    // Convert Eigen to CGAL
    clock_t start, stop;
    start = clock();

    SurfaceMesh mesh;
    std::vector<SurfaceMesh::Vertex_index> vtx_indices;
    for (int i = 0; i < V.rows(); ++i)
        vtx_indices.push_back(mesh.add_vertex(Point(V(i, 0), V(i, 1), V(i, 2))));

    for (int i = 0; i < F.rows(); ++i)
        mesh.add_face(vtx_indices[F(i, 0)], vtx_indices[F(i, 1)], vtx_indices[F(i, 2)]);
    stop = clock();

    double time = double (stop - start) / CLOCKS_PER_SEC;

    std::cout << "time for converting mesh to CGAL : " <<  time << " seconds" << std::endl;

    // Iterative remeshing with a vertex count cap
    start = clock();
    for (int i = 0; i < iterations; ++i) {
        PMP::isotropic_remeshing(
            faces(mesh),
            target_edge_length,
            mesh,
            PMP::parameters::number_of_iterations(1).protect_constraints(false)
        );

        if (mesh.number_of_vertices() >= max_vertices) {
            std::cout << "Reached vertex limit (" << mesh.number_of_vertices() << "), stopping remeshing early.\n";
            break;
        }
    }
    stop = clock();
    time = double (stop - start) / CLOCKS_PER_SEC;
    std::cout << "time for remeshing : " <<  time << " seconds" << std::endl;



    // Convert back to Eigen
    start = clock();
    V.resize(mesh.number_of_vertices(), 3);
    std::map<SurfaceMesh::Vertex_index, int> vi_map;
    int vi = 0;
    for (auto v : mesh.vertices()) {
        Point p = mesh.point(v);
        V.row(vi) = Eigen::RowVector3d(p.x(), p.y(), p.z());
        vi_map[v] = vi++;
    }

    F.resize(mesh.number_of_faces(), 3);
    int fi = 0;
    for (auto f : mesh.faces()) {
        int j = 0;
        for (auto v : vertices_around_face(mesh.halfedge(f), mesh))
            F(fi, j++) = vi_map[v];
        fi++;
    }
    stop = clock();
    time = double (stop - start) / CLOCKS_PER_SEC;
    std::cout << "time for converting mesh back from CGAL : " << time << " seconds" << std::endl;


}


void sparseLSsolve(Eigen::SparseMatrix<double>& A, Eigen::MatrixXd& b, Eigen::MatrixXd& X, bool factorize) {
    static Eigen::MatrixXd previousSolution;
    static bool firstSolving = true;

    std::cout << "Rows in matrix : " << A.rows() << std::endl;
    #define LSCG
    #ifndef LSCG
        static Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
        // Solve the constraints system with least-squares
        clock_t start = clock();
        if (factorize) {
            A.makeCompressed();
            solver.compute(A);
        }
        clock_t end = clock();
        double time = double (end - start) / CLOCKS_PER_SEC;
        std::cout << "Factorization time : " << time << " seconds" << std::endl;

        start = clock();
        X = solver.solve(b);
        end = clock();
        time = double (end - start) / CLOCKS_PER_SEC;
        std::cout << "Solving time : " << time << " seconds" << std::endl;
        std::cout << std::endl;

    #else
        Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver;

        clock_t start = clock();
        solver.compute(A);
        clock_t end = clock();
        double time = double (end - start) / CLOCKS_PER_SEC;
        std::cout << "Preconditioning time : " << time << " seconds" << std::endl;

        // solver.setTolerance(0.0001);
        solver.setTolerance(0.001);
        // solver.setMaxIterations( 100);

        start = clock();
        firstSolving = true;
        if (firstSolving) {
            printf("=!=!=!=!=!=!=!=!=!=!=!=!= RECALCUL DEPUIS DEBUT\n");
            // X = solver.solve(b);
            X = solver.solveWithGuess( b, X);
            firstSolving = false;
            // firstSolving = true;
            // std::cout << X << std::endl;
        }
        previousSolution = X;
        end = clock();
        time = double (end - start) / CLOCKS_PER_SEC;
        std::cout << "time for solving : " << time << " seconds" << std::endl;

        std::cout << "Number of iterations for solving : " << solver.iterations() << std::endl;
        std::cout << std::endl;
    #endif
}



void sparseVerticalConcat(
    Eigen::SparseMatrix<double> M1,
    Eigen::SparseMatrix<double> M2,
    Eigen::SparseMatrix<double>& M) {

    if (M1.cols() != M2.cols()) {
        std::cerr << "Concatenated matrices must have the same number of columns."
        << std::endl;
        return;
    }

    // Sparse matrix concatenation tip taken from 
    // https://stackoverflow.com/questions/41756428/concatenate-sparse-matrix-eigen
    M = Eigen::SparseMatrix<double>(M1.rows() + M2.rows(), M1.cols());
    // Create a list of triplets storing non zero coordinates of the matrix A
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(M1.nonZeros() + M2.nonZeros());

    // Fill with M1 part
    for (int k = 0; k < M1.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(M1, k); it; ++it) {
            triplets.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
        }
    }
    // Fill with M2 part
    for (int k = 0; k < M2.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(M2, k); it; ++it) {
            triplets.push_back(Eigen::Triplet<double>(M1.rows() + it.row(), it.col(), it.value()));
        }
    }
    M.setFromTriplets(triplets.begin(), triplets.end());
}

void computeFromOctahedronLapAndBary(double w, const Eigen::MatrixXd& inputPts,
                                      const Eigen::MatrixXd& VP,
                                      const Eigen::MatrixXi& FP,
                                      Eigen::MatrixXd& Vcl) {
    // Step 1: Project input points onto the mesh
    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    clock_t start = clock();
    igl::point_mesh_squared_distance(inputPts, VP, FP, sqrD, I, C);
    clock_t stop = clock();
    double time = double (stop - start) / CLOCKS_PER_SEC;
    std::cout << "time for projecting points on proxy : " <<  time << " seconds" << std::endl;

    const int n = inputPts.rows();
    const int n_vertices = VP.rows();

    start = clock();
    // Step 2: Build barycentric constraint matrix
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n * 3); // Each point contributes 3 non-zero entries
    Eigen::MatrixXd bc = inputPts; // Our positional constraints

    for (int i = 0; i < n; ++i) {
        const int face_idx = I[i];
        const Eigen::RowVector3i face = FP.row(face_idx);
    
        // Get triangle vertex positions
        Eigen::RowVector3d A = VP.row(face(0));
        Eigen::RowVector3d B = VP.row(face(1));
        Eigen::RowVector3d C_vert = VP.row(face(2));
    
        Eigen::RowVector3d p = C.row(i); // Closest point on the triangle
    
        // Compute barycentric coordinates
        Eigen::RowVector3d bary_coords;
        igl::barycentric_coordinates(p, A, B, C_vert, bary_coords);
    
        // Optionally normalize if numerical instability arises
        bary_coords /= bary_coords.sum();
    
        // Fill constraint matrix triplets
        for (int j = 0; j < 3; ++j) {
            triplets.emplace_back(i, face(j), bary_coords(j));
        }
    }
    

    // Build sparse barycentric constraint matrix
    Eigen::SparseMatrix<double> Ac(n, n_vertices);
    Ac.setFromTriplets(triplets.begin(), triplets.end());
    Ac.makeCompressed();

    // Step 3: Compute cotangent Laplacian
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(VP, FP, L);


    // Debug output
    std::cout << "Barycentric constraint matrix: " << Ac.rows() << "x" << Ac.cols() << std::endl;
    std::cout << "Laplacian matrix: " << L.rows() << "x" << L.cols() << std::endl;

    // Step 4: Build and solve linear system
    Eigen::SparseMatrix<double> Acl_w;
    sparseVerticalConcat(sqrt(1.0/w) * Ac, L, Acl_w);

    Eigen::MatrixXd bcl_w(bc.rows() + L.rows(), 3);
    bcl_w << sqrt(1.0/w) * bc,
             Eigen::MatrixXd::Zero(L.rows(), 3);

    stop = clock();
    time = double (stop - start) / CLOCKS_PER_SEC;
    std::cout << "time for computing all barycentric coordinates, cotangent matrix and all required matrices : " <<  time << " seconds" << std::endl;

    // Solve the linear system
    sparseLSsolve(Acl_w, bcl_w, Vcl, true);
}

extern "C" DLL_EXPORT int computeSurfaceMesh(
    const Point3D* points, const Point3D* normals, int numPoints,
    double* outVertices, int* outFaces,
    int* outVertexCount, int* outFaceCount,
    int maxVertices, int maxFaces) {

    if (numPoints <= 0 || !points || !normals || !outVertices || !outFaces || !outVertexCount || !outFaceCount) {
        return 0;
    }

    // Convert to CGAL/Eigen
    std::vector<Point> cgal_points;
    std::vector<Vector> cgal_normals;
    std::unordered_set<Point> unique_vertices;

    Eigen::MatrixXd points_mat(3, numPoints);

    for (int i = 0; i < numPoints; ++i) {
        Point p(points[i].x, points[i].y, points[i].z);
        Vector n(normals[i].x, normals[i].y, normals[i].z);
        cgal_points.push_back(p);
        cgal_normals.push_back(n);
        points_mat.col(i) = Eigen::Vector3d(p.x(), p.y(), p.z());
    }

    // Delaunay triangulation
    Delaunay dt;
    std::vector<std::pair<Point, int>> points_with_info;
    for (int i = 0; i < numPoints; ++i) {
        points_with_info.emplace_back(cgal_points[i], i);
    }
    dt.insert(points_with_info.begin(), points_with_info.end());

    // Compute circumspheres
    std::vector<std::pair<Point, double>> valid_spheres;

    for (auto it = dt.finite_cells_begin(); it != dt.finite_cells_end(); ++it) {
        Tetrahedron tetra(it->vertex(0)->point(), it->vertex(1)->point(),
                          it->vertex(2)->point(), it->vertex(3)->point());
        Point c = CGAL::circumcenter(tetra);
        double r = std::sqrt(CGAL::squared_distance(c, tetra.vertex(0)));

        if (dt.is_infinite(dt.locate(c))) continue;

        bool valid = true;
        for (int i = 0; i < 4; ++i) {
            int idx = it->vertex(i)->info();
            Vector to_v = it->vertex(i)->point() - c;
            if (compute_angle(to_v, cgal_normals[idx]) >= 1.0) {
                valid = false;
                break;
            }
        }

        if (valid) 
        {   valid_spheres.emplace_back(c, r);
            for (int i = 0; i < 4; ++i) {
                unique_vertices.insert(it->vertex(i)->point());
        }
        }
    }

    // Voxel grid setup
    int res = 128;
    Eigen::Vector3d min_corner = points_mat.rowwise().minCoeff();
    Eigen::Vector3d max_corner = points_mat.rowwise().maxCoeff();
    Eigen::Vector3d diag = max_corner - min_corner;
    double mesh_size = diag.norm();
    double voxel_size = mesh_size / res;
    double padding = 0.1 * mesh_size;

    min_corner -= Eigen::Vector3d::Constant(padding);
    max_corner += Eigen::Vector3d::Constant(padding);
    double step = mesh_size / static_cast<double>(res);

    int N = res * res * res;
    Eigen::VectorXd S = Eigen::VectorXd::Constant(N, std::numeric_limits<double>::max());
    Eigen::MatrixXd GV(N, 3);

    for (int i = 0; i < res; ++i) {
        for (int j = 0; j < res; ++j) {
            for (int k = 0; k < res; ++k) {
                int index = i * res * res + j * res + k;
                GV.row(index) = min_corner + Eigen::Vector3d(i, j, k) * step;
            }
        }
    }

    // Fill SDF
    for (const auto& [center_pt, radius] : valid_spheres) {
        Eigen::Vector3d center(center_pt.x(), center_pt.y(), center_pt.z());

        Eigen::Vector3d min_bb = (center - Eigen::Vector3d::Constant(radius) - min_corner) / step;
        Eigen::Vector3d max_bb = (center + Eigen::Vector3d::Constant(radius) - min_corner) / step;

        int i_min = std::max(0, (int)std::floor(min_bb.x()));
        int j_min = std::max(0, (int)std::floor(min_bb.y()));
        int k_min = std::max(0, (int)std::floor(min_bb.z()));

        int i_max = std::min(res - 1, (int)std::ceil(max_bb.x()));
        int j_max = std::min(res - 1, (int)std::ceil(max_bb.y()));
        int k_max = std::min(res - 1, (int)std::ceil(max_bb.z()));

        for (int i = i_min; i <= i_max; ++i) {
            for (int j = j_min; j <= j_max; ++j) {
                for (int k = k_min; k <= k_max; ++k) {
                    int index = i * res * res + j * res + k;
                    Eigen::Vector3d p = min_corner + Eigen::Vector3d(i, j, k) * step;
                    double d = (p - center).norm() - radius;
                    S(index) = std::min(S(index), d);
                }
            }
        }
    }

    // Marching cubes
    Eigen::MatrixXd mcV;
    Eigen::MatrixXi mcF;
    igl::copyleft::marching_cubes(S, GV, res, res, res, mcV, mcF);

    Eigen::MatrixXd filte_v = unorderedSetToEigenMatrix(unique_vertices);

    uniform_remesh(mcV, mcF, voxel_size);


    computeFromOctahedronLapAndBary(1.0, filte_v, mcV, mcF, mcV);

    int vCount = mcV.rows();
    int fCount = mcF.rows();

    *outVertexCount = std::min(vCount, maxVertices);
    *outFaceCount = std::min(fCount, maxFaces);

    for (int i = 0; i < *outVertexCount; ++i) {
        outVertices[i * 3 + 0] = mcV(i, 0);
        outVertices[i * 3 + 1] = mcV(i, 1);
        outVertices[i * 3 + 2] = mcV(i, 2);
    }

    for (int i = 0; i < *outFaceCount; ++i) {
        outFaces[i * 3 + 0] = mcF(i, 0);
        outFaces[i * 3 + 1] = mcF(i, 1);
        outFaces[i * 3 + 2] = mcF(i, 2);
    }

     // Invert face orientation by swapping two indices in each face
    for (int i = 0; i < mcF.rows(); ++i) {
        std::swap(mcF(i, 0), mcF(i, 1));  // Flip orientation
    }



    return 1;
}
