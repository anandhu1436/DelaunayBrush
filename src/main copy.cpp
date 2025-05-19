


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <igl/harmonic.h>

#include <sstream>
#include <vector>
#include <string>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/IO/OBJ.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include "circumspheredll.h"

#include <igl/copyleft/marching_cubes.h>
#include <Eigen/Dense>
#include <vector>
#include <limits>

#include <igl/writeOBJ.h>   
#include <igl/readPLY.h>
#include <igl/loop.h>

#include <set>
#include <stdexcept>

#include <igl/readOBJ.h>
#include <igl/biharmonic_coordinates.h>
#include <igl/point_mesh_squared_distance.h>
#include <chrono>

#include <igl/barycentric_coordinates.h>
#include <igl/cotmatrix.h>
#include <Eigen/Sparse>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

namespace PMP = CGAL::Polygon_mesh_processing;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> SurfaceMesh;

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


void generate_sphere_mesh(
    double radius, int resolution,
    Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    std::vector<Eigen::RowVector3d> vertices;
    std::vector<Eigen::RowVector3i> faces;

    // Create vertices
    for (int i = 0; i <= resolution; ++i)
    {
        double theta = M_PI * i / resolution; // Angle from top to bottom
        for (int j = 0; j <= resolution; ++j)
        {
            double phi = 2 * M_PI * j / resolution; // Angle around the sphere
            double x = radius * sin(theta) * cos(phi);
            double y = radius * sin(theta) * sin(phi);
            double z = radius * cos(theta);
            vertices.emplace_back(x, y, z);
        }
    }

    // Create faces
    for (int i = 0; i < resolution; ++i)
    {
        for (int j = 0; j < resolution; ++j)
        {
            int v0 = i * (resolution + 1) + j;
            int v1 = v0 + 1;
            int v2 = v0 + resolution + 1;
            int v3 = v2 + 1;

            faces.emplace_back(v2, v1, v0);
            faces.emplace_back(v2, v3, v1);
        }
    }

    // Convert to Eigen matrices
    V.resize(vertices.size(), 3);
    F.resize(faces.size(), 3);

    for (size_t i = 0; i < vertices.size(); ++i)
        V.row(i) = vertices[i];

    for (size_t i = 0; i < faces.size(); ++i)
        F.row(i) = faces[i];
}


void save_combined_spheres(
    const std::vector<std::pair<Eigen::RowVector3d, double>>& circumspheres, // List of circumspheres (center, radius)
    const std::string& output_file)
{
    int resolution = 20;  // Sphere resolution
    Eigen::MatrixXd V_base;
    Eigen::MatrixXi F_base;
    generate_sphere_mesh(1.0, resolution, V_base, F_base); // Generate unit sphere

    std::vector<Eigen::MatrixXd> vertices_list;
    std::vector<Eigen::MatrixXi> faces_list;
    int vertex_offset = 0;

    // Iterate through all circumspheres
    for (size_t i = 0; i < circumspheres.size(); ++i) {
        const Eigen::RowVector3d& center = circumspheres[i].first;
        double radius = circumspheres[i].second;

        // Scale and translate the base sphere
        Eigen::MatrixXd V_scaled = V_base * radius; 
        Eigen::MatrixXd V_translated = V_scaled.rowwise() + center;

        // Adjust face indices
        Eigen::MatrixXi F_translated = F_base.array() + vertex_offset;


        // Store the results
        vertices_list.push_back(V_translated);
        faces_list.push_back(F_translated);
        vertex_offset += V_translated.rows();
    }


    // Combine all vertices and faces into final matrices
    int total_vertices = 0;
    int total_faces = 0;

    for (size_t i = 0; i < vertices_list.size(); ++i) {
        total_vertices += vertices_list[i].rows();
        total_faces += faces_list[i].rows();
    }

    Eigen::MatrixXd V_final(total_vertices, 3);
    Eigen::MatrixXi F_final(total_faces, 3);

    int v_offset = 0, f_offset = 0;
    for (size_t i = 0; i < vertices_list.size(); ++i) {
        V_final.block(v_offset, 0, vertices_list[i].rows(), 3) = vertices_list[i];
        F_final.block(f_offset, 0, faces_list[i].rows(), 3) = faces_list[i];
        v_offset += vertices_list[i].rows();
        f_offset += faces_list[i].rows();
    }

    // Save the final mesh to an OBJ file
    igl::writeOBJ(output_file, V_final, F_final);
}

void saveSpheresToFile(const std::vector<std::pair<Eigen::RowVector3d, double>>& spheres, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file for writing: " << filename << std::endl;
        return;
    }

    for (const auto& [center, radius] : spheres) {
        outFile << center.x() << " " << center.y() << " " << center.z() << " " << radius << "\n";
    }

    outFile.close();
}

// typedef CGAL::Nef_polyhedron_3<K> Nef_polyhedron;

typedef std::tuple<int, int, int> Voxel;
struct VoxelHash {
    std::size_t operator()(const Voxel& mcV) const {
        auto [x, y, z] = mcV;
        return std::hash<int>()(x) ^ (std::hash<int>()(y) << 1) ^ (std::hash<int>()(z) << 2);
    }
};





// Define CGAL types
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<int, K> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Tetrahedron_3 Tetrahedron;
typedef CGAL::Surface_mesh<Point> Mesh;
struct PointHash {
    std::size_t operator()(const Point& p) const {
        std::size_t h1 = std::hash<double>()(CGAL::to_double(p.x()));
        std::size_t h2 = std::hash<double>()(CGAL::to_double(p.y()));
        std::size_t h3 = std::hash<double>()(CGAL::to_double(p.z()));
        return h1 ^ (h2 << 1) ^ (h3 << 2);  // Combine the hash values for x, y, and z
    }
};



#include <igl/dual_contouring.h>

 void dualContourSpheres( const std::vector<std::pair<Eigen::RowVector3d,double>>& spheres,
                    Eigen::MatrixXd&  V,          // ← vertices
                    Eigen::MatrixXi&  F,          // ← quads (or triangles)
                    bool              triangles   = false, // false→quads
                    int               grid_res    = 128)    // nx=ny=nz
{
    using Row3d = Eigen::RowVector3d;

    // ------------------------------------------------------------
    // 1.  Signed‑distance of union of spheres
    // ------------------------------------------------------------
    auto f = [&](const Row3d& x)->double
    {
        double d = std::numeric_limits<double>::infinity();
        for(const auto& s : spheres)
        d = std::min(d, (x-s.first).norm() - s.second);
        return d;
    };

    // Finite‑difference gradient (normalised)
    auto f_grad = [&](const Row3d& x)->Row3d
    {
        const double eps = 1e-6;
        Row3d g;
        for(int c=0;c<3;++c)
        {
        Row3d e = Row3d::Zero(); e(c)=eps;
        g(c) = (f(x+e) - f(x-e)) / (2*eps);
        }
        return g.normalized();
    };

    // ------------------------------------------------------------
    // 2.  Bounding box (+ 5 % padding)
    // ------------------------------------------------------------
    Row3d bb_min( 1e30, 1e30, 1e30),
            bb_max(-1e30,-1e30,-1e30);
    for(const auto& s : spheres){
        // Option A: work in Array domain
bb_min = bb_min.array().min(s.first.array() - s.second).matrix();
bb_max = bb_max.array().max(s.first.array() + s.second).matrix();

    }
    double pad = 0.05 * (bb_max - bb_min).norm();
    bb_min.array() -= pad;  bb_max.array() += pad;

    // ------------------------------------------------------------
    // 3.  Run dual contouring
    // ------------------------------------------------------------
    bool constrained = false;   // unconstrained QEF
    bool root_find   = true;    // Newton refinement
    igl::dual_contouring(
        f, f_grad,
        bb_min, bb_max,
        grid_res, grid_res, grid_res,
        constrained,
        triangles,
        root_find,
        V, F);
}


// Convert Point3D to Eigen::Vector3d
Eigen::Vector3d toEigen(const Point3D &p) {
    return Eigen::Vector3d(p.x, p.y, p.z);
}

// Define the implicit function for the union of spheres
double union_of_spheres(const Eigen::Vector3d &x, const std::vector<Point3D> &centers, const std::vector<double> &radii) {
    double min_val = std::numeric_limits<double>::max();
    for (size_t i = 0; i < centers.size(); ++i) {
        double val = (x - toEigen(centers[i])).squaredNorm() - radii[i] * radii[i];
        min_val = std::min(min_val, val);
    }
    return min_val;
}



// Function to voxelize the spheres
void voxelize_spheres(const std::vector<std::pair<Point, double>>& spheres, double voxel_size,
    std::unordered_set<Voxel, VoxelHash>& voxel_grid) {
    for (const auto& [center, radius] : spheres) {
    int min_x = std::floor((center.x() - radius) / voxel_size);
    int max_x = std::ceil((center.x() + radius) / voxel_size);
    int min_y = std::floor((center.y() - radius) / voxel_size);
    int max_y = std::ceil((center.y() + radius) / voxel_size);
    int min_z = std::floor((center.z() - radius) / voxel_size);
    int max_z = std::ceil((center.z() + radius) / voxel_size);

    for (int x = min_x; x <= max_x; ++x) {
    for (int y = min_y; y <= max_y; ++y) {
    for (int z = min_z; z <= max_z; ++z) {
    Point voxel_center(x * voxel_size, y * voxel_size, z * voxel_size);
    if (CGAL::squared_distance(voxel_center, center) <= radius * radius) {
        voxel_grid.insert({x, y, z});
    }
    }
    }
    }
    }
    }

    // Function to extract boundary faces and save as OBJ
    void extract_and_save_boundary_faces(const std::unordered_set<Voxel, VoxelHash>& voxel_grid,
        double voxel_size, const std::string& filename) {
    std::vector<std::array<int, 3>> offsets = {
    {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}
    };

    std::vector<std::array<Point, 3>> faces;
    std::unordered_map<Point, int, PointHash> point_index;
    int index = 1;

    std::ofstream out(filename);
    if (!out) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
    }

    for (const auto& voxel : voxel_grid) {
    auto [x, y, z] = voxel;
    for (size_t i = 0; i < offsets.size(); ++i) {
    auto& offset = offsets[i];
    Voxel neighbor = {x + offset[0], y + offset[1], z + offset[2]};
    if (voxel_grid.find(neighbor) == voxel_grid.end()) {
    Point base(x * voxel_size, y * voxel_size, z * voxel_size);
    double half = voxel_size / 2.0;

    std::array<Point, 8> corners = {
    Point(base.x() - half, base.y() - half, base.z() - half),
    Point(base.x() + half, base.y() - half, base.z() - half),
    Point(base.x() + half, base.y() + half, base.z() - half),
    Point(base.x() - half, base.y() + half, base.z() - half),
    Point(base.x() - half, base.y() - half, base.z() + half),
    Point(base.x() + half, base.y() - half, base.z() + half),
    Point(base.x() + half, base.y() + half, base.z() + half),
    Point(base.x() - half, base.y() + half, base.z() + half)
    };

    std::array<std::array<int, 3>, 2> faces_indices;
    switch (i) {
    case 0: faces_indices = {{{0, 4, 7}, {7, 3, 0}}}; break;
    case 1: faces_indices = {{{1, 2, 6}, {6, 5, 1}}}; break;
    case 2: faces_indices = {{{0, 1, 5}, {5, 4, 0}}}; break;
    case 3: faces_indices = {{{2, 3, 7}, {7, 6, 2}}}; break;
    case 4: faces_indices = {{{0, 3, 2}, {2, 1, 0}}}; break;
    case 5: faces_indices = {{{4, 5, 6}, {6, 7, 4}}}; break;
    }

    for (const auto& indices : faces_indices) {
    faces.push_back({corners[indices[0]], corners[indices[1]], corners[indices[2]]});
    }
    }
    }
    }

    for (const auto& face : faces) {
    for (const auto& vertex : face) {
    if (point_index.find(vertex) == point_index.end()) {
    point_index[vertex] = index++;
    out << "mcV " << CGAL::to_double(vertex.x()) << " "
    << CGAL::to_double(vertex.y()) << " "
    << CGAL::to_double(vertex.z()) << "\n";
    }
    }
    }

    for (const auto& face : faces) {
    out << "mcF " << point_index[face[0]] << " "
    << point_index[face[1]] << " "
    << point_index[face[2]] << "\n";
    }

    out.close();
    std::cout << "Saved boundary voxelized mesh to " << filename << "\n";
}
// Function to generate a single sphere mesh

// Function to read vertices and normals from an OBJ file
bool read_obj(const std::string& filename, Eigen::MatrixXd& points, Eigen::MatrixXd& normals) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return false;
    }

    std::vector<Eigen::Vector3d> temp_points;
    std::vector<Eigen::Vector3d> temp_normals;

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string type;
        iss >> type;

        if (type == "v") {
            double x, y, z;
            iss >> x >> y >> z;
            temp_points.emplace_back(x, y, z);
        } else if (type == "vn") {
            double nx, ny, nz;
            iss >> nx >> ny >> nz;
            // temp_normals.emplace_back(nx, ny, nz);
            temp_normals.emplace_back(-1*nx, -1*ny, -1*nz);
        }
    }

    file.close();

    if (temp_points.size() != temp_normals.size()) {
        std::cerr << "Mismatch between the number of vertices and normals in the OBJ file." << std::endl;
        return false;
    }

    // Resize matrices to 3 x N
    size_t N = temp_points.size();
    points.resize(N, 3);
    normals.resize(N, 3);

    for (size_t i = 0; i < N; ++i) {
        points.row(i) = temp_points[i];
        normals.row(i) = temp_normals[i];
    }

    return true;
}
double distance(const Point3D& p1, const Point3D& p2) {
    return std::sqrt((p1.x - p2.x) * (p1.x - p2.x) +
                     (p1.y - p2.y) * (p1.y - p2.y) +
                     (p1.z - p2.z) * (p1.z - p2.z));
}

double findShortestDistance(const std::vector<Point3D>& points) {
    double min_dist = std::numeric_limits<double>::max();

    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = i + 1; j < points.size(); ++j) {
            double d = distance(points[i], points[j]);
            if (d < min_dist) {
                min_dist = d;
            }
        }
    }

    return min_dist;
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
            X = solver.solve(b);
            // X = solver.solveWithGuess( b, previousSolution);
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
    igl::point_mesh_squared_distance(inputPts, VP, FP, sqrD, I, C);

    const int n = inputPts.rows();
    const int n_vertices = VP.rows();

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

    // Solve the linear system
    sparseLSsolve(Acl_w, bcl_w, Vcl, true);
}

int main(int argc, char** argv) {
    double remsh_factor = 1;
    if (argc >= 4) {
        remsh_factor = std::stod(argv[3]);
    }

    Eigen::MatrixXd points;
    Eigen::MatrixXd normals;

    std::string obj_file = argv[1];
    int res=std::stoi(argv[2]);


    Eigen::MatrixXd mcV,dcV;
    Eigen::MatrixXi mcF,dcF;

    // Save results
    std::string output_file = obj_file.substr(0, obj_file.find_last_of(".")) + "_marching.obj";
    std::string dual_contour = obj_file.substr(0, obj_file.find_last_of(".")) + "_dual_contouring.obj";
    std::string sphere_file = obj_file.substr(0, obj_file.find_last_of(".")) + "_circumspheres.txt";
    std::string sphere_file_obj = obj_file.substr(0, obj_file.find_last_of(".")) + "_circumspheres.obj";
    std::string filtered = obj_file.substr(0, obj_file.find_last_of(".")) + "_filtered.ply";
    std::string final_mc = obj_file.substr(0, obj_file.find_last_of(".")) + "_final_mc.obj";
    std::string final_dc = obj_file.substr(0, obj_file.find_last_of(".")) + "_final_dc.obj";
    std::string timing = obj_file.substr(0, obj_file.find_last_of(".")) + "timing.txt";
    std::string uniorm_mesh = obj_file.substr(0, obj_file.find_last_of(".")) + "_uniform.obj";

    std::ofstream outfile(timing, std::ios::app); // open in append mode

    double  voxel_size ;
    double total_time = 0.0;
    {
        if (!read_obj(obj_file, points, normals)) {
            std::cerr << "Failed to read OBJ file: " << obj_file << std::endl;
            return 1;
        }

        std::cout << "Read " << points.size() << " vertices and normals " << normals.size() << std::endl;

        // Convert to Point3D format for DLL
        std::vector<Point3D> input_points(points.rows());
        std::vector<Point3D> input_normals(normals.rows());

        for (size_t i = 0; i < points.rows(); ++i) {
            input_points[i] = {CGAL::to_double(points(i, 0)), CGAL::to_double(points(i, 1)), CGAL::to_double(points(i, 2))};      
            input_normals[i] = {CGAL::to_double(normals(i, 0)), CGAL::to_double(normals(i, 1)), CGAL::to_double(normals(i, 2))};
        }

        // Allocate output arrays
        int maxOutputSize = 90000;  // Adjust as needed
        std::vector<Point3D> outCenters(maxOutputSize);
        std::vector<double> outRadii(maxOutputSize);

        clock_t start = clock();
        // Extract mesh using marching cubes
        
        int num_spheres = computeCircumspheres(
            input_points.data(), input_normals.data(),
            static_cast<int>(input_points.size()),
            outCenters.data(), outRadii.data(), maxOutputSize,filtered);

        clock_t stop = clock();
        double time = double (stop - start) / CLOCKS_PER_SEC;
        total_time += time;
            // Call the DLL function

        std::cout << "Computed " << num_spheres << " valid circumspheres\n";
        outfile << "Proxy time: " << time << " seconds" << std::endl;
        // outfile.close();
        // Store results in a vector
        std::vector<std::pair<Eigen::RowVector3d, double>> spheres;
        for (int i = 0; i < num_spheres; ++i) {
            spheres.emplace_back(Eigen::RowVector3d(outCenters[i].x, outCenters[i].y, outCenters[i].z), outRadii[i]);
        }

        // std::unordered_set<Voxel, VoxelHash> voxel_grid;

        // // Timer start for voxelization
        // voxel_size = 2*findShortestDistance(input_points);
        start = clock();
        saveSpheresToFile(spheres, sphere_file);
        save_combined_spheres(spheres, sphere_file_obj);
        stop = clock();
        time = double (stop - start) / CLOCKS_PER_SEC;
        std::cout << "time for saving spheres : " <<  time << " seconds" << std::endl;

        // std::cout << "Saved circumspheres to " << sphere_file << std::endl;

        // std::cout << "Voxel size: " << voxel_size << std::endl;
        // auto t1 = std::chrono::high_resolution_clock::now();
        // voxelize_spheres(spheres, voxel_size, voxel_grid);
        // auto t2 = std::chrono::high_resolution_clock::now();
        // std::cout << "Voxelization time: " 
        //         << std::chrono::duration<double>(t2 - t1).count() 
        //         << " seconds" << std::endl;

        // auto t3 = std::chrono::high_resolution_clock::now();
        // extract_and_save_boundary_faces(voxel_grid, voxel_size, voxel_mesh);
        // auto t4 = std::chrono::high_resolution_clock::now();
        // std::cout << "bpondary time: " 
        //         << std::chrono::duration<double>(t4 - t3).count() 
        //         << " seconds" << std::endl;

        // return 0;

        start = clock();
        // float voxel_size = findShortestDistance(input_points);
        // Compute bounding box
        Eigen::Vector3d min_corner = points.rowwise().minCoeff();
        Eigen::Vector3d max_corner = points.rowwise().maxCoeff();

        // Compute diagonal length of the bounding box (optional, for scale)
        Eigen::Vector3d diag = max_corner - min_corner;
        double mesh_size = 10*diag.norm();  // or diag.norm(), depending on your goal

        voxel_size = diag.norm() / static_cast<double>(res);

        // Add padding (optional)
        double padding =  0.3*mesh_size;  // 10% padding  
        min_corner -= Eigen::Vector3d::Constant(padding);
        max_corner += Eigen::Vector3d::Constant(padding);

        // Recompute step and resolution info
        double step = mesh_size / static_cast<double>(res);
        int N = res * res * res;

        // Allocate
        Eigen::VectorXd S = Eigen::VectorXd::Constant(N, std::numeric_limits<double>::max());
        Eigen::MatrixXd GV(N, 3);

        // Fill voxel grid positions
        for (int i = 0; i < res; ++i) {
            for (int j = 0; j < res; ++j) {
                for (int k = 0; k < res; ++k) {
                    int index = i * res * res + j * res + k;
                    Eigen::Vector3d p = min_corner + Eigen::Vector3d(i, j, k) * step;
                    GV.row(index) = p;
                }
            }
        }


        // Per-sphere iteration
        for (int s = 0; s < outCenters.size(); ++s) {
            Eigen::Vector3d center(outCenters[s].x, outCenters[s].y, outCenters[s].z);

            double radius = outRadii[s];

            Eigen::Vector3d min_bb = (center - Eigen::Vector3d::Constant(radius) - min_corner) / step;
            Eigen::Vector3d max_bb = (center + Eigen::Vector3d::Constant(radius) - min_corner) / step;

            int i_min = std::max(0, static_cast<int>(std::floor(min_bb.x())));
            int j_min = std::max(0, static_cast<int>(std::floor(min_bb.y())));
            int k_min = std::max(0, static_cast<int>(std::floor(min_bb.z())));

            int i_max = std::min(res - 1, static_cast<int>(std::ceil(max_bb.x())));
            int j_max = std::min(res - 1, static_cast<int>(std::ceil(max_bb.y())));
            int k_max = std::min(res - 1, static_cast<int>(std::ceil(max_bb.z())));

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

        stop = clock();
        time = double (stop - start) / CLOCKS_PER_SEC;
        total_time += time;
        std::cout << "grid time : " << time << " seconds" << std::endl;
        outfile << "Grid time: " << time << " seconds" << std::endl;
        start = clock();
        // Extract mesh using marching cubes
        igl::copyleft::marching_cubes(S, GV, res, res, res, mcV, mcF);
        
        stop = clock();
        time = double (stop - start) / CLOCKS_PER_SEC;
        total_time += time;
        outfile << "Marching time: " << time << " seconds" << std::endl;
        std::cout << "marching time : " << time << " seconds" << std::endl;

        // start = clock();
        // // Extract mesh using marching cubes
        // dualContourSpheres(spheres, dcV, dcF, true, res);
        // end = clock();
        // time = double (end - start) / CLOCKS_PER_SEC;
        // std::cout << "dual contouring time : " << time << " seconds" << std::endl;
        // igl::writeOBJ(dual_contour, dcV, dcF);

        // Save the result
        igl::writeOBJ(output_file, mcV, mcF);

        std::cout << "Saved circumspheres to " << output_file << "\n";

}

    Eigen::MatrixXd filte_v;
    Eigen::MatrixXi filter_f;
    if(res<10)
    {   std::cout << output_file << std::endl;
        igl::readOBJ(output_file, mcV, mcF);
        igl::readOBJ(filtered, filte_v, filter_f);
        // final = obj_file.substr(0, obj_file.find_last_of(".")) + "_final2.obj";
    }
    else
    {
        igl::readPLY(filtered, filte_v, filter_f);
    }
    std::cout << "Read circumspheres: " << mcV.rows() << " vertices and " << mcF.rows() << " faces" << std::endl;
    // igl::loop(mcV,mcF,mcV,mcF);
    clock_t st = clock();
    uniform_remesh(mcV, mcF, remsh_factor*voxel_size); // adjust edge length to your threshold
    
    clock_t en = clock();
    double utime = double (en - st) / CLOCKS_PER_SEC;
    total_time += utime;
    outfile << "Uniform mesh time: " << utime << " seconds" << std::endl;
    if(!igl::writeOBJ(uniorm_mesh, mcV, mcF)) {
        std::cerr << "Failed to save deformed mesh!" << std::endl;
        return -1;
    }

    // Eigen::VectorXi I;
    // Eigen::MatrixXd C, sqrD;
    // std::cout << "Closest vertex: " << mcF.rows() << std::endl;
    // igl::point_mesh_squared_distance(filte_v, mcV, mcF, sqrD, I, C);

    // std::set<int> used_indices;
    // std::vector<int> unique_indices;
    // std::vector<Eigen::RowVector3d> unique_positions;

    // for (int i = 0; i < filte_v.rows(); ++i) {
    //     int face_idx = I(i);
    //     Eigen::RowVector3d point = filte_v.row(i);

    //     // Get vertex indices of the closest triangle
    //     int v0 = mcF(face_idx, 0);
    //     int v1 = mcF(face_idx, 1);
    //     int v2 = mcF(face_idx, 2);

    //     // Find closest vertex among the 3
    //     std::vector<int> verts = {v0, v1, v2};
    //     int closest_vid = -1;
    //     double min_dist = std::numeric_limits<double>::max();

    //     for (int vid : verts) {
    //         double dist = (point - mcV.row(vid)).squaredNorm();
    //         if (dist < min_dist) {
    //             min_dist = dist;
    //             closest_vid = vid;
                
    //         }

    //     }

    //     // Skip if already used
    //     if (used_indices.count(closest_vid)) continue;

    //     used_indices.insert(closest_vid);
    //     unique_indices.push_back(closest_vid);
    //     unique_positions.push_back(point);
    // }
    // std::cout << "Unique indices: " << unique_indices.size() << std::endl;

    // // Convert to Eigen matrices
    // Eigen::VectorXi b(unique_indices.size());
    // Eigen::MatrixXd bc(unique_positions.size(), 3);
    // std::cout << "Unique indices: " << unique_indices.size() << std::endl;
    // for (int i = 0; i < unique_indices.size(); ++i) {
    //     b(i) = unique_indices[i];
    //     bc.row(i) = mcV.row(unique_indices[i]);
    // }

    // Compute weights
    // Eigen::MatrixXd W;
    // igl::harmonic(mcV, mcF, b, bc, 2, W);


    // Compute new vertex positions
    // Eigen::MatrixXd V_new = W ;

    // Step 5: Save the deformed mesh
    Eigen::MatrixXd mcV_new,dcV_new;

    clock_t start = clock();
    computeFromOctahedronLapAndBary(1.0, filte_v, mcV, mcF, mcV_new);
    // computeFromOctahedronLapAndBary(1.0, filte_v, dcV, dcF, dcV_new);
    // if (!igl::writeOBJ(final_dc, dcV_new, dcF)) {
    //     std::cerr << "Failed to save deformed mesh!" << std::endl;
    //     return -1;
    // }
    // Extract mesh using marching cubes
    clock_t end = clock();
    double time = double (end - start) / CLOCKS_PER_SEC;
    total_time += time;
    outfile << "Smoothing time: " << time << " seconds" << std::endl;
    outfile << "Total time: " << total_time << " seconds" << std::endl;
    outfile.close();

    if (!igl::writeOBJ(final_mc, mcV_new, mcF)) {
        std::cerr << "Failed to save deformed mesh!" << std::endl;
        return -1;
    }


    std::cout << "Deformation complete! Saved as 'deformed_mesh.obj'"<<final_mc<< std::endl;

    return 0;
}