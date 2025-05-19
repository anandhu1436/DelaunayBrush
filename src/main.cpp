#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include "CircumsphereDLL.h"

// Wrapper function to call the DLL
bool computeSurfaceMeshWrapper(
    const std::vector<Point3D>& input_points,
    const std::vector<Point3D>& input_normals,
    int res,
    Eigen::MatrixXd& outV,
    Eigen::MatrixXi& outF) 
{
    const int maxVertices = 90000;
    const int maxFaces = 180000;

    std::vector<double> outVertices(maxVertices * 3);
    std::vector<int> outFaces(maxFaces * 3);
    int outVertexCount = 0;
    int outFaceCount = 0;

    int success = computeSurfaceMesh(
        input_points.data(), input_normals.data(), static_cast<int>(input_points.size()),
        outVertices.data(), outFaces.data(),
        &outVertexCount, &outFaceCount,
        maxVertices, maxFaces
    );

    if (!success) {
        std::cerr << "computeSurfaceMesh failed\n";
        return false;
    }

    outV.resize(outVertexCount, 3);
    for (int i = 0; i < outVertexCount; ++i) {
        outV(i, 0) = outVertices[i * 3 + 0];
        outV(i, 1) = outVertices[i * 3 + 1];
        outV(i, 2) = outVertices[i * 3 + 2];
    }

    outF.resize(outFaceCount, 3);
    for (int i = 0; i < outFaceCount; ++i) {
        outF(i, 0) = outFaces[i * 3 + 0];
        outF(i, 1) = outFaces[i * 3 + 1];
        outF(i, 2) = outFaces[i * 3 + 2];
    }

    return true;
}


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
            temp_normals.emplace_back(nx, ny, nz);
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

// Main function
int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " input.obj resolution [remesh_factor]\n";
        return 1;
    }

    std::string obj_file = argv[1];
    int res = std::stoi(argv[2]);

    Eigen::MatrixXd points, normals;
    if (!read_obj(obj_file, points, normals)) {
        std::cerr << "Failed to read OBJ file: " << obj_file << std::endl;
        return 1;
    }

    std::cout << "Loaded " << points.rows() << " vertices\n";

    // Convert to Point3D format
    std::vector<Point3D> input_points(points.rows());
    std::vector<Point3D> input_normals(points.rows());

    for (int i = 0; i < points.rows(); ++i) {
        input_points[i] = { points(i, 0), points(i, 1), points(i, 2) };
        if (normals.rows() == points.rows()) {
            input_normals[i] = { -1*normals(i, 0), -1*normals(i, 1), -1*normals(i, 2) };
        } else {
            input_normals[i] = { 0, 0, 1 }; // default normal if not available
        }
    }

    // Compute surface mesh
    Eigen::MatrixXd mcV;
    Eigen::MatrixXi mcF;
    std::cout << "Computing surface mesh...\n";

    clock_t start = clock();
    bool success = computeSurfaceMeshWrapper(input_points, input_normals, res, mcV, mcF);
    double elapsed = double(clock() - start) / CLOCKS_PER_SEC;

    if (!success) {
        return 1;
    }

    std::cout << "Mesh generated with " << mcV.rows() << " vertices and " << mcF.rows() << " faces in " << elapsed << " sec\n";

    // Save result
    std::string out_path = obj_file.substr(0, obj_file.find_last_of(".")) + "_surface_mesh.obj";
    igl::writeOBJ(out_path, mcV, mcF);
    std::cout << "Saved mesh to: " << out_path << std::endl;

    return 0;
}
