#ifndef CIRCUMSPHERE_DLL_H
#define CIRCUMSPHERE_DLL_H

#ifdef _WIN32
  #define DLL_EXPORT __declspec(dllexport)
#else
  #define DLL_EXPORT
#endif

// Struct representing a 3D point or normal
struct Point3D {
    double x, y, z;
};

extern "C" {

// Computes a surface mesh from input points and normals via circumspheres + marching cubes.
// - points: input point cloud (array of Point3D)
// - normals: corresponding normals (same length as points)
// - numPoints: number of input points
// - outVertices: output array of doubles (size maxVertices * 3)
// - outFaces: output array of ints (size maxFaces * 3)
// - outVertexCount: pointer to number of vertices written
// - outFaceCount: pointer to number of faces written
// - maxVertices: maximum number of vertices to write
// - maxFaces: maximum number of faces to write
//
// Returns 1 on success, 0 on failure.
DLL_EXPORT int computeSurfaceMesh(
    const Point3D* points,
    const Point3D* normals,
    int numPoints,
    double* outVertices,
    int* outFaces,
    int* outVertexCount,
    int* outFaceCount,
    int maxVertices,
    int maxFaces
);

}

#endif // CIRCUMSPHERE_DLL_H
