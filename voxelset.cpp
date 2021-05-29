#include <QFile>
#include <QDataStream>
#include "voxelset.h"
#include "Triangulation.h"

using namespace R3Graph;

const int NEIGHBOURS6[6][3] = {         // slice, x, y
    { 0, -1, 0 },
    { 0, 1, 0 },
    { 0, 0, -1 },
    { 0, 0, 1 },    // First 4 element MUST BE on the same slice!!!
    //-------------------------------------------------------------
    { -1, 0, 0 },
    { 1, 0, 0 }
};

const int NEIGHBOURS4_3D[4][3] = {      // slice, x, y
    { 0, -1, 0 },
    { 0, 1, 0 },
    { 0, 0, -1 },
    { 0, 0, 1 }
};

const int Voxel::FACE_DIRECTIONS[6][3] = {
    {  0, -1,  0 }, // FACE_LEFT    slice, x, y
    {  0,  1,  0 }, // FACE_RIGHT
    {  0,  0, -1 }, // FACE_FRONT
    {  0,  0,  1 }, // FACE_BACK
    { -1,  0,  0 }, // FACE_BOTTOM
    {  1,  0,  0 }  // FACE_TOP
};

void VoxelSet::initialize(
   int x_max, int y_max, int max_slices
) {
    xMax = x_max;;
    yMax = y_max;
    maxSlices = max_slices;
    bitmasks.resize(maxSlices);
    numVoxels = 0.;
    for (int i = 0; i < maxSlices; ++i) {
        bitmasks.at(i).resizeBitmask(x_max, y_max);
        bitmasks.at(i).clear();
    }
    detected3D = false;
}

void VoxelSet::clear() {
    for (int i = 0; i < maxSlices; ++i) {
        bitmasks.at(i).clear();
    }
    numVoxels = 0.;
}

void VoxelSet::computeVoxelBox() {
    numVoxels = 0.;
    for (int slice = 0; slice < maxSlices; ++slice) {
        const Bitmask& bitmask = bitmaskAt(slice);
        for (int y = 0; y < yMax; ++y) {
            for (int x = 0; x < xMax; ++x) {
                if (bitmask.pixelAt(x, y) != 0) {
                    if (numVoxels == 0) {
                        // Initialize a voxel box
                        voxelBox.origin = Voxel(slice, x, y);
                        voxelBox.width = 0;
                        voxelBox.depth = 0;
                        voxelBox.height = 0;
                    } else {
                        if (x < voxelBox.xMin()) {
                            voxelBox.width += voxelBox.xMin() - x;
                            voxelBox.setXMin(x);
                        } else if (x > voxelBox.xMax()) {
                            voxelBox.width += x - voxelBox.xMax();
                        }

                        if (y < voxelBox.yMin()) {
                            voxelBox.depth += voxelBox.yMin() - y;
                            voxelBox.setYMin(y);
                        } else if (y > voxelBox.yMax()) {
                            voxelBox.depth += y - voxelBox.yMax();
                        }

                        if (slice < voxelBox.sliceMin()) {
                            voxelBox.height += voxelBox.sliceMin() - slice;
                            voxelBox.setSliceMin(slice);
                        } else if (slice > voxelBox.sliceMax()) {
                            voxelBox.height += slice - voxelBox.sliceMax();
                        }
                    }
                    numVoxels += 1.;
                } // end if
            } // end for (x...
        } // end for (y...
    } // end for (slice...
}

bool VoxelSet::save(QString path) const {
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly))
        return false;
    QDataStream out(&f);

    out << threshold;
    out << xMax;
    out << yMax;
    out << maxSlices;
    if (out.status() == QDataStream::WriteFailed)
        return false;
    for (int s = 0; s < maxSlices; ++s) {
        const Bitmask& bitmask = bitmaskAt(s);
        if (!bitmask.write(out))
            return false;
    }
    return true;
}

bool VoxelSet::load(QString path) {
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly))
        return false;
    QDataStream in(&f);
    in >> threshold;
    in >> xMax;
    in >> yMax;
    in >> maxSlices;
    if (in.status() != QDataStream::Ok)
        return false;
    for (int s = 0; s < maxSlices; ++s) {
        Bitmask& bitmask = bitmaskAt(s);
        if (!bitmask.read(in))
            return false;
    }
    return true;
}

PackedVoxelSet::PackedVoxelSet(const VoxelSet& vs) {
    operator=(vs);
}

PackedVoxelSet& PackedVoxelSet::operator=(const VoxelSet& vs) {
    xMax = vs.xMax;
    yMax = vs.yMax;
    maxSlices = vs.maxSlices;
    numVoxels = vs.numVoxels;
    voxelBox = vs.voxelBox;
    threshold = vs.threshold;
    detected3D = vs.detected3D;

    packedBitmasks.resize(vs.bitmasks.size());
    assert(maxSlices == (int) vs.bitmasks.size());
    for (int slice = 0; slice < maxSlices; ++slice) {
        packedBitmasks.at(slice) = vs.bitmasks.at(slice);
    }

    return *this;
}

VoxelSet& PackedVoxelSet::unpack(VoxelSet& vs) const {
    vs.xMax = xMax;
    vs.yMax = yMax;
    vs.maxSlices = maxSlices;
    vs.numVoxels = numVoxels;
    vs.voxelBox = voxelBox;
    vs.threshold = threshold;
    vs.detected3D = detected3D;

    vs.bitmasks.resize(packedBitmasks.size());
    assert(maxSlices == (int) packedBitmasks.size());
    for (int slice = 0; slice < maxSlices; ++slice) {
        packedBitmasks.at(slice).unpack(vs.bitmasks.at(slice));
    }

    return vs;
}

bool VoxelSet::faceOpen(const Voxel& v, int face) const {
    assert(voxelAt(v) != 0);
    Voxel neighbour(
        v.slice + Voxel::FACE_DIRECTIONS[face][0],
        v.point.x + Voxel::FACE_DIRECTIONS[face][1],
        v.point.y + Voxel::FACE_DIRECTIONS[face][2]
    );
    return (
        !voxelInVolume(neighbour) ||
        voxelAt(neighbour) == 0
    );
}

void computeTriangulationOfVoxelSet(
    Triangulation& triangulation,
    const VoxelSet& voxelSet,
    const R3Point& origin,
    double dx, double dy, double dz
) {
    const VoxelBox& voxelBox = voxelSet.voxelBox;
    int sliceStart = voxelBox.origin.slice - 2;
    if (sliceStart < 0)
        sliceStart = 0;

    int sliceFinish = voxelBox.origin.slice + voxelBox.height + 1;
    if (sliceFinish >= voxelSet.maxSlices)
        sliceFinish = voxelSet.maxSlices - 1;

    int ixmin = voxelBox.origin.point.x - 2;
    if (ixmin < 0)
        ixmin = 0;

    int ixmax = voxelBox.origin.point.x + voxelBox.width + 1;
    if (ixmax >= voxelSet.xMax)
        ixmax = voxelSet.xMax - 1;

    int iymin = voxelBox.origin.point.y - 2;
    if (iymin < 0)
        iymin = 0;

    int iymax = voxelBox.origin.point.y + voxelBox.depth + 1;
    if (iymax >= voxelSet.yMax)
        iymax = voxelSet.yMax - 1;

    triangulation.clear();

    // Indices of extended voxels vertices in array
    // Each voxel produce 8 vertices == extended voxels
    // voxel (slice, x, y) -> 8 extended voxels:
    //       (2*slice + s0, 2*x + s1, 2*y + s2), where
    //                                           si = +-1
    std::map<Voxel, int> vertexIndices;

    for (int slice = sliceStart; slice <= sliceFinish; ++slice) {
        for (int iy = iymin; iy <= iymax; ++iy) {
            for (int ix = ixmin; ix <= ixmax; ++ix) {
                if (voxelSet.voxelAt(slice, ix, iy) == 0)
                    continue;

                // Enumeration of cube vertices and faces:
                //        7         6
                //       +---------+          z
                //      /   top   /|        ^
                //   4 / |       / |        |
                //    +---------+5 |        |
                //Left|         |  | Right  |    ^ y
                //    |  |  back|  |        |   /
                //    |  + - - -|- +        |  /
                //    |  3      | /2        | /
                //    |/  Front |/          |/
                //    +---------+           ------> x
                //   0  bottom   1

                R3Point cubeVertices[8];
                //... R3Vector cubeGradients[8];

                Voxel cube(slice, ix, iy);
                R3Point cubeCenter = voxel3DCoord(
                    cube, origin, dx, dy, dz
                );

                // Bottom
                Voxel neighborVoxel =
                    Voxel(cube.slice - 1, cube.point.x, cube.point.y);
                R3Point bottomCenter = voxel3DCoord(
                    neighborVoxel,
                    origin, dx, dy, dz
                );

                // Top
                neighborVoxel.slice += 2;
                R3Point topCenter = voxel3DCoord(
                    neighborVoxel,
                    origin, dx, dy, dz
                );

                // Left
                neighborVoxel =
                    Voxel(cube.slice, cube.point.x - 1, cube.point.y);
                R3Point leftCenter = voxel3DCoord(
                    neighborVoxel,
                    origin, dx, dy, dz
                );

                // Right
                neighborVoxel.point.x += 2;
                R3Point rightCenter = voxel3DCoord(
                    neighborVoxel,
                    origin, dx, dy, dz
                );

                // Front
                neighborVoxel =
                    Voxel(cube.slice, cube.point.x, cube.point.y - 1);
                R3Point frontCenter = voxel3DCoord(
                    neighborVoxel,
                    origin, dx, dy, dz
                );

                // Back
                neighborVoxel.point.y += 2;
                R3Point backCenter = voxel3DCoord(
                    neighborVoxel,
                    origin, dx, dy, dz
                );

                cubeVertices[0] = cubeCenter +
                    (frontCenter - cubeCenter) * 0.5 +
                    (leftCenter - cubeCenter) * 0.5 +
                    (bottomCenter - cubeCenter) * 0.5;
                cubeVertices[1] = cubeCenter +
                    (frontCenter - cubeCenter) * 0.5 +
                    (rightCenter - cubeCenter) * 0.5 +
                    (bottomCenter - cubeCenter) * 0.5;
                cubeVertices[2] = cubeCenter +
                    (backCenter - cubeCenter) * 0.5 +
                    (rightCenter - cubeCenter) * 0.5 +
                    (bottomCenter - cubeCenter) * 0.5;
                cubeVertices[3] = cubeCenter +
                    (backCenter - cubeCenter) * 0.5 +
                    (leftCenter - cubeCenter) * 0.5 +
                    (bottomCenter - cubeCenter) * 0.5;

                cubeVertices[4] = cubeCenter +
                    (frontCenter - cubeCenter) * 0.5 +
                    (leftCenter - cubeCenter) * 0.5 +
                    (topCenter - cubeCenter) * 0.5;
                cubeVertices[5] = cubeCenter +
                    (frontCenter - cubeCenter) * 0.5 +
                    (rightCenter - cubeCenter) * 0.5 +
                    (topCenter - cubeCenter) * 0.5;
                cubeVertices[6] = cubeCenter +
                    (backCenter - cubeCenter) * 0.5 +
                    (rightCenter - cubeCenter) * 0.5 +
                    (topCenter - cubeCenter) * 0.5;
                cubeVertices[7] = cubeCenter +
                    (backCenter - cubeCenter) * 0.5 +
                    (leftCenter - cubeCenter) * 0.5 +
                    (topCenter - cubeCenter) * 0.5;

                Voxel extendedVoxels[8];
                for (int iv = 0; iv < 8; ++iv) {
                    extendedVoxels[iv] = Voxel(
                        cube.slice*2, cube.point.x*2, cube.point.y*2
                    );
                }
                --(extendedVoxels[0].slice);
                --(extendedVoxels[1].slice);
                --(extendedVoxels[2].slice);
                --(extendedVoxels[3].slice);
                ++(extendedVoxels[4].slice);
                ++(extendedVoxels[5].slice);
                ++(extendedVoxels[6].slice);
                ++(extendedVoxels[7].slice);

                --(extendedVoxels[0].point.x);
                --(extendedVoxels[3].point.x);
                --(extendedVoxels[4].point.x);
                --(extendedVoxels[7].point.x);
                ++(extendedVoxels[1].point.x);
                ++(extendedVoxels[2].point.x);
                ++(extendedVoxels[5].point.x);
                ++(extendedVoxels[6].point.x);

                --(extendedVoxels[0].point.y);
                --(extendedVoxels[1].point.y);
                --(extendedVoxels[4].point.y);
                --(extendedVoxels[5].point.y);
                ++(extendedVoxels[2].point.y);
                ++(extendedVoxels[3].point.y);
                ++(extendedVoxels[6].point.y);
                ++(extendedVoxels[7].point.y);

                int indices[8]; // Indices of vertices in array
                for (int iv = 0; iv < 8; ++iv)
                    indices[iv] = (-1);

                // Front face
                if (voxelSet.faceOpen(cube, Voxel::FACE_FRONT)) {
                    // Add vertices
                    if (
                        vertexIndices.count(
                            extendedVoxels[0]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(
                            cubeVertices[0]
                        );
                        indices[0] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[0]] = indices[0];
                    } else {
                        // Point is already in the array
                        indices[0] = vertexIndices[extendedVoxels[0]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[1]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(cubeVertices[1]);
                        indices[1] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[1]] = indices[1];
                    } else {
                        // Point is already in the array
                        indices[1] = vertexIndices[extendedVoxels[1]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[4]
                        ) == 0
                    )
                    {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(cubeVertices[4]);
                        indices[4] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[4]] = indices[4];
                    } else {
                        // Point is already in the array
                        indices[4] = vertexIndices[extendedVoxels[4]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[5]
                        ) == 0
                    )
                    {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(cubeVertices[5]);
                        indices[5] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[5]] = indices[5];
                    } else {
                        // Point is already in the array
                        indices[5] = vertexIndices[extendedVoxels[5]];
                    }

                    // Add triangles for this face
                    assert(
                        indices[0] >= 0 &&
                        indices[1] >= 0 &&
                        indices[4] >= 0 &&
                        indices[5] >= 0
                    );

                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[0], indices[1], indices[5]
                        )
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[0], indices[5], indices[4]
                        )
                    );
                } // end if (... FRONT_FACE

                // Back face
                if (voxelSet.faceOpen(cube, Voxel::FACE_BACK)) {
                    // Add vertices
                    if (
                        vertexIndices.count(
                            extendedVoxels[2]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(
                            cubeVertices[2]
                        );
                        indices[2] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[2]] = indices[2];
                    } else {
                        // Point is already in the array
                        indices[2] = vertexIndices[extendedVoxels[2]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[3]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(
                            cubeVertices[3]
                        );
                        indices[3] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[3]] = indices[3];
                    } else {
                        // Point is already in the array
                        indices[3] = vertexIndices[extendedVoxels[3]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[6]
                        ) == 0
                    )
                    {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(
                            cubeVertices[6]
                        );
                        indices[6] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[6]] = indices[6];
                    } else {
                        // Point is already in the array
                        indices[6] = vertexIndices[extendedVoxels[6]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[7]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(
                            cubeVertices[7]
                        );
                        indices[7] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[7]] = indices[7];
                    } else {
                        // Point is already in the array
                        indices[7] = vertexIndices[extendedVoxels[7]];
                    }

                    // Add triangles for this face
                    assert(
                        indices[2] >= 0 &&
                        indices[3] >= 0 &&
                        indices[6] >= 0 &&
                        indices[7] >= 0
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[2], indices[3], indices[6]
                        )
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[6], indices[3], indices[7]
                        )
                    );
                } // end if (... BACK_FACE

                // Left face
                if (voxelSet.faceOpen(cube, Voxel::FACE_LEFT)) {
                    // Add vertices
                    if (
                        vertexIndices.count(
                            extendedVoxels[0]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(cubeVertices[0]);
                        indices[0] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[0]] = indices[0];
                    } else {
                        // Point is already in the array
                        indices[0] = vertexIndices[extendedVoxels[0]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[3]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(cubeVertices[3]);
                        indices[3] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[3]] = indices[3];
                    } else {
                        // Point is already in the array
                        indices[3] = vertexIndices[extendedVoxels[3]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[4]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(cubeVertices[4]);
                        indices[4] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[4]] = indices[4];
                    } else {
                        // Point is already in the array
                        indices[4] = vertexIndices[extendedVoxels[4]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[7]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(cubeVertices[7]);
                        indices[7] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[7]] = indices[7];
                    } else {
                        // Point is already in the array
                        indices[7] = vertexIndices[extendedVoxels[7]];
                    }

                    // Add triangles for this face
                    assert(
                        indices[0] >= 0 &&
                        indices[3] >= 0 &&
                        indices[4] >= 0 &&
                        indices[7] >= 0
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[0], indices[7], indices[3]
                        )
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[7], indices[0], indices[4]
                        )
                    );
                } // end if (... LEFT_FACE

                // Right face
                if (voxelSet.faceOpen(cube, Voxel::FACE_RIGHT)) {
                    // Add vertices
                    if (
                        vertexIndices.count(
                            extendedVoxels[1]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(cubeVertices[1]);
                        indices[1] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[1]] = indices[1];
                    } else {
                        // Point is already in the array
                        indices[1] = vertexIndices[extendedVoxels[1]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[2]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(cubeVertices[2]);
                        indices[2] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[2]] = indices[2];
                    } else {
                        // Point is already in the array
                        indices[2] = vertexIndices[extendedVoxels[2]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[5]
                        ) == 0
                    ) {
                        // Add vertex to triangulation
                        triangulation.vertices.push_back(cubeVertices[5]);
                        indices[5] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[5]] = indices[5];
                    } else {
                        // Point is already in the array
                        indices[5] = vertexIndices[extendedVoxels[5]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[6]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(cubeVertices[6]);
                        indices[6] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[6]] = indices[6];
                    } else {
                        // Point is already in the array
                        indices[6] = vertexIndices[extendedVoxels[6]];
                    }

                    // Add triangles for this face
                    assert(
                        indices[1] >= 0 &&
                        indices[2] >= 0 &&
                        indices[5] >= 0 &&
                        indices[6] >= 0
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[1], indices[2], indices[6]
                        )
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[6], indices[5], indices[1]
                        )
                    );
                } // end if (... RIGHT_FACE

                // Bottom face
                if (voxelSet.faceOpen(cube, Voxel::FACE_BOTTOM)) {
                    // Add vertices
                    if (
                        vertexIndices.count(
                            extendedVoxels[0]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(cubeVertices[0]);
                        indices[0] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[0]] = indices[0];
                    } else {
                        // Point is already in the array
                        indices[0] = vertexIndices[extendedVoxels[0]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[1]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(cubeVertices[1]);
                        indices[1] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[1]] = indices[1];
                    } else {
                        // Point is already in the array
                        indices[1] = vertexIndices[extendedVoxels[1]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[2]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(cubeVertices[2]);
                        indices[2] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[2]] = indices[2];
                    } else {
                        // Point is already in the array
                        indices[2] = vertexIndices[extendedVoxels[2]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[3]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(cubeVertices[3]);
                        indices[3] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[3]] = indices[3];
                    } else {
                        // Point is already in the array
                        indices[3] = vertexIndices[extendedVoxels[3]];
                    }

                    // Add triangles for this face
                    assert(
                        indices[0] >= 0 &&
                        indices[1] >= 0 &&
                        indices[2] >= 0 &&
                        indices[3] >= 0
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[0], indices[2], indices[1]
                        )
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[0], indices[3], indices[2]
                        )
                    );
                } // end if (... BOTTOM_FACE

                // Top face
                if (voxelSet.faceOpen(cube, Voxel::FACE_TOP)) {
                    // Add vertices
                    if (
                        vertexIndices.count(
                            extendedVoxels[4]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(
                            cubeVertices[4]
                        );
                        indices[4] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[4]] = indices[4];
                    } else {
                        // Point is already in the array
                        indices[4] = vertexIndices[extendedVoxels[4]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[5]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(
                            cubeVertices[5]
                        );
                        indices[5] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[5]] = indices[5];
                    } else {
                        // Point is already in the array
                        indices[5] = vertexIndices[extendedVoxels[5]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[6]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(
                            cubeVertices[6]
                        );
                        indices[6] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[6]] = indices[6];
                    } else {
                        // Point is already in the array
                        indices[6] = vertexIndices[extendedVoxels[6]];
                    }

                    if (
                        vertexIndices.count(
                            extendedVoxels[7]
                        ) == 0
                    ) {
                        triangulation.vertices.push_back(
                            cubeVertices[7]
                        );
                        indices[7] = (int) triangulation.vertices.size() - 1;
                        vertexIndices[extendedVoxels[7]] = indices[7];
                    } else {
                        // Point is already in the array
                        indices[7] = vertexIndices[extendedVoxels[7]];
                    }

                    // Add triangles for this face
                    assert(
                        indices[4] >= 0 &&
                        indices[5] >= 0 &&
                        indices[6] >= 0 &&
                        indices[7] >= 0
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[4], indices[5], indices[6]
                        )
                    );
                    triangulation.triangles.push_back(
                        Triangulation::Triangle(
                            indices[4], indices[6], indices[7]
                        )
                    );
                } // end if (... BOTTOM_FACE
            } // end for (ix...
        } // end for (iy...
    } // end for (slice...
}
