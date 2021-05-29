#include <QApplication>
#include "mainwindow.h"
#include "drawarea.h"
#include "voxelset.h"
#include "voxfunc.h"
#include "R3Graph.h"

using namespace R3Graph;

void MainWindow::on_regionGrowButton_clicked()
{
    QApplication::setOverrideCursor(Qt::WaitCursor);
    //... setCursor(QCursor(Qt::WaitCursor));
    double dx = (XMAX - XMIN)/double(NUM_VOXELS_X);
    double dy = (YMAX - YMIN)/double(NUM_VOXELS_Y);
    double dz = (ZMAX - ZMIN)/double(NUM_VOXELS_Z);

    double dx2 = dx/2.;
    double dy2 = dy/2.;
    double dz2 = dz/2.;

    Voxel seed;
    bool seedFound = false;

    R3Graph::R3Point p(XMIN, YMIN, ZMIN);
    for (int ix = 0; ix < NUM_VOXELS_X; ++ix) {
        for (int iy = 0; iy < NUM_VOXELS_Y; ++iy) {
            for (int iz = 0; iz < NUM_VOXELS_Z; ++iz) {
                R3Graph::R3Point p(
                    XMIN + dx2 + ix*dx,
                    YMIN + dy2 + iy*dy,
                    ZMIN + dz2 + iz*dz
                );
                if (spaceFunc(p) > 0.) {
                    seed = Voxel(iz, I2Point(ix, iy));
                    seedFound = true;
                }
            }
        }
    }
    if (!seedFound)
        return;

    VoxelBox voxelBox(
        Voxel(0, I2Point(0, 0)),
        NUM_VOXELS_X, NUM_VOXELS_Y, NUM_VOXELS_Z
    );

    if (voxelSet == 0) {
        voxelSet = new VoxelSet();
    }
    voxelSet->initialize(
        NUM_VOXELS_X, NUM_VOXELS_Y, NUM_VOXELS_Z
    );

    detectVoxelSet(
        &voxelFunc,
        0.,                     // threshold,
        voxelBox,
        seed,
        *voxelSet
    );

    voxelSet->computeVoxelBox();
    computeTriangulationOfVoxelSet(
        drawArea->triangulation,
        *voxelSet,
        R3Point(XMIN, YMIN, ZMIN),
        dx, dy, dz
    );
    drawArea->update();
    QApplication::restoreOverrideCursor();
}
