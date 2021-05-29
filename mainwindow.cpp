#include <cmath>
#include <string>
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "drawarea.h"
#include <QFileDialog>
#include "triFunc.h"
#include "skala.h"
#include "voxelset.h"
#include "voxfunc.h"
#include "R3Graph.h"
#include "QtLL1Calc/calc.h"
#include "QtLL1Calc/lparser.h"
#include "exp3ddlg.h"
#include <QtDebug>
#include <QElapsedTimer>

using namespace R3Graph;

MainWindow* mainWindow = 0;
DrawArea* drawArea = 0;

static double f(const R2Point&);
static double f3D(const R3Point&);
static double torus3D(const R3Point& p);

/*
static double f(const R2Point& p) {
    double r = R2Point(0., 0.).distance(p);
    return 2.*cos(4.*r)/(1. + 0.2*r*r);
}
*/

static double f(const R2Point& p) {
    if (mainWindow == 0 || !mainWindow->parsing2DDone)
        return 0.;
    bool defined;
    double v = mainWindow->calc2D.calculate(defined, p.x, p.y);
    if (!defined)
        v = 0.;
    return v;
}

static double f3D(const R3Point& p) {
    if (mainWindow == 0 || !mainWindow->parsing3DDone)
        return 0.;
    bool defined;
    double v = mainWindow->calc3D.calculate(defined, p.x, p.y, p.z);
    if (!defined)
        v = 0.;
    return v;
}

static double torus3D(const R3Point& p) {
    /*
    double R = R3Point(0., 0., 0.).distance(p);
    return R - 3.;
    // return (p.x*p.y - p.z);
    */

    /*
    R2Point q(p.x, p.y);
    double r = R2Point(0., 0.).distance(q);
    return 2.*cos(4.*r)/(1. + 0.2*r*r) - p.z;
    */

    // Slanted torus
    R3Vector n(-2., -1., 1);
    n.normalize();
    R3Vector ez = n;

    R3Vector ey(1., 2., 4.);
    ey.normalize();

    R3Vector ex = ey.vectorProduct(ez);
    ex.normalize();

    R3Point o(0., 0., 0.);  // Center
    double R = 4.;
    double r = 1.;

    R3Vector v = p - o;
    R3Vector proj_v = ex*(v*ex) + ey*(v*ey);
    double l = v.length();
    if (l < (R - r)/2.)
        return (-(R - r)/2.);
    proj_v.normalize();
    R3Point t = o + proj_v*R;
    double d = t.distance(p);
    return (r - d);
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    triangulationPath(),
    xMin(-8.),
    xMax(8.),
    yMin(-8.),
    yMax(8.),
    zMin(-8.),
    zMax(8.),
    voxelSet(0),
    shading(FLAT_SHADING),
    calc2D(),
    calc3D(),
    parsing2DDone(false),
    parsing3DDone(false),
    smooth_lambda(0.30),
    smooth_mu(0.31),
    smooth_iterations(1),
    ui(new Ui::MainWindow)
{
    mainWindow = this;
    ui->setupUi(this);

    setDimensionsText();

    QString tmp;
    ui->smoothLambda->setText(tmp.sprintf("%.3f", smooth_lambda));
    ui->smoothMu->setText(tmp.sprintf("%.3f", smooth_mu));
    ui->iterations->setText(tmp.sprintf("%d", smooth_iterations));
}

MainWindow::~MainWindow()
{
    if (voxelSet != 0)
        delete voxelSet;
    delete ui;
}

void MainWindow::setDimensionsText() {
    QString tmp;
    ui->xMin->setText(tmp.sprintf("%.2f", xMin));
    ui->xMax->setText(tmp.sprintf("%.2f", xMax));
    ui->yMin->setText(tmp.sprintf("%.2f", yMin));
    ui->yMax->setText(tmp.sprintf("%.2f", yMax));
    ui->zMin->setText(tmp.sprintf("%.2f", zMin));
    ui->zMax->setText(tmp.sprintf("%.2f", zMax));
}

void MainWindow::on_path_returnPressed()
{
    triangulationPath = ui->path->text();
    on_loadButton_clicked();
}

void MainWindow::on_browseButton_clicked()
{
    QFileDialog openFileDlg;
    QString fileName = QFileDialog::getOpenFileName(
        this, QString("Triangulation File")
    );
    if (fileName == QString())
        return;
    ui->path->setText(fileName);
    on_loadButton_clicked();
}

void MainWindow::on_loadButton_clicked() {
    setCursor(QCursor(Qt::WaitCursor));

    triangulationPath = ui->path->text();
    bool res = drawArea->triangulation.load(
        triangulationPath.toStdString().c_str()
    );

    if (res) {
        // Set dimensions
        double extendFactor = 1.75;
        const R3Vector& boxSize = drawArea->triangulation.box.size;
        R3Point boxCenter = drawArea->triangulation.box.origin +
            drawArea->triangulation.box.size * 0.5;
        double s = boxSize.x;
        if (boxSize.y > s)
            s = boxSize.y;
        if (boxSize.z > s)
            s = boxSize.z;
        s *= extendFactor;

        xMin = boxCenter.x - s/2.; xMax = boxCenter.x + s/2.;
        yMin = boxCenter.y - s/2.; yMax = boxCenter.y + s/2.;
        zMin = boxCenter.z - s/2.; zMax = boxCenter.z + s/2.;
        setDimensionsText();

        drawArea->xmin = xMin;
        drawArea->xmax = xMax;
        drawArea->ymin = yMin;
        drawArea->ymax = yMax;
        drawArea->onResizeGL();
    }

    setCursor(QCursor(Qt::ArrowCursor));

    drawArea->update();
}

void MainWindow::on_plot2DButton_clicked()
{
    setCursor(QCursor(Qt::WaitCursor));

    double xmin = XMIN;
    double xmax = XMAX;
    double ymin = YMIN;
    double ymax = YMAX;
    double zmin = ZMIN;
    double zmax = ZMAX;
    xMin = xmin; xMax = xmax;
    yMin = ymin; yMax = ymax;
    zMin = zmin; zMax = zMax;
    setDimensionsText();
    setCoordinates(false);

    plot2DFunction(
        &f,
        R2Rectangle(
            R2Point(xMin, yMin),
            xMax - xMin, yMax - yMin
        ),
        128, 128,
        drawArea->triangulation
    );

    drawArea->xmin = xMin;
    drawArea->xmax = xMax;
    drawArea->ymin = yMin;
    drawArea->ymax = yMax;
    drawArea->onResizeGL();

    setCursor(QCursor(Qt::ArrowCursor));

    drawArea->update();
}

void MainWindow::on_saveButton_clicked()
{
    triangulationPath = ui->path->text();
    if (triangulationPath == "")
        return;
    drawArea->triangulation.save(
        triangulationPath.toStdString().c_str()
    );
}

void MainWindow::on_plot3DButton_clicked()
{
    setCursor(QCursor(Qt::WaitCursor));

    // double ymin = xMin;
    // double ymax = xMax;
    // double zmin = xMin;
    // double zmax = xMax;
    double xmin = XMIN;
    double xmax = XMAX;
    double ymin = YMIN;
    double ymax = YMAX;
    double zmin = ZMIN;
    double zmax = ZMAX;
    xMin = xmin; xMax = xmax;
    yMin = ymin; yMax = ymax;
    zMin = zmin; zMax = zMax;
    setCoordinates(false);

    setDimensionsText();

    skalaMethod(
        &f3D,
        // &spaceFunc,
        R3Box(
            R3Point(xmin, ymin, zmin),
            R3Vector(
                xmax-xmin, ymax-ymin, zmax-zmin
            )
        ),
        NUM_VOXELS_X, NUM_VOXELS_Y, NUM_VOXELS_Z,
        drawArea->triangulation
    );

    drawArea->xmin = xMin;
    drawArea->xmax = xMax;
    drawArea->ymin = yMin;
    drawArea->ymax = yMax;
    drawArea->onResizeGL();

    setCursor(QCursor(Qt::ArrowCursor));
    drawArea->update();
}

void MainWindow::on_redrawButton_clicked()
{
    setCoordinates();
    drawArea->update();
}

void MainWindow::setCoordinates(bool copyFromTextboxes /* = true */) {
    if (copyFromTextboxes) {
        xMin = ui->xMin->text().toDouble();
        xMax = ui->xMax->text().toDouble();
        yMin = ui->yMin->text().toDouble();
        yMax = ui->yMax->text().toDouble();
        zMin = ui->zMin->text().toDouble();
        zMax = ui->zMax->text().toDouble();
    }
    double dx = fabs(xMax - xMin)*1.2;
    double xCenter = (xMax + xMin)/2.;
    drawArea->xmin = xCenter - dx/2.;
    drawArea->xmax = xCenter + dx/2.;
    drawArea->onResizeGL();
}

void MainWindow::on_radioWireframe_clicked()
{
    if (ui->radioWireframe->isChecked()) {
        shading = WIREFRAME_SHADING;
        ui->radioFlatShading->setChecked(false);
        ui->radioSmoothShading->setChecked(false);
    } else if (ui->radioFlatShading->isChecked()) {
        shading = FLAT_SHADING;
        ui->radioWireframe->setChecked(false);
        ui->radioSmoothShading->setChecked(false);
    } else if (ui->radioSmoothShading->isChecked()) {
        shading = SMOOTH_SHADING;
        ui->radioWireframe->setChecked(false);
        ui->radioFlatShading->setChecked(false);
    }
    drawArea->update();
}

void MainWindow::on_radioFlatShading_clicked()
{
    on_radioWireframe_clicked();
}

void MainWindow::on_radioSmoothShading_clicked()
{
    on_radioWireframe_clicked();
}

void MainWindow::on_initialViewButton_clicked()
{
    drawArea->initializeView();
}

void MainWindow::on_fxy_editingFinished()
{
    char scanLine[1024];
    QString txt = ui->fxy->text();
    const char* str = txt.toUtf8().constData();
    strncpy(scanLine, str, 1022);
    scanLine[1022] = 0;

    // Call the parser
    std::string errorMessage;

    parsing2DDone = lparser::parseFormula(
        scanLine, calc2D, errorMessage
    );

    if (!parsing2DDone) {
        ui->twoDLabel->setText(errorMessage.c_str());
    } else {
        ui->twoDLabel->setText("Parsing done: OK");
    }

    // drawArea->onDraw();
}

void MainWindow::on_fxyz_editingFinished()
{
    char scanLine[1024];
    QString txt = ui->fxyz->text();
    const char* str = txt.toUtf8().constData();
    strncpy(scanLine, str, 1022);
    scanLine[1022] = 0;

    // Call the parser
    std::string errorMessage;

    parsing3DDone = lparser::parseFormula(
        scanLine, calc3D, errorMessage
    );

    if (!parsing3DDone) {
        ui->threeDLabel->setText(errorMessage.c_str());
    } else {
        ui->threeDLabel->setText("Parsing done: OK");
    }
}

void MainWindow::on_torusTriangulationButton_clicked()
{
    setCursor(QCursor(Qt::WaitCursor));

    // double ymin = xMin;
    // double ymax = xMax;
    // double zmin = xMin;
    // double zmax = xMax;
    double xmin = XMIN;
    double xmax = XMAX;
    double ymin = YMIN;
    double ymax = YMAX;
    double zmin = ZMIN;
    double zmax = ZMAX;
    xMin = xmin; xMax = xmax;
    yMin = ymin; yMax = ymax;
    zMin = zmin; zMax = zMax;
    setCoordinates(false);
    setDimensionsText();

    skalaMethod(
        // &f3D,
        &spaceFunc,
        R3Box(
            R3Point(xmin, ymin, zmin),
            R3Vector(
                xmax-xmin, ymax-ymin, zmax-zmin
            )
        ),
        NUM_VOXELS_X, NUM_VOXELS_Y, NUM_VOXELS_Z,
        drawArea->triangulation
    );

    drawArea->xmin = xMin;
    drawArea->xmax = xMax;
    drawArea->ymin = yMin;
    drawArea->ymax = yMax;
    drawArea->onResizeGL();

    setCursor(QCursor(Qt::ArrowCursor));
    drawArea->update();
}

void MainWindow::on_uniformLaplaceButton_clicked()
{
    setCursor(QCursor(Qt::WaitCursor));

    for (int i = 0; i < smooth_iterations; ++i) {
        drawArea->triangulation.uniformLaplaceSmoothing(
            smooth_lambda
        );
    }

    setCursor(QCursor(Qt::ArrowCursor));
    drawArea->update();
}

void MainWindow::on_cotangentLaplaceButton_clicked()
{
    setCursor(QCursor(Qt::WaitCursor));

    for (int i = 0; i < smooth_iterations; ++i) {
        drawArea->triangulation.cotangentLaplaceSmoothing(
            smooth_lambda
        );
    }

    setCursor(QCursor(Qt::ArrowCursor));
    drawArea->update();
}

void MainWindow::on_taubinUniformButton_clicked()
{
    setCursor(QCursor(Qt::WaitCursor));

    drawArea->triangulation.taubinSmoothing(
        smooth_iterations,
        smooth_lambda, smooth_mu,
        false                   // Not using cotangent smoothing
    );

    setCursor(QCursor(Qt::ArrowCursor));
    drawArea->update();
}

void MainWindow::on_taubinCotangentButton_clicked()
{
    setCursor(QCursor(Qt::WaitCursor));

    drawArea->triangulation.taubinSmoothing(
        smooth_iterations,
        smooth_lambda, smooth_mu,
        true                    // Using cotangent smoothing
    );

    setCursor(QCursor(Qt::ArrowCursor));
    drawArea->update();
}

void MainWindow::on_checkEdgesButton_clicked()
{
    setCursor(QCursor(Qt::WaitCursor));

    int numManifoldEdges;
    int numBorderEdges;
    int numNonManifoldEdges;

    drawArea->triangulation.checkEdges(
        numManifoldEdges,
        numBorderEdges,
        numNonManifoldEdges
    );

    qDebug() << "Manifold edges num: " << numManifoldEdges;
    qDebug() << "Border edges num: " << numBorderEdges;
    qDebug() << "NonManifold edges num: " << numNonManifoldEdges;
    qDebug() << "Triangles num: " << drawArea->triangulation.triangles.size();
    qDebug() << "Verices num: " << drawArea->triangulation.vertices.size();

    setCursor(QCursor(Qt::ArrowCursor));
}

void MainWindow::on_anniTriButton_clicked()
{
    setCursor(QCursor(Qt::WaitCursor));

    int numRemovedEdges;
    int numRemovedCompleteTriangles;
    int numRemovedTriangles;

    drawArea->triangulation.annihilateSmallTriangles(
        numRemovedEdges,
        numRemovedCompleteTriangles,
        numRemovedTriangles
    );

    setCursor(QCursor(Qt::ArrowCursor));
    drawArea->update();
}

void MainWindow::on_simplifyButton_clicked()
{
    setCursor(QCursor(Qt::WaitCursor));

    int removedVerticesNum = 0;
    int removedEdgesNum = 0;
    int removedTrianglesNum = 0;

    ui->removedVerticesNumLine->setText( QString::number(removedVerticesNum) );
    ui->removedEdgesNumLine->setText( QString::number(removedEdgesNum) );
    ui->removedTrianglesNumLine->setText( QString::number(removedTrianglesNum) );

    int iterationsTotal = ui->simplifyIterNumLine->text().toInt();
    int iterationNum = 0;

    QElapsedTimer timer;
    timer.start();

    while(iterationsTotal == 0 || iterationNum < iterationsTotal)
    {
        removedVerticesNum = 0;
        removedEdgesNum = 0;
        removedTrianglesNum = 0;

        drawArea->triangulation.simplify(
            removedVerticesNum,
            removedEdgesNum,
            removedTrianglesNum,
            ui->epsLine->text().toDouble()
        );

        ui->removedVerticesNumLine->setText(
            QString::number(ui->removedVerticesNumLine->text().toInt() + removedVerticesNum)
        );
        ui->removedEdgesNumLine->setText(
            QString::number(ui->removedEdgesNumLine->text().toInt() + removedEdgesNum)
        );
        ui->removedTrianglesNumLine->setText(
            QString::number(ui->removedTrianglesNumLine->text().toInt() + removedTrianglesNum)
        );

        ++iterationNum;

        if( removedVerticesNum == 0 && removedEdgesNum == 0 && removedTrianglesNum == 0 )
            break;
    }

    qDebug() << "Time: " << timer.elapsed() << "milliseconds";

    setCursor(QCursor(Qt::ArrowCursor));
    drawArea->update();
}

void MainWindow::on_mySimplifyButton_clicked()
{
    setCursor(QCursor(Qt::WaitCursor));

    int removedVerticesNum = 0;
    int removedEdgesNum = 0;
    int removedTrianglesNum = 0;

    ui->removedVerticesNumLine->setText( QString::number(removedVerticesNum) );
    ui->removedEdgesNumLine->setText( QString::number(removedEdgesNum) );
    ui->removedTrianglesNumLine->setText( QString::number(removedTrianglesNum) );

    int iterationsTotal = ui->simplifyIterNumLine->text().toInt();
    int iterationNum = 0;

    QElapsedTimer timer;
    timer.start();

    while(iterationsTotal == 0 || iterationNum < iterationsTotal)
    {
        removedVerticesNum = 0;
        removedEdgesNum = 0;
        removedTrianglesNum = 0;

        drawArea->triangulation.mySimplify(
            removedVerticesNum,
            removedEdgesNum,
            removedTrianglesNum,
            ui->epsLine->text().toDouble()
        );

        ui->removedVerticesNumLine->setText( QString::number(ui->removedVerticesNumLine->text().toInt() + removedVerticesNum) );
        ui->removedEdgesNumLine->setText( QString::number(ui->removedEdgesNumLine->text().toInt() + removedEdgesNum) );
        ui->removedTrianglesNumLine->setText( QString::number(ui->removedTrianglesNumLine->text().toInt() + removedTrianglesNum) );

        iterationNum++;

        if( removedVerticesNum == 0 && removedEdgesNum == 0 && removedTrianglesNum == 0 )
            break;
    }

    qDebug() << "Time: " << timer.elapsed() << "milliseconds";

    setCursor(QCursor(Qt::ArrowCursor));
    drawArea->update();
}

void MainWindow::on_smoothLambda_editingFinished()
{
    double l = ui->smoothLambda->text().toDouble();
    if (l > 0.)
        smooth_lambda = l;
}

void MainWindow::on_smoothMu_editingFinished()
{
    double m = ui->smoothMu->text().toDouble();
    if (m > 0.)
        smooth_mu = m;
}

void MainWindow::on_iterations_editingFinished()
{
    int n = ui->iterations->text().toInt();
    if (n > 0)
        smooth_iterations = n;
}

void MainWindow::on_expotrButton_clicked()
{
    Exp3DDlg export3DDlg(this);
    export3DDlg.exec();
}
