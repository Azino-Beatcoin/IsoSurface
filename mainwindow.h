#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
// #include "voxelset.h"

#include "QtLL1Calc/calc.h"
#include "QtLL1Calc/lparser.h"

class DrawArea;
class MainWindow;

extern MainWindow* mainWindow;
extern DrawArea* drawArea;

class VoxelSet;

namespace Ui {
class MainWindow;
}

const int WIREFRAME_SHADING = 0;
const int FLAT_SHADING = 1;
const int SMOOTH_SHADING = 2;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    QString triangulationPath;

    double xMin, xMax, yMin, yMax, zMin, zMax;
    VoxelSet* voxelSet;
    int shading; // WIREFRAME_SHADING / FLAT_SHADING / SMOOTH_SHADING

    Calc calc2D;
    Calc calc3D;
    bool parsing2DDone;
    bool parsing3DDone;

    double smooth_lambda;
    double smooth_mu;
    int smooth_iterations;

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    void setCoordinates(bool copyFromTextboxes = true);
    void setDimensionsText();
    bool loadTriangulation(QString path);

private slots:
    void on_path_returnPressed();

    void on_browseButton_clicked();

    void on_loadButton_clicked();

    void on_plot2DButton_clicked();

    void on_saveButton_clicked();

    void on_plot3DButton_clicked();

    void on_redrawButton_clicked();

    void on_regionGrowButton_clicked();

    void on_radioWireframe_clicked();

    void on_radioFlatShading_clicked();

    void on_radioSmoothShading_clicked();

    void on_initialViewButton_clicked();

    void on_fxy_editingFinished();

    void on_fxyz_editingFinished();

    void on_torusTriangulationButton_clicked();

    void on_uniformLaplaceButton_clicked();

    void on_cotangentLaplaceButton_clicked();

    void on_taubinUniformButton_clicked();

    void on_taubinCotangentButton_clicked();

    void on_anniTriButton_clicked();

    void on_simplifyButton_clicked();

    void on_mySimplifyButton_clicked();

    void on_checkEdgesButton_clicked();

    void on_smoothLambda_editingFinished();

    void on_smoothMu_editingFinished();

    void on_iterations_editingFinished();

    void on_expotrButton_clicked();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
