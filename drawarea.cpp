#include <QMouseEvent>
#include <QWheelEvent>
#include <cassert>
#include <cmath>
#include "drawarea.h"
#include "mainwindow.h"
//... #include "GL/glu.h"

using namespace R3Graph;

DrawArea::DrawArea(QWidget *parent):
    QGLWidget(parent),
    triangulation(),
    xmin(-8.),
    xmax(8.),
    ymin(-8.),
    ymax(8.),
    zmin(-8.),
    zmax(8.),
    alpha(ALPHA0),
    beta(BETA0),
    // shiftX((xmin + xmax)/2.),
    // shiftY((ymin + ymax)/2.),
    // shiftZ((zmin + zmax)/2.),
    shiftX(0.),
    shiftY(0.),
    shiftZ(0.),
    mouseX(-1),     // Undefined
    mouseY(-1)      // Undefined
{
    setMouseTracking(true);
    drawArea = this;
}

void DrawArea::initializeGL() {
    /*
    glDisable(GL_COLOR_MATERIAL);
    glEnable(GL_BLEND);
    glEnable(GL_POLYGON_SMOOTH);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClearColor(0., 0.1, 0.2, 1.);
    */

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    glEnable(GL_NORMALIZE);
    glDepthFunc(GL_LEQUAL);

    // Lighting
    glEnable(GL_LIGHT0);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

    GLfloat pos[4];
    // pos[0] = 0.; pos[1] = 0.; pos[2] = 1.; pos[3] = 0.; // Directional light
    pos[0] = -0.5; pos[1] = 1.; pos[2] = 1.; pos[3] = 0.; // Directional light
    glLightfv(GL_LIGHT0, GL_POSITION, pos);

    GLfloat light[4];
    light[0] = 0.25; light[1] = 0.25; light[2] = 0.25; light[3] = 1.;
    glLightfv(GL_LIGHT0, GL_AMBIENT, light);

    light[0] = 1.; light[1] = 1.; light[2] = 1.; light[3] = 1.;
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void DrawArea::resizeGL(int w, int h) {
    onResizeGL(w, h);
}

void DrawArea::onResizeGL(int w /* = (-1) */, int h /* = (-1) */) {
    if (w < 0) {
        w = width();
        h = height();
    }
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if (w == 0)
        w = 1;
    if (h == 0)
        h = 1;
    double aspect = (double) w / (double) h;
    assert(aspect > 0.);

    // R3Point center = triangulation.box.origin + triangulation.box.size*0.5;
    double dy = (xmax - xmin)/aspect;
    // ymin = center.y - dy/2.;
    // ymax = center.y + dy/2.;
    double centerY = (ymin + ymax)/2.;
    ymin = centerY - dy/2.;
    ymax = centerY + dy/2.;
    mainWindow->yMin = ymin;
    mainWindow->yMax = ymax;
    mainWindow->setDimensionsText();

    double depth = ymax - ymin;
    if (xmax - xmin > depth)
        depth = xmax - xmin;
    double centerZ = (zmin + zmax)/2.;
    zmin = centerZ - depth*10.;
    zmax = centerZ + depth*10.;

    glOrtho(
        xmin, xmax,
        ymin, ymax,
        // -10. * depth, 10. * depth
        zmin, zmax
    );
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void DrawArea::paintGL() {
    // Rotate the scene
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glEnable(GL_LIGHT0);        // Turn on the first light source

    // Shift a model center
    glTranslatef(shiftX, shiftY, shiftZ);

    // Rotate a model around its center
    R3Point center = triangulation.box.origin +
        triangulation.box.size*0.5;
    glTranslatef(center.x, center.y, center.z);
    glRotatef(beta, 1., 0., 0.);
    glRotatef(alpha, 0., 1., 0.);
    glTranslatef(-center.x, -center.y, -center.z);

    // glTranslatef(-shiftX, -shiftY, 0.);

    render();           // Draw a scene graph

    // swapBuffers();   // Not needed - done automatically
}

void DrawArea::mousePressEvent(QMouseEvent * /* event */) {
    mouseX = (-1);      // Undefined
    mouseY = (-1);      // Undefined
}

void DrawArea::mouseReleaseEvent(QMouseEvent * /* event */) {
    mouseX = (-1);      // Undefined
    mouseY = (-1);      // Undefined
}

void DrawArea::mouseMoveEvent(QMouseEvent *event) {
    //... printf("%d, %d\n", event->x(), event->y());
    int x = event->x();
    int y = event->y();
    if ((event->buttons() & Qt::LeftButton) != 0) {
        if (mouseX >= 0) {
            int dx = x - mouseX;
            int dy = y - mouseY;
            if (beta > 100. || beta < (-100.))
                dx = -dx;

            alpha += dx * 0.1;
            if (fabs(alpha) > 360.)
                alpha -= ((int) alpha / 360) * 360.;

            beta += dy * 0.1;
            if (fabs(beta) > 360.)
                beta -= ((int) beta / 360) * 360.;

            repaint();
        }
        mouseX = x; mouseY = y;

    } else if ((event->buttons() & Qt::RightButton) != 0) {
        if (mouseX >= 0) {
            int dx = x - mouseX;
            int dy = mouseY - y;
            shiftX += dx * 0.05;
            shiftY += dy * 0.05;

            repaint();
        }
        mouseX = x; mouseY = y;

    } else {
        mouseX = (-1);      // Undefined
        mouseY = (-1);      // Undefined
    }
}

void DrawArea::keyPressEvent(QKeyEvent* event) {
    switch (event->key()) {
    case Qt::Key_Escape:
    case Qt::Key_Q:
        close();
        break;
    case Qt::Key_Space:
        alpha = ALPHA0;
        beta = BETA0;
        repaint();
        break;
    default:
        event->ignore();
        break;
    }
}

void DrawArea::render() {
    GLfloat color[4];

    // glClearColor(0.1, 0.2, 0.4, 1.); // Background color: dark blue
    glClearColor(1., 1., 1., 1.); // Background color: white
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // color[0] = 0.5; color[1] = 1.; color[2] = 0.5;  // Green
    GLfloat gray = 0.6;
    color[0] = gray; color[1] = gray; color[2] = gray;  // Gray scale
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);

    // Draw the triangualtion
    if (
        mainWindow->shading == FLAT_SHADING ||
        mainWindow->shading == SMOOTH_SHADING
    ) {
        glBegin(GL_TRIANGLES);
        for (
            size_t t = 0;
            t < triangulation.triangles.size();
            ++t
        ) {
            const Triangulation::Triangle& tr =
                triangulation.triangles[t];

            for (int i = 0; i < 3; ++i) {
                int vertexIdx = tr[i];
                const Triangulation::Vertex& v =
                    triangulation.vertices[vertexIdx];

                if (mainWindow->shading == FLAT_SHADING) {
                    // Set normal to triangle for flat shading
                    if (i == 0) {
                        R3Point p0 = triangulation.vertices[tr[0]].point;
                        R3Point p1 = triangulation.vertices[tr[1]].point;
                        R3Point p2 = triangulation.vertices[tr[2]].point;
                        R3Vector n = (p1 - p0).vectorProduct(p2 - p0);
                        n.normalize();
                        glNormal3f(
                            n.x, n.y, n.z
                        );
                    }
                } else {
                    assert(mainWindow->shading == SMOOTH_SHADING);
                    // Set a normal to the current point
                    glNormal3f(
                        v.normal.x, v.normal.y, v.normal.z
                    );
                }
                glVertex3f(
                    v.point.x, v.point.y, v.point.z
                );
            }
        }
        glEnd();
    } else {
        assert(mainWindow->shading == WIREFRAME_SHADING);
        gray = 0.3;
        color[0] = gray; color[1] = gray; color[2] = gray;  // Gray scale
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);

        glBegin(GL_LINES);
        for (
            size_t t = 0;
            t < triangulation.triangles.size();
            ++t
        ) {
            const Triangulation::Triangle& tr =
                triangulation.triangles[t];

            for (int i = 0; i < 3; ++i) {
                int j = i + 1;
                if (j == 3)
                    j = 0;
                int vertexIdx0 = tr[i];
                int vertexIdx1 = tr[j];
                const Triangulation::Vertex& v0 =
                    triangulation.vertices[vertexIdx0];
                const Triangulation::Vertex& v1 =
                    triangulation.vertices[vertexIdx1];

                // Set normal to triangle for flat shading
                if (i == 0) {
                    R3Point p0 = triangulation.vertices[tr[0]].point;
                    R3Point p1 = triangulation.vertices[tr[1]].point;
                    R3Point p2 = triangulation.vertices[tr[2]].point;
                    R3Vector n = (p1 - p0).vectorProduct(p2 - p0);
                    n.normalize();
                    glNormal3f(
                        n.x, n.y, n.z
                    );
                }

                glVertex3f(
                    v0.point.x, v0.point.y, v0.point.z
                );
                glVertex3f(
                    v1.point.x, v1.point.y, v1.point.z
                );
            }
        } // end for

        glEnd();
    }
}

void DrawArea::initializeView() {
    alpha = ALPHA0;
    beta = BETA0;
    // shiftX = (xmin + xmax)/2.;
    // shiftY = (ymin + ymax)/2.;
    shiftX = 0.;
    shiftY = 0.;
    shiftZ = 0.;

    update();
}

static const double WHEEL_FACTOR = 1.05;

void DrawArea::wheelEvent(QWheelEvent *event) {
    QPoint delta = event->angleDelta();
    double centerX =
        (mainWindow->xMin + mainWindow->xMax)/2.;
    double centerY =
        (mainWindow->yMin + mainWindow->yMax)/2.;
    double centerZ =
        (mainWindow->zMin + mainWindow->zMax)/2.;
    double zoomFactor = 1.;
    if (delta.y() > 0) {
        zoomFactor *= WHEEL_FACTOR;
    } else if (delta.y() < 0) {
        zoomFactor /= WHEEL_FACTOR;
    } else {
        event->ignore();
    }
    event->accept();
    mainWindow->xMin = centerX +
        (mainWindow->xMin - centerX)*zoomFactor;
    mainWindow->xMax = centerX +
        (mainWindow->xMax - centerX)*zoomFactor;

    mainWindow->yMin = centerY +
        (mainWindow->yMin - centerY)*zoomFactor;
    mainWindow->yMax = centerY +
        (mainWindow->yMax - centerY)*zoomFactor;

    mainWindow->zMin = centerZ +
        (mainWindow->zMin - centerZ)*zoomFactor;
    mainWindow->zMax = centerZ +
        (mainWindow->zMax - centerZ)*zoomFactor;
    mainWindow->setDimensionsText();
    mainWindow->setCoordinates();
    update();
}

