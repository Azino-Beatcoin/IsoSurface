#include <QFileDialog>
#include <QMessageBox>
#include "mainwindow.h"
#include "drawarea.h"
#include "exp3ddlg.h"
#include "ui_exp3ddlg.h"

Exp3DDlg::Exp3DDlg(QWidget *parent) :
    QDialog(parent),
    dddFormatGroup(),
    stlFormatGroup(),
    exportPath(),
    exportFormat(EXPORT_VRML),
    stlFormat(STL_FORMAT_BINARY),
    saved(false),
    ui(new Ui::Exp3DDlg)
{
    ui->setupUi(this);
    dddFormatGroup.addButton(ui->radioVRML);
    dddFormatGroup.addButton(ui->radioSTL);
    dddFormatGroup.addButton(ui->radioObj);
    stlFormatGroup.addButton(ui->radioSTLBinary);
    stlFormatGroup.addButton(ui->radioSTLAscii);
    updateState();
}

Exp3DDlg::~Exp3DDlg()
{
    delete ui;
}

void Exp3DDlg::on_browseButton_clicked()
{
    QFileDialog fd(this);
    fd.setFileMode(QFileDialog::AnyFile);

    //... fd.setNameFilter("VRML Files (*.wrl);;STL Files (*.stl)");
    if (exportFormat == EXPORT_VRML) {
        LExportVRML: ;
        fd.setNameFilter("VRML Files (*.wrl)");
        fd.setDefaultSuffix(".wrl");
    } else if (exportFormat == EXPORT_STL) {
        fd.setNameFilter("STL Files (*.stl)");
        fd.setDefaultSuffix(".stl");
    } else if (exportFormat == EXPORT_OBJ) {
        fd.setNameFilter("Wavefront OBJ Files (*.obj)");
        fd.setDefaultSuffix(".obj");
    } else {
        Q_ASSERT(false);
        goto LExportVRML;
    }
    fd.setAcceptMode(QFileDialog::AcceptSave);
    QStringList fileNames;
    QString path;
    if (fd.exec()) {
        fileNames = fd.selectedFiles();
        if (fileNames.size() > 0) {
            ui->path->setText(fileNames[0]);
            path = fileNames[0];
        }
        if (exportFormat == EXPORT_VRML) {
            if (!path.isEmpty()) {
                QApplication::setOverrideCursor(Qt::WaitCursor);
                if (
                    drawArea->triangulation.saveVRML(
                        path, 1.
                    )
                ) {
                    QApplication::restoreOverrideCursor();
                    saved = true;
                    accept();
                    return;
                } else {
                    QApplication::restoreOverrideCursor();
                    QMessageBox::warning(
                        this, "Save Error", "Error in triangulation saving."
                    );
                    saved = false;
                }
            }
        } else if (exportFormat == EXPORT_STL) {
            if (!path.isEmpty()) {
                QApplication::setOverrideCursor(Qt::WaitCursor);
                if (stlFormat == STL_FORMAT_BINARY) {
                    if (
                        drawArea->triangulation.saveSTLbinary(
                            path, 1.
                        )
                    ) {
                        QApplication::restoreOverrideCursor();
                        saved = true;
                        accept();
                        return;
                    } else {
                        QApplication::restoreOverrideCursor();
                        QMessageBox::warning(
                            this, "Save Error", "Error in triangulation saving."
                        );
                        saved = false;
                    }
                } else {
                    if (
                        drawArea->triangulation.saveSTLascii(
                            path, 1.
                        )
                    ) {
                        QApplication::restoreOverrideCursor();
                        saved = true;
                        accept();
                        return;
                    } else {
                        QApplication::restoreOverrideCursor();
                        QMessageBox::warning(
                            this, "Save Error", "Error in triangulation saving."
                        );
                        saved = false;
                    }
                }
            }
        } else if (exportFormat == EXPORT_OBJ) {
            if (!path.isEmpty()) {
                QApplication::setOverrideCursor(Qt::WaitCursor);
                if (
                    !drawArea->triangulation.saveOBJ(
                        path, 1.
                    )
                ) {
                    QApplication::restoreOverrideCursor();
                    QMessageBox::warning(
                        this, "Save Error", "Error in triangulation saving."
                    );
                    saved = false;
                } else {
                    QApplication::restoreOverrideCursor();
                    saved = true;
                    accept();
                    return;
                }
            }
        } else {
            Q_ASSERT(false);
        }
    }
}

void Exp3DDlg::on_radioVRML_clicked()
{
    exportFormat = EXPORT_VRML;
    ui->radioVRML->setChecked(true);
    ui->radioSTL->setChecked(false);
    ui->radioObj->setChecked(false);
    updateState();
}

void Exp3DDlg::on_radioSTL_clicked()
{
    exportFormat = EXPORT_STL;
    ui->radioSTL->setChecked(true);
    ui->radioVRML->setChecked(false);
    ui->radioObj->setChecked(false);
    updateState();
}

void Exp3DDlg::on_radioSTLBinary_clicked()
{
    if (ui->radioSTLBinary->isChecked())
        stlFormat = STL_FORMAT_BINARY;
    else
        stlFormat = STL_FORMAT_ASCII;
    updateState();
}

void Exp3DDlg::on_radioSTLAscii_clicked()
{
    if (ui->radioSTLAscii->isChecked())
        stlFormat = STL_FORMAT_ASCII;
    else
        stlFormat = STL_FORMAT_BINARY;
    updateState();
}

void Exp3DDlg::on_buttonBox_accepted()
{
    QString path = ui->path->text();
    if (!path.isEmpty()) {
        if (
            drawArea->triangulation.saveVRML(path)
        ) {
            saved = true;
            accept();
            return;
        }
    }
    QMessageBox::warning(
        this, "Save Error", "Error in triangulation saving."
    );
    saved = false;
}

void Exp3DDlg::on_buttonBox_rejected()
{
}

void Exp3DDlg::updateState() {
    bool stl = (exportFormat == EXPORT_STL);
    // ui->radioVRML->setChecked(vrml);
    // ui->radioSTL->setChecked(!vrml);

    bool stlBinary = (stlFormat == STL_FORMAT_BINARY);
    ui->radioSTLBinary->setChecked(stlBinary);
    ui->radioSTLAscii->setChecked(!stlBinary);

    ui->radioSTLBinary->setEnabled(stl);
    ui->radioSTLAscii->setEnabled(stl);
}

void Exp3DDlg::accept() {
    if (!saved) {
        return;
    }
    QDialog::accept();
}

void Exp3DDlg::on_radioObj_clicked()
{
    exportFormat = EXPORT_OBJ;
    ui->radioObj->setChecked(true);
    ui->radioVRML->setChecked(false);
    ui->radioSTL->setChecked(false);
    updateState();
}
