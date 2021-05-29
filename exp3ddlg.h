#ifndef EXP3DDLG_H
#define EXP3DDLG_H

#include <QDialog>
#include <QButtonGroup>

const int EXPORT_VRML = 0;
const int EXPORT_STL = 1;
const int EXPORT_OBJ = 2;       // Wavefront OBJ file

const int STL_FORMAT_BINARY = 0;
const int STL_FORMAT_ASCII = 1;

namespace Ui {
class Exp3DDlg;
}

class Exp3DDlg : public QDialog
{
    Q_OBJECT

public:
    explicit Exp3DDlg(QWidget *parent = 0);
    ~Exp3DDlg();
    void updateState();

public slots:
    void accept();      // Virtual method

public:
    QButtonGroup dddFormatGroup;
    QButtonGroup stlFormatGroup;
    QString exportPath;
    int exportFormat;   // EXPORT_VRML / EXPORT_STL / EXPORT_OBJ
    int stlFormat;      // STL_FORMAT_BINARY / STL_FORMAT_ASCII
    bool saved;

private slots:
    void on_browseButton_clicked();

    void on_radioVRML_clicked();

    void on_radioSTL_clicked();

    void on_radioSTLBinary_clicked();

    void on_radioSTLAscii_clicked();

    void on_buttonBox_accepted();

    void on_buttonBox_rejected();

    void on_radioObj_clicked();

private:
    Ui::Exp3DDlg *ui;
};

#endif // EXP3DDLG_H
