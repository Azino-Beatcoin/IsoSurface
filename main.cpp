#include "mainwindow.h"
#include <QApplication>
#include <locale.h>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    setlocale(LC_ALL, "C");
    MainWindow w;
    w.show();

    return a.exec();
}
