#include "calibwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    CalibWindow w;
    w.show();

    return a.exec();
}
