#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "calibration.h"

namespace Ui {
class CalibWindow;
}

class CalibWindow : public QMainWindow
{
    Q_OBJECT

    void setupPeaksTable();

    Q_SLOT void setPeaksNumber();

    Q_SLOT void calibrate();

    Q_SLOT void msg(QString strMsg);
public:
    explicit CalibWindow(QWidget *parent = 0);
    ~CalibWindow();

private:

    Ui::CalibWindow *ui;

    QScopedPointer<Calibration> m_pCalibration;
};

#endif // MAINWINDOW_H
