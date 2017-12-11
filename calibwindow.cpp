#include "calibwindow.h"
#include "ui_mainwindow.h"
#include <map>
#include <QMessageBox>

CalibWindow::CalibWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::CalibWindow)
{
    setFixedSize(400, 300);
    ui->setupUi(this);
    setCentralWidget(ui->gridLayoutWidget);
    setupPeaksTable();

    std::vector<std::string> vStrCalTypes = Calibration::supportedTypes();
    for(size_t i = 0; i < vStrCalTypes.size(); ++i)
    {
        ui->comboBox->addItem(QString::fromStdString(vStrCalTypes[i]));
    }
}

CalibWindow::~CalibWindow()
{
    delete ui;
}

void CalibWindow::setupPeaksTable()
{
    ui->tableWidget->setColumnCount(4);
    ui->tableWidget->setHorizontalHeaderItem(0, new QTableWidgetItem("Time"));
    ui->tableWidget->setHorizontalHeaderItem(1, new QTableWidgetItem("Theor. mass"));
    ui->tableWidget->setHorizontalHeaderItem(2, new QTableWidgetItem("Calc. mass"));
    ui->tableWidget->setHorizontalHeaderItem(3, new QTableWidgetItem("delta (ppm)"));

    connect(ui->pushButton, SIGNAL(clicked()), this, SLOT(setPeaksNumber()));
    connect(ui->pushButton_2, SIGNAL(clicked()), this, SLOT(calibrate()));
}

void CalibWindow::setPeaksNumber()
{
    int nPeaks = ui->lineEdit->text().toInt();
    ui->tableWidget->setRowCount(nPeaks);
}

void CalibWindow::calibrate()
{
    std::set<double> ms;
    std::map<double, double> ps;
    for(int i = 0; i < ui->tableWidget->rowCount(); ++i)
    {
        double fTime = ui->tableWidget->item(i, 0)->text().toDouble();
        double fMass = ui->tableWidget->item(i, 1)->text().toDouble();
        ps.emplace(std::make_pair(fTime, 1.0));
        ms.emplace(fMass);
    }
    m_pCalibration.reset
    (
        Calibration::create
        (
            Calibration::strToType(ui->comboBox->currentText().toStdString()),
            ms, ps
        )
    );
    if(m_pCalibration->isValid())
    {
        for(int i = 0; i < ui->tableWidget->rowCount(); ++i)
        {
            double fMass = ui->tableWidget->item(i, 1)->text().toDouble();
            double fTime = ui->tableWidget->item(i, 0)->text().toDouble();
            double fMassCalc = m_pCalibration->timeToMass(fTime);
            ui->tableWidget->setItem(i, 2, new QTableWidgetItem(QString::number(fMassCalc)));
            double fPrecision = (fMassCalc - fMass) / fMass * 1e6;
            ui->tableWidget->setItem(i, 3, new QTableWidgetItem(QString::number((int)fPrecision)));
        }
    }
    else
    {
        msg("Errors occured during the calibration");
    }
}

void CalibWindow::msg(QString strMsg)
{
    QMessageBox::warning(Q_NULLPTR, "Calibration message", strMsg);
}
