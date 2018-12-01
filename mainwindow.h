#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "dlt.hpp"
#include "qcustomplot.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    void addPoint(double x, double y);
    void clearData();
    void plotPolygon();
    void plot();

private slots:
    void on_btAdd_clicked();

    void on_btClear_clicked();

    void clickedGraph(QMouseEvent *event);


    void on_btDLT_clicked();

    void on_btNaive_clicked();

    void on_btNormalizedDLT_clicked();

    void on_btTransform_clicked();

    void on_btNewPoint_clicked();

    void on_btUndo_clicked();

private:
    Ui::MainWindow *ui;

    QVector<double> qv_x, qv_y;
    QCPCurve *polyIn, *polyOut;

    Eigen::Matrix3f P;
    Eigen::UniformScaling<float> S;

    std::vector<Eigen::Vector3f> inPts;
    std::vector<Eigen::Vector3f> outPts;

    void initPoints(std::vector<Eigen::Vector3f> &inPts, std::vector<Eigen::Vector3f> &outPts);
    void printEigenMatrix(Eigen::Matrix3f P);
};

#endif // MAINWINDOW_H
