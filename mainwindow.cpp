#include "mainwindow.h"
#include "ui_mainwindow.h"


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow){
        ui->setupUi(this);

        ui->plot->addGraph();
        ui->plot->graph(0)->setScatterStyle(QCPScatterStyle::ssCircle);
        ui->plot->graph(0)->setLineStyle(QCPGraph::lsNone);
        ui->plot->xAxis->setRange(-10, 10);
        ui->plot->yAxis->setRange(-10, 10);
        ui->plot->xAxis->setLabel("x");
        ui->plot->yAxis->setLabel("y");

        connect(ui->plot, SIGNAL(mousePress(QMouseEvent*)), SLOT(clickedGraph(QMouseEvent*)));

        QBrush shadowBrush (QColor(0,0,0), Qt::Dense7Pattern);
        polyIn = new QCPCurve(ui->plot->xAxis, ui->plot->yAxis);
        polyIn->setBrush(shadowBrush);

        polyOut = new QCPCurve(ui->plot->xAxis, ui->plot->yAxis);
        polyOut->setBrush(shadowBrush);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::addPoint(double x, double y)
{
    qv_x.append(x);
    qv_y.append(y);
}

void MainWindow::clearData()
{
    qv_x.clear();
    qv_y.clear();

    polyIn->setData(qv_x, qv_y);
    polyOut->setData(qv_x, qv_y);

    ui->ptOut->clear();
    ui->lbInX->clear();
    ui->lbInY->clear();
    ui->lbOut->clear();
    ui->lbNewPoint1->clear();
    ui->lbNewPoint2->clear();
}

void MainWindow::plotPolygon() {

    QVector<double> qv_xFirst, qv_yFirst;
    qv_xFirst.clear();
    qv_yFirst.clear();

    //first polygon
    for (int i = 0; i < qv_x.size(); i+=2) {
        qv_xFirst.append(qv_x[i]);
        qv_yFirst.append(qv_y[i]);
    }
    qv_xFirst.append(qv_x[0]);
    qv_yFirst.append(qv_y[0]);

    polyIn->setData(qv_xFirst, qv_yFirst);

    qv_xFirst.clear();
    qv_yFirst.clear();

    //second polygon
    for (int i = 1; i < qv_x.size(); i+=2) {
        qv_xFirst.append(qv_x[i]);
        qv_yFirst.append(qv_y[i]);
    }
    qv_xFirst.append(qv_x[1]);
    qv_yFirst.append(qv_y[1]);

    polyOut->setData(qv_xFirst, qv_yFirst);
}

void MainWindow::plot()
{
    ui->plot->graph(0)->setData(qv_x, qv_y);
    ui->plot->replot();
    ui->plot->update();
}

void MainWindow::on_btAdd_clicked()
{
    addPoint(ui->dsX->value(), ui->dsY->value());
    plot();
}

void MainWindow::on_btClear_clicked()
{
    clearData();
    plot();
}

void MainWindow::on_btUndo_clicked()
{
    if(qv_x.size() > 0) {
        qv_x.pop_back();
        qv_y.pop_back();
        plot();
    }
}

void MainWindow::clickedGraph(QMouseEvent *event)
{
   QPoint point = event->pos();
   double x = ui->plot->xAxis->pixelToCoord(point.x());
   double y = ui->plot->yAxis->pixelToCoord(point.y());
   ui->lbInX->setText("x: " + QString::number(x));
   ui->lbInY->setText("y: " + QString::number(y));

   addPoint(x, y);
   plot();

}

void MainWindow::initPoints(std::vector<Eigen::Vector3f>& inPts, std::vector<Eigen::Vector3f>& outPts) {
    inPts.clear();
    outPts.clear();

    float scale = float(ui->vsScale->value()) / 5;
    Eigen::UniformScaling<float>S1(scale);

    S = S1;

    double *tx = qv_x.data();
    double *ty = qv_y.data();
    for (int i = 0; i < qv_x.size(); i+=2) {
        Eigen::Vector3f point = {float(qv_x.at(i)), float(qv_y.at(i)), 1};
        point = S * point;
        inPts.push_back(point);
        tx[i] = double(point.coeff(0));
        ty[i] = double(point.coeff(1));
    }

    for (int i = 1 ; i < qv_x.size(); i+=2) {
        Eigen::Vector3f point = {float(qv_x.at(i)), float(qv_y.at(i)), 1};
        point = S * point;
        outPts.push_back(point);
        tx[i] = double(point.coeff(0));
        ty[i] = double(point.coeff(1));
    }
}

void MainWindow::printEigenMatrix(Eigen::Matrix3f P) {
    ui->ptOut->clear();
    ui->lbOut->setText("Projective matrix P:");

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ui->ptOut->insertPlainText(QString::number(double(P.coeff(i, j))) + "\t\t");
        }
        ui->ptOut->insertPlainText("\n");
    }
}

void MainWindow::on_btDLT_clicked()
{
    ui->lbNewPoint2->clear();
    plotPolygon();

    if(qv_x.size() < 8) {
        ui->ptOut->setPlainText("Enter 8 points");
        return;
    }

    if(qv_x.size() % 2 != 0) {
        ui->ptOut->setPlainText("Some points are missing. Number of points must be even");
        return;
    }

    initPoints(inPts, outPts);

    P = projectiveDLT(inPts, outPts);

    plot();

    printEigenMatrix(P);
    ui->lbNewPoint1->setText("Insert new point:");
}

void MainWindow::on_btNaive_clicked()
{
    ui->lbNewPoint2->clear();
    plotPolygon();

    if(qv_x.size() != 8) {
        ui->ptOut->setPlainText("Enter 8 points");
        return;
    }

    initPoints(inPts, outPts);

    P = naiveProjective4pts(inPts, outPts);

    plot();

    printEigenMatrix(P);
    ui->lbNewPoint1->setText("Insert new point:");
}

void MainWindow::on_btNormalizedDLT_clicked()
{
    ui->lbNewPoint2->clear();

    if(qv_x.size() < 8) {
        ui->ptOut->setPlainText("Enter 8 points");
        return;
    }

    if(qv_x.size() % 2 != 0) {
        ui->ptOut->setPlainText("Some points are missing. Number of points must be even");
        return;
    }

    initPoints(inPts, outPts);
    plotPolygon();

    plot();

    P = normalizedDLT(inPts, outPts);

    printEigenMatrix(P);
    ui->lbNewPoint1->setText("Insert new point:");
}

void MainWindow::on_btTransform_clicked()
{
    double *tx = qv_x.data();
    double *ty = qv_y.data();
    Eigen::Vector3f tmp;
    for (unsigned i = 0; i < outPts.size(); i++) {
        tmp = P * inPts[i];

        tx[2*i + 1] = double(tmp(0)/tmp(2));
        ty[2*i + 1] = double(tmp(1)/tmp(2));
    }

    plotPolygon();
    plot();

    //restore
    for (unsigned i = 0; i < outPts.size(); i++) {
        tx[2*i + 1] = double(outPts[i](0));
        ty[2*i + 1] = double(outPts[i](1));
    }
}

void MainWindow::on_btNewPoint_clicked()
{
    if (ui->ptOut->toPlainText().length() == 0) {
        ui->lbNewPoint2->setText("Matrix P hasn't been calculated");
        return;
    }

    //every point have their image point
    //only last point can be project with P
    if (qv_x.size() % 2 == 0) {
        return;
    }

    //take last point
    double lastX = qv_x.last();
    double lastY = qv_y.last();

    Eigen::Vector3f newPoint = P * Eigen::Vector3f(float(lastX), float(lastY), 1);
    double newX = double(newPoint(0)/newPoint(2));
    double newY = double(newPoint(1)/newPoint(2));

    //new point = P * old point
    ui->lbNewPoint2->setText("New point: " + QString::number(newX) + ", " + QString::number(newY));

    //add to q
    qv_x.push_back(newX);
    qv_y.push_back(newY);

    plot();

    //remove last two points from input
    qv_x.pop_back();
    qv_y.pop_back();
    qv_x.pop_back();
    qv_y.pop_back();
}
