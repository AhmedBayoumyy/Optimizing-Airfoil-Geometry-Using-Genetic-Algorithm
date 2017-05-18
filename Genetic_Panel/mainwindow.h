#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include<Eigen/Dense>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <Eigen/Dense>
#include <iterator>
#include <vector>
#include <algorithm>
#include <QtDebug>
#include <ctime>
using namespace std;



namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    void Shape_Function(std::vector<double>a,std::vector<double>b,double X,double &Z_u0,double &Z_d0);
    void Shape_Function_PARSIK(std::vector<double>a,double X,double &Z_u0,double &Z_d0);
    void Parsec(std::vector<double> p, std::vector<double> &a);
    void Panel_Method_Solver(std::vector<double> P0, double AOA, double uinf, int NPanel, double &Cl, double &Cdp, double &max_Thickness);
    void Genetic(int genNo, std::vector<double> P0, std::vector<double>Range, double uinf, double AOA, double Npanels, double &Cloriginal, double & Cl_fittest,std::vector<double>(&fittest_fittest));
    void RandP(std::vector<double>P, std::vector<double>Range, std::vector<double> &P1);

private slots:
    void on_pushButton_clicked();



private:
    Ui::MainWindow *ui;

};

#endif // MAINWINDOW_H
