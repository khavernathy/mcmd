#include <QGuiApplication>
//#include <QApplication>
#include <QObject>
#include <QQmlApplicationEngine>
#include <QDebug>
#include <QtCharts/QChartView>
#include <QtCharts/QPieSeries>
#include <QtCharts/QPieSlice>
#include <QtCharts/QAbstractBarSeries>
#include <QtCharts/QPercentBarSeries>
#include <QtCharts/QStackedBarSeries>
#include <QtCharts/QBarSeries>
#include <QtCharts/QBarSet>
#include <QtCharts/QLineSeries>
#include <QtCharts/QSplineSeries>
#include <QtCharts/QScatterSeries>
#include <QtCharts/QAreaSeries>
#include <QtCharts/QLegend>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QFormLayout>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QLabel>
#include <QtCore/QTime>
#include <QtCharts/QBarCategoryAxis>



#include "backend.h"
#include "fileio.h"

#include <unistd.h>
#include <iostream>
#include <string>
#include <fstream>



using namespace std;



int main(int argc, char *argv[])
{
    system("pwd");

    string homeDir;
    string linuxcheck="/proc/cpuinfo";
        //linux
        if (std::ifstream(linuxcheck.c_str())) {
            string usfspace = "/home/dfranz";
            string homebox = "/home/khavernathy";
            if (std::ifstream(usfspace.c_str())) {
                homeDir = usfspace;
            } else {
                homeDir = homebox;
            }
        } else {
            // mac
            homeDir = "/Users/douglasfranz";
        }


    //string commandString = homeDir+"/mcmd/mcmd "+homeDir+"/mcmd/testzone/mcmd.inp | tee "+homeDir+"/mcmd/testzone/runlog.log &";
    // dev/null is important to supress MCMD output in Qt Application Output
    string commandString = homeDir+"/mcmd/mcmd "+homeDir+"/mcmd/testzone/mcmd.inp | tee "+homeDir+"/mcmd/testzone/runlog.log >/dev/null &";
    system(commandString.c_str());
    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);

    QGuiApplication app(argc, argv);
    //QApplication app(argc, argv);

    // c++ classes and corresponding headers for C++ <--> QML comm.
    qmlRegisterType<BackEnd>("io.qt.examples.backend", 1, 0, "BackEnd");
    qmlRegisterType<FileIO, 1>("FileIO", 1, 0, "FileIO");

    QQmlApplicationEngine engine;
    engine.load(QUrl(QLatin1String("qrc:/main.qml")));
    if (engine.rootObjects().isEmpty())
        return -1;

    return app.exec();
}
