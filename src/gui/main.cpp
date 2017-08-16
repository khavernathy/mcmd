#include <QApplication>
//#include <QGuiApplication>
#include <QObject>
#include <QQmlApplicationEngine>
#include <QDebug>
#include <QStandardPaths>
//#include <QFileDialog>

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
    // RUN MCMD
    //system(commandString.c_str());
    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);

    //QGuiApplication app(argc, argv);
    QApplication app(argc, argv);

    // c++ classes and corresponding headers for C++ <--> QML comm.
    qmlRegisterType<BackEnd>("io.qt.examples.backend", 1, 0, "BackEnd");
    qmlRegisterType<FileIO, 1>("FileIO", 1, 0, "FileIO");
    //qmlRegisterType<Graphs, 1>("Graphs", 1, 0, "Graphs" );

    QQmlApplicationEngine engine;



    // go
    engine.load(QUrl(QLatin1String("qrc:/main.qml")));
    if (engine.rootObjects().isEmpty())
        return -1;
    return app.exec();
}
