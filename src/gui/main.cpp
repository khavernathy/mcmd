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
    printf("Result of pwd: ");
    system("pwd");

    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);

    //QGuiApplication app(argc, argv);
    QApplication app(argc, argv);
    qDebug() << "Result of qApp->applicationDirPath(): " << qApp->applicationDirPath();


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
