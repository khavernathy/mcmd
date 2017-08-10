#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <QDebug>

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
    // dev/null is important to supress MCMD output in Qt Application Output
    system("$HOME/mcmd/mcmd $HOME/mcmd/testzone/mcmd.inp | tee $HOME/mcmd/testzone/runlog.log  >/dev/null & ");
    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);

    QGuiApplication app(argc, argv);

    // c++ classes and corresponding headers for C++ <--> QML comm.
    qmlRegisterType<BackEnd>("io.qt.examples.backend", 1, 0, "BackEnd");
    qmlRegisterType<FileIO, 1>("FileIO", 1, 0, "FileIO");

    QQmlApplicationEngine engine;
    engine.load(QUrl(QLatin1String("qrc:/main.qml")));
    if (engine.rootObjects().isEmpty())
        return -1;

    return app.exec();
}
