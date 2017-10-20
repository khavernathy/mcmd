#include "viz.h"
#include <QFile>
#include <QTextStream>
#include <QDebug>
#include <QProcess>
#include "unistd.h"

#include <iostream>
#include <string>
#include <fstream>
using namespace std;

Viz::Viz(QObject *parent) :
    QObject(parent)
{

}

void Viz::openVMD(const QString& traj_name)
{

    // try to open VMD
    QProcess *process = new QProcess();
    QString file = "/usr/local/bin/vmd";
    process->start(file);
    qDebug() << "running openVMD() with traj_name = " + traj_name;
}
