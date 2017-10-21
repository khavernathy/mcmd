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

void Viz::openVMD(const QString& traj_name, const QString& exe_path)
{

    // try to open VMD
    string command = "/usr/local/bin/vmd";
    command += " -pdb " + exe_path.toStdString() + "/" + traj_name.toStdString();
    system(command.c_str());
    qDebug() << "running openVMD() with CLI command = " + QString::fromStdString(command);

    /* other way
    QProcess *process = new QProcess();
    QString file = "/usr/local/bin/vmd -pdb" + traj_name;
    process->start(file);
    */
}
