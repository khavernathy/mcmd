#include "fileio.h"
#include <QFile>
#include <QTextStream>
#include <QDebug>
#include "unistd.h"

#include <iostream>
#include <string>
#include <fstream>
using namespace std;

FileIO::FileIO(QObject *parent) :
    QObject(parent)
{

}

QString FileIO::read()
{
    // dynamically find the home directory (based on Linux or Mac)
    QString homeDir;
    string linuxcheck="/proc/cpuinfo";
        //linux
        if (std::ifstream(linuxcheck.c_str())) {
            string usfspace = "/home/dfranz";
            string homebox = "/home/khavernathy";
            if (std::ifstream(usfspace.c_str())) {
                homeDir = "/home/dfranz";
            } else {
                homeDir = "/home/khavernathy";
            }
        } else {
            // mac
            homeDir = "/Users/douglasfranz";
        }

        mSource = homeDir + "/mcmd/testzone/runlog.log";

    if (mSource.isEmpty()){
        emit error("source is empty");
        return QString();
    }

    QFile file(mSource);
    QString fileContent;
    if ( file.open(QIODevice::ReadOnly) ) {
        QString line;
        QTextStream t( &file );
        int recounter=0;
        int limit=mLinecount;
        bool write = false;
        bool firstCheck = false;
        if (recounter == 0 && limit == 0) firstCheck = true;
        //printf("firstCheck is %i and mLinecount is currently %i\n",firstCheck,mLinecount);
        do {
            line = t.readLine();
            //if (line.indexOf("MCMD") != -1) printf("recounter %i; lim %i; %s\n", recounter, limit, line.toLatin1().data());
            if (firstCheck) {
                write = true;
            } else {
                if (limit < recounter+1) {
                    write = true;
                }
            }
            if (write) {

                QString spacer="";
                if (recounter<9) spacer =      "      | ";
                else if (recounter<99) spacer = "     | ";
                else if (recounter<999) spacer = "    | ";
                else if (recounter<9999) spacer = "   | ";
                else if (recounter<99999) spacer = "  | ";
                else if (recounter<999999) spacer = " | ";
                else if (recounter<9999999) spacer = "| ";

                // make sure it's a real line.
                if (!line.isNull())
                    fileContent += "\n"+QString::number(recounter+1)+spacer+line;

                //fileContent += "\n"+line;
            }
            recounter++;

        } while (!line.isNull());
        //printf(" ---> TOTAL LINES at %i\n",recounter);

        mLinecount = recounter-1;

        file.close();
    } else {
        emit error("Unable to open the file");
        return QString();
    }
    return fileContent;
}

QString FileIO::read_gr() {

    mgrSource = "/home/khavernathy/mcmd/src/build-gui_mcmd-Desktop_Qt_5_9_1_GCC_64bit-Debug/radial_distribution.dat0";
    //mgrSource = mworkingDir + "/radial_distribution.dat0";
    //mgrSource = qApp->applicationDirPath() + "/radial_distribution.dat0";
    if (mgrSource.isEmpty()){
        emit error("source is empty");
        return QString();
    }

    QFile file(mgrSource);
    QString fileContent;
    if ( file.open(QIODevice::ReadOnly) ) {
        QString line;
        QTextStream t( &file );
        do {
            line = t.readLine();
            fileContent += line+"\n";
        } while (!line.isNull());
        file.close();
    } else {
        emit error("Unable to open the file");
        return QString();
    }
    return fileContent;
}

QString FileIO::read_other(const QString& data) {
    motherSource = data;
    motherSource.remove(0,7); // truncate "file://"

    if (motherSource.isEmpty()){
        emit error("source is empty");
        return QString();
    }

    QFile file(motherSource);
    QString fileContent;
    if ( file.open(QIODevice::ReadOnly) ) {
        QString line;
        QTextStream t( &file );
        do {
            line = t.readLine();
            fileContent += line+"\n";
        } while (!line.isNull());
        file.close();
    } else {
        emit error("Unable to open the file"+motherSource);
        return QString();
    }
    return fileContent;
}


bool FileIO::write(const QString& data)
{
    if (motherSource.isEmpty())
        return false;

    QFile file(motherSource);
    if (!file.open(QFile::WriteOnly | QFile::Truncate))
        return false;

    QTextStream out(&file);
    out << data;

    file.close();

    return true;
}

int FileIO::startSimulation(const QString& data)
{
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
    string filename = data.toStdString();
    string placeholder = filename;

    // determine working directory at location of input file.
    string delimeter = "/";
    string token;
    string workingDir; // = delimeter; // start at root
    size_t pos = 0;
    while ((pos = placeholder.find(delimeter)) != std::string::npos) {
        token = placeholder.substr(0, pos);
        //std::cout << token << endl;
        workingDir = workingDir + token + delimeter; // build the working dir piece by piece
        placeholder.erase(0, pos + delimeter.length());
    }
    //std::cout << "final token: " << token << endl;
    std::cout << "working dir is " << workingDir << endl;
    mworkingDir = QString::fromStdString(workingDir);

    // write it to a runlog in working dir.
    string cdString = "cd "+workingDir;
    system(cdString.c_str());
    string commandString = homeDir+"/mcmd/mcmd "+filename+" | tee "+workingDir+"runlog.log >/dev/null &";
    // RUN MCMD
    system(commandString.c_str());
    return 0;
}

int FileIO::killSimulation(const QString& data)
{
    string filename = data.toStdString();
    string commandString = "kill -9 $(ps aux | grep '"+filename+"' | awk {'print $2'})";
    system(commandString.c_str());
    qDebug() << "Tried to kill simulation via: " << QString::fromStdString(commandString);
    return 0;
}
