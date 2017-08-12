#include "fileio.h"
#include <QFile>
#include <QTextStream>
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

bool FileIO::write(const QString& data)
{
    if (mSource.isEmpty())
        return false;

    QFile file(mSource);
    if (!file.open(QFile::WriteOnly | QFile::Truncate))
        return false;

    QTextStream out(&file);
    out << data;

    file.close();

    return true;
}
