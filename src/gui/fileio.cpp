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
        printf("firstCheck is %i and mLinecount is currently %i\n",firstCheck,mLinecount);
        do {
            line = t.readLine();

            if (firstCheck) {
                write = true;
            } else {
                if (limit < recounter) {
                    write = true;

                }
            }
            if (write) {
                /*
                QString spacer="";
                if (recounter<9) spacer =      "      | ";
                else if (recounter<99) spacer = "     | ";
                else if (recounter<999) spacer = "    | ";
                else if (recounter<9999) spacer = "   | ";
                else if (recounter<99999) spacer = "  | ";
                else if (recounter<999999) spacer = " | ";
                else if (recounter<9999999) spacer = "| ";
                fileContent += "\n"+QString::number(recounter+1)+spacer+line;
                */
                fileContent += "\n"+line;
            }
            recounter++;

        } while (!line.isNull());
        //printf("HEYYY %i",count);

        mLinecount = recounter;

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
