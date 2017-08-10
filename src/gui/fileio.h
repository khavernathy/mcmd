#ifndef FILEIO_H
#define FILEIO_H

#include <QObject>

class FileIO : public QObject
{
    Q_OBJECT

public:
    Q_PROPERTY(QString source
               READ source
               WRITE setSource
               NOTIFY sourceChanged)
    Q_PROPERTY(int linecount
               READ linecount
               WRITE setLinecount
               NOTIFY linecountChanged)
    explicit FileIO(QObject *parent = 0);

    Q_INVOKABLE QString read();
    Q_INVOKABLE bool write(const QString& data);

    QString source() { return mSource; };
    int linecount() { return mLinecount; };

public slots:
    void setSource(const QString& source) { mSource = source; };
    void setLinecount(int& linecount) { mLinecount = linecount; };

signals:
    void sourceChanged(const QString& source);
    void linecountChanged(int& linecount);
    void error(const QString& msg);

private:
    QString mSource;

public:
    int mLinecount=0;
};

#endif // FILEIO_H
