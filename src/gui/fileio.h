#ifndef FILEIO_H
#define FILEIO_H

#include <QObject>

class FileIO : public QObject
{
    Q_OBJECT

public:
    // objects
    Q_PROPERTY(QString source
               READ source
               WRITE setSource
               NOTIFY sourceChanged)
    Q_PROPERTY(QString grSource
               READ grSource
               WRITE setgrSource
               NOTIFY grSourceChanged);
    Q_PROPERTY(int linecount
               READ linecount
               WRITE setLinecount
               NOTIFY linecountChanged)
    Q_PROPERTY(QString homeDir
               READ homeDir
               WRITE setHomeDir
               NOTIFY homeDirChanged)
    explicit FileIO(QObject *parent = 0);

    // functions
    Q_INVOKABLE QString read();
    Q_INVOKABLE QString read_gr();
    Q_INVOKABLE bool write(const QString& data);

    // local objects
    QString source() { return mSource; };
    QString grSource() { return mgrSource; };
    int linecount() { return mLinecount; };
    QString homeDir() { return mHomedir; };

    // setting functions
public slots:
    void setSource(const QString& source) { mSource = source; };
    void setLinecount(int& linecount) { mLinecount = linecount; };
    void setHomeDir(const QString& homeDir) { mHomedir = homeDir; };
    void setgrSource(const QString& grSource) { mgrSource = grSource; };

    // change signals
signals:
    void sourceChanged(const QString& source);
    void grSourceChanged(const QString& grSource);
    void linecountChanged(int& linecount);
    void homeDirChanged(const QString& homeDir);
    void error(const QString& msg);

private:
    QString mSource;
    QString mgrSource;

public:
    int mLinecount=0;
    QString mHomedir;
};

#endif // FILEIO_H
