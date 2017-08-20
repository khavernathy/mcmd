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
               NOTIFY grSourceChanged)
    Q_PROPERTY(QString inputSource
               READ inputSource
               WRITE setinputSource
               NOTIFY inputSourceChanged)
    Q_PROPERTY(int linecount
               READ linecount
               WRITE setLinecount
               NOTIFY linecountChanged)
    Q_PROPERTY(QString homeDir
               READ homeDir
               WRITE setHomeDir
               NOTIFY homeDirChanged)
    Q_PROPERTY(QString workingDir
               READ workingDir
               WRITE setWorkingDir
               NOTIFY workingDirChanged)
    Q_PROPERTY(QString mode
               READ mode
               WRITE setMode
               NOTIFY modeChanged)
    explicit FileIO(QObject *parent = 0);

    // functions
    Q_INVOKABLE QString read(const QString& data);
    Q_INVOKABLE QString read_gr(const QString& data);
    Q_INVOKABLE QString read_input(const QString& data);
    Q_INVOKABLE int startSimulation(const QString& data, const QString& data2);
    Q_INVOKABLE int killSimulation(const QString& data);
    Q_INVOKABLE bool write(const QString& data);

    // local objects
    QString source() { return mSource; };
    QString grSource() { return mgrSource; };
    QString inputSource() { return minputSource; };
    int linecount() { return mLinecount; };
    QString homeDir() { return mHomedir; };
    QString workingDir() { return mworkingDir; };
    QString mode() { return mmode; };

    // setting functions
public slots:
    void setSource(const QString& source) { mSource = source; };
    void setLinecount(int& linecount) { mLinecount = linecount; };
    void setHomeDir(const QString& homeDir) { mHomedir = homeDir; };
    void setgrSource(const QString& grSource) { mgrSource = grSource; };
    void setinputSource(const QString& inputSource) { minputSource = inputSource; };
    void setWorkingDir(const QString& workingDir) { mworkingDir = workingDir; };
    void setMode(const QString& mode) { mmode = mode; };

    // change signals
signals:
    void sourceChanged(const QString& source);
    void grSourceChanged(const QString& grSource);
    void inputSourceChanged(const QString& inputSource);
    void linecountChanged(int& linecount);
    void homeDirChanged(const QString& homeDir);
    void workingDirChanged(const QString& workingDir);
    void modeChanged(const QString& mode);
    void error(const QString& msg);

private:
    QString mSource;
    QString mgrSource;
    QString minputSource;
    QString mworkingDir;
    QString mmode; // "mc" or "md"

public:
    int mLinecount=0;
    QString mHomedir;
};

#endif // FILEIO_H
