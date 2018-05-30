#ifndef VIZ_H
#define VIZ_H

#include <QObject>

class Viz : public QObject
{
    Q_OBJECT

public:
    // objects
    Q_PROPERTY(QString traj_name
               READ traj_name
               WRITE setTrajName
               NOTIFY trajNameChanged)

    explicit Viz(QObject *parent = 0);

    // functions
    Q_INVOKABLE void openVMD(const QString& data, const QString& data2);


    // local objects
    QString traj_name() { return mTrajName; };


    // setting functions
public slots:
    void setTrajName(const QString& traj_name) { mTrajName = traj_name; };


    // change signals
signals:
    void trajNameChanged(const QString& traj_name);
    void error(const QString& msg);

private:


public:
    QString mTrajName;

};






#endif // VIZ_H
