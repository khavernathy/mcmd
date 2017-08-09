#ifndef BACKEND_H
#define BACKEND_H

#include <QObject>
#include <QString>

class BackEnd : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QString userName READ userName WRITE setUserName NOTIFY userNameChanged)
    Q_PROPERTY(int outputLineNumber READ outputLineNumber WRITE setOutputLineNumber NOTIFY outputLineNumberChanged)

public:
    explicit BackEnd(QObject *parent = nullptr);

    QString userName();
    int outputLineNumber();
    void setUserName(const QString &userName);
    void setOutputLineNumber(const int &outputLineNumber);

signals:
    void userNameChanged();
    void outputLineNumberChanged();

private:
    QString m_userName;
    int m_outputLineNumber;
};

#endif // BACKEND_H
