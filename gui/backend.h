#ifndef BACKEND_H
#define BACKEND_H

#include <QObject>
#include <QString>
#include <QList>

class BackEnd : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QString userName READ userName WRITE setUserName NOTIFY userNameChanged)
    Q_PROPERTY(int outputLineNumber READ outputLineNumber WRITE setOutputLineNumber NOTIFY outputLineNumberChanged)
    Q_PROPERTY(QList<qreal> Qst READ Qst WRITE setQst NOTIFY QstChanged)

public:
    explicit BackEnd(QObject *parent = nullptr);

    QString userName();
    int outputLineNumber();
    QList<qreal> Qst();

    void setUserName(const QString &userName);
    void setOutputLineNumber(const int &outputLineNumber);
    void setQst(QList<qreal> &Qst);

signals:
    void userNameChanged();
    void outputLineNumberChanged();
    void QstChanged();

private:
    QString m_userName;
    int m_outputLineNumber;


public:
    QList<qreal> m_Qst;
};

#endif // BACKEND_H
