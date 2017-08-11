#include "backend.h"

BackEnd::BackEnd(QObject *parent) :
    QObject(parent)
{
}

QString BackEnd::userName()
{
    return m_userName;
}

void BackEnd::setUserName(const QString &userName)
{
    if (userName == m_userName)
        return;

    m_userName = userName;
    emit userNameChanged();
}

int BackEnd::outputLineNumber()
{
    return m_outputLineNumber;
}

void BackEnd::setOutputLineNumber(const int &outputLineNumber)
{
    if (outputLineNumber == m_outputLineNumber)
        return;

    m_outputLineNumber = outputLineNumber;
    emit outputLineNumberChanged();
}

QList<qreal> BackEnd::Qst()
{
    return m_Qst;
}

void BackEnd::setQst(QList<qreal> &Qst)
{
    if (Qst == m_Qst)
        return;

    m_Qst = Qst;
    emit QstChanged();
}
