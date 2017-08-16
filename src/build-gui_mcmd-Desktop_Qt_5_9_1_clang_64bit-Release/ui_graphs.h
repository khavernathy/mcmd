/********************************************************************************
** Form generated from reading UI file 'graphs.ui'
**
** Created by: Qt User Interface Compiler version 5.9.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_GRAPHS_H
#define UI_GRAPHS_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QWidget>
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_Graphs
{
public:
    QCustomPlot *plot;
    QDoubleSpinBox *bx_x;
    QDoubleSpinBox *bx_y;
    QPushButton *btn_add;
    QPushButton *btn_clear;
    QLabel *label;
    QLabel *label_2;

    void setupUi(QWidget *Graphs)
    {
        if (Graphs->objectName().isEmpty())
            Graphs->setObjectName(QStringLiteral("Graphs"));
        Graphs->resize(666, 659);
        plot = new QCustomPlot(Graphs);
        plot->setObjectName(QStringLiteral("plot"));
        plot->setGeometry(QRect(40, 10, 581, 381));
        bx_x = new QDoubleSpinBox(Graphs);
        bx_x->setObjectName(QStringLiteral("bx_x"));
        bx_x->setGeometry(QRect(140, 470, 69, 26));
        bx_y = new QDoubleSpinBox(Graphs);
        bx_y->setObjectName(QStringLiteral("bx_y"));
        bx_y->setGeometry(QRect(350, 480, 69, 26));
        btn_add = new QPushButton(Graphs);
        btn_add->setObjectName(QStringLiteral("btn_add"));
        btn_add->setGeometry(QRect(440, 480, 89, 25));
        btn_clear = new QPushButton(Graphs);
        btn_clear->setObjectName(QStringLiteral("btn_clear"));
        btn_clear->setGeometry(QRect(540, 480, 89, 25));
        label = new QLabel(Graphs);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(60, 480, 67, 17));
        label_2 = new QLabel(Graphs);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setGeometry(QRect(260, 480, 67, 17));

        retranslateUi(Graphs);

        QMetaObject::connectSlotsByName(Graphs);
    } // setupUi

    void retranslateUi(QWidget *Graphs)
    {
        Graphs->setWindowTitle(QApplication::translate("Graphs", "Form", Q_NULLPTR));
        btn_add->setText(QApplication::translate("Graphs", "Add", Q_NULLPTR));
        btn_clear->setText(QApplication::translate("Graphs", "Clear", Q_NULLPTR));
        label->setText(QApplication::translate("Graphs", "TextLabel", Q_NULLPTR));
        label_2->setText(QApplication::translate("Graphs", "TextLabel", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class Graphs: public Ui_Graphs {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_GRAPHS_H
