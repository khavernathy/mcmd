/********************************************************************************
** Form generated from reading UI file 'page1.ui'
**
** Created by: Qt User Interface Compiler version 5.9.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_PAGE1_H
#define UI_PAGE1_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QHeaderView>

QT_BEGIN_NAMESPACE

class Ui_Page1
{
public:
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *Page1)
    {
        if (Page1->objectName().isEmpty())
            Page1->setObjectName(QStringLiteral("Page1"));
        Page1->resize(613, 408);
        buttonBox = new QDialogButtonBox(Page1);
        buttonBox->setObjectName(QStringLiteral("buttonBox"));
        buttonBox->setGeometry(QRect(30, 240, 341, 32));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        retranslateUi(Page1);
        QObject::connect(buttonBox, SIGNAL(accepted()), Page1, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), Page1, SLOT(reject()));

        QMetaObject::connectSlotsByName(Page1);
    } // setupUi

    void retranslateUi(QDialog *Page1)
    {
        Page1->setWindowTitle(QApplication::translate("Page1", "Dialog", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class Page1: public Ui_Page1 {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PAGE1_H
