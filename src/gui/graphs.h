#ifndef GRAPHS_H
#define GRAPHS_H

#include <QWidget>

namespace Ui {
class Graphs;
}

class Graphs : public QWidget
{
    Q_OBJECT

public:
    explicit Graphs(QWidget *parent = 0);
    ~Graphs();

    void addPoint(double x, double y);
    void clearData();
    void plot();

private slots:
    void on_btn_add_clicked();

    void on_btn_clear_clicked();

private:
    Ui::Graphs *ui;

    QVector<double> qv_x, qv_y;
};

#endif // GRAPHS_H
