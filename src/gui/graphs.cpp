#include "graphs.h"
#include "ui_graphs.h"

Graphs::Graphs(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Graphs)
{
    ui->setupUi(this);

    ui->plot->addGraph();
    ui->plot->graph(0)->setScatterStyle(QCPScatterStyle::ssCircle); // circles on scatter plot
    ui->plot->graph(0)->setLineStyle(QCPGraph::lsNone); // no line
}

Graphs::~Graphs()
{
    delete ui;
}

void Graphs::addPoint(double x, double y)
{
    qv_x.append(x);
    qv_y.append(y);
}

void Graphs::clearData()
{

}

void Graphs::plot()
{
    ui->plot->graph(0)->setData(qv_x, qv_y);
    ui->plot->replot();
    ui->plot->update();
}

void Graphs::on_btn_add_clicked()
{
    addPoint(ui->bx_x->value(), ui->bx_y->value());
    plot();
}

void Graphs::on_btn_clear_clicked()
{

}
