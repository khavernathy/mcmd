#include <QGuiApplication>
#include <QQmlApplicationEngine>
using namespace std;

/*
void myMessageOutput(QtMsgType type, const QMessageLogContext &context, const QString &msg)
{
    QByteArray localMsg = msg.toLocal8Bit();
    switch (type) {
    case QtDebugMsg:
        fprintf(stderr, "Debug: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        break;
    case QtInfoMsg:
        fprintf(stderr, "Info: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        break;
    case QtWarningMsg:
        fprintf(stderr, "Warning: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        break;
    case QtCriticalMsg:
        fprintf(stderr, "Critical: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        break;
    case QtFatalMsg:
        fprintf(stderr, "Fatal: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        abort();
    }
}
*/

int main(int argc, char *argv[])
{
    //std::system("cd /home/dfranz/mcmd/examples/md/nve");
    std::system("pwd");
    std::system("/home/dfranz/mcmd/mcmd /home/dfranz/mcmd/testzone/mcmd.inp /home/dfranz/mcmd/testzone/input.pdb | tee runlog.log &");
    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
    //qInstallMessageHandler(myMessageOutput); // Install the handler


    QGuiApplication app(argc, argv);

    QQmlApplicationEngine engine;
    engine.load(QUrl(QLatin1String("qrc:/main.qml")));
    if (engine.rootObjects().isEmpty())
        return -1;

    return app.exec();
}
