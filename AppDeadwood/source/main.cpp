#include <QtWidgets/QApplication>
#include <QtWidgets/QStyleFactory>
#include <QtCore/QLibraryInfo>
#include "mainwindow.h"

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);

	QStringList styles = QStyleFactory::keys();
	if (styles.contains("windowsvista", Qt::CaseInsensitive)) {
		QApplication::setStyle("windowsvista");
	}
	else if (styles.contains("fusion", Qt::CaseInsensitive)) {
		QApplication::setStyle("fusion");
	}
	else {
		qDebug() << "Styles 'windowsvista' and 'fusion' not found!";
		qDebug() << "Available Qt Styles:";
		for (const QString& style : styles) {
			qDebug() << " " << style.toLower();
		}
		qDebug() << "The application's appearance may differ from the intended design.";
		qDebug() << "Plugins path:" << QLibraryInfo::path(QLibraryInfo::PluginsPath);
	}

	QSurfaceFormat f;
	f.setVersion(4, 5);
	f.setProfile(QSurfaceFormat::CoreProfile);
	QSurfaceFormat::setDefaultFormat(f);

	MainWindow mainWin;
	mainWin.showMaximized();

	return app.exec();
}
