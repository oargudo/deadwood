#pragma once

#include <QtCharts/QChart>
#include <QtCharts/QValueAxis>

class EventChart : public QChart
{
public:

		void initAxes(const QString& yaxisTitle = "");
	  void addEvent(int year, const QString& str);
	  void reset();

		QValueAxis* getAxisX() const { return axisX; }
		QValueAxis* getAxisY() const { return axisY; }

		QList<QGraphicsTextItem*> events;

protected:
	  QMap<QGraphicsTextItem*, QList<QMetaObject::Connection>> textItemConnections;

		QValueAxis* axisX = nullptr;
		QValueAxis* axisY = nullptr;
};
