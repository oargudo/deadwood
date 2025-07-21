#include "EventChart.h"
#include <QtCharts/QValueAxis>
#include <QtCharts/qlineseries.h>
#include <QtWidgets/qgraphicsscene.h>
#include <QtCore/QPointer>
#include <iostream>

void EventChart::initAxes(const QString& yaxisTitle)
{
    axisX = new QValueAxis();
    axisY = new QValueAxis();
    axisX->setRange(0, 1);
    axisY->setRange(0, 1);
    this->addAxis(axisX, Qt::AlignBottom);
    this->addAxis(axisY, Qt::AlignLeft);
    this->legend()->hide();
    if (yaxisTitle.length() > 0) {
        axisY->setTitleText(yaxisTitle);
    }
}


void EventChart::addEvent(int year, const QString& str)
{
    // Vertical dashed line
    QLineSeries* lineSeries = new QLineSeries();
    lineSeries->append(year, axisY->min());
    lineSeries->append(year, axisY->max());

    QPen pen(Qt::DashLine);
    pen.setColor(Qt::red);
    lineSeries->setPen(pen);

    this->addSeries(lineSeries);
    //this->addAxis(this->axes(Qt::Vertical).first(), Qt::AlignBottom);
    //this->addAxis(this->axes(Qt::Horizontal).first(), Qt::AlignLeft);
    lineSeries->attachAxis(axisX);
    lineSeries->attachAxis(axisY);

    // Text
    QGraphicsTextItem* textItem = new QGraphicsTextItem(str);
    QGraphicsScene* scene = this->scene();
    if (scene != nullptr)
      scene->addItem(textItem);

    QFont font = textItem->font();
    font.setPointSize(8);
    textItem->setFont(font);
    textItem->setDefaultTextColor(Qt::blue);

    QPointer<QGraphicsTextItem> safeTextItem = textItem;
    QPointer<QChart> safeChart = this;

    auto updateTextPosition = [=]() {
      if (!safeChart || !safeTextItem) return;
        qreal maxX = axisX->max();
        safeTextItem->setPos(safeChart->plotArea().left() + safeChart->plotArea().width() * year / maxX,
                             safeChart->plotArea().top() - safeTextItem->boundingRect().height());
        safeTextItem->setZValue(safeChart->zValue() + 1);
        };

    // Connect signals to update text position when: range changed or chart "resized"

    // Save connections to disconect later if textItems are deleted.
    QList<QMetaObject::Connection> connections;
    connections.push_back(QObject::connect(axisY, &QValueAxis::rangeChanged, [=](qreal min, qreal max) {
        lineSeries->replace(1, QPointF(year, max)); // update vertical line "max" value
        updateTextPosition();
        }));
    connections.push_back(QObject::connect(axisX, &QValueAxis::rangeChanged, updateTextPosition));
    connections.push_back(QObject::connect(this, &QChart::geometryChanged, updateTextPosition));

    events.push_back(textItem);
    textItemConnections.insert(textItem, connections);

    // Init text pos
    updateTextPosition();
}

void EventChart::reset()
{
    if (this->scene() != nullptr) {
        for (QGraphicsTextItem* item : events) {
            this->scene()->removeItem(item);
        }
    }
    events.clear();

    this->removeAllSeries();

    for (QAbstractAxis* ax : this->axes()) {
        this->removeAxis(ax);
        delete ax;
    }
}
