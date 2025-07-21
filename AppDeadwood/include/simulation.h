#pragma once

#include <QtCore/qobject.h>

class Simulation : public QObject
{
	Q_OBJECT

public:
	Simulation();
	~Simulation();
public slots:
	void simulationStep();
signals:
	void _finishedSimStep();
};

