#pragma once

#include <QtCore/qobject.h>
#include "sim.h"

class SimWorker : public QObject
{
	Q_OBJECT

public:
	SimWorker(ecosim::Simulation* sim);
	~SimWorker();
public slots:
	void simulationStep();
signals:
	void _finishedSimStep();
private:
	ecosim::Simulation* sim;

	float dtime = 0.0f; // time required for death and decay

};

