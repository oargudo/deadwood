#include "simWorker.h"
#include <Windows.h>
#include <iostream>

SimWorker::SimWorker(ecosim::Simulation* sim)
{
	this->sim = sim;
}

SimWorker::~SimWorker()
{
}

void SimWorker::simulationStep() {
	sim->simStep(dtime);
	emit _finishedSimStep();
}