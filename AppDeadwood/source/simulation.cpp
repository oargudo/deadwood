#include "simulation.h"
#include <Windows.h>
#include <iostream>

Simulation::Simulation()
{
}

Simulation::~Simulation()
{
}

void Simulation::simulationStep() {
	std::cout << "Simulation work..." << std::endl;
	Sleep(1000);
	emit _finishedSimStep();
}