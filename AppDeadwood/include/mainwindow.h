#ifndef __MainWindow__
#define __MainWindow__

#include <QtWidgets/QMainWindow>
#include "ui_main.h"
#include "ecosimwidget.h"
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCore/QPropertyAnimation>
#include "EventChart.h"

class MainWindow : public QMainWindow
{
	Q_OBJECT
private:
	Ui::MainWindow ui; //!< Interface : QtDesigner.
	EcosimWidget* ecosimWidget;

	QString sceneDir;

	// chart
	int chartLastYear;
	std::vector<QLineSeries*> speciesCntSeries, speciesAreaSeries, typeCntSeries, typeAreaSeries;
	EventChart speciesCntChart, speciesAreaChart, typeCntChart, typeAreaChart;
	float maxSpeciesCnt, maxSpeciesArea, maxTypeCnt, maxTypeArea;

	std::vector<ecosim::ScriptedDisturbanceEvent*> eventlist;
	int editingEventIndex;

	QPropertyAnimation* saveEventListAnim;

public:
	MainWindow();
	~MainWindow();

private:
	void createActions();

public slots:
	void sceneLoad();
	void exportSim();
	void importSim();

	void simRun();
	void simStop();
	void simStep();
	void simReset();

	void changeTexture();

	void setInfoText(const QString& s);
	void setPickedPoint(const QPoint& p);

	void resetEventsUI(bool unselectEvent);
	void loadEvent();
	void setEvent(ecosim::DisturbanceEventCategory);
	void removeEvent(int index);
	void saveEventList();
	void loadEventList();

	void updateFilters();

	void changeChart(EventChart* c);

protected:
	void fillSpeciesCombobox();
	void createSpeciesCheckbox();

	void initCharts();
	void updateCharts();
	void resetCharts();

	void initEventsList();
	void updateEventsList();
	//void addVerticalLineAndText(QChart* chart, qreal xPosition, const QString& text);

};

#endif