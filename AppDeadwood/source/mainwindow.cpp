#include "mainwindow.h"
#include <QtWidgets/QFileDialog>
#include <QtCore/QCoreApplication>
#include <QtCore/QSettings>
#include <QtCharts/QValueAxis>
#include <QtWidgets/qgraphicseffect.h>
#include <fstream>

/*!
\brief Create the main window.
*/
MainWindow::MainWindow()
{
    // Loading interface
    ui.setupUi(this);
    this->statusBar()->hide();

    // Loading GLWidget
    ecosimWidget = new EcosimWidget();
    //ecosimWidget->SetLight(Vector(0.0, 0.0, 100.0));
    QGridLayout* GLlayout = new QGridLayout;
    GLlayout->addWidget(ecosimWidget, 0, 0);
    GLlayout->setContentsMargins(0, 0, 0, 0);
    ui.widget_GL->setLayout(GLlayout);

    // Connecting QT Actions
    createActions();

    initEventsList();

    editingEventIndex = -1;

    QGraphicsColorizeEffect* effect = new QGraphicsColorizeEffect(ui.groupBox_8);
    effect->setColor(Qt::green);
    effect->setStrength(0);
    ui.groupBox_8->setGraphicsEffect(effect);

    saveEventListAnim = new QPropertyAnimation(effect, "strength");
    saveEventListAnim->setDuration(1000);
    saveEventListAnim->setStartValue(1.0);
    saveEventListAnim->setEndValue(0.0);
}

/*!
\brief Destructor.
*/
MainWindow::~MainWindow()
{
    delete ecosimWidget;
}


/*!
\brief Create callbacks between member slots and user interface.
*/
void MainWindow::createActions()
{
    connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));

    connect(ui.btnSimRun, SIGNAL(clicked()), this, SLOT(simRun()));
    connect(ui.btnSimStop, SIGNAL(clicked()), this, SLOT(simStop()));
    connect(ui.btnSimStep, SIGNAL(clicked()), this, SLOT(simStep()));
    connect(ui.btnSimReset, SIGNAL(clicked()), this, SLOT(simReset()));

    connect(ui.btnSceneLoad, SIGNAL(clicked()), this, SLOT(sceneLoad()));
    connect(ui.btnExportSim, SIGNAL(clicked()), this, SLOT(exportSim()));
    connect(ui.btnImportSim, SIGNAL(clicked()), this, SLOT(importSim()));

    connect(ui.rbTexTerrain, SIGNAL(clicked()), this, SLOT(changeTexture()));
    connect(ui.rbTexMonthly, SIGNAL(clicked()), this, SLOT(changeTexture()));
    connect(ui.rbTexSpecies, SIGNAL(clicked()), this, SLOT(changeTexture()));
    connect(ui.rbTexMostViab, SIGNAL(clicked()), this, SLOT(changeTexture()));
    connect(ui.cbTexTerrain, SIGNAL(currentIndexChanged(int)), this, SLOT(changeTexture()));
    connect(ui.cbTexMonthly, SIGNAL(currentIndexChanged(int)), this, SLOT(changeTexture()));
    connect(ui.cbTexSpecies, SIGNAL(currentIndexChanged(int)), this, SLOT(changeTexture()));
    connect(ui.cbPlantColor, SIGNAL(currentIndexChanged(int)), ecosimWidget, SLOT(setPlantColor(int)));
    connect(ui.checkShading, SIGNAL(toggled(bool)), ecosimWidget, SLOT(setRenderShadeTerrain(bool)));
    connect(ui.checkRenderStats, SIGNAL(toggled(bool)), ecosimWidget, SLOT(setRenderShowStats(bool)));

    connect(ui.checkPlt, &QCheckBox::toggled, ecosimWidget, &EcosimWidget::toggleRenderPlts);
    connect(ui.checkSnag, &QCheckBox::toggled, ecosimWidget, &EcosimWidget::toggleRenderSnags);
    connect(ui.checkLog, &QCheckBox::toggled, ecosimWidget, &EcosimWidget::toggleRenderLogs);

    //filters 
    connect(ui.cbMinHeight, &QCheckBox::toggled, ui.dsbMinHeight, &QDoubleSpinBox::setEnabled);
    connect(ui.cbMaxHeight, &QCheckBox::toggled, ui.dsbMaxHeight, &QDoubleSpinBox::setEnabled);
    connect(ui.cbMinAge, &QCheckBox::toggled, ui.sbMinAge, &QSpinBox::setEnabled);
    connect(ui.cbMaxAge, &QCheckBox::toggled, ui.sbMaxAge, &QSpinBox::setEnabled);
    connect(ui.cbMinHeight, &QCheckBox::toggled, this, &MainWindow::updateFilters);
    connect(ui.cbMaxHeight, &QCheckBox::toggled, this, &MainWindow::updateFilters);
    connect(ui.cbMinAge, &QCheckBox::toggled, this, &MainWindow::updateFilters);
    connect(ui.cbMaxAge, &QCheckBox::toggled, this, &MainWindow::updateFilters);
    connect(ui.dsbMinHeight, &QDoubleSpinBox::valueChanged, this, &MainWindow::updateFilters);
    connect(ui.dsbMaxHeight, &QDoubleSpinBox::valueChanged, this, &MainWindow::updateFilters);
    connect(ui.sbMinAge, &QSpinBox::valueChanged, this, &MainWindow::updateFilters);
    connect(ui.sbMaxAge, &QSpinBox::valueChanged, this, &MainWindow::updateFilters);

    connect(ecosimWidget, &EcosimWidget::_signalFinishedYear, this, &MainWindow::updateCharts);
    connect(ecosimWidget, &EcosimWidget::_signalTextInfo, this, &MainWindow::setInfoText);
    connect(ecosimWidget, &EcosimWidget::_signalClickedTerrainCell, this, &MainWindow::setPickedPoint);
    connect(ecosimWidget, &EcosimWidget::_signalResetSim, this, [this]() { resetCharts(); initCharts(); });

    connect(ui.dockWidget, &QDockWidget::topLevelChanged, this, [this](bool toplevel) {
        if (toplevel)
            ui.dockWidget->resize(800, 600);
        });

    connect(ui.cbChartSelection, &QComboBox::currentIndexChanged, this, [this](int idx) {
      switch (idx)
      {
      case 0: changeChart(&speciesCntChart); break;
      case 1: changeChart(&speciesAreaChart); break;
      case 2: changeChart(&typeCntChart); break;
      case 3: changeChart(&typeAreaChart); break;
      default: break;
      }
      });

    //connect(ui.tableWidget, &QTableWidget::itemSelectionChanged, this, [this]() { ui.btnEditEvent->setEnabled(not ui.tableWidget->selectedItems().isEmpty());});
    //connect(ui.btnEditEvent, &QPushButton::clicked, this, &MainWindow::loadEvent);
    connect(ui.tableWidget, &QTableWidget::itemSelectionChanged, this, &MainWindow::loadEvent);

    connect(ui.tableWidget, &QTableWidget::itemSelectionChanged, this, [this]() { ui.btnRmEvent->setEnabled(not ui.tableWidget->selectedItems().isEmpty());});

    connect(ui.btnRmEvent, &QPushButton::clicked, this, [this]() { removeEvent(ui.tableWidget->row(ui.tableWidget->currentItem())); updateEventsList(); });
    connect(ui.btnNewEvent, &QPushButton::clicked, this, [this]() { resetEventsUI(true); });
    connect(ui.btnSaveEventList, &QPushButton::clicked, this, &MainWindow::saveEventList);
    connect(ui.btnLoadEventList, &QPushButton::clicked, this, &MainWindow::loadEventList);

    connect(ui.btnSelectPtWf, &QPushButton::clicked, this, [this]() { ui.btnSelectPtDs->setChecked(false); 
                                                                      ecosimWidget->toggleTerrainPicking(ui.btnSelectPtWf->isChecked()); });
    connect(ui.btnSelectPtDs, &QPushButton::clicked, this, [this]() { ui.btnSelectPtWf->setChecked(false);
                                                                      ecosimWidget->toggleTerrainPicking(ui.btnSelectPtDs->isChecked());  });

    // wildfire
    connect(ui.lwOrgPtsWf, &QListWidget::itemSelectionChanged, this, [this]() { ui.btnRmOrgPtWf->setEnabled(not ui.lwOrgPtsWf->selectedItems().isEmpty()); });
    connect(ui.btnRmOrgPtWf, &QPushButton::clicked, this, [this]() { delete ui.lwOrgPtsWf->takeItem(ui.lwOrgPtsWf->row(ui.lwOrgPtsWf->currentItem())); });
    connect(ui.btnAddOrgPtWf, &QPushButton::clicked, this, [this]() {
        QString itemText = QString("(%1, %2) r:%3")
            .arg(ui.dsbWfX->value())
            .arg(ui.dsbWfY->value())
            .arg(ui.dsbWfR->value());
        QListWidgetItem* item = new QListWidgetItem(itemText, ui.lwOrgPtsWf);
        item->setData(Qt::UserRole, QVariant(ui.dsbWfX->value()));
        item->setData(Qt::UserRole + 1, QVariant(ui.dsbWfY->value()));
        item->setData(Qt::UserRole + 2, QVariant(ui.dsbWfR->value()));
        });
    connect(ui.pbSetFireEvent, &QPushButton::clicked, this, [this]() { setEvent(ecosim::DisturbanceEventCategory::FIRE); });

    // windstorm
    connect(ui.pbSetWindEvent, &QPushButton::clicked, this, [this]() { setEvent(ecosim::DisturbanceEventCategory::WIND); });

    // disease
    connect(ui.pbSetDiseaseEvent, &QPushButton::clicked, this, [this]() { setEvent(ecosim::DisturbanceEventCategory::DISEASE); });

    // drought
    connect(ui.lwDrValues, &QListWidget::itemSelectionChanged, this, [this]() { ui.btnDrRmOrgPt->setEnabled(not ui.lwDrValues->selectedItems().isEmpty()); });
    connect(ui.btnDrRmOrgPt, &QPushButton::clicked, this, [this]() { delete ui.lwDrValues->takeItem(ui.lwDrValues->row(ui.lwDrValues->currentItem())); });
    connect(ui.btnDrAddOrgPt, &QPushButton::clicked, this, [this]() {
        QListWidgetItem* item = new QListWidgetItem(QString::number(ui.dsbDrOrgPt->value()), ui.lwDrValues);
        item->setData(Qt::UserRole, QVariant(ui.dsbDrOrgPt->value()));
        });
    connect(ui.pbSetDroughtEvent, &QPushButton::clicked, this, [this]() { setEvent(ecosim::DisturbanceEventCategory::DROUGHT); });

}


void MainWindow::sceneLoad()
{
    // Load the last used directory from settings
    QSettings settings("deadwood", "viewer");
    QString lastDirectory = settings.value("LastDirectory").toString();

    // If lastDirectory is empty, use the current directory
    if (lastDirectory.isEmpty())
        lastDirectory = QDir::currentPath();

    // Open file dialog with the last used directory
    sceneDir = QFileDialog::getExistingDirectory(
        this, "Open scene directory", lastDirectory
    );

    if (!sceneDir.isEmpty()) {

        // Store the selected directory as the last used directory
        settings.setValue("LastDirectory", sceneDir);

        // TODO: replace this loading code for a call to ecosim Terrain::loadElv
        // and then use ecosim Terrain data to create the HeightField
        QString qs = sceneDir + QString("/heightfield.elv");
        std::string filename = qs.toStdString();
        
        std::ifstream infile;
        infile.open(filename, std::ios_base::in);

        if (infile.is_open()) {
            infile.close();

            // update widget
            if (!ecosimWidget->loadScene(sceneDir.toStdString())) return;
            ecosimWidget->resetCamera();
            ecosimWidget->updatePlantInstances();

            // update charts
            resetCharts();
            initCharts();
            editingEventIndex = -1;
        }
        else
        {
            std::cerr << "Error, unable to open file " << filename << std::endl;
        }

        // Load other files (species database, sun/wet/temp maps, etc.)

        fillSpeciesCombobox();
        createSpeciesCheckbox();
        updateEventsList();
        
    }

}

void MainWindow::exportSim()
{
    QString selectedFilter;
    QString qfilename = QFileDialog::getSaveFileName(
        this,
        tr("Export Simulation"),
        QString(),
        "JSON (*.json);;PDB Terse (*.pdb);;PDB Full (*.pdb)",
        &selectedFilter
    );
    if (qfilename.isEmpty()) return;

    bool terse = false;
    if (selectedFilter.contains("JSON")) {
        if (!qfilename.endsWith(".json", Qt::CaseInsensitive))
            qfilename += ".json";
        terse = true; // Not relevant for JSON
    }
    else if (selectedFilter.contains("Terse")) {
        if (!qfilename.endsWith(".pdb", Qt::CaseInsensitive))
            qfilename += ".pdb";
        terse = true;
    }
    else if (selectedFilter.contains("Full")) {
        if (!qfilename.endsWith(".pdb", Qt::CaseInsensitive))
            qfilename += ".pdb";
        terse = false;
    }
    if (!ecosimWidget->exportSim(qfilename.toStdString(), terse)) {
        std::cerr << "ERROR: export failed" << std::endl;
    }
}

void MainWindow::importSim()
{
    // TODO
}


void MainWindow::simRun()
{
    ecosimWidget->runSim();
}

void MainWindow::simStop()
{
    ecosimWidget->pauseSim();
}

void MainWindow::simStep()
{
    ecosimWidget->stepSim(ui.sbSimNumSteps->value());
}

void MainWindow::simReset()
{
    ecosimWidget->resetSim();
    ui.textTree->clear();
    //resetCharts();
    //initCharts();
}


void MainWindow::changeTexture()
{
    if (ui.rbTexTerrain->isChecked()) {
        switch (ui.cbTexTerrain->currentIndex()) {
            case 0: ecosimWidget->setRenderTexture(EcosimWidget::PLAIN); break;
            case 1: ecosimWidget->setRenderTexture(EcosimWidget::ELEVATION); break;
            case 2: ecosimWidget->setRenderTexture(EcosimWidget::NORMALS); break;
            case 3: ecosimWidget->setRenderTexture(EcosimWidget::ASPECT_E); break;
            case 4: ecosimWidget->setRenderTexture(EcosimWidget::ASPECT_N); break;
            case 5: ecosimWidget->setRenderTexture(EcosimWidget::SLOPE); break;
        }
    }
    else if (ui.rbTexMonthly->isChecked()) {
        switch (ui.cbTexMonthly->currentIndex()) {
            case 0: ecosimWidget->setRenderTexture(EcosimWidget::SUNLIGHT); break;
            case 1: ecosimWidget->setRenderTexture(EcosimWidget::MOISTURE); break;
            case 2: ecosimWidget->setRenderTexture(EcosimWidget::TEMPERATURE); break;
        }
    }
    else if (ui.rbTexSpecies->isChecked()) {
        ecosimWidget->setRenderTexture(EcosimWidget::VIABILITY, ui.cbTexSpecies->currentIndex());
    }
    else if (ui.rbTexMostViab->isChecked()) {
        ecosimWidget->setRenderTexture(EcosimWidget::MOSTVIABLE);
    }
}

void MainWindow::initCharts()
{
    int numSpecies = static_cast<int>(ecosimWidget->getSpeciesNames().size());
    
    speciesCntSeries.resize(numSpecies);
    speciesAreaSeries.resize(numSpecies);
    typeCntSeries.resize(3);
    typeAreaSeries.resize(3);

    speciesCntChart.initAxes("#plants / ha");
    speciesAreaChart.initAxes("m2 / ha");

    for (int i = 0; i < speciesCntSeries.size(); ++i) {
        QColor c = ecosimWidget->getColorForSpecies(i);
        speciesCntSeries[i] = new QLineSeries();
        speciesCntSeries[i]->setColor(c);
        speciesCntChart.addSeries(speciesCntSeries[i]);
        speciesCntSeries[i]->attachAxis(speciesCntChart.getAxisX());
        speciesCntSeries[i]->attachAxis(speciesCntChart.getAxisY());

        speciesAreaSeries[i] = new QLineSeries();
        speciesAreaSeries[i]->setColor(c);
        speciesAreaChart.addSeries(speciesAreaSeries[i]);
        speciesAreaSeries[i]->attachAxis(speciesAreaChart.getAxisX());
        speciesAreaSeries[i]->attachAxis(speciesAreaChart.getAxisY());
    }

    typeCntChart.initAxes("#plants / ha");
    typeAreaChart.initAxes("m2 / ha");

    for (int i = 0; i < 3; i++) {
        QColor c;
        switch (i) {
        case 0: c = QColor(232, 0, 0); break;
        case 1: c = QColor(160, 0, 0); break;
        case 2: c = QColor(74, 34, 6); break;
        }

        typeCntSeries[i] = new QLineSeries();
        typeCntSeries[i]->setColor(c);
        typeCntChart.addSeries(typeCntSeries[i]);
        typeCntSeries[i]->attachAxis(typeCntChart.getAxisX());
        typeCntSeries[i]->attachAxis(typeCntChart.getAxisY());

        typeAreaSeries[i] = new QLineSeries();
        typeAreaSeries[i]->setColor(c);
        typeAreaChart.addSeries(typeAreaSeries[i]);
        typeAreaSeries[i]->attachAxis(typeAreaChart.getAxisX());
        typeAreaSeries[i]->attachAxis(typeAreaChart.getAxisY());
    }


    ui.graphicsView->setChart(&speciesCntChart);
    ui.graphicsView->setRenderHint(QPainter::Antialiasing);

    chartLastYear = -1;
    maxSpeciesCnt = 0;
    maxSpeciesArea = 0;
    maxTypeCnt = 0;
    maxTypeArea = 0;
}

void MainWindow::updateCharts()
{
    int currentYear = ecosimWidget->getSimulatedYears(); // int(speciesCntHistory[0].size()) - 1;
    int numPFT = ecosimWidget->getNumSpecies(); //int(speciesCntHistory.size());

    for (int year = chartLastYear + 1; year <= currentYear; ++year) {
        const std::vector<float>& datacounts = ecosimWidget->getDataCount(year, 0);
        const std::vector<float>& databasalarea = ecosimWidget->getDataBasal(year, 0);
        for (int pft = 0; pft < numPFT; ++pft) {
            speciesCntSeries[pft]->append(year, datacounts[pft]);
            if (datacounts[pft] > maxSpeciesCnt) maxSpeciesCnt = datacounts[pft];
            speciesAreaSeries[pft]->append(year, databasalarea[pft]);
            if (databasalarea[pft] > maxSpeciesArea) maxSpeciesArea = databasalarea[pft];
        }

        for (int i = 0; i < 3; i++) {
            const std::vector<float>& typecounts = ecosimWidget->getDataCount(year, i);
            const std::vector<float>& typebasalarea = ecosimWidget->getDataBasal(year, i);

            float sumCount = 0;
            float sumArea = 0;
            for (int pft = 0; pft < numPFT; ++pft) {
                sumCount += typecounts[pft];
                sumArea += typebasalarea[pft];
            }

            typeCntSeries[i]->append(year, sumCount);
            if (sumCount > maxTypeCnt) maxTypeCnt = sumCount;
            typeAreaSeries[i]->append(year, sumArea);
            if (sumArea > maxTypeArea) maxTypeArea = sumArea;
        }
    }

    std::vector<ecosim::DisturbanceEvent> newEvents = ecosimWidget->getPastNewEvents();
    int i = 0;
    while (i < newEvents.size()) {
        int year = newEvents[i].year;
        QString eventsStr;
        while (i < newEvents.size() && year == newEvents[i].year) {
            if (!eventsStr.isEmpty())
                eventsStr += ", ";
            switch (newEvents[i].cat) {
            case ecosim::DisturbanceEventCategory::FIRE:
                eventsStr += "Wildfire";
                break;
            case ecosim::DisturbanceEventCategory::WIND:
                eventsStr += "Windstorm";
                break;
            case ecosim::DisturbanceEventCategory::DISEASE:
                eventsStr += "Disease";
                break;
            case ecosim::DisturbanceEventCategory::DROUGHT:
                eventsStr += "Drought";
                break;
            default:
                eventsStr += "Error";
                break;
            }
            i++;
        }
        speciesCntChart.addEvent(newEvents[i - 1].year, eventsStr);
        speciesAreaChart.addEvent(newEvents[i - 1].year, eventsStr);
        typeCntChart.addEvent(newEvents[i - 1].year, eventsStr);
        typeAreaChart.addEvent(newEvents[i - 1].year, eventsStr);
        i++;
    }

    chartLastYear = currentYear;

    speciesCntChart.getAxisX()->setRange(0, currentYear + 1);
    speciesAreaChart.getAxisX()->setRange(0, currentYear + 1);
    typeCntChart.getAxisX()->setRange(0, currentYear + 1);
    typeAreaChart.getAxisX()->setRange(0, currentYear + 1);

    speciesCntChart.getAxisY()->setRange(0, std::max(1.0f, 1.1f * maxSpeciesCnt));
    speciesAreaChart.getAxisY()->setRange(0, std::max(0.01f, 1.1f * maxSpeciesArea));
    typeCntChart.getAxisY()->setRange(0, std::max(1.0f, 1.1f * maxTypeCnt));
    typeAreaChart.getAxisY()->setRange(0, std::max(0.01f, 1.1f * maxTypeArea));

    speciesCntChart.update();
    speciesAreaChart.update();
    typeCntChart.update();
    typeAreaChart.update();
}

void MainWindow::resetCharts()
{
    speciesCntChart.reset();
    speciesAreaChart.reset();
    typeCntChart.reset();
    typeAreaChart.reset();

    speciesCntSeries.clear();
    speciesAreaSeries.clear();
    typeCntSeries.clear();
    typeAreaSeries.clear();

    ui.cbChartSelection->setCurrentIndex(0);
}

void MainWindow::initEventsList()
{
    ui.tableWidget->setColumnCount(3);
    ui.tableWidget->setHorizontalHeaderLabels({"Year", "Month", "Event Type"});
    ui.tableWidget->resizeColumnsToContents();
    ui.tableWidget->setSelectionBehavior(QAbstractItemView::SelectionBehavior::SelectRows);
    ui.tableWidget->setSelectionMode(QAbstractItemView::SelectionMode::SingleSelection);
    ui.tableWidget->setEditTriggers(QAbstractItemView::EditTrigger::NoEditTriggers);
}

void MainWindow::updateEventsList()
{
    ui.tableWidget->clearContents();
    ui.tableWidget->setRowCount(0);
    eventlist = ecosimWidget->getEvents();
    for (int i = 0; i < eventlist.size(); ++i) {
        ecosim::ScriptedDisturbanceEvent* e = eventlist.at(i);
        QString year = QString::number(e->year);
        QString month = QString::number(e->month);
        QString etype;
        switch (e->distEvent->type())
        {
        case ecosim::DisturbanceEventCategory::FIRE:
            etype = "Wildfire";
            break;
        case ecosim::DisturbanceEventCategory::WIND:
            etype = "Windstorm";
            break;
        case ecosim::DisturbanceEventCategory::DISEASE:
            etype = "Disease";
            break;
        case ecosim::DisturbanceEventCategory::DROUGHT:
            etype = "Drought";
            break;
        default:
            etype = "Error";
            break;
        }
        ui.tableWidget->insertRow(i);
        ui.tableWidget->setItem(i, 0, new QTableWidgetItem(year));
        ui.tableWidget->setItem(i, 1, new QTableWidgetItem(month));
        ui.tableWidget->setItem(i, 2, new QTableWidgetItem(etype));
    }
    ui.tableWidget->resizeColumnsToContents();
}

void MainWindow::setInfoText(const QString& s)
{
    ui.textTree->setText(s);
}

void MainWindow::setPickedPoint(const QPoint& p)
{
    int dx, dy;
    ecosimWidget->getTerrainGridDim(dx, dy);

    if (ui.btnSelectPtWf->isChecked()) {
        ui.btnSelectPtWf->setChecked(false);
        ui.dsbWfX->setValue(p.x() / float(dx));
        ui.dsbWfY->setValue(p.y() / float(dy));
    }
    else if (ui.btnSelectPtDs->isChecked()) {
        ui.btnSelectPtDs->setChecked(false);
        ui.dsbDsX->setValue(p.x() / float(dx));
        ui.dsbDsY->setValue(p.y() / float(dy));
    }
}

void MainWindow::resetEventsUI(bool unselectEvent)
{
    if (unselectEvent) {
        ui.tableWidget->clearSelection();
        editingEventIndex = -1;
    }

    ui.sbYear->setEnabled(true);
    ui.sbMonth->setEnabled(true);

    // TODO: get current year and month
    ui.sbYear->setValue(1);
    ui.sbMonth->setValue(1);

    // wildfire
    ui.dsbWfWindX->setValue(0);
    ui.dsbWfWindY->setValue(0);
    ui.lwOrgPtsWf->clear();

    //windstorm
    ui.dsbWsX->setValue(0);
    ui.dsbWsY->setValue(0);

    //disease
    ui.dsbDsX->setValue(0);
    ui.dsbDsY->setValue(0);
    ui.sbTargetSpecies->setValue(0);
    ui.dsbSpreadRadius->setValue(0);
    ui.dsbSpreadProb->setValue(0);
    ui.dsbSeverity->setValue(0);
    ui.sbMContagious->setValue(0);
    ui.sbMRecovery->setValue(0);
    ui.sbMImmune->setValue(0);

    //drought
    ui.dsbDrOrgPt->setValue(0);
    ui.lwDrValues->clear();

}

void MainWindow::loadEvent()
{
    if (ui.tableWidget->selectedItems().isEmpty()) {
        resetEventsUI(true);
        return;
    }

    resetEventsUI(false);
    int index = ui.tableWidget->selectedItems()[0]->row();
    ecosim::ScriptedDisturbanceEvent* e = eventlist.at(index);

    editingEventIndex = index;

    ui.sbYear->setValue(e->year);
    ui.sbMonth->setValue(e->month);

    ui.sbYear->setEnabled(false);
    ui.sbMonth->setEnabled(false);

    switch (e->distEvent->type())
    {
    case ecosim::DisturbanceEventCategory::FIRE:
    {
        //editingEvent = e;
        ecosim::Wildfire* wfev = dynamic_cast<ecosim::Wildfire*>(e->distEvent);
        float windX, windY;
        wfev->getWind(windX, windY);
        ui.dsbWfWindX->setValue(windX);
        ui.dsbWfWindY->setValue(windY);

        ui.lwOrgPtsWf->clear();
        std::vector<float> origins = wfev->getOrigins();
        for (int i = 0; i < origins.size(); i += 3) {
            QString itemText = QString("(%1, %2) r:%3")
                .arg(origins[i])
                .arg(origins[i + 1])
                .arg(origins[i + 2]);
            QListWidgetItem* item = new QListWidgetItem(itemText, ui.lwOrgPtsWf);
            item->setData(Qt::UserRole, QVariant(origins[i]));
            item->setData(Qt::UserRole + 1, QVariant(origins[i + 1]));
            item->setData(Qt::UserRole + 2, QVariant(origins[i + 2]));
            //ui.lwOrgPtsWf->addItem(itemText);
        }

        ui.tabEventCreator->setCurrentWidget(ui.tabWf);
    }
    break;
    case ecosim::DisturbanceEventCategory::WIND:
    {
        ecosim::Windstorm* wsev = dynamic_cast<ecosim::Windstorm*>(e->distEvent);

        float windX, windY;
        wsev->getWind(windX, windY);
        ui.dsbWsX->setValue(windX);
        ui.dsbWsY->setValue(windY);

        ui.tabEventCreator->setCurrentWidget(ui.tabWs);
    }
        break;
    case ecosim::DisturbanceEventCategory::DISEASE:
    {
        ecosim::DiseaseOutbreak* dsev = dynamic_cast<ecosim::DiseaseOutbreak*>(e->distEvent);

        float x, y;
        dsev->getOutbreakLocation(x, y);
        ui.dsbDsX->setValue(x);
        ui.dsbDsY->setValue(y);

        const ecosim::Disease* dparams = dsev->getDiseaseParams();
        ui.sbTargetSpecies->setValue(dparams->targetSpecies);
        ui.dsbSpreadRadius->setValue(dparams->spreadRadius);
        ui.dsbSpreadProb->setValue(dparams->spreadProbability);
        ui.dsbSeverity->setValue(dparams->severity);
        ui.sbMContagious->setValue(dparams->monthsContagious);
        ui.sbMRecovery->setValue(dparams->monthsRecovery);
        ui.sbMImmune->setValue(dparams->monthsImmune);

        ui.tabEventCreator->setCurrentWidget(ui.tabDs);
    }
        break;
    case ecosim::DisturbanceEventCategory::DROUGHT:
    {
        ecosim::Drought* drev = dynamic_cast<ecosim::Drought*>(e->distEvent);

        ui.lwDrValues->clear();

        std::vector<float> drValues = drev->getDroughtValues();
        for (int i = 0; i < drValues.size(); ++i) {
            QListWidgetItem* item = new QListWidgetItem(QString::number(drValues[i]), ui.lwDrValues);
            item->setData(Qt::UserRole, QVariant(drValues[i]));
        }

        ui.tabEventCreator->setCurrentWidget(ui.tabDr);
    }
    }
}

void MainWindow::setEvent(ecosim::DisturbanceEventCategory dec)
{

    ecosim::DisturbanceAction* da = nullptr;
    switch (dec)
    {
    case ecosim::DisturbanceEventCategory::FIRE:
    {
        ecosim::Wildfire* devt = new ecosim::Wildfire();

        float windx, windy;
        windx = ui.dsbWfWindX->value();
        windy = ui.dsbWfWindY->value();
        devt->setWind(windx, windy);

        for (int i = 0; i < ui.lwOrgPtsWf->count(); ++i)
        {
            QListWidgetItem* item = ui.lwOrgPtsWf->item(i);
            float x = item->data(Qt::UserRole).toFloat();
            float y = item->data(Qt::UserRole + 1).toFloat();
            float r = item->data(Qt::UserRole + 2).toFloat();
            devt->addOrigin(x, y, r);
        }
        da = devt;
    }
    break;
    case ecosim::DisturbanceEventCategory::WIND:
    {
        ecosim::Windstorm* devt = new ecosim::Windstorm();
        float windx, windy;
        windx = ui.dsbWsX->value();
        windy = ui.dsbWsY->value();
        devt->setWind(windx, windy);
        da = devt;
    }
        break;
    case ecosim::DisturbanceEventCategory::DISEASE:
    {
        ecosim::DiseaseOutbreak* devt = new ecosim::DiseaseOutbreak();
        float ox, oy;
        ox = ui.dsbDsX->value();
        oy = ui.dsbDsY->value();
        devt->setOutbreakLocation(ox, oy);

        ecosim::Disease* params = new ecosim::Disease();
        params->diseaseId = ++ecosim::ScriptedDisturbanceEvent::diseaseNum;
        params->targetSpecies = ui.sbTargetSpecies->value();
        params->spreadRadius = ui.dsbSpreadRadius->value();
        params->spreadProbability = ui.dsbSpreadProb->value();
        params->severity = ui.dsbSeverity->value();
        params->monthsContagious = ui.sbMContagious->value();
        params->monthsRecovery = ui.sbMRecovery->value();
        params->monthsImmune = ui.sbMImmune->value();
        params->totalInfected = 0;
        params->currentInfected = 0;
        devt->setDiseaseParams(params);

        da = devt;
    }
        break;
    case ecosim::DisturbanceEventCategory::DROUGHT:
    {
        ecosim::Drought* devt = new ecosim::Drought();

        std::vector<float> vals(ui.lwDrValues->count());
        for (int i = 0; i < ui.lwDrValues->count(); ++i)
        {
            QListWidgetItem* item = ui.lwDrValues->item(i);
            vals[i] = item->data(Qt::UserRole).toFloat();
        }
        devt->setDroughtValues(vals);

        da = devt;
    }
        break;
    default:
        break;
    }

    ecosim::ScriptedDisturbanceEvent* sdev = new ecosim::ScriptedDisturbanceEvent();
    sdev->year = ui.sbYear->value();
    sdev->month = ui.sbMonth->value();
    sdev->distEvent = da;
    // edit event
    if (editingEventIndex >= 0) {
        ecosimWidget->editEvent(editingEventIndex, sdev);
    }
    // new event
    else {
        ecosimWidget->addEvent(sdev);
    }

    updateEventsList();
}

void MainWindow::removeEvent(int index)
{
    ecosimWidget->removeEvent(index);
}

void MainWindow::saveEventList()
{
    QString filepath = sceneDir + QString("/scriptedEvents.txt");
    ecosimWidget->saveEventList(eventlist, filepath.toStdString());
    saveEventListAnim->start();
}

void MainWindow::loadEventList()
{
    QString event_filepath = QFileDialog::getOpenFileName(
        this, "Open events file", sceneDir, tr("text files (*.txt)")
    );

    ecosimWidget->setEventList(event_filepath.toStdString());

    updateEventsList();
}

void MainWindow::updateFilters()
{
    float min_h = (ui.cbMinHeight->isChecked()) ? ui.dsbMinHeight->value() : 0;
    float max_h = (ui.cbMaxHeight->isChecked()) ? ui.dsbMaxHeight->value() : -1;
    int min_age = (ui.cbMinAge->isChecked()) ? ui.sbMinAge->value() : 0;
    int max_age = (ui.cbMaxAge->isChecked()) ? ui.sbMaxAge->value() : -1;
    ecosimWidget->setFilters(min_h, max_h, min_age, max_age);
}

void MainWindow::changeChart(EventChart* c)
{
    EventChart* oldChart = dynamic_cast<EventChart*>(ui.graphicsView->chart());
    QGraphicsScene* scene = oldChart->scene();
    for (auto item : oldChart->events) {
        scene->removeItem(item);
    }

    ui.graphicsView->setChart(c);
    for (auto item : c->events) {
        c->scene()->addItem(item);
    }
}



void MainWindow::fillSpeciesCombobox()
{
    std::vector<QString> names = ecosimWidget->getSpeciesNames();

    ui.cbTexSpecies->clear();
    for (const QString& s : names) {
        ui.cbTexSpecies->addItem(s);
    }
}

void MainWindow::createSpeciesCheckbox()
{
    std::vector<QString> names = ecosimWidget->getSpeciesNames();
    QVBoxLayout* layout = new QVBoxLayout();
    for (int i = 0; i < names.size(); ++i) {
        QString s = names[i];
        QCheckBox* checkBox = new QCheckBox(s);
        checkBox->setChecked(true);

        // QLabel with background color
        QLabel* colorLabel = new QLabel();
        QColor color = ecosimWidget->getColorForSpecies(i);
        colorLabel->setStyleSheet("background-color: " + color.name() + "; border: 1px solid black;");
        colorLabel->setFixedWidth(20);
        colorLabel->setFixedHeight(20);

        QHBoxLayout* hLayout = new QHBoxLayout();
        hLayout->addWidget(colorLabel);
        hLayout->addWidget(checkBox);

        layout->addLayout(hLayout);

        connect(checkBox, &QCheckBox::clicked, this, [this, i] (bool b) {
            ecosimWidget->toggleRenderSpecie(i, b); 
        });
    }
    ui.groupBox_7->setLayout(layout);
}
