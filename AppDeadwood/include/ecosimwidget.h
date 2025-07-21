#ifndef __EcosimWidget__
#define __EcosimWidget__

#include "glheaders.h"
#include <QtCore/QElapsedTimer>
#include <QtCore/QThread>
#include <QtGui/QKeyEvent>
#include <QtGui/QMouseEvent>
#include <QtOpenGLWidgets/QOpenGLWidget>
#include <QtOpenGL/QOpenGLShaderProgram>
#include "core.h"
#include "simWorker.h"
#include "sim.h"
#include "sun.h"
#include "terrain.h"
#include "disturbance_event.h"


class EcosimWidget : public QOpenGLWidget
{
	Q_OBJECT

public:
	enum TerrainTexture {
		PLAIN,
		ELEVATION,
		NORMALS,
		ASPECT,
		ASPECT_E,
		ASPECT_N,
		SLOPE,
		SUNLIGHT,
		MOISTURE,
		TEMPERATURE,
		VIABILITY,
		MOSTVIABLE
	};

	EcosimWidget(QWidget* parent = nullptr);
	virtual ~EcosimWidget();

	void runSim();
	void pauseSim();
	void resetSim();
	void stepSim(int months = 1);

	bool loadScene(const std::string& base_data_dir);

	int getNumSpecies() const;
	QColor getColorForSpecies(int index);
	std::vector<QString> getSpeciesNames() const;

	void getTerrainGridDim(int& dx, int& dy) const;
	void getTerrainDim(float& tx, float& ty) const;
	float getTerrainHeight(int x, int y) const;
	bool getRayTerrainIntersection(const Ray& ray, QPoint& hitCell) const;
	Camera getCamera() const;

	const std::vector<ecosim::ScriptedDisturbanceEvent*> getEvents() const;
	ecosim::ScriptedDisturbanceEvent* getEvent(int index) const;
	void addEvent(ecosim::ScriptedDisturbanceEvent* e) const;
	void editEvent(int index, ecosim::ScriptedDisturbanceEvent* e) const;
	void removeEvent(int index) const;
	void setEventList(const std::string& dir);
	std::vector<ecosim::DisturbanceEvent> getPastNewEvents();
	bool saveEventList(const std::vector<ecosim::ScriptedDisturbanceEvent*>& events, const std::string& filepath);

	int getSimulatedYears() const;
	std::vector<float> getDataCount(int year, int type) const;
	std::vector<float> getDataBasal(int year, int type) const;

	bool exportSim(const std::string& filename, bool terse=true) const;

public slots:
	virtual void mousePressEvent(QMouseEvent* e);
	virtual void mouseMoveEvent(QMouseEvent* e);
	virtual void mouseReleaseEvent(QMouseEvent* e);
	virtual void wheelEvent(QWheelEvent* e);

	void updatePlantInstances(bool falseColor = false);

	void resetCamera();
	void setRenderShadeTerrain(bool b);
	void setRenderShowStats(bool b);
	void setPlantColor(int t);
	void toggleRenderPlts(bool);
	void toggleRenderSnags(bool);
	void toggleRenderLogs(bool);
	void toggleRenderSpecie(int, bool);
	void toggleTerrainPicking(bool);
	void setRenderTexture(TerrainTexture t, int subtype = -1);
	void setFilters(float min_h, float max_h, int min_age, int max_age);

signals:
	void _signalMousePress(QMouseEvent*);
	void _signalMouseMoved(QMouseEvent*);
	void _signalClickedTerrainCell(QPoint);
	void _signalOperate();
	void _signalFinishedYear();
	void _signalTextInfo(const QString&);
	void _signalResetSim();

protected:

	class PlantInstanceData {
		  GLfloat color[4];
		  GLfloat transform[16];
	public:
			PlantInstanceData(const QColor& c, const QMatrix4x4& m);
	};

	// paint
	virtual void initializeGL();
	virtual void paintGL();

	void initPlantModels();
	void deletePlantModelBuffers();
	void initTerrainGeometry();
	void updateTexture(TerrainTexture texType, int subtype = -1);
	void updateTexture();
	void updateCameraPlanes();
	void paintVegetationIds();
	void paintTextOverlay();

	// simulation step
	void callSimulationStep();
	void simulationStepFinished();

	// helpers
	Vector3 ecosimToWidgetCoords(const ecosim::vpPoint& c, bool clampToGround = true) const;
	QColor encodeId(int type, int id) const;
	unsigned int decodeId(const QColor& c, int& type) const;

	std::string get_dirprefix();

	ecosim::Terrain* getTerrain();
	ecosim::EcoSystem* getEcoSys();
	ecosim::Simulation* getSim();
	MapFloat* getSunlight(int month);
	MapFloat* getSlope();
	MapFloat* getMoisture(int month);
	ecosim::Biome* getBiome();
	ecosim::GLSun* getGLSun();

	void loadPlantInfo(ecosim::SimLivePlant* plt);
	void loadPlantInfo(ecosim::SimDeadPlant* plt);
	void loadPlantInfo(ecosim::SimLog* plt);
	std::string PlantSimState_toString(ecosim::PlantSimState pss);
	std::string SnagDecayClass_toString(ecosim::SnagDecayClass sdc);

	bool filtersOk(ecosim::SimLivePlant* plt);


protected:
	bool isSimRunning = false;
	int  doSimSteps = 0;
	int  simulatedMonths = 0;
	double simTime = 0;
	double simPerf = 0;
	QElapsedTimer timer;

	QThread simThread;
	SimWorker* simWorker;
	bool isSimWorking = false;
	bool waitingForReset = false;
	bool pickingActive = false;
	bool paintingFalseColor = false;

	bool simvalid; //< whether or not a simulation can be run
	bool firstsim; //< is this the first simulation run

	std::string datadir;
	std::string biome_db_filepath; //< biome database

	// Scene
	ecosim::Terrain* terrain;
	Box3 terrainBBox;
	ecosim::EcoSystem* eco;
	ecosim::Biome* biome;
	ecosim::Simulation* sim;
	ecosim::GLSun* glsun;

	QOpenGLShaderProgram shaderTerrain, shaderSkybox, shaderVegInstances, shaderSelection;
	GLuint meshVAO = 0;
	GLuint bufferVerts = 0, bufferIndices = 0;
	GLuint numTerrainTriangles = 0;
	QImage texImg;
	GLuint texId = 0;
	GLuint skyboxVAO = 0;
	GLuint trunkVAO = 0;
	GLuint trunkVAOSize = 0;
	std::vector<GLuint> pftVAO;
	std::vector<GLuint> pftVAOSize;
	std::vector<GLuint> pftBuffers;
	GLuint instanceVBO = 0;

	Camera camera;
	int x0 = 0, y0 = 0;

	TerrainTexture renderTexType = ELEVATION;
	int plantColorType = 0;
	int renderTexSubtype = -1;
	bool renderShadeTerrain = true;
	bool renderShowStats = true;	
	bool renderPlts, renderSnags, renderLogs;	
	std::vector<bool> renderSpecies;
	Vector3 lightPosition;
	int selectedPlantId = -1, selectedPlantType = -1;
	std::vector<std::vector<PlantInstanceData> > instancesCanopy;
	std::vector<PlantInstanceData> instancesTrunks;
	std::vector<PlantInstanceData> instancesLogs;

	QImage viabilitiesTex;
	bool viabilitiesTexComputed = false;
	bool updatedMonthlyTexture = false;

	// filters
	float minHeight = 0;
	float maxHeight = FLT_MAX;
	int minAge = 0;
	int maxAge = INT32_MAX;
};


inline void EcosimWidget::setRenderTexture(TerrainTexture t, int sub)
{
	renderTexType = t;
	renderTexSubtype = sub;
	updateTexture(t, sub);
	update();
}

inline void EcosimWidget::setRenderShadeTerrain(bool b)
{
	renderShadeTerrain = b;
	updateTexture();
	update();
}

inline void EcosimWidget::setRenderShowStats(bool b)
{
	renderShowStats = b;
	update();
}

inline void EcosimWidget::setPlantColor(int t)
{
	plantColorType = t;
	update();
}

inline Vector3 EcosimWidget::ecosimToWidgetCoords(const ecosim::vpPoint& c, bool clampToGround) const
{	
	const double x = c.x;
	const double y = c.z;
	const double z = clampToGround ? terrain->getHeightAtCoords(c.x, c.z, true) : c.y;
	return Vector3(x, y, z);
}

inline int EcosimWidget::getNumSpecies() const {
	if (sim && sim->getBiome()) return sim->getBiome()->numPFTypes();
	else return 0;
}

inline std::vector<float> EcosimWidget::getDataCount(int year, int type) const {
	return sim->getDataCount(year, type);
}

inline std::vector<float> EcosimWidget::getDataBasal(int year, int type) const {
	return sim->getDataBasal(year, type);
}

inline int EcosimWidget::getSimulatedYears() const {
	return simulatedMonths / 12;
}

inline void EcosimWidget::toggleRenderPlts(bool b)
{
	this->renderPlts = b; updatePlantInstances();
}

inline void EcosimWidget::toggleRenderSnags(bool b)
{
	this->renderSnags = b; updatePlantInstances();
}

inline void EcosimWidget::toggleRenderLogs(bool b)
{
	this->renderLogs = b; updatePlantInstances();
}

inline void EcosimWidget::toggleRenderSpecie(int i, bool b)
{
	this->renderSpecies[i] = b; updatePlantInstances();
}

inline void EcosimWidget::toggleTerrainPicking(bool b)
{
	pickingActive = b;
}

inline std::vector<QString> EcosimWidget::getSpeciesNames() const
{
	std::vector<std::string> names_std = sim->getSpeciesNames();
	std::vector<QString> names(names_std.size());
	for (int i = 0; i < names_std.size(); ++i) {
		names[i] = QString::fromStdString(names_std[i]);
	}
	return names;
}

inline ecosim::Terrain* EcosimWidget::getTerrain()
{
	return terrain;
}

inline ecosim::EcoSystem* EcosimWidget::getEcoSys()
{
	return eco;
}

inline ecosim::Simulation* EcosimWidget::getSim()
{
	return sim;
}

inline MapFloat* EcosimWidget::getSunlight(int month)
{
	return sim->getSunlightMap(month);
}

inline MapFloat* EcosimWidget::getSlope()
{
	return sim->getSlopeMap();
}

inline MapFloat* EcosimWidget::getMoisture(int month)
{
	return sim->getMoistureMap(month);
}

inline ecosim::Biome* EcosimWidget::getBiome()
{
	return biome;
}

inline ecosim::GLSun* EcosimWidget::getGLSun()
{
	return glsun;
}

inline const std::vector<ecosim::ScriptedDisturbanceEvent*> EcosimWidget::getEvents() const {
	if (sim == nullptr) return std::vector<ecosim::ScriptedDisturbanceEvent*>();
	return sim->getEvents();
}

inline ecosim::ScriptedDisturbanceEvent* EcosimWidget::getEvent(int index) const {
	if (sim == nullptr) return nullptr;
	return sim->getEvent(index);
}

inline void EcosimWidget::addEvent(ecosim::ScriptedDisturbanceEvent* e) const {
	if (sim != nullptr)
		sim->addEvent(e);
}

inline void EcosimWidget::editEvent(int index, ecosim::ScriptedDisturbanceEvent* e) const {
	if (sim != nullptr)
		sim->editEvent(index, e);
}

inline void EcosimWidget::removeEvent(int index) const {
	if (sim != nullptr)
		sim->removeEvent(index);
}

inline void EcosimWidget::setEventList(const std::string& dir)
{
	getSim()->setEvents(ecosim::ScriptedDisturbanceEvent::loadFromFile(dir));
}

inline std::vector<ecosim::DisturbanceEvent> EcosimWidget::getPastNewEvents()
{
	return sim->getPastNewEvents();
}

inline QColor EcosimWidget::getColorForSpecies(int index)
{
	GLfloat* col = sim->getBiome()->getPFType(index)->basecol;
	return toQColor(Vector3(col[0], col[1], col[2]));
}

inline bool EcosimWidget::saveEventList(const std::vector<ecosim::ScriptedDisturbanceEvent*>& events, const std::string& filepath)
{
	return ecosim::ScriptedDisturbanceEvent::saveToFile(events, filepath);
}

inline void EcosimWidget::getTerrainGridDim(int& dx, int& dy) const
{
	terrain->getGridDim(dx, dy);
}

inline void EcosimWidget::getTerrainDim(float& tx, float& ty) const
{
	terrain->getTerrainDim(tx, ty);
}

inline float EcosimWidget::getTerrainHeight(int x, int y) const
{
	return terrain->getHeight(x, y);
}

inline Camera EcosimWidget::getCamera() const
{
	return camera;
}

inline bool EcosimWidget::exportSim(const std::string& filename, bool terse) const
{
	if (sim != nullptr) {
		if (filename.ends_with("json")) return sim->exportSimJSON(filename);
		else if (filename.ends_with("pdb")) {
			if (terse) return sim->exportSimPDBTerse(filename);
			else return sim->exportSimPDBFull(filename);
		}
		else return false;
	}
	else return false;
}


#endif