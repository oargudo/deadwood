#include "ecosimwidget.h"
#include <QtGui/QPainter>
#include <QtGui/QPen>
#include <QtGui/QOpenGLContext>
#include "disturbance_event.h"
#include "model.h"


EcosimWidget::EcosimWidget(QWidget* parent) : QOpenGLWidget(parent)
{
    sim = nullptr;
    eco = nullptr;
    terrain = nullptr;
    biome = nullptr;
    glsun = nullptr;

    simWorker = nullptr;

    simvalid = false;
    firstsim = true;

    renderPlts = true;
    renderSnags = true;
    renderLogs = true;
}

EcosimWidget::~EcosimWidget()
{
    simThread.quit();
    simThread.wait();
    if (simWorker != nullptr) delete simWorker;
    if (sim != nullptr) delete sim;
    if (eco != nullptr) delete eco;
    if (biome != nullptr) delete biome;
    if (terrain != nullptr) delete terrain;
    deletePlantModelBuffers();
}


void EcosimWidget::runSim()
{
    this->makeCurrent();
	  isSimRunning = true;
    update();
}

void EcosimWidget::pauseSim()
{
    this->makeCurrent();
	  isSimRunning = false;
    update();
}

void EcosimWidget::resetSim()
{
    if (sim == nullptr)
        return;

    if (isSimWorking) {
        waitingForReset = true;
        return;
    }
    waitingForReset = false;

    this->makeCurrent();

	simulatedMonths = 0;
    doSimSteps = 0;
    simTime = 0;
    simPerf = 0;
    isSimRunning = false;
    updatedMonthlyTexture = false;

    selectedPlantId = -1;
    selectedPlantType = 0;
    getSim()->setEvents(ecosim::ScriptedDisturbanceEvent::loadFromFile(datadir + "/scriptedEvents.txt"));
    getSim()->setupSim(getEcoSys(), false);

    updatePlantInstances();

    emit _signalResetSim();
}

void EcosimWidget::stepSim(int numSteps)
{
    this->makeCurrent();
	  isSimRunning = true;
    doSimSteps = numSteps;
	  update();
}

void EcosimWidget::initializeGL()
{
    if (!context()->isValid()) {
        std::cerr << "OpenGL context is not valid!" << std::endl;
        return;
    }

    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW!" << std::endl;
        return;
    }

    shaderTerrain.addShaderFromSourceFile(QOpenGLShader::Vertex, ":/shaders/terrain.vert");
    shaderTerrain.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/shaders/terrain.frag");
    shaderTerrain.link();
    shaderSkybox.addShaderFromSourceFile(QOpenGLShader::Vertex, ":/shaders/skybox.vert");
    shaderSkybox.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/shaders/skybox.frag");
    shaderSkybox.link();
    shaderVegInstances.addShaderFromSourceFile(QOpenGLShader::Vertex, ":/shaders/vegInstances.vert");
    shaderVegInstances.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/shaders/vegInstances.frag");
    shaderVegInstances.link();
    shaderSelection.addShaderFromSourceFile(QOpenGLShader::Vertex, ":/shaders/vegInstancesId.vert");
    shaderSelection.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/shaders/vegInstancesId.frag");
    shaderSelection.link();
    CE();

    glGenVertexArrays(1, &skyboxVAO);
    CE();

    unsigned char foo = 255;
    glGenTextures(1, &texId);
    glBindTexture(GL_TEXTURE_2D, texId);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 1, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, (void*)(&foo));
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);
    CE();

    glsun = new ecosim::GLSun();
    glsun->initializeGL();
    CE();
}

void EcosimWidget::paintGL()
{
    if (paintingFalseColor) return;

    this->makeCurrent();
    CE();

    // For sunlight and moisture monthly maps, update according to curr month
    if ((renderTexType == SUNLIGHT || renderTexType == MOISTURE || renderTexType == TEMPERATURE) && !updatedMonthlyTexture) {
        updateTexture(renderTexType);
    }

    QMatrix4x4 matPerspective;
    matPerspective.perspective(camera.getAngleOfViewV(width(), height()) * 180 / M_PI,
                              (GLdouble)width() / (GLdouble)height(),
                              camera.getNearPlane(), camera.getFarPlane());

    QMatrix4x4 matView;
    matView.lookAt(QVector3D(camera.getEye()[0], camera.getEye()[1], camera.getEye()[2]),
                   QVector3D(camera.getAt()[0], camera.getAt()[1], camera.getAt()[2]),
                   QVector3D(camera.getUp()[0], camera.getUp()[1], camera.getUp()[2]));


    // Clear
    glClearColor(0.62f, 0.74f, 0.85f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Sky
    CE();
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glUseProgram(shaderSkybox.programId());
    glBindVertexArray(skyboxVAO);
    glUniform3f(0, camera.getEye()[0], camera.getEye()[1], camera.getEye()[2]);
    glUniform3f(1, camera.getAt()[0], camera.getAt()[1], camera.getAt()[2]);
    glUniform3f(2, camera.getUp()[0], camera.getUp()[1], camera.getUp()[2]);
    glUniform2f(3, width(), height());
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    CE();

    // Terrain
    if (meshVAO > 0) {
        GLuint shader = shaderTerrain.programId();
        glUseProgram(shader);

        glUniformMatrix4fv(glGetUniformLocation(shader, "ProjectionMatrix"), 1, GL_FALSE, matPerspective.constData());
        glUniformMatrix4fv(glGetUniformLocation(shader, "ModelViewMatrix"), 1, GL_FALSE, matView.constData());
        glUniform2f(glGetUniformLocation(shader, "u_worldMin"), terrainBBox.getMin()[0], terrainBBox.getMin()[1]);
        glUniform2f(glGetUniformLocation(shader, "u_worldSize"), terrainBBox.width(), terrainBBox.height());
        
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texId);
        glUniform1i(glGetUniformLocation(shader, "u_texture"), 0);

        glBindVertexArray(meshVAO); CE();
        glDrawElements(GL_TRIANGLES, numTerrainTriangles * 3, GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);

        glUseProgram(0);
        CE();
    }

    // Draw vegetation using instancing
    if (instanceVBO > 0) {
        GLuint shader = shaderVegInstances.programId();
        glUseProgram(shader);

        QVector3D lightPos = matView.map(QVector3D(2.0, 2.0, 8.0)*terrainBBox.radius());
        glUniform3f(glGetUniformLocation(shader, "lightPos"), lightPos.x(), lightPos.y(), lightPos.z());
        glUniformMatrix4fv(glGetUniformLocation(shader, "ProjectionMatrix"), 1, GL_FALSE, matPerspective.constData());
        glUniformMatrix4fv(glGetUniformLocation(shader, "ModelViewMatrix"), 1, GL_FALSE, matView.constData());

        // trunks and logs
        glBindVertexArray(trunkVAO);
        glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, instancesTrunks.size() * sizeof(PlantInstanceData), instancesTrunks.data());
        glDrawElementsInstanced(GL_TRIANGLES, trunkVAOSize, GL_UNSIGNED_INT, 0, static_cast<GLsizei>(instancesTrunks.size()));
        glBufferSubData(GL_ARRAY_BUFFER, 0, instancesLogs.size() * sizeof(PlantInstanceData), instancesLogs.data());
        glDrawElementsInstanced(GL_TRIANGLES, trunkVAOSize, GL_UNSIGNED_INT, 0, static_cast<GLsizei>(instancesLogs.size()));
    
        // per-species canopies
        for (int pft = 0; pft < pftVAO.size(); pft++) {
            glBindVertexArray(pftVAO[pft]);
            glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
            glBufferSubData(GL_ARRAY_BUFFER, 0, instancesCanopy[pft].size() * sizeof(PlantInstanceData), instancesCanopy[pft].data());
            glDrawElementsInstanced(GL_TRIANGLES, pftVAOSize[pft], GL_UNSIGNED_INT, 0, static_cast<GLsizei>(instancesCanopy[pft].size()));
        }

        CE();
    }

    glBindVertexArray(0);
    glUseProgram(0);
    CE();


    // Then draw our text overlay
    if (renderShowStats)
        paintTextOverlay();


    // Simulate and call next frame draw
    if (isSimRunning && !isSimWorking) {
        // do one simulation step
        callSimulationStep();
	  }

    CE();
}

void EcosimWidget::paintVegetationIds()
{
    QMatrix4x4 matPerspective;
    matPerspective.perspective(camera.getAngleOfViewV(width(), height()) * 180 / M_PI,
      (GLdouble)width() / (GLdouble)height(),
      camera.getNearPlane(), camera.getFarPlane());

    QMatrix4x4 matView;
    matView.lookAt(QVector3D(camera.getEye()[0], camera.getEye()[1], camera.getEye()[2]),
      QVector3D(camera.getAt()[0], camera.getAt()[1], camera.getAt()[2]),
      QVector3D(camera.getUp()[0], camera.getUp()[1], camera.getUp()[2]));

    // Save the current blending state
    GLboolean blendEnabled = glIsEnabled(GL_BLEND);
    GLint srcRGB, dstRGB, srcAlpha, dstAlpha;
    GLint blendEquationRGB, blendEquationAlpha;
    glGetIntegerv(GL_BLEND_SRC_RGB, &srcRGB);
    glGetIntegerv(GL_BLEND_DST_RGB, &dstRGB);
    glGetIntegerv(GL_BLEND_SRC_ALPHA, &srcAlpha);
    glGetIntegerv(GL_BLEND_DST_ALPHA, &dstAlpha);
    glGetIntegerv(GL_BLEND_EQUATION_RGB, &blendEquationRGB);
    glGetIntegerv(GL_BLEND_EQUATION_ALPHA, &blendEquationAlpha);

    glDisable(GL_BLEND);

    // Clear
    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    GLuint shader = shaderSelection.programId();
    glUseProgram(shader);

    glUniformMatrix4fv(glGetUniformLocation(shader, "ProjectionMatrix"), 1, GL_FALSE, matPerspective.constData());
    glUniformMatrix4fv(glGetUniformLocation(shader, "ModelViewMatrix"), 1, GL_FALSE, matView.constData());

    glBindVertexArray(trunkVAO);
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);

    glBufferSubData(GL_ARRAY_BUFFER, 0, instancesTrunks.size() * sizeof(PlantInstanceData), instancesTrunks.data());
    glDrawElementsInstanced(GL_TRIANGLES, trunkVAOSize, GL_UNSIGNED_INT, 0, static_cast<GLsizei>(instancesTrunks.size()));

    glBufferSubData(GL_ARRAY_BUFFER, 0, instancesLogs.size() * sizeof(PlantInstanceData), instancesLogs.data());
    glDrawElementsInstanced(GL_TRIANGLES, trunkVAOSize, GL_UNSIGNED_INT, 0, static_cast<GLsizei>(instancesLogs.size()));

    for (int pft = 0; pft < pftVAO.size(); pft++) {
      glBindVertexArray(pftVAO[pft]);
      glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
      glBufferSubData(GL_ARRAY_BUFFER, 0, instancesCanopy[pft].size() * sizeof(PlantInstanceData), instancesCanopy[pft].data());
      glDrawElementsInstanced(GL_TRIANGLES, pftVAOSize[pft], GL_UNSIGNED_INT, 0, static_cast<GLsizei>(instancesCanopy[pft].size()));
    }

    glBindVertexArray(0);
    glUseProgram(0);

    if (blendEnabled) glEnable(GL_BLEND);
    glBlendFuncSeparate(srcRGB, dstRGB, srcAlpha, dstAlpha);
    glBlendEquationSeparate(blendEquationRGB, blendEquationAlpha);
}

void EcosimWidget::paintTextOverlay()
{
    // We need to unbind VAO and program for now because it causes problem with commands below
    glUseProgram(0);
    glBindVertexArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    QPainter painter;
    painter.begin(this);

    painter.setRenderHint(QPainter::Antialiasing);
    QPen penLineGrey(QColor(50, 50, 50));
    QPen penLineWhite(QColor(250, 250, 250));

    constexpr int bX = 10;
    constexpr int bY = 10;
    constexpr int sizeX = 170;
    constexpr int sizeY = 130;

    // Background
    painter.setPen(penLineGrey);
    painter.fillRect(QRect(bX, bY, sizeX, sizeY), QColor(0, 0, 255, 25));
    painter.drawRect(bX, bY, sizeX, sizeY);

    // Text
    QFont f;
    f.setPointSize(10);
    f.setFamily("Arial");
    painter.setPen(penLineWhite);
    painter.setFont(f);

    int livePlantsCnt = 0;
    int deadPlantsCnt = 0;
    int logsCnt = 0;
    if (sim != nullptr) {
        livePlantsCnt = int(sim->getLivePlants().size());
        deadPlantsCnt = int(sim->getDeadPlants().size());
        logsCnt = int(sim->getLogs().size());
    }

    painter.drawText(bX + 5, bY + 10 + 10, "Year " + QString::number((simulatedMonths - 1)/12)
                                        + " month " + QString::number(simulatedMonths > 0 && simulatedMonths%12 == 0 ? 12 : (simulatedMonths) % 12));
    painter.drawText(bX + 5, bY + 10 + 30, "#Trees: " + QString::number(livePlantsCnt));
    painter.drawText(bX + 5, bY + 10 + 50, "#Snags: " + QString::number(deadPlantsCnt));
    painter.drawText(bX + 5, bY + 10 + 70, "#Logs:  " + QString::number(logsCnt));
    painter.drawText(bX + 5, bY + 10 + 90, "Curr perf: " + QString::number(simPerf, 'f', 1) + " s/month");
    painter.drawText(bX + 5, bY + 10 + 110,"Avg  perf: " + QString::number(simTime / double(simulatedMonths), 'f', 1) + " s/month");
    painter.end();

    // Reset GL depth test and alpha
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void EcosimWidget::resetCamera()
{
    // adapt camera far plane to the terrain
    if (terrain) {
        double R = terrainBBox.radius();
        camera = Camera(
            terrainBBox.center() + 4 * R * Normalized(Vector3(0.75, -1.0, 1.0)),
            terrainBBox.center(),
            Vector3(0, 0, 1)
        );
        updateCameraPlanes();
    }
    update();
}

void EcosimWidget::updateCameraPlanes()
{
    // update camera planes
    if (terrain) {
        double camDist = Norm(camera.getEye() - terrainBBox.center());
        double tradius = terrainBBox.radius();
        if (camDist <= tradius) {
            camera.setPlanes(std::max(0.01, 0.0001*tradius), 2 * tradius);
        }
        else {
            camera.setPlanes(camDist - tradius, 2 * tradius + camDist);
        }
    }
    else {
        camera.setPlanes(0.01, 1);
    }
}


void EcosimWidget::loadPlantInfo(ecosim::SimLivePlant* plt)
{
    //plt->state
    string info = std::format(
        "LIVE PLANT: {}\n"
        "Species: {}\n"
        "Age: {} years\n"
        "Nominal age: {:.3f}\n"
        "Height: {:.3f} m\n"
        "Canopy: {:.3f} m\n"
        "DBH: {:.3f} m\n"
        "Vigour: {:.3f}",
      plt->id, biome->getPFType(plt->pft)->code, plt->age, plt->nomage,
      plt->height, plt->canopy, plt->dbh, plt->vigourhistory.back());

    emit _signalTextInfo(QString::fromStdString(info));
}

void EcosimWidget::loadPlantInfo(ecosim::SimDeadPlant* plt)
{
    //plt->state
    string info = std::format(
        "SNAG: {}\n"
        "Decay class: {}\n"
        "Species: {}\n"
        "Age: {} years ({} at death)\n"
        "Height: {:.3f} m\n"
        "Canopy: {:.3f} m\n"
        "Decay: {:.3f}\n",
        plt->id, SnagDecayClass_toString(plt->dc), biome->getPFType(plt->parent->pft)->code, plt->age, plt->parent->age - plt->age,
        plt->parent->height * plt->trunkremains, plt->parent->canopy, plt->decay);

    emit _signalTextInfo(QString::fromStdString(info));
}

void EcosimWidget::loadPlantInfo(ecosim::SimLog* plt)
{
    //plt->state
    string info = std::format(
        "LOG: {}\n"
        "Age: {} years\n"
        "Diameter: {:.3f} m\n"
        "Length: {:.3f} m\n"
        , plt->id, plt->age, plt->diam, plt->len);

    emit _signalTextInfo(QString::fromStdString(info));
}

std::string EcosimWidget::PlantSimState_toString(ecosim::PlantSimState pss) {
    switch (pss) {
        case ecosim::PlantSimState::ALIVE: return "Alive";
        case ecosim::PlantSimState::SPROUT: return "Sprout";
        case ecosim::PlantSimState::DEAD: return "Dead";
        case ecosim::PlantSimState::STATIC: return "Static";
        default: return "Unknown";
    }
}

std::string EcosimWidget::SnagDecayClass_toString(ecosim::SnagDecayClass pss) {
    switch (pss) {
        case ecosim::SnagDecayClass::INTACT: return "Intact";
        case ecosim::SnagDecayClass::BARE: return "Bare";
        case ecosim::SnagDecayClass::DEBRANCHED: return "Debranched";
        case ecosim::SnagDecayClass::SNAPPED: return "Snapped";
        case ecosim::SnagDecayClass::STUMP: return "Stump";
        default: return "Unknown";
    }
}

void EcosimWidget::mousePressEvent(QMouseEvent* e)
{
    x0 = e->globalPosition().x();
    y0 = e->globalPosition().y();

    if (sim != nullptr && !e->modifiers() && e->buttons() & Qt::RightButton) {

        this->makeCurrent();
        paintingFalseColor = true;
        updatePlantInstances(true);
        paintVegetationIds();

        qreal dpr = this->devicePixelRatioF();
        int x = int(e->position().x() * dpr);
        int y = int((height() - e->position().y()) * dpr);
        unsigned char pixelColor[4];
        glReadPixels(x, y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, pixelColor);
        paintingFalseColor = false;

        int type = -1;
        unsigned int index = decodeId(QColor(pixelColor[0], pixelColor[1], pixelColor[2], pixelColor[3]), type);

        selectedPlantId = -1;
        selectedPlantType = 0;
        if (index > 0 && type >= 1 && type <= 3) {
            selectedPlantType = type;
            switch (type) {
                case 1: if (index < sim->getLivePlants().size()) selectedPlantId = sim->getLivePlants()[index]->id; break;
                case 2: if (index < sim->getDeadPlants().size()) selectedPlantId = sim->getDeadPlants()[index]->id; break;
                case 3: if (index < sim->getLogs().size()) selectedPlantId = sim->getLogs()[index]->id; break;
                default: break;
            }
        }
        else {
            emit _signalTextInfo("");
        }

        updatePlantInstances(false);
        update();
    }
}

void EcosimWidget::mouseMoveEvent(QMouseEvent* e) 
{
    int x = e->globalPosition().x();
    int y = e->globalPosition().y();

    double moveScale = Norm(camera.getEye() - camera.getAt()) * 0.015 * 0.1;

    if (e->buttons() & Qt::LeftButton) {
        if (e->modifiers() & Qt::ShiftModifier) {
            camera.backForth((y - y0) * moveScale);
            updateCameraPlanes();
        }
        else if (e->modifiers() & Qt::AltModifier) {
            camera.leftRightPlane((x - x0) * moveScale);
            camera.upDownPlane((y - y0) * moveScale);
        }
        else {
            camera.leftRightRound((x0 - x) * 0.01);
            camera.upDownRound((y0 - y) * 0.005);
        }
    }
    else if (e->buttons() & Qt::MiddleButton) {
        camera.leftRightPlane((x - x0) * moveScale);
        camera.upDownPlane((y - y0) * moveScale);
    }

    x0 = e->globalPosition().x();
    y0 = e->globalPosition().y();

    update();
}

void EcosimWidget::mouseReleaseEvent(QMouseEvent* e)
{
    int x = e->globalPosition().x();
    int y = e->globalPosition().y();

    if (e->button() & Qt::LeftButton) {
        if (pickingActive) {
            Ray ray = camera.pixelToRay(e->position().x(), e->position().y() - 1, width(), height());
            QPoint phit;
            if (getRayTerrainIntersection(ray, phit)) {
            emit _signalClickedTerrainCell(phit);
            }
        }
    }

    x0 = e->globalPosition().x();
    y0 = e->globalPosition().y();

    update();
}


void EcosimWidget::wheelEvent(QWheelEvent* e) 
{
    int numDegrees = e->angleDelta().y() / 8;
    int numSteps = numDegrees / 15;

    if (!(e->modifiers() & Qt::ShiftModifier)) {
        Vector3 direction = camera.getAt() - camera.getEye();
        double currentDist = Norm(direction);
        direction = Normalized(direction);
        double speed = 0.1 * currentDist;

        camera.setEye(camera.getEye() + direction * speed * numSteps);

        updateCameraPlanes();
    }

    update();
}


void EcosimWidget::callSimulationStep()
{
    timer.start();

    if (simvalid) {
        if (firstsim) {
            firstsim = false;
        }
        isSimWorking = true;
        emit _signalOperate();
    }
    else {
        std::cerr << "Simulation is invalidated, so not executed" << std::endl;
        isSimWorking = false;
        isSimRunning = false;
        doSimSteps = 0;
    }
}

void EcosimWidget::simulationStepFinished()
{
    isSimWorking = false;

    if (waitingForReset) {
        resetSim();
        return;
    }

    const double dt = 1e-3 * double(timer.elapsed()); // ms to s
    simPerf = dt;
    simTime += dt;

    simulatedMonths++;

    updatedMonthlyTexture = false;

    if (doSimSteps > 0) {
        doSimSteps--;
        isSimRunning = doSimSteps > 0;
    }

    // year simulated
    if (simulatedMonths % 12 == 0) {
        emit _signalFinishedYear();
    }

    // synchronize models between sim data and graphics
    // this call indirectly causes another paintGL()
    updatePlantInstances();
}



std::string EcosimWidget::get_dirprefix()
{
    // std::cout << "Datadir before fixing: " << datadir << std::endl;
    while (datadir.back() == '/')
        datadir.pop_back();

    // std::cout << "Datadir after fixing: " << datadir << std::endl;
    int slash_idx = int(datadir.find_last_of("/"));
    std::string setname = datadir.substr(slash_idx + 1);
    std::string dirprefix = datadir + "/" + setname;
    return dirprefix;
}

bool EcosimWidget::loadScene(const std::string& base_data_dir)
{
    this->makeCurrent();

    if (eco) delete eco;
    eco = new ecosim::EcoSystem();
    if (terrain) delete terrain;
    terrain = new ecosim::Terrain();
    if (biome) delete biome;
    biome = new ecosim::Biome(); 

    datadir = base_data_dir;

    simvalid = true;
    std::string terfile  = base_data_dir + "/heightfield.elv";
    std::string sunfile  = base_data_dir + "/sun.txt";
    std::string wetfile  = base_data_dir + "/wet.txt";
    std::string climfile = base_data_dir + "/clim.txt";
    std::string tempfile = base_data_dir + "/temp.txt";
    std::string catfile  = base_data_dir + "/plt.png";
    std::string almfile  = base_data_dir + "/alm";
    biome_db_filepath    = base_data_dir + "/biome.db";

    // load terrain
    if (!getTerrain()->loadElv(terfile)) return false;
    std::cerr << "Elevation file loaded" << std::endl;
    getTerrain()->calcMeanHeight(); 
    initTerrainGeometry(); CE();
    updateTexture(renderTexType, renderTexSubtype); CE();


    // load biome
    if (getBiome()->read_dataimporter(biome_db_filepath))
    {
        getBiome()->readAllometries(almfile);

        std::cerr << "Biome file load" << std::endl;

        if (sim != nullptr) {
            delete sim;
        }
        const int subcellFactor = 4;
        sim = new ecosim::Simulation(getTerrain(), getBiome(), subcellFactor);
        CE();

        selectedPlantId = -1;
        simulatedMonths = 0;
        doSimSteps = 0;
        simTime = 0;
        simPerf = 0;
        isSimRunning = false;
        updatedMonthlyTexture = false;

        // TEMP
        if (simWorker != nullptr) {
            delete simWorker;
        }
        simWorker = new SimWorker(sim);
        simWorker->moveToThread(&simThread);
        connect(this, &EcosimWidget::_signalOperate, simWorker, &SimWorker::simulationStep);
        connect(simWorker, &SimWorker::_finishedSimStep, this, &EcosimWidget::simulationStepFinished);
        simThread.start();
        // ----

        // read climate parameters
        if (!getSim()->readClimate(climfile))
        {
            simvalid = false;
            std::cerr << "No climate file " << climfile << " found. Simulation invalidated" << std::endl;
        }
        else // write temperature, to be used by grass simulator
        {
            getSim()->writeTemperature(tempfile, getTerrain());
        }

        // read soil moisture, and if that fails simulate it and store results
        if (!getSim()->readMoisture(wetfile))
        {
            std::cerr << "No Soil moisture file " << wetfile << " found, so simulating soil moisture" << std::endl;
            getSim()->calcMoisture();
            getSim()->writeMoisture(wetfile);
        }
        else
        {
            std::cerr << "Soil moisture file loaded" << std::endl;
        }

        // loading plant distribution
        getEcoSys()->setBiome(getBiome());
        
        std::string simsun_file = sunfile;

        // read sunlight, and if that fails simulate it and store results
        if (!getSim()->readSun(simsun_file))
        {
            std::cerr << "No Sunlight file " << simsun_file << " found, so (not) simulating sunlight" << std::endl;
            CE();
            getSim()->calcSunlight(glsun, 15, 64);  // direct light every 15min + 64 hemisphere samples as AO
            getSim()->writeSunlight(simsun_file);
        }
        else
        {
            std::cerr << "Sunlight file loaded" << std::endl;
        }

        getSim()->setEvents(ecosim::ScriptedDisturbanceEvent::loadFromFile(datadir + "/scriptedEvents.txt"));
        getSim()->setupSim(getEcoSys(), false);

        renderSpecies.resize(biome->numPFTypes());
        std::fill(renderSpecies.begin(), renderSpecies.end(), true);
    }
    else
    {
        std::cerr << "loadScene: Unable to load Biome from database." << endl;
        return false;
    }

    initPlantModels();

    if (simvalid)
        firstsim = true;
    return simvalid;
}


void EcosimWidget::initTerrainGeometry()
{
    // scene
    int nx, ny;
    float sizex, sizey;
    float zmin, zmax;    
    terrain->getGridDim(nx, ny);
    terrain->getTerrainDim(sizex, sizey);
    terrain->getHeightBounds(zmin, zmax);
    terrainBBox = Box3(Vector3(0, 0, zmin), Vector3(sizex, sizey, zmax));
    float dx = sizex / (nx - 1);
    float dy = sizey / (ny - 1);

    // verts
    int idx = 0;
    std::vector<GLfloat> verts(3 * nx * ny);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            verts[idx++] = GLfloat(i * dx);
            verts[idx++] = GLfloat(j * dy);
            verts[idx++] = GLfloat(terrain->getHeight(i, j));
        }
    }

    // tris
    numTerrainTriangles = (nx - 1) * (ny - 1) * 2;
    idx = 0;
    std::vector<GLuint> indices(numTerrainTriangles * 3);
    for (int i = 1; i < nx; i++) {
        for (int j = 1; j < ny; j++) {
            GLuint v00 = (i - 1) * ny + j - 1;
            GLuint v01 = (i - 1) * ny + j;
            GLuint v10 = i * ny + j - 1;
            GLuint v11 = i * ny + j;

            indices[idx++] = v00;
            indices[idx++] = v01;
            indices[idx++] = v10;

            indices[idx++] = v10;
            indices[idx++] = v01;
            indices[idx++] = v11;
        }
    }

    // update buffers
    CE();
    if (bufferVerts > 0) glDeleteBuffers(1, &bufferVerts); 
    if (bufferIndices > 0) glDeleteBuffers(1, &bufferIndices); 
    if (meshVAO > 0) glDeleteVertexArrays(1, &meshVAO);

    glGenVertexArrays(1, &meshVAO); CE();
    glBindVertexArray(meshVAO); CE();

    glUseProgram(shaderTerrain.programId()); CE();
    GLuint attribVertexLoc = glGetAttribLocation(shaderTerrain.programId(), "a_position"); CE();

    glGenBuffers(1, &bufferVerts); CE();
    glBindBuffer(GL_ARRAY_BUFFER, bufferVerts); CE();
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * verts.size(), &verts[0], GL_STATIC_DRAW); CE();
    glVertexAttribPointer(attribVertexLoc, 3, GL_FLOAT, GL_FALSE, 0, 0); CE();
    glEnableVertexAttribArray(attribVertexLoc); CE();

    glGenBuffers(1, &bufferIndices); CE();
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufferIndices); CE();
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * indices.size(), &indices[0], GL_STATIC_DRAW); CE();

    glBindVertexArray(0); CE();
    glUseProgram(0); CE();


    // light position
    double boxRad = terrainBBox.radius();
    lightPosition = terrainBBox.center() + 2 * boxRad * Vector3(0.5, 0.5, 1);        
}


void EcosimWidget::updateTexture(TerrainTexture t, int subtype)
{
    renderTexType = t;

    if (terrain == nullptr) return;

    if (subtype == -1) subtype = (simulatedMonths - 1) % 12;
    if (subtype < 0) subtype = 0;

    updatedMonthlyTexture = true;
    int nx, ny;
    terrain->getGridDim(nx, ny);

    texImg = QImage(nx, ny, QImage::Format_RGBA8888);

    switch (t) {
    case PLAIN:
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                Vector3 c(0.9,0.9,0.9);
                texImg.setPixelColor(i, j, toQColor(c));
            }
        }
        break;

    case ELEVATION: 
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                float h = terrain->getHeight(i, j);
                Vector3 c = ColorPalette::Relief().getColor((h-terrainBBox.getMin()[2]) / (terrainBBox.depth()));
                texImg.setPixelColor(i, j, toQColor(c));
            }
        }
        break;

    case NORMALS: 
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                ecosim::Vector n;
                terrain->getNormal(i, j, n);
                Vector3 n2 = 0.5 * Vector3(1 + n.i, 1 + n.k, 1 + n.j);
                texImg.setPixelColor(i, j, toQColor(n2));
            }
        }
        break;

    case ASPECT:
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                ecosim::vpPoint2D g = terrain->getGradient(i, j);
                float a = M_PI + std::atan2(g.x, g.y);
                float u = a / (2 * M_PI);
                Vector3 c = ColorPalette::CoolWarm().getColor(u);
                texImg.setPixelColor(i, j, toQColor(c));
            }
        }
        break;

    case ASPECT_E:
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                ecosim::vpPoint2D g = terrain->getGradient(i, j);
                float n = std::sqrt(g.x * g.x + g.y * g.y);
                float u = 0.5 - 0.5 * g.x / n;
                Vector3 c = ColorPalette::CoolWarm().getColor(u);
                texImg.setPixelColor(i, j, toQColor(c));
            }
        }
        break;

    case ASPECT_N:
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                ecosim::vpPoint2D g = terrain->getGradient(i, j);
                float n = std::sqrt(g.x * g.x + g.y * g.y);
                float u = 0.5 - 0.5 * g.y / n;
                Vector3 c = ColorPalette::CoolWarm().getColor(u);
                texImg.setPixelColor(i, j, toQColor(c));
            }
        }
        break;

    case SLOPE:
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                ecosim::vpPoint2D g = terrain->getGradient(i, j);
                float s = std::sqrt(g.x*g.x + g.y*g.y);
                Vector3 c = ColorPalette::CoolWarm().getColor(s/3);
                texImg.setPixelColor(i, j, toQColor(c));
            }
        }
        break;

    case SUNLIGHT: {
        MapFloat* sunlight = getSunlight(subtype);
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                float s = sunlight->get(i, j);
                Vector3 c = ColorPalette::CoolWarm().getColor((s - 4)/12);
                texImg.setPixelColor(i, j, toQColor(c));
            }
        }
    } break;

    case MOISTURE: {
        float vmin = FLT_MAX;
        float vmax = FLT_MIN;
        MapFloat* moisture = getMoisture(subtype);
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                float m = std::log(moisture->get(i, j));
                if (m > vmax) vmax = m;
                if (m < vmin) vmin = m;
            }
        }
        if (vmin < 0) vmin = 0;

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                float m = std::log(moisture->get(i, j));
                Vector3 c = ColorPalette::Blues().getColor((m - vmin) / (vmax - vmin));
                texImg.setPixelColor(i, j, toQColor(c));
            }
        }
    } break;

    case TEMPERATURE: {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                float t = sim->getTemperature(i, j, subtype);
                Vector3 c = ColorPalette::CoolWarm().getColor((t - 5) / (25 - 5));
                texImg.setPixelColor(i, j, toQColor(c));
            }
        }
    } break;

    case VIABILITY: {
        MapFloat* viability = sim->getViabilityMap(subtype);
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                float v = viability->get(i, j);
                Vector3 c = ColorPalette::CoolWarm().getColor((v + 1) / 2);
                texImg.setPixelColor(i, j, toQColor(c));
            }
        }
        break;
    }

    case MOSTVIABLE:
        if (!viabilitiesTexComputed) {
            viabilitiesTex = QImage(nx, ny, QImage::Format_RGBA8888);
            std::vector<MapFloat*> viabilities(getBiome()->numPFTypes());
            for (int k = 0; k < getBiome()->numPFTypes(); ++k) {
                viabilities[k] = sim->getViabilityMap(k);
            }
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    float maxV = 0;
                    int maxK = 0;
                    for (int k = 0; k < getBiome()->numPFTypes(); ++k) {
                        float v = viabilities[k]->get(i, j);
                        if (v > maxV) {
                            maxV = v;
                            maxK = k;
                        }
                    }
                    viabilitiesTex.setPixelColor(i, j, getColorForSpecies(maxK));
                }
            }
            viabilitiesTexComputed = true;
        }
        texImg = viabilitiesTex;
        break;
    }

    updateTexture();
}

void EcosimWidget::updateTexture()
{
    QImage tex = texImg.copy();
    if (renderShadeTerrain) {
        for (int i = 0; i < tex.width(); i++) {
            for (int j = 0; j < tex.height(); j++) {
                QColor qc = tex.pixelColor(i, j);
                Vector3 cz = Vector3(qc.redF(), qc.greenF(), qc.blueF());

                ecosim::Vector n;
                terrain->getNormal(i, j, n);
                Vector3 normal(n.i, n.k, n.j);
                double s = normal * Normalized(Vector3(-2.0, 1.0, 4.0));
                s = 0.5 * (1.0 + s);
                s *= s;
                Vector3 c = 0.2 * Vector3(1.0, 1.0, 1.0) + 0.6 * s * cz + 0.2 * s * Vector3(1.0, 1.0, 1.0);

                tex.setPixelColor(i, j, toQColor(c));
            }
        }
    }

    // update texture
    glBindTexture(GL_TEXTURE_2D, texId);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, tex.width(), tex.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, tex.bits());
    glGenerateMipmap(GL_TEXTURE_2D);
}


GLuint createVAO(Model& model, std::vector<GLuint>& buffersList, const GLuint instanceBuffer)
{
    GLuint vao, vboVertices, vboNormals, ebo; 

    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    glGenBuffers(1, &vboVertices);
    glBindBuffer(GL_ARRAY_BUFFER, vboVertices);
    glBufferData(GL_ARRAY_BUFFER, model.getVertexCoords().size() * sizeof(float), model.getCoordsPtr(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0); 
    glEnableVertexAttribArray(0);

    glGenBuffers(1, &vboNormals);
    glBindBuffer(GL_ARRAY_BUFFER, vboNormals);
    glBufferData(GL_ARRAY_BUFFER, model.getNormals().size() * sizeof(float), model.getNormalsPtr(), GL_STATIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
    glEnableVertexAttribArray(1);

    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, instanceBuffer);
    glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 20 * sizeof(GLfloat), 0);
    glVertexAttribDivisor(2, 1);
    for (int i = 0; i < 4; i++) {
      glEnableVertexAttribArray(3 + i);
      glVertexAttribPointer(3 + i, 4, GL_FLOAT, GL_FALSE, 20 * sizeof(GLfloat), (void*)((i + 1) * 4 * sizeof(GLfloat)));
      glVertexAttribDivisor(3 + i, 1);
    }

    glGenBuffers(1, &ebo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, model.getIndices().size() * sizeof(unsigned int), model.getIndicesPtr(), GL_STATIC_DRAW);
    
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    buffersList.push_back(vboVertices);
    buffersList.push_back(vboNormals);
    buffersList.push_back(ebo);
    return vao;
}


void EcosimWidget::initPlantModels()
{
    deletePlantModelBuffers();

    // instance data buffer
    glGenBuffers(1, &instanceVBO);
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glBufferData(GL_ARRAY_BUFFER, 10000000 * sizeof(PlantInstanceData), nullptr, GL_DYNAMIC_DRAW);


    Model trunk = Model::createCylinder(6, false);
    trunkVAO = createVAO(trunk, pftBuffers, instanceVBO);
    trunkVAOSize = trunk.numFaces() * 3;

    for (int i = 0; i < sim->getBiome()->numPFTypes(); i++) {
        const ecosim::PFType* pft = sim->getBiome()->getPFType(i);
        if (pft->isTree) {
            Model canopy = pft->isConifer ? Model::createCone(16) : Model::createIcosphere(1);
            pftVAO.push_back(createVAO(canopy, pftBuffers, instanceVBO));
            pftVAOSize.push_back(canopy.numFaces()*3);
        }
        else {
            Model canopy = Model::createCube();
            pftVAO.push_back(createVAO(canopy, pftBuffers, instanceVBO));
            pftVAOSize.push_back(canopy.numFaces() * 3);
        }
    }


}

void EcosimWidget::deletePlantModelBuffers()
{
    for (GLuint bufId : pftBuffers)
        glDeleteBuffers(1, &bufId);
    for (GLuint vaoId : pftVAO)
        glDeleteVertexArrays(1, &vaoId);
    if (trunkVAO > 0)
        glDeleteVertexArrays(1, &trunkVAO);
    if (instanceVBO > 0)
        glDeleteBuffers(1, &instanceVBO);
}

bool EcosimWidget::filtersOk(ecosim::SimLivePlant* plt)
{
    return minHeight <= plt->height && plt->height <= maxHeight
        && minAge <= plt->age && plt->age <= maxAge;
}

void EcosimWidget::setFilters(float min_h, float max_h, int min_age, int max_age)
{
    minHeight = min_h;
    maxHeight = (max_h >= 0) ? max_h : FLT_MAX;

    minAge = min_age;
    maxAge = (max_age >= 0) ? max_age : INT32_MAX;
    updatePlantInstances();
}


EcosimWidget::PlantInstanceData::PlantInstanceData(const QColor& c, const QMatrix4x4& mat)
{     
    this->color[0] = c.red() / 255.0f;
    this->color[1] = c.green() / 255.0f;
    this->color[2] = c.blue() / 255.0f;
    this->color[3] = c.alpha() / 255.0f;
    for (int i = 0; i < 16; i++)
        this->transform[i] = mat.constData()[i];
}

QColor EcosimWidget::encodeId(int type, int id) const
{
    unsigned char r = (id >> 24) & 0xFF | (type << 6);
    unsigned char g = (id >> 16) & 0xFF;
    unsigned char b = (id >> 8) & 0xFF;
    unsigned char a = id & 0xFF;
    return QColor(r, g, b, a);
}

unsigned int EcosimWidget::decodeId(const QColor& c, int& type) const
{
    type = (c.red() & 0xC0) >> 6;
    return ((c.red() & 0x3F) << 24) | (c.green() << 16) | (c.blue() << 8) | c.alpha();
}

inline void transformsConifer(double x, double y, double z, float h, float canopy, float rad, QMatrix4x4& matC, QMatrix4x4& matT)
{
    const double baseline = 0.25;
    matC.setToIdentity();
    matC.translate(x, y, z + h * baseline);
    matC.scale(canopy, canopy, h * (1 - baseline));
    matC.rotate(90, 1, 0, 0);

    matT.setToIdentity();
    matT.translate(x, y, z + h * baseline / 2.0);
    matT.scale(rad, rad, 1.1 * h * baseline / 2.0);
    matT.rotate(90, 1, 0, 0);
}

inline void transformsBroadleaf(double x, double y, double z, float h, float canopy, float rad, QMatrix4x4& matC, QMatrix4x4& matT)
{
    const double baseline = 1.0 / 3.0;
    matC.setToIdentity();
    matC.translate(x, y, z + 0.5 * (baseline + 1) * h);
    matC.scale(canopy, canopy, h * (1 - baseline) / 2.0);
  
    matT.setToIdentity();
    matT.translate(x, y, z + h * baseline / 2.0);
    matT.scale(rad, rad, 1.1 * h * baseline / 2.0);
    matT.rotate(90, 1, 0, 0);
}

inline void transformNontree(double x, double y, double z, float h, float canopy, QMatrix4x4& mat)
{
    mat.setToIdentity();
    mat.translate(x, y, z + h / 2.0);
    mat.scale(canopy, canopy, h / 2.0);
}

void EcosimWidget::updatePlantInstances(bool falseColor)
{
    if (terrain == nullptr)
        return;

    if (falseColor != paintingFalseColor) return;

    instancesCanopy.clear();
    instancesCanopy.resize(sim->getBiome()->numPFTypes());
    instancesTrunks.clear();
    instancesLogs.clear();

    QColor selectionColor(255, 92, 0, 255);
    bool foundSelected = false;

    QMatrix4x4 matC, matT;

    if (this->renderPlts) {

        const std::vector<ecosim::SimLivePlant*>& livePlants = sim->getLivePlants();
        for (int i = 0; i < livePlants.size(); i++) {
            Vector3 vc = ecosimToWidgetCoords(livePlants[i]->pos);
            const double x = vc[0];
            const double y = vc[1];
            const double z = vc[2];
            const float h = livePlants[i]->height;
            const float canopy = livePlants[i]->canopy;
            const float rad = livePlants[i]->dbh / 2;
            const float vigour = livePlants[i]->vigourhistory.size() > 0 ? livePlants[i]->vigourhistory.back() : 0;
            const int pft = livePlants[i]->pft;
            const ecosim::PFType* pftype = getBiome()->getPFType(pft);

            if (livePlants[i]->state != ecosim::PlantSimState::ALIVE) continue;
            if (!renderSpecies[pft]) continue;
            if (!filtersOk(livePlants[i])) continue;

            QColor color;
            if (plantColorType == 0) // species
                color = toQColor(Vector3(pftype->basecol[0], pftype->basecol[1], pftype->basecol[2]));
            else if (plantColorType == 1) // vigour
                color = toQColor(ColorPalette::GreenGreyBrown().getColor(0.5*vigour + 0.5));
            if (selectedPlantType == 1 && livePlants[i]->id == selectedPlantId) {
                color = selectionColor;
                loadPlantInfo(livePlants[i]);
                foundSelected = true;
            }
            if (falseColor) color = encodeId(1, i);

            if (pftype->isTree) {
                if (pftype->isConifer)
                    transformsConifer(x, y, z, h, canopy, rad, matC, matT);
                else
                    transformsBroadleaf(x, y, z, h, canopy, rad, matC, matT);
                instancesCanopy[pft].push_back(PlantInstanceData(color, matC));
                instancesTrunks.push_back(PlantInstanceData(color, matT));
            }
            else {
                transformNontree(x, y, z, h, canopy, matC);
                instancesCanopy[pft].push_back(PlantInstanceData(color, matC));
            }
        }
    }


    // dead plants
    if (this->renderSnags) {
        const std::vector<ecosim::SimDeadPlant*>& deadPlants = sim->getDeadPlants();
        for (int i = 0; i < deadPlants.size(); i++) {
            Vector3 vc = ecosimToWidgetCoords(deadPlants[i]->parent->pos);
            const double x = vc[0];
            const double y = vc[1];
            const double z = vc[2];
            const double canopy = deadPlants[i]->parent->canopy;            
            const int pft = deadPlants[i]->parent->pft;
            const float rad = deadPlants[i]->parent->dbh / 2;
            const float decay = deadPlants[i]->decay;
            const ecosim::PFType* pftype = getBiome()->getPFType(pft);
            const ecosim::SnagDecayClass decayClass = deadPlants[i]->dc;

            if (!renderSpecies[pft]) continue;
            if (not filtersOk(deadPlants[i]->parent)) continue;

            QColor color;
            if (plantColorType == 0) // species
                color = toQColor(0.5*Vector3(pftype->basecol[0], pftype->basecol[1], pftype->basecol[2]));
            else if (plantColorType == 1) // vigour
                color = toQColor(ColorPalette::GreenGreyBrown().getColor(0.5*decay));
            if (deadPlants[i]->eventhistory.size() > 0 && deadPlants[i]->eventhistory.back().cat == ecosim::DisturbanceEventCategory::FIRE ||
                deadPlants[i]->parent->eventhistory.size() > 0 && deadPlants[i]->parent->eventhistory.back().cat == ecosim::DisturbanceEventCategory::FIRE)
            {
                color = QColor(0, 0, 0);
            }
            else {
                if (decayClass == ecosim::SnagDecayClass::BARE) {
                    float tt = std::min(1.0f, (def_snagDC1 - decay) / (def_snagDC1 - def_snagDC2));
                    color = toQColor((1 - tt) * fromQColor(color) + tt * Vector3(0.5, 0.5, 0.5));
                }
            }
            if (selectedPlantType == 2 && deadPlants[i]->id == selectedPlantId) {
                color = selectionColor;
                loadPlantInfo(deadPlants[i]);
                foundSelected = true;
            }
            if (falseColor) color = encodeId(2, i);

            double h = deadPlants[i]->parent->height;
            if (decayClass > ecosim::SnagDecayClass::BARE) h *= deadPlants[i]->trunkremains;

            if (pftype->isTree) {
                if (pftype->isConifer)
                    transformsConifer(x, y, z, h, canopy, rad, matC, matT);
                else
                    transformsBroadleaf(x, y, z, h, canopy, rad, matC, matT);

                if (decayClass < ecosim::SnagDecayClass::DEBRANCHED) {
                    instancesCanopy[pft].push_back(PlantInstanceData(color, matC));
                }
                instancesTrunks.push_back(PlantInstanceData(color, matT));
            }
            else {
                transformNontree(x, y, z, h, canopy, matC);
                instancesCanopy[pft].push_back(PlantInstanceData(color, matC));
            }
        }
    }


    if (this->renderLogs) {
        const std::vector<ecosim::SimLog*>& logs = sim->getLogs();    
        for (int i = 0; i < logs.size(); i++) {
            const double x = logs[i]->liveparent->pos.x;
            const double y = logs[i]->liveparent->pos.z;
            const double z = terrain->getHeightAtCoords(x, y);
            const double rad = logs[i]->diam/2;
            const double h = logs[i]->len;
            const double dirn = qRadiansToDegrees(logs[i]->dirn);

            QColor color(52, 33, 0);
            if (selectedPlantType == 3 && logs[i]->id == selectedPlantId) {
              color = selectionColor;
              loadPlantInfo(logs[i]);
              foundSelected = true;
            }
            if (falseColor) color = encodeId(3, i);

            int ti, tj;
            terrain->toGrid(logs[i]->liveparent->pos, ti, tj);
            ecosim::vpPoint2D grad = terrain->getGradient(ti, tj);
            Vector3 terrNormal = Normalized(Vector3(-grad.x, -grad.y, 1.0));
            Vector3 terrTangent = cross(Vector3(0, 0, 1), terrNormal);
            Vector3 slopeDir;
            if (Norm(terrTangent) < 1e-3) slopeDir = Vector3(0, 1, 0);
            else slopeDir = Normalized(cross(terrNormal, Normalized(terrTangent)));
            Vector3 rotAxis = Normalized(cross(Vector3(0, 1, 0), slopeDir));
            double rotAngle = qRadiansToDegrees(std::acos(slopeDir[1]));

            matT.setToIdentity();
            matT.translate(x, y, z + rad);
            matT.rotate(dirn, terrNormal[0], terrNormal[1], terrNormal[2]);
            matT.rotate(rotAngle, rotAxis[0], rotAxis[1], rotAxis[2]);
            matT.scale(rad, h / 2.0, rad);

            instancesTrunks.push_back(PlantInstanceData(color, matT));
        }
    }

    if (!foundSelected) {
      selectedPlantId = -1;
      selectedPlantType = 0;
      emit _signalTextInfo("");
    }

    update();
}

bool EcosimWidget::getRayTerrainIntersection(const Ray& ray, QPoint& hitCell) const
{
    ecosim::vpPoint hitp;
    bool intersection = terrain->rayIntersect(ecosim::vpPoint(ray.origin()[0], ray.origin()[2], ray.origin()[1]),
                                              ecosim::Vector(ray.direction()[0], ray.direction()[2], ray.direction()[1]),
                                              hitp);
    if (intersection) {
        int hx, hy;
        terrain->toGrid(hitp, hx, hy);
        hitCell.setX(hx);
        hitCell.setY(hy);
        return true;
    }
    return false;
}
