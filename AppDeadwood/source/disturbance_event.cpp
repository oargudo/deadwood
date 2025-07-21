#include "disturbance_event.h"
#include <fstream>
#include <algorithm>
#include <set>
//#include <QtGui/QImage>

using namespace ecosim;

int ScriptedDisturbanceEvent::diseaseNum = 0;

DiceRoller DisturbanceAction::dice = DiceRoller(0, 100000);

float DisturbanceAction::randomUniform() {
    return dice.generate() / 100000.0;
}

float DisturbanceAction::randomUniform(float rmin, float rmax) {
    return randomUniform() * (rmax - rmin) + rmin;
}

std::vector<ScriptedDisturbanceEvent*> ScriptedDisturbanceEvent::loadFromFile(const std::string& filepath)
{
    std::vector<ScriptedDisturbanceEvent*> v;

    std::fstream fin(filepath, std::fstream::in);
    if (fin.good()) {

        int year, month, evtType;
        while (fin >> year >> month >> evtType) {
            
            ScriptedDisturbanceEvent* evt = new ScriptedDisturbanceEvent();
            evt->year = year;
            evt->month = month;

            switch (evtType) {
            case 1: { // wildfire
                Wildfire* devt = new Wildfire();

                float windx, windy;
                fin >> windx >> windy;
                devt->setWind(windx, windy);

                int numOrigins;
                fin >> numOrigins;
                for (int i = 0; i < numOrigins; i++) {
                    float x, y, r;
                    fin >> x >> y >> r;
                    devt->addOrigin(x, y, r);
                }

                evt->distEvent = devt;
                std::cout << "Loaded wildfire event on year/month = " << year << "/" << month << std::endl;
            }
            break;
            case 2: { // windstorm
                Windstorm* devt = new Windstorm();

                float windx, windy;
                fin >> windx >> windy;
                devt->setWind(windx, windy);

                evt->distEvent = devt;
                std::cout << "Loaded windstorm event on year/month = " << year << "/" << month << std::endl;
            }
            break;
            case 3: { // disease
                DiseaseOutbreak* devt = new DiseaseOutbreak();
                float ox, oy;
                fin >> ox >> oy;
                devt->setOutbreakLocation(ox, oy);
                
                Disease* params = new Disease();
                params->diseaseId = ++diseaseNum;
                fin >> params->targetSpecies;
                fin >> params->spreadRadius;
                fin >> params->spreadProbability;
                fin >> params->severity;
                fin >> params->monthsContagious;
                fin >> params->monthsRecovery;
                fin >> params->monthsImmune;
                params->totalInfected = 0;
                params->currentInfected = 0;
                devt->setDiseaseParams(params);

                evt->distEvent = devt;
                std::cout << "Loaded disease event on year/month = " << year << "/" << month << std::endl;
            }
            break;
            case 4: { // drought
                Drought* devt = new Drought();

                int duration;
                fin >> duration;
                std::vector<float> vals(duration);
                for (int i = 0; i < duration; i++) {
                    fin >> vals[i];
                }
                devt->setDroughtValues(vals);

                evt->distEvent = devt;
                std::cout << "Loaded drought event on year/month = " << year << "/" << month << std::endl;
            }
            break;
            default:
                std::cerr << "Unknown event type" << std::endl;
                continue;
            }

            v.push_back(evt);
        }
    }
    else {
        std::cout << "Could not find disturbance events file: " << filepath << std::endl;
    }
    fin.close();

    // sort by date 
    std::sort(v.begin(), v.end(), 
             [](ScriptedDisturbanceEvent* a, ScriptedDisturbanceEvent* b) {return *a < *b;});

    std::vector<ScriptedDisturbanceEvent*> q;
    for (int i = 0; i < v.size(); i++) {
        q.push_back(v[i]);
        //std::cerr << "   EVENT IN QUEUE: " << v[i]->year << " " << v[i]->month << std::endl;
    }

    return q;
}

bool ecosim::ScriptedDisturbanceEvent::saveToFile(const std::vector<ScriptedDisturbanceEvent*>& events, const std::string& filepath)
{
    std::fstream fout(filepath, std::fstream::out);
    if (fout.good()) {

        int year, month;
        for (int k = 0; k < events.size(); ++k) {

            ScriptedDisturbanceEvent* evt = events[k];
            year = evt->year;
            month = evt->month;

            fout << year << " " << month;
            fout << " " << static_cast<int>(evt->distEvent->type()) + 1;
            //cout << year << " " << month << " " << static_cast<int>(evt->distEvent->type()) + 1 << endl;

            switch (evt->distEvent->type()) {
            case DisturbanceEventCategory::FIRE: { // wildfire
                Wildfire* wfev = dynamic_cast<Wildfire*>(evt->distEvent);

                float windX, windY;
                wfev->getWind(windX, windY);
                fout << " " << windX << " " << windY;

                std::vector<float> origins = wfev->getOrigins();
                fout << " " << origins.size() / 3;
                for (int i = 0; i < origins.size(); i += 3) {
                    fout << " " << origins[i] << " " << origins[i+1] << " " << origins[i+2];
                }
            }
            break;
            case DisturbanceEventCategory::WIND: { // windstorm
                Windstorm* wsev = dynamic_cast<Windstorm*>(evt->distEvent);

                float windX, windY;
                wsev->getWind(windX, windY);
                fout << " " << windX << " " << windY;
            }
            break;
            case DisturbanceEventCategory::DISEASE: { // disease

                DiseaseOutbreak* dsev = dynamic_cast<DiseaseOutbreak*>(evt->distEvent);

                float x, y;
                dsev->getOutbreakLocation(x, y);
                fout << " " << x << " " << y;

                const Disease* dparams = dsev->getDiseaseParams();
                fout << " " << dparams->targetSpecies << " " << dparams->spreadRadius << " " << dparams->spreadProbability <<
                    " " << dparams->severity << " " << dparams->monthsContagious << " " << dparams->monthsRecovery << " " << dparams->monthsImmune;
            }
            break;
            case DisturbanceEventCategory::DROUGHT: { // drought
                Drought* drev = dynamic_cast<Drought*>(evt->distEvent);

                std::vector<float> drValues = drev->getDroughtValues();
                fout << " " << drValues.size();
                for (int i = 0; i < drValues.size(); ++i) {
                    fout << " " << drValues[i];
                }
            }
            break;
            default:
                std::cerr << "Unknown event type" << std::endl;
                continue;
            }

            fout << endl;
        }
        fout.close();
        return true;
    }
    else {
        std::cout << "Could not find disturbance events file: " << filepath << std::endl;
        return false;
    }
}



/****************************************************************************************
 *                                       WILDFIRE                                       *
 ****************************************************************************************/

int computeScaleFactor(int numCells, float cellSize, float desiredCellSize) {
    // fast and dirty, there's probably a more clever method
    int scale = 1;
    int k = 2;
    while (cellSize*scale < desiredCellSize) {
        if (numCells%k == 0) {
            numCells /= k;
            scale *= k;
        }
        else {
            k++;
        }
    }
    return scale;
}


void Wildfire::execute(Simulation* sim)
{
    // read necessary info about simulation   
    Terrain* terrain = sim->getTerrain();
    float terDimX, terDimY;
    terrain->getTerrainDim(terDimX, terDimY);

    MapFloat simHeight, simDensity, simMoisture;
    sim->deriveFireMaps(simHeight, simDensity, simMoisture);

    // compute number of cells to make them ~1m side
    int ncX, ncY;
    simDensity.getDim(ncX, ncY);
    
    int scaleFactorX = computeScaleFactor(ncX, terDimX/ncX, 1.0);
    int scaleFactorY = computeScaleFactor(ncY, terDimY/ncY, 1.0);
    numCellsX = ncX / scaleFactorX;
    numCellsY = ncY / scaleFactorY;
    dimCellX = terDimX / numCellsX;
    dimCellY = terDimY / numCellsY;

    cells.resize(numCellsX*numCellsY);

    std::cout << "FIRE: Terrain dimensions " << terDimX << " " << terDimY << std::endl;
    std::cout << "FIRE: Maps dimensions " << ncX << " " << ncY << std::endl;
    std::cout << "FIRE: Scale factors " << scaleFactorX << " " << scaleFactorY << std::endl;
    std::cout << "FIRE: Num cells " << numCellsX << " " << numCellsY << std::endl;
    std::cout << "FIRE: Cell resolution " << dimCellX << " " << dimCellY << std::endl;


    // redimension maps by average filtering
    MapFloat height, density, moisture;
    height.setDim(numCellsX, numCellsY);
    density.setDim(numCellsX, numCellsY);
    moisture.setDim(numCellsX, numCellsY);
    for (int i = 0; i < numCellsX; i++) {
        for (int j = 0; j < numCellsY; j++) {
            int n = 0;
            float hsum = 0;
            float dsum = 0;
            float msum = 0;
            for (int ii = 0; ii < scaleFactorX; ii++) {
                for (int jj = 0; jj < scaleFactorY; jj++) {
                    hsum += simHeight.get(i*scaleFactorX + ii, j*scaleFactorY + jj);
                    dsum += simDensity.get(i*scaleFactorX + ii, j*scaleFactorY + jj);
                    msum += simMoisture.get(i*scaleFactorX + ii, j*scaleFactorY + jj);
                    n++;
                }
            }
            height.set(i, j, hsum/n);
            density.set(i, j, dsum/n);
            moisture.set(i, j, msum/n);
        }
    }

    // initialize cell grid as flammable where there is vegetation
    // use this loop at the same time to compute avg and std of density map
    float densSum = 0;
    float densSqSum = 0;
    int densCnt = 0;
    for (int i = 0; i < numCellsX; i++) {
        for (int j = 0; j < numCellsY; j++) {
            float d = density.get(i, j);
            cell(i, j) = d > 0 ? FireState::UNBURNED : FireState::NONFLAMMABLE;
            if (d > 0) {
                densSum += d;
                densSqSum += d*d;
                densCnt++;
            }
        }
    }

    // remap densities to 0..1
    float densMean = densSum/densCnt;
    float densStd = std::sqrt(densSqSum/densCnt - densMean*densMean);
    float densMin = std::max(0.0f, densMean - densStd);
    float densMax = densMean + densStd;
    float densRng = densMax - densMin;
    for (int i = 0; i < numCellsX; i++) {
        for (int j = 0; j < numCellsY; j++) {
            float d = density.get(i, j);
            d = (d - densMin)/densRng;
            d = std::max(0.0f, std::min(1.0f, d));
            density.set(i, j, d);
        }
    }

    // blur maps
    double blurSigma = 6;
    height.gaussianBlur(6);
    density.gaussianBlur(6);
    moisture.gaussianBlur(6);

    
    /*QImage imgH(numCellsX, numCellsY, QImage::Format_RGBA8888);
    QImage imgD(numCellsX, numCellsY, QImage::Format_RGBA8888);
    QImage imgM(numCellsX, numCellsY, QImage::Format_RGBA8888);
    for (int i = 0; i < numCellsX; i++) {
        for (int j = 0; j < numCellsY; j++) {
            int cH = int(255 * std::min(1.0f, height.get(i, j) / 10.0f));
            int cM = int(255 * moisture.get(i, j));
            int cD = int(255 * density.get(i, j));
            imgH.setPixelColor(i, j, QColor(cH, cH, cH));
            imgM.setPixelColor(i, j, QColor(cM, cM, cM));
            imgD.setPixelColor(i, j, QColor(cD, cD, cD));
        }
    }
    imgH.mirrored(false, true).save("fire_" + QString::number(sim->year) + "_heights.png");
    imgD.mirrored(false, true).save("fire_" + QString::number(sim->year) + "_density.png");
    imgM.mirrored(false, true).save("fire_" + QString::number(sim->year) + "_moisture.png");*/
    

 
    // initial fire positions
    if (fireOrigin.size() > 0) {
        for (int k = 0; k < int(fireOrigin.size()); k += 3) {
            int i = int(fireOrigin[k] * numCellsX);
            int j = int(fireOrigin[k + 1] * numCellsY);
            float r = fireOrigin[k + 2];

            int dri = int(r/dimCellX) + 1;
            int drj = int(r/dimCellY) + 1;
            for (int di = -dri; di <= dri; di++) {
                for (int dj = -drj; dj <= drj; dj++) {
                    if (i + di < 0 || i + di >= numCellsX) continue;
                    if (j + dj < 0 || j + dj >= numCellsX) continue;

                    float d2 = di*di*dimCellX*dimCellX + dj*dj*dimCellY*dimCellY;
                    float pf = 1.0 - d2/(r*r);
                    float u = randomUniform();
                    if (u*u < pf) {
                        cell(i + di, j + dj) = FireState::BURNING;                        
                    }
                }
            }
            cell(i, j) = FireState::BURNING;
        }
    }
    else {
        int i, j;
        do {
            i = randomUniform(0, numCellsX-1);
            j = randomUniform(0, numCellsY-1);
        } while (cell(i, j) == FireState::NONFLAMMABLE);
        cell(i, j) = FireState::BURNING;
    }


    // simulation
    int numSteps = 0;
    int burnedCells = 0;
    int totalBurnedCells = 0;
    do {
        numSteps++;
        burnedCells = propagationStep(terrain, height, density, moisture);
        totalBurnedCells += burnedCells;
        //std::cerr << "FIRE: propagation step " << numSteps << " burned " << burnedCells << " cells" << std::endl;

        // DEBUG
        //saveDebugImage("fire_" + QString::number(sim->year) + "_" + QString::number(numSteps) + QString(".png"));
    } while (burnedCells > 0);

    std::cout << "Fire event burned " << totalBurnedCells << " cells = " 
      << 0.0001 * totalBurnedCells * dimCellX * dimCellY << " ha" << std::endl;


    // modify ecosystem according to burned cells
    MapFloat fireMap;
    fireMap.setDim(numCellsX, numCellsY);
    for (int i = 0; i < numCellsX; i++) {
        for (int j = 0; j < numCellsY; j++) {
            fireMap.set(i, j, cell(i, j) == FireState::BURNED ? 1.0f : 0.0f);
        }
    }
    sim->applyFireMap(fireMap);
}

int Wildfire::propagationStep(Terrain* terrain, const MapFloat& vegHeight, const MapFloat& vegDensity, const MapFloat& vegMoisture)
{
    // constants: see Table 4 in [Alexandridis et al. 2011]
    const double p0 = 0.6;
    const double moist_a = 3.258;
    const double moist_b = 0.111;
    const double wind_c1 = 0.045;
    const double wind_c2 = 0.191;
    const double height_d = 0.932;
    const double slope_a = 0.063;

    int burnedCells = 0;
    std::vector<FireState> next(cells);

    for (int i = 0; i < numCellsX; i++) {
        for (int j = 0; j < numCellsY; j++) {

            if (cell(i, j) == FireState::BURNING) {

                // cells burn after one step
                next[index(i, j)] = FireState::BURNED;
                burnedCells++;

                // burning cells can propagate fire to neighbors
                for (int di = -1; di <= 1; di++) {
                    for (int dj = -1; dj <= 1; dj++) {
                        if (di == 0 && dj == 0) continue;
                        if (!insideIndex(i + di, j + dj)) continue;

                        int ni = i + di;
                        int nj = j + dj;

                        if (next[index(ni, nj)] == FireState::UNBURNED) {

                            // get origin and neighbor points
                            vpPoint pi(i * dimCellX, 0, j * dimCellY);
                            vpPoint pn(ni * dimCellY, 0, nj * dimCellY);
                            int tii, tij, tni, tnj;
                            terrain->toGrid(pi, tii, tij);
                            terrain->toGrid(pn, tni, tnj);
                            pi.z = terrain->getHeight(tii, tij);
                            pn.z = terrain->getHeight(tni, tni);                            
                            float dist = std::sqrt(di*di*dimCellX*dimCellX + dj*dj*dimCellY*dimCellY);
                            float propDirX = di / dist;
                            float propDirY = dj / dist;

                            // density coefficient (sparse -0.4, normal 0, dense 0.3)
                            float d = vegDensity.get(ni, nj);
                            double pDen = 1 - 0.4 + d*0.7;

                            // moisture term 
                            float m = vegMoisture.get(ni, nj);
                            double pMoist = moist_a * std::exp(-moist_b * m);

                            // height term
                            float h = vegHeight.get(ni, nj);
                            double pHeight = std::pow(h, height_d);

                            // slope term
                            float slope = std::atan((pn.z - pi.z) / dist);
                            double pSlope = std::exp(slope_a * slope);

                            // wind term
                            float cosWindProp = windDirX * propDirX + windDirY * propDirY;
                            double pWind = std::exp(windSpeed * (wind_c1 + wind_c2 * (cosWindProp - 1)));

                            // merge them together to obtain burn probability
                            double pBurn = p0 * pDen * pMoist * pHeight * pSlope * pWind;

                            /*
                            std::cerr << "FIRE SPREAD PROB FROM " << i << " " << j << " -> " << ni << " " << nj 
                                << " = " << pBurn << std::endl;
                            std::cerr << "    p0      = " << p0 << std::endl;
                            std::cerr << "    pDen    = " << pDen << std::endl;
                            std::cerr << "    pMoist  = " << pMoist << std::endl;
                            std::cerr << "    pHeight = " << pHeight << " (" << h << ")" << std::endl;
                            std::cerr << "    pSlope  = " << pSlope << " (" << slope << ")" << std::endl;
                            std::cerr << "    pWind   = " << pWind << " (" << cosWindProp << ")" << std::endl;
                            */

                            // do we propagate?
                            if (randomUniform() < pBurn) {
                                next[index(ni, nj)] = FireState::BURNING;
                            }
                        }
                    }
                }

                // random spotting to a further location
                // TODO
            }
        }
    }

    cells = next;

    return burnedCells;
}

void Wildfire::saveDebugImage(const QString& imgName) const
{
    QImage img = QImage(numCellsX, numCellsY, QImage::Format_RGBA8888);
    for (int i = 0; i < numCellsX; i++) {
        for (int j = 0; j < numCellsY; j++) {
            switch (cell(i, j)) {
            case FireState::NONFLAMMABLE:
                img.setPixelColor(i, j, QColor(64, 128, 255)); break;
            case FireState::UNBURNED:
                img.setPixelColor(i, j, QColor(0, 128, 64)); break;
            case FireState::BURNING:
                img.setPixelColor(i, j, QColor(255, 128, 64)); break;
            case FireState::BURNED:
                img.setPixelColor(i, j, QColor(64, 64, 64)); break;
            default:
                img.setPixelColor(i, j, QColor(255, 255, 255)); break;
            }
        }
    }
    img.mirrored(false, true).save(imgName);
}




/****************************************************************************************
 *                                       WINDSTORM                                      *
 ****************************************************************************************/

/*
float getCreg(int type) {
    switch (type) {
        case 3: return 162.0; // pinus uncinata (using pinus nigra coefficient)
        case 4: return 145.0; // abies alba
        case 5: return 162.0; // betula pendula
        case 6: return 162.0; // quercus petraea
        case 7: return 162.0; // fagus sylvatica
        default: return 162.0; // seems to be pretty common in iLand DB
    }
}

float getMOR(int type) {
    switch (type) {
        case 3: return 44.0; // pinus uncinata (using pinus nigra coefficient)
        case 4: return 42.0; // abies alba
        case 5: return 63.0; // betula pendula
        case 6: return 59.0; // quercus petraea
        case 7: return 65.0; // fagus sylvatica
        default: return 45.0; // something in between?
    }
}

float getWetBiomassFactor(int type) {
    return 1.85; // iLand DB used same value everywhere, but I guess could be per-species
}

// values from: https://www.wood-database.com
float getWoodDensity(int type) {  // kg/m3
    switch (type) {
        case 3: return 505.0; // pinus uncinata: from https://agritrop.cirad.fr/581489/1/2016%2008%20WCTE%20Eduard%20Pinus%20uncinata.pdf
        case 4: return 415.0; // abies alba
        case 5: return 640.0; // betula pendula
        case 6: return 710.0; // quercus petraea
        case 7: return 710.0; // fagus sylvatica
        default: return 600.0; // something average
    }
}

float getStemMass(int type, float dbh) 
{
    float a = 1;
    float b = 1;
    switch (type) 
    {
        case 3: a = 0.035558;    b = 2.71865;    break; // pinus uncinata (using pinus nigra coefficient)
        case 4: a = 0.03007189;  b = 2.74001367; break; // abies alba
        case 5: a = 0.114223297; b = 2.2919;     break; // betula pendula
        case 6: a = 0.180585;    b = 2.20123;    break; // quercus petraea
        case 7: a = 0.22225;     b = 2.2503739;  break; // fagus sylvatica
        default: break;
    }
    return a * std::pow(dbh, b);
}

*/

float getFormFactor(int type)
{
    switch (type) {
        case 3: return 0.531f; // pinus uncinata (using pinus nigra coefficient)
        case 4: return 0.472f; // abies alba
        case 5: return 0.338f; // betula pendula
        case 6: return 0.454f; // quercus petraea
        case 7: return 0.433f; // fagus sylvatica
        default: return 1.0f;
    }
}


void Windstorm::execute(Simulation* sim)
{
    std::cerr << "WINDSTORM: " << windSpeed << " m/s winds towards " << windDirX << " " << windDirY << std::endl;

    // compute wind map
    MapFloat canopyHeight, density;
    sim->deriveWindMaps(canopyHeight, density);
    canopyHeight.gaussianBlur(6);
    density.gaussianBlur(6);
    computeWindSpeedMap(sim, canopyHeight, density);
    std::cerr << "Computed wind speed map" << std::endl;
               
    // apply wind effect
    sim->rebuildNeighborsAccelerator();
    applyWindEffect(sim, canopyHeight);
    
    // debug
    /*QImage windFracMap(numCellsX, numCellsY, QImage::Format_RGBA8888);
    for (int i = 0; i < numCellsX; i++) {
        for (int j = 0; j < numCellsY; j++) {
            int c = min(255, max(0, int(255*(windSpeedMap.get(i, j)/windSpeed))));
            windFracMap.setPixelColor(i, j, QColor(c,c,c, 255));
        }
    }
    windFracMap.save("wind_" + QString::number(sim->year) + "_exposure.png");
    
    QImage imgH(numCellsX, numCellsY, QImage::Format_RGBA8888);
    QImage imgD(numCellsX, numCellsY, QImage::Format_RGBA8888);
    for (int j = 0; j < numCellsY; j++) {
        for (int i = 0; i < numCellsX; i++) {
            int cH = int(255 * std::min(1.0f, canopyHeight.get(i, j)/5.0f));
            int cD = int(255 * std::max(0.0f, std::min(1.0f, density.get(i, j))));
            imgH.setPixelColor(i, j, QColor(cH, cH, cH));
            imgD.setPixelColor(i, j, QColor(cD, cD, cD));
        }
    }
    imgH.save("wind_" + QString::number(sim->year) + "_heights.png");
    imgD.save("wind_" + QString::number(sim->year) + "_density.png");*/
}


void Windstorm::computeWindSpeedMap(Simulation* sim, MapFloat& simHeight, MapFloat& simDensity)
{
    const int NUM_SAMPLES = 64;
    const float MAX_SAMPLE_DIST = 50;
    const float MAX_SAMPLE_ANGLE = float(M_PI_2);

    // read necessary info about simulation   
    float terDimX, terDimY;
    int terCellsX, terCellsY;
    Terrain* terrain = sim->getTerrain();
    terrain->getTerrainDim(terDimX, terDimY);
    terrain->getGridDim(terCellsX, terCellsY);
    float MIN_SAMPLE_DIST = std::max(terDimX/terCellsX, terDimY/terCellsY);

    // remap densities to 0..1
    float densSum = 0;
    float densSqSum = 0;
    int densCnt = 0;
    for (int i = 0; i < numCellsX; i++) {
        for (int j = 0; j < numCellsY; j++) {
            float d = simDensity.get(i, j);
            if (d > 0) {
                densSum += d;
                densSqSum += d*d;
                densCnt++;
            }
        }
    }
    float densMean = densSum/densCnt;
    float densStd = std::sqrt(densSqSum/densCnt - densMean*densMean);
    float densMin = std::max(0.0f, densMean - densStd);
    float densMax = densMean + densStd;
    float densRng = densMax - densMin;
    for (int i = 0; i < numCellsX; i++) {
        for (int j = 0; j < numCellsY; j++) {
            float d = simDensity.get(i, j);
            d = (d - densMin)/densRng;
            d = std::max(0.0f, std::min(1.0f, d));
            simDensity.set(i, j, d);
        }
    }

    sim->getSimMap()->getDim(numCellsX, numCellsY);
    cellSizeX = terDimX/(numCellsX - 1);
    cellSizeY = terDimY/(numCellsY - 1);

    windSpeedMap.setDim(numCellsX, numCellsY);
    windSpeedMap.initMap();

    // compute per cell wind damage
    for (int j = 0; j < numCellsY; j++) {
        for (int i = 0; i < numCellsX; i++) {

            // current point p(x, y)
            float x = (i + 0.5f)*cellSizeX;
            float y = (j + 0.5f)*cellSizeY;

            float p_hterrain = terrain->getHeightAtCoords(x, y, true);
            float p_hcanopy = simHeight.get(i, j);
            float p_topheight = p_hterrain + p_hcanopy;

            // approximate occlusion factor...
            float windFractionSum = 0;
            float windFractionWeight = 0;

            for (int k = 0; k < NUM_SAMPLES; k++) {

                // sample random direction and distance
                float a = randomUniform(-1, 1);
                float dAngle = MAX_SAMPLE_ANGLE * a*a * (a > 0 ? 1 : -1);
                float dirAngle = windAngle + dAngle;
                float dx = -std::cos(dirAngle);
                float dy = -std::sin(dirAngle);
                float weightDir = std::cos(dAngle);
                
                float r = randomUniform();
                float dist = MIN_SAMPLE_DIST + (MAX_SAMPLE_DIST - MIN_SAMPLE_DIST)*r*r;
                float sx = x + dist*dx;
                float sy = y + dist*dy;
                float weightDist = dist / MAX_SAMPLE_DIST;
                
                float windFraction = 1.0;

                if (terrain->inWorldBounds(vpPoint(sx, 0, sy))) {
                    
                    // get sample location height and density
                    int si = int(sx / cellSizeX);
                    int sj = int(sy / cellSizeY);

                    float s_hterrain  = terrain->getHeightAtCoords(sx, sy, true);
                    float s_hcanopy   = simHeight.get(si, sj);
                    float s_topheight = s_hterrain + s_hcanopy;
                    float s_density   = simDensity.get(si, sj);

                    // p canopy below sampled ground
                    if (p_topheight <= s_hterrain) {
                        windFraction = 0;
                    }
                    // p canopy above sampled ground but below sampled canopy
                    else if (p_topheight <= s_topheight) {
                        windFraction = 1.0 - s_density;
                    }
                    // else: p canopy above sample canopy (no occlusion)
                }

                float weight = weightDir * weightDist;
                windFractionSum += windFraction * weight;
                windFractionWeight += weight;
            }

            // compute and set damage
            float windFraction = windFractionSum / windFractionWeight;
            float windInCell = windFraction * windSpeed;
            windSpeedMap.set(i, j, windInCell);
        }
    }
}


float Windstorm::criticalWindSpeed(Biome* biome, int type, float h, float dbh, bool alive, float& cwsUproot, float& cwsBreak)
{
    float stemVol  = getFormFactor(type) * h * M_PI_4 * dbh * dbh;
    float stemMass  = stemVol * biome->getWoodDensity(type) * (alive ? biome->getWetBiomassFactor(type) : 1);
    //float stemMass2 = getStemMass(type, 100*dbh)            * (alive ? biome->getWetBiomassFactor(type) : 1);

    //float Chegyi = 0.5;
    //float Tc = 4.42 + 122.1*dbh*dbh*h - 0.141*Chegyi - 14.6*dbh*dbh*h*Chegyi;
    float Tc   = 111.7 * dbh * dbh * h;  // Hale et al. 2015, eq 10 
    float Creg = biome->getCreg(type);
    float MOR  = biome->getMOR(type);

    cwsUproot = std::sqrt((Creg * stemMass)/Tc);
    cwsBreak  = std::sqrt((MOR * 1e6 * dbh*dbh*dbh * M_PI)/(32*Tc)); // 1e6 to convert dbh to cm

    return std::min(cwsBreak, cwsUproot);
}


void Windstorm::applyWindEffect(Simulation* sim, const MapFloat& simHeight)
{
    std::vector<SimLivePlant*>& livePlants = sim->getLivePlants();
    std::vector<SimDeadPlant*>& deadPlants = sim->getDeadPlants();
    Biome* biome = sim->getBiome();

    // wind effect on live individuals
    int numUprooted = 0;
    int numTrunksplit = 0;
    int numShaded = 0;
    for (auto it = livePlants.begin(); it != livePlants.end(); it++) {

        SimLivePlant* plnt = *it;
        SimDeadPlant* snag = nullptr;
        std::vector<SimDeadPlant*>::iterator itsnag = deadPlants.end();
        float h   = plnt->height;
        float dbh = plnt->dbh;

        // wind affects only trees and after they've reached a minimum height
        if (!biome->isTree(plnt->pft) || h < Windstorm::WIND_EFFECT_MIN_HEIGHT) {
            continue;
        }

        // snag?
        if (plnt->state == PlantSimState::DEAD) {
            for (auto dit = deadPlants.begin(); dit != deadPlants.end(); dit++) {
                if ((*dit)->parent == plnt) {
                    snag = *dit;
                    itsnag = dit;
                    h *= snag->trunkremains;
                    break;
                }
            }
        }

        // critical wind speed for this individual
        float cws_uproot, cws_break;
        float cws = criticalWindSpeed(biome, plnt->pft, h, dbh, snag != nullptr, cws_uproot, cws_break);

        // wind speed at its location, at canopy top
        int mapX = int(plnt->pos.x / cellSizeX);
        int mapY = int(plnt->pos.z / cellSizeY);
        float ws = windSpeedMap.get(mapX, mapY);

        // wind speed at actual location
        float hMap = simHeight.get(mapX, mapY);
        float hRatio = hMap > 0 ? h/hMap : 1.0f;
        hRatio = std::min(hRatio, 1.0f); // precision errors in cell location might give wrong hMap
        ws *= hRatio;

        // small wind effect, no damage
        if (ws < Windstorm::WIND_SPEED_TWIGS_BREAK) {
            continue;
        }
        // wind < CWS, break branches with a certain probability
        else if (ws < cws) {
            float prob = (ws - Windstorm::WIND_SPEED_TWIGS_BREAK)/(cws - Windstorm::WIND_SPEED_TWIGS_BREAK);
            if (randomUniform() < prob) {
                sim->windDamage(it, DropEventCategory::PRIMARY, windDirX, windDirY);
            }   
        }
        // wind > CWS, only uproot/split trees larger than a certain size
        else if (h > Windstorm::WIND_EFFECT_MIN_HEIGHT) {
            // linear probability of surviving below highest tree in area
            if (randomUniform() > hRatio) {
                numShaded++;
                continue;
            }
            
            DropEventCategory dropCat = DropEventCategory::DPCEND;

            if (ws > cws_uproot && ws > cws_break) {                
                // if ws larger than both critical speeds, decide randomly what happened
                float p_up = cws_uproot/(cws_uproot + cws_break);
                if (randomUniform() < p_up) {
                    dropCat = DropEventCategory::TRUNKSPLIT;
                    numTrunksplit++;
                }
                else {
                    dropCat = DropEventCategory::UPROOTING;
                    numUprooted++;
                }
            }
            else if (ws > cws_break) {
                dropCat = DropEventCategory::TRUNKSPLIT;
                numTrunksplit++;
            }
            else if (ws > cws_uproot) {
                dropCat = DropEventCategory::UPROOTING;
                numUprooted++;
            }

            if (dropCat != DropEventCategory::DPCEND) {
                if (snag) sim->windDamage(itsnag, dropCat, windDirX, windDirY);
                else      sim->windDamage(it, dropCat, windDirX, windDirY);
            }
        }
    }

    std::cout << "Wind event on trees:" << std::endl;
    std::cout << "\t" << numUprooted << " uprooted" << std::endl;
    std::cout << "\t" << numTrunksplit << " trunksplit" << std::endl;
    std::cout << "\t" << numShaded << " shaded by taller trees" << std::endl;
}


/****************************************************************************************
 *                                       DISEASE                                        *
 ****************************************************************************************/


void DiseaseOutbreak::execute(Simulation* sim)
{
    cout << "DISEASE OUTBREAK " << " for species " << disease->targetSpecies << endl;

    float terDimX, terDimY;
    sim->getTerrain()->getTerrainDim(terDimX, terDimY);
    float outX = outbreakX * terDimX;
    float outY = outbreakY * terDimY;
    float r2 = disease->spreadRadius * disease->spreadRadius;

    std::vector<SimLivePlant*>& livePlants = sim->getLivePlants();
    std::vector<SimLivePlant*>::iterator pind;

    // get all candidates to be infected: same species in outbreak location
    std::vector<std::vector<SimLivePlant*>::iterator> susceptibles;
    susceptibles.reserve(livePlants.size());
    for (pind = livePlants.begin(); pind != livePlants.end(); pind++) {
        SimLivePlant* p = *pind;
        if (p->pft == disease->targetSpecies) {            
            float dx = p->pos.x - outX;
            float dy = p->pos.z - outY;
            if (dx*dx + dy*dy < r2)
                susceptibles.push_back(pind);
        }
    }

    if (susceptibles.size() > 0) {
        // infect patient 0
        pind = susceptibles[int(randomUniform()*susceptibles.size())];
        sim->infectPlant(pind, disease);

        // record the "new" disease to keep track of its metrics
        sim->recordDisease(disease);
    }
    else {
        std::cerr << "No plant infected" << std::endl;
    }
}


/****************************************************************************************
 *                                       DROUGHT                                        *
 ****************************************************************************************/


void Drought::execute(Simulation* sim)
{
    cout << "DROUGHT START" << endl;

    sim->beginDrought(droughtValues);
}
