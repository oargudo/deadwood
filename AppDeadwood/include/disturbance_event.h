#ifndef _DISTURBANCE_EVENTS_
#define _DISTURBANCE_EVENTS_

#include "dice_roller.h"
#include "grass.h"
#include "pft.h"
#include "sim.h"
#include <vector>
#include <queue>
#include <set>
#include <sstream>
//#include <QtCore/QString>


namespace ecosim {


class DisturbanceAction 
{
public:
    virtual DisturbanceEventCategory type() const = 0;
    virtual void execute(Simulation* sim) = 0;

    static void setSeed(int seed) { dice.setSeed(seed); }

protected:
    static DiceRoller dice;
    static float randomUniform();
    static float randomUniform(float rmin, float rmax);
};

class ScriptedDisturbanceEvent {
public:
    
    static std::vector<ScriptedDisturbanceEvent*> loadFromFile(const std::string& filepath);

    static bool saveToFile(const std::vector<ScriptedDisturbanceEvent*>& events, const std::string& filepath);

    static int diseaseNum;
    
    int year = 0;
    int month = 0;
    DisturbanceAction* distEvent = nullptr;

    ScriptedDisturbanceEvent() {
        year = month = 0;
        distEvent = nullptr;
    }
    ~ScriptedDisturbanceEvent() {
        delete distEvent;
    }
    
    bool operator<(const ScriptedDisturbanceEvent& e) {
        return year < e.year || (year == e.year && month < e.month);
    }
    bool operator>(const ScriptedDisturbanceEvent& e) {
        return year > e.year || (year == e.year && month > e.month);
    }
};


/**
 * Fire propagation using a Cellular Automata model based on
 * [Alexandridis et al. 2011] Wildland fire spread modelling using cellular automata
 */
class Wildfire : public DisturbanceAction {
protected:

    enum class FireState { NONFLAMMABLE, UNBURNED, BURNING, BURNED };

    std::vector<FireState> cells;    //< cellular automata grid
    int numCellsX, numCellsY;        //< dimensions of grid
    float dimCellX, dimCellY;        //< world dimensions of grid

    std::vector<float> fireOrigin;   //< fire origin positions as list of (x, y, r). If empty, will be one at random
    float windDirX;                  //< wind direction x
    float windDirY;                  //< wind direction y
    float windSpeed;                 //< wind speed in m/s

    bool insideIndex(int i, int j) { return i > 0 && j > 0 && i < numCellsX && j < numCellsY; }
    int index(int i, int j) const { return j * numCellsX + i; }
    FireState& cell(int i, int j) { return cells[index(i, j)]; }
    FireState cell(int i, int j) const { return cells[index(i, j)]; }

    int propagationStep(Terrain* terrain, const MapFloat& vegHeight, const MapFloat& vegDensity, const MapFloat& vegMoisture);
    //void saveDebugImage(const QString& imgName) const;

public:

    Wildfire() { 
        numCellsX = numCellsY = 0;
        dimCellX = dimCellY = 1;
        windDirX = windDirY = windSpeed = 0;
    }

    virtual DisturbanceEventCategory type() const {
        return DisturbanceEventCategory::FIRE; 
    }
    
    /**
     * @brief execute simulates the propagation of fire using the CA model and vegetation information derived from sim
     * @param sim simulation
     */
    virtual void execute(Simulation* sim);

    /**
     * @brief setWind sets wind speed as (v_x, v_y)
     */
    void setWind(float x, float y) { 
        windSpeed = std::sqrt(x*x + y*y);
        windDirX = x/windSpeed; 
        windDirY = y/windSpeed; 
    }

    void getWind(float& windX, float& windY) const {
        windX = windDirX * windSpeed;
        windY = windDirY * windSpeed;
    }

    /**
     * @brief addOrigin adds a fire point of origin at (x, y)
     * @param x,y assumed to be in normalized range [0,1] representing whole domain
     * @param r radius, assumed to be in m
     */
    void addOrigin(float x, float y, float r) {
        fireOrigin.push_back(x);
        fireOrigin.push_back(y);
        fireOrigin.push_back(r);
    }

    std::vector<float> getOrigins() const {
        return fireOrigin;
    }

    void saveDebugImage(const QString& imgName) const;
};

/**
 * Windstorm simulation that can damage trees
 */
class Windstorm : public DisturbanceAction {
protected:
    float windDirX;                 //< wind direction x
    float windDirY;                 //< wind direction y
    float windSpeed;                //< wind speed in m/s
    float windAngle;                //< wind direction angle (0 = +x)

    MapFloat windSpeedMap;
    int numCellsX = 0, numCellsY = 0;
    float cellSizeX = 1, cellSizeY = 1;

public:
    virtual DisturbanceEventCategory type() const {
        return DisturbanceEventCategory::WIND; 
    }

    virtual void execute(Simulation* sim);

    /**
     * @brief setWind sets wind speed as (v_x, v_y)
     */
    void setWind(float x, float y) {
        windSpeed = std::sqrt(x * x + y * y);
        windDirX  = x / windSpeed;
        windDirY  = y / windSpeed;
        windAngle = std::atan2(y, x);
    }

    void getWind(float& windX, float& windY) const {
        windX = windDirX * windSpeed;
        windY = windDirY * windSpeed;
    }

protected:

    void computeWindSpeedMap(Simulation* sim, MapFloat& canopyHeights, MapFloat& vegDensity);
    void applyWindEffect(Simulation* sim, const MapFloat& canopyHeights);

    static float criticalWindSpeed(Biome* biome, int type, float h, float dbh, bool alive, float& cwsUproot, float& cwsBreak);

    static constexpr float WIND_EFFECT_MIN_HEIGHT = 4.0;  // m
    static constexpr float WIND_SPEED_TWIGS_BREAK = 18.0; // m/s according to Beaufort scale (~65 km/h)
};

/**
 * Disease spread simulation
 */
class DiseaseOutbreak : public DisturbanceAction {

protected:
    Disease* disease;
    float outbreakX, outbreakY;  // outbreak location

public:
    virtual DisturbanceEventCategory type() const {
        return DisturbanceEventCategory::DISEASE;
    }

    void setOutbreakLocation(float x, float y) {
        outbreakX = x;
        outbreakY = y;
    }

    void getOutbreakLocation(float& x, float& y) const {
        x = outbreakX;
        y = outbreakY;
    }

    void setDiseaseParams(Disease* dis) {
        disease = dis;
    }

    const Disease* getDiseaseParams() const {
        return disease;
    }

    virtual void execute(Simulation* sim);
};


/**
 * Drought
 */
class Drought : public DisturbanceAction {

protected:
    std::vector<float> droughtValues;

public:
    virtual DisturbanceEventCategory type() const {
        return DisturbanceEventCategory::DROUGHT;
    }

    void setDroughtValues(const std::vector<float>& vals) {
        droughtValues = vals;
    }

    std::vector<float> getDroughtValues() const {
        return droughtValues;
    }

    virtual void execute(Simulation* sim);
};

}

#endif
