#pragma once

#include <iostream>
#include <queue>
#include "pft.h"
#include "terrain.h"
#include "dice_roller.h"
#include "eco.h"
#include "sun.h"
#include "moisture.h"

#define DICE_RNG_SEED 421337

// simulation parameters
// #define STEPFILE
#define lapserate 6.5f // in moderately wet air conditions decrease is on average 6.5C per 1000m
#define def_reservecapacity 3.0f
#define def_moisturedemand 4.0f
#define def_seeddistmult 5.0f
#define def_seedprob 0.005f // assuming a 20x20cm cell size = 0.003
// #define def_logseedprob 0.0174f // average of seeding probablilty multiplier over all species is 5.8X more likely
#define def_logseedprob 0.0279f
#define def_logadjust 0.1f
#define def_mortalitybase 0.05f
#define def_viabilityneg -0.2f
#define def_stresswght 1.0f
#define def_hghtthreshold 3.0f // kill plants above this height because they are too tall for undergrowth
#define def_snagdecaymult 2.25f  // multiplier on basic rate of decay for snags
#define def_logdecaymult 2.75f   // multiplier on basic rate of decay for logs
#define def_singleprimary 0.1f;     // single primary branch as proportion of the total
#define def_singlesecondary 0.1f;   // single secondary branch as proportion of the total
#define def_trunkdensity 0.8f;      // wood mass density of trunk section as a proportion of overall plant height
#define def_canopydensity 0.3f;     // wood mass desnity of canopy section as a proportion of overall plant height
#define def_firedmgscale 6.0f;  // scale factor to control how much fire damage is transmitted to vitality reserve
#define def_winddmg 2.0f;      // amount of damage to vitality induced by a wind-based trunk split
#define def_collisiondmg 0.3f; // amount of damage induced by collision
#define def_loggingdbh 0.0f; //0.075f; // only include trees with dbh above this meter threshold in logging
#define def_branchratio 0.8f; // ratio of branch to trunk diameter

namespace ecosim {


    enum class PlantCompartment {
        ROOT,   //< address root compartment of plant
        CANOPY, //< address canopy compartment of plant
        PCEND
    };

    enum class PlantSimState {
        ALIVE, //< plant is active and can grow
        SPROUT, //< growth from an existing tree that is otherwise dead or damaged. Will be co-located
        DEAD, //< plant has died and will transition to the dead list
        STATIC, //< canopy plant that is fixed in size over the course of the simulation
        PSSEND
    };

    enum class DisturbanceEventCategory {
        FIRE,
        WIND,
        DISEASE,
        DROUGHT,
        LOGGING,
        DECAY,      //< not a direct disturbance but used as a cause for drop events
        DEATH,
        COLLISION,
        DECEND
    };

    enum class DropEventCategory {
        UPROOTING,  //< entire tree falls with roots exposed, usually as a result of wind
        TRUNKSPLIT, //< break in the main trunk
        PRIMARY,    //< primary branch and all dependents breaks off
        DPCEND
    };

    enum class SnagDecayClass {
        INTACT,     //< DC1: starts fully intact with leaves. Gradually loses leaves, twigs and secondary branches.
        BARE,       //< DC2: starts with not leaves or twigs and few secondary branches. Gradually loses primary branches.
        DEBRANCHED, //< DC3: Stumps and broken primary branches remain, still full height.
        SNAPPED,    //< DC4: Main trunk reduced in height.
        STUMP,      //< DC5: less than 2m of trunk remaining.
        SDCEND
    };

    // data structure for post-sim distribution analysis
    struct SimLogging {
        int numspecies;     //< number of species in simulation
        int numyears;       //< length of simulation
        //int nummonths;       //< length of simulation in months
        float terrainArea;         //< terrain hect area
        std::vector<string> species; //< species names
        std::vector<std::vector<float>> colour; // species colours
        std::vector<std::vector<float>> datacnt; //< species counts per hectare for each timestep
        std::vector<std::vector<float>> databasal; //< species basal area per hectare for each timestep
        std::vector<int> treemap; //< map species idx to tree idx
    };

    struct DisturbanceEvent {
        DisturbanceEventCategory cat; //< nature of the disturbance event
        int year;   //< simulation year in which the event took place
        float severity; //< severity of the event in terms of tree impact [0,1]
    };

    struct DropEvent {
        DropEventCategory cat; //< category in terms of how much of the tree falls
        DisturbanceEventCategory cause; //< cause of the split
        int year;   //< simulation year in which the drop occured
        long pieceID;   //< unique ID of the resulting log
    };

    struct Disease {
        int diseaseId;              //< id of this disease
        int targetSpecies;          //< id of disease target species
        float spreadRadius;         //< spread radius to search for neighbors to infect
        float spreadProbability;    //< probability of infecting a neighbor in radius
        int monthsContagious;       //< months until plant is no longer contagious, -1 for permanent
        int monthsRecovery;         //< months until disease disappears, -1 for permanent
        int monthsImmune;           //< months until plant can get reinfected, -1 for permanent
        float severity;             //< amount that will be deduced from vitality every month

        int totalInfected;          //< counters for tracking disease stats
        int currentInfected;
    };

    enum class InfectionStage {
        UNINFECTED,
        INFECTED_CONTAGIOUS,
        INFECTED_NONCONTAGIOUS,
        IMMUNE,
        ISTEND
    };

    struct Infection {
        Disease* disease;
        int periodContagious;   //< amount of months this infection can transmit disease
        int periodSick;         //< total amount of months of sickness (including contagious phase)
        InfectionStage stage;   //< current stage along the infection
        int monthsSick;         //< increasing counter of months this infection has lasted
        int remainingImmunity;  //< decreasing counter of months of immunity to this disease, -1 for permanent
    };

    struct SimLivePlant
    {
        long id;        //< unique plant ID
        long parent;    //< source ID for sprouts and dead trees
        PlantSimState state;     //< tree status
        int age;        //< age in years
        float nomage;   //< nominal age accounting for stunted growth in [0,1]
        int numsub;     //< number of subsidiary plant forms, such as shoots and decaying branches
        vpPoint pos;    //< position on terrain in world units
        float gx, gy;   //< position on terrain in terrain grid units
        float height;   //< plant height in metres
        float canopy;   //< canopy radius in metres
        float root;     //< roots radius in metres
        float dbh;      //< diameter at breast height
        float avgblen;  //< average length of primary branches
        int numb;       //< number of primary branches
        float reserves; //< plant biomass pool that sees it through short periods of poor conditions
        float stress;   //< stress factor on plant due to poor conditions
        float vigour;   //< reflection of general health [-1, 1]
        glm::vec4 col;  //< colour variation randomly assigned to plant, but keep consistent to track growth
        int pft;        //< plant functional type to which tree belongs
        float water;    //< accumulator for water contribution
        float sunlight; //< accumulator for sunlight contribution
        float filterlight; //< accumulator for extent of shade
        int watercnt;   //< number of cells intersected by roots
        int sunlightcnt; //< number of cells intersected by canopy
        std::map<int, Infection*> infections; //< current infections, map disease id -> infection
        std::vector<float> hghthistory; //< record of height over the plant's lifespan
        std::vector<float> radhistory; //< record of canopy radius over the plant's lifespan
        std::vector<float> vigourhistory; //< record of plant health over its lifespan
        std::vector<DisturbanceEvent> eventhistory; //< record of disturbance events over the plant's lifespan
        std::vector<DropEvent> drophistory; //< record of log fall over the plant's lifespan
    };

    struct SimDeadPlant
    {
        long id;        //< unique plant ID
        SimLivePlant* parent;    //< source of dead biomass, points to the plant from which this came
        int age;        //< time in years since death
        int deathspan;  //< nominal time from death to total decay (ignoring abiotic decay factors)
        float decay;    //< current level in progression of decay (1.0 at moment of death, 0.0 fully converted to humus)
        float decayage; //< nominal age on decay span, adjusted by abiotic factors
        float deldecay; //< decay increment used for loss of branches
        std::vector<float> decayhistory; //< record of decay over the plant's deathspan
        float trunkremains; //< proportion of original trunk height that remains in the snag (1 = full trunk height, 0 = nothing, not even a stump)
        int initialbranches; //< number of initial primary branches
        int remainingbranches; //< number of primary branches remaining
        std::vector<DisturbanceEvent> eventhistory; //< record of disturbance events over the snag's deathspan
        std::vector<DropEvent> drophistory; //< record of log fall over the snag's deathspan
        SnagDecayClass dc;      //< decay class for snag appearance
        float decaystr;     //< year multiplier from abiotics
        // pft, position, original height and canopy attributes can be derived from parent
    };

    struct SimLog
    {
        long id;        //< unique plant ID
        SimLivePlant* liveparent; //< original live source of the log
        SimDeadPlant* deadparent; //< dead source of the log, NULL if it arose from a live plant
        float dirn;     //< direction of fall relative to the source parent
        float len;      //< length of the log
        float diam;     //< original diameter of the log
        int age;        //< time in years since creation of log
        int deathspan;  //< nominal time from death to total decay (ignoring abiotic decay factors)
        float decay;    //< current level in progression of decay (1.0 at moment of death, 0.0 fully converted to humus)
        float decayage; //< nominal age on decay span, adjusted by abiotic factors
        DropEventCategory cat;  //< category of drop from which log size can be derived

    };

    struct PlntInCell
    {
        SimLivePlant* plnt;            //< index of plant in plant population
        float hght;         //< height of plant in metres, used for sorting
    };

    struct SimParams
    {
        int subcellfactor;      //< no. of times a landscape cell is subdivided along an edge to create the a simulation cell
        float reservecapacity;  //< no. of months in plants carbohydrate reserves
        float moisturedemand;   //< minimum amount of moisture taken by dominant plant
        float seeddistmult;     //< number of times the canopy radius to which seeding extends - not currently used?
        float seedprob;         //< probability in [0,1] of a seed sprouting in a given cell per year (originally 0.00001, or 1e-5)
        float logseedprob;      //< probablility of a seed sprouting in a decaying log seedbed
        float logadjust;        //< adjustment to vitality for saplings below 0.25m to represent protection offered by logs
        float stresswght;       //< weighting of stress as a role in mortality
        float mortalitybase;    //< base term in mortality function
        float viabilityneg;     //< cap on how negative viability can be in a given month
        float firedmgscale;     //< proportion by which fire damage affects carbohydrate reserves
        float winddmg;          //< flat depletion of reserves due to being uprooted by wind
        float collisiondmg;     //< when a tree topples into another
        float loggingdbh;       //< lower limit on dbh for recording
        float branchratio;      //< average ratio of branch to trunk diameter
        // also viability c&r values for each species
    };

    float cmpHeight(PlntInCell i, PlntInCell j);

    struct SimCell
    {
        std::vector<PlntInCell> canopies;    //< all plants whose canopies intersect the cell
        std::vector<PlntInCell> roots;       //< all plnts whose roots intersect the cell
        std::vector<SimDeadPlant*> snags;       //< all snags whose initial canopies intersect the cell, unsorted
        std::vector<SimLog*> logs;        //< all logs that intersect the cell, unsorted
        bool growing;               //< cell is in a growing region
        bool available;				//< cell does not intersect with any canopy or undergrowth tree trunks
        float leftoversun;          //< avg. remaining sunlight per growing month for seeding
        float leftoverwet;          //< avg. remaining moisture per growing month for seeding
        bool canopysorted;          //< is the list of canopies sorted by decreasing tree height
        bool rootsorted;            //< is the list of roots sorted by decreasing tree height
        std::vector<int> seedbank;  //< list of subbiomes that can provide seeds for this cell
    };


    class Simulation;

    class MapSimCell
    {
    private:
        int step;                       //< multiple for conversion from terrain grid to simulation grid
        int gx, gy;                     //< grid dimensions
        std::vector<SimCell> smap;      //< grid of simcells
        DiceRoller* dice;              //< random number generator
        const float radius_mult = 6.0f;
        const float seed_radius_mult = 6.0f;

        /// return the row-major linearized value of a grid position
        inline int flatten(int dx, int dy) { return dx * gy + dy; }

        /// check if coordinates are within the grid
        inline bool ingrid(int x, int y) { return x >= 0 && x < gx && y >= 0 && y < gy; }

        /**
         * @brief notInSeedbank check if a particular subbiome is NOT found in the seedbank for a particular cell
         * @param sbidx subbiome index
         * @param x     x-coord grid position (sub cell accuracy)
         * @param y     y-coord grid position (sub cell accuracy)
         * @return      true if the subbiome is not present, false otherwise
         */
        bool notInSeedbank(int sbidx, int x, int y);

        /**
         * @brief inscribeSeeding  Add a seeding area into the simulation grid
         * @param sbidx     index of subbiome being placed
         * @param x         x-coord grid position (sub cell accuracy)
         * @param y         y-coord grid position (sub cell accuracy)
         * @param rcanopy   radius of seeding
         */
        void inscribeSeeding(int sbidx, int spec_idx, float px, float py, float rcanopy, Terrain* ter);

        /**
         * @brief createLivePlant   Generate a new live plant
         * @param x         x-coord grid position
         * @param y         y-coord grid position
         * @param pft       species of newly created plant
         * @param viability initial seedling viability, which determines size
         * @param simstate  initial plant state (ALIVE, SPROUT, etc.)
         * @param sim       Ongoing simulation state
         * @param ter       Underlying terrain
         * @param biome     Biome statistics for plant species
         * @return pointer to newly allocated plant
         */
        SimLivePlant* createLivePlant(int x, int y, int pft, float viability, PlantSimState simstate, Simulation* sim, Terrain* ter, Biome* biome);

        /**
         * @brief singleSeed    Instantiate a single sapling randomly from the seed bank
         * @param x         x-coord grid position
         * @param y         y-coord grid position
         * @param plntpop   Plant population
         * @param sim       Ongoing simulation state
         * @param ter       Underlying terrain
         * @param biome     Biome statistics for plant species
         */
        bool singleSeed(int x, int y, std::vector<SimLivePlant*>* plntpop, Simulation* sim, Terrain* ter, Biome* biome, std::vector<int>& noseed_count, SimLogging* simlog);

        /**
         * @brief deathSpan Calculate the expected duration in year over which plant remains will decay
         * @param pft       Plant species
         * @param diameter  Diameter of the plant remains in meters
         * @param biome     Biome statistics for plant species
         * @return          Expected duration in years
         */
        int deathSpan(int pft, float diameter, Biome* biome);

        /**
         * @brief inTrunk   Determine whether a cell is occupied by the trunk of a plant
         * @param h     plant height
         * @param dbh   plant diameter at breast height
         * @param pos   plant position on terrain
         * @param x     grid cell x-pos
         * @param y     grid cell y-pos
         * @param ter   underlying terrain
         * @param biome biome statistics for plant species
         * @return      true if cell center falls within the plants trunk, false otherwise
         */
        bool inTrunk(float h, float dbh, vpPoint pos, int x, int y, Terrain* ter, Biome* biome);

    public:

        MapSimCell() { gx = 0; gy = 0; initMap(); dice = new DiceRoller(0, 10000); dice->setSeed(DICE_RNG_SEED); }

        ~MapSimCell() { delMap(); }

        /// getter for grid dimensions
        void getDim(int& dx, int& dy) { dx = gx; dy = gy; }

        /// setter for grid dimensions
        void setDim(int dx, int dy, int mult) { gx = dx; gy = dy; step = mult; initMap(); }

        /// gett for step
        int getStep() { return step; }

        /**
         * @brief toTerGrid Convert from simulation grid coordinates to terrain grid coordinates, taking care to check for out of bounds issues
         * @param mx    simgrid x
         * @param my    simgrid y
         * @param tx    terrain grid x
         * @param ty    terrain grid y
         * @param ter   underlying terrain
         */
        void toTerGrid(int mx, int my, int& tx, int& ty, Terrain* ter) const;

        /// clear the contents of the grid to empty
        void initMap();

        /// completely delete map
        void delMap();

        /// reset per year seeding suitability to zero
        void resetSeeding();

        /// convert from terrain grid position to simulation grid position
        inline float convert(float pos) { return (float)step * pos; }

        /// convert from simulation grid position to terrain grid position
        inline float convert_to_tergrid(float pos) { return (float)pos / step; }

        /// clamp to simgrid domain
        void clamp(int& x, int& y) const;

        /// getter and setter for map elements
        SimCell* get(int x, int y) { return &smap[flatten(x, y)]; }
        void set(int x, int y, SimCell& val) { smap[flatten(x, y)] = val; }

        /**
         * @brief onLog Check whether location coincides with at least one log
         * @param x x-coord terrain grid position
         * @param y y-coord terrain grid position
         * @return  True if there is a log at that location, False otherwise
         */
        bool onLog(float x, float y)
        {
            float grx = convert(x);
            float gry = convert(y);
            return ((int)get(grx, gry)->logs.size() > 0);
        }

        const std::vector<SimCell>& get_smap()
        {
            return smap;
        }

        /**
         * @brief inscribe  Add a new plant into the simulation grid
         * @param plntidx   pointer to plant being placed
         * @param px         x-coord grid position (sub cell accuracy)
         * @param py         y-coord grid position (sub cell accuracy)
         * @param rcanopy   diameter of plant canopy
         * @param rroot     diameter of plant root
         * @param isStatic  true if a canopy plant is being inscribed
         * @param ter       underlying terrain
         */
        void inscribe(SimLivePlant* plntidx, float px, float py, float rcanopy, float rroot, bool isStatic, Terrain* ter, Simulation* sim);

        /**
         * @brief inscribeSnag  Add a new snag into the simulation grid
         * @param snag      pointer to previously created snag
         * @param px        x-coord grid position (sub cell accuracy)
         * @param py        y-coord grid position (sub cell accuracy)
         * @param canopy    diameter of plant canopy
         * @param ter       underlying terrain
         */
        void inscribeSnag(SimDeadPlant* snag, float px, float py, float canopy, Terrain* ter);

        /**
         * @brief inscribeLog  Add a new log into the simulation grid
         * @param log       pointer to previously created log
         * @param px        x-coord grid position (sub cell accuracy)
         * @param py        y-coord grid position (sub cell accuracy)
         * @param ter       underlying terrain
         */
        void inscribeLog(SimLog* log, float px, float py, Terrain* ter);

        /**
         * @brief expand  increase radius of plant in the simulation grid
         * @param plntidx   pointer plant being placed
         * @param x         x-coord grid position (sub cell accuracy)
         * @param y         y-coord grid position (sub cell accuracy)
         * @param prevrcanopy   previous diameter of plant canopy
         * @param prevrroot     previous diameter of plant root
         * @param newrcanopy    new expanded diameter of plant canopy
         * @param newrroot      new expanded diameter of plant root
         */
        void expand(SimLivePlant* plntidx, float px, float py, float prevrcanopy, float prevrroot, float newrcanopy, float newrroot);

        /**
         * @brief uproot Remove a plant from the simulation grid
         * @param plntidx   pointer to plant being uprooted
         * @param px        x-coord grid position (sub cell accuracy)
         * @param py        y-coord grid position (sub cell accuracy)
         * @param rcanopy   diameter of plant canopy
         * @param rroot     diameter of plant root
         * @param ter       terrain containing ecosystem
         */
        void uproot(SimLivePlant* plntidx, float px, float py, float rcanopy, float rroot, Terrain* ter);

        /**
         * @brief uproot Remove a snag from the simulation grid
         * @param snag      pointer to snag being uprooted
         * @param px        x-coord grid position (sub cell accuracy)
         * @param py        y-coord grid position (sub cell accuracy)
         * @param canopy    diameter of snag canopy
         */
        void uprootSnag(SimDeadPlant* snag, float px, float py, float canopy);

        /** simcells.inscribe(plntpop.begin(), x, y, ter->toGrid(sp->canopy), ter->toGrid(sp->root), true, ter, this);
         * @brief uproot Remove a log from the simulation grid
         * @param log       pointer to log being uprooted
         * @param px        x-coord grid position (sub cell accuracy)
         * @param py        y-coord grid position (sub cell accuracy)
         */
        void uprootLog(SimLog* log, float px, float py);

        /**
         * @brief traverse  Traverse the simulation grid and collate sunlight and soil moisture contributions to each plant
         * @param seedable  Does this traversal contribute to seeding determination due to growing season
         * @param precip    Precipitation ratio [0..1] for this particular cycle, used during drought events
         */
        void traverse(std::vector<SimLivePlant*>* plntpop, Simulation* sim, Biome* biome, MapFloat* sun, MapFloat* wet, bool seedable, float precip);

        /**
         * @brief establishSeedBank Render canopy plants into the seedbank
         * @param plntpop   Plant population
         * @param biome     Biome statistics for plant species
         * @param ter       Terrain onto which the plants are placed
         */
        void establishSeedBank(std::vector<SimLivePlant*>* plntpop, int plntpopsize, Biome* biome, Terrain* ter);

        /**
         * @brief uniformSeedBank Add a single subbiomes worth of seeds to the entire growing region of the grid
         * @param sub_biome Subbiome index to be added
         */
        void uniformSeedBank(int sub_biome);

        /**
         * @brief seeding Execute once-yearly seeding test
         * @param plntpop   Existing list of plants to which new seedlings will be added
         * @param sim       Ongoing simulation state
         * @param ter       Terrain onto which plants are placed
         * @param biome     Biome statistics for plant species
         */
        void seeding(std::vector<SimLivePlant*>* plntpop, int plntpopsize, Simulation* sim, Terrain* ter, Biome* biome, SimLogging* simlog);

        /**
         * @brief createDeadPlant   Generate a new standing dead plant
         * @param parent    pointer to the source plant from which the dead plant arises
         * @param ter       Terrain onto which plants are placed
         * @param sim       Ongoing simulation state
         * @param biome     Biome statistics for plant species
         * @return          Pointer to newly created dead plant
         */
        SimDeadPlant* createDeadPlant(SimLivePlant* parent, Terrain* ter, Simulation* sim, Biome* biome);

        /**
         * @brief createLogFromLive Generate a new log from a live plant
         * @param parent    Pointer to the source plant from which this log has fallen
         * @param dirn      angular direction in radians of logfall
         * @param len       the length in meters of the resulting log
         * @param category  class of log
         * @param ter       Terrain onto which plants are placed
         * @param sim       Ongoing simulation state
         * @param biome     Biome statistics for plant species
         * @return          Pointer to newly created log
         */
        SimLog* createLogFromLive(SimLivePlant* parent, float dirn, float len, DropEventCategory category, Terrain* ter, Simulation* sim, Biome* biome);

        /**
         * @brief createLogFromSnag Generate a new log from a snag
         * @param#fig.update_layout(showlegend=True) parent    Pointer to the source snag from which this log has fallen
         * @param dirn      angular direction in radians of logfall
         * @param len       the length in meters of the resulting log
         * @param category  class of log
         * @param ter       Terrain onto which plants are placed
         * @param sim       Ongoing simulation state
         * @return          Pointer to newly created log
         */
        SimLog* createLogFromSnag(SimDeadPlant* parent, float dirn, float len, DropEventCategory category, Terrain* ter, Simulation* sim, Biome* biome);

        /**
         * @brief testSeeding   Test seeding of a single plant from mixed sub-biomes
         * @param pos       position on terrain to force plant seeding for test purposes
         * @param sim       Ongoing simulation state
         * @param ter       Terrain onto which plants are placed
         * @param biome     Biome statistics for plant species
         */
         // void testSeeding(vpPoint pos, Simulation * sim, Terrain * ter, Biome * biome);

         /**
         * @brief visualize Display the canopy grid as an image, with plants coded by colour
         * @param visimg    visualisation image
         * @param plnts     list of plants in simulation
         */
         // void visualize(QImage * visimg, std::vector<SimPlant* > *plnts);

         /**
         * @brief unitTests A set of basic unit tests for the MapSimCell class
         * @param visimg     image into which visualisation of unit tests is written
         * @return  true if the unit tests pass, false otherwise
         */
         // bool unitTests(QImage * visimg, Terrain *ter, Simulation * sim);

         /**
         * @brief validate  check validity of mapsim structure
         * @param plnts     list of plants in simulation
         * @return          true if the mapsim is valid, false otherwise
         */
         // bool validate(std::vector<SimPlant *> *plnts);

         /**
         * @brief extractCanopyHeight   Derive the height above ground level for the canopy level in each grid cell and store as a map
         *                              Note that snags are not included in the canopy
         * @param canopyHeight Map of canopy heights
         */
        void extractCanopyHeight(MapFloat& canopyHeight);

        /**
         * @brief extractDensity    Derive the density of plants in each grid cell and store as a map. Quantized values are: low = 0.0f, medium = 0.5f, high = 1.0f
         * @param density   Map of (low, medium, high) plant densities
         * @param ter       Underlying terrain
         * @param biome     Biome statistics for plant species
         */
        void extractDensity(MapFloat& density, Terrain* ter, Biome* biome);

        /**
         * @brief extractWoodMoisture   Derive the normalized moisture of plants in each grid cell and store as a map. Value are in the range: 0 = completely dry to 1 = completely saturated
         * @param woodMoisture  Map of wood moisture
         * @param ter       Underlying terrain
         * @param sim       Ongoing simulation state
         * @param biome     Biome statistics for plant species
         */
        void extractWoodMoisture(MapFloat& woodMoisture, Terrain* ter, Simulation* sim, Biome* biome);
    };

    // other map inputs: slope, temperature
    class ScriptedDisturbanceEvent;

    class Simulation
    {

    private:
        std::vector<SimLivePlant*> plntpop;  //< all alive plants in the simulation (dead ones are deleted, hence the list data structure)
        std::vector<SimDeadPlant*> snagpop;  //< all decaying snags in the simulation
        std::vector<SimLog*> logpop;         //< all decaying fallen logs in the simulation
        int plntpopsize;                //< count of number of live plants in the simualation
        int snagpopsize;                //< count of number of decaying plants in the simulation
        int logpopsize;                 //< count of number of fallen logs in the simulation
        long newid;                     //< for allocating plant IDs

        MapSimCell simcells;            //< flattened sub-cell grid structure storing intersecting roots and canopies
        std::vector<MapFloat> sunlight; //< local per cell illumination for each month
        std::vector<MapFloat> moisture; //< local per cell moisture for each month
        MapFloat slope;                 //< local per cell slope derived from terrain
        std::vector<float> temperature; //< average monthly temperature
        std::vector<float> cloudiness;  //< average monthly cloud cover
        std::vector<float> rainfall;    //< average monthly rainfall
        float totalrainfall;            //< total rainfall in mm over the course of the year
        Terrain* ter;                   //< terrain over which ecosystem is being simulated
        Biome* biome;                   //< biome within which plants are simulated

        float time;                     //< internal ecosystem clock in years
        DiceRoller* dice;               //< random number generator
        SimLogging simlog;               //< simulation distribution log
        bool suppressDecay;              //< if true then decay processes are ignored
        std::vector<int> speciesCnt;     //< number of individuals per species
        std::vector<int> speciesDeadCnt; //< number of snags per species
        std::vector<MapFloat> viability; //< viability map for each pft

        std::vector<ScriptedDisturbanceEvent*> eventsQueue; //< list of disturbance events
        int eventIndex;                 //< index of the first not yet executed event
        int pastEventIndex;             //< index of the first not yet reported executed event
        std::vector<Disease*> diseases; //< list of diseases occured throughout the simulation
        std::list<float> droughtModifiers; //< precipitation modifiers during drought period
        float totalrainfallDeficit;     //< current year rainfall deficit due to drought

        std::vector<std::list<std::vector<SimLivePlant*>::iterator>> coarseGridLive; //< coarse grid for live plant neighbor queries
        std::vector<std::list<std::vector<SimDeadPlant*>::iterator>> coarseGridDead; //< coarse grid for dead plant neighbor queries
        float coarseGridSpacing;                 //< coarse grid spacing
        int coarseGridSizeX, coarseGridSizeZ;    //< coarse grid dimensions


        /** @brief initSim    Initialize internal data for the simulation
        * @param dx, dy       Main void Simulation::checkpointedSim(EcoSystem * eco, std::string seedbank_file, std::string seedchance_filename, int delYears, std::string out_filename)terrain dimensions
        * @param subcellfactor  number of sub-cells per main simulation cell
        */
        void initSim(int dx, int dy, int subcellfactor);

        /**
         * @brief readMonthlyMap Read a set of 12 maps from file
         * @param filename   name of file to be read
         * @param monthly    content will be loaded into this vector of maps
         */
        bool readMonthlyMap(std::string filename, std::vector<MapFloat>& monthly);

        /**
         * @brief writeMonthlyMap    Write a set of 12 maps to file
         * @param filename   name of file to be written
         * @param monthly    map content to be written
         */
        bool writeMonthlyMap(std::string filename, std::vector<MapFloat>& monthly);

        /// clear internal simulation data
        void delSim();

        /**
         * @brief decayfn Inverse sigmoidal smoothstep decay function from decay to nominal age
         * @param d decay value in [0,1]
         * @return nominal age in [0,1]
         */
        float invdecayfn(float d);

        /**
         * @brief decayfn Sigmoidal smoothstep decay function accounting for nominal age of decayed wood
         * @param t nominal age
         * @return decay value in [0,1]
         */
        float decayfn(float t);

        /**
         * @brief decayeval Use flattened gaussian to evluate decay response to abiotic factors
         * @param val   input abiotic value
         * @param r     abiotic range
         * @param c     abiotic center
         * @return abiotic response value in [0,1]
         */
        float decayeval(float val, float r, float c);

        /**
         * @brief decay Calculate the increment to the decay age of a plant based on its location
         * @param x, y      abiotic map location
         * @param gx, gy    simulation grid location
         * @return      increment to decay age
         */
        float decay(int x, int y, int gx, int gy);

        /**
         * @brief logFall   Create a log that falls from a live plant (primarily due to wind)
         * @param pind      live source plant
         * @param devent    event category for the fall
         * @param dirn      direction in which the log is oriented
         * @param len       length of the log
         */
        void logFall(std::vector<SimLivePlant*>::iterator pind, DropEventCategory devent, float dirn, float len);

        /**
         * @brief logFall
         * @param dind
         * @param devent
         * @param dirn
         * @param len
         */
        void logFall(std::vector<SimDeadPlant*>::iterator dind, DropEventCategory devent, float dirn, float len);

        /**
         * @brief burnSnag  Apply burn damage to a snag as a result of a fire
         * @param dind      Index of the dead plant being burned
         * @param dmg       Extent of the burn damage
         * @return          True if the snag is completely reduced to ash and destroyed, False otherwise
         */
        bool burnSnag(std::vector<SimDeadPlant*>::iterator dind, float dmg);

        /**
         * @brief burnSnag  Apply burn damage to a log as a result of a fire
         *SnagDecayClass Simulation::decayClassLookup(float decay) @param dind      Index of the log being burned
         * @param dmg       Extent of the burn damage
         * @return          True if the log is completely reduced to ash and destroyed, False otherwise
         */
        bool burnLog(std::vector<SimLog*>::iterator lind, float dmg);

        /**
         * @brief decayClassLookup  Look up the decay class corresponding to a decay value
         * @param decay Current decay value
         * @return Corresponding decay class
         */
        SnagDecayClass decayClassLookup(float decay);

        /**
         * @brief decay Update the decay value of a snag to reflect gradual decomposition
         * @param dind  Index of the dead plant being decayed
         * @return      True if dead plant has completely decayed away and needs deleting, False otherwise
         */
        bool decaySnag(std::vector<SimDeadPlant*>::iterator dind);

        /**
         * @brief decay Update the decay value of a log to reflect gradual decomposition
         * @param lind  Index of the log being decayed
         * @return      True if log has completely decayed away and needs deleting, False otherwise
         */
        bool decayLog(std::vector<SimLog*>::iterator lind);

        /**
         * @brief checkDecayedAway Determine whether a plant has completely decayed away and if so, remove
         * @param pind  Index of plant being checked
         * @return whether or not the plant can be completely deleted
         */
        bool checkDecayedAwayPlant(std::vector<SimLivePlant*>::iterator& pind);

        /**
         * @brief checkDecayedAway Remove plants with no extant logs or snags from the simulation
         * @param needsAging If true and snags and/or logs still remain, advance the age of the source tree
         */
        void checkDecayedAway(bool needsAging);

        /**
         * @brief damage    Deplete reserves and increase stress to reflect plant damage
         * @param pind  Index of plant being damaged
         * @param dmg   Extent of damage (relative to reserves)
         */
        void damage(std::vector<SimLivePlant*>::iterator pind, float dmg);

        /**
         * @brief kill Transition plant from live population to dead.
         * @param pind Index of plant being transitioned
         */
        void kill(std::vector<SimLivePlant*>::iterator& pind);

        /**
         * @brief death Test to see whether a plant dies, based on current viability
         * @param plntind   Index of the plant being tested
         * @param stress    Plants stress from environmental factors
         */
        bool death(std::vector<SimLivePlant*>::iterator pind, float stress);

        /**
         * @brief growth    Apply monthly growth equation to plant moderated by current vitality
         * @param pind      Index of the plant being grown
         * @param vitality  Plant vitality from environmental factors
         * @param sunstr    measure of how much shade case by adjacent trees, 1 = none, 0 = total
         */
        void growth(std::vector<SimLivePlant*>::iterator pind, float vitality, float sunstr);


        /**
         * @brief recordLiveLogFall Add an event to the list of log falls for a live tree
         * @param pind  Iterator for the plant undergoing the log fall
         * @param cat   what part of the tree falls
         * @param cause cause of the log fall
         * @param newLogID  unique ID of the resulting log
         */
        void recordLiveLogFall(std::vector<SimLivePlant*>::iterator pind, DropEventCategory cat, DisturbanceEventCategory cause, long newLogID);

        /**
         * @brief recordDeadLogFall Add and event to the list of log falls for a snag
         * @param pind  Iterator for the snag undergoing the log fall
         * @param cat   what part of the tree falls
         * @param cause cause of the log fall
         * @param newLogID  unique ID of the resulting log
         */
        void recordDeadLogFall(std::vector<SimDeadPlant*>::iterator dind, DropEventCategory cat, DisturbanceEventCategory cause, long newLogID);

        void computeViabilityMaps();

        /**
         * @brief averageViability  Calculate and report the average viability of each species over growing areas of the current terrain.
         *                          This is the same basis that is used for seeding probability.
         * @param targetnums        Average number of expected undergrowth plants per canopy tree for each species
         */
        void averageViability(std::vector<float>& targetnums);

        /// Simstep: a single simulation step, equal to one month
        /**
         * @brief simStep   Undertake a single month of ecosystem simulation
         * @param month     Current month
         * @param dtime     Accumulated computation time spent on death and decay
         */
        void simStep(int month, float& dtime);

        /**
         * @brief recordLiveEvent   Add an event to the list of events for a live plant
         * @param pind      Iterator for the plant undergoing the event
         * @param cat       The form of the disturbance event
         * @param severity  Severity of the event on a scale [0,1] with 1 being maximally severe
         */
        void recordLiveEvent(std::vector<SimLivePlant*>::iterator pind, DisturbanceEventCategory cat, float severity);

        /**
         * @brief recordLiveEvent   Add an event to the list of events for a snag
         * @param pind      Iterator for the snag undergoing the event
         * @param cat       The form of the disturbance event
         * @param severity  Severity of the event on a scale [0,1] with 1 being maximally severe
         */
        void recordDeadEvent(std::vector<SimDeadPlant*>::iterator dind, DisturbanceEventCategory cat, float severity);

        void reportYear();

        /**
         * @brief simlogSetup   Initialize simulation log
         */
        void simlogSetup();
        
        /**
         * @brief simlogYear    Record data for the current simulation year to the log
         */
        void simlogYear();


    public:
        /**
         * @brief basalArea calculate the cross-section area of a tree given the dbh
         * @param diam  diameter at breast height
         * @return basal area
         */
        float basalArea(float diam);

        const std::vector<int>& getSpeciesCnt() const;
        const std::vector<int>& getSpeciesDeadCnt() const;

        std::vector<std::string> getSpeciesNames() const;


        void simStep(float& dtime);


        SimParams sparams;              //< simulation parameters for testing outcomes

        int year;
        int month;

        Simulation() { ter = NULL; initSim(1000, 1000, 5); }

        ~Simulation() { ter = NULL; delete dice; delSim(); }

        /** @brief Simulation    Initialize the simulation using an input terrain
        *                        position on the earth, average cloud cover per month, a per cell elevation, canopy height, and
        *                        canopy density. Store output in 12 sunlight maps, one for each month
        * @param terrain        Terrain onto which the plants are placed
        * @param subcellfactor  number of sub-cells per main simulation cell
        */
        Simulation(Terrain* terrain, Biome* simbiome, int subcellfactor);

        void incrPlantPop(int pft) { plntpopsize++; speciesCnt[pft]++; }

        /**
         * @brief getNewPlantID Provide a unique ID for plant creation
         * @return new plant ID
         */
        int getNewPlantID()
        {
            newid++;
            return (newid - 1);
        }

        Terrain* getTerrain() { return ter; }
        Biome* getBiome() { return biome; }

        std::vector<SimLivePlant*>& getLivePlants() { return plntpop; }
        std::vector<SimDeadPlant*>& getDeadPlants() { return snagpop; }
        const std::vector<SimLog*>& getLogs() const { return logpop; }

        std::vector<float> getDataCount(int year, int type) const;
        std::vector<float> getDataBasal(int year, int type) const;


        /**
         * @brief calcSlope    Calculate per cell ground slope
         */
        void calcSlope();

        /**
         * @brief calcMoisture    Calculate per cell and per month moisture values, based on rainfall values and calculation of stream power
         */
        void calcMoisture();

        /**
         * @brief calcMoisture    Calculate per cell and per month sunlight values, based on sun position and cloudiness values
         */
        void calcSunlight(GLSun* glsun, int minstep, int nsamples);
        void reportSunAverages();

        /**
         * @brief setupSim    Set up the simulation on the first run
         * @param        Ecosystem storing canopy trees
         */
        void setupSim(ecosim::EcoSystem* eco, bool nodecay);

        /**
           * @brief setEvents    Set the disturbance events queue
           * @param q       Queue of disturbance events with year/month, is assumed to be ordered
           */
        void setEvents(const std::vector<ScriptedDisturbanceEvent*>& q) {
            eventsQueue = q;
            eventIndex = 0;
            pastEventIndex = 0;
        }

        const std::vector<ScriptedDisturbanceEvent*>& getEvents() const {
            return eventsQueue;
        }
        ScriptedDisturbanceEvent* getEvent(int index) const {
            return eventsQueue.at(index);
        }

        void addEvent(ScriptedDisturbanceEvent* e);
        void editEvent(int index, ScriptedDisturbanceEvent* e);
        void removeEvent(int index);
        std::vector<DisturbanceEvent> getPastNewEvents();



        /*
         * @brief rangeSearch find plants inside a circle defined by center p and radius rad
         * #return the total number of plants + snags found in query
         */
        int rangeSearch(const vpPoint2D& p, float rad, 
            std::list<std::vector<SimLivePlant*>::iterator>& rLivePlants,
            std::list<std::vector<SimDeadPlant*>::iterator>& rDeadPlants);

        /**
          * @brief rebuildNeighborsAccelerator    Recomputes the acceleration structure for NN and range queries
          */
        void rebuildNeighborsAccelerator();

        bool readSun(const std::string& filename)       { return readMonthlyMap(filename, sunlight); }
        bool readMoisture(const std::string& filename)  { return readMonthlyMap(filename, moisture); }
        bool writeSunlight(const std::string& filename) { return writeMonthlyMap(filename, sunlight); }
        bool writeMoisture(const std::string& filename) { return writeMonthlyMap(filename, moisture); }
        bool writeTemperature(const std::string& filename, Terrain* ter);

        /**
         * @brief readClimate    read in average monthly climatic variables
         * @param filename   name of the file containing temperature values
         * @return   false, if the file does not exist, true otherwise
         */
        bool readClimate(std::string filename);


        /**
         * @brief getTemperature Return the temperature at a given position on the terrain
         * @param x      x grid location
         * @param y      y grid location
         * @param mth    month of the year
         * @return       temperature in degrees celsius
         */
        float getTemperature(int x, int y, int mth);

        /// getters for condition maps
        MapFloat* getSunlightMap(int mth) { return &sunlight[mth]; }
        MapFloat* getMoistureMap(int mth) { return &moisture[mth]; }
        MapFloat* getSlopeMap() { return &slope; }
        //SunLight* getSun() { return sunsim; }
        MapFloat* getViabilityMap(int pft) { return &viability[pft]; }

        MapSimCell* getSimMap() { return &simcells; }

        /**
        * @brief simulate   Run the simulation for a certain number of iterations
        * @param eco        ecosystem for storing plant outputs
        * @param delYears   number of years to simulate
        * @param nodecay    if true then decay processes are not included in the simulation
        */
        void simulate(EcoSystem* eco, int delYears, bool nodecay, int exportPeriod, const std::string& exportPrefix);

        /// Hook for Fire Simulator

        /**
        * @brief deriveFireMaps Extract biotic maps from the current state of the simulation for use in a fire simulation.
        *                       Resolution is at the granularity of the simulation
        * @param canopyHeight   Height of the top of the canopy in meters from the ground
        * @param density        A vertical measure (similar to rainfall) of the height of solid wood in meters if it was compacted in a cell (i.e., it would have celldim x celldim x density volume)
        * @param woodMoisture   A measure of how dry (and hence flammable) the plants and ground are in [0,1], where 0 = completely dry and 1 = completely saturated.
        * @return               The length of a simulation cell in meters
        */
        float deriveFireMaps(MapFloat& canopyHeight, MapFloat& density, MapFloat& woodMoisture);

        /**
        * @brief deriveWindMaps Extract biotic maps from the current state of the simulation for use in a wind simulation.
        *                       Resolution is at the granularity of the simulation
        * @param canopyHeight   Height of the top of the canopy in meters from the ground
        * @param density        A vertical measure (similar to rainfall) of the height of solid wood in meters if it was compacted in a cell (i.e., it would have celldim x celldim x density volume)
        * @return               The length of a simulation cell in meters
        */
        float deriveWindMaps(MapFloat& canopyHeight, MapFloat& density);
        void extractCanopyHeight(MapFloat& canopyHeight);

        ///

        // Responses to disturbance events
        void applyFireMap(const MapFloat& fireMap);

        void windDamage(const std::vector<SimLivePlant*>::iterator pind, DropEventCategory type, float wdirX, float wdirY);
        void windDamage(const std::vector<SimDeadPlant*>::iterator pind, DropEventCategory type, float wdirX, float wdirY);

        void  recordDisease(Disease* disease);
        void  infectPlant(std::vector<SimLivePlant*>::iterator pind, Disease* d);
        float updateInfections(std::vector<SimLivePlant*>::iterator pind);
        int   spreadDiseases();
        //void  spreadDiseaseFromInfected(std::vector<SimLivePlant*>::iterator pind);
        void  debugLogDiseaseInfo(const Disease* disease, int y, int m) const;

        void beginDrought(const std::vector<float>& precipRates);


        // Export simulation as different formats
        bool exportSimJSON(const std::string& filename) const;
        bool exportSimPDBTerse(const std::string& filename) const;
        bool exportSimPDBFull(const std::string& filename) const;
    };
}