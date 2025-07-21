/*******************************************************************************
 *
 * EcoSynth - Data-driven Authoring of Large-Scale Ecosystems (Undergrowth simulator)
 * Copyright (C) 2020  J.E. Gain  (jgain@cs.uct.ac.za)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 ********************************************************************************/


#ifndef PFT
#define PFT
/* file: pft.h
   author: (c) James Gain, 2018
   notes: plant functional type database
*/

#include "terrain.h"
#include "data_importer.h"
// temp
#include "glheaders.h"

namespace ecosim {

    // SIMULATION PARAMETERS

    // growth parameters
    const int shortgrowstart = 11; // month in which short growing season starts
    const int shortgrowend = 4; // and ends
    const int shortgrowmonths = 6; // length of short growing season
    const int longgrowstart = 8; // month in which long growing season starts
    const int longgrowend = 4; // and ends
    const int longgrowmonths = 9; // length of long growing season

    // sun viability
    const float shadeLow = 4.0f; // shade is low end, in hours of sunlight
    const float shadeMedium = 1.5f;
    const float shadeHigh = -2.0f;
    const float scorchLow = 5.0f; // scorch is high end
    const float scorchMedium = 9.0f;
    const float scorchHigh = 14.0f;

    // moisture viability
    const float droughtLow = 30.0f; // drought is low end, in mm
    const float droughtMedium = 20.0f;
    const float droughtHigh = 10.0f;
    const float floodLow = 100.0f; // flood is high end
    const float floodMedium = 120.0f;
    const float floodHigh = 140.0f;

    // temperature viability
    const float coldLow = 6.0f; // tolerance to cold in celsius
    const float coldMedium = -4.0f;
    const float coldHigh = -15.0f;
    const float coldVeryHigh = -25.0f;

    // slope viability
    const float steepMedium = 45.0f; // steepness in degrees
    const float steepHigh = 60.0f;

    enum TreeShapeType
    {
        SPHR,     //< sphere shape, potentially elongated
        BOX,      //< cuboid shape
        CONE,     //< cone shape
        INVCONE,  //< inverted cone, capped by hemisphere
        HEMISPHR, //< hemisphere
        CYL       //< cylinder
    };

    class Viability
    {
    public:
        float c; //< center of viability range
        float r; //< extent of range
        char cmin, cmax; //< codes for upper and lower bounds: 'L', 'M', 'H' for low, medium, high

        Viability() { c = r = 0.0f; }
        ~Viability() {}

        /// Set critical point values for the viability function
        void setValues(char minc, char maxc, float center, float range)
        {
            c = center;
            r = range;
            cmin = minc;
            cmax = maxc;
        }

        /// Set critical point values for the viability function
        void getValues(char& minc, char& maxc, float& center, float& range)
        {
            center = c;
            range = r;
            minc = cmin;
            maxc = cmax;
        }

        /**
         * @brief eval  Evaluate the viability function at a particular value
         * @param val   input value to be evaluated
         * @param neg   range of negative viability (e.g., -0.2)
         * @return      viability result
         */
        float eval(float val, float neg);
    };

    struct AllometryMeasures
    {
        std::vector<float> hght;    //< tree height, all units in metres
        std::vector<float> rad;     //< canopy radius
        std::vector<float> dbh;     //< diameter at breast height
        std::vector<float> branchlen;   //< average length of primary branches
        std::vector<int> numbranches;   //< number of primary branches
    };

    struct PFType
    {
        string code;        //< a mnemonic for the plant type
        //< rendering parameters
        GLfloat basecol[4];  //< base colour for the PFT, individual plants will vary
        float draw_hght;    //< canopy height scaling
        float draw_radius;  //< canopy radius scaling
        float draw_box1;    //< box aspect ratio scaling
        float draw_box2;    //< box aspect ration scaling
        TreeShapeType shapetype; //< shape for canopy: sphere, box, cone
        //< simulation parameters
        Viability sun, wet, temp, slope; //< viability functions for different environmental factors
        float alpha;        //< sunlight attenuation multiplier
        float logseeding;   //< seeding probability multiplier if on logs
        int maxage;         //< expected maximum age of the species in years
        int maxdeadage;     //< expected duration of decay after death in years
        int maxallomage;    //< maximum available age of sampled allometry data, extrapolate beyond this
        float maxhght;      //< expected maximum height in meters
        float maxdbh;       //< expected maximum diameter in meters
        float minhght;      //< height below which trees are culled post simulation for display purposes
        bool isTree;        //< whether or not the species represents a tree
        bool isConifer;     //< whether or not the species represents a conifer
        bool isHerb;        //< whether or not the species represents a herb or forbe
        AllometryMeasures denseAllometries; //< various allometries up to 100 years for plants in a dense setting
        AllometryMeasures isolateAllometries; //< various allometries up to 100 years for plants in an isolated setting
        float qval;         //< logistic function parameter for growth

        //< disturbance event parameters
        float creg;
        float mor;
        float wetbiomass;
        float wooddensity;

        //< growth parameters
        char growth_period; //< long (L) or short (S) growing season
        int grow_months;    //< number of months of the year in the plants growing season
        int grow_start, grow_end; //< start and stop months for growth (in range 0 to 11)
        float grow_m, grow_c1, grow_c2; // terms in growth equation: m * (c1 + exp(c2))

        //< allometry parameters
        char allometry_code; //< allometry patterns: A, B, C, D, or E
        float alm_a, alm_b; //< terms in allometry equation to convert height to canopy radius: r = e ** (a + b ln(h))
        float alm_rootmult; //< multiplier to convert canopy radius to root radius
        float hghtdense;    //< maximum height of species in a dense forest setting
        float hghtisolate;  //< maximum height of species when isolated
        float raddense;     //< maximum radius of species in a dense forest setting
        float radisolate;   //< maximum radius of species when isolated
        float dbhdense;     //< maximum diameter at breast height in a dense forest setting
        float dbhisolate;   //< maximum diameter at breast height when isolated
        float blendense;    //< maximum average primary branch length in a dense forest setting
        float blenisolate;  //< maximum average primary branch length when isolated
        int numbdense;      //< maximum number of primary branches in a dense forest setting
        int numbisolate;    //< maximum number of primary branches when isolated
    };

    struct SubBiome
    {
        string code;    //< a name for the sub-biome
        std::vector<int> canopies; //< indices of plants that occur in the canopy of the sub-biome
        std::vector<int> understorey; //< indices of plants that occur in the understorey of the sub-biome
    };

    class Biome
    {
    private:
        std::vector<PFType> pftypes; //< vector of species in the biome
        std::vector<std::string> catTable;  //< lookup of category names corresponding to category numbers
        std::vector<SubBiome> subTable; //< access to sub-biomes
        std::vector<int> subLookup; //< subbiome number from canopy plant index
        std::string name; //< biome name
        data_importer::common_data* cdata; // access to pdb database

        /**
         * @brief sunmapping    Map sunlight letter codes to response parameters
         * @param pft   plant functional type to store mapping
         * @param cmin  code for shade tolerance
         * @param cmax  code for sunlight tolerance
         */
        void sunmapping(PFType& pft, char cmin, char cmax);

        /** // cerr << "total sun = " << plntpop[p].sunlight << " lightcnt = " << plntpop[p].sunlightcnt  << endl;
                            // cerr << "total water = " << plntpop[p].water << " watercnt = " << plntpop[p].watercnt << endl;
                            // cerr << "simcell occupancy = " << plntpop[p].sunlightcnt << endl;
         * @brief wetmapping    Map moisture letter codes to response parameters
         * @param pft   plant functional type to store mapping
         * @param cmin  code for drought tolerance
         * @param cmax  code for flood tolerance
         */
        void wetmapping(PFType& pft, char cmin, char cmax);

        /**
         * @brief tempmapping    Map temperature letter codes to response parameters
         * @param pft   plant functional type to store mapping
         * @param cmin  code for cold tolerance
         */
        void tempmapping(PFType& pft, char cmin);

        /**
         * @brief slopemapping    Map slope letter codes to response parameters
         * @param pft   plant functional type to store mapping
         * @param cmax  code for steepness tolerance
         */
        void slopemapping(PFType& pft, char cmax);

        /**
         * @brief growthfn  logistic growth function that maps normalized age to normalized growth
         * @param pft   species identifier
         * @param t     age as a proportion of maximum lifespan (in [0,1])
         * @return      corresponding growth value (in [0,1])
         */
        float growthfn(int pft, float t);

        /**
         * @brief cullHght return the culling height corresponding to a particular pft
         * @param id  species identifier
         * @return minumum height for the correponding pft
         */
        float cullHght(int id);

        /**
         * @brief readSpeciesAllometries Read in the measured allometries for a given species
         * @param pft       species identifier
         * @param filename  name of the text file containing the allometry data
         * @param allom     pointer to the allometry into which the data is placed
         * @return true if the file is successfully parsed, false otherwise
         */
        bool readSpeciesAllometries(int pft, const std::string& filename, AllometryMeasures* allom);

        /**
         * @brief allometryDeltas   Use measured allometries and change in age to derive change in other attributes for plants <= 100 years old
         * @param pft           species identifier
         * @param nominalage    nominal age accounting for stunted growth
         * @param sunstr        sun strength for determining degree of crowding
         * @param delage        change in nominal age
         * @param delhght       derived change in height
         * @param delrad        derived change in canopy radius
         * @param deldbh        derived change in diameter at breast height
         * @param delblen       derived change in average branch length
         * @param numb          absolute number of branches
         */
        void allometryDeltas(int pft, float nominalage, float sunstr, float delage, float& delhght, float& delrad, float& deldbh, float& delblen, int& numb);

        /**
         * @brief allometryExtrapolation Use growth curve to extrapolate allometries beyond 100 years
         * @param pft           species identifier
         * @param nominalage    nominal age accounting for stunted growth
         * @param delage        change in nominal age
         * @param delhght       derived change in height
         * @param delrad        derived change in canopy radius
         * @param deldbh        derived change in diameter at breast height
         * @param delblen       derived change in average branch length
         * @param numb          absolute number of branches
         */
        void allometryExtrapolation(int pft, float nominalage, float sunstr, float delage, float& delhght, float& delrad, float& deldbh, float& delblen, int& numb);

        /**
         * @brief underAllomMax Test whether a plant is less than the sampled maximum allometry age
         * @param pft           species identifier
         * @param nominalage    Nominal age in lifespan range [0,1]
         * @return True if less than one hundred years old, false otherwise
         */
        bool underAllomMax(int pft, float nominalage);

    public:

        // soil moisture infiltration parameters
        float slopethresh;       //< slope at which runoff starts increasing linearly
        float slopemax;          //< slope at which water stops being absorbed altogether
        float evaporation;       //< proportion of rainfall that is evaporated
        float runofflevel;       //< cap on mm per month that can be asorbed by the soil before runoff
        float soilsaturation;    //< cap in mm on amount of water that soil can hold
        float waterlevel;        //< surface water level above which a location is marked as a river

        Biome() { cdata = nullptr; }

        ~Biome() { pftypes.clear(); delete cdata; }

        /// numPFTypes: returns the number of plant functional types in the biome
        int numPFTypes() { return (int)pftypes.size(); }

        /// getPFType: get the ith plant functional type in the biome
        PFType* getPFType(int i) { return &pftypes[i]; }

        /// numTrees: returns the number of species that are trees
        int numTrees();

        /// getAlpha: return the alpha canopy factor for PFT type i
        float getAlpha(int i) { return pftypes[i].alpha; }

        /// getMinIdealMoisture: return the start of the ideal moisture range for PFT type i
        float getMinIdealMoisture(int i);

        /// categoryNameLookup: get the category name corresponding to a category number
        void categoryNameLookup(int idx, std::string& catName);

        /// return the sub-biome index for which the species is canopy dominant or co-dominant
        int getSubBiomeFromSpecies(int species);

        /// getSubBiome: get a particular subbiome from the subtable
        SubBiome* getSubBiome(int sbidx) { return &subTable[sbidx]; }

        /// getcreg: return creg, a regression coefficient between stem weight and resistance to overturning (N m/kg)
        float getCreg(int i) { return pftypes[i].creg; }

        /// getMOR: return MOR, Modulus Of Rupture (Pa)
        float getMOR(int i) { return pftypes[i].mor; }

        /// getLogSeeding: return logseeding value for species
        float getLogSeeding(int pft) { return pftypes[pft].logseeding; }

        /// getWetBiomassFactor: return wet biomass factor for disturbance events
        float getWetBiomassFactor(int i) { return pftypes[i].wetbiomass; }

        /// getWoodDensity: return wood density for disturbance events
        float getWoodDensity(int i) { return pftypes[i].wooddensity; }

        /// getGrowMonths: return length of growing season
        int getGrowMonths(int pft) { return pftypes[pft].grow_months; }

        void add_pftype(const PFType& pft)
        {
            pftypes.push_back(pft);
        }

        void add_subbiome(const SubBiome& bme, int bme_idx)
        {
            if ((int)subTable.size() <= bme_idx)		// it seems silly that we are doing this resizing instead of push_back, but I do not want
                // to assume a particular order of the bme_idx indices
                subTable.resize(bme_idx + 1);
            subTable[bme_idx] = bme;
            std::cout << "number of canopies, understorey species in subbiome: " << bme.canopies.size() << ", " << bme.understorey.size() << std::endl;
            for (auto& sp_idx : bme.canopies)
            {
                if ((int)subLookup.size() <= sp_idx)
                    subLookup.resize(sp_idx + 1, -1);
            }
            for (auto& sp_idx : bme.understorey)
            {
                if ((int)subLookup.size() <= sp_idx)
                    subLookup.resize(sp_idx + 1, -1);
            }
            for (int cp_species_idx = 0; cp_species_idx < (int)bme.canopies.size(); cp_species_idx++)
            {
                int sp_idx = bme.canopies[cp_species_idx];
                subLookup[sp_idx] = bme_idx;
            }
        }

        /**
         * @brief linearAllometryDBH  Use species allometry to derive diameter from height (in meters)
         * @param pft       the species index of the plant
         * @param height    plant height in meters
         * @return          diameter at breast height in meters
         */
        float linearAllometryDBH(int pft, float height);

        /**
         * @brief maxGrowthTest Return the maximum possible height achieved by a species
         * @param pft       the species index of the plant
         * @return maximum height achieved by tree under ideal conditions
         */
        float maxGrowthTest(int pft);
        void testgrowthfn();

        /**
         * @brief growth    Provide a 1-month growth increment for a particular plant species relative to its life-cycle
         * @param pft       the species index of the plant
         * @param currhght  current height of the plant in meters
         * @param viability strength of the growth in [0, 1]
         * @return          growth increment in meters
         */
         // float growth(int pft, float currhght, float viability);

         /**
          * @brief growth    Provide a 1-month growth increment for a particular plant species relative to its life-cycle
          * @param pft           the species index of the plant
          * @param nominalage    current nominal age (which accounts for stunting and is not the actual age)
          * @param viability     strength of the growth in [0, 1]
          * @param sunstrength   amount of incident sunlight in [0, 1]
          * @param delage        change in nominal age
          * @param delhght       change in plant height
          * @param delrad        change in plant canopy radius
          * @param deldbh        derived change in diameter at breast height
          * @param delblen       derived change in average branch length
          * @param numb          absolute number of branches
          */
        void growth(int pft, float nominalage, float viability, float sunstrength, float& delage, float& delhght, float& delrad, float& deldbh, float& delblen, int& numb);

        /**
         * @brief viability Calculate the viability score of a particular plant
         * @param pft       the plant functional type index of the plant
         * @param sunlight  avg. hours of sunlight per day in the given month
         * @param moisture  total water available in mm
         * @param temperature   avg. temperature over the given month
         * @param slope     incline in degrees where the plant is located
         * @param negval    how deep negative viability can go
         * @return          viability in [-1, 1]
         */
        float viability(int pft, float sunlight, float moisture, float temperature, float slope, std::vector<float>& adaptvals, float negval, bool print);

        /**
         * @brief Biome::overWet Assuming a moisture check is negative, check whether it is over or under
         * @param pft   the plant functional type index of the plant
         * @param moisture total water available in mm
         * @return True if the test is overwet, false otherwise
         */
        bool overWet(int pft, float moisture);

        /**
         * @brief read  read plant functional type data for a biome from file
         * @param filename  name of file to be read
         * @return true if the file is found and has the correct format, false otherwise
         */
        bool read(const std::string& filename);

        /**
         * @brief write  write plant functional type data for a biome to a file
         * @param filename  name of file to be read
         * @return true if the file is found and has the correct format, false otherwise
         */
        bool write(const std::string& filename);
        bool read_dataimporter(std::string cdata_fpath);
        bool read_dataimporter(data_importer::common_data& cdata);

        /**
         * @brief readAllometries Read in the measured allometries for a all species
         * @param rootfilename  name of the text file prefix common to all allometries
         * @return true if all files successfully parsed, false otherwise
         */
        bool readAllometries(const std::string& rootfilename);

        /**
         * @brief printParams Print the biome parameters for viability to console
         */
        void printParams();

        /**
         * @brief isTree Find out if a particular species is classified as a tree
         * @param pft   species identifier
         * @return      return true if the species represents a tree, false otherwise
         */
        bool isTree(int pft) { return pftypes[pft].isTree; }
        bool isConifer(int pft) { return pftypes[pft].isConifer; }
        bool isHerb(int pft) { return pftypes[pft].isHerb; }
    };
}

#endif
