#include "sim.h"
#include "timer.h"
#include "disturbance_event.h"
#include "glheaders.h"

using namespace ecosim;

void Simulation::initSim(int dx, int dy, int subcellfactor)
{
    // -
    speciesCnt.resize(biome->numPFTypes());
    speciesDeadCnt.resize(biome->numPFTypes());

    simcells.setDim(dx * subcellfactor, dy * subcellfactor, subcellfactor);
    for (int m = 0; m < 12; m++)
    {
        MapFloat sun, mst;
        sun.setDim(dx, dy);
        mst.setDim(dx, dy);
        sunlight.push_back(sun);
        sunlight.back().fill(0.0f);
        moisture.push_back(mst);
        moisture.back().fill(0.0f);
        temperature.push_back(0.0f);
        cloudiness.push_back(0.0f);
        rainfall.push_back(0.0f);
    }
    slope.setDim(dx, dy);
    //sunsim = new SunLight();
    dice = new DiceRoller(0, 1000);
    dice->setSeed(DICE_RNG_SEED);
    DisturbanceAction::setSeed(DICE_RNG_SEED);
    time = 0.0f;
    plntpopsize = 0;
    snagpopsize = 0;
    logpopsize = 0;
    newid = 0;

    // set simulation parameters to default values
    sparams.subcellfactor = subcellfactor;
    sparams.reservecapacity = def_reservecapacity;
    sparams.moisturedemand = def_moisturedemand;
    sparams.seeddistmult = def_seeddistmult;
    sparams.seedprob = def_seedprob;
    sparams.logseedprob = def_logseedprob;
    sparams.logadjust = def_logadjust;
    sparams.stresswght = def_stresswght;
    sparams.mortalitybase = def_mortalitybase;
    sparams.viabilityneg = def_viabilityneg;
    sparams.firedmgscale = def_firedmgscale;
    sparams.winddmg = def_winddmg;
    sparams.collisiondmg = def_collisiondmg;
    sparams.loggingdbh = def_loggingdbh;
    sparams.branchratio = def_branchratio;
}

bool Simulation::readMonthlyMap(std::string filename, std::vector<MapFloat>& monthly)
{
    float val;
    ifstream infile;
    int gx, gy, dx, dy;

    infile.open(filename.c_str(), ios_base::in);
    if (infile.is_open())
    {
        infile >> gx >> gy;
#ifdef STEPFILE
        float step;
        infile >> step; // new format
#endif
        ter->getGridDim(dx, dy);
        if ((gx != dx) || (gy != dy))
            cerr << "Error Simulation::readMonthlyMap: map dimensions do not match terrain" << endl;

        for (int m = 0; m < 12; m++)
            monthly[m].setDim(gx, gy);

        for (int y = 0; y < gy; y++)
            for (int x = 0; x < gx; x++)
                for (int m = 0; m < 12; m++)
                {
                    infile >> val;
                    monthly[m].set(x, y, val);
                }
        infile.close();
        return true;
    }
    else
    {
        cerr << "Error Simulation::readMonthlyMap: unable to open file " << filename << endl;
        return false;
    }
}

bool Simulation::writeMonthlyMap(std::string filename, std::vector<MapFloat>& monthly)
{
    int gx, gy;
    ofstream outfile;
    monthly[0].getDim(gx, gy);

    outfile.open(filename.c_str(), ios_base::out);
    if (outfile.is_open())
    {
        outfile << gx << " " << gy;
#ifdef STEPFILE
        outfile << " 0.9144"; // hardcoded step
#endif
        outfile << endl;
        for (int y = 0; y < gy; y++)
            for (int x = 0; x < gx; x++)
                for (int m = 0; m < 12; m++)
                    outfile << monthly[m].get(x, y) << " ";

        outfile << endl;
        outfile.close();
        return true;
    }
    else
    {
        cerr << "Error Simulation::writeMonthlyMap:unable to open file " << filename << endl;
        return false;
    }

}

void Simulation::delSim()
{
    simcells.delMap();
    for (int m = 0; m < 12; m++)
    {
        sunlight[m].delMap();
        moisture[m].delMap();
    }
    std::vector<SimLivePlant*>::iterator plnt;
    for (plnt = plntpop.begin(); plnt != plntpop.end(); plnt++)
        delete (*plnt);
    plntpop.clear();
    std::vector<SimDeadPlant*>::iterator dead;
    for (dead = snagpop.begin(); dead != snagpop.end(); dead++)
        delete (*dead);
    snagpop.clear();
    std::vector<SimLog*>::iterator log;
    for (log = logpop.begin(); log != logpop.end(); log++)
        delete (*log);
    logpop.clear();
    temperature.clear();
    cloudiness.clear();
    rainfall.clear();
    slope.delMap();
    //if (sunsim) delete sunsim;
    time = 0.0f;
    plntpopsize = 0;
    snagpopsize = 0;
    logpopsize = 0;
}

float Simulation::invdecayfn(float d)
{
    // inverse of smoothstep function
    // inverse of sigmoidal-like (smoothstep) function

    // bounds checks
    if (d < -0.0f)
    {
        cerr << "Error Simulation::invdecayfn: decay is below zero" << endl;
        return 0.0f;
    }

    if (d > 1.0f)
    {
        return  1.0f;
        cerr << "Error Simulation::invdecayfn: decay is above maximum" << endl;
    }
    return 0.5f - sinf(asinf(1.0f - 2.0f * d) / 3.0f);
}

float Simulation::decayfn(float t)
{
    // sigmoidal-like (smoothstep) function for decay

    // bounds checks
    if (t < 0.0f)
    {
        // cerr << "Error Simulation::decayfn: age is below zero" << endl;
        return 1.0f;
    }

    if (t > 1.0f)
    {
        // cerr << "Error Simulation::decayfn: age is above maximum" << endl;
        return  0.0f;

    }
    return 1.0f - t * t * (3.0f - 2.0f * t);
}
float Simulation::decayeval(float val, float r, float c)
{
    // evaluate exponential. Power term introduces flat top
    float s = log(0.2f) / pow(r / 2.0, 4.5f);
    float v = pow(M_E, s * pow(fabs(val - c), 4.5f));

    return v + 1.4f; // range [0.5, 1.5]
}

float Simulation::decay(int x, int y, int gx, int gy)
{
    float tresponse, mresponse, wet, hot, response;

    wet = moisture[7].get(gx, gy); // totalrainfall - totalrainfallDeficit + simcells.get(x,y)->leftoverwet;
    hot = getTemperature(gx, gy, 6);

    // moderate decay step according to temperature and moisture
    tresponse = decayeval(hot, 35.0f, 20.0f); // center at 30 deg celsius in July
    mresponse = decayeval(wet, 100.0f, 60.0f); // center at soil moisture in August
    response = tresponse * mresponse;
    response = std::max(response, 0.8f); // clamp to [0.8, 1.2]
    response = std::min(response, 1.2f);
    return response;
}

void Simulation::logFall(std::vector<SimLivePlant*>::iterator pind, DropEventCategory devent, float dirn, float len)
{
    SimLog* log;
    long brid;

    log = simcells.createLogFromLive((*pind), dirn, len, devent, ter, this, biome);
    logpop.push_back(log);
    logpopsize++;
    brid = log->id;
    recordLiveLogFall(pind, devent, DisturbanceEventCategory::WIND, brid);
}

void Simulation::logFall(std::vector<SimDeadPlant*>::iterator dind, DropEventCategory devent, float dirn, float len)
{
    SimLog* log;
    long brid;

    log = simcells.createLogFromSnag((*dind), dirn, len, devent, ter, this, biome);
    logpop.push_back(log);
    logpopsize++;
    brid = log->id;
    recordDeadLogFall(dind, devent, DisturbanceEventCategory::DECAY, brid);
}

bool Simulation::burnSnag(std::vector<SimDeadPlant*>::iterator dind, float dmg)
{
    (*dind)->decay -= dmg;

    if ((*dind)->decay <= 0.01f) // snag is completely destroyed
    {
        (*dind)->decay = 0.0f;
        return true;
    }

    // adjust related decay parameters
    (*dind)->decayage = invdecayfn((*dind)->decay) * (float)(*dind)->deathspan;
    (*dind)->dc = decayClassLookup((*dind)->decay);
    (*dind)->deldecay = 0.0f; // used to control branch fall. Since branches are assumed burned away, it is not needed

    return false;
}

bool Simulation::burnLog(std::vector<SimLog*>::iterator lind, float dmg)
{
    (*lind)->decay -= dmg;

    if ((*lind)->decay <= 0.01f) // snag is completely destroyed
    {
        (*lind)->decay = 0.0f;
        return true;
    }

    // adjust related decay parameters
    (*lind)->decayage = invdecayfn((*lind)->decay) * (float)(*lind)->deathspan;

    return false;
}

SnagDecayClass Simulation::decayClassLookup(float decay)
{
    if (decay < def_snagDC4)
        return SnagDecayClass::STUMP;
    else if (decay < def_snagDC3)
        return SnagDecayClass::SNAPPED;
    else if (decay < def_snagDC2)
        return SnagDecayClass::DEBRANCHED;
    else if (decay < def_snagDC1)
        return SnagDecayClass::BARE;
    else
        return SnagDecayClass::INTACT;
}

bool Simulation::decaySnag(std::vector<SimDeadPlant*>::iterator dind)
{
    int x, y;
    float deldecayage;
    float dirn, len, del, lenvar;

    // find temperature and moisture values at the snags location
    x = (int)simcells.convert((*dind)->parent->gx);
    y = (int)simcells.convert((*dind)->parent->gy);
    simcells.clamp(x, y);
    deldecayage = decay(x, y, (*dind)->parent->gx, (*dind)->parent->gy);
    //deldecayage = 1.0f;
    (*dind)->decayage += deldecayage * def_snagdecaymult;
    (*dind)->decaystr = deldecayage;

    // smoothstep on decayage as proportion of total deathspan
    del = (*dind)->decay;


    /*
    if((* dind)->decayage / (float) (* dind)->deathspan < 0.0f)
    {
        cerr << "decaySnag error:" << endl;
        cerr << "decayage = " << (* dind)->decayage << " deathspan = " << (* dind)->deathspan << endl;
    }*/

    if ((*dind)->deathspan > 0)
        (*dind)->decay = decayfn((*dind)->decayage / (float)(*dind)->deathspan);
    else
        (*dind)->decay = 0.0f;
    del -= (*dind)->decay;
    (*dind)->deldecay += del;

    // test for break off of branches or trunk to form logs
    // depends on the decay class of the snag
    switch ((*dind)->dc)
    {
    case SnagDecayClass::INTACT: // gradual loss of secondary branches, but not modelled as logs because they are too small
        if ((*dind)->decay < def_snagDC1) // transition point
        {
            (*dind)->dc = SnagDecayClass::BARE;
            (*dind)->deldecay = def_snagDC1 - (*dind)->decay;
        }
        break;
    case SnagDecayClass::BARE:   // gradual loss of primary branches
        // drop a primary branch proportionately over the duration of DC2
        // that is evenly spread transition points on branch drops
        if ((*dind)->deldecay > (def_snagDC1 - def_snagDC2) / (*dind)->initialbranches)
        {
            (*dind)->deldecay = 0.0f;
            (*dind)->remainingbranches--;
            dirn = (dice->generate() / 1000.0f) * PI2;
            lenvar = 0.9f + 0.2f * (dice->generate() / 1000.0f);
            len = lenvar * (*dind)->parent->avgblen;
            logFall(dind, DropEventCategory::PRIMARY, dirn, len);
        }

        if ((*dind)->decay < def_snagDC2) // transition point
        {
            (*dind)->dc = SnagDecayClass::DEBRANCHED;
            /*
            if((* dind)->remainingbranches > 0)
            {
                 cerr << "ERROR decaySnag: " << (* dind)->remainingbranches << " remaining in DC2"  << endl;
                 cerr << "for pft " << (* dind)->parent->pft << endl;
            }*/
            while ((*dind)->remainingbranches > 0)
            {
                (*dind)->remainingbranches--;
                dirn = (dice->generate() / 1000.0f) * PI2;
                lenvar = 0.9f + 0.2f * (dice->generate() / 1000.0f);
                len = lenvar * (*dind)->parent->avgblen;
                logFall(dind, DropEventCategory::PRIMARY, dirn, len);
            }
        }
        break;
    case SnagDecayClass::DEBRANCHED: // chance of high trunk snap
        //if ((*dind)->remainingbranches > 0)
        //    cerr << "WARNING: Debranched Snag still has brances remaining" << endl;
        if ((*dind)->decay < def_snagDC3) // transition point
        {
            // snap point in [0.2, 0.6] of height
            float snap = 0.4f * (dice->generate() / 1000.0f) + 0.2f;
            dirn = (dice->generate() / 1000.0f) * PI2;
            len = (*dind)->parent->avgblen;
            len *= snap;
            (*dind)->trunkremains = snap;
            logFall(dind, DropEventCategory::TRUNKSPLIT, dirn, len);

            (*dind)->dc = SnagDecayClass::SNAPPED;
        }
        break;
    case SnagDecayClass::SNAPPED: // chance of low trunk snap
        //if ((*dind)->remainingbranches > 0)
        //    cerr << "WARNING: Snapped Snag still has brances remaining" << endl;
        if ((*dind)->decay < def_snagDC4) // transition point
        {
            // snap point in [0.01, 0.05] of height
            float snap = 0.04f * (dice->generate() / 1000.0f) + 0.01f;
            dirn = (dice->generate() / 1000.0f) * PI2;
            len = (*dind)->trunkremains * (*dind)->parent->height;
            len -= snap * (*dind)->parent->avgblen;
            (*dind)->trunkremains = snap;
            logFall(dind, DropEventCategory::TRUNKSPLIT, dirn, len);

            (*dind)->dc = SnagDecayClass::STUMP;
        }
        break;
    case SnagDecayClass::STUMP: // nothing more to lose
        //if ((*dind)->remainingbranches > 0)
        //    cerr << "WARNING: Stump Snag still has brances remaining" << endl;
        break;
    default:
        break;
    }
    return ((*dind)->decay <= 0.01f);
}

bool Simulation::decayLog(std::vector<SimLog*>::iterator lind)
{
    int x, y;
    float deldecayage;
    SimLivePlant* parent;

    parent = (*lind)->liveparent;

    // find temperature and moisture values at the snags location
    x = (int)simcells.convert(parent->gx);
    y = (int)simcells.convert(parent->gy);
    simcells.clamp(x, y);
    deldecayage = decay(x, y, parent->gx, parent->gy);
    (*lind)->decayage += deldecayage * def_logdecaymult;

    // smoothstep on decayage as proportion of total deathspan
    if ((*lind)->decayage / (float)(*lind)->deathspan < 0.0f)
    {
        cerr << "decayLog error:" << endl;
        cerr << "decayage = " << (*lind)->decayage << " deathspan = " << (*lind)->deathspan << endl;
    }

    (*lind)->decay = decayfn((*lind)->decayage / (float)(*lind)->deathspan);

    // adjust dbh of log as it decays

    return ((*lind)->decay <= 0.01f);
}

bool Simulation::checkDecayedAwayPlant(std::vector<SimLivePlant*>::iterator& pind)
{
    if ((*pind)->state == PlantSimState::DEAD && (*pind)->numsub == 0)
    {
        int sx = (int)simcells.convert((*pind)->gx);
        int sy = (int)simcells.convert((*pind)->gy);
        simcells.clamp(sx, sy);
        simcells.get(sx, sy)->growing = true; // flag mapsimcell location as growable so that a new plant can be placed there
        simcells.get(sx, sy)->available = true;
        return true;
    }
    else
    {
        return false;
    }
}

void Simulation::checkDecayedAway(bool needsAging)
{
    std::vector<SimLivePlant*>::iterator pind;

    // check that dead plants still have subsidiary snags and logs
    pind = plntpop.begin();
    while (pind != plntpop.end())
    {
        if (checkDecayedAwayPlant(pind)) // all subsidiaries are dead and plant can be removed
        {
            // cnt
            speciesCnt[(*pind)->pft]--;

            //int indx = simlog.treemap[(*pind)->pft];
            //if (indx > 0)
            //    // remove basal area
            //    simlog.databasal[simlog.numyears][3 * indx] -= basalArea((*pind)->dbh);

            // remove from plant population, but make sure iterator remains valid
            delete (*pind);
            pind = plntpop.erase(pind); // require c++11
            plntpopsize--;

        }
        else
        {
            if (needsAging)
                (*pind)->age += 1;
            pind++;
        }
    }
}

void Simulation::simStep(int month, float& dtime)
{
    std::vector<float> adaptvals(4);
    bool shortgrow;
    float x, y, h;
    std::vector<SimLivePlant*>::iterator plant;
    std::vector<SimDeadPlant*>::iterator zombie;
    std::vector<SimDeadPlant*>::iterator zcheck;
    std::vector<SimLog*>::iterator log;

    // drought modifier
    float precipRatio = droughtModifiers.size() > 0 ? droughtModifiers.front() : 1.0f;
    totalrainfallDeficit += std::max(0.0f, (1.0f - precipRatio) * rainfall[month]);

    // traverse all cells contributing moisture and sunlight to plant pool
    shortgrow = (month >= shortgrowstart || month <= shortgrowend);
    simcells.traverse(&plntpop, this, biome, &sunlight[month], &moisture[month], shortgrow, precipRatio);

    int ndied = 0;
    int ntotal_died = 0;
    int p = 0;
    int ngonelogs = 0;
    int nnewlogs = 0;

    Timer simd, simp;

    simd.start();

    // decay all snags and logs at the end of the year
    if (month == 11)
    {
        // decay logs
        log = logpop.begin();
        while (log != logpop.end())
        {
            if (decayLog(log)) // true if log has decayed away
            {
                // reduce subsidiary count in parent
                (*log)->liveparent->numsub--;

                // remove from simulation grid
                ter->toGrid((*log)->liveparent->pos, x, y, h);
                simcells.uprootLog((*log), x, y);

                delete (*log);
                log = logpop.erase(log); // require c++11
                logpopsize--;
                ngonelogs++;
            }
            else
            {
                (*log)->age += 1;
                log++;
            }
        }

        // decay snags
        zombie = snagpop.begin();
        nnewlogs = (int)logpop.size();
        while (zombie != snagpop.end())
        {
            if (decaySnag(zombie))
            {
                // cnt
                speciesDeadCnt[(*zombie)->parent->pft]--;


                // reduce subsidiary count in parent
                (*zombie)->parent->numsub--;

                // remove from simulation grid
                ter->toGrid((*zombie)->parent->pos, x, y, h);
                simcells.uprootSnag((*zombie), x, y, ter->toGrid((*zombie)->parent->canopy));

                delete (*zombie);
                zombie = snagpop.erase(zombie); // require c++11
                snagpopsize--;
            }
            else
            {
                (*zombie)->age += 1;
                (*zombie)->decayhistory.push_back((*zombie)->decay);
                zombie++;
            }
        }
        nnewlogs = int(logpop.size()) - nnewlogs;

        checkDecayedAway(true);
    }

    // check for death and seeding once a year and store height in history
    if (month == 11)
        for (plant = plntpop.begin(); plant != plntpop.end(); plant++)
        {
            if ((*plant)->state == PlantSimState::ALIVE)
            {
                bool died;

                // check for death from old age or stress
                died = death(plant, (*plant)->stress);

                // calculate average vigour, currently stores acumulated vitality
                (*plant)->vigour = (*plant)->vigour / 12.0f;
                (*plant)->vigour -= (*plant)->stress / 12.0f;
                (*plant)->stress = 0.0f;
                if (died)
                {
                    // remove from simulation grid
                    ter->toGrid((*plant)->pos, x, y, h);
                    simcells.uproot((*plant), x, y, ter->toGrid((*plant)->canopy), ter->toGrid((*plant)->root), ter);

                    // transition plant to deadlist unless decay is supressed
                    // do not remove from live population until the snag and all logs have decayed away
                    kill(plant);
                    ntotal_died++;

                    // update currently infected counters
                    for (auto itinf = (*plant)->infections.begin(); itinf != (*plant)->infections.end(); itinf++) {
                        if (itinf->second->stage == InfectionStage::INFECTED_CONTAGIOUS ||
                            itinf->second->stage == InfectionStage::INFECTED_NONCONTAGIOUS)
                            itinf->second->disease->currentInfected--;
                    }
                }
                else // update plant history
                {
                    (*plant)->hghthistory.push_back((*plant)->height);
                    (*plant)->radhistory.push_back((*plant)->canopy / 2.0f);
                    (*plant)->vigourhistory.push_back((*plant)->vigour);
                }
            }
        }

    simd.stop();
    dtime += simd.peek();


    simp.start();

    // for(plnt = plntpop.begin(); plnt != plntpop.end(); plnt++)
    int step = (int)plntpop.size() / 16;
    #pragma omp parallel for
    for (int pp = 0; pp < (int)plntpop.size(); pp += step)
    {
        int qstep = (int)plntpop.size() / 16;
        int qstart = pp;
        int qend = pp + qstep;
        if (qend > (int) plntpop.size())
            qend = (int)plntpop.size();
        for (int q = qstart; q < qend; q++)
        {

            std::vector<SimLivePlant*>::iterator plnt;
            plnt = plntpop.begin() + q;

            float sun, wet, temp, slope, str, shade;
            float pool, stress = 0.0f, vitality = 0.0f;

            if ((*plnt)->state == PlantSimState::ALIVE)
            {
                // cerr << endl << "PLANT " << p << " of " << plntpopsize << endl;
                if ((*plnt)->sunlightcnt > 0)
                {
                    sun = (*plnt)->sunlight / (float)(*plnt)->sunlightcnt; // average sunlight
                    shade = (*plnt)->filterlight / (float)(*plnt)->sunlightcnt; // average filtering
                }
                else
                {
                    sun = (*plnt)->sunlight;
                    shade = (*plnt)->filterlight;
                }
                wet = (*plnt)->water / (float)(*plnt)->watercnt; // average moisture
                // cerr << "total sun = " << plntpop[p].sunlight << " lightcnt = " << plntpop[p].sunlightcnt  << endl;
                // cerr << "total water = " << plntpop[p].water << " watercnt = " << plntpop[p].watercnt << endl;
                // cerr << "simcell occupancy = " << plntpop[p].sunlightcnt << endl;
                temp = getTemperature((*plnt)->gx, (*plnt)->gy, month);
                slope = getSlopeMap()->get((*plnt)->gx, (*plnt)->gy);
                str = biome->viability((*plnt)->pft, sun, wet, temp, slope, adaptvals, sparams.viabilityneg, false);

                // adjustment for logs since they are a favourable environment for tree seedlings
                float x, y, h;
                ter->toGrid((*plnt)->pos, x, y, h);
                if (simcells.onLog(x, y))
                    if (biome->isTree((*plnt)->pft) && (*plnt)->height < 0.25f)
                        str += sparams.logadjust;

                // account for plant reserve pool
                pool = (*plnt)->reserves + str;

                // account for infections
                pool -= updateInfections(plnt);

                // cerr << "pool = " << pool << " ";
                if (pool < 0.0f) // potential death due to stress
                {
                    (*plnt)->stress += -1.0f * pool;
                    pool = 0.0f;
                }
                else if (pool > sparams.reservecapacity) // reserves are full so growth is possible
                {
                    vitality = pool - sparams.reservecapacity;
                    pool = sparams.reservecapacity;
                }
                (*plnt)->reserves = pool;

                // cerr << "vitality = " << vitality << " reserves = " << plntpop[p].reserves << " ";
                // cerr << "stress = " << plntpop[p].stress << " ";


                // cerr << "vitality = " << vitality << " reserves = " << (* plnt)->reserves << " ";
                // cerr << "stress = " << (* plnt)->stress << " ";

                // use vitality to determine growth based on allometry
                // but check if this falls in the growing season
                // cerr << "pre-growth" << endl;

                int gs = biome->getPFType((*plnt)->pft)->grow_start;
                int ge = biome->getPFType((*plnt)->pft)->grow_end;
                bool inseason;
                if (gs > ge)
                    inseason = month >= gs || month <= ge;
                else
                    inseason = month >= gs && month <= ge;
                if (inseason)
                {
                    // growth depends on density of surrounding plants for which we use sunlight as a proxy
                    shade = 2.0f * std::min(shade, 0.5f);
                    float sunstr = 1.0f - shade;
                    growth(plnt, vitality, sunstr);
                }
                (*plnt)->vigour += vitality;
                // cerr << "post-growth" << endl;
            }
        }
        // p++;
    }
    simp.stop();

    // spread diseases
    spreadDiseases();
    //for (Disease* d : diseases)
    //    debugLogDiseaseInfo(d, year, month);

    // advance drought
    if (droughtModifiers.size() > 0)
        droughtModifiers.pop_front();
    totalrainfallDeficit = 0; // put back to 0 to accumulate yearly cycles    

}

void Simulation::simStep(float& dtime)
{
    if (this->month == 0) {
        simcells.resetSeeding(); // clear seeding

        // biome->testgrowthfn();
        dtime = 0.0f; // time required for death and decay

        //// new basal vector for new year
        //std::vector<float> lastBasal(biome->numTrees() * 3);
        //if (simlog.numyears > 0) { // copy last basal data
        //    lastBasal = simlog.databasal[simlog.numyears - 1];
        //}

        //simlog.databasal.push_back(lastBasal);
    }

    // disturbance events
    while (eventIndex < eventsQueue.size() && this->year == eventsQueue.at(eventIndex)->year && this->month + 1 == eventsQueue.at(eventIndex)->month)
    {
        ScriptedDisturbanceEvent* currEvt = eventsQueue.at(eventIndex);
        currEvt->distEvent->execute(this);
        eventIndex++;
    }


    simStep(this->month, dtime);


    if (this->month == 11) {
        Timer sims;
        sims.start();
        simcells.seeding(&plntpop, plntpopsize, this, ter, biome, &simlog);
        sims.stop();
        //cerr << "SEEDING TOOK " << sims.peek() << "s" << endl;
        //reportYear();
        year++;

        simlogYear();

        simlog.numyears++;
    }

    this->month = (this->month + 1) % 12;
}

const std::vector<int>& Simulation::getSpeciesCnt() const
{
    //std::vector<int> _speciesCnt(biome->numPFTypes(), 0);
    //for (auto it = plntpop.begin(); it != plntpop.end(); it++) {
    //    _speciesCnt[(*it)->pft]++;
    //}

    //for (int i = 0; i < _speciesCnt.size(); ++i) {
    //    assert(_speciesCnt[i] == speciesCnt[i]);
    //}

    return speciesCnt;
}

const std::vector<int>& Simulation::getSpeciesDeadCnt() const
{
    //std::vector<int> _speciesCnt(biome->numPFTypes(), 0);
    //for (auto it = snagpop.begin(); it != snagpop.end(); it++) {
    //    _speciesCnt[(*it)->parent->pft]++;
    //}

    //for (int i = 0; i < _speciesCnt.size(); ++i) {
    //    assert(_speciesCnt[i] == speciesDeadCnt[i]);
    //}
    return speciesDeadCnt;
}

std::vector<std::string> Simulation::getSpeciesNames() const
{
    std::vector<std::string> names(biome->numPFTypes());
    for (int i = 0; i < biome->numPFTypes(); ++i) {
        names[i] = biome->getPFType(i)->code;
    }
    return names;
}

void Simulation::recordLiveEvent(std::vector<SimLivePlant*>::iterator pind, DisturbanceEventCategory cat, float severity)
{
    DisturbanceEvent devent;
    devent.cat = cat;
    devent.severity = severity;
    devent.year = (*pind)->age;
    (*pind)->eventhistory.push_back(devent);
}

void Simulation::recordDeadEvent(std::vector<SimDeadPlant*>::iterator dind, DisturbanceEventCategory cat, float severity)
{
    DisturbanceEvent devent;
    devent.cat = cat;
    devent.severity = severity;
    devent.year = (*dind)->age;
    (*dind)->eventhistory.push_back(devent);
}

void Simulation::reportYear()
{
    cerr << "\t plants/snags/logs: " << plntpop.size() << "/" << snagpop.size() << "/" << logpop.size() << endl;
    cerr << "pfts: " << endl;
    // calculate species statistics
    std::vector<int> speciesCnt(biome->numPFTypes(), 0);
    std::vector<float> maxHght(biome->numPFTypes(), 0.0);
    std::vector<float> avgHght(biome->numPFTypes(), 0.0);
    std::vector<int> snagCnt(biome->numPFTypes(), 0);
    std::vector<float> avgDeathAge(biome->numPFTypes(), 0.0);
    std::vector<int> numDeathZero(biome->numPFTypes(), 0);
    std::vector<float> avgDecayStr(biome->numPFTypes(), 0.0);
    std::vector<float> avgDeathDiameter(biome->numPFTypes(), 0.0);

    for (auto it = plntpop.begin(); it != plntpop.end(); it++) {
        speciesCnt[(*it)->pft]++;
        if ((*it)->height > maxHght[(*it)->pft])
            maxHght[(*it)->pft] = (*it)->height;
        avgHght[(*it)->pft] += (*it)->height;
    }

    for (auto it = snagpop.begin(); it != snagpop.end(); it++) {
        snagCnt[(*it)->parent->pft]++;
        avgDeathAge[(*it)->parent->pft] += (float)(*it)->deathspan;
        avgDecayStr[(*it)->parent->pft] += (*it)->decaystr;
        avgDeathDiameter[(*it)->parent->pft] += ((*it)->parent->dbh);
        if ((*it)->deathspan == 0)
            numDeathZero[(*it)->parent->pft]++;
    }
    for (int i = 0; i < biome->numPFTypes(); i++) {
        cerr << i << ": " << speciesCnt[i] << " [avgh = " << avgHght[i] / (float)speciesCnt[i] << " maxh = " << maxHght[i];
        if (snagCnt[i] > 0) {
            cerr << " avgdeathspan = " << avgDeathAge[i] / (float)snagCnt[i] << " avg diameter = " << avgDeathDiameter[i] / (float)snagCnt[i] << " zero proportion = " << (float)numDeathZero[i] / (float)snagCnt[i] << " avgdecaystr = " << avgDecayStr[i] / (float)snagCnt[i] << "]" << endl;
        }
        else {
            cerr << " avgdecayage = N/A avgdecaystr = N/A]" << endl;
        }
    }
    cerr << endl;

}

void Simulation::simlogSetup()
{
    int tidx = 0;
    simlog.numyears = 0;
    //simlog.nummonths = 0;
    simlog.numspecies = 3 * biome->numTrees(); // record live plants, snags, and log counts for each species

    cerr << "Terrain area in hectares is " << getTerrain()->getTerrainHectArea() << endl;
    // clear any existing leftover data
    simlog.species.clear();
    for (int i = 0; i < (int)simlog.colour.size(); i++)
        simlog.colour[i].clear();
    simlog.colour.clear();
    for (int i = 0; i < (int)simlog.datacnt.size(); i++)
        simlog.datacnt[i].clear();
    simlog.datacnt.clear();
    for (int i = 0; i < (int)simlog.databasal.size(); i++)
        simlog.databasal[i].clear();
    simlog.databasal.clear();

    // create species names
    for (int s = 0; s < biome->numPFTypes(); s++)
    {
        if (biome->isTree(s))
        {
            simlog.treemap.push_back(tidx); tidx++;
            float r, g, b;
            std::vector<float> col, colsnag, collog;
            string name = biome->getPFType(s)->code;
            r = biome->getPFType(s)->basecol[0]; g = biome->getPFType(s)->basecol[1]; b = biome->getPFType(s)->basecol[2];
            col.resize(3, 0.0f); colsnag.resize(3, 0.0f); collog.resize(3, 0.0f);
            col[0] = r; col[1] = g; col[2] = b;
            colsnag[0] = 0.75f * r; colsnag[1] = 0.75f * g; colsnag[2] = 0.75f * b;
            collog[0] = 0.5f * r; collog[1] = 0.5f * g; collog[2] = 0.5f * b;
            simlog.colour.push_back(col);
            simlog.colour.push_back(colsnag);
            simlog.colour.push_back(collog);
            simlog.species.push_back(name);
            simlog.species.push_back(name + " snags");
            simlog.species.push_back(name + " logs");
        }
        else
        {
            simlog.treemap.push_back(-1); // not recorded, since it isn't a tree
        }
    }

    simlog.terrainArea = getTerrain()->getTerrainHectArea();

    simlogYear();
}

float Simulation::basalArea(float diam)
{
    float r = 0.5f * diam;

    return PI * r * r;
}

void Simulation::simlogYear()
{
    int tidx;
    std::vector<float> spcdata;
    std::vector<float> spcbasal;
    spcdata.resize(biome->numTrees() * 3, 0.0f);
    spcbasal.resize(biome->numTrees() * 3, 0.0f);

    // record species counts for live plants
    for (auto it = plntpop.begin(); it != plntpop.end(); it++) {
        tidx = simlog.treemap[(*it)->pft];
        if (tidx != -1 && (*it)->dbh >= sparams.loggingdbh && (*it)->state == PlantSimState::ALIVE)
        {
            spcdata[3 * tidx] += 1.0f;
            spcbasal[3 * tidx] += basalArea((*it)->dbh);
        }
    }

    // record species counts for snags
    for (auto it = snagpop.begin(); it != snagpop.end(); it++) {
        int pft = (*it)->parent->pft; // snag pft
        tidx = simlog.treemap[pft];
        if (tidx != -1 && (*it)->parent->dbh >= sparams.loggingdbh)
        {
            spcdata[3 * tidx + 1] += 1.0f;
            spcbasal[3 * tidx + 1] += basalArea((*it)->parent->dbh);
        }
    }

    // record species counts for logs
    for (auto it = logpop.begin(); it != logpop.end(); it++) {
        int pft = (*it)->liveparent->pft; // log pft
        tidx = simlog.treemap[pft];
        if (tidx != -1 && (*it)->diam >= sparams.loggingdbh)
        {
            spcdata[3 * tidx + 2] += 1.0f;
            spcbasal[3 * tidx + 2] += basalArea((*it)->diam);
        }
    }

    // normalize by area in hectares
    float area = getTerrain()->getTerrainHectArea();
    for (auto& it : spcdata)
        it = it / area;
    for (auto& it : spcbasal)
        it = it / area;

    simlog.datacnt.push_back(spcdata);
    simlog.databasal.push_back(spcbasal);
}

void Simulation::simulate(EcoSystem* eco, int delYears, bool nodecay, int exportPeriod, const std::string& exportPrefix)
{
    year = 0;
    for (int y = 0; y < delYears; y++)
    {
        Timer sims;
        
        simcells.resetSeeding(); // clear seeding

        // biome->testgrowthfn();
        float dtime = 0.0f; // time required for death and decay

        for (int m = 0; m < 12; m++)
        {
            simStep(m, dtime);
        }

        //if (y == 0)
        //    averageViability(targetnums);

        sims.start();
        simcells.seeding(&plntpop, plntpopsize, this, ter, biome, &simlog);
        sims.stop();
        cerr << "SEEDING TOOK " << sims.peek() << "s" << endl;
        //simy.stop();

        //simlogYear(); // add statistics for the year to the log

        float stime = sims.peek(); // time required for seeding
        //cerr << "YEAR " << y << " took " << ytime << "s" << "\t";
        cerr << "\t plants/snags/logs: " << plntpop.size() << "/" << snagpop.size() << "/" << logpop.size() << endl;

        cerr << "pfts: " << endl;
        // calculate species statistics
        std::vector<int> speciesCnt(biome->numPFTypes(), 0);
        std::vector<float> maxHght(biome->numPFTypes(), 0.0);
        std::vector<float> avgHght(biome->numPFTypes(), 0.0);
        std::vector<int> snagCnt(biome->numPFTypes(), 0);
        std::vector<float> avgDeathAge(biome->numPFTypes(), 0.0);
        std::vector<int> numDeathZero(biome->numPFTypes(), 0);
        std::vector<float> avgDecayStr(biome->numPFTypes(), 0.0);
        std::vector<float> avgDeathDiameter(biome->numPFTypes(), 0.0);

        for (auto it = plntpop.begin(); it != plntpop.end(); it++) {
            speciesCnt[(*it)->pft]++;
            if ((*it)->height > maxHght[(*it)->pft])
                maxHght[(*it)->pft] = (*it)->height;
            avgHght[(*it)->pft] += (*it)->height;
        }

        for (auto it = snagpop.begin(); it != snagpop.end(); it++) {
            snagCnt[(*it)->parent->pft]++;
            avgDeathAge[(*it)->parent->pft] += (float)(*it)->deathspan;
            avgDecayStr[(*it)->parent->pft] += (*it)->decaystr;
            avgDeathDiameter[(*it)->parent->pft] += ((*it)->parent->dbh);
            if ((*it)->deathspan == 0)
                numDeathZero[(*it)->parent->pft]++;
        }
        for (int i = 0; i < biome->numPFTypes(); i++) {
            cerr << i << ": " << speciesCnt[i] << " [avgh = " << avgHght[i] / (float)speciesCnt[i] << " maxh = " << maxHght[i];
            if (snagCnt[i] > 0) {
                cerr << " avgdeathspan = " << avgDeathAge[i] / (float)snagCnt[i] << " avg diameter = " << avgDeathDiameter[i] / (float)snagCnt[i] << " zero proportion = " << (float)numDeathZero[i] / (float)snagCnt[i] << " avgdecaystr = " << avgDecayStr[i] / (float)snagCnt[i] << "]" << endl;
            }
            else {
                cerr << " avgdecayage = N/A avgdecaystr = N/A]" << endl;
            }
        }
        cerr << endl;

        year++;

        if (exportPeriod > 0 && year % exportPeriod == 0) {
            std::string fname = exportPrefix;
            fname += std::to_string(year);
            fname += ".pdb";
            //writeSimTerse(fname, getBiome());
        }

    }
}

void Simulation::recordLiveLogFall(std::vector<SimLivePlant*>::iterator pind, DropEventCategory cat, DisturbanceEventCategory cause, long newLogID)
{
    DropEvent drp;
    drp.cat = cat;
    drp.cause = cause;
    drp.year = (*pind)->age;
    drp.pieceID = newLogID;
    (*pind)->drophistory.push_back(drp);
}

void Simulation::recordDeadLogFall(std::vector<SimDeadPlant*>::iterator dind, DropEventCategory cat, DisturbanceEventCategory cause, long newLogID)
{
    DropEvent drp;
    drp.cat = cat;
    drp.cause = cause;
    drp.year = (*dind)->age;
    drp.pieceID = newLogID;
    (*dind)->drophistory.push_back(drp);
}

void Simulation::computeViabilityMaps() {

    std::vector<float> adaptvals(4);

    // drought modifier
    const float precipRatio = droughtModifiers.size() > 0 ? droughtModifiers.front() : 1.0f;

    // traverse all cells contributing moisture and sunlight to plant pool

    //simcells.traverse(&plntpop, this, biome, &sunlight[month], &moisture[month], shortgrow, precipRatio);
    int tx, ty;
    //int gx, gy;

    //simcells.getDim(gx, gy);
    ter->getGridDim(tx, ty);

    float v, sun, wet, temp, slope;

    viability.resize(biome->numPFTypes());
    for (MapFloat& m : viability)
        m.setDim(tx, ty);


    for (int x = 0; x < tx; x++)
        for (int y = 0; y < ty; y++) {
            float leftoversun = 0;
            float leftoverwet = 0;
            for (int month = 0; month < 12; ++month)
            {
                const bool seedable = (month >= shortgrowstart || month <= shortgrowend);

                const float sunlight_v = sunlight[month].get(x, y);

                // add value to sunlight accumulation for later seeding check
                if (seedable)
                    leftoversun += sunlight_v;

                const float moisture_v = moisture[month].get(x, y) * precipRatio;

                // remainder spread equally
                if (moisture_v > 0.0f)
                {
                    // add share of leftover water to water accumulator for later seeding check
                    if (seedable)
                        leftoverwet += moisture_v; // potential water for seed access
                }
                else {
                    //cout << "moisture_v <= 0.0f" << moisture_v << endl;
                }
            }

            // average sunlight and moisture over growing season
            sun = leftoversun / (float)shortgrowmonths;
            wet = leftoverwet / (float)shortgrowmonths;
            temp = getTemperature(x, y, shortgrowend); // use end of growing season
            slope = getSlopeMap()->get(x, y);

            for (int pft = 0; pft < biome->numPFTypes(); ++pft) {
                v = biome->viability(pft, sun, wet, temp, slope, adaptvals, sparams.viabilityneg, false);
                viability[pft].set(x, y, v);
            }
        }
}

void Simulation::averageViability(std::vector<float>& targetnums)
{
    std::vector<float> adaptvals(4);
    std::vector<double> totadaptvals(4);
    double totstr;
    float sun, wet, temp, slope;
    int gx, gy, step, avgcount, wetcount, numover, nonviable, lowviable, highviable;

    simcells.getDim(gx, gy);
    step = simcells.getStep();

    cerr << "Average Viability by Species: " << endl;

    // sum viability per species over growing areas of map
    for (int s = 0; s < biome->numPFTypes(); s++)
    {
        totstr = 0.0f;
        avgcount = 0;
        wetcount = 0; numover = 0; nonviable = 0; lowviable = 0; highviable = 0;
        for (int i = 0; i < 4; i++)
            totadaptvals[i] = 0.0f;

        if (targetnums[s] >= 0.0f)
        {
            for (int x = 0; x < gx; x++) {
                for (int y = 0; y < gy; y++) {
                    if (simcells.get(x, y)->growing)
                    {
                        float str;
                        int tgx, tgy;

                        // convert from x, y simcell to tx, ty terrain grid
                        simcells.toTerGrid(x, y, tgx, tgy, ter);

                        // average sunlight and moisture over growing season
                        sun = simcells.get(x, y)->leftoversun / (float)shortgrowmonths;
                        wet = simcells.get(x, y)->leftoverwet / (float)shortgrowmonths;
                        temp = getTemperature(tgx, tgy, shortgrowend); // use end of growing season
                        slope = getSlopeMap()->get(tgx, tgy);
                        str = biome->viability(s, sun, wet, temp, slope, adaptvals, sparams.viabilityneg, false);
                        totstr += (double)str;
                        for (int i = 0; i < 4; i++)
                            totadaptvals[i] += (double)adaptvals[i];
                        avgcount++;
                        if (adaptvals[1] < 0.0f) // negative wetness
                        {
                            if (biome->overWet(s, wet))
                                numover++;
                            wetcount++;
                        }
                        if (str >= 0.3f) // high viability here
                        {
                            highviable++;
                        }
                        if (str < 0.3f && str > 0.0f) // low viability here
                        {
                            lowviable++;
                        }
                        if (str <= 0.0f) // non-viable here
                        {
                            nonviable++;
                        }
                    }
                }
            }

            cerr << s << ": (avg. viability, sun, wet, temp, slope) = " << totstr / (double)avgcount << ", " << totadaptvals[0] / (double)avgcount << ", " << totadaptvals[1] / (double)avgcount << ", " << totadaptvals[2] / (double)avgcount << ", " << totadaptvals[3] / (double)avgcount;
            cerr << endl;
            cerr << "overwet proportion = " << (float)numover / (float)wetcount << " of " << (float)wetcount / (float)avgcount << " proportion of cells" << endl;
            cerr << "high viability proportion = " << (float)highviable / (float)avgcount << endl;
            cerr << "low viability proportion = " << (float)lowviable / (float)avgcount << endl;
            cerr << "non-viabile proportion = " << (float)nonviable / (float)avgcount << endl;
            float c, r, hmin, hmax;
            c = biome->getPFType(s)->slope.c;
            r = biome->getPFType(s)->slope.r;
            hmin = biome->getPFType(s)->slope.cmin;
            hmax = biome->getPFType(s)->slope.cmax;
            cerr << "Slope Response: c = " << c << " r = " << r << " hmin = " << hmin << " hmax = " << hmax << endl;

            // str = biome->viability(s, sun, wet, temp, 90.0f, adaptvals, sparams.viabilityneg, false);
            // cerr << "90 def slope test = " << adaptvals[3] << endl;
        }
        cerr << endl;
    }
}

void Simulation::damage(std::vector<SimLivePlant*>::iterator pind, float dmg)
{
    (*pind)->reserves -= dmg;
    if ((*pind)->reserves < 0.0f)
    {
        (*pind)->stress -= (*pind)->reserves;
        (*pind)->reserves = 0.0f;
    }
}

void Simulation::kill(std::vector<SimLivePlant*>::iterator& pind)
{
    if (!suppressDecay)
    {
        SimDeadPlant* dp = simcells.createDeadPlant((*pind), ter, this, biome);
        snagpop.push_back(dp);
        snagpopsize++;

        // cnt
        speciesDeadCnt[(*pind)->pft]++;

    }

    // add death event to the plants eventlist TODO
    //recordLiveEvent(pind, DisturbanceEventCategory::DEATH, 1.0f);
}

bool Simulation::death(std::vector<SimLivePlant*>::iterator pind, float stress)
{
    bool dead = false;

    // age and stress factors are combined using a probabilistic apprach
    // use constant background mortaility to model age effects
    // use bioclamatic envelope and a carbohydrate pool to model stress effects

    // background mortality
    float ageFactor = 1.0f / (biome->getPFType((*pind)->pft)->maxage);
    float age = 1.0f - pow(sparams.mortalitybase, ageFactor);
    // cerr << "maxage = " << biome->getPFType((* currind)->pft)->maxage << " ";

    // stress mortality. stress score is remainder after effects of carbohydrate pool accounted for
    float p = age + sparams.stresswght * stress;

    // test against a uniform random variable in [0,1000]
    int r = dice->generate();
    dead = r/1000.0f < p;

    if (dead) // if plant has died change its state
        (*pind)->state = PlantSimState::DEAD;

    return dead;
}

void Simulation::growth(std::vector<SimLivePlant*>::iterator pind, float vitality, float sunstr)
{
    PFType* pft;
    int numb;
    float x, y, h, prevcanopy, prevroot, prevhght, delh, delr, dela, deldbh, delblen;

    assert((*pind)->state != PlantSimState::DEAD);

    if (vitality > 0.0f)
    {
        // apply growth equation for particular pft moderated by vitality
        pft = biome->getPFType((*pind)->pft);
        prevhght = (*pind)->height;
        prevcanopy = (*pind)->canopy;
        biome->growth((*pind)->pft, (*pind)->nomage, min(vitality, 1.0f), sunstr, dela, delh, delr, deldbh, delblen, numb);

        //float oldBasalArea = basalArea((*pind)->dbh);

        (*pind)->height += delh;
        (*pind)->canopy += 2.0f * delr;
        (*pind)->nomage += dela;
        (*pind)->dbh += deldbh;
        (*pind)->avgblen += delblen;
        (*pind)->numb = numb;
        // cerr << (* pind)->height << " ";

        //float newBasalArea = basalArea((*pind)->dbh);

        //int indx = simlog.treemap[(*pind)->pft];
        //if (indx > 0)
        //    // add new basal area
        //    simlog.databasal[simlog.numyears][3 * indx] += newBasalArea - oldBasalArea;

        /*
        if(isnan((* pind)->height) && !isnan(prehght))
        {
            cerr << "Simulation::growth error: negative height = " << (* pind)->height << endl;
            cerr << "canopy = " << (* pind)->canopy << " nomage = " << (* pind)->nomage << " vitality = " << vitality << endl;
            cerr << "prehght = " << prevhght << " sunstr = " << sunstr << endl;
        }*/

        // use allometries for canopy and root
        prevroot = (*pind)->root;
        (*pind)->root = (*pind)->canopy * pft->alm_rootmult;

        // adjust coverage in simcells accordingly
        ter->toGrid((*pind)->pos, x, y, h);
        simcells.expand((*pind), x, y, ter->toGrid(prevcanopy), ter->toGrid(prevroot), ter->toGrid((*pind)->canopy), ter->toGrid((*pind)->root));

        // cerr << " * GROWTH by " << vitality << " to new hght " << plntpop[pind].height << " * ";
    }
}

Simulation::Simulation(Terrain* terrain, Biome* simbiome, int subcellfactor)
{
    int dx, dy;

    ter = terrain;
    biome = simbiome;
    ter->getGridDim(dx, dy);
    initSim(dx, dy, subcellfactor);
    calcSlope();
}

void Simulation::setupSim(ecosim::EcoSystem* eco, bool nodecay)
{
    year = 0;
    month = 0;
    suppressDecay = nodecay;

    std::fill(speciesCnt.begin(), speciesCnt.end(), 0);
    std::fill(speciesDeadCnt.begin(), speciesDeadCnt.end(), 0);


    //////
    
    auto start_tp = std::chrono::steady_clock::now().time_since_epoch();

    std::vector<Plant> plnts;
    int dx, dy;

    eco->pickAllPlants(ter); // must gather plants before vectorizing

    // make sure list of plants is cleared
    std::vector<SimLivePlant*>::iterator plnt;
    for (plnt = plntpop.begin(); plnt != plntpop.end(); plnt++)
        delete (*plnt);
    plntpop.clear();
    std::vector<SimDeadPlant*>::iterator dead;
    for (dead = snagpop.begin(); dead != snagpop.end(); dead++)
        delete (*dead);
    snagpop.clear();
    std::vector<SimLog*>::iterator log;
    for (log = logpop.begin(); log != logpop.end(); log++)
        delete (*log);
    logpop.clear();
    simcells.initMap();
    ter->getGridDim(dx, dy);
    // cerr << "ter dimensions = " << dx << " X " << dy << endl;
    simcells.getDim(dx, dy);
    // cerr << "sim map dimensions = " << dx << " X " << dy << endl;

    std::set<int> all_specidxes;

    // If there are existing plants in the ecosystem copy them across
    // This allows us to continue where a previous simulation left off
    // iterate over plant functional types
    /*
    for(int pft = 0; pft < biome->numPFTypes(); pft++)
    {
        eco->getPlants()->vectoriseByPFT(pft, plnts);
        // cerr << "PFT = " << pft << "numplants = " << (int) plnts.size() << endl;
        for(int p = 0; p < (int) plnts.size(); p++)
        {
            SimLivePlant * sp = new SimLivePlant;
            sp->state = PlantSimState::ALIVE;
            sp->age = plnts[p].height / biome->getPFType(pft)->maxhght * (float) biome->getPFType(pft)->maxage;
            sp->pos = plnts[p].pos;
            // cerr << "POSITION = " << sp.pos.x << ", " << sp.pos.y << endl;
            sp->height = plnts[p].height;
            hghts.push_back(sp->height);
            // cerr << "HEIGHT = " << sp.height << endl;
            sp->canopy = plnts[p].canopy;
            // cerr << "CANOPY = " << sp.canopy << endl;
            sp->root = sp->canopy;
            sp->reserves = sparams.reservecapacity;
            sp->stress = 0.0f;
            sp->col = plnts[p].col;
            sp->pft = pft;
            sp->water = 0.0f;
            sp->sunlight = 0.0f;
            plntpop.push_back(sp);
            plntpopsize++;

            all_specidxes.insert(sp->pft);

            // inscribe plant into simcells
            // cerr << "TERRAIN CELL SIZE = " << ter->getCellExtent() << endl;

            ter->toGrid(sp->pos, x, y, h);
            if(x > dx-1)
                x = dx-1;
            if(y > dy-1)
                y = dy-1;
            sp->gx = x;
            sp->gy = y;
            simcells.inscribe(plntpop.begin(), x, y, ter->toGrid(sp->canopy), ter->toGrid(sp->root), true, ter, this);
        }
    }*/

    // create seedbank
    float tw, th;
    uint gw, gh;
    ter->getGridDim(gw, gh);
    ter->getTerrainDim(tw, th);

    simcells.uniformSeedBank(0); // arbitrarily choose sub-biome 0 for seeding
    // for authoring we will need to do something more sophisticated, such as a map

    computeViabilityMaps();
    simlogSetup();

    auto end_tp = std::chrono::steady_clock::now().time_since_epoch();
}

void Simulation::calcSlope()
{
    int dx, dy;
    Vector up, n;
    std::vector<int> hist;

    hist.resize(10, 0);

    // slope is dot product of terrain normal and up vector
    up = Vector(0.0f, 1.0f, 0.0f);
    ter->getGridDim(dx, dy);
    slope.setDim(dx, dy);
    slope.fill(0.0f);
    for (int x = 0; x < dx; x++)
        for (int y = 0; y < dy; y++)
        {
            ter->getNormal(x, y, n);
            float rad = acos(up.dot(n));
            float deg = RAD2DEG * rad;
            slope.set(x, y, deg);
            hist[(int)(deg / 10.0f)]++;
        }

    cerr << "SLOPE HISTOGRAM" << endl;
    for (int h = 0; h < 10; h++)
        cerr << h << ": " << (float)hist[h] / (float)(dx * dy) << endl;
}

void Simulation::calcMoisture()
{
    Timer t;
    MoistureSim wet;

    t.start();

    cerr << "MOISTURE SIM PARAMETERS:" << endl;
    cerr << " slope threshold = " << biome->slopethresh << " slope max = " << biome->slopemax << " evaporation = " << biome->evaporation;
    cerr << " runofflevel = " << biome->runofflevel << " soilsaturation = " << biome->soilsaturation << " water level = " << biome->waterlevel << endl;
    wet.simSoilCycle(ter, &slope, rainfall, biome->slopethresh, biome->slopemax, biome->evaporation,
        biome->runofflevel, biome->soilsaturation, biome->waterlevel, moisture);
    t.stop();
    cerr << "MOISTURE SIMULATION TOOK " << t.peek() << "s IN TOTAL" << endl;
}

void Simulation::calcSunlight(GLSun* glsun, int minstep, int nsamples)
{
    Timer t;
    t.start();

    SunLight sunsim;
    sunsim.setLatitude(ter->getLatitude());
    sunsim.setNorthOrientation(Vector(0.0f, 0.0f, -1.0f));
    sunsim.setTerrainDimensions(ter);
    glsun->setScene(getTerrain());

    bool projectenable = true;
    int startm = 1, endm = 12, mincr = 1;
    if (minstep == 0)
    {
        projectenable = false;
        minstep = 60;
    }
    bool diffuseenable = nsamples > 0;

    MapFloat diffusemap;
    std::vector<float> sunhours;

    std::cout << "projecting sun..." << std::endl;
    sunsim.projectSun(ter, sunlight, glsun, sunhours, minstep, startm, endm, mincr, projectenable);

    std::cout << "diffusing sun..." << std::endl;
    sunsim.diffuseSun(ter, &diffusemap, glsun, nsamples);

    std::cout << "merging sun..." << std::endl;
    sunsim.mergeSun(sunlight, &diffusemap, cloudiness, sunhours);

    t.stop();

    cerr << "SUNLIGHT SIMULATION TOOK " << t.peek() << "s IN TOTAL" << endl;

   /* int dw, dh;
    diffusemap.getDim(dw, dh);
    QImage dsun(dw, dh, QImage::Format_RGB888);
    for (int i = 0; i < dw; i++)
        for (int j = 0; j < dh; j++) {
            float v = diffusemap.get(i, j);
            unsigned int c = 255 * v;
            dsun.setPixelColor(i, j, QColor(c, c, c));
        }
    dsun.save("diffusesun.png");*/


    delete glsun;
}


void Simulation::reportSunAverages()
{
    int dx, dy;
    float sum;
    ter->getGridDim(dx, dy);

    for (int m = 0; m < 12; m++)
    {
        sum = 0.0f;
        for (int x = 0; x < dx; x++)
            for (int y = 0; y < dy; y++)
            {
                sum += sunlight[m].get(x, y);
            }
        cerr << "sun month " << m << " avg = " << sum / (float)(dx * dy) << endl;
    }
}

void Simulation::addEvent(ScriptedDisturbanceEvent* e)
{
    int l = 0;
    int r = int(eventsQueue.size()) - 1;
    int m = 0;
    while (l <= r) {
        m = (l + r) / 2;
        ScriptedDisturbanceEvent* de = eventsQueue.at(m);

        if (*e < *de) {
            r = m - 1;
        }
        else if (*e > *de) {
            l = m + 1;
        }
        else {
            break;
        }
    }

    if (*e > *eventsQueue.at(m)) ++m; // Insert next to if the element is greater

    if (e->year < this->year || (e->year == this->year && e->month <= e->month + 1))
        eventIndex++;

    eventsQueue.insert(eventsQueue.begin() + m, e);
}

void Simulation::editEvent(int index, ScriptedDisturbanceEvent* e)
{
    if (index < eventsQueue.size()) {
        delete eventsQueue.at(index);
        eventsQueue.at(index) = e;
    }
}

void Simulation::removeEvent(int index)
{
    if (index < eventsQueue.size())
        eventsQueue.erase(eventsQueue.begin() + index);
}

std::vector<DisturbanceEvent> Simulation::getPastNewEvents()
{
    std::vector<DisturbanceEvent> vde;
    while (pastEventIndex < eventIndex) {
        ScriptedDisturbanceEvent* e = eventsQueue[pastEventIndex];
        vde.push_back({e->distEvent->type(), e->year, 0.f});
        pastEventIndex++;
    }
    return vde;
}

int Simulation::rangeSearch(const vpPoint2D& p, float rad,
                            std::list<std::vector<SimLivePlant*>::iterator>& rLivePlants,
                            std::list<std::vector<SimDeadPlant*>::iterator>& rDeadPlants)
{
    int xmin = std::max(int((p.x - rad) / coarseGridSpacing), 0);
    int xmax = std::min(int((p.x + rad) / coarseGridSpacing), coarseGridSizeX - 1);
    int zmin = std::max(int((p.y - rad) / coarseGridSpacing), 0);
    int zmax = std::min(int((p.y + rad) / coarseGridSpacing), coarseGridSizeZ - 1);
    float r2 = rad * rad;
    int n = 0;

    for (int i = xmin; i <= xmax; i++) {
        for (int j = zmin; j <= zmax; j++) {
            int cellId = i * coarseGridSizeZ + j;
            for (auto it : coarseGridLive[cellId]) {
                vpPoint2D nPos((*it)->pos.x, (*it)->pos.z);
                float d2 = (p.x - nPos.x) * (p.x - nPos.x) + (p.y - nPos.y) * (p.y - nPos.y);
                if (d2 <= r2) {
                    rLivePlants.push_back(it);
                    n++;
                }
            }
            for (auto it : coarseGridDead[cellId]) {
                vpPoint2D nPos((*it)->parent->pos.x, (*it)->parent->pos.z);
                float d2 = (p.x - nPos.x) * (p.x - nPos.x) + (p.y - nPos.y) * (p.y - nPos.y);
                if (d2 <= r2) {
                    rDeadPlants.push_back(it);
                    n++;
                }
            }
        }
    }
    return n;
}

void Simulation::rebuildNeighborsAccelerator()
{
    // coarse grid for neighbor queries
    coarseGridSpacing = 10.0f;
    float tx, tz;
    getTerrain()->getTerrainDim(tx, tz);
    coarseGridSizeX = int(std::ceil(tx / coarseGridSpacing));
    coarseGridSizeZ = int(std::ceil(tz / coarseGridSpacing));
    coarseGridLive = std::vector<std::list<std::vector<SimLivePlant*>::iterator>>(coarseGridSizeX * coarseGridSizeZ);
    coarseGridDead = std::vector<std::list<std::vector<SimDeadPlant*>::iterator>>(coarseGridSizeX * coarseGridSizeZ);

    for (auto it = plntpop.begin(); it != plntpop.end(); it++) {
        float x = (*it)->pos.x;
        float z = (*it)->pos.z;
        int i = std::min(int(x / coarseGridSpacing), coarseGridSizeX - 1);
        int j = std::min(int(z / coarseGridSpacing), coarseGridSizeZ - 1);
        coarseGridLive[i * coarseGridSizeZ + j].push_back(it);
    }
    for (auto it = snagpop.begin(); it != snagpop.end(); it++) {
        float x = (*it)->parent->pos.x;
        float z = (*it)->parent->pos.z;
        int i = std::min(int(x / coarseGridSpacing), coarseGridSizeX - 1);
        int j = std::min(int(z / coarseGridSpacing), coarseGridSizeZ - 1);
        coarseGridDead[i * coarseGridSizeZ + j].push_back(it);
    }
}

bool Simulation::writeTemperature(const std::string& filename, Terrain* ter)
{
    int gx, gy;
    ofstream outfile;
    ter->getGridDim(gx, gy);

    outfile.open(filename.c_str(), ios_base::out);
    if (outfile.is_open())
    {
        outfile << gx << " " << gy;
#ifdef STEPFILE
        outfile << " 0.9144"; // hardcoded step
#endif
        outfile << endl;
        for (int y = 0; y < gy; y++)
            for (int x = 0; x < gx; x++)
                for (int m = 0; m < 12; m++)
                    outfile << getTemperature(x, y, m) << " ";
        outfile << endl;
        outfile.close();
        return true;
    }
    else
    {
        cerr << "Error Simulation::writeTemperature:unable to open file " << filename << endl;
        return false;
    }

}

bool Simulation::readClimate(std::string filename)
{
    float elv, val;
    ifstream infile;

    infile.open(filename.c_str(), ios_base::in);
    if (infile.is_open())
    {
        totalrainfall = 0.0f;
        infile >> elv; // elevation at which the temperature readings were taken. Used to calculate a sea-level temperature equivalent

        // temperature values
        for (int m = 0; m < 12; m++)
        {
            infile >> val;
            temperature[m] = val + (elv / 1000.0f * lapserate);
        }

        // sky clarity
        for (int m = 0; m < 12; m++)
        {
            infile >> val;
            cloudiness[m] = 1.0f - val;
        }

        // rainfall
        for (int m = 0; m < 12; m++)
        {
            infile >> val;
            rainfall[m] = val;
            totalrainfall += rainfall[m];
        }

        infile.close();
        return true;
    }
    else
    {
        cerr << "Error Simulation::readClimate: unable to open file" << filename << endl;
        return false;
    }
}

float Simulation::getTemperature(int x, int y, int mth)
{
    return temperature[mth] - ter->getHeight(x, y) / 1000.0f * lapserate;
}


float Simulation::deriveFireMaps(MapFloat& canopyHeight, MapFloat& density, MapFloat& woodMoisture)
{
    float celllen;
    int gx, gy;

    // clear maps and set to simulation grid dimensions
    simcells.getDim(gx, gy);
    canopyHeight.setDim(gx, gy);
    canopyHeight.initMap();
    density.setDim(gx, gy);
    density.initMap();
    woodMoisture.setDim(gx, gy);
    woodMoisture.initMap();

    // calculate length of a sim cell
    celllen = ter->getCellExtent() / (float)sparams.subcellfactor;

    // extract each map in turn
    simcells.extractCanopyHeight(canopyHeight);
    simcells.extractDensity(density, ter, biome);
    simcells.extractWoodMoisture(woodMoisture, ter, this, biome);

    return celllen;
    return 0.0f;
}


float Simulation::deriveWindMaps(MapFloat& canopyHeight, MapFloat& density)
{
    float celllen;
    int gx, gy;

    // clear maps and set to simulation grid dimensions
    simcells.getDim(gx, gy);
    canopyHeight.setDim(gx, gy);
    canopyHeight.initMap();
    density.setDim(gx, gy);
    density.initMap();

    // calculate length of a sim cell
    celllen = ter->getCellExtent() / (float)sparams.subcellfactor;

    // extract each map in turn
    simcells.extractCanopyHeight(canopyHeight);
    simcells.extractDensity(density, ter, biome);

    return celllen;
}

void Simulation::extractCanopyHeight(MapFloat& canopyHeight)
{
    int gx, gy;
    canopyHeight.getDim(gx, gy);
    float tx, ty;
    getTerrain()->getTerrainDim(tx, ty);
    float cellSizeX = tx / gx;
    float cellSizeY = ty / gy;

    for (int x = 0; x < gx; x++)
        for (int y = 0; y < gy; y++)
            canopyHeight.set(x, y, 0.0f);

    for (auto pind = plntpop.begin(); pind != plntpop.end(); pind++) {
        float x = (*pind)->pos.x;
        float y = (*pind)->pos.z;
        float r = (*pind)->canopy;
        float h = (*pind)->height;
        float r2 = r * r;

        int imin = std::max(0, int((x - r) / cellSizeX));
        int jmin = std::max(0, int((y - r) / cellSizeY));
        int imax = std::min(gx - 1, int((x + r) / cellSizeX));
        int jmax = std::min(gy - 1, int((y + r) / cellSizeY));
        for (int i = imin; i <= imax; i++) {
            for (int j = jmin; j <= jmax; j++) {
                float cx = (i + 0.5) * cellSizeX;
                float cy = (j + 0.5) * cellSizeY;
                float d2 = (cx - x) * (cx - x) + (cy - y) * (cy - y);
                if (d2 < r2) {
                    canopyHeight.set(i, j, std::max(h, canopyHeight.get(i, j)));
                }
            }
        }

        int i = std::min(gx - 1, int(x / cellSizeX));
        int j = std::min(gy - 1, int(y / cellSizeY));
        canopyHeight.set(i, j, std::max(h, canopyHeight.get(i, j)));
    }
}

void Simulation::applyFireMap(const MapFloat& fireMap)
{
    std::vector<SimDeadPlant*>::iterator sind;
    std::vector<SimLog*>::iterator lind;
    float terSizeX, terSizeY, h;
    ter->getTerrainDim(terSizeX, terSizeY);
    int ncx, ncy;
    fireMap.getDim(ncx, ncy);
    float dimCellX = terSizeX / ncx;
    float dimCellY = terSizeY / ncy;

    // burn live plants
    for (auto pind = plntpop.begin(); pind != plntpop.end(); pind++) {
        float x = (*pind)->pos.x;
        float y = (*pind)->pos.z;
        float i = x / dimCellX;
        float j = y / dimCellY;
        float dmg = fireMap.getBilinear(i, j);

        if (dmg > 0.0f) {
            recordLiveEvent(pind, DisturbanceEventCategory::FIRE, dmg);

            // a healthy tree with strong carbohydrate reserves has some change of regrowing
            // reduce carbohydrate reserves first then induce stress in proportion to damage
            damage(pind, dmg * sparams.firedmgscale);
        }
    }

    // burn snags by advancing decay
    sind = snagpop.begin();
    while (sind != snagpop.end())
    {
        float x = (*sind)->parent->pos.x;
        float y = (*sind)->parent->pos.z;
        float i = x / dimCellX;
        float j = y / dimCellY;
        float damage = fireMap.getBilinear(i, j);

        if (damage > 0)
        {
            recordDeadEvent(sind, DisturbanceEventCategory::FIRE, damage);

            // adjust decay stage according to damage
            // note that there is no log fall as part of this process, because subsidiary branches are simply assumed to have burnt away
            if (burnSnag(sind, damage)) // delete if totally burnt away
            {
                // reduce subsidiary count in parent
                (*sind)->parent->numsub--;

                // remove from simulation grid
                ter->toGrid((*sind)->parent->pos, x, y, h);
                simcells.uprootSnag((*sind), x, y, ter->toGrid((*sind)->parent->canopy));

                delete (*sind);
                sind = snagpop.erase(sind); // require c++11
                snagpopsize--;
            }
            else
            {
                sind++;
            }
        }
        else
        {
            sind++;
        }
    }

    // burn logs
    lind = logpop.begin();
    while (lind != logpop.end())
    {
        float x = (*lind)->liveparent->pos.x;
        float y = (*lind)->liveparent->pos.z;
        float i = x / dimCellX;
        float j = y / dimCellY;
        float damage = fireMap.getBilinear(i, j);

        if (damage > 0.0f)
        {
            if (burnLog(lind, damage)) // delete if totally burnt away
            {
                // reduce subsidiary count in parent
                (*lind)->liveparent->numsub--;

                // remove from simulation grid
                ter->toGrid((*lind)->liveparent->pos, x, y, h);
                simcells.uprootLog((*lind), x, y);

                delete (*lind);
                lind = logpop.erase(lind); // require c++11
                logpopsize--;
            }
            else
            {
                lind++;
            }
        }
        else
        {
            lind++;
        }
    }

    checkDecayedAway(false);
}

float distanceToSegment(float ax, float ay, float bx, float by, float px, float py)
{
    float sx = bx - ax;
    float sy = by - ay;
    float len = std::sqrt(sx * sx + sy * sy);
    float dx = sx / len;
    float dy = sy / len;

    float proj = (px - ax) * dx + (py - ay) * dy;
    if (proj < 0 || proj > len) {
        // projection outside segment, minimum distance from p to either endpoint
        float da2 = (px - ax) * (px - ax) + (py - ay) * (py - ay);
        float db2 = (px - bx) * (px - bx) + (py - by) * (py - by);
        return std::sqrt(std::min(da2, db2));
    }
    else {
        // projection inside segment, distance from p to projected point q
        float qx = ax + proj * dx;
        float qy = ay + proj * dy;
        return std::sqrt((qx - px) * (qx - px) + (qy - py) * (qy - py));
    }
}

void Simulation::windDamage(const std::vector<SimLivePlant*>::iterator pind, DropEventCategory type, float wdirX, float wdirY)
{
    // TODO: account for terrain?
    float falldirX = wdirX;
    float falldirY = wdirY;
    float dirn = std::atan2(falldirY, falldirX);
    float norm = std::sqrt(falldirX * falldirX + falldirY * falldirY);
    falldirX /= norm;
    falldirY /= norm;

    // plant data
    SimLivePlant* plnt = *pind;
    float p_x = plnt->pos.x;
    float p_y = plnt->pos.z;
    float p_h = plnt->height;
    float p_r = 0.5f * plnt->canopy;
    float p_dbh = plnt->dbh;

    if (type == DropEventCategory::UPROOTING || type == DropEventCategory::TRUNKSPLIT) {

        // impact on main plant
        if (type == DropEventCategory::TRUNKSPLIT)
        {
            // randomly split in first quarter of height
            float splithght = (float)dice->generate() / 1000.0f * (*pind)->height * 0.25f;
            // create trunk log
            logFall(pind, DropEventCategory::TRUNKSPLIT, dirn, (*pind)->height - splithght);

            (*pind)->height = splithght;

            // adjust carbohydrate store to reflect stress. This allows possibility of tree surviving
            damage(pind, sparams.winddmg);
            recordLiveEvent(pind, DisturbanceEventCategory::WIND, 0.75f);
        }
        else if (type == DropEventCategory::UPROOTING)
        {
            // create massive log representing the trunk of the tree
            // TODO account for primary branches as extra separate logs

            // create trunk log
            logFall(pind, DropEventCategory::UPROOTING, dirn, (*pind)->height);

            // kill parent and remove from sim grid
            float x, y, h;
            ter->toGrid((*pind)->pos, x, y, h);
            simcells.uproot((*pind), x, y, ter->toGrid((*pind)->canopy), ter->toGrid((*pind)->root), ter);
            if ((*pind)->height < 0.0f)
                cerr << "UPROOTING neg height " << (*pind)->height << endl;
            // note: this is not added to the snags but goes directly to logs
            // add death event to the plants eventlist
            recordLiveEvent(pind, DisturbanceEventCategory::WIND, 1.0f);
            recordLiveEvent(pind, DisturbanceEventCategory::DEATH, 1.0f);
        }

        // search neighbors for potential damage
        std::list<std::vector<SimLivePlant*>::iterator> neighPlants;
        std::list<std::vector<SimDeadPlant*>::iterator> neighSnags;
        float fallTopX = p_x + p_h * falldirX;
        float fallTopY = p_y + p_h * falldirY;
        float fallCtrX = p_x + 0.5 * p_h * falldirX;
        float fallCtrY = p_y + 0.5 * p_h * falldirY;
        float fallRad = std::max(p_r, 0.5f * p_h);
        rangeSearch(vpPoint2D(fallCtrX, fallCtrY), fallRad, neighPlants, neighSnags);

        for (std::vector<SimLivePlant*>::iterator neigh : neighPlants) {
            // self, ignore
            if ((*neigh)->id == plnt->id) continue;

            // neighbor data
            float n_x = (*neigh)->pos.x;
            float n_y = (*neigh)->pos.z;
            float n_h = (*neigh)->height;
            float n_r = 0.5f * (*neigh)->canopy;

            // distance between neighbor and fall trajectory (a segment from p to p+h)
            float distToFall = distanceToSegment(p_x, p_y, fallTopX, fallTopY, n_x, n_y);

            // compute overlap ratio of intersection diameter with neighbor diameter
            float l = std::max(distToFall - n_r, -p_r);
            float r = std::min(distToFall + n_r, p_r);
            float overlap = std::max(0.0f, (r - l) / (2 * n_r));

            // do damage if overlap
            if (overlap > 0) {
                // we will always break a primary branch if hit by trunk
                // otherwise, it will depend on the amount of overlap (x2 since overlap accounted for diameter)
                bool hitByTrunk = distToFall < p_dbh + n_r;
                if (hitByTrunk || dice->generate() / 1000.0f < 2 * overlap) {
                    float logDir = dirn;
                    float logLen = std::max(std::min(1.0f, 2 * overlap) * n_r, 1.0f);
                    logFall(neigh, DropEventCategory::PRIMARY, logDir, logLen);
                    damage(neigh, sparams.collisiondmg * overlap);
                }

                // record collision
                recordLiveEvent(neigh, DisturbanceEventCategory::COLLISION, overlap);
            }
        }

        for (std::vector<SimDeadPlant*>::iterator neigh : neighSnags) {

            // neighbor data
            float n_x = (*neigh)->parent->pos.x;
            float n_y = (*neigh)->parent->pos.z;
            float n_h = (*neigh)->parent->height * (*neigh)->trunkremains;

            // distance between neighbor and fall trajectory (a segment from p to p+h)
            float distToFall = distanceToSegment(p_x, p_y, fallTopX, fallTopY, n_x, n_y);

            // if snag still has branches, we could break a primary branch
            if ((*neigh)->dc == SnagDecayClass::BARE || (*neigh)->dc == SnagDecayClass::INTACT) {
                float n_r = 0.5f * (*neigh)->parent->canopy;
                logFall(neigh, DropEventCategory::PRIMARY, dirn, n_r); // create log
            }
            // if just the trunk, we could force a transition to split
            else if ((*neigh)->dc == SnagDecayClass::DEBRANCHED) {
                // snap point in [0.2, 0.6] of height
                float snap = 0.4f * (dice->generate() / 1000.0f) + 0.2f;
                float len = (*neigh)->parent->height;
                len *= snap;
                (*neigh)->trunkremains = snap;
                logFall(neigh, DropEventCategory::TRUNKSPLIT, dirn, len);

                (*neigh)->decay = def_snagDC3;
                (*neigh)->decayage = invdecayfn((*neigh)->decay) * (float)(*neigh)->deathspan;
                (*neigh)->deldecay = 0.0f;
                (*neigh)->dc = SnagDecayClass::SNAPPED;
            }
        }
    }
    else if (type == DropEventCategory::PRIMARY) {
        logFall(pind, DropEventCategory::PRIMARY, dirn, p_r); // create log
        damage(pind, sparams.collisiondmg);
        recordLiveEvent(pind, DisturbanceEventCategory::WIND, 0.25f);
    }
}

void Simulation::windDamage(const std::vector<SimDeadPlant*>::iterator dind, DropEventCategory type, float wdirX, float wdirY)
{
      float dirn = std::atan2(wdirY, wdirX);
      float h = (*dind)->parent->height * (*dind)->trunkremains;
      float len;
      bool logfalling = true;

      if (type == DropEventCategory::TRUNKSPLIT) {
          if ((*dind)->dc > SnagDecayClass::DEBRANCHED) // don't split it if it is already split
          {
              logfalling = false;
          }
          else
          {
              len = std::max(h - 1.3f, 0.2f);
              (*dind)->decay = def_snagDC3;
              (*dind)->decayage = invdecayfn((*dind)->decay) * (float)(*dind)->deathspan;
              (*dind)->deldecay = 0.0f;
              (*dind)->dc = SnagDecayClass::SNAPPED;
              recordDeadEvent(dind, DisturbanceEventCategory::WIND, 0.75f);
          }
      }
      else if (type == DropEventCategory::UPROOTING) {
          len = h;
          // set snag to 0 decay, so that it will be subsequently removed
          (*dind)->decay = 0.0f;
          (*dind)->decayage = invdecayfn(0.0f) * (float)(*dind)->deathspan;
          (*dind)->deldecay = 0.0f;
          (*dind)->dc = SnagDecayClass::STUMP;
          recordDeadEvent(dind, DisturbanceEventCategory::WIND, 1.0f);
      }
      else if (type == DropEventCategory::PRIMARY) {
          if ((*dind)->dc > SnagDecayClass::DEBRANCHED) // cannot lose a primary if there are none left
              logfalling = false;
          else {
              len = (*dind)->parent->canopy;
              recordDeadEvent(dind, DisturbanceEventCategory::WIND, 0.25f);
          }
      }
      else {
          return;
      }

      if (logfalling)
      {
          SimLog* log = simcells.createLogFromSnag((*dind), dirn, len, type, ter, this, biome);
          logpop.push_back(log);
          logpopsize++;
          recordDeadLogFall(dind, type, DisturbanceEventCategory::WIND, log->id);
      }
}

void Simulation::recordDisease(Disease* dis)
{
    diseases.push_back(dis);
}

void Simulation::infectPlant(std::vector<SimLivePlant*>::iterator pind, Disease* d)
{
    const int devC = 3;
    const int devS = 4;
    const int devI = 6;

    SimLivePlant* p = *pind;

    Infection* inf = new Infection();
    inf->disease = d;
    inf->stage = InfectionStage::INFECTED_CONTAGIOUS;
    inf->monthsSick = 0;

    if (d->monthsContagious > 0)
        inf->periodContagious = std::max(1, d->monthsContagious + int(2 * devC * (dice->generate() / 1000.0f - 0.5)));
    else
        inf->periodContagious = d->monthsContagious;

    if (d->monthsRecovery > d->monthsContagious)
        inf->periodSick = std::max(inf->periodContagious, d->monthsRecovery + int(2 * devS * (dice->generate() / 1000.0f - 0.5)));
    else if (d->monthsRecovery >= 0)
        inf->periodSick = inf->periodContagious;
    else
        inf->periodSick = d->monthsRecovery;

    if (d->monthsImmune > 0)
        inf->remainingImmunity = std::max(0, d->monthsImmune + int(2 * devI * (dice->generate() / 1000.0f - 0.5)));
    else
        inf->remainingImmunity = d->monthsImmune;

    p->infections[d->diseaseId] = inf;

    d->totalInfected++;
    d->currentInfected++;

    recordLiveEvent(pind, DisturbanceEventCategory::DISEASE, 1);
}

float Simulation::updateInfections(std::vector<SimLivePlant*>::iterator plnt)
{
    SimLivePlant* p = *plnt;

    float damage = 0;
    std::list<int> curedInfections;

    for (auto it = p->infections.begin(); it != p->infections.end(); it++) {

        Infection* inf = it->second;
        inf->monthsSick++;

        switch (inf->stage) {
        case InfectionStage::INFECTED_CONTAGIOUS:
        case InfectionStage::INFECTED_NONCONTAGIOUS:
            damage += inf->disease->severity;

            // plant stops being contagious
            if (inf->stage == InfectionStage::INFECTED_CONTAGIOUS &&
                inf->periodContagious >= 0 && inf->monthsSick >= inf->periodContagious) {
                inf->stage = InfectionStage::INFECTED_NONCONTAGIOUS;
            }

            // plant recovers
            if (inf->stage == InfectionStage::INFECTED_NONCONTAGIOUS &&
                inf->periodSick >= 0 && inf->monthsSick >= inf->periodSick) {
                recordLiveEvent(plnt, DisturbanceEventCategory::DISEASE, 0);
                inf->disease->currentInfected--;
                // immunity or cured
                if (inf->remainingImmunity != 0)
                    inf->stage = InfectionStage::IMMUNE;
                else
                    curedInfections.push_back(it->first);
            }
            break;
        case InfectionStage::IMMUNE:
            if (inf->remainingImmunity > 0)
                inf->remainingImmunity--;
            if (inf->remainingImmunity == 0)
                curedInfections.push_back(it->first);
            break;
        }
    }

    // erase cured infections
    for (int disId : curedInfections) {
        p->infections.erase(p->infections.find(disId));
    }

    return damage;
}

int Simulation::spreadDiseases()
{
    // which species are potentially affected (r > 0)
    std::vector<float> maxSpreadRadius(biome->numPFTypes(), 0.0f);

    // count infected
    int numInfected = 0;
    float largestRadius = 1;
    for (const Disease* d : diseases) {
        numInfected += d->currentInfected;
        maxSpreadRadius[d->targetSpecies] = std::max(maxSpreadRadius[d->targetSpecies], d->spreadRadius);
        largestRadius = std::max(largestRadius, d->spreadRadius);

        //std::cerr << "Disease " << d->diseaseId << ": " << d->currentInfected << std::endl;
    }

    // early exit if no diseases to spread
    if (numInfected <= 0)
        return 0;

    // build grid accelerator for diseases
    float gridSize = largestRadius;
    float tx, tz;
    this->getTerrain()->getTerrainDim(tx, tz);
    int numCellsX = int(std::ceil(tx/gridSize));
    int numCellsZ = int(std::ceil(tz/gridSize));
    int n = 0;

    std::map<int,int> diseaseIdToIndex;
    for (int i = 0; i < diseases.size(); i++) {
        diseaseIdToIndex[diseases[i]->diseaseId] = i;
    }

    std::vector<std::vector<std::vector<const SimLivePlant*> > > infectedGrids(diseases.size(),
        std::vector<std::vector<const SimLivePlant*> > (numCellsX*numCellsZ));
    for (auto it = plntpop.begin(); it != plntpop.end(); it++) {
        for (auto itinf = (*it)->infections.begin(); itinf != (*it)->infections.end(); itinf++) {
            Infection* inf = itinf->second;
            if (inf->stage == InfectionStage::INFECTED_CONTAGIOUS) {
                float x = (*it)->pos.x;
                float z = (*it)->pos.z;
                int i = std::min(int(x/gridSize), numCellsX - 1);
                int j = std::min(int(z/gridSize), numCellsZ - 1);
                infectedGrids[diseaseIdToIndex[inf->disease->diseaseId]][i*numCellsZ + j].push_back(*it);
                n++;
            }
        }
    }
    if (n <= 0) {
        return 0;
    }

    // check if plant gets infected
    int cntInfected = 0;
    for (auto plnt = plntpop.begin(); plnt != plntpop.end(); plnt++) {

        // no diseases targetting this plant species
        SimLivePlant *p = *plnt;
        if (maxSpreadRadius[p->pft] <= 0) continue;

        // for each disease
        for (int di = 0; di < diseases.size(); di++) {

            // disease target not of same species
            Disease* d = diseases[di];
            if (d->targetSpecies != p->pft) 
                continue;

            // check plant is not already infected nor immune to this disease
            if (p->infections.find(d->diseaseId) != p->infections.end()) 
                continue;
                        
            // count infected neighbors of this disease
            float spread2 = d->spreadRadius*d->spreadRadius;
            float x = p->pos.x;
            float z = p->pos.z;
            int imin = std::max(0, int((x - d->spreadRadius)/gridSize));
            int imax = std::min(numCellsX - 1, int((x + d->spreadRadius)/gridSize));
            int jmin = std::max(0, int((z - d->spreadRadius)/gridSize));
            int jmax = std::min(numCellsZ - 1, int((z + d->spreadRadius)/gridSize));
            float avgDist = 0;
            float sumArea = 0;
            int neighCnt = 0;
            for (int i = imin; i <= imax; i++) {
                for (int j = jmin; j <= jmax; j++) {
                    for (const SimLivePlant* infNeigh : infectedGrids[di][i*numCellsZ + j]) {                        
                        // distance
                        float dx = x - infNeigh->pos.x;
                        float dz = z - infNeigh->pos.z;
                        float dist2 = dx*dx + dz*dz;
                        if (dist2 < spread2) {
                            avgDist += std::sqrt(dist2);
                            sumArea += infNeigh->canopy * infNeigh->canopy * PI;
                            neighCnt++;
                        }                        
                    }
                }
            }

            // check for infection
            if (neighCnt > 0) {

                // avg dist modifier to prob
                avgDist /= neighCnt;
                float tdist = 2*avgDist/d->spreadRadius - 1;
                tdist = 1.0f - std::min(1.0f, std::max(0.0f, tdist))*0.5;

                // coverage modifier to prob
                float cover = sumArea / (d->spreadRadius * d->spreadRadius * PI);
                float tdens = 0.1 + 0.9 * cover * cover * (3 - 2 * cover);

                float reduction = tdens;
                float r = dice->generate() / 1000.0f;
                if (r < reduction*d->spreadProbability) {
                    infectPlant(plnt, d);
                    cntInfected++;
                    break;
                }
            }
        }
    }

    return cntInfected;
}

void Simulation::debugLogDiseaseInfo(const Disease* disease, int y, int m) const
{
    std::string fname = "diseaseFiles/disease_" + std::to_string(disease->diseaseId) + "_" +
        std::to_string(y) + "-" + std::to_string(m + 1) + ".txt";
    std::fstream fout;
    fout.open(fname, std::fstream::out);

    for (auto plnt = plntpop.begin(); plnt != plntpop.end(); plnt++) {
        if ((*plnt)->pft == disease->targetSpecies) {
            fout << (*plnt)->id << " ";
            fout << (* plnt)->pos.x << " ";
            fout << (* plnt)->pos.y << " ";
            fout << (* plnt)->pos.z << " ";
            fout << (*plnt)->canopy/2.0 << " ";
            fout << (*plnt)->vigour << " ";            
            fout << ((* plnt)->state == PlantSimState::DEAD ? 0 : 1) << " ";

            if ((*plnt)->infections.find(disease->diseaseId) != (*plnt)->infections.end()) {
                fout << int((*plnt)->infections[disease->diseaseId]->stage);
            }
            else {
                fout << 0;
            }

            fout << std::endl;
        }
    }

    fout.close();
}

void Simulation::beginDrought(const std::vector<float>& precipRates)
{
    std::list<float>::iterator it = droughtModifiers.begin();
    for (int i = 0; i < precipRates.size(); i++) {
        if (it != droughtModifiers.end()) {
            *it = std::min(*it, precipRates[i]);
            it++;
        }
        else {
            droughtModifiers.push_back(precipRates[i]);
        }
    }

    std::cerr << "  monthly drought modifiers:";
    for (float f : droughtModifiers) std::cerr << " " << f;
    std::cerr << std::endl;
}

float ecosim::cmpHeight(PlntInCell i, PlntInCell j)
{
    return i.hght > j.hght;;
}

bool MapSimCell::notInSeedbank(int sbidx, int x, int y)
{
    bool found = false;
    int i = 0;

    while (!found && i < (int)get(x, y)->seedbank.size())
    {
        found = (get(x, y)->seedbank[i] == sbidx);
        i++;
    }

    return !found;
}

SimLivePlant* MapSimCell::createLivePlant(int x, int y, int pft, float viability, PlantSimState simstate, Simulation* sim, Terrain* ter, Biome* biome)
{
    SimLivePlant* sp = new SimLivePlant;
    float sgx, sgy, rndoff;
    float wx, wy, wz;
    int dx, dy, numb;
    float delh, delr, dela, deldbh, delblen;

    ter->getGridDim(dx, dy);
    sp->id = sim->getNewPlantID();
    sp->parent = -1;
    sp->state = simstate;
    sp->age = 0;
    sp->nomage = 0.0f;
    sp->numsub = 0;
    sp->height = 0.0f;
    sp->canopy = 0.0f;
    sp->dbh = 0.0f;
    sp->avgblen = 0.0f;
    sp->numb = 0;

    // random subgrid position within simcell
    rndoff = (float)dice->generate() / 10000.0f;
    sgx = (float)x + rndoff;
    rndoff = (float)dice->generate() / 10000.0f;
    sgy = (float)y + rndoff;
    // convert to terrain position
    sgx /= (float)step; sgy /= (float)step; // terrain grid

    // clamp to terrain grid if necessary
    if (sgx > (float)(dx - 1))
        sgx = (float)(dx - 1);
    if (sgy > (float)(dy - 1))
        sgy = (float)(dy - 1);
    sgx = std::max(0.0f, sgx);
    sgy = std::max(0.0f, sgy);

    sp->pos = ter->toWorld(sgx, sgy, ter->getHeight(sgx, sgy));
    // one years growth in proportion to viability
    for (int m = 0; m < biome->getGrowMonths(pft); m++)
    {
        biome->growth(pft, sp->nomage, std::max(0.1f, viability), 0.5f, dela, delh, delr, deldbh, delblen, numb); // one years growth in proportion to viability
        sp->height += delh;
        sp->nomage += dela;
        sp->canopy += 2.0f * delr;
        sp->dbh += deldbh;
        sp->avgblen += delblen;
        sp->numb = numb;
    }
    // cerr << "hght = " << sp->height << endl;
    if (sp->height > 10.0f || sp->height < 0.0f || sp->dbh < 0.f || sp->dbh > 1.0f)
    {
        cerr << "CREATELIVEPLANT error height or dbh out of bounds: pft = " << pft << " height = " << sp->height << " dbh = " << sp->dbh << endl;
    }

    sp->root = sp->canopy;
    sp->reserves = sim->sparams.reservecapacity;
    sp->stress = 0.0f;
    sp->vigour = 0.0f;
    rndoff = (float)dice->generate() / 10000.0f * 0.4f;
    sp->col = glm::vec4(-0.2f + rndoff, -0.2f + rndoff, -0.2f + rndoff, 1.0f);
    sp->pft = pft;
    sp->water = 0.0f;
    sp->sunlight = 0.0f;
    ter->toGrid(sp->pos, wx, wy, wz);
    if (wx > dx - 1)
        wx = dx - 1;
    if (wy > dy - 1)
        wy = dy - 1;
    sp->gx = wx;
    sp->gy = wy;

    // get(x, y)->growing = false; // no other plants can intersect the trunk
    // get(x, y)->available = false;
    inscribe(sp, sgx, sgy, ter->toGrid(sp->canopy), ter->toGrid(sp->root), false, ter, sim);
    return sp;
}

bool MapSimCell::singleSeed(int x, int y, std::vector<SimLivePlant*>* plntpop, Simulation* sim, Terrain* ter, Biome* biome, std::vector<int>& noseed_count, SimLogging* simlog)
{
    // Build a roulette wheel distribution according to the viability of all individual species with seeds at this position.
    // Select the Sapling accordingly. Note that the total viability is factored in so it is possible for locations with low
    // overall viability that no sapling will sprout.

    std::vector<float> adaptvals(4);
    std::vector<float> seedviability;
    std::vector<float> cumviability;
    std::vector<int> seedspecies;
    float totalstr, maxstr, cumstr, sun, wet, temp, slope, logmult = 1.0f;
    int sidx, dx, dy, tgx, tgy;
    bool found;
    bool seeded = false, onlog = false;

    ter->getGridDim(dx, dy);

    // gather understorey plants
    for (int b = 0; b < (int)get(x, y)->seedbank.size(); b++)
    {
        SubBiome* sb = biome->getSubBiome(get(x, y)->seedbank[b]);
        for (int u = 0; u < (int)sb->understorey.size(); u++)
            seedspecies.push_back(sb->understorey[u]);
        for (int o = 0; o < (int)sb->canopies.size(); o++)
            seedspecies.push_back(sb->canopies[o]);
    }

    // determine viability of plants
    totalstr = 0.0f;
    maxstr = 0.0f;

    // convert from x, y simcell to tx, ty terrain grid
    toTerGrid(x, y, tgx, tgy, sim->getTerrain());

    // average sunlight and moisture over growing season
    sun = get(x, y)->leftoversun / (float)shortgrowmonths;
    wet = get(x, y)->leftoverwet / (float)shortgrowmonths;

    temp = sim->getTemperature(tgx, tgy, shortgrowend); // use end of growing season
    slope = sim->getSlopeMap()->get(tgx, tgy);

    // calculate base probability adjustment if seeding in or near a decaying log
    if ((int)get(x, y)->logs.size() > 0)
    {
        bool aconifer = false;
        int lcount = 0;
        float decayavg = 0.0f;

        // find average log decay
        for (auto log : get(x, y)->logs)
        {
            if (biome->isConifer(log->liveparent->pft))
                aconifer = true;
            decayavg += log->decay;
            lcount++;
        }
        if (lcount > 0)
            decayavg /= (float)lcount;

        if (aconifer) // conifers twice as likely as seed beds
            logmult = 2.0f;
        logmult *= (1.5f - decayavg); // more decayed logs are better as seedbeds
        onlog = true;
    }

    for (int s = 0; s < (int)seedspecies.size(); s++)
    {
        float str, spcadj, spclog = 1.0f;

        str = max(0.0f, biome->viability(seedspecies[s], sun, wet, temp, slope, adaptvals, sim->sparams.viabilityneg, false));
        // The species adjustment modifier (spcadj) is crucial since it factors in the higher seed density of short-lived plants, which
        // mature and therefore seed more rapidly.
        if (str > maxstr)
            maxstr = str;
        spcadj = pow((float)biome->getPFType(seedspecies[s])->maxage, 1.2f);
        str /= spcadj;

        if (onlog) // incorporate species specific on log seeding multiplier
            spclog = logmult * biome->getLogSeeding(seedspecies[s]);
        str *= spclog;
        totalstr += str;
        seedviability.push_back(str);
        cumviability.push_back(totalstr);

        for (int i = 0; i < adaptvals.size(); i++)
        {
            if (adaptvals.at(i) <= 0.0f)
                noseed_count.at(i) += 1;
        }
    }

#ifdef TESTING
    cerr << "SEEDBANK" << endl;
    for (int s = 0; s < (int)seedspecies.size(); s++)
    {
        cerr << "s: " << seedspecies[s] << " v: " << seedviability[s] << endl;
    }
    cerr << "total viability = " << totalstr << endl;
#endif

    // choose particular seed randomly according to viability by roulette wheel selection
    if (totalstr > 0.0f)
    {
        float select;
        // first check random suitability of location using totalstr compared to maxstr
        // plants will thus be more likely to seed in generally more favourable areas

        // select = (float) dice->generate() / 10000.0f * maxstr;
        select = (float)dice->generate() / 10000.0f;
        if (select < maxstr) // passes initial suiability check so seedling will be placed
        {
            select = (float)dice->generate() / 10000.0f * totalstr;
            cumstr = 0.0f; found = false; sidx = 0;

            // find particular seed according to random selection
            while (!found)
            {
                cumstr = cumviability[sidx];
                if (select <= cumstr)
                {
                    found = true;
                }
                else
                {
                    sidx++;
                    if (sidx >= (int)cumviability.size())
                        found = true;
                }
            }

#ifdef TESTING
            cerr << "RND GEN" << endl;
            cerr << "rnd number = " << select << " for entry = " << sidx << endl;
#endif

            // cerr << "   sidx = " << sidx << endl;
            // assign height, random position within cell, and other attribues
            if (sidx < (int)cumviability.size())
            {
                SimLivePlant* sp = createLivePlant(x, y, seedspecies[sidx], seedviability[sidx], PlantSimState::ALIVE, sim, ter, biome);
                plntpop->push_back(sp);
                sim->incrPlantPop(seedspecies[sidx]);

                //int indx = simlog->treemap[seedspecies[sidx]];
                //if (indx > 0)
                //    // add new basal area
                //    simlog->databasal[simlog->numyears][3 * indx] += sim->basalArea(sp->dbh);
                seeded = true;
            }
        }
    }
    seedviability.clear();
    seedspecies.clear();

    return seeded;
}

int MapSimCell::deathSpan(int pft, float diameter, Biome* biome)
{
    int deathspan;
    float dwght, diamcm;

    // capture species specific decay period and modify by diameter
    diamcm = diameter * 100.0f;
    if (diamcm < 1.0f)
        dwght = diamcm * 0.6f;
    else
        dwght = 0.00816f * diamcm + 0.592f;
    deathspan = (int)(dwght * (float)biome->getPFType(pft)->maxdeadage);
    return deathspan;
}

bool MapSimCell::inTrunk(float h, float dbh, vpPoint pos, int x, int y, Terrain* ter, Biome* biome)
{
    // is inside trunk if distance of cell center from plant center is less than diameter at breast height
    vpPoint cellpos = ter->toWorld(convert_to_tergrid((float)x), convert_to_tergrid((float)y), 0.0f);
    vpPoint plntpos = pos;
    plntpos.y = 0.0f; cellpos.y = 0.0f;
    return (plntpos.dist(cellpos) < 0.5f * dbh);
}

void MapSimCell::toTerGrid(int mx, int my, int& tx, int& ty, Terrain* ter) const
{
    int dx, dy;

    tx = mx / step;
    ty = my / step;

    // upper bounds check
    ter->getGridDim(dx, dy);
    tx = std::min(tx, dx - 1);
    ty = std::min(ty, dy - 1);
}

void MapSimCell::initMap()
{
    smap.clear();
    smap.resize(gx * gy);
    for (int c = 0; c < gx * gy; c++) // empty cells are always sorted
    {
        smap[c].canopysorted = true;
        smap[c].rootsorted = true;
        smap[c].growing = true;
        smap[c].available = true;
    }
    resetSeeding();

}

void MapSimCell::resetSeeding()
{
    for (int c = 0; c < gx * gy; c++) // empty cells are always sorted
    {
        smap[c].leftoversun = 0.0f;
        smap[c].leftoverwet = 0.0f;
    }
}

void MapSimCell::delMap()
{
    if (smap.empty()) return;
    for (int c = 0; c < gx * gy; c++) // empty cells are always sorted
    {
        smap[c].canopies.clear();
        smap[c].roots.clear();
        smap[c].seedbank.clear();
    }
    smap.clear();
}



void MapSimCell::clamp(int& x, int& y) const
{
    if (x < 0)
        x = 0;
    if (x > gx - 1)
        x = gx - 1;
    if (y < 0)
        y = 0;
    if (y > gy - 1)
        y = gy - 1;
}

void MapSimCell::inscribe(SimLivePlant* plntidx, float px, float py, float rcanopy, float rroot, bool isStatic, Terrain* ter, Simulation* sim)
{
    float grx, gry, grcanopy, grroot, rmax, distsq, grcanopysq, grrootsq;
    int dx, dy, cx, cy, cmax;

    // assumes that plant index is not already inscribed in the simulation grid
    ter->getGridDim(dx, dy);

    // also convert from diameter to radius

    // convert to sim grid coordinates
    grx = convert(px);
    gry = convert(py);

    // also convert from diameter to radius
    grcanopy = convert(rcanopy / 2.0f);
    grroot = convert(rroot / 2.0f);

    /*
    // square of diagonal extent of half a single cell in terrain coordinates
    float singlecell = ter->longEdgeDist() / (float) gx;
    singlecell = 2.0 * singlecell * singlecell;
    */

    grcanopysq = grcanopy * grcanopy;
    grrootsq = grroot * grroot;
    rmax = fmax(grroot * radius_mult, grcanopy * radius_mult);

    vpPoint pt = ter->toWorld(rcanopy / 2.0f, rcanopy / 2.0f, rcanopy / 2.0f);

    float real_rcanopysq = pt.x * pt.x;
    pt = ter->toWorld(rroot / 2.0f, rroot / 2.0f, rroot / 2.0f);
    float real_rrootsq = pt.x * pt.x;

    float tercellsize = ter->getCellExtent();

    float distsqmax = fmax(real_rcanopysq, real_rrootsq) * radius_mult * radius_mult;
    float distmax = sqrt(distsqmax);
    float distsqmax_error = 2.0f * (distmax + tercellsize) * (distmax + tercellsize);

    cx = (int)grx;
    cy = (int)gry;
    cmax = (int)(rmax + 0.5f);

    if (!ingrid(cx, cy))
    {
        cerr << "plant out of bounds" << grx << ", " << gry << endl;
        cerr << "terrain location = " << px << ", " << py << " with terrain bounds = " << dx << ", " << dy << endl;
        cerr << "max bounds = " << gx << ", " << gy << endl;
    }

    get(cx, cy)->available = false;
    get(cx, cy)->growing = false;

    for (int x = cx - cmax; x <= cx + cmax; x++)
        for (int y = cy - cmax; y <= cy + cmax; y++)
        {
            PlntInCell pic;

            if (ingrid(x, y))		// cell is in grid
            {
                pic.hght = 0.0f; // height will be instantiated later on demand
                pic.plnt = plntidx;

                if (x == cx && y == cy) // always include center
                {
                    get(x, y)->canopies.push_back(pic);
                    if ((int)get(x, y)->canopies.size() > 1)
                        get(x, y)->canopysorted = false;
                    get(x, y)->roots.push_back(pic);
                    if ((int)get(x, y)->roots.size() > 1)
                        get(x, y)->rootsorted = false;
                }
                else
                {
                    float delx = convert_to_tergrid((float)x - grx);
                    float dely = convert_to_tergrid((float)y - gry);
                    vpPoint temp = ter->toWorld(delx, dely, 0.0f);
                    delx = temp.x;
                    dely = temp.z;
                    distsq = delx * delx + dely * dely;

                    assert(distsq < distsqmax_error);
                    // if (distsq > distsqmax) distsq = distsqmax;

                    if (!get(x, y)->available)
                        assert(!get(x, y)->growing);

                    if (distsq <= real_rcanopysq)
                    {
                        get(x, y)->canopies.push_back(pic);
                        if ((int)get(x, y)->canopies.size() > 1)
                            get(x, y)->canopysorted = false;
                    }
                    if (distsq <= real_rrootsq)
                    {
                        get(x, y)->roots.push_back(pic);
                        if ((int)get(x, y)->roots.size() > 1)
                            get(x, y)->rootsorted = false;
                    }

                }
            }
        }
}

void MapSimCell::inscribeSnag(SimDeadPlant* snag, float px, float py, float canopy, Terrain* ter)
{
    float grx, gry, grcanopy, rmax, distsq, grcanopysq;

    // assumes that plant index is not already inscribed in the simulation grid

    // convert to sim grid coordinates
    grx = convert(px);
    gry = convert(py);

    // also convert from diameter to radius
    grcanopy = convert(canopy / 2.0f);
    grcanopysq = grcanopy * grcanopy;
    rmax = grcanopy * radius_mult;

    vpPoint pt = ter->toWorld(canopy / 2.0f, canopy / 2.0f, canopy / 2.0f);
    float real_rcanopysq = pt.x * pt.x;

    float tercellsize = ter->getCellExtent();

    float distsqmax = real_rcanopysq * radius_mult * radius_mult;
    float distmax = sqrt(distsqmax);
    float distsqmax_error = 2.0f * (distmax + tercellsize) * (distmax + tercellsize);

    for (int x = (int)(grx - rmax); x <= (int)(grx + rmax); x++)
        for (int y = (int)(gry - rmax); y <= (int)(gry + rmax); y++)
        {
            if (ingrid(x, y))		// cell is in grid
            {
                float delx = convert_to_tergrid((float)x - grx);
                float dely = convert_to_tergrid((float)y - gry);
                vpPoint temp = ter->toWorld(delx, dely, 0.0f);
                delx = temp.x;
                dely = temp.z;
                distsq = delx * delx + dely * dely;

                assert(distsq < distsqmax_error);
                if (distsq > distsqmax) distsq = distsqmax;

                if (distsq <= real_rcanopysq)
                {
                    get(x, y)->snags.push_back(snag);
                }
            }
        }
}

void MapSimCell::inscribeLog(SimLog* log, float px, float py, Terrain* ter)
{
    bool draw = false;
    vpPoint rect[4];
    Vector dirn;

    // yes, I know this isn't an efficient way of doing it, but sometimes simplicity trumps other considerations

    // convert to sim grid coordinates
    float grx = convert(px);
    float gry = convert(py);
    float grlen = convert(log->len);

    // generate corners of the rectangle representing the projecting of the log onto the terrain
    dirn = Vector(1.0f, 0.0f, 0.0f);
    dirn.mult(log->len);
    dirn.rotate(log->dirn);
    formRect(vpPoint(px, py, 0.0f), log->diam / 2.0f, dirn, rect);

    // test point in the simulation grid to see if they fall inside the rectangle
    for (int x = (int)(grx - grlen); x <= (int)(grx + grlen); x++)
        for (int y = (int)(gry - grlen); y <= (int)(gry + grlen); y++)
        {
            if (ingrid(x, y))		// cell is in grid
            {
                float terx = convert_to_tergrid((float)x);
                float tery = convert_to_tergrid((float)y);
                vpPoint test = ter->toWorld(terx, tery, 0.0f);
                if (insideRect(rect, test))
                {
                    get(x, y)->logs.push_back(log);
                }

            }
        }
}

void MapSimCell::expand(SimLivePlant* plntidx, float px, float py, float prevrcanopy, float prevrroot, float newrcanopy, float newrroot)
{
    float grx, gry, gprevrcanopy, gprevrroot, gnewrcanopy, gnewrroot, rmax, distsq, gprevrcanopysq, gprevrrootsq, gnewrcanopysq, gnewrrootsq;
    int cx, cy, cmax;

    if (prevrcanopy > newrcanopy)
        cerr << "EXPAND: CANOPY INCORRECTLY SHRUNK" << endl;
    // convert to sim grid coordinates
    grx = convert(px);
    gry = convert(py);
    gprevrcanopy = convert(prevrcanopy / 2.0f);
    gprevrroot = convert(prevrroot / 2.0f);
    gnewrcanopy = convert(newrcanopy / 2.0f);
    gnewrroot = convert(newrroot / 2.0f);

    gprevrcanopysq = gprevrcanopy * gprevrcanopy;
    gprevrrootsq = gprevrroot * gprevrroot;
    gnewrcanopysq = gnewrcanopy * gnewrcanopy;
    gnewrrootsq = gnewrroot * gnewrroot;
    rmax = fmax(gnewrroot, gnewrcanopy);

    cx = (int)grx;
    cy = (int)gry;
    cmax = (int)(rmax + 0.5f);

    for (int x = cx - cmax; x <= cx + cmax; x++)
        for (int y = cy - cmax; y <= cy + cmax; y++)
        {
            PlntInCell pic;

            if (ingrid(x, y))
            {
                if (!(x == cx && y == cy)) // always exclude center
                {
                    pic.hght = 0.0f; // height will be instantiated later on demand
                    pic.plnt = plntidx;

                    float delx = (float)x - grx;
                    float dely = (float)y - gry;
                    distsq = delx * delx + dely * dely;

                    if (distsq <= gnewrcanopysq && distsq > gprevrcanopysq)
                    {
                        // do a sanity check to ensure we are not re-adding
                        bool found = false;
                        for (auto it = get(x, y)->canopies.begin(); it != get(x, y)->canopies.end(); it++)
                        {
                            if (it->plnt->id == plntidx->id)
                                found = true;
                        }
                        if (!found)
                        {
                            get(x, y)->canopies.push_back(pic);
                            if ((int)get(x, y)->canopies.size() > 1)
                                get(x, y)->canopysorted = false;
                        }
                    }
                    if (distsq <= gnewrrootsq && distsq > gprevrrootsq)
                    {
                        // do a sanity check to ensure we are not re-adding
                        bool found = false;
                        for (auto it = get(x, y)->roots.begin(); it != get(x, y)->roots.end(); it++)
                        {
                            if (it->plnt->id == plntidx->id)
                                found = true;
                        }
                        if (!found)
                        {
                            get(x, y)->roots.push_back(pic);
                            if ((int)get(x, y)->roots.size() > 1)
                                get(x, y)->rootsorted = false;
                        }
                    }
                }
            }
        }
}

void MapSimCell::uproot(SimLivePlant* plntidx, float px, float py, float rcanopy, float rroot, Terrain* ter)
{
    float grx, gry, grcanopy, grroot, rmax, grcanopysq, grrootsq;
    std::vector<PlntInCell>::iterator cidx, delidx;

    // convert to sim grid coordinates
    grx = convert(px);
    gry = convert(py);
    grcanopy = convert(rcanopy / 2.0f);
    grroot = convert(rroot / 2.0f);

    grcanopysq = grcanopy * grcanopy;
    grrootsq = grroot * grroot;
    rmax = fmax(grroot * radius_mult, grcanopy * radius_mult) + 2.0f;

    for (int x = (int)(grx - rmax); x <= (int)(grx + rmax); x++)
        for (int y = (int)(gry - rmax); y <= (int)(gry + rmax); y++)
        {
            if (ingrid(x, y))		// cell is in grid
            {

                // find and remove plntidx from canopies
                // allow for it to appear multiple times even though it shouldn't
                bool fin = false;
                if (!fin)
                {
                    bool found = false;
                    for (cidx = get(x, y)->canopies.begin(); cidx != get(x, y)->canopies.end(); cidx++)
                        if (cidx->plnt == plntidx)
                        {
                            if (found)
                                cerr << "REMOVE CANOPY: MULTIPLE FINDS ON PLANT CANOPY" << endl;
                            found = true;
                            delidx = cidx;
                        }
                    if (found)
                        get(x, y)->canopies.erase(delidx);
                    else
                        fin = true;
                }


                // find and remove plntidx from roots
                fin = false;
                if (!fin)
                {
                    bool found = false;
                    for (cidx = get(x, y)->roots.begin(); cidx != get(x, y)->roots.end(); cidx++)
                        if (cidx->plnt == plntidx)
                        {
                            if (found)
                                cerr << "REMOVE ROOT: MULTIPLE FINDS ON PLANT CANOPY" << endl;
                            found = true;
                            delidx = cidx;
                        }
                    if (found)
                        get(x, y)->roots.erase(delidx);
                    else
                        fin = true;
                }
            }
        }
}

void MapSimCell::uprootSnag(SimDeadPlant* plntidx, float px, float py, float canopy)
{
    float grx, gry, grcanopy, rmax;
    std::vector<SimDeadPlant*>::iterator cidx, delidx;

    // convert to sim grid coordinates
    grx = convert(px);
    gry = convert(py);
    grcanopy = convert(canopy / 2.0f);
    rmax = grcanopy * radius_mult + 2.0f;

    for (int x = (int)(grx - rmax); x <= (int)(grx + rmax); x++)
        for (int y = (int)(gry - rmax); y <= (int)(gry + rmax); y++)
        {
            if (ingrid(x, y))		// cell is in grid
            {
                // find and remove plntidx from snags
                // assumes it appears only once
                bool found = false;
                for (cidx = get(x, y)->snags.begin(); cidx != get(x, y)->snags.end(); cidx++)
                    if ((*cidx) == plntidx)
                    {
                        if (found)
                            cerr << "REMOVE SNAGS: MULTIPLE FINDS ON SNAGS" << endl;
                        found = true;
                        delidx = cidx;
                    }
                if (found)
                    get(x, y)->snags.erase(delidx);
            }
        }
}

void MapSimCell::uprootLog(SimLog* log, float px, float py)
{
    std::vector<SimLog*>::iterator cidx, delidx;

    // convert to sim grid coordinates
    float grx = convert(px);
    float gry = convert(py);
    float grlen = convert(log->len);

    // check and remove log from inside enclosing bounding box
    for (int x = (int)(grx - grlen); x <= (int)(grx + grlen); x++)
        for (int y = (int)(gry - grlen); y <= (int)(gry + grlen); y++)
        {
            if (ingrid(x, y))		// cell is in grid
            {
                bool found = false;
                for (cidx = get(x, y)->logs.begin(); cidx != get(x, y)->logs.end(); cidx++)
                    if ((*cidx) == log)
                    {
                        if (found)
                            cerr << "REMOVE LOGS: MULTIPLE FINDS ON LOGS" << endl;
                        found = true;
                        delidx = cidx;
                    }
                if (found)
                    get(x, y)->logs.erase(delidx);
            }
        }
}

void MapSimCell::traverse(std::vector<SimLivePlant*>* plntpop, Simulation* sim, Biome* biome, MapFloat* sun, MapFloat* wet, bool seedable, float precipRatio)
{
    int tx, ty;
    int canopycnt = 0;

    // clear sunlight and moisture plant pools
    std::vector<SimLivePlant*>::iterator plnt;
    for (plnt = plntpop->begin(); plnt != plntpop->end(); plnt++)
    {
        if ((*plnt)->state == PlantSimState::ALIVE)
        {
            (*plnt)->sunlight = 0.0f;
            (*plnt)->filterlight = 0.0f;
            (*plnt)->water = 0.0f;
            (*plnt)->sunlightcnt = 0;
            (*plnt)->watercnt = 0;
        }
    }

    for (int x = 0; x < gx; x++)
        for (int y = 0; y < gy; y++)
        {
            // sort canopy trees by height if necessary
            if (!get(x, y)->canopysorted)
            {
                // update heights
                for (int p = 0; p < (int)get(x, y)->canopies.size(); p++)
                    get(x, y)->canopies[p].hght = get(x, y)->canopies[p].plnt->height;
                std::sort(get(x, y)->canopies.begin(), get(x, y)->canopies.end(), cmpHeight);
                get(x, y)->canopysorted = true;
            }

            toTerGrid(x, y, tx, ty, sim->getTerrain());
            float sunlight = sun->get(tx, ty);
            // traverse trees by height supplying and reducing sunlight
            canopycnt += (int)smap[flatten(x, y)].canopies.size();
            for (int p = 0; p < (int)smap[flatten(x, y)].canopies.size(); p++)
            {
                PlantSimState pstate = get(x, y)->canopies[p].plnt->state;

                assert(pstate != PlantSimState::DEAD); // dead plants should already be removed from grid
                if (pstate == PlantSimState::DEAD)
                    cerr << "WARNING: DEAD PLANTS IN CANOPY GRID" << endl;

                if (pstate != PlantSimState::STATIC)
                {
                    // cerr << "sunlight = " << sunlight << " at " << (int) x/step << ", " << (int) y/step << endl;
                    get(x, y)->canopies[p].plnt->sunlight += sunlight;
                    if (sun->get(tx, ty) > 0.0f)
                        get(x, y)->canopies[p].plnt->filterlight += (sun->get(tx, ty) - sunlight) / sun->get(tx, ty);
                    else
                        get(x, y)->canopies[p].plnt->filterlight += 0.0f;
                    get(x, y)->canopies[p].plnt->sunlightcnt++;
                    // reduce sunlight by alpha of current plant. note reciprocal
                    sunlight *= (1.0f - biome->getAlpha(get(x, y)->canopies[p].plnt->pft));
                    // cerr << "sunlight incr " << sunlight << " for " << get(x,y)->canopies[p].plnt->id << endl;
                }
            }

            // add value to sunlight accumulation for later seeding check
            if (seedable)
                get(x, y)->leftoversun += sunlight;

            // sort root trees by height if necessary
            if (!get(x, y)->rootsorted)
            {
                // update heights
                for (int p = 0; p < (int)get(x, y)->roots.size(); p++)
                    get(x, y)->roots[p].hght = get(x, y)->roots[p].plnt->height;
                std::sort(get(x, y)->roots.begin(), get(x, y)->roots.end(), cmpHeight);
                get(x, y)->rootsorted = true;
            }

            float moisture = wet->get(tx, ty) * precipRatio;
            float watershare;
            int livecnt = 0;
            // traverse trees by height (proxy for root depth) supplying and reducing moisture
            for (int p = 0; p < (int)smap[flatten(x, y)].roots.size(); p++)
            {
                PlantSimState pstate = get(x, y)->roots[p].plnt->state;

                if (pstate == PlantSimState::DEAD)
                    cerr << "WARNING: DEAD PLANTS IN ROOT GRID" << endl;

                // plants grab flat min
                watershare = fmin(moisture, sim->sparams.moisturedemand);
                //if(pstate == PlantSimState::STATIC)
                //     watershare *= 0.0f;
                // watershare = fmin(moisture, biome->getMinIdeal   dice = new DiceRoller(0,1000);Moisture(get(x,y)->roots[p].plnt->pft));
                get(x, y)->roots[p].plnt->water += watershare;
                get(x, y)->roots[p].plnt->watercnt++;
                moisture -= watershare;
                if (pstate != PlantSimState::STATIC)
                    livecnt++;
            }
            // remainder spread equally
            if (moisture > 0.0f)
            {
                if (livecnt > 0)
                    watershare = moisture / (float)livecnt;
                else
                    watershare = moisture;

                // add share of leftover water to water accumulator for later seeding check
                if (seedable)
                    // get(x,y)->leftoverwet += watershare;
                    get(x, y)->leftoverwet += moisture; // potential water for seed access

                for (int p = 0; p < (int)smap[flatten(x, y)].roots.size(); p++)
                    if (get(x, y)->roots[p].plnt->state != PlantSimState::DEAD)
                        get(x, y)->roots[p].plnt->water += watershare;
            }
        }
    // cerr << "CANOPY COUNT = " << canopycnt << endl;
}

void MapSimCell::uniformSeedBank(int sub_biome)
{
    // cover entire grid
    for (int x = 0; x < gx; x++)
        for (int y = 0; y < gy; y++)
        {
            // not already in seedbank and in growing region
            if (notInSeedbank(sub_biome, x, y) && get(x, y)->growing)
                get(x, y)->seedbank.push_back(sub_biome);
        }
}

void MapSimCell::seeding(std::vector<SimLivePlant*>* plntpop, int plntpopsize, Simulation* sim, Terrain* ter, Biome* biome, SimLogging* simlog)
{
    int nseed_attempts = 0;
    int nseeded = 0;
    int nseeded_total = 0;
    int size_before = plntpopsize;
    std::vector<int> noseed_count(4, 0);
    std::vector<int> noseed_count_outside(4, 0);

    for (int x = 0; x < gx; x++)
        for (int y = 0; y < gy; y++)
            if (get(x, y)->growing) // is cell in a growth zone
            {
                bool seed;
                // check chance to seed

                if (get(x, y)->logs.size() > 0) // seeding on or near a log
                    seed = dice->generate() < (int)(sim->sparams.logseedprob * 10000.0f);
                else
                    seed = dice->generate() < (int)(sim->sparams.seedprob * 10000.0f);

                if (seed) // Note that this does not absolutely guarantee seeding. Abiotic factors still play a part.
                {
                    bool seeded = singleSeed(x, y, plntpop, sim, ter, biome, noseed_count_outside, simlog);
                    if (seeded)
                        nseeded_total++;
                }
            }
    nseeded_total += nseeded;

    /*
    cerr << "num plants before seeding = " << size_before << std::endl;
    cerr << "num plants after seeding = " << plntpopsize << endl;
    cerr << "Number seeded in target location, attempted: " << nseeded << ", " << nseed_attempts << std::endl;
    cerr << "Total number seeded: " << nseeded_total << std::endl;

    cerr << "No seed count due to sun, moisture, temp slope in target area:  ";
    for (auto &v : noseed_count)
    {
        cerr << v << " ";
    }
    cerr << std::endl;
    */
}

SimDeadPlant* MapSimCell::createDeadPlant(SimLivePlant* parent, Terrain* ter, Simulation* sim, Biome* biome)
{
    SimDeadPlant* dp = new SimDeadPlant;
    float dbh, px, py, h;

    dp->id = sim->getNewPlantID();
    dp->parent = parent;
    dp->age = 0;
    dp->decayage = 0.0f;
    dbh = parent->dbh;
    dp->deathspan = deathSpan(parent->pft, dbh, biome);

    if (dp->deathspan < 0)
    {
        cerr << "NEG DEATHSPAN CreateDeadPlant " << dp->deathspan << endl;
        cerr << "pft = " << parent->pft << " dbh = " << dbh << " parent height = " << parent->height << endl;
        cerr << "id = " << parent->id << " parent nomage = " << parent->nomage << endl;
    }
    dp->decay = 1.0f;
    dp->trunkremains = 1.0f;
    dp->decaystr = 1.0f;
    dp->initialbranches = dp->parent->numb; // same number of branches as parent
    dp->remainingbranches = dp->initialbranches;
    dp->dc = SnagDecayClass::INTACT;
    dp->parent->numsub++;
    ter->toGrid(parent->pos, px, py, h);
    inscribeSnag(dp, px, py, parent->canopy, ter);
    return dp;
}

SimLog* MapSimCell::createLogFromLive(SimLivePlant* parent, float dirn, float len, DropEventCategory category, Terrain* ter, Simulation* sim, Biome* biome)
{
    SimLog* lg = new SimLog;
    float dbh, px, py, h;

    lg->id = sim->getNewPlantID();
    lg->liveparent = parent;
    lg->deadparent = nullptr;
    lg->dirn = dirn;
    lg->len = len;
    lg->age = 0;
    dbh = parent->dbh;
    // now modify according to category
    switch (category)
    {
    case DropEventCategory::UPROOTING:
        // use trunk dbh as is
        break;
    case DropEventCategory::TRUNKSPLIT:
        // use trunk dbh as is
        break;
    case DropEventCategory::PRIMARY:
        dbh *= def_branchratio;
        break;
    default:
        break;
    }

    lg->deathspan = deathSpan(parent->pft, dbh, biome);
    if (lg->deathspan < 0)
        cerr << "NEG DEATHSPAN" << lg->deathspan << endl;
    lg->diam = dbh;
    lg->decay = 1.0f;
    lg->decayage = 0.0f;
    lg->cat = category;
    lg->liveparent->numsub++;
    ter->toGrid(parent->pos, px, py, h); // TO FIX - not subcell accurate but terrain grid accurate
    inscribeLog(lg, px, py, ter);

    return lg;
}

SimLog* MapSimCell::createLogFromSnag(SimDeadPlant* parent, float dirn, float len, DropEventCategory category, Terrain* ter, Simulation* sim, Biome* biome)
{
    SimLog* lg = new SimLog;
    float px, py, h;

    lg->id = sim->getNewPlantID();
    lg->liveparent = parent->parent;
    lg->deadparent = parent;
    lg->dirn = dirn;
    lg->len = len;
    lg->age = 0;
    lg->liveparent->numsub++;

    float dbh = parent->parent->dbh; // diameter of source tree
    // now modify according to category
    switch (category)
    {
    case DropEventCategory::UPROOTING:
    case DropEventCategory::TRUNKSPLIT:
        // use trunk dbh as is
        break;
    case DropEventCategory::PRIMARY:
        dbh *= def_branchratio;
        break;
    default:
        break;
    }
    lg->diam = dbh;

    // log has already partially decayed while it was part of the snag
    lg->deathspan = parent->deathspan;
    if (lg->deathspan < 0)
        cerr << "NEG DEATHSPAN FROM SNAG" << lg->deathspan << endl;
    lg->decay = parent->decay;
    lg->decayage = parent->decayage;
    lg->cat = category;
    ter->toGrid(parent->parent->pos, px, py, h);
    inscribeLog(lg, px, py, ter); // TO FIX - not subcell accurate but terrain grid accurate

    return lg;
}

void MapSimCell::extractCanopyHeight(MapFloat& canopyHeight)
{
    for (int x = 0; x < gx; x++)
        for (int y = 0; y < gy; y++)
        {
            // sort canopy trees by height if necessary
            if (!get(x, y)->canopysorted)
            {
                // update heights
                for (int p = 0; p < (int)get(x, y)->canopies.size(); p++)
                    get(x, y)->canopies[p].hght = get(x, y)->canopies[p].plnt->height;
                std::sort(get(x, y)->canopies.begin(), get(x, y)->canopies.end(), cmpHeight);
                get(x, y)->canopysorted = true;
            }

            if ((int)get(x, y)->canopies.size() > 0) // at least one plant
                canopyHeight.set(x, y, get(x, y)->canopies[0].hght); // first plant is tallest
            else // no plants in cell
                canopyHeight.set(x, y, 0.0f);
        }
}

void MapSimCell::extractDensity(MapFloat& density, Terrain* ter, Biome* biome)
{
    for (int x = 0; x < gx; x++)
        for (int y = 0; y < gy; y++)
        {
            float d = 0.0f;

            // density from live plants
            for (int p = 0; p < (int)smap[flatten(x, y)].canopies.size(); p++)
            {
                SimLivePlant* plnt = get(x, y)->canopies[p].plnt;

                if (plnt->state == PlantSimState::ALIVE)
                {
                    // trunk or canopy area, depending on distance from plant center
                    float h = plnt->height;
                    float dbh = plnt->dbh;

                    if (inTrunk(h, dbh, plnt->pos, x, y, ter, biome))
                    {
                        d += h * def_trunkdensity;
                    }
                    else
                    {
                        d += h * def_canopydensity;
                    }
                }
            }

            // density from snags
            for (int p = 0; p < (int)smap[flatten(x, y)].snags.size(); p++)
            {
                SimDeadPlant* snag = get(x, y)->snags[p];

                // trunk or canopy area, depending on distance from plant center
                float h = snag->parent->height;
                float dbh = snag->parent->dbh;

                if (inTrunk(h, dbh, snag->parent->pos, x, y, ter, biome))
                {
                    d += h * snag->trunkremains * def_trunkdensity;
                }
                else
                {
                    float w;
                    // density depends on decay class
                    switch (snag->dc)
                    {
                    case SnagDecayClass::INTACT:
                        w = 1.0f; // as per normal tree, since the snag has only recently died and twigs and leaf detritus lie below
                        break;
                    case SnagDecayClass::BARE:
                        w = 0.6f; // no leaves and twigs, loss of some secondary branches
                        break;
                    case SnagDecayClass::DEBRANCHED:
                        w = 0.1f; // broken primary brances counted as logs, so avoid double counting
                    case SnagDecayClass::SNAPPED: // no canopy remains for any of the rest of the decay classes
                    case SnagDecayClass::STUMP:
                    case SnagDecayClass::SDCEND:
                        w = 0.0f;
                        break;
                    }
                    d += w * h * def_canopydensity;
                }
            }

            // density from logs
            for (int p = 0; p < (int)smap[flatten(x, y)].logs.size(); p++)
            {
                SimLog* log = get(x, y)->logs[p];
                // use diameter as proxy for height since the log is lying on the ground
                d += log->diam; // could also simulate collapse using decay class but not going to bother for now
            }
            density.set(x, y, d);
        }
}

void MapSimCell::extractWoodMoisture(MapFloat& woodMoisture, Terrain* ter, Simulation* sim, Biome* biome)
{
    for (int x = 0; x < gx; x++)
        for (int y = 0; y < gy; y++)
        {
            float d, ad = 0.0f, wd = 0.0f; // d is actual density, wd is fully wet density

            // moisture contribution from live plants
            for (int p = 0; p < (int)smap[flatten(x, y)].canopies.size(); p++)
            {
                SimLivePlant* plnt = get(x, y)->canopies[p].plnt;

                if (plnt->state == PlantSimState::ALIVE)
                {
                    // trunk or canopy area, depending on distance from plant center
                    float h = plnt->height;
                    float dbh = plnt->dbh;
                    if (inTrunk(h, dbh, plnt->pos, x, y, ter, biome))
                    {
                        d = h * def_trunkdensity;
                    }
                    else
                    {
                        d = h * def_canopydensity;
                    }
                    wd += d;
                    if (plnt->reserves > sim->sparams.reservecapacity || plnt->reserves < 0.0f)
                        cerr << "Error MapSimCell::ExtractWoodMoisture - PLANT RESERVES OUT OF RANGE = " << plnt->reserves << endl;
                    ad += d * (0.5f + 0.5f * plnt->reserves / sim->sparams.reservecapacity); // current health of plant influences dryness, can't use stress because that resets to zero at end of year
                }
            }

            // moisture contribution from snags
            for (int p = 0; p < (int)smap[flatten(x, y)].snags.size(); p++)
            {
                SimDeadPlant* snag = get(x, y)->snags[p];

                // trunk or canopy area, depending on distance from plant center
                float h = snag->parent->height;
                float dbh = snag->parent->dbh;
                if (inTrunk(h, dbh, snag->parent->pos, x, y, ter, biome))
                {
                    d = h * snag->trunkremains * def_trunkdensity;
                }
                else
                {
                    float w;
                    // density depends on decay class
                    switch (snag->dc)
                    {
                    case SnagDecayClass::INTACT:
                        w = 1.0f; // as per normal tree, since the snag has only recently died and twigs and leaf detritus lie below
                        break;
                    case SnagDecayClass::BARE:
                        w = 0.6f; // no leaves and twigs, loss of some secondary branches
                        break;
                    case SnagDecayClass::DEBRANCHED:
                        w = 0.1f; // broken primary brances counted as logs, so avoid double counting
                    case SnagDecayClass::SNAPPED: // no canopy remains for any of the rest of the decay classes
                    case SnagDecayClass::STUMP:
                    case SnagDecayClass::SDCEND:
                        w = 0.0f;
                        break;
                    }
                    d = w * h * def_canopydensity;
                }
                wd += d;
                if (snag->decay > 1.0f || snag->decay < 0.0f)
                    cerr << "ERROR MapSimCell::ExtractWoodMoisture: decay value out of range = " << snag->decay << endl;
                ad += d * 0.5f * snag->decay; // snags get dryer as they decay
            }

            // moisture contribution from logs
            for (int p = 0; p < (int)smap[flatten(x, y)].logs.size(); p++)
            {
                SimLog* log = get(x, y)->logs[p];
                // use diameter as proxy for height since the log is lying on the ground
                d = log->diam; // could also simulate collapse using decay class but not going to bother for now

                wd += d;
                ad += d * std::min(1.0f, get(x, y)->leftoverwet / 120.0f); // amount of left over moisture influences how wet logs are
            }

            if (wd > 0.0f)
                d = ad / wd;
            else
                d = 0.0f; // no wood at all
            woodMoisture.set(x, y, d); // moist density / total density
        }
}

std::vector<float> Simulation::getDataCount(int year, int type) const {
  std::vector<float> v(biome->numPFTypes(), 0.f);

  if (year >= simlog.databasal.size()) return v;

  for (int pft = 0; pft < biome->numPFTypes(); ++pft) {
    int indx = simlog.treemap[pft];
    if (indx != -1) {
      v[pft] = simlog.datacnt[year][3 * indx + type];
    }
  }
  return v;
}

std::vector<float> Simulation::getDataBasal(int year, int type) const {
  std::vector<float> v(biome->numPFTypes(), 0.f);

  if (year >= simlog.databasal.size()) return v;

  for (int pft = 0; pft < biome->numPFTypes(); ++pft) {
    int indx = simlog.treemap[pft];
    if (indx != -1) {
      v[pft] = simlog.databasal[year][3 * indx + type];
    }
  }
  return v;
}


bool Simulation::exportSimJSON(const std::string& filename) const
{
  ofstream outfile;
  int writtenPlants = 0;

  std::vector<int> activeSpecies, liveCount, deadCount, logCount;

  std::cerr << "write sim as JSON (plants) " << filename << endl;
  outfile.open(filename.c_str(), ios_base::out | ios_base::trunc);
  if (not outfile.is_open()) return false;

  outfile.setf(ios::fixed, ios::floatfield);
  outfile.precision(2);  // for fixed format, two decimal places!

  outfile << "{\"pfts\": [";

  for (int pft = 0; pft < biome->numPFTypes(); ++pft) {
    PFType* pftype = biome->getPFType(pft);
    if (pft > 0) {
      outfile << ",";
    }
    outfile << "{" << std::format("\"isTree\":{},\"isConifer\":{}", pftype->isTree, pftype->isConifer) << "}";
  }

  outfile << "],\"plants\":[";
  for (auto plnt = plntpop.begin(); plnt != plntpop.end(); ++plnt) {
    const SimLivePlant* p = *plnt;

    if (plnt != plntpop.begin()) {
      outfile << ",";
    }
    outfile << "{";

    const float x = p->pos.x;
    const float y = p->pos.y;
    const float z = p->pos.z;
    const float height = p->height;
    const float dbh = p->dbh;
    const float canopy = p->canopy;
    const float vigour = p->vigour;
    const int age = p->age;
    const int pft = p->pft;

    outfile << std::format("\"pft\":{},\"age\":{},\"pos\":{{\"x\":{:.3f},\"y\":{:.3f},\"z\":{:.3f}}},\"height\":{:.4f},\"dbh\":{:.4f},\"canopy\":{:.4f},\"vigour\":{:.4f}",
      pft, age, x, y, z, height, dbh, canopy, vigour);

    outfile << "}";
    writtenPlants++;
  }

  outfile << "],\"snags\":[";
  for (auto plnt = snagpop.begin(); plnt != snagpop.end(); ++plnt) {
    const SimDeadPlant* p = *plnt;

    if (plnt != snagpop.begin()) {
      outfile << ",";
    }
    outfile << "{";

    const float x = p->parent->pos.x;
    const float y = p->parent->pos.y;
    const float z = p->parent->pos.z;
    const float height = p->parent->height * p->trunkremains;
    const float dbh = p->parent->dbh;
    const float canopy = p->parent->canopy;
    const float decay = std::min(std::max(p->decay, 0.0f), 1.0f);
    const int age = p->age;
    const int pft = p->parent->pft;

    outfile << std::format("\"pft\":{},\"age\":{},\"pos\":{{\"x\":{:.3f},\"y\":{:.3f},\"z\":{:.3f}}},\"height\":{:.4f},\"dbh\":{:.4f},\"canopy\":{:.4f},\"decay\":{:.2f}",
      pft, age, x, y, z, height, dbh, canopy, decay);

    outfile << "}";
    writtenPlants++;
  }

  outfile << "],\"logs\":[";
  for (auto plnt = logpop.begin(); plnt != logpop.end(); ++plnt) {
    const SimLog* p = *plnt;

    if (plnt != logpop.begin()) {
      outfile << ",";
    }
    outfile << "{";

    const float x = p->liveparent->pos.x;
    const float y = p->liveparent->pos.y;
    const float z = p->liveparent->pos.z;
    const int pft = p->liveparent->pft;
    const int age = p->age;
    const float diam = p->diam;
    const float len = p->len;
    const float dirn = p->dirn;
    const float decay = std::min(std::max(p->decay, 0.0f), 1.0f);
    const int cat = int(p->cat);

    outfile << std::format("\"pft\":{},\"cat\":{},\"age\":{},\"pos\":{{\"x\":{:.3f},\"y\":{:.3f},\"z\":{:.3f}}},\"length\":{:.4f},\"diam\":{:.4f},\"dirn\":{:.4f},\"decay\":{:.2f}",
      pft, cat, age, x, y, z, len, diam, dirn, decay);

    outfile << "}";
    writtenPlants++;
  }

  outfile << "]}";
  outfile.close();
  std::cerr << "num written plants = " << writtenPlants << endl;
  return true;
}

bool Simulation::exportSimPDBTerse(const std::string& filename) const
{
  ofstream outfile;
  int writtenPlants = 0;
  std::vector<int> activeSpecies, liveCount, deadCount, logCount;

  cerr << "write sim as PDB (terse) " << filename << endl;
  outfile.open((char*)filename.c_str(), ios_base::out | ios_base::trunc);

  if (outfile.is_open())
  {
    outfile.setf(ios::fixed, ios::floatfield);
    outfile.precision(3);  

    outfile << "v2.3" << endl;

    // collate species that have a presence in the simulation
    for (int s = 0; s < biome->numPFTypes(); s++)
    {
      int liveNum = 0, deadNum = 0, logNum = 0;

      // find if species has a presence
      for (auto plnt = plntpop.begin(); plnt != plntpop.end(); plnt++)
        if ((*plnt)->pft == s)
          liveNum++;

      for (auto zombie = snagpop.begin(); zombie != snagpop.end(); zombie++)
      {
        if ((*zombie)->parent->pft == s)
          deadNum++;
      }

      for (auto log = logpop.begin(); log != logpop.end(); log++)
      {
        if ((*log)->liveparent->pft == s)
          logNum++;
      }

      if (liveNum > 0)
      {
        activeSpecies.push_back(s);
        liveCount.push_back(liveNum);
        deadCount.push_back(deadNum);
        logCount.push_back(logNum);
      }
    }

    outfile << (int)activeSpecies.size() << endl;

    for (int k = 0; k < (int)activeSpecies.size(); k++)
    {
      int s;
      std::vector<Plant> tpop; tpop.clear();

      s = activeSpecies[k];
      cerr << "activeSpecies = " << s << " ";
      outfile << s << " N/A " << endl;
      outfile << liveCount[k] << " " << deadCount[k] << " " << logCount[k] << endl;
      cerr << "liveCount = " << liveCount[k] << " deadCount = " << deadCount[k] << " logCount = " << logCount[k] << endl;

      for (auto plnt = plntpop.begin(); plnt != plntpop.end(); plnt++)
        if ((*plnt)->pft == s)
        {
          long id = (*plnt)->id;
          long pid = (*plnt)->parent;
          int a = (*plnt)->age;
          float x = (*plnt)->pos.x;
          float y = (*plnt)->pos.z;
          float z = (*plnt)->pos.y;

          outfile << id << " " << pid << " " << a << " " << x << " " << y << " " << z << " ";

          float h = (*plnt)->height;
          float r = (*plnt)->canopy / 2.0f;
          float v = (*plnt)->vigour;
          float d = (*plnt)->dbh;
          outfile << h << " " << r << " " << v << " "; // << d << " ";  Oscar: used during debug

          // dead or alive?
          if ((*plnt)->state == PlantSimState::DEAD)
            outfile << 0 << " ";
          else
            outfile << 1 << " ";

          // last entry of growth history
          // int hidx = (int) (* plnt)->hghthistory.size()-1;
          // outfile << (* plnt)->vigourhistory[hidx] << " " << (* plnt)->hghthistory[hidx] << " " << (* plnt)->radhistory[hidx] << endl;

          // disturbance events
          outfile << (int)(*plnt)->eventhistory.size() << " ";
          for (int e = 0; e < (int)(*plnt)->eventhistory.size(); e++)
            outfile << (int)(*plnt)->eventhistory[e].cat << " " << (*plnt)->eventhistory[e].year << " " << (*plnt)->eventhistory[e].severity << " ";
          outfile << endl;

          writtenPlants++;
        }

      for (auto zombie = snagpop.begin(); zombie != snagpop.end(); zombie++)
        if ((*zombie)->parent->pft == s)
        {
          long id = (*zombie)->id;
          long pid = (*zombie)->parent->id;
          int a = (*zombie)->age;
          float x = (*zombie)->parent->pos.x;
          float y = (*zombie)->parent->pos.z;
          float z = (*zombie)->parent->pos.y;

          outfile << id << " " << pid << " " << a << " " << x << " " << y << " " << z << " ";

          float h = (*zombie)->parent->height * (*zombie)->trunkremains;
          float r = (*zombie)->parent->canopy / 2.0f;
          float d = (*zombie)->decay;
          d = std::min(d, 1.0f);
          d = std::max(d, 0.0f);
          outfile << h << " " << r << " " << d << " ";

          // maximum height that tree reached when it was alive
          float maxh = (*zombie)->parent->height;
          outfile << maxh << " ";

          // disturbance events
          outfile << (int)(*zombie)->eventhistory.size() << " ";
          for (int e = 0; e < (int)(*zombie)->eventhistory.size(); e++)
            outfile << (int)(*zombie)->eventhistory[e].cat << " " << (*zombie)->eventhistory[e].year << " " << (*zombie)->eventhistory[e].severity << " ";
          outfile << endl;

          writtenPlants++;
        }

      for (auto log = logpop.begin(); log != logpop.end(); log++)
      {
        bool isSpecies = false;
        long pid;

        if ((*log)->liveparent->pft == s)
        {
          isSpecies = true;
          if ((*log)->deadparent != nullptr)
            pid = (*log)->deadparent->id;
          else
            pid = (*log)->liveparent->id;
        }

        if (isSpecies)
        {
          outfile << (*log)->id << " " << pid << " " << (*log)->age << " " << (*log)->decay << " "
            << (*log)->dirn << " " << (*log)->len << " " << (*log)->diam << " ";
          float x = (*log)->liveparent->pos.x;
          float y = (*log)->liveparent->pos.z;
          float z = (*log)->liveparent->pos.y;
          outfile << x << " " << y << " " << z << " " << (int)(*log)->cat << endl;
          writtenPlants++;
        }
      }
    }
    outfile.close();
  }
  else
  {
    cerr << "Error Simulation::writeSim: unable to open " << filename << endl;
    return false;
  }
  cerr << "num written plants = " << writtenPlants << endl;
  return true;
}

bool Simulation::exportSimPDBFull(const std::string& filename) const
{
  ofstream outfile;
  int writtenPlants = 0;
  std::vector<int> activeSpecies, liveCount, deadCount, logCount;

  cerr << "write sim as PDB (full) " << filename << endl;
  outfile.open((char*)filename.c_str(), ios_base::out | ios_base::trunc);
  if (outfile.is_open())
  {
    outfile.setf(ios::fixed, ios::floatfield);
    outfile.precision(3);

    outfile << "v2.1" << endl;

    // collate species that have a presence in the simulation
    for (int s = 0; s < biome->numPFTypes(); s++)
    {
      int liveNum = 0, deadNum = 0, logNum = 0;

      // find if species has a presence
      for (auto plnt = plntpop.begin(); plnt != plntpop.end(); plnt++)
        if ((*plnt)->pft == s)
          liveNum++;

      for (auto zombie = snagpop.begin(); zombie != snagpop.end(); zombie++)
      {
        if ((*zombie)->parent->pft == s)
          deadNum++;
      }

      for (auto log = logpop.begin(); log != logpop.end(); log++)
      {
        if ((*log)->liveparent->pft == s)
          logNum++;
      }

      if (liveNum > 0)
      {
        activeSpecies.push_back(s);
        liveCount.push_back(liveNum);
        deadCount.push_back(deadNum);
        logCount.push_back(logNum);
      }
    }

    outfile << (int)activeSpecies.size() << endl;

    for (int k = 0; k < (int)activeSpecies.size(); k++)
    {
      int s;
      std::vector<Plant> tpop; tpop.clear();

      s = activeSpecies[k];
      cerr << "activeSpecies = " << s << " ";
      outfile << s << " N/A " << endl;
      outfile << liveCount[k] << " " << deadCount[k] << " " << logCount[k] << endl;
      cerr << "liveCount = " << liveCount[k] << " deadCount = " << deadCount[k] << " logCount = " << logCount[k] << endl;

      for (auto plnt = plntpop.begin(); plnt != plntpop.end(); plnt++)
        if ((*plnt)->pft == s)
        {
          long id = (*plnt)->id;
          long pid = (*plnt)->parent;
          int a = (*plnt)->age;
          float x = (*plnt)->pos.x;
          float y = (*plnt)->pos.z;
          float z = (*plnt)->pos.y;

          outfile << id << " " << pid << " " << a << " " << x << " " << y << " " << z << " ";

          float h = (*plnt)->height;
          float r = (*plnt)->canopy / 2.0f;
          float v = (*plnt)->vigour;
          outfile << h << " " << r << " " << v << " ";

          // plants growth history
          outfile << (int)(*plnt)->hghthistory.size() << " ";
          for (int h = 0; h < (int)(*plnt)->hghthistory.size(); h++)
            outfile << (*plnt)->vigourhistory[h] << " " << (*plnt)->hghthistory[h] << " " << (*plnt)->radhistory[h] << " ";
          outfile << endl;

          // disturbance events
          outfile << (int)(*plnt)->eventhistory.size() << " ";
          for (int e = 0; e < (int)(*plnt)->eventhistory.size(); e++)
            outfile << (int)(*plnt)->eventhistory[e].cat << " " << (*plnt)->eventhistory[e].year << " " << (*plnt)->eventhistory[e].severity << " ";
          outfile << endl;

          // log events
          outfile << (int)(*plnt)->drophistory.size() << " ";
          for (int d = 0; d < (int)(*plnt)->drophistory.size(); d++)
            outfile << (int)(*plnt)->drophistory[d].cat << " " << (int)(*plnt)->drophistory[d].cause << " " << (*plnt)->drophistory[d].year << " " << (*plnt)->drophistory[d].pieceID << " ";
          outfile << endl;

          writtenPlants++;
        }

      for (auto zombie = snagpop.begin(); zombie != snagpop.end(); zombie++)
        if ((*zombie)->parent->pft == s)
        {
          long id = (*zombie)->id;
          long pid = (*zombie)->parent->id;
          int a = (*zombie)->age;
          float x = (*zombie)->parent->pos.x;
          float y = (*zombie)->parent->pos.z;
          float z = (*zombie)->parent->pos.y;

          outfile << id << " " << pid << " " << a << " " << x << " " << y << " " << z << " ";

          float h = (*zombie)->parent->height * (*zombie)->trunkremains;
          float r = (*zombie)->parent->canopy / 2.0f;
          float d = (*zombie)->decay;
          d = std::min(d, 1.0f);
          d = std::max(d, 0.0f);
          outfile << h << " " << r << " " << d << " ";

          // if(h > 2.0f)
          // {
          outfile << (int)(*zombie)->decayhistory.size() << " ";
          for (int h = 0; h < (int)(*zombie)->decayhistory.size(); h++)
          {
            d = (*zombie)->decayhistory[h];
            d = std::min(d, 1.0f);
            d = std::max(d, 0.0f);
            outfile << d << " ";
            // << (* zombie)->hghthistory[h] << " " << (* zombie)->radhistory[h] << " ";
          }
          outfile << endl;

          // disturbance events
          outfile << (int)(*zombie)->eventhistory.size() << " ";
          for (int e = 0; e < (int)(*zombie)->eventhistory.size(); e++)
            outfile << (int)(*zombie)->eventhistory[e].cat << " " << (*zombie)->eventhistory[e].year << " " << (*zombie)->eventhistory[e].severity << " ";
          outfile << endl;

          // log events
          outfile << (int)(*zombie)->drophistory.size() << " ";
          for (int d = 0; d < (int)(*zombie)->drophistory.size(); d++)
            outfile << (int)(*zombie)->drophistory[d].cat << " " << (int)(*zombie)->drophistory[d].cause << " " << (*zombie)->drophistory[d].year << " " << (*zombie)->drophistory[d].pieceID << " ";
          outfile << endl;

          writtenPlants++;
        }

      for (auto log = logpop.begin(); log != logpop.end(); log++)
      {
        bool isSpecies = false;
        long pid;

        if ((*log)->liveparent->pft == s)
        {
          isSpecies = true;
          if ((*log)->deadparent != nullptr)
            pid = (*log)->deadparent->id;
          else
            pid = (*log)->liveparent->id;
        }

        if (isSpecies)
        {
          outfile << (*log)->id << " " << pid << " " << (*log)->age << " " << (*log)->decay << " " << (*log)->dirn << " " << (*log)->len << " ";
          outfile << endl;
          writtenPlants++;
        }
      }
    }
    outfile.close();
  }
  else
  {
    cerr << "Error Simulation::writeSim: unable to open " << filename << endl;
    return false;
  }
  cerr << "num written plants = " << writtenPlants << endl;
  return true;
}
