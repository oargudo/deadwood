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


#include "grass.h"
#include "dice_roller.h"

bool GrassFloat::read(std::string filename)
{
    float val;
    ifstream infile;

    infile.open((char *) filename.c_str(), ios_base::in);
    if(infile.is_open())
    {
        infile >> gx >> gy;
        initMap();

        for (int x = 0; x < gx; x++)
        {
            for (int y = 0; y < gy; y++)
            {
                infile >> val;
                set(x, y, val);
            }
        }
        infile.close();
        return true;
    }
    else
    {
        cerr << "Error TypeMap::loadTxt: unable to open file" << filename << endl;
        return false;
    }
}

float GrassFloat::getBilinear(float x, float y) const
{
    int i0 = std::max(0, std::min(gx-1, int(x)));
    int j0 = std::max(0, std::min(gy-1, int(y)));
    int i1 = std::max(0, std::min(gx-1, i0+1));
    int j1 = std::max(0, std::min(gy-1, j0+1));
    float fx = std::max(0.0f, x - float(i0));
    float fy = std::max(0.0f, y - float(j0));

    float v00 = get(i0, j0);
    float v01 = get(i0, j1);
    float v10 = get(i1, j0);
    float v11 = get(i1, j1);
    float v0  = fy*v01 + (1 - fy)*v00;
    float v1  = fy*v11 + (1 - fy)*v10;
    float v   = fx*v1  + (1 - fx)*v1;
    return v;
}
