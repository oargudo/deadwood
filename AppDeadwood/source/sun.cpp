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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sun.h"
#include "dice_roller.h"

#define GLM_FORCE_RADIANS
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/rotate_vector.hpp>
using namespace std;
using namespace ecosim;

#define AXIS_TILT 0.408407f
const float SunLight::_axis_tilt = AXIS_TILT;
const float SunLight::_monthly_axis_tilt = (AXIS_TILT)/3.f;
const float SunLight::_half_day_in_minutes = 720.0f;
const float SunLight::_quarter_day_in_minutes = SunLight::_half_day_in_minutes/2.0f;
const float SunLight::_3_quarters_day_in_minutes = SunLight::_half_day_in_minutes + SunLight::_quarter_day_in_minutes;
const float avgtransmit = 0.45f; // proportion of light blocked by leaves


////
// SunScene
////

SunScene::SunScene()
{
    view = new View();
    terrain = new Terrain();
    terrain->initGrid(1024, 1024, 10000.0f, 10000.0f);
    view->setForcedFocus(terrain->getFocus());
    view->setViewScale(terrain->longEdgeDist());
    view->setViewType(ViewState::ORTHOGONAL);

    float tx, ty;
    terrain->getTerrainDim(tx, ty);
    view->setOrthoViewExtent(tx,ty);
}

SunScene::~SunScene()
{
    delete view;
    //if (alpha) delete alpha;
}

////
// GLSun
////
//
GLSun::GLSun() 
{
    scene = new SunScene();
}

GLSun::~GLSun()
{
    if (scene) delete scene;
    if (fbo) delete fbo;
    deleteTerrainOpenGLbuffers();
}

void GLSun::colToCoord(QColor col, int &gx, int &gy)
{
    int r, g, b, idx, dx, dy;

    getTerrain()->getGridDim(dx, dy);
    col.getRgb(&r, &g, &b);

    // assumes 8-bits per colour channel
    idx = (r * 65536) + (g * 256) + b;

    // then derive grid coordinates
    //gx = (int) (idx / dy);
    //gy = idx - gx * dy;
    gx = idx % dx;
    gy = idx / dx;
}

View * GLSun::getView()
{
    return scene->view;
}

Terrain * GLSun::getTerrain()
{
    return scene->terrain;
}

//MapFloat * GLSun::getCanopyHeight()
//{
//    return scene->chght;
//}
//
//MapFloat * GLSun::getCanopyDensity()
//{
//    return scene->cdense;
//}
//
//MapFloat * GLSun::getAlpha()
//{
//    return scene->alpha;
//}

//PMrender::TRenderer * GLSun::getRenderer()
//{
//    return renderer;
//}

void GLSun::setScene(Terrain * ter/*, MapFloat* ch, MapFloat* cd*/)
{
    float tx, ty;
    int dx, dy;

    scene->terrain = ter;
    //scene->chght = ch;
    //scene->cdense = cd;

    getView()->setForcedFocus(getTerrain()->getFocus());
    getView()->setViewScale(getTerrain()->longEdgeDist());
    getTerrain()->calcMeanHeight();
    getTerrain()->getTerrainDim(tx, ty);
    getView()->setOrthoViewExtent(tx,ty);
    getTerrain()->getGridDim(dx, dy);
    //if (scene->alpha == nullptr)
    //    scene->alpha = new MapFloat();
    //scene->alpha->setDim(dx, dy);

    updateHeightMap(dx, dy, tx, ty, getTerrain()->getGridPtr(), true);
}
//
//void GLSun::deriveAlpha(EcoSystem *eco, Biome * biome)
//{
//    getAlpha()->fill(0.0f);
//    eco->sunSeeding(getTerrain(), biome, getAlpha());
//}

//void GLSun::alphaMapStats()
//{
//    int dx, dy;
//    float avgalpha = 0.0f, a;
//    int cntalpha = 0;
//
//    cerr << "**** ALPHA MAP STATS ***" << endl;
//    if(getAlpha() != nullptr)
//    {
//        getAlpha()->getDim(dx, dy);
//        cerr << "Alpha dim = " << dx << ", " << dy << endl;
//        for(int x = 0; x < dx; x++)
//            for(int y = 0; y < dy; y++)
//            {
//                a = getAlpha()->get(x, y);
//                if(a > 0.0f)
//                {
//                    avgalpha += a;
//                    cntalpha++;
//                }
//            }
//        cerr << "Average alpha = " << avgalpha / (float) cntalpha << " with nonzero of " << cntalpha << " from " << dx*dy << endl;
//        cerr << "with percentage coverage = " << (float) cntalpha / (float) (dx*dy) << endl;
//    }
//    else
//    {
//        cerr << "alpha map does not exist" << endl;
//    }
//}

//void GLSun::bind()
//{
//    if(sun)
//        delete sun;
//    sun = new CanopyShape(0.75, getTerrain()->getCellExtent());
//}

void GLSun::initializeGL()
{
    QColor col = QColor::fromCmykF(0.0, 0.0, 0.0, 0.0).lighter(); //TODO: check
    glClearColor(0, 0, 0, 0);

    //int mu;
    //glGetIntegerv(GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS, &mu);
    //cerr << "max texture units = " << mu << endl;

    // set up light
    Vector dl = Vector(0.6f, 1.0f, 0.6f);
    dl.normalize();

    // initialise renderer/compile shaders
    initShaders();

    htmapTexUnit = GL_TEXTURE0;
    normalMapTexUnit = GL_TEXTURE1;

    fbo = new QOpenGLFramebufferObject(fboSize, QOpenGLFramebufferObject::CombinedDepthStencil);
}

void GLSun::initShaders()
{
    if (shadersReady) return; // already compiled

    QString shadersPath = ":/shaders/";

    QOpenGLShaderProgram* s = new QOpenGLShaderProgram();
    s->addShaderFromSourceFile(QOpenGLShader::Vertex,   shadersPath + "sun.vert");
    s->addShaderFromSourceFile(QOpenGLShader::Fragment, shadersPath + "sun.frag");
    s->link();
    shaders["sunShader"] = s;

    s = new QOpenGLShaderProgram();
    s->addShaderFromSourceFile(QOpenGLShader::Vertex, shadersPath + "genNormal.vert");
    s->addShaderFromSourceFile(QOpenGLShader::Fragment, shadersPath + "genNormal.frag");
    s->link();
    shaders["normalShader"] = s;

    s = new QOpenGLShaderProgram();
    s->addShaderFromSourceFile(QOpenGLShader::Vertex, shadersPath + "canopy.vert");
    s->addShaderFromSourceFile(QOpenGLShader::Fragment, shadersPath + "canopy.frag");
    s->link();
    shaders["canopyShader"] = s;

    shadersReady = true;
}

void GLSun::drawSun(View* view, int renderPass)
{
    // Save the current state
    GLboolean depthTestEnabled = glIsEnabled(GL_DEPTH_TEST);
    GLboolean depthMaskEnabled;
    glGetBooleanv(GL_DEPTH_WRITEMASK, &depthMaskEnabled);
    GLint depthFunc;
    glGetIntegerv(GL_DEPTH_FUNC, &depthFunc);
    GLfloat depthRange[2];
    glGetFloatv(GL_DEPTH_RANGE, depthRange);
    GLboolean cullFaceEnabled = glIsEnabled(GL_CULL_FACE);

    // render at 2X resolution for later linear downsampling (basic anti-aliasing)
    //updateRadianceScalingBuffers(2 * viewport[2], 2 * viewport[3]);

    // Set the clear color to white
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f); 

    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE); 
    glDepthFunc(GL_LEQUAL); 
    glDepthRange(0.0f, 1.0f);
    glEnable(GL_CULL_FACE);

    if (renderPass == 2) // enable blending
    {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    CE();

    // configure shading params.


    // ************************* Terrain Setup & render code *****************************************
    glm::mat3x3 normalMatrix = view->getNormalMtx();
    glm::mat4x4 MVP = view->getMatrix();
    glm::mat4x4 MVmx = view->getViewMtx();
    glm::mat4x4 projMx = view->getProjMtx();

    std::string shaderName;
    shaderName = "sunShader";
    GLuint programID = (*shaders[shaderName]).programId();
    glUseProgram(programID); CE();

    glUniform1i(glGetUniformLocation(programID, "drawCanopies"), renderPass - 1); CE(); // whether or not to use indexed colours

    // PHASE 1: Draw Terrain
    glUniformMatrix3fv(glGetUniformLocation(programID, "normMx"), 1, GL_FALSE, &normalMatrix[0][0]); CE();
    glUniformMatrix4fv(glGetUniformLocation(programID, "MV"), 1, GL_FALSE, &MVmx[0][0]); CE();
    glUniformMatrix4fv(glGetUniformLocation(programID, "MVproj"), 1, GL_FALSE, &MVP[0][0]); CE();

    glUniform2f(glGetUniformLocation(programID, "terdim"), (float)(twidth), (float)(theight)); CE();

    // pass height and normal map to shader
    GLuint textur = glGetUniformLocation(programID, "htMap"); CE();
    glUniform1i(textur, (GLint)(htmapTexUnit - GL_TEXTURE0));  CE(); // assumes texture unit 0 is bound to heightmap texture

    GLuint textur2 = glGetUniformLocation(programID, "normalMap"); CE();
    glUniform1i(textur2, (GLint)(normalMapTexUnit - GL_TEXTURE0)); CE(); // assumes texture unit 1 is bound to normal map texture

    // draw terrain:
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); CE();
    glUniform1i(glGetUniformLocation(programID, "drawWalls"), 0); CE(); // do NOT draw walls


    glBindVertexArray(vaoTerrain); CE();
    glDrawElements(GL_TRIANGLE_STRIP, indexSize, GL_UNSIGNED_INT, 0); CE();

    // PHASE 2: Draw capping walls
    // --- move base up/down to avoid edits to surface which 'punch through' base

    glUniform1f(glGetUniformLocation(programID, "terrainBase"), terrainBase); CE();
    glUniform1f(glGetUniformLocation(programID, "terrainBasePad"), terrainBasePad); CE();
    glUniform1i(glGetUniformLocation(programID, "drawWalls"), 1); // draw walls - ignore normal map lookup
    GLuint loc = glGetUniformLocation(programID, "normalWall");

    // turn off region layer texturing for walls:
    glUniform1i(glGetUniformLocation(programID, "useRegionTexture"), 0); CE(); // always off for walls!

    for (int i = 0; i < 5; i++)
    {
        glUniform3fv(loc, 1, &normalWalls[i][0]); CE();
        glBindVertexArray(vaoWalls[i]); CE();
        glDrawElements(GL_TRIANGLE_STRIP, wallDrawEls[i], GL_UNSIGNED_INT, 0); CE();
    }

    // PHASE 3: Draw canopies
    // draw alpha blended canopies
    if (renderPass == 2)
    {
        // disable depth writes so that more than one canopy can block sunlight
        glDepthMask(GL_FALSE); CE();

        programID = (*shaders["canopyShader"]).programId();
        //drawManipulators(programID);
        // revert to previous settings
        glDepthMask(GL_TRUE); CE();
    }

    glUseProgram(0);  CE();

    // Restore the previous state
    if (!depthTestEnabled) glDisable(GL_DEPTH_TEST);
    glDepthMask(depthMaskEnabled);
    glDepthFunc(depthFunc);
    glDepthRange(depthRange[0], depthRange[1]);
    if (!cullFaceEnabled) glDisable(GL_CULL_FACE);

    CE();
}


void GLSun::calcVisibility(MapFloat* sunvis, float timestep)
{
    static int exportCounter = 0;

    QImage baseImg, canImg;
    int dx, dy, gx, gy;
    MapFloat mask; // check off wether a gridpoint is visible

    // north = (0,0,-1), west = (-1,0,0), east = (1,0,0), south = (0, 0, 1)
    sunvis->getDim(dx, dy);
    mask.setDim(dx, dy);
    mask.fill(0.0f);

    fbo->bind();
    glViewport(0, 0, fboSize.width(), fboSize.height());

    getView()->setOrthoQuadrant(0, 0);

    // first pass: terrain indices
    drawSun(getView(), 1);
    glFlush();
    baseImg = fbo->toImage();

    // second pass: // delete terrain; alpha-blended canopies
    drawSun(getView(), 2);
    glFlush();
    canImg = fbo->toImage();

    // first use exact location for incrementing
    for (int x = 0; x < baseImg.width(); x++)
    {
        for (int y = 0; y < baseImg.height(); y++)
        {
            QColor col = baseImg.pixelColor(x, y);
            colToCoord(col, gx, gy);

            if (gx < dx && gy < dy) // not the background
            {
                if (mask.get(gx, gy) < 0.5f) // not already incremented
                {
                    float r, g, b;
                    QColor viscol = canImg.pixelColor(x, y);
                    viscol.getRgbF(&r, &g, &b); // all channels store the same info so just use red

                    sunvis->set(gx, gy, sunvis->get(gx, gy) + (float)r * timestep);
                    mask.set(gx, gy, 1.0f);
                }
            }
        }
    }

    // now do a pass on the neighbours for hole filling
    for (int x = 0; x < baseImg.width(); x++)
    {
        for (int y = 0; y < baseImg.height(); y++)
        {
            QColor col = baseImg.pixelColor(x, y);
            colToCoord(col, gx, gy);

            if (gx < dx && gy < dy) // not the background
            {
                float r, g, b;
                QColor viscol = canImg.pixelColor(x, y);
                viscol.getRgbF(&r, &g, &b); // all channels store the same info so just use red

                for (int i = std::max(0, gx - 2); i <= std::min(dx - 1, gx + 2); i++)
                    for (int j = std::max(0, gy - 2); j <= std::min(dy - 1, gy + 2); j++)
                    {
                        if (mask.get(i, j) < 0.5f)
                        {
                            sunvis->set(i, j, sunvis->get(i, j) + (float)r * timestep);
                            mask.set(i, j, 1.0f);
                        }
                    }
            }
        }
    }

    fbo->release();
}

//void GLSun::updateRadianceScalingBuffers(int vwd, int vht)
//{
    // Technically, for sun calc we would always return here
    //if (shadModel != RADIANCE_SCALING)
    //    return;

    //if (_w == 0 || _h == 0 || vwd != _w || vht != _h)
    //{
    //    // std::cerr << "Delete old rad scaling buffer\n";

    //    // clean old Radiance scaling FBO data if required
    //    deleteFBOrscalingBuffers();

    //    //std::cerr << "Calling Rad scaling...\n";
    //    initRadianceScalingBuffers(vwd, vht);
    //    //std::cerr << "Done\n";
    //}
//}

void GLSun::updateHeightMap(int wd, int ht, float scx, float scy, const float* data, bool force)
{
    if (data == NULL)
    {
        std::cerr << "TRenderer::updateHeightMapTexture - NULL pointer for heightmap?";
        return;
    }

    assert(wd != 0 && ht != 0);

    if (twidth == wd && theight == ht && !force) // nothing to do - make sure binding is intact
    {
        //std::cerr << "rebind heightmap texture\n";
        if (heightmapTexture == 0)
        {
            std::cerr << "Error! Heighmap texture undefined!\n";
        }
        glActiveTexture(htmapTexUnit); CE();
        glBindTexture(GL_TEXTURE_2D, heightmapTexture); CE();
        return;
    }

    // if grid dimensions have changed:
    if (heightmapTexture != 0 && (twidth != wd || theight != ht))
    {
        // std::cerr << "- Delete texture\n";
        glDeleteTextures(1, &heightmapTexture);  CE();
        heightmapTexture = 0;
    }

    if (heightmapTexture == 0) // create texture if it does not exist
    {
        // std::cerr << "- Create texture\n";
        glGenTextures(1, &heightmapTexture); CE();
        glActiveTexture(htmapTexUnit); CE();
        glBindTexture(GL_TEXTURE_2D, heightmapTexture); CE();

        glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, wd, ht, 0, GL_RED, GL_FLOAT, (GLfloat*)data); CE();
        // no filtering
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); CE();
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); CE();
        // deal with out of array access
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
    }
    else // otherwise sub in new texture data
    {
        // std::cerr << " - sub texture\n";
        glActiveTexture(htmapTexUnit); CE();
        glBindTexture(GL_TEXTURE_2D, heightmapTexture); CE();
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, wd, ht, GL_RED, GL_FLOAT, (GLfloat*)data); CE();
    }

    // test all values for lowest terrain height:
    terrainBase = 1000000.0; // +infinity
    for (int i = 0; i < wd * ht; i++)
        terrainBase = std::min(terrainBase, data[i]);

    twidth = wd;
    theight = ht;
    scalex = scx;
    scaley = scy;

    // rebuild VAO and everything else if this was first image or new
    // dimensions or data in heightmap has changed
    deleteTerrainOpenGLbuffers();
    prepareTerrainGeometry(); // set up VBO, IBO, FBO etc
    prepareWalls(); // build capping walls for terrain
    generateNormalTexture(); // generate new normals and set up texture unit
}

void GLSun::deleteTerrainOpenGLbuffers(void)
{
    // delete old VAO etc
    if (vboTerrain != 0) glDeleteBuffers(1, &vboTerrain);  CE();
    if (iboTerrain != 0) glDeleteBuffers(1, &iboTerrain);  CE();
    if (vaoTerrain != 0) glDeleteVertexArrays(1, &vaoTerrain);  CE();
    if (normalTexture != 0) glDeleteTextures(1, &normalTexture);  CE();
    if (fboNormalMap != 0) glDeleteFramebuffers(1, &fboNormalMap);  CE();
    if (vaoScreenQuad != 0) glDeleteVertexArrays(1, &vaoScreenQuad);  CE();
    if (vboScreenQuad != 0) glDeleteBuffers(1, &vboScreenQuad);  CE();

    for (int i = 0; i < 5; i++)
    {
        if (vboWalls[i] != 0) glDeleteBuffers(1, &vboWalls[i]); CE();
        if (iboWalls[i] != 0) glDeleteBuffers(1, &iboWalls[i]); CE();
        if (vaoWalls[i] != 0) glDeleteVertexArrays(1, &vaoWalls[i]); CE();
    }
}

bool GLSun::prepareTerrainGeometry(void)
{
    GLfloat* vertexStorage;
    GLuint* indexStorage;

    vertexStorage = new GLfloat[twidth * theight * 5];
    if (vertexStorage == NULL)
    {
        std::cerr << "prepareTerrainGeometry: vertex allocation failed for " << (twidth * theight) <<
            " vertices\n";
        return false;
    }

    int vidx = 0;

    for (int y = 0; y < theight; y++)
    {
        for (int x = 0; x < twidth; x++, vidx++)
        {
            // positions: z wil be offset from height texture in the shader
            vertexStorage[5 * vidx] = (float)x / (float)(twidth - 1) * scalex;
            //vertexStorage[5*vidx+1] = 1.0f - ((float) y / (float) (height - 1)) - 0.5f;
            //vertexStorage[5*vidx+2] = 0.0f;
            vertexStorage[5 * vidx + 1] = 0.0f;
            vertexStorage[5 * vidx + 2] = /* 1.0f - */ (float)y / (float)(theight - 1) * scaley;

            // texture coordinates
            vertexStorage[5 * vidx + 3] = (float)x / (float)(twidth - 1);
            vertexStorage[5 * vidx + 4] = (float)y / (float)(theight - 1);
        }
    }


    int numStripsRequired = theight - 1;
    int numDegensRequired = 2 * (numStripsRequired - 1);
    int verticesPerStrip = 2 * twidth;

    indexSize = verticesPerStrip * numStripsRequired + numDegensRequired;

    indexStorage = new GLuint[indexSize];
    if (indexStorage == NULL)
    {
        std::cerr << "prepareTerrainGeometry: index buffer  allocation failed\n";
        return false;
    }


    int offset = 0;
    for (int y = 0; y < theight - 1; y++)
    {
        if (y > 0) // Degenerate begin: repeat first vertex
            indexStorage[offset++] = (GLuint)(y * twidth);

        for (int x = 0; x < twidth; x++)
        {
            // One part of the strip
            indexStorage[offset++] = (GLuint)((y * twidth) + x);
            indexStorage[offset++] = (GLuint)(((y + 1) * twidth) + x);
        }

        if (y < theight - 2)   // Degenerate end: repeat last vertex
            indexStorage[offset++] = (GLuint)(((y + 1) * twidth) + (twidth - 1));

    }

    // generate index array: set up for triangle strips with degeneraret tris linking strips

    glGenVertexArrays(1, &vaoTerrain); CE();
    glBindVertexArray(vaoTerrain); CE();

    // set up vertex buffer an copy in data
    glGenBuffers(1, &vboTerrain); CE();
    glBindBuffer(GL_ARRAY_BUFFER, vboTerrain); CE();
    glBufferData(GL_ARRAY_BUFFER, 5 * sizeof(GLfloat) * twidth * theight, vertexStorage, GL_STATIC_DRAW); CE();

    // enable position attribute
    glEnableVertexAttribArray(0); CE();
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), 0); CE();
    // enable texture coord attribute
    const int sz = 3 * sizeof(GLfloat);
    glEnableVertexAttribArray(1); CE();
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), reinterpret_cast<void*>(static_cast<uintptr_t>(sz))); CE();

    // set up index buffer
    glGenBuffers(1, &iboTerrain);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iboTerrain); CE();
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * indexSize, indexStorage, GL_STATIC_DRAW); CE();

    // unbind everything and clean up

    delete[] vertexStorage;
    delete[] indexStorage;

    glBindBuffer(GL_ARRAY_BUFFER, 0); CE();
    //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); CE();
    glBindVertexArray(0); CE();

    //set up screen quad for screen rendering
    glGenVertexArrays(1, &vaoScreenQuad); CE();
    glBindVertexArray(vaoScreenQuad); CE();

    glGenBuffers(1, &vboScreenQuad); CE();
    glBindBuffer(GL_ARRAY_BUFFER, vboScreenQuad); CE();
    glBufferData(GL_ARRAY_BUFFER, sizeof(screenQuad), screenQuad, GL_STATIC_DRAW); CE();

    // enable position attribute
    //glEnableVertexAttribArray(0); CE();
    //glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (void*)(0)); CE();
    // enable position attribute
    glEnableVertexAttribArray(0); CE();
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0); CE();
    // enable texture coord attribute
    const int sz2 = 2 * sizeof(GLfloat);
    glEnableVertexAttribArray(1); CE();
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), reinterpret_cast<void*>(static_cast<uintptr_t>(sz2))); CE();

    glBindBuffer(GL_ARRAY_BUFFER, 0); CE();
    glBindVertexArray(0); CE();


    // create an FBO for normal map rendering
    glGenFramebuffers(1, &fboNormalMap); CE();
    glBindFramebuffer(GL_FRAMEBUFFER, fboNormalMap); CE();

    // create texture target for normal map computation
    glActiveTexture(normalMapTexUnit); // normal map is bound to this TIU
    glGenTextures(1, &normalTexture); CE();
    glBindTexture(GL_TEXTURE_2D, normalTexture); CE();

    // set up texture state.
    // Give an empty image to OpenGL ( the last "0" )
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, twidth, theight, 0, GL_RGBA, GL_FLOAT, 0); CE();

    // no filtering
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); CE();
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); CE();
    // deal with out of array access
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);

    // configure FBO
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, normalTexture, 0); CE();

    // Set the list of draw buffers.
    GLenum DrawBuffers[1] = { GL_COLOR_ATTACHMENT0 };
    glDrawBuffers(1, DrawBuffers);  CE(); // "1" is the size of DrawBuffers

    // Always check that our framebuffer is ok
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    {
        std::cerr << "Normal map FBO initialisation failed\n";
        return false;
    }

    // unbind FBO

    glBindFramebuffer(GL_FRAMEBUFFER, 0); CE();
    return true;
}

bool GLSun::prepareWalls(void)
{
    GLfloat* vertexStorage;
    GLuint* indexStorage;

    // max space required across both vert/horiz walls - we'll re-use this array

    int vertexStorageSize = std::max(twidth, theight) * 2 * 5; // 5 float attribs per vertex
    int indexStorageSize = std::max(twidth, theight) * 2;

    vertexStorage = new GLfloat[vertexStorageSize]; // 4 walls: top, bottom, left, right
    if (vertexStorage == NULL)
    {
        std::cerr << "prepareWalls: vertex allocation failed for " << vertexStorageSize <<
            " vertices\n";
        return false;
    }

    indexStorage = new GLuint[indexStorageSize];
    if (indexStorage == NULL)
    {
        std::cerr << "prepareWalls: index buffer  allocation failed\n";
        return false;
    }

    makeXwall(theight - 1, vaoWalls[0], vboWalls[0], iboWalls[0], vertexStorage, indexStorage, false);
    normalWalls[0] = glm::vec3(0.0f, 0.0f, 1.0f); wallDrawEls[0] = 2 * twidth;
    makeXwall(0, vaoWalls[1], vboWalls[1], iboWalls[1], vertexStorage, indexStorage, true);
    normalWalls[1] = glm::vec3(0.0f, 0.0f, -1.0f); wallDrawEls[1] = 2 * twidth;
    makeYwall(0, vaoWalls[2], vboWalls[2], iboWalls[2], vertexStorage, indexStorage, false);
    normalWalls[2] = glm::vec3(-1.0f, 0.0f, 0.0f); wallDrawEls[2] = 2 * theight;
    makeYwall(twidth - 1, vaoWalls[3], vboWalls[3], iboWalls[3], vertexStorage, indexStorage, true);
    normalWalls[3] = glm::vec3(1.0f, 0.0f, 0.0f); wallDrawEls[3] = 2 * theight;
    makeBase(vaoWalls[4], vboWalls[4], iboWalls[4], vertexStorage, indexStorage);
    normalWalls[4] = glm::vec3(0.0f, -1.0f, 0.0f); wallDrawEls[4] = 4;

    delete[] vertexStorage;
    delete[] indexStorage;

    return true;
}

void GLSun::makeXwall(int atY, GLuint& vao, GLuint& vbo, GLuint& ibo, GLfloat* verts, GLuint* indices, bool reverse)
{
    int vidx = 0, x, y;

    y = atY;  // create wall at y = atY;

    for (x = 0; x < twidth; x++, vidx++)
    {
        // positions: z wil be offset from height texture in the shader
        verts[5 * vidx] = (float)x / (float)(twidth - 1) * scalex;
        //verts[5*vidx+1] = 1.0f - ((float) y / (float) (height - 1)) - 0.5f;
        //verts[5*vidx+2] = 0.0f;
        verts[5 * vidx + 1] = 0.0f;
        verts[5 * vidx + 2] = /* 1.0f - */ ((float)y / (float)(theight - 1)) * scaley;
        // texture coordinates
        verts[5 * vidx + 3] = (float)(x + 0.5f) / (float)(twidth);
        verts[5 * vidx + 4] = (float)(y + 0.5f) / (float)(theight);

        // positions: z will remain 0  - bottom edge of bottom wall, fixed to base plane
        // use negative text coords to signal this in shader
        verts[5 * vidx + 5 * twidth] = verts[5 * vidx];
        verts[5 * vidx + 1 + 5 * twidth] = verts[5 * vidx + 1];
        verts[5 * vidx + 2 + 5 * twidth] = verts[5 * vidx + 2];
        // texture coordinates
        verts[5 * vidx + 3 + 5 * twidth] = -1.0f;
        verts[5 * vidx + 4 + 5 * twidth] = -1.0f;

        indices[2 * vidx] = (reverse ? twidth - 1 - x : x);
        indices[2 * vidx + 1] = (reverse ? 2 * twidth - 1 - x : x + twidth);

    }

    glGenVertexArrays(1, &vao); CE();
    glBindVertexArray(vao); CE();

    // set up vertex buffer an copy in data
    glGenBuffers(1, &vbo); CE();
    glBindBuffer(GL_ARRAY_BUFFER, vbo); CE();
    glBufferData(GL_ARRAY_BUFFER, 5 * sizeof(GLfloat) * twidth * 2, verts, GL_STATIC_DRAW); CE();

    // enable position attribute
    glEnableVertexAttribArray(0); CE();
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), 0); CE();
    // enable texture coord attribute
    const int sz = 3 * sizeof(GLfloat);
    glEnableVertexAttribArray(1); CE();
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), reinterpret_cast<void*>(static_cast<uintptr_t>(sz))); CE();

    // set up index buffer
    glGenBuffers(1, &ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo); CE();
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * 2 * twidth, indices, GL_STATIC_DRAW); CE();

    // unbind everything and clean up

    glBindBuffer(GL_ARRAY_BUFFER, 0); CE();
    glBindVertexArray(0); CE();
}

void GLSun::makeYwall(int atX, GLuint& vao, GLuint& vbo, GLuint& ibo, GLfloat* verts, GLuint* indices, bool reverse)
{
    int vidx = 0, x, y;

    x = atX;  // create vertical wall at x = atX;

    for (y = 0; y < theight; y++, vidx++)
    {
        // positions: z  wil be offset from height texture in the shader
        verts[5 * vidx] = (float)x / (float)(twidth - 1) * scalex;
        //verts[5*vidx+1] = 1.0f - ((float) y / (float) (height - 1)) - 0.5f;
        //verts[5*vidx+2] = 0.0f;
        verts[5 * vidx + 1] = 0.0f;
        verts[5 * vidx + 2] = /*1.0f -*/ ((float)y / (float)(theight - 1)) * scaley;
        // texture coordinates
        verts[5 * vidx + 3] = (float)(x + 0.5f) / (float)(twidth);
        verts[5 * vidx + 4] = (float)(y + 0.5f) / (float)(theight);

        // positions: z will remain 0  - bottome edge of bottom wall, fixed to base plane
        // use negative text coords to signal this in shader
        verts[5 * vidx + 5 * theight] = verts[5 * vidx];
        verts[5 * vidx + 1 + 5 * theight] = verts[5 * vidx + 1];
        verts[5 * vidx + 2 + 5 * theight] = verts[5 * vidx + 2];
        // texture coordinates
        verts[5 * vidx + 3 + 5 * twidth] = -1.0f;
        verts[5 * vidx + 4 + 5 * twidth] = -1.0f;

        indices[2 * vidx] = (reverse ? theight - 1 - y : y);
        indices[2 * vidx + 1] = (reverse ? 2 * theight - 1 - y : y + theight);

    }

    glGenVertexArrays(1, &vao); CE();
    glBindVertexArray(vao); CE();

    // set up vertex buffer an copy in data
    glGenBuffers(1, &vbo); CE();
    glBindBuffer(GL_ARRAY_BUFFER, vbo); CE();
    glBufferData(GL_ARRAY_BUFFER, 5 * sizeof(GLfloat) * theight * 2, verts, GL_STATIC_DRAW); CE();

    // enable position attribute
    glEnableVertexAttribArray(0); CE();
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), 0); CE();
    // enable texture coord attribute
    const int sz = 3 * sizeof(GLfloat);
    glEnableVertexAttribArray(1); CE();
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), reinterpret_cast<void*>(static_cast<uintptr_t>(sz))); CE();

    // set up index buffer
    glGenBuffers(1, &ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo); CE();
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * 2 * theight, indices, GL_STATIC_DRAW); CE();

    // unbind everything and clean up

    glBindBuffer(GL_ARRAY_BUFFER, 0); CE();
    glBindVertexArray(0); CE();
}

void GLSun::makeBase(GLuint& vao, GLuint& vbo, GLuint& ibo, GLfloat* verts, GLuint* indices)
{
    float coords[4][2] = { {0.0f, 0.0f}, {scalex, 0.0f}, {0.0f, scaley}, {scalex, scaley} };
    for (int vidx = 0; vidx < 4; vidx++)
    {
        verts[5 * vidx] = coords[vidx][0];
        verts[5 * vidx + 1] = 0.0f; // y=0, ground plane
        verts[5 * vidx + 2] = coords[vidx][1];

        verts[5 * vidx + 3] = -1.0f; // tex coords: these verts are NOT displaced in shader...
        verts[5 * vidx + 4] = -1.0f;
    }
    // ensure correct strip winding: normal must point down into ground plane, y=0: N=(0,-1,0)
    indices[0] = 1; indices[1] = 3; indices[2] = 0; indices[3] = 2;

    // standard VAO/VBO setup
    glGenVertexArrays(1, &vao); CE();
    glBindVertexArray(vao); CE();

    // set up vertex buffer an copy in data
    glGenBuffers(1, &vbo); CE();
    glBindBuffer(GL_ARRAY_BUFFER, vbo); CE();
    glBufferData(GL_ARRAY_BUFFER, 5 * sizeof(GLfloat) * 4, verts, GL_STATIC_DRAW); CE();

    // enable position attribute
    glEnableVertexAttribArray(0); CE();
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), 0); CE();
    // enable texture coord attribute
    const int sz = 3 * sizeof(GLfloat);
    glEnableVertexAttribArray(1); CE();
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), reinterpret_cast<void*>(static_cast<uintptr_t>(sz))); CE();

    // set up index buffer
    glGenBuffers(1, &ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo); CE();
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * 4, indices, GL_STATIC_DRAW); CE();

    // unbind everything and clean up

    glBindBuffer(GL_ARRAY_BUFFER, 0); CE();
    glBindVertexArray(0); CE();
}

void GLSun::generateNormalTexture(void)
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); CE();

    glDisable(GL_DEPTH_TEST); CE();
    glDisable(GL_CULL_FACE); CE();


    glClear(GL_COLOR_BUFFER_BIT); CE();

    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport); // save current viewport

    // reset current viewport
    glViewport(0, 0, twidth, theight); CE();

    std::string shaderName = "normalShader";

    GLuint programID = (*shaders[shaderName]).programId();
    glUseProgram(programID); CE();


    GLfloat imgDims[2] = { float(twidth), float(theight) };
    GLuint locDims = glGetUniformLocation(programID, "imgSize");  CE();
    glUniform2fv(locDims, 1, imgDims); CE();

    GLuint textur = glGetUniformLocation(programID, "htMap");  CE();
    glUniform1i(textur, (GLint)(htmapTexUnit - GL_TEXTURE0)); CE(); // assumes heightmap texture is bound to this TIU

    // pass in scale
    GLfloat terDims[2] = { float(scalex), float(scaley) };
    GLuint scDims = glGetUniformLocation(programID, "scale");  CE();
    glUniform2fv(scDims, 1, terDims); CE();


    // Render to our framebuffer
    glBindFramebuffer(GL_FRAMEBUFFER, fboNormalMap); CE();

    // set shader program to normal map gen

    glBindVertexArray(vaoScreenQuad); CE();

    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);  CE();

    // unbind everthing

    glBindFramebuffer(GL_FRAMEBUFFER, 0);  CE();
    glBindVertexArray(0);  CE();
    glUseProgram(0);  CE();

    // reset viewport
    glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
}



////
// CanopyShape
///
/*
void CanopyShape::genCanopyBox(float trunkratio, float cellscale)
{

    glm::mat4 idt, tfm;
    glm::vec3 trs, rotx;
    float canopyheight;

    rotx = glm::vec3(1.0f, 0.0f, 0.0f);
    canopyheight = 1.0f - trunkratio;

    GLfloat basecol[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    canopybox->setColour(basecol);
    // canopy - tapered box
    idt = glm::mat4(1.0f);
    trs = glm::vec3(0.0f, trunkratio, 0.0f);
    tfm = glm::translate(idt, trs);
    tfm = glm::rotate(tfm, glm::radians(-90.0f), rotx);
    canopybox->genPyramid(cellscale * 1.0f, cellscale * 1.0f, canopyheight, tfm);
}

void CanopyShape::bindCanopy(Terrain * ter, View * view, MapFloat * hght, MapFloat *dnsty)
{
    int dx, dy, bndplants = 0;
    std::vector<glm::mat4> xform; // transformation to be applied to each instance
    std::vector<glm::vec4> colvar; // colour variation to be applied to each instance
    float rndoff;

    xform.clear();
    colvar.clear();

    ter->getGridDim(dx, dy);
    DiceRoller * dice = new DiceRoller(-400,400);

    for(int x = 0; x < dx; x++)
        for(int y = 0; y < dy; y++)
        {
            float h = hght->get(x, y);
            h = h * 0.3048f; // convert feet to metres
            // float d = dnsty->get(x, y);
            // stop using density because it is not available from pipeline
            // mean density is 0.91, add random variation of +- 0.4
            // rndoff = (float) dice->generate() / 10000.0f;
            float d = (0.91f + rndoff) * dnsty->get(x, y);

            // previously used #define avgtransmit

            if(h > 1.0f)
            {
                // setup transformation for individual plant, including scaling and translation
                glm::mat4 idt, tfm;
                glm::vec3 trs, sc;
                vpPoint loc = ter->toWorld(y, x, ter->getHeight(x, y)); // center of cell
                idt = glm::mat4(1.0f);
                trs = glm::vec3(loc.x, loc.y, loc.z);
                tfm = glm::translate(idt, trs); // translate to correct position
                sc = glm::vec3(1.0f, h, 1.0f); // scale to correct tree height
                tfm = glm::scale(tfm, sc);
                xform.push_back(tfm);

                colvar.push_back(glm::vec4(0.0f, 0.0f, 0.0f, d));
                // d is the blocking density of the canopy, adjusted to allow some light through beyond what can be seen from the ground
                // avgtransmit is calculated as the average transmission factor for sonoma tree species.
                bndplants++;
            }
        }

    std::cout << "binding " << bndplants << " instances..." << std::endl;
    if(!canopybox->bindInstances(nullptr, &xform, &colvar))
        cerr << "CanopyShape::bindCanopies: binding failed" << endl;
    std::cout << "finished binding instances" << std::endl;
    delete dice;
}

void CanopyShape::drawCanopy(std::vector<ShapeDrawData> &drawParams)
{
    ShapeDrawData sdd;

    sdd = canopybox->getDrawParameters();
    sdd.current = false;
    drawParams.push_back(sdd);
}
*/

////
// HemSample
///

// Efficient sampling of the hemisphere based on the Hammersley Point Set
// http://holger.dammertz.org/stuff/notes_HammersleyOnHemisphere.html

float HemSample::radicalInverse_VdC(uint bits)
{
    bits = (bits << 16u) | (bits >> 16u);
    bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
    bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
    bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
    bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
    return float(bits) * 2.3283064365386963e-10; // / 0x100000000
}

void HemSample::hammersley2d(uint i, uint N, float &u, float &v)
{
    u = (float) i / (float) N;
    v = radicalInverse_VdC(i);
}

void HemSample::convertUniform(float u, float v, Vector &dirn)
{
    float phi = v * 2.0 * PI;
    float cosTheta = 1.0 - u;
    float sinTheta = sqrt(1.0 - cosTheta * cosTheta);
    dirn = Vector(cos(phi) * sinTheta, cosTheta, sin(phi) * sinTheta);
}

void HemSample::getSample(int s, int totSamples, Vector &dirn)
{
    float u, v;
    hammersley2d((uint) s, (uint) totSamples, u, v);
    convertUniform(u, v, dirn);
}


////
// SunLight
////

SunLight::SunLight()
{
    tx = 0; ty = 0;
    month = 0; time = 0;
    latitude = 0;
    setNorthOrientation(Vector(0.0, 0.0, 1.0f)); // note: original was 0,0,-1
    refreshSun();
}

SunLight::~SunLight()
{
}

// Various conversion helper routines
float SunLight::minutesToAngle(float minutes)
{
    return ((_half_day_in_minutes - minutes) / _half_day_in_minutes) * M_PI;
}

/*
float SunLight::getAxisTiltAngle(int mnth)
{
    return -_axis_tilt + ((float) std::abs(6 - mnth) * _monthly_axis_tilt);
}*/

float SunLight::getAxisTiltAngle(float mnth)
{
    return -_axis_tilt + (std::fabs(6.0 - mnth) * _monthly_axis_tilt);
}

float SunLight::latitudeToAngle(float lat)
{
    return -lat / (180.0f / M_PI);
}

void SunLight::splitTilt(int time_of_day, float & pitch, float & roll)
{
    float f_time_of_day ((float) time_of_day);
    // Pitch
    {
        if(time_of_day <= _half_day_in_minutes) // Before midday
            pitch = 1.0f - ((f_time_of_day/_half_day_in_minutes) * 2);
        else // After midday
            pitch = -1 + (((f_time_of_day-_half_day_in_minutes)/_half_day_in_minutes) * 2);
    }

    // Roll
    {
        if(time_of_day < (_quarter_day_in_minutes))
            roll = (f_time_of_day/_quarter_day_in_minutes) * 1.0f;
        else if(f_time_of_day >= _quarter_day_in_minutes && f_time_of_day <= _3_quarters_day_in_minutes)
            roll = 1 - (((f_time_of_day-_quarter_day_in_minutes)/_half_day_in_minutes)*2.0f);
        else // 6 pm -> midnight
            roll = -1 + ((f_time_of_day - _3_quarters_day_in_minutes) / _quarter_day_in_minutes);
    }
}

void SunLight::refreshSun()
{   
    // First calculate some orientations we need values we need
    glm::vec3 east_orientation = glm::rotateY(north, (float)M_PI_2);
    glm::vec3 true_north_orientation = glm::rotate(north, glm::radians((float) (-latitude)), east_orientation); // is sign correct here?

    int sun_trajectory_radius(500000);
    float max_axis_tilt(SunLight::getAxisTiltAngle(month));
    float day_angle(SunLight::minutesToAngle((float) time));
    glm::vec3 cp_tn_and_east (glm::normalize(glm::cross(true_north_orientation, east_orientation)));

    // First calculate the sun position at midday during the equinox
    glm::vec3 sun_position ( ((float)sun_trajectory_radius) * cp_tn_and_east );

    // Now take into consideration axis tilt based on the month
    sun_position = glm::rotate( sun_position, -max_axis_tilt, east_orientation ); // PITCH

    // Now rotate around true north for the day
    sun_position = glm::rotate(sun_position, day_angle, true_north_orientation);

    // Now align to the center of the terrain (i.e the center of the terrain is at the latitude specified)
    // sun_position += glm::vec3(center.x, center.y, center.z);

    sunpos = vpPoint(sun_position[0], sun_position[1], sun_position[2]);
}

//void SunLight::genSphereSun()
//{
//    glm::mat4 idt;
//
//    // simple unit diameter sphere
//    idt = glm::mat4(1.0f);
//    GLfloat suncol[] = {1.0f, 1.0f, 0.0f, 1.0f};
//    //sunRender.genSphere(20.0f, 20, 20, idt);
//    //sunRender.setColour(suncol);
//}

//void SunLight::bindDiffuseSun(View * view, Terrain * ter)
//{
//    std::vector<glm::mat4> xform; // transformation to be applied to each instance
//    std::vector<glm::vec4> colvar; // colour variation to be applied to each instance
//    HemSample hem;
//
//    xform.clear();
//    colvar.clear();
//
//    for(int t = 0; t < 100; t++)
//    {
//        // setup transformation for individual sun
//        glm::mat4 idt, tfm;
//        glm::vec3 trs, sc;
//        Vector sundir;
//        vpPoint termid, loc;
//
//        hem.getSample(t, 100, sundir);
//        sundir.mult(400.0f);
//        ter->getMidPoint(termid);
//        loc = vpPoint(termid.x+sundir.i, termid.y+sundir.j, termid.z+sundir.k);
//        //cerr << "sundir = " << sundir.i << " " << sundir.j << " " << sundir.k << endl;
//
//        idt = glm::mat4(1.0f);
//        trs = glm::vec3(loc.x+center.x, loc.y+center.y, loc.z+center.z);
//        tfm = glm::translate(idt, trs);
//        xform.push_back(tfm);
//        glm::vec4 col = glm::vec4(0.0f, 0.0f, 0.0f, 0.0f); // no colour variation
//        colvar.push_back(col); // colour variation
//     }
//    //if(!sunRender.bindInstances(view, &xform, &colvar))
//    //    cerr << "SUN BINDING FAILED" << endl;
//}
//
//void SunLight::bindSun(View * view)
//{
//    std::vector<glm::mat4> xform; // transformation to be applied to each instance
//    std::vector<glm::vec4> colvar; // colour variation to be applied to each instance
//
//    xform.clear();
//    colvar.clear();
//    float tothours = 0.0f;
//
//    for(int t = 0; t < 1440; t+= 1)
//    {
//        // setup transformation for individual sun
//        glm::mat4 idt, tfm;
//        glm::vec3 trs, sc;
//
//        setTime(t);
//        Vector sundir = Vector(sunpos.x, sunpos.y, sunpos.z);
//        sundir.normalize();
//        sundir.mult(400.0f);
//        vpPoint loc = vpPoint(sundir.i, sundir.j, sundir.k);
//        //cerr << "sundir = " << sundir.i << " " << sundir.j << " " << sundir.k << endl;
//        if(loc.y > 0.0f)
//            tothours += (1.0 / 60.0f);
//
//        idt = glm::mat4(1.0f);
//        trs = glm::vec3(loc.x+center.x, loc.y+center.y, loc.z+center.z);
//        tfm = glm::translate(idt, trs);
//        xform.push_back(tfm);
//        glm::vec4 col = glm::vec4(abs((float) t - 720.0f) / -1400.0f, abs((float) t - 720.0f) / -1400.0f, abs((float) t - 720.0f) / -1400.0f, 0.0f); // midday is peak colour, all others are darker
//        colvar.push_back(col); // colour variation
//     }
//    //if(!sunRender.bindInstances(view, &xform, &colvar))
//    //    cerr << "SUN BINDING FAILED" << endl;
//    cerr << "TOTAL HOURS OF SUNLIGHT = " << tothours << endl;
//}

//void SunLight::drawSun(std::vector<ShapeDrawData>& drawParams)
//{
//    ShapeDrawData sdd;
//
//    sdd = sunRender.getDrawParameters();
//    sdd.current = false;
//    drawParams.push_back(sdd);
//}

void SunLight::projectSun(Terrain * ter, std::vector<MapFloat> &sunmaps, GLSun * glsun, std::vector<float> &sunhours, int minutestep, int startm, int endm, int mincr, bool enable)
{
    float tothours;
    int dx, dy;
    vpPoint gridpos, xsectpos;
    Vector sunvec;
    float hourstep = (float) minutestep / 60.0f;
    Timer tmm;
    float tstep, tx, ty;

    ter->getGridDim(dx, dy);
    ter->getTerrainDim(tx, ty);
    tstep = 1.0f * tx / (float) dx;
    sunhours.clear();

    int month_incr = mincr;
    int start_month = startm;
    int end_month = endm;

    for (int m = start_month; m <= end_month; m += month_incr) // months of the year
    {
        tothours = 0.0f;

        sunmaps[m-1].fill(0.0f);

        // we choose one day of the month
        int d = 21;
        setMonth((float) m-1 + ((float) (d-1) / 30.0f));

        tmm.start();
        for (int t = 0; t < 1440; t += minutestep) // 24 hours
        {
            setTime(t);

            // cerr << "zoom dist = " << glsun->getView()->getZoom() << endl;
            // cerr << "sunpos " << sunpos.x << ", " << sunpos.y << ", " << sunpos.z << endl;
            if(sunpos.y > 0.0f) // above the horizon
            {
                if (enable)
                {
                    sunvec = Vector(sunpos.x, sunpos.y, -sunpos.z);
                    sunvec.normalize();
                    glsun->getView()->sundir(sunvec);
                    glsun->calcVisibility(&sunmaps[m-1], hourstep);
                }
                tothours += hourstep; // sun above horizon
            }
        }
        tmm.stop();
        cerr << "Month " << m << " Sunlight Pass Complete in " << tmm.peek() << " seconds" << endl;
        cerr << "with " << tothours << " of sunlight" << endl;
        sunhours.push_back(tothours);
    }
}

void SunLight::diffuseSun(Terrain* ter, MapFloat* diffusemap, GLSun* glsun, int numSamples, bool enable)
{
    int dx, dy;
    Timer tmm;
    HemSample hem;

    ter->getGridDim(dx, dy);
    diffusemap->setDim(dx, dy);
    diffusemap->fill(0.0f);
    tmm.start();
    if (enable)
    {
        for (int s = 0; s < numSamples; s++) // hemisphere samples
        {
            Vector sunvec;
            hem.getSample(s, numSamples, sunvec);
            if (sunvec.j > 0.0f) // above the horizon
            {
                glsun->getView()->sundir(sunvec);
                glsun->calcVisibility(diffusemap, 1.0f);
            }
        }

        // normalization
        for (int x = 0; x < dx; x++)
            for (int y = 0; y < dy; y++)
            {
                diffusemap->set(x, y, diffusemap->get(x, y) / (float)numSamples);
            }
    }
    tmm.stop();
    cerr << "Diffuse sampling complete in " << tmm.peek() << " seconds" << endl;
}

void SunLight::mergeSun(std::vector<MapFloat> &sunmaps, MapFloat * diffusemap, std::vector<float> cloudiness, std::vector<float> sunhours)
{
    int dx, dy;
    float direct, diffuse;

    diffusemap->getDim(dx, dy);
    for(int m = 0; m < 12; m++) // months of the year
    {
        blurSun(sunmaps[m], 2, 2);
        blurSun((* diffusemap), 2, 2);
        for(int x = 0; x < dx; x++)
            for(int y = 0; y < dy; y++)
            {
                // direct sunlight portion
                direct = sunmaps[m].get(x, y) * (1.0f - cloudiness[m]);
                // diffuse sunlight portion
                diffuse = diffusemap->get(x, y) * cloudiness[m] * sunhours[m];

                // includes re-orientation
                //sunmaps[m].set(y, x, direct+diffuse);
                sunmaps[m].set(x, y, direct+diffuse);
            }

        // cerr << "Cloudiness for month " << m << " is " << cloudiness[m] << " and sun hours are " << sunhours[m] << endl;
    }
}

//void SunLight::applyAlpha(GLSun * sun, std::vector<MapFloat> &sunmaps)
//{
//    int dx, dy;
//
//    // apply previously derived alpha map
//    sun->getAlpha()->getDim(dx, dy);
//
//    // do radial smoothing. Note this changes the alpha map
//    radialBlur((* sun->getAlpha()), 15);
//
//    for(int m = 0; m < 12; m++) // months of the year
//    {
//        for(int x = 0; x < dx; x++)
//            for(int y = 0; y < dy; y++)
//                sunmaps[m].set(x, y, sunmaps[m].get(x, y) * (1.0f - sun->getAlpha()->get(x, y)));
// //               if(sun->getAlpha()->get(x, y) > 0.1f)
// //                   sunmaps[m].set(x, y, 0.0f);
//    }
//}
//
//void SunLight::radialBlur(MapFloat &map, int radius)
//{
//    MapFloat newmap;
//    int filterwidth = radius+1;
//    int sqrrad = radius * radius;
//
//    int dx, dy;
//    map.getDim(dx, dy);
//    newmap.setDim(dx, dy);
//
//    for(int x = 0; x < dx; x++)
//        for(int y = 0; y < dy; y++)
//        {
//            float avg = 0.0f;
//            int cnt = 0;
//
//            for(int cx = x-filterwidth; cx <= x+filterwidth; cx++)
//                for(int cy = y-filterwidth; cy <= y+filterwidth; cy++)
//                {
//                    int r = (cx-x) * (cx-x) + (cy-y) * (cy-y);
//                    // within radius and within bounds
//                    if(r <= sqrrad && cx >= 0 && cx < dx && cy >= 0 && cy < dy)
//                    {
//                        avg += map.get(cx, cy);
//                        cnt++;
//                    }
//                }
//            newmap.set(x,y, avg / (float) cnt);
//        }
//
//    for(int x = 0; x < dx; x++)
//        for(int y = 0; y < dy; y++)
//            map.set(x, y, newmap.get(x, y));
//}

void SunLight::blurSun(MapFloat &map, int filterwidth, int passes)
{
    float filterarea;
    MapFloat newmap;

    filterarea = (float) ((filterwidth*2+1)*(filterwidth*2+1));

    int dx, dy;
    map.getDim(dx, dy);
    newmap.setDim(dx, dy);

    for(int i = 0; i < passes; i++)
    {
        for(int x = 0; x < dx; x++)
            for(int y = 0; y < dy; y++)
            {
                float avg = 0.0f;

                for(int cx = x-filterwidth; cx <= x+filterwidth; cx++)
                    for(int cy = y-filterwidth; cy <= y+filterwidth; cy++)
                    {
                            if(cx < 0 || cx >= dx || cy < 0 || cy >= dy)
                                avg += map.get(x, y);
                            else
                                avg += map.get(cx, cy);
                    }
                    newmap.set(x,y, avg / filterarea);
            }

        for(int x = 0; x < dx; x++)
            for(int y = 0; y < dy; y++)
                map.set(x, y, newmap.get(x, y));
    }
}
