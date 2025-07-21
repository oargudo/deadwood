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


#ifndef GLHEADERS_H_
#define GLHEADERS_H_

#if defined(_WIN32)
#include <glew.h>
#elif defined(__APPLE__)
#include <GL/glew.h>
#include <OpenGL/gl3.h>
#include <OpenGL/glu.h>
#else
#define GL3_PROTOTYPES
#include <GL/glew.h>
#include <GL/gl3.h>
#include <GL/glu.h>
#endif

#include <QtCore/QDebug>
#include <QtCore/QString>

/*
 * Append the CE() macro on the same line after an OpenGL function call to check for errors.
 * Use IF_FAIL("Error message") = ErrorCodeReturningFunction(); To print out a message on failure.
 */

class CheckSucc{
 public:
    CheckSucc(QString msg){
        succ = true;
        msg_ = msg;
    }

    CheckSucc operator = (bool rhs){
        succ = rhs;
        if(!rhs)
            qDebug() << msg_;
        return *this;
    }

    CheckSucc operator = (int rhs){
        if(rhs == -1){
            succ = false;
            qDebug() << msg_;
        }
        return *this;
    }

    bool operator == (bool rhs){
        return this->succ == rhs;
    }


 private:
    bool succ;
    QString msg_;
};

#ifndef _DEBUG_GL_FUNCS_
#define _DEBUG_GL_FUNCS_

#define IF_FAIL(msg) CheckSucc(QString(QString(msg) + QString(" @ Line ") + QString::number(__LINE__) +  QString(" of ") +  QString(__FILE__)))

#ifdef CC_GL_DEBUG_
    static GLenum gl_err = 0;
    #define CE() gl_err = glGetError();\
    if (gl_err != GL_NO_ERROR) {\
        const char* err_str = reinterpret_cast<const char *>(gluErrorString(gl_err));\
        QString errString(err_str);\
        qDebug() << "GL Error: " << gl_err << " (" << errString << ") on line" << __LINE__ << "of" << __FILE__;\
    }
#else
    #define CE()
#endif  // CC_GL_DEBUG_

#endif

#endif  // GLHEADERS_H_
