#ifndef _CORE_H_
#define _CORE_H_

#include <cmath>
#include <vector>
#include <random>
#include <QtGui/QColor>

class Vector3 {
protected:
    double c[3];

public:
    Vector3() {
        c[0] = c[1] = c[2] = 0;
    }
    explicit Vector3(double d) {
        c[0] = c[1] = c[2] = d;
    }
    explicit Vector3(double d0, double d1, double d2) {
        c[0] = d0;
        c[1] = d1;
        c[2] = d2;
    }

    double& operator[] (int i) { return c[i]; };
    double operator[] (int i) const { return c[i]; };
    Vector3 operator-() const { return Vector3(-c[0], -c[1], -c[2]); }

    friend Vector3 operator+(const Vector3& u, const Vector3& v) { return Vector3(u[0]+v[0], u[1]+v[1], u[2]+v[2]); };
    friend Vector3 operator-(const Vector3& u, const Vector3& v) { return Vector3(u[0]-v[0], u[1]-v[1], u[2]-v[2]); };
    friend Vector3 operator*(const Vector3& u, double a) { return Vector3(u[0]*a, u[1]*a, u[2]*a); }
    friend Vector3 operator*(double a, const Vector3& v) { return v * a; }
    friend Vector3 operator/(const Vector3& u, double a) { return Vector3(u[0]/a, u[1]/a, u[2]/a); }
    friend double operator* (const Vector3& u, const Vector3& v)  { return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]; }

    friend double Norm(const Vector3& u) { return sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]); }
    friend double SquaredNorm(const Vector3& u) { return u[0]*u[0] + u[1]*u[1] + u[2]*u[2]; }
    friend Vector3 Normalized(const Vector3& u) { return u/Norm(u); }

    friend double dot(const Vector3& u, const Vector3& v) { return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]; }
    friend Vector3 cross(const Vector3& u, const Vector3& v) {
        return Vector3(u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]);
    }

    friend Vector3 Bilinear(const Vector3& a00, const Vector3& a10, const Vector3& a11, const Vector3& a01, const double& u, const double& v) {
        return (1 - u) * (1 - v) * a00 + (1 - u) * (v)*a01 + (u) * (1 - v) * a10 + (u) * (v)*a11;
    }
};

class Vector2 {
protected:
    double c[2];

public:
    Vector2() {
        c[0] = c[1] = 0;
    }
    explicit Vector2(double d) {
        c[0] = c[1] = d;
    }
    explicit Vector2(double d0, double d1) {
        c[0] = d0;
        c[1] = d1;
    }
    Vector2(const Vector3& v) {
        c[0] = v[0];
        c[1] = v[1];
    }

    double& operator[] (int i) { return c[i]; };
    double operator[] (int i) const { return c[i]; };

    Vector2 operator- () const { return Vector2(-c[0], -c[1]); };
    friend Vector2 operator+(const Vector2& u, const Vector2& v) { return Vector2(u[0]+v[0], u[1]+v[1]); };
    friend Vector2 operator-(const Vector2& u, const Vector2& v) { return Vector2(u[0]-v[0], u[1]-v[1]); };
    friend Vector2 operator*(const Vector2& u, double a) { return Vector2(u[0]*a, u[1]*a); }
    friend Vector2 operator*(double a, const Vector2& v) { return v * a; }
    friend Vector2 operator/(const Vector2& u, double a) { return Vector2(u[0]/a, u[1]/a); }
    friend double operator* (const Vector2& u, const Vector2& v)  { return u[0]*v[0] + u[1]*v[1]; }

    friend double Norm(const Vector2& u) { return sqrt(u[0]*u[0] + u[1]*u[1]); }
    friend double SquaredNorm(const Vector2& u) { return u[0]*u[0] + u[1]*u[1]; }
    friend Vector2 Normalized(const Vector2& u) { return u/Norm(u); }
    double Angle() const { return std::atan2(c[1], c[0]); }

    Vector3 toVector3(double d) const { return Vector3(c[0], c[1], d); }
};


class Ray
{
protected:
    Vector3 p; // Origin of the ray.
    Vector3 d; // Direction.

public:
    Ray() {}
    explicit Ray(const Vector3& p, const Vector3& d) : p(p), d(d) {}

    Vector3 origin() const { return p; }
    Vector3 direction() const { return d; }

    Vector3 operator()(double t) const { return p + t * d; }

    Ray reflect(const Vector3& p, const Vector3& n) { return Ray(p, n - 2 * n * dot(d, n)); }
};


class Box2 {
protected:
    Vector2 bmin; // min
    Vector2 bmax; // max

public:
    Box2() : bmin(0), bmax(0) {};
    Box2(const double a, const double b) : bmin(-Vector2(a/2, b/2)), bmax(Vector2(a/2,b/2)) {}
    Box2(const Vector2& pmin, const Vector2& pmax) : bmin(pmin), bmax(pmax) {}
    Box2(const Vector2& c, double r) : bmin(c - Vector2(r)), bmax(c + Vector2(r)) {}

    Vector2 getMin() const { return bmin; }
    Vector2 getMax() const { return bmax; }
    Vector2 center() const { return 0.5*(bmin + bmax); }
    double radius() const { return 0.5 * Norm(bmax - bmin); }
    double width() const { return bmax[0] - bmin[0]; }
    double height() const { return bmax[1] - bmin[1]; }

    bool Inside(const Vector2& p) const {
      if ((p[0] < bmin[0]) || (p[0] > bmax[0]) || (p[1] < bmin[1]) || (p[1] > bmax[1]))
        return false;
      return true;
    }
    bool Intersect(const Vector2& s0, const Vector2& s1, double& tmin, double& tmax);
};


class Box3 {
protected:
    Vector3 bmin; // min
    Vector3 bmax; // max

public:
    Box3() : bmin(0), bmax(0) {};
    Box3(const Vector3& pmin, const Vector3& pmax) : bmin(pmin), bmax(pmax) {}
    Box3(const Box2& box2, double zmin, double zmax) {
        bmin = Vector3(box2.getMin()[0], box2.getMin()[1], zmin);
        bmax = Vector3(box2.getMax()[0], box2.getMax()[1], zmax);
    }

    Vector3 getMin() const { return bmin; }
    Vector3 getMax() const { return bmax; }
    Vector3 center() const { return 0.5*(bmin + bmax); }
    Vector3 size() const { return bmax - bmin; }
    double radius() const { return 0.5 * Norm(bmax - bmin); }
    double width() const { return bmax[0] - bmin[0]; }
    double height() const { return bmax[1] - bmin[1]; }
    double depth() const { return bmax[2] - bmin[2]; }

    int Intersect(const Ray& ray, double& tmin, double& tmax) const;
};


class Camera {
protected:
    Vector3 eye;      // Eye
    Vector3 at;       // Look at point
    Vector3 up;       // Up vector
    double cah;       // Camera aperture horizontal
    double cav;       // Camera aperture vertical
    double nearplane; // Near plane
    double farplane;  // Far plane
    double fl;        // Focal length

public:
    Camera();
    Camera(const Vector3& eye, const Vector3& at, const Vector3& up = Vector3(0,0,1), double near = 1.0, double far = 100000.0);

    Vector3 getEye() const { return eye; }
    Vector3 getAt() const { return at; }
    Vector3 getUp() const { return up; }
    Vector3 getViewDir() const { return Normalized(at - eye); }
    double getNearPlane() const { return nearplane; }
    double getFarPlane() const { return farplane; }
    double getAngleOfViewH(double, double) const;
    double getAngleOfViewV(double, double) const;

    void setAt(const Vector3& p) { at = p; up = Vector3(0,0,1); }
    void setEye(const Vector3& p) { eye = p; }
    void setPlanes(double n, double f) { nearplane = n; farplane = f;}

    // Move camera around eye
    void upDownRound(double a);
    void leftRightRound(double a);
    void backForth(double a, bool moveAt = false);

    // Move camera in a plane
    void upDownPlane(double);
    void leftRightPlane(double);

    Ray pixelToRay(int px, int py, int w, int h) const;

    static Camera View(const Box3& box);
};


inline Vector3 fromQColor(const QColor& c) {
    return Vector3(c.red() / 255.0, c.green() / 255.0, c.blue() / 255.0);
}

inline QColor toQColor(const Vector3& c) {
    return QColor(int(255.0 * std::max(std::min(c[0], 1.0), 0.0)),
                  int(255.0 * std::max(std::min(c[1], 1.0), 0.0)),
                  int(255.0 * std::max(std::min(c[2], 1.0), 0.0)),
                  255);
}

class ColorPalette
{
protected:
    std::vector<Vector3> colors;
    std::vector<double> anchors;

public:
    ColorPalette() : colors({Vector3(1)}), anchors({0}) {}
    ColorPalette(const std::vector<Vector3>& c, const std::vector<double>& a) : colors(c), anchors(a) {}

    Vector3 getColor(double u) const;

    static ColorPalette CoolWarm() {
        const Vector3 Cool = Vector3(97, 130, 234) / 255.0;
        const Vector3 White = Vector3(221, 221, 221) / 255.0;
        const Vector3 Warm = Vector3(220, 94, 75) / 255.0;
        return ColorPalette({Cool, White, Warm}, {0, 0.5, 1});
    }

    static ColorPalette GreenGreyBrown() {
        const Vector3 Brown = Vector3(74, 41, 13)/255.0;
        const Vector3 Green = Vector3(20, 156, 38)/255.0;
        const Vector3 Grey = Vector3(232, 232, 232)/255.0;
        return ColorPalette({ Brown, Grey, Green }, { 0, 0.5, 1 });
    }

    static ColorPalette Relief() {
        const std::vector<Vector3> c = { Vector3(0.627, 0.863, 0.411),
                                         Vector3(1.00, 0.90, 0.45),
                                         Vector3(0.659, 0.607, 0.541),
                                         Vector3(0.95, 0.95, 0.95) };
        return ColorPalette(c, {0, 0.375, 0.625, 1.0});
    }

    static ColorPalette Blues() {
        const std::vector<Vector3> c = {
            Vector3(222,235,247)/255.0,
            Vector3(198,219,239)/255.0,
            Vector3(158,202,225)/255.0,
            Vector3(107,174,214)/255.0,
            Vector3(66,146,198)/255.0,
            Vector3(33,113,181)/255.0,
            Vector3(8,81,156)/255.0
        };
        return ColorPalette(c, { 0, 0.167, 0.33, 0.5, 0.667, 0.833, 1 });
    }

};


class LookupPalette
{
protected:
    std::vector<Vector3> c; //!< Set of colors.
public:
    LookupPalette(const std::vector<Vector3>& v) : c(v) {};
    virtual Vector3 getColor(int) const;
};

#endif
