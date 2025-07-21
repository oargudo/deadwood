#include "core.h"


bool Box2::Intersect(const Vector2 &s0, const Vector2 &s1, double &tmin, double &tmax)
{
    const double epsilon = 1.0e-5;

    tmin = -1e16;
    tmax = 1e16;

    const Vector2& a = bmin;
    const Vector2& b = bmax;
    Vector2 p = s0;
    Vector2 d = s1 - s0;

    double t;
    // Ox
    if (d[0] < -epsilon) {
        t = (a[0] - p[0]) / d[0];
        if (t < tmin)
            return false;
        if (t <= tmax)
            tmax = t;
        t = (b[0] - p[0]) / d[0];
        if (t >= tmin) {
            if (t > tmax)
                return false;
            tmin = t;
        }
    }
    else if (d[0] > epsilon) {
        t = (b[0] - p[0]) / d[0];
        if (t < tmin)
            return false;
        if (t <= tmax)
            tmax = t;
        t = (a[0] - p[0]) / d[0];
        if (t >= tmin) {
            if (t > tmax)
                return false;
            tmin = t;
        }
    }
    else if (p[0]<a[0] || p[0]>b[0])
        return false;

    // Oy
    if (d[1] < -epsilon) {
        t = (a[1] - p[1]) / d[1];
        if (t < tmin)
            return false;
        if (t <= tmax)
            tmax = t;
        t = (b[1] - p[1]) / d[1];
        if (t >= tmin) {
            if (t > tmax)
                return false;
            tmin = t;
        }
    }
    else if (d[1] > epsilon) {
        t = (b[1] - p[1]) / d[1];
        if (t < tmin)
            return false;
        if (t <= tmax)
            tmax = t;
        t = (a[1] - p[1]) / d[1];
        if (t >= tmin) {
            if (t > tmax)
                return false;
            tmin = t;
        }
    }
    else if (p[1]<a[1] || p[1]>b[1])
        return false;

    return true;
}


int Box3::Intersect(const Ray& ray, double& tmin, double& tmax) const
{
  const double epsilon = 1.0e-5;

  tmin = -1e16;
  tmax = 1e16;

  Vector3 p = ray.origin();
  Vector3 d = ray.direction();

  double t;
  // Ox
  if (d[0] < -epsilon)
  {
    t = (bmin[0] - p[0]) / d[0];
    if (t < tmin)
      return 0;
    if (t <= tmax)
      tmax = t;
    t = (bmax[0] - p[0]) / d[0];
    if (t >= tmin)
    {
      if (t > tmax)
        return 0;
      tmin = t;
    }
  }
  else if (d[0] > epsilon)
  {
    t = (bmax[0] - p[0]) / d[0];
    if (t < tmin)
      return 0;
    if (t <= tmax)
      tmax = t;
    t = (bmin[0] - p[0]) / d[0];
    if (t >= tmin)
    {
      if (t > tmax)
        return 0;
      tmin = t;
    }
  }
  else if (p[0]<bmin[0] || p[0]>bmax[0])
    return 0;

  // Oy
  if (d[1] < -epsilon)
  {
    t = (bmin[1] - p[1]) / d[1];
    if (t < tmin)
      return 0;
    if (t <= tmax)
      tmax = t;
    t = (bmax[1] - p[1]) / d[1];
    if (t >= tmin)
    {
      if (t > tmax)
        return 0;
      tmin = t;
    }
  }
  else if (d[1] > epsilon)
  {
    t = (bmax[1] - p[1]) / d[1];
    if (t < tmin)
      return 0;
    if (t <= tmax)
      tmax = t;
    t = (bmin[1] - p[1]) / d[1];
    if (t >= tmin)
    {
      if (t > tmax)
        return 0;
      tmin = t;
    }
  }
  else if (p[1]<bmin[1] || p[1]>bmax[1])
    return 0;

  // Oz
  if (d[2] < -epsilon)
  {
    t = (bmin[2] - p[2]) / d[2];
    if (t < tmin)
      return 0;
    if (t <= tmax)
      tmax = t;
    t = (bmax[2] - p[2]) / d[2];
    if (t >= tmin)
    {
      if (t > tmax)
        return 0;
      tmin = t;
    }
  }
  else if (d[2] > epsilon)
  {
    t = (bmax[2] - p[2]) / d[2];
    if (t < tmin)
      return 0;
    if (t <= tmax)
      tmax = t;
    t = (bmin[2] - p[2]) / d[2];
    if (t >= tmin)
    {
      if (t > tmax)
        return 0;
      tmin = t;
    }
  }
  else if (p[2]<bmin[2] || p[2]>bmax[2])
    return 0;

  return 1;
}



Camera::Camera()
{
    Camera::eye = Vector3(0.0);
    Camera::at = Vector3(0.0, 1.0, 0.0);
    Camera::up = Vector3(0.0, 0.0, 1.0);

    // Near and far planes
    Camera::nearplane = 1.0;
    Camera::farplane = 1000.0;

    // Aperture
    Camera::cah = 0.980;
    Camera::cav = 0.735;
    Camera::fl = 35.0;
}

Camera::Camera(const Vector3& eye, const Vector3& at, const Vector3& up, double near, double far)
{
    Camera::eye = eye;
    Camera::at = at;
    Camera::up = up;

    // Near and far planes
    Camera::nearplane = near;
    Camera::farplane = far;

    // Aperture
    Camera::cah = 0.980;
    Camera::cav = 0.735;
    Camera::fl = 35.0;
}

double Camera::getAngleOfViewH(double, double) const
{
    return 2.0 * atan(cah * 25.4 * 0.5 / fl);
}

double Camera::getAngleOfViewV(double w, double h) const
{
    double avh = getAngleOfViewH(w, h);
    return 2.0 * atan(tan(avh / 2.0) * double(h) / double(w));
}

void Camera::upDownRound(double a)
{
    Vector3 z = at - eye;
    double length = Norm(z);
    z = z/length;
    Vector3 left = Normalized(cross(up, z));

    // Rotate
    z = z * cos(a) + up * sin(a);

    // Update Vector
    up = cross(z, left);
    eye = at - z * length;
}

void Camera::leftRightRound(double a)
{
    Vector3 e = eye - at;
    Vector3 left = cross(up, e);
    e = Vector3(e[0] * cos(a) - e[1] * sin(a), e[0] * sin(a) + e[1] * cos(a), e[2]);
    left = Vector3(left[0] * cos(a) - left[1] * sin(a), left[0] * sin(a) + left[1] * cos(a), 0.0);
    up = Normalized(cross(left, -e));
    eye = at + e;
}

void Camera::backForth(double a, bool moveAt)
{
    Vector3 z = at - eye;
    double length = Norm(z);
    z = z/length;
    eye = eye + a * z;
    if (moveAt) {
        at = at + a * z;
    }
}

void Camera::upDownPlane(double a)
{
    Vector3 z = at - eye;
    double length = Norm(z);
    z = z/length;
    Vector3 left = Normalized(cross(Vector3(0, 0, 1), z));

    eye = eye + a * cross(z, left);
    at = at + a * cross(z, left);
}

void Camera::leftRightPlane(double a)
{
    Vector3 z = at - eye;
    z[2] = 0.0;
    double length = Norm(z);
    z = z/length;
    Vector3 left = Normalized(cross(Vector3(0, 0, 1), z));

    eye = eye + a * left;
    at = at + a * left;
}

Ray Camera::pixelToRay(int px, int py, int w, int h) const
{
    // Get coordinates
    Vector3 view = getViewDir();
    Vector3 horizontal = Normalized(cross(view, up));
    Vector3 vertical = Normalized(cross(horizontal, view));

    double length = 1.0;

    // Convert to radians
    double rad = getAngleOfViewV(w, h);  // fov

    double vLength = tan(rad / 2.0) * length;
    double hLength = vLength * (double(w) / double(h));
    vertical = vertical*vLength;
    horizontal = horizontal*hLength;

    // Translate mouse coordinates so that the origin lies in the center of the view port
    double x = px - w / 2.0;
    double y = h / 2.0 - py;

    // Scale mouse coordinates so that half the view port width and height becomes 1.0
    x /= w / 2.0;
    y /= h / 2.0;

    // Direction is a linear combination to compute intersection of picking ray with view port plane
    return Ray(eye, Normalized(view * length + horizontal * x + vertical * y));
}

Camera Camera::View(const Box3& box)
{
    Vector3 v = 0.5*(box.getMax() - box.getMin());
	v[2] = 0;
	double r = Norm(v);
	v = 2*v;
	v[2] = -r;
    return Camera(box.center() - v, box.center(), Vector3(0.0, 0.0, 1.0), r, 3*r);
}


Vector3 ColorPalette::getColor(double t) const {
    if (colors.size() == 0) return Vector3(1);
    if (colors.size() == 1) return colors[0];

    if (t <= anchors.front()) return colors.front();
    if (t >= anchors.back()) return colors.back();
    for (int i = 0; i < int(colors.size() - 1); i++) {
        if (t < anchors[i+1]) {
            double s = (t - anchors[i]) / (anchors[i + 1] - anchors[i]);
            return (1-s)*colors[i] + s*colors[i+1];
        }
    }

    return Vector3(0);
}

Vector3 LookupPalette::getColor(int i) const {
    if ((i < 0) || (i >= int(c.size()))) return Vector3(0,0,0);
    return c.at(i);
}
