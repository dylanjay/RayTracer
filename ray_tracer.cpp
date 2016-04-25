/**
 * ray_tracer.cpp
 * CS230
 * -------------------------------
 * Implement ray tracer here.
 */

#define SET_RED(P, C)   (P = (((P) & 0x00ffffff) | ((C) << 24)))
#define SET_GREEN(P, C)  (P = (((P) & 0xff00ffff) | ((C) << 16)))
#define SET_BLUE(P, C) (P = (((P) & 0xffff00ff) | ((C) << 8)))

#include "ray_tracer.h"

using namespace std;

const double Object::small_t=1e-6;
//--------------------------------------------------------------------------------
// utility functions
//--------------------------------------------------------------------------------
double sqr(const double x)
{
    return x*x;
}

Pixel Pixel_Color(const Vector_3D<double>& color)
{
    Pixel pixel=0;
    SET_RED(pixel,(unsigned char)(min(color.x,1.0)*255));
    SET_GREEN(pixel,(unsigned char)(min(color.y,1.0)*255));
    SET_BLUE(pixel,(unsigned char)(min(color.z,1.0)*255));
    return pixel;
}

double dist(const Vector_3D<double> obj0, const Vector_3D<double> obj1)
{
  double dx = obj0.x - obj1.x;
  double dy = obj0.y - obj1.y;
  double dz = obj0.z - obj1.z;
  return sqr(dx) + sqr(dy) + sqr(dz);
}
//--------------------------------------------------------------------------------
// Shader
//--------------------------------------------------------------------------------
Vector_3D<double> Phong_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,const Vector_3D<double>& intersection_point,const Vector_3D<double>& same_side_normal) const
{
    Vector_3D<double> color;

    // TODO: determine the color   
    for(int i = 0; i < world.lights.size(); i++)
    {
      Vector_3D<double> k_a = color_ambient;
      Vector_3D<double> col = world.lights[i] -> Emitted_Light(ray);
      Vector_3D<double> L_m = (world.lights[i] -> position - intersection_point).Normalized();
      Vector_3D<double> N = same_side_normal;
      Vector_3D<double> k_d = color_diffuse;
      Vector_3D<double> k_s = color_specular;
      Vector_3D<double> R_m = (N * 2 * Vector_3D<double>::Dot_Product(L_m,N)) - L_m;
      Vector_3D<double> V = (world.camera.position - intersection_point).Normalized();
      Ray shadow(intersection_point, world.lights[i] -> position);
      color += k_a*col*1.0;
      if(world.enable_shadows)
      {
	for(int j = 0; j < world.objects.size(); j++)
	{
	  if(!world.objects[j] -> Intersection(shadow))
	  {
	    color += k_d * Vector_3D<double>::Dot_Product(L_m,N) * col + k_s * pow(Vector_3D<double>::Dot_Product(R_m,V),specular_power) * col; 
	  }
	}
      }
      else
      {
	color += k_d * Vector_3D<double>::Dot_Product(L_m,N) * col + k_s * pow(Vector_3D<double>::Dot_Product(R_m,V),specular_power) * col;
      }
    }
    

    return color;
}

Vector_3D<double> Reflective_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,const Vector_3D<double>& intersection_point,const Vector_3D<double>& same_side_normal) const
{
    Vector_3D<double> color;

    // TODO: determine the color

    return color;
}

Vector_3D<double> Flat_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,const Vector_3D<double>& intersection_point,const Vector_3D<double>& same_side_normal) const
{  
  return color;
}

//--------------------------------------------------------------------------------
// Objects
//--------------------------------------------------------------------------------
// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and semi_infinite as appropriate and return true
bool Sphere::
Intersection(Ray& ray) const
{
  // TODO
  
  Vector_3D<double> u_vec = ray.direction.Normalized();
  Vector_3D<double> s_vec = ray.endpoint - center;
  double intersections = sqr(Vector_3D<double>::Dot_Product(s_vec, u_vec)) - s_vec.Length_Squared() + sqr(radius);
  
  if(intersections >= 0)
  {
    double length = -(Vector_3D<double>::Dot_Product(s_vec, u_vec) + sqrt(intersections));
    if(-(Vector_3D<double>::Dot_Product(s_vec, u_vec) - sqrt(intersections)) < length)
    {
      length = (Vector_3D<double>::Dot_Product(s_vec, u_vec) - sqrt(intersections));
    }
    cout << length << endl;
    if(length < small_t && length >= -small_t)
    {
      return false;
    }
    ray.t_max = length;
    ray.current_object = this;
    ray.semi_infinite = false;  
    return true;
  }

  return false;
}

Vector_3D<double> Sphere::
Normal(const Vector_3D<double>& location) const
{
    Vector_3D<double> normal;
    // TODO: set the normal
    normal = (location - center);
    normal.Normalize();
    return normal;
}

// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and semi_infinite as appropriate and return true
bool Plane::
Intersection(Ray& ray) const
{
  // TODO
 
  if(Vector_3D<double>::Dot_Product(ray.direction, normal) != 0)
  {
    ray.t_max = Vector_3D<double>::Dot_Product((x1-ray.endpoint),normal) / Vector_3D<double>::Dot_Product(ray.direction, normal);
    ray.current_object = this;
    ray.semi_infinite = false;
    return true;
    }
  
  return false;
}

Vector_3D<double> Plane::
Normal(const Vector_3D<double>& location) const
{
    return normal;
}
//--------------------------------------------------------------------------------
// Camera
//--------------------------------------------------------------------------------
// Find the world position of the input pixel
Vector_3D<double> Camera::
World_Position(const Vector_2D<int>& pixel_index)
{
  Vector_3D<double> result;
  // TODO 
  Vector_2D<double> grid2D = film.pixel_grid.X(pixel_index);
  result = focal_point + (vertical_vector * grid2D.y) + (horizontal_vector * grid2D.x);
  return result;
}
//--------------------------------------------------------------------------------
// Render_World
//--------------------------------------------------------------------------------
// Find the closest object of intersection and return a pointer to it
//   if the ray intersects with an object, then ray.t_max, ray.current_object, and ray.semi_infinite will be set appropriately
//   if there is no intersection do not modify the ray and return 0
int cnt = 0;
const Object* Render_World::
Closest_Intersection(Ray& ray)
{
  // TODO
  double min = 0.0;
  int closest_index = 0;
  for(int i = 0; i < objects.size(); i++)
  {
    if(objects[i] -> Intersection(ray))
    {
      /*if(i==1)
	{
	cnt++;
	cout << cnt << endl;
	}*/
	
      if(min == 0.0)
      {
	min = ray.t_max;
	closest_index = i;
      }
      else
      {
	if(ray.t_max < min)
	{
	  min = ray.t_max;
	  closest_index = i;
	}
      }
    }
  }
  
  if(min != 0.0)
  {
    ray.t_max = min;
    ray.current_object = objects[closest_index];
    ray.semi_infinite = false;
    return objects[closest_index];
  }


  return 0;
}

// set up the initial view ray and call 
void Render_World::
Render_Pixel(const Vector_2D<int>& pixel_index)
{
    // TODO
    Ray ray; // TODO: set up the initial view ray here
    ray.endpoint = camera.position;
    //ray.direction = camera.World_Position(pixel_index).Normalized();
    ray.direction = camera.World_Position(pixel_index)-camera.position;

    Ray dummy_root;
    Vector_3D<double> color=Cast_Ray(ray,dummy_root);
    camera.film.Set_Pixel(pixel_index,Pixel_Color(color));
}

// cast ray and return the color of the closest intersected surface point, 
// or the background color if there is no object intersection
Vector_3D<double> Render_World::
Cast_Ray(Ray& ray,const Ray& parent_ray)
{
  // TODO
  Vector_3D<double> color;
  // determine the color here
  if(Closest_Intersection(ray))
  {
    const Object * obj = Closest_Intersection(ray);
    color = obj -> material_shader -> Shade_Surface(ray, *obj, ray.Point(ray.t_max), obj -> Normal(ray.Point(ray.t_max)));
  }

  return color;
}
