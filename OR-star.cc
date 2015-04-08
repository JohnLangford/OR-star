#include<iostream>
#include<math.h>
#include<stdlib.h>

using namespace std;

double sim_scale = 0.00001;
int clocks = 100000000;

struct ab {
  double A;
  double B;
};

struct rAB {
  double radius;
  ab geom;
};

rAB* geometry;

size_t max_geom;

double schwarz_r;
double schwarz_k;

double sign(double s)
{
  if (s < 0)
    return -1;
  else
    return 1;
}

ab find_nearest(double radius, size_t& hint)
{
  ab ret;

  if (radius < 1)
    {
      ret.A = 1;
      ret.B = 1;
      return ret;
    }
  if (radius > geometry[max_geom].radius)
    {//schwarzchild
      ret.A = schwarz_k * (1 - schwarz_r / radius);
      ret.B = 1. / (1. - schwarz_r / radius);
      return ret;
    }
  
  double delta_radius = fabs(geometry[hint].radius - radius);
  while(hint < max_geom) {
    double new_delta_radius = fabs(geometry[hint+1].radius - radius);
    if (new_delta_radius > delta_radius)
      break;
    else
      {
	delta_radius = new_delta_radius;
	hint = hint+1;
      }
  } 

  while(hint > 0) {
    double new_delta_radius = fabs(geometry[hint-1].radius - radius);
    if (new_delta_radius > delta_radius)
      break;
    else
      {
	delta_radius = new_delta_radius;
	hint = hint-1;
      }
  }

  size_t second_best;
  if (hint == 0)
    second_best = 1;
  else if (hint == max_geom)
    second_best = max_geom-1;
  else
    if (sign(geometry[hint-1].radius - radius) ==  sign(radius - geometry[hint].radius))
      second_best = hint-1;
    else
      second_best = hint+1;

  //radius is a mixure of hint radius and second_best radius.  What are the coefficients?
  // a * hint_radius + (1-a)*second_radius = radius
  // => a * (hint_radius - second_radius) = radius - second_radius
  // => a = (radius - second_radius) / (hint_radius - second_radius)
  float hint_mix = (radius - geometry[second_best].radius) / (geometry[hint].radius - geometry[second_best].radius);
  
  ret.A = hint_mix * geometry[hint].geom.A + (1 - hint_mix) * geometry[second_best].geom.A;
  ret.B = hint_mix * geometry[hint].geom.B + (1 - hint_mix) * geometry[second_best].geom.B;
  
  return ret;
}

double fourpi = 4. * 3.1415;

double m(rAB current)
{
  return current.radius*(current.geom.B - 1) / (2 * current.geom.B);
}

double dm(double K, rAB current)
{
  double A = current.geom.A;
  double mass = m(current);
  
  double numerator = fourpi * K * sqrt(current.radius - 2 * mass);
  double denominator = A * sqrt(current.radius * A);
  return numerator / denominator;
}

double dphi(double K, rAB current, double impact)
{
  double A = current.geom.A;
  double mass = m(current);
  double r = current.radius;

  double first_term = A * mass + fourpi * K * sqrt(r*r - impact * impact * A);
  double second_term = sqrt(current.geom.B) * sqrt(r*r - impact * impact * A);
  double denominator = r * r * r * A * sqrt(A);
  return (first_term * second_term) / denominator;
}

double dr(rAB current, double impact)
{
  double r = current.radius;  

  double numerator = sqrt(r*r - current.geom.A * impact * impact);
  double denominator = r * sqrt(current.geom.A * current.geom.B);
  return numerator / denominator;
}

double phi(rAB current)
{
  return 0.5 * log(current.geom.A);
}

struct photon {
  double radius;
  double origin_angle;
  double vr;
  double vperp;
};

int main(int argc, char* argv[])
{
  double impact = 1.;
  if (argc > 1)
    impact = atof(argv[1]);
  
  double e_at_origin = 0.81 / 2. / fourpi;
  if (argc > 2)
    e_at_origin = atof(argv[2]) / 2. / fourpi;
  
  rAB current = {impact + sim_scale/10000, {1., 1.}};
  
  geometry = (rAB*)calloc(clocks, sizeof(rAB));
  
  int i;
  for (i = 0; i < clocks; i++)
    {
      geometry[i] = current;
      
      double new_dm = sim_scale * dm(e_at_origin, current);
      double new_dphi = sim_scale * dphi(e_at_origin, current, impact);
      double new_dr = sim_scale * dr(current, impact);

      double new_r = current.radius + new_dr;
      double new_m = m(current) + new_dm;
      double new_B = new_r / (new_r - 2 * new_m);
      double new_phi = phi(current) + new_dphi;
      double new_A = exp(2. * new_phi);
      
      rAB new_rAB = {new_r, {new_A, new_B}};

      current = new_rAB;
      if ( impact * impact * current.geom.A / (current.radius * current.radius) > 1)
	{
	  cerr << "reached outer radius, total mass = " << m(current) << " total flux = " << i * sim_scale * e_at_origin << endl;
	  break;
	}
    }

  max_geom = min(i, clocks-1);
  schwarz_r = geometry[max_geom].radius * (1. - 1. / geometry[max_geom].geom.B);
  schwarz_k = geometry[max_geom].geom.A / (1. - schwarz_r / geometry[max_geom].radius);

  if (i == clocks)
    cerr << "ran out of clock ticks" << endl;

  photon ray = {impact, 0., 0., 1.};

  if (argc>3)
    ray.radius = atof(argv[3]);

  if (argc>4)
    {
      ray.vr = atof(argv[4]);
      ray.vperp = sqrt(1-ray.vr*ray.vr);
    }

  bool breaker = true;
  if (argc>5)
    breaker = false;

  size_t hint = 0;
  cout << "radius\t" << "angle\t" << "vr\t" << "vperp\t" << "x\t" << "y\t" << "A\t" << "B\t" << "hint\t" << endl;

  for (int i = 0; i < clocks; i++)
    {
      ab geom = find_nearest(ray.radius, hint);

      if (i % 100 == 0)
	cout << ray.radius << "\t" << ray.origin_angle << "\t" << ray.vr << "\t" << ray.vperp << 
	  "\t" << ray.radius * cos(ray.origin_angle) << 
	  "\t" << ray.radius * sin(ray.origin_angle) << 
	  "\t" << geom.A << 
	  "\t" << geom.B << 
	  "\t" << hint << endl;   
      
      double sqrt_A = sqrt(geom.A);

      double derived_impact = ray.vperp * ray.radius / sqrt_A;

      double delta_radius = ray.vr * sim_scale * sqrt_A;
      double delta_perp = ray.vperp * sim_scale * sqrt_A;
      
      double new_radius = sqrt((ray.radius+delta_radius)*(ray.radius+delta_radius) 
			       + delta_perp*delta_perp);
      ab geom_new = find_nearest(new_radius, hint);
      
      double new_vperp = derived_impact * sqrt(geom_new.A) / new_radius;
      double s_radius = sign(new_radius - ray.radius);

      double new_vr;
      if (1. < new_vperp * new_vperp) // This is not possible with infinite precision.  It can only happen when vr is extremely small.  We assume that this is a turnaround point.
	{
	  if (breaker)
	    {	  
	      cerr << "reached outer radius" << endl;
	      break;
	    }
	  new_vr = - ray.vr;
	  delta_radius = new_vr * sim_scale * sqrt_A;
	  new_radius = sqrt((ray.radius+delta_radius)*(ray.radius+delta_radius) 
			    + delta_perp*delta_perp);
	  geom_new = find_nearest(new_radius, hint);
	  new_vperp = derived_impact * sqrt(geom_new.A) / new_radius;
	}
      else
	new_vr = s_radius * sqrt (1. - new_vperp * new_vperp);

      double new_origin_angle = ray.origin_angle + asin(delta_perp / ray.radius);

      photon new_ray = {new_radius, new_origin_angle, new_vr, new_vperp};

      if (new_radius > 400)
	break;
      ray = new_ray;
      }
}
