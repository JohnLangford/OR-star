# OR-star
A relativistic solution to Einstein's equations and a ray tracer for that solution

Usage:
make
./OR-star &lt;impact&gt; &lt;flux&gt; &lt;photon_radius&gt; &lt;photon_radial_velocity&gt; [nobreak]

&lt;impact&gt; and &lt;flux&gt; define the OR-star while
&lt;photon_radius&gt; and &lt;photon_radial_velocity&gt; define the
starting radius and angle of a photon interacting with the OR-star.
Photons which deflect inwards either halt the simulation (as per
characteristic geodesic plots) or continue if there is a 5th argument.

Example:
./OR-star 1 10 200 -0.999

A photon star of radius 1 (typical) and flux 10 (high) interacts with
a photon starting at radius 200 aimed near the origin.  
