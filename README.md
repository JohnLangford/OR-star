# OR-star
A relativistic solution to Einstein's equations and a ray tracer for that solution

Usage:
./OR-star &lt;impact&gt; <flux> <photon_radius> <photon_radial_velocity> [nobreak]

<impact> and <flux> define the OR-star while <photon_radius> and
<photon_radial_velocity> define the starting radius and angle of a
photon interacting with the OR-star.  Photons on the characteristic
geodesic must be handled carefully---the [nobreak] argument makes them
reflect inside while by default the simulation ends.

Example:
./OR-star 1 10 200 -0.999

A photon star of radius 1 (typical) and flux 10 (high) interacts with
a photon starting at radius 200 aimed near the origin.  
