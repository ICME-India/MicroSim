#ifndef FILLING_H_
#define FILLING_H_

#include <time.h>
#include <gsl/gsl_rng.h>

int checkOverlap(long centx, long centy, long centz, double rad, long shield_dist, int *occupancy) {
  long x,y,z;
    if (!shield_dist)
        return 0;

    for (x = 0; x < rows_x; x++)
        for (y = 0; y < rows_y; y++)
            for (z = 0; z < rows_z; z++)
                if (occupancy[x*layer_size + z*rows_y + y] == 1)
                    if ((x - centx)*(x - centx) + (y - centy)*(y - centy)
                        + (z - centz)*(z - centz) <= (shield_dist*rad)*(shield_dist*rad))
                    return 1;
    return 0;
}

void fill_phase_cube (struct fill_cube fill_cube_parameters, struct fields* gridinfo, long b) {
  long x, y, z, index;
  long a;
  double sum;

  for(x=0;x < rows_x; x++) {
    for(z=0; z < rows_z; z++) {
      for (y=0; y < rows_y; y++) {
        index = x*layer_size + z*rows_y + y;

        gridinfo[index].deltaphi[b] = 0.0;

        if (b < (NUMPHASES-1)) {
          if ((x >= fill_cube_parameters.x_start) && (x <= fill_cube_parameters.x_end)
               && (y >=fill_cube_parameters.y_start) && (y<=fill_cube_parameters.y_end)
               && (z >=fill_cube_parameters.z_start) && (z<=fill_cube_parameters.z_end)) {

            gridinfo[index].phia[b] = 1.00000;
            for (a=0; a < NUMPHASES; a++) {
              if (b!=a) {
              gridinfo[index].phia[a] = 0.00000;
              }
            }
          } else {
            if (gridinfo[index].phia[b] != 1.0000) {
              gridinfo[index].phia[b] = 0.00000;
            }
          }
        } else {
          sum = 0.0;
          for (a=0; a < NUMPHASES-1; a++) {
            sum += gridinfo[index].phia[a];
          }
          if (sum > 1.0) {
            printf("Wrong filling operation, will fill it with liquid\n");
            gridinfo[index].phia[b] = 1.0;
            for (a=0; a < NUMPHASES-1; a++) {
              if (a!=b) {
                gridinfo[index].phia[a] = 0.00000;
              }
            }
          } else {
            gridinfo[index].phia[b]   = 1.0-sum;
          }
        }
      }
    }
  }
}


void fill_cube_pattern(long variants, long sx, long sy, long sz,
                       double sfrac, long gap, double gfrac)
{
    /* Randomly fill multiple variants of square/cubic particles
       in an array.
       NOTE: the last phase is assumed to be matrix. */
    ldiv_t resx, resy, resz;
    long i, j, k, index, sgn, nparticles;
    long xlo=0, ylo=0, zlo=0;
    long xhi=1, yhi=1, zhi=1;
    double r;
    gsl_rng *rng;

    // Set up GSL RNG.
    gsl_rng_env_setup();
    rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, time(0));

    if ( NUMPHASES > 1 )
    {
        // First, make matrix = 1 everywhere.
        for ( long x=0; x<rows_x; x++ )
        {
            for ( long y=0; y<rows_y; y++ )
            {
                for ( long z=0; z<rows_z; z++ )
                {
                    index = x*layer_size + z*rows_y + y;
                    gridinfo[index].phia[NUMPHASES-1] = 1.0;
                }
            }
        }
    }

    resx = ldiv(MESH_X, sx+gap);
    resy = ldiv(MESH_Y, sy+gap);
    resz = ldiv(MESH_Z, sz+gap);
    if ( resx.quot == 0 )
        resx.quot = 1;
    if ( resy.quot == 0 )
        resy.quot = 1;
    if ( resz.quot == 0 )
        resz.quot = 1;
    nparticles = resx.quot * resy.quot * resz.quot;
    printf("Filling %ld particles.\n", nparticles);
    for ( i=0; i<resx.quot; i++ )
    {
        if ( MESH_X > 1 )
        {
            r = gsl_rng_uniform(rng);
            sgn = 2*lround(r) - 1;
            xlo = i*(gap+sx) + gap + sgn*gfrac*gap - (2*r-1)*sfrac*sx;
            xhi = (i+1) * (gap+sx) + sgn*gfrac*gap + (2*r-1)*sfrac*sx;
        }
        for ( j=0; j<resy.quot; j++ )
        {
            if ( MESH_Y > 1 )
            {
                r = gsl_rng_uniform(rng);
                sgn = 2*lround(r) - 1;
                ylo = j*(gap+sy) + gap + sgn*gfrac*gap - (2*r-1)*sfrac*sy;
                yhi = (j+1) * (gap+sy) + sgn*gfrac*gap + (2*r-1)*sfrac*sy;
            }
            for (k=0; k<resz.quot; k++ )
            {
                if ( MESH_Z > 1 )
                {
                    r = gsl_rng_uniform(rng);
                    sgn = 2*lround(r) - 1;
                    zlo = k*(gap+sz) + gap + sgn*gfrac*gap - (2*r-1)*sfrac*sz;
                    zhi = (k+1) * (gap+sz) + sgn*gfrac*gap + (2*r-1)*sfrac*sz;
                }
                fill_cube_parameters.x_start = xlo;
                fill_cube_parameters.x_end   = xhi;
                fill_cube_parameters.y_start = ylo;
                fill_cube_parameters.y_end   = yhi;
                fill_cube_parameters.z_start = zlo;
                fill_cube_parameters.z_end   = zhi;
                if ( variants >= NUMPHASES )
                    variants = NUMPHASES - 1;
                r = gsl_rng_uniform(rng);
                long phase = r * variants;
                printf("phase now = %ld\n", phase);
                fill_phase_cube(fill_cube_parameters, gridinfo, phase);
            }
        }
    }
    // Free GSL RNG.
    gsl_rng_free(rng);
}


int check_overlap_sq(struct fill_cube cube, long shield, long vol)
{
    /* Check if the to-be created precipitate is overlapping an existing one. */
    long index;
    for ( long x=(cube.x_start-shield); x<(cube.x_end+shield); x++ )
    {
        for ( long y=(cube.y_start-shield); y<(cube.y_end+shield); y++ )
        {
            for ( long z=(cube.z_start-shield); z<(cube.z_end+shield); z++ )
            {
                index = x*layer_size + z*rows_y + y;
                if ( (0 <= index) && (index <= vol) )
                {
                    if ( NUMPHASES > 1 )
                    {
                        // For multiple phases, check if phi(matrix) is zero.
                        if ( !(gridinfo[index].phia[NUMPHASES-1] > 0) )
                            return 1;
                    }
                    else
                    {
                        // Else, perform usual check, i.e. if phi is one.
                        if ( gridinfo[index].phia[NUMPHASES-1] > 0 )
                            return 1;
                    }
                }
            }
        }
    }
    return 0;
}


void fill_phase_cube_random_variants(long variants, long sx, long sy, long sz,
                                     double sfrac, double vf, long shield)
{
    /* Randomly fill multiple variants of square/cubic particles
       with a given volume fraction and size.
       NOTE: (1) The last phase is assumed to be matrix.
             (2) `vf` is the total volume fraction. */
    long index, sgn;
    double r;
    gsl_rng *rng;

    // Set up GSL RNG.
    gsl_rng_env_setup();
    rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, time(0));

    if ( NUMPHASES > 1 )
    {
        // First, make matrix = 1 everywhere.
        for ( long x=0; x<rows_x; x++ )
        {
            for ( long y=0; y<rows_y; y++ )
            {
                for ( long z=0; z<rows_z; z++ )
                {
                    index = x*layer_size + z*rows_y + y;
                    gridinfo[index].phia[NUMPHASES-1] = 1.0;
                }
            }
        }
    }

    long volume_domain = MESH_X * MESH_Y * MESH_Z;
    double volume_per_particle = sx * sy * sz;
    int num_particles = ceil(volume_domain*vf/volume_per_particle);
    printf("Domain volume = %ld, sx = %ld, sy = %ld, sz = %ld,"
           " volume_per_particle = %lf, number of particles = %d\n",
           volume_domain, sx, sy, sz, volume_per_particle, num_particles);

    long random_count = 0, random_limit = 1e+4;
    while ( num_particles )
    {
        printf("num_particles = %d\n", num_particles);
        r = gsl_rng_uniform(rng);
        long sxnow = sx + (2*r - 1) * sfrac*sx;
        long synow = sy + (2*r - 1) * sfrac*sy;
        long sznow = sz + (2*r - 1) * sfrac*sz;
        long xnow = MESH_X * gsl_rng_uniform(rng);
        long ynow = MESH_Y * gsl_rng_uniform(rng);
        long znow = MESH_Z * gsl_rng_uniform(rng);
        long xlo = xnow - sxnow/2;
        long xhi = xnow + sxnow/2;
        long ylo = ynow - synow/2;
        long yhi = ynow + synow/2;
        long zlo = znow - sznow/2;
        long zhi = znow + sznow/2;
        printf("xlo = %ld, xhi = %ld\n", xlo, xhi);
        printf("ylo = %ld, yhi = %ld\n", ylo, yhi);
        printf("zlo = %ld, zhi = %ld\n", zlo, zhi);
        fill_cube_parameters.x_start = xlo;
        fill_cube_parameters.x_end = xhi;
        fill_cube_parameters.y_start = ylo;
        fill_cube_parameters.y_end = yhi;
        fill_cube_parameters.z_start = zlo;
        fill_cube_parameters.z_end = zhi;
        // Skip if this particle is overlapping an existing particle.
        if ( check_overlap_sq(fill_cube_parameters, shield, volume_domain) )
            continue;
        // Else, create the particle.
        if ( variants >= NUMPHASES )
            variants = NUMPHASES - 1;
        r = gsl_rng_uniform(rng);
        long phase = r * variants;
        printf("phase now = %ld\n", phase);
        fill_phase_cube(fill_cube_parameters, gridinfo, phase);
        --num_particles;
        ++random_count;
        if (random_count > random_limit)
        {
            printf("Warning: random filling limit = %ld reached."
                   " %d particles could not be filled.\n",
                   random_limit, num_particles);
            break;
        }
    }
    // Free GSL RNG.
    gsl_rng_free(rng);
}


void fill_phase_ellipse (struct fill_ellipse fill_ellipse_parameters, struct fields* gridinfo, long b) {
  long x, y, z, gidy;
  long a;
  double x_;
  double y_;
  double sum;
  long x_center, y_center, z_center;
  double l, e, angle;
  double angle_rad = angle*M_PI/180.0;

  x_center =  fill_ellipse_parameters.x_center;
  y_center =  fill_ellipse_parameters.y_center;
  z_center =  fill_ellipse_parameters.z_center;
  l        =  fill_ellipse_parameters.major_axis;
  e        =  fill_ellipse_parameters.eccentricity;
  angle    =  fill_ellipse_parameters.rot_angle;

  for(x=0;x<=(MESH_X-1);x++) {
    for (y=0; y<=(MESH_Y-1); y++) {
      if (b < (NUMPHASES-1)) {
        if ((double)(x-x_center)*(x-x_center)/(l*l) + ((double)(y-y_center)*(y-y_center))/(e*e*l*l) <= 1.0) {
          x_ =  x_center + floor((double)(x-x_center)*cos(angle_rad)  + (double)(y-y_center)*sin(angle_rad));
          y_ =  y_center + floor(-(double)(x-x_center)*sin(angle_rad) + (double)(y-y_center)*cos(angle_rad));
        } else {
          x_ = x;
          y_ = y;
        }
        gidy = x_*MESH_Y + y_;
        if ((gidy < (MESH_X*MESH_Y)) && (gidy > 0)) {
          if ((double)(x-x_center)*(x-x_center)/(l*l) + ((double)(y-y_center)*(y-y_center))/(e*e*l*l) <= 1.0) {
            gridinfo[gidy].phia[b] = 1.00000;
            for (a=0; a < NUMPHASES; a++) {
              if (b!=a) {
              gridinfo[gidy].phia[a] = 0.00000;
              }
            }
          } else {
            if (gridinfo[gidy].phia[b] != 1.0000) {
              gridinfo[gidy].phia[b] = 0.00000;
            }
          }
        }
      } else {
        gidy = x*MESH_Y + y;
        sum = 0.0;
        for (a=0; a < NUMPHASES-1; a++) {
          sum += gridinfo[gidy].phia[a];
        }
        if (sum > 1.0) {
          printf("Wrong filling operation, will fill it with liquid\n");
          gridinfo[gidy].phia[b] = 1.0;
          for (a=0; a < NUMPHASES-1; a++) {
            if (a!=b) {
             gridinfo[gidy].phia[a] = 0.00000;
            }
          }
        } else {
          gridinfo[gidy].phia[b]   = 1.0-sum;
        }
      }
    }
  }
}
void fill_phase_cylinder(struct fill_cylinder fill_cylinder_parameters, struct fields* gridinfo, long b) {
  long x, y, z, index;
  long a;
  double sum;

  long x_center, y_center, z_start, z_end;
  double radius;

  x_center = fill_cylinder_parameters.x_center;
  y_center = fill_cylinder_parameters.y_center;
  z_start  = fill_cylinder_parameters.z_start;
  z_end    = fill_cylinder_parameters.z_end;
  radius   = fill_cylinder_parameters.radius;

  for(x=0;x < rows_x; x++) {
    for (z=0; z < rows_z; z++) {
      for (y=0; y < rows_y; y++) {
        index = x*layer_size + z*rows_y + y;
        if (b < (NUMPHASES-1)) {
          if (((x-x_center)*(x-x_center) + (y-y_center)*(y-y_center) <= radius*radius) && (z>=z_start) && (z<=z_end)) {
            gridinfo[index].phia[b] = 1.00000;
            for (a=0; a < NUMPHASES; a++) {
              if (b!=a) {
              gridinfo[index].phia[a] = 0.00000;
              }
            }
          } else {
            if (gridinfo[index].phia[b] != 1.0000) {
              gridinfo[index].phia[b] = 0.00000;
            }
          }
        } else {
          sum = 0.0;
          for (a=0; a < NUMPHASES-1; a++) {
            sum += gridinfo[index].phia[a];
          }
          if (sum > 1.0) {
            printf("Wrong filling operation, will fill it with liquid\n");
            gridinfo[index].phia[b] = 1.0;
            for (a=0; a < NUMPHASES-1; a++) {
              if (a!=b) {
              gridinfo[index].phia[a] = 0.00000;
              }
            }
          } else {
            gridinfo[index].phia[b]   = 1.0-sum;
          }
        }
      }
    }
  }
}
void fill_phase_sphere(struct fill_sphere fill_sphere_parameters, struct fields* gridinfo, long b) {
  long x, y, z, index;
  long a;
  double sum;

  long x_center, y_center, z_center;
  double radius;

  x_center = fill_sphere_parameters.x_center;
  y_center = fill_sphere_parameters.y_center;
  z_center = fill_sphere_parameters.z_center;
  radius   = fill_sphere_parameters.radius;

  for(x=0;x < rows_x; x++) {
    for (z=0; z < rows_z; z++) {
      for (y=0; y < rows_y; y++) {
        index = x*layer_size + z*rows_y + y;
        if (b < (NUMPHASES-1)) {
          if (((x-x_center)*(x-x_center) + (y-y_center)*(y-y_center) + (z- z_center)*(z-z_center) <= radius*radius)) {
            gridinfo[index].phia[b] = 1.00000;
            for (a=0; a < NUMPHASES; a++) {
              if (b!=a) {
              gridinfo[index].phia[a] = 0.00000;
              }
            }
          } else {
            if (gridinfo[index].phia[b] != 1.0000) {
              gridinfo[index].phia[b] = 0.00000;
            }
          }
        } else {
          sum = 0.0;
          for (a=0; a < NUMPHASES-1; a++) {
            sum += gridinfo[index].phia[a];
          }
          if (sum > 1.0) {
            printf("Wrong filling operation, will fill it with liquid\n");
            gridinfo[index].phia[b] = 1.0;
            for (a=0; a < NUMPHASES-1; a++) {
              if (a!=b) {
              gridinfo[index].phia[a] = 0.00000;
              }
            }
          } else {
            gridinfo[index].phia[b]   = 1.0-sum;
          }
        }
      }
    }
  }
}
void fill_phase_sphere_occupancy(struct fill_sphere fill_sphere_parameters, struct fields* gridinfo, long b, int* occupancy) {
  long x, y, z, index;
  long a;
  double sum;

  long x_center, y_center, z_center;
  double radius;

  x_center = fill_sphere_parameters.x_center;
  y_center = fill_sphere_parameters.y_center;
  z_center = fill_sphere_parameters.z_center;
  radius   = fill_sphere_parameters.radius;

  for(x=0;x < rows_x; x++) {
    for (z=0; z < rows_z; z++) {
      for (y=0; y < rows_y; y++) {
        index = x*layer_size + z*rows_y + y;
        if (b < (NUMPHASES-1)) {
          if (((x-x_center)*(x-x_center) + (y-y_center)*(y-y_center) + (z- z_center)*(z-z_center) <= radius*radius)) {
            gridinfo[index].phia[b] = 1.00000;
             occupancy[index] = 1;
            for (a=0; a < NUMPHASES; a++) {
              if (b!=a) {
              gridinfo[index].phia[a] = 0.00000;
              }
            }
          } else {
            if (gridinfo[index].phia[b] != 1.0000) {
              gridinfo[index].phia[b] = 0.00000;
            }
          }
        } else {
          sum = 0.0;
          for (a=0; a < NUMPHASES-1; a++) {
            sum += gridinfo[index].phia[a];
          }
          if (sum > 1.0) {
            printf("Wrong filling operation, will fill it with liquid\n");
            gridinfo[index].phia[b] = 1.0;
            for (a=0; a < NUMPHASES-1; a++) {
              if (a!=b) {
              gridinfo[index].phia[a] = 0.00000;
              }
            }
          } else {
            gridinfo[index].phia[b]   = 1.0-sum;
          }
        }
      }
    }
  }
}

void fill_phase_sphere_random(long phase, double ppt_radius, double volume_fraction, long shield_dist, double spread) {
  int overlap = 1;
  int *occupancy;
  long i;
  occupancy = (int*)malloc(sizeof(int)*rows_x*rows_y*rows_z);

  for (i = 0; i < rows_x*rows_y*rows_z; i++) {
    occupancy[i] = 0;
  }
  clock_t SEED = clock();
  double volume_domain = (double) MESH_X*MESH_Y*MESH_Z;
  double volume_per_particle = (double)(ppt_radius*ppt_radius*ppt_radius)*(4.0/3.0)*M_PI;

  int num_particles = ceil(volume_domain*volume_fraction/volume_per_particle);
  printf("Domain volume = %lf, ppt_radius = %le, volume_per_particle = %lf, Number of particles = %d\n", volume_domain, ppt_radius, volume_per_particle, num_particles);

  long particle_index = 1, random_count = 0, random_limit = 1e+5;

  struct fill_sphere temp_sph;

  while (particle_index <= num_particles) {

      long centx = MESH_X*drand48();
      long centy = MESH_Y*drand48();
      long centz = MESH_Z*drand48();

      temp_sph.x_center = centx;
      temp_sph.y_center = centy;
      temp_sph.z_center = centz;

      double mdev = spread*ppt_radius;
      double rad  = (double)ppt_radius + (1.0*drand48() - 0.5)*mdev;

      temp_sph.radius = rad;
      //printf("particle_index = %d, x_center = %ld, y_center = %ld, z_start = %ld, z_end = %ld, radius = %lf\n", particle_index, centx, centy, temp_cyl.z_start, temp_cyl.z_end, rad);

      random_count++;
      if (random_count > random_limit) {
        printf("Random filling attempt count limit = %ld reached. Only %ld particles could be filled.\n", random_limit, particle_index-1);
        break;
      }

      if (((centx - rad) < 0.0) || ((centx + rad) > (MESH_X-1)) ||
          ((centy - rad) < 0.0) || ((centy + rad) > (MESH_Y-1)) ||
          ((centz - rad) < 0.0) || ((centz + rad) > (MESH_Z-1)) )
          continue;

      overlap = checkOverlap(centx, centy, centz, rad, shield_dist, occupancy);

      if(overlap == 0) {
        fill_phase_sphere_occupancy(temp_sph, gridinfo, phase, occupancy);
        //fill_phase_cylinder(temp_cyl, gridinfo, b);
      }
      else {
        overlap = 1;
        continue;
      }
      particle_index++;
  }
  fill_phase_sphere_occupancy(temp_sph, gridinfo, NUMPHASES-1, occupancy);
  free(occupancy);
}

void fill_phase_cylinder_occupancy(struct fill_cylinder fill_cylinder_parameters, struct fields* gridinfo, long b, int* occupancy) {
    long x, y, z, index;
    long a;
    double sum;

    long x_center, y_center, z_start, z_end;
    double radius;

    x_center = fill_cylinder_parameters.x_center;
    y_center = fill_cylinder_parameters.y_center;
    z_start  = fill_cylinder_parameters.z_start;
    z_end    = fill_cylinder_parameters.z_end;
    radius   = fill_cylinder_parameters.radius;

    for(x = 0; x < rows_x; x++) {
      for (y = 0; y < rows_y; y++) {
        for (z = 0; z < rows_z; z++) {
            index = x*layer_size + z*rows_y + y;
            if ((z >= z_start) && (z <= z_end)) {
              if (b < (NUMPHASES-1)) {
                if (((x-x_center)*(x-x_center) + (y-y_center)*(y-y_center) <= radius*radius)) {
                  gridinfo[index].phia[b] = 1.00000;
                  occupancy[index] = 1;
                  for (a = 0; a < NUMPHASES; a++) {
                    if (b != a) {
                        gridinfo[index].phia[a] = 0.00000;
                    }
                  }
                } else if (gridinfo[index].phia[b] != 1.0000) {
                    gridinfo[index].phia[b] = 0.00000;
                }
            } else {
              sum = 0.0;
              for (a=0; a < NUMPHASES-1; a++) {
                sum += gridinfo[index].phia[a];
              }
              if (sum > 1.0) {
                printf("Wrong filling operation, will fill it with liquid\n");
                gridinfo[index].phia[b] = 1.0;
                for (a=0; a < NUMPHASES-1; a++) {
                  if (a!=b) {
                    gridinfo[index].phia[a] = 0.00000;
                  }
                }
              } else {
                gridinfo[index].phia[b]   = 1.0-sum;
              }
           }
        }
      }
    }
  }
}
void fill_phase_cylinder_random(long phase, double ppt_radius, double volume_fraction, long shield_dist, double spread) {
    int overlap = 1;
    int *occupancy;
    long i;
    occupancy = (int*)malloc(sizeof(int)*rows_x*rows_y*rows_z);

    for (i = 0; i < rows_x*rows_y*rows_z; i++) {
        occupancy[i] = 0;
    }

    clock_t SEED = clock();
    double volume_domain = (double) MESH_X*MESH_Y*MESH_Z;
    double volume_per_particle = (double)(ppt_radius*ppt_radius)*M_PI;

    int num_particles = ceil(volume_domain*volume_fraction/volume_per_particle);
    printf("Domain volume = %lf, ppt_radius = %le, volume_per_particle = %lf, Number of particles = %d\n", volume_domain, ppt_radius, volume_per_particle, num_particles);

    long particle_index = 1, random_count = 0, random_limit = 1e+4;

    struct fill_cylinder temp_cyl;
    temp_cyl.z_start  = 0;
    temp_cyl.z_end    = 0;

    while (particle_index <= num_particles)
    {
        long centx = MESH_X*drand48();
        long centy = MESH_Y*drand48();
        long centz = 0;

        temp_cyl.x_center = centx;
        temp_cyl.y_center = centy;

        double mdev = spread*ppt_radius;
        double rad  = (double)ppt_radius + (1.0*drand48() - 0.5)*mdev;

        temp_cyl.radius = rad;
//         printf("particle_index = %d, x_center = %ld, y_center = %ld, z_start = %ld, z_end = %ld, radius = %lf\n", particle_index, centx, centy, temp_cyl.z_start, temp_cyl.z_end, rad);

        random_count++;
        if (random_count > random_limit)
        {
            printf("Random filling attempt count limit = %ld reached. Only %ld particles could be filled.\n", random_limit, particle_index-1);
            break;
        }

        if (((centx - rad) < 0.0) || ((centx + rad) > (MESH_X-1)) ||
            ((centy - rad) < 0.0) || ((centy + rad) > (MESH_Y-1)) )
            continue;

        overlap = checkOverlap(centx, centy, centz, rad, shield_dist, occupancy);

        if(overlap == 0)
        {
            fill_phase_cylinder_occupancy(temp_cyl, gridinfo, phase, occupancy);
//             fill_phase_cylinder_occupancy(temp_cyl, gridinfo, 0, occupancy);
//             fill_phase_cylinder(temp_cyl, gridinfo, b);
        }
        else
        {
            overlap = 1;
            continue;
        }

        particle_index++;
    }
    fill_phase_cylinder_occupancy(temp_cyl, gridinfo, NUMPHASES-1, occupancy);
    free(occupancy);
}
void fill_phase_voronoi_2D(struct fill_cube fill_cube_parameters, struct fields* gridinfo, long NUMPOINTS_VORONOI, double size_min) {
  long x, y, gidy, gidy1;
  long limit_x, limit_y;
//   static long A = 1000;
  int k, location, s;
//   double n[NUMPOINTS_VORONOI], m[NUMPOINTS_VORONOI], l[NUMPOINTS_VORONOI], minimum;
//   long phase[NUMPOINTS_VORONOI];
//   int FLAG[rows_x*rows_y];
  double *n, *m, *p, *l, minimum;
  long *phase;
  int *FLAG;

  FLAG  = (int*)malloc(sizeof(int)*rows_x*rows_y);
  n     = (double *)malloc(sizeof(double)*NUMPOINTS_VORONOI);
  m     = (double *)malloc(sizeof(double)*NUMPOINTS_VORONOI);
  p     = (double *)malloc(sizeof(double)*NUMPOINTS_VORONOI);
  l     = (double *)malloc(sizeof(double)*NUMPOINTS_VORONOI);
  phase = (long *)malloc(sizeof(long)*NUMPOINTS_VORONOI);

  long rand_x, rand_y;
  //static long size_min=400;
  int PHASE_FILLED=0;

  limit_x =  fill_cube_parameters.x_end - fill_cube_parameters.x_start;
  limit_y =  fill_cube_parameters.y_end - fill_cube_parameters.y_start;

  for(k=0;k<NUMPOINTS_VORONOI;k++) {
    while (PHASE_FILLED!=1) {
      rand_x = fill_cube_parameters.x_start + (drand48())*limit_x;
      rand_y = fill_cube_parameters.y_start + (drand48())*limit_y;
      gidy = rand_x*rows_y + rand_y;
      if (FLAG[gidy]!=1) {
        n[k] = rand_x;
        m[k] = rand_y;
        for(x=fill_cube_parameters.x_start;x<fill_cube_parameters.x_end;x++){
         for (y=fill_cube_parameters.y_start; y<fill_cube_parameters.y_end; y++) {
          if ((m[k]-y)*(m[k]-y) +(n[k]-x)*(n[k]-x) <= size_min*size_min) {
            gidy1 = x*rows_y + y;
            FLAG[gidy1] = 1;
          }
        }
      }
      PHASE_FILLED=1;
    }
   }
   PHASE_FILLED=0;
  }
  for(k=0;k<NUMPOINTS_VORONOI;k++) {
    phase[k] = lrand48()%(NUMPHASES-1);
  }
  for(x=fill_cube_parameters.x_start;x<fill_cube_parameters.x_end;x++) {
    for (y=fill_cube_parameters.y_start; y<fill_cube_parameters.y_end; y++) {
      gidy   = x*rows_y + y;
      for(k=0; k<NUMPOINTS_VORONOI; k++) {
        l[k] = (n[k]-x)*(n[k]-x) + (m[k]-y)*(m[k]-y);
      }
      minimum  = l[0];
      location = 0;

      for(k=0;k<NUMPOINTS_VORONOI;k++) {
        if(l[k]<minimum) {
          minimum  = l[k];
          location = phase[k];
        }
      }
      for(s=0;s<NUMPHASES-1;s++) {
        if(location == s) {
          gridinfo[gidy].phia[s] = 1.0;
        }
        if(location != s) {
          gridinfo[gidy].phia[s] = 0.0;
        }
      }
    }
  }
  fill_phase_cube(fill_cube_parameters, gridinfo, NUMPHASES-1);
  free(FLAG);
  free(n);
  free(m);
  free(p);
  free(l);
  free(phase);
}
void fill_phase_voronoi_3D(struct fill_cube fill_cube_parameters, struct fields* gridinfo, long NUMPOINTS_VORONOI, double size_min) {
  long x, y, z, gidy, gidy1;
  long limit_x, limit_y, limit_z;
//   static long A = 1000;
  int k, location, s;
//   double n[NUMPOINTS_VORONOI], m[NUMPOINTS_VORONOI], p[NUMPOINTS_VORONOI], l[NUMPOINTS_VORONOI], minimum;
//   long phase[NUMPOINTS_VORONOI];

  double *n, *m, *p, *l, minimum;
  long *phase;

//   int FLAG[rows_x*rows_y*rows_z];
  int *FLAG;
  FLAG  = (int*)malloc(sizeof(int)*rows_x*rows_y*rows_z);
  n     = (double *)malloc(sizeof(double)*NUMPOINTS_VORONOI);
  m     = (double *)malloc(sizeof(double)*NUMPOINTS_VORONOI);
  p     = (double *)malloc(sizeof(double)*NUMPOINTS_VORONOI);
  l     = (double *)malloc(sizeof(double)*NUMPOINTS_VORONOI);
  phase = (long *)malloc(sizeof(long)*NUMPOINTS_VORONOI);
  long rand_x, rand_y, rand_z;
  //static long size_min=400;
  int PHASE_FILLED=0;

  limit_x =  fill_cube_parameters.x_end - fill_cube_parameters.x_start;
  limit_y =  fill_cube_parameters.y_end - fill_cube_parameters.y_start;
  limit_z =  fill_cube_parameters.z_end - fill_cube_parameters.z_start;

  for(k=0;k<NUMPOINTS_VORONOI;k++) {
    while (PHASE_FILLED!=1) {
      rand_x = fill_cube_parameters.x_start + (drand48())*(limit_x-1);
      rand_y = fill_cube_parameters.y_start + (drand48())*(limit_y-1);
      rand_z = fill_cube_parameters.z_start + (drand48())*(limit_z-1);

      gidy = rand_x*layer_size + rand_z*rows_y + rand_y;

      if (FLAG[gidy]!=1) {
        n[k] = rand_x;
        m[k] = rand_y;
        p[k] = rand_z;
        for(x=fill_cube_parameters.x_start;x < fill_cube_parameters.x_end;x++) {
         for (y=fill_cube_parameters.y_start; y < fill_cube_parameters.y_end; y++) {
           for (z=fill_cube_parameters.z_start; z < fill_cube_parameters.z_end; z++) {
            if ((m[k]-y)*(m[k]-y) + (n[k]-x)*(n[k]-x) + (p[k]-z)*(p[k]-z) <= size_min*size_min) {
              gidy1 = x*layer_size + z*rows_y + y;
              FLAG[gidy1] = 1;
            }
          }
        }
      }
      PHASE_FILLED=1;
    }
   }
   PHASE_FILLED=0;
  }

  for(k=0;k<NUMPOINTS_VORONOI;k++) {
    phase[k] = lrand48()%(NUMPHASES-1);
  }

  for(x=fill_cube_parameters.x_start; x < fill_cube_parameters.x_end; x++) {
    for (y=fill_cube_parameters.y_start; y < fill_cube_parameters.y_end; y++) {
      for (z=fill_cube_parameters.z_start; z < fill_cube_parameters.z_end; z++) {
        gidy   = x*layer_size + z*rows_y + y;
        for(k=0; k<NUMPOINTS_VORONOI; k++) {
          l[k] = (n[k]-x)*(n[k]-x) + (m[k]-y)*(m[k]-y) + (p[k]-z)*(p[k]-z);
        }
        minimum  = l[0];
        location = 0;

        for(k=0;k<NUMPOINTS_VORONOI;k++) {
          if(l[k]<minimum) {
            minimum  = l[k];
            location = phase[k];
          }
        }
        for(s=0;s<NUMPHASES-1;s++) {
          if(location == s) {
            gridinfo[gidy].phia[s] = 1.0;
          }
          if(location != s) {
            gridinfo[gidy].phia[s] = 0.0;
          }
        }
      }
    }
  }
  fill_phase_cube(fill_cube_parameters, gridinfo, NUMPHASES-1);
  free(FLAG);
  free(n);
  free(m);
  free(p);
  free(l);
  free(phase);
}
void fill_composition_cube(struct fields* gridinfo) {
  long x, y, z, index;
  long k;
  long b;
  double c[NUMCOMPONENTS-1];
  double chemical_potential;
  long PHASE_FILLED=0;

  for(x=0;x<rows_x;x++) {
    for(z=0; z<rows_z; z++) {
      for (y=0; y<rows_y; y++) {
        index = x*layer_size + z*rows_y + y;
        PHASE_FILLED = 0;
        for (b=0; b < NUMPHASES-1; b++) {
          if (gridinfo[index].phia[b] == 1.0) {
            for (k=0; k < NUMCOMPONENTS-1; k++) {
              c[k] = ceq[b][b][k];
            }
//             init_propertymatrices(Teq);
            Mu(c, Teq, b, gridinfo[index].compi);
            for (k=0; k < NUMCOMPONENTS-1; k++) {
//               chemical_potential = Mu(c, Teq, b, k);
//               gridinfo[index].compi[k] = chemical_potential;
              gridinfo[index].composition[k] = c[k];
  //             printf("solid_chemical_potential=%le\n",gridinfo[gidy].compi[k]);
            }
            PHASE_FILLED =1;
            break;
          }
        }
        if (!PHASE_FILLED) {
          //Fill with liquid
          for (k=0; k < NUMCOMPONENTS-1; k++) {
            c[k] = cfill[NUMPHASES-1][NUMPHASES-1][k];
          }
  //         c[0] = 0.160413;
  //         c[1] = 0.245857;
//           init_propertymatrices(Teq);
          Mu(c, Teq, NUMPHASES-1, gridinfo[index].compi);
          for (k=0; k < NUMCOMPONENTS-1; k++) {
//             chemical_potential = Mu(c, Teq, NUMPHASES-1, k);
// //             printf("chemical_potential =%le\n", chemical_potential);
//             gridinfo[index].compi[k] = chemical_potential;
             gridinfo[index].composition[k] = c[k];
          }
        }
      }
    }
  }
//   exit(0);
//   init_propertymatrices(T);
}
#endif





// void initdomain_struct(struct variables *gridinfo) {
//   long x,y;
//   double chemical_potential;
//   FILE *fp;
//   long b, k;
//
//   char name[1000];
//
//   //Filling algorithm for the phases---------------------------------------------------------------
//   long cubes;
//   long start_cube;
//
//   char pattern[10000];
//   sprintf(pattern, "01");
//   long pattern_length = 2;
//   long period = (MESH_X-1);
//   long numperiods=1;
//   //Filling cubes of phases
//   //Filling abg configuration for half the domain--------------------------------------------------
//   start_cube = 0;
//   //Determine the number of cubes of each phase in the pattern
//   find_numcubes(pattern, filling_type_phase);
//   long lambda = period/numperiods;
//
//   for (a=0; a < NUMPHASES; a++) {
//     filling_type_phase[a].length = ((filling_type_phase[a].volume_fraction/filling_type_phase[a].NUMCUBES)*lambda);
//   }
//
// //   fillpattern(pattern, filling_type_phase, start_cube, numperiods, gridinfo);
//
// //You need to write your filling function here
// //   fill_random();
//   fill_cube_parameters.x_start = 0;
//   fill_cube_parameters.x_end   = MESH_X-1;
//   fill_cube_parameters.y_start = 0;
//   fill_cube_parameters.y_end   = 20;
//  fill_phase_cube(fill_cube_parameters, gridinfo, 0);
//  fill_phase_cube(fill_cube_parameters, gridinfo, NUMPHASES-1);
//
//
// //   fill_phase_circle(10, 0,               0, gridinfo, 0);
// //   fill_phase_circle(10, 0.25*(MESH_X-1), 0, gridinfo, 0);
// //   fill_phase_circle(10, 0.5*(MESH_X-1),  0, gridinfo, 0);
// //   fill_phase_circle(10, 0.75*(MESH_X-1), 0, gridinfo, 0);
// //   fill_phase_circle(10, (MESH_X-1),      0, gridinfo, 0);
//
//
// //   fill_phase_circle(10, 0, 0, gridinfo, NUMPHASES-1);
//
// //   fill_phase_ellipse(50, 0.1,  -0.0, 0.5*(MESH_X-1), 0.5*(MESH_Y-1), gridinfo, 0);
// //   fill_phase_ellipse(50, 0.1,  -0.0, 0.5*(MESH_X-1), 0.5*(MESH_Y-1), gridinfo, NUMPHASES-1);
//
//   fill_composition_cube(gridinfo);
//   writetofile_struct(gridinfo, 0);
// }
//
// long assign_phase (double rand_phase) {
//   long i;
//   for (i=0; i < NUMPHASES-1; i++) {
//     if ((rand_phase >= i) && (rand_phase < (i+1)) ) {
//       return i;
//     }
//   }
// }
// int check_FLAG(int *FLAG_FILLED, long xo, long yo, long radius) {
//   long x,  y;
//   long gidy;
//   for (x =0; x < MESH_X; x++) {
//     for (y =0; y < MESH_Y; y++) {
//       gidy = x*MESH_Y + y;
//       if ((((x-xo)*(x-xo) + (y-yo)*(y-yo)) <= ((radius+2.0*epsilon)*(radius+2.0*epsilon)))) {
//         if (FLAG_FILLED[gidy] == 1) {
//           return 1;
//         }
//       }
//     }
//   }
//   return 0;
// }
// void assign_FLAG(int *FLAG_FILLED, long xo, long yo, long radius) {
//   long x,  y;
//   long gidy;
//   for (x =0; x < MESH_X; x++) {
//     for (y =0; y < MESH_Y; y++) {
//       gidy = x*MESH_Y + y;
//       if (((x-xo)*(x-xo) + (y-yo)*(y-yo)) <= ((radius+2.0*epsilon)*(radius+2.0*epsilon))) {
//        FLAG_FILLED[gidy] = 1;
//       }
//     }
//   }
// }
// void fill_random() {
//   long particlesize_bin[NUMBINS];
//   long numparticles_bin[NUMBINS];
//   int FLAG, *FLAG_FILLED;
//
//   long fill_count;
//   long x, y, b;
//   long i;
//   long gidy;
//
//   double rand_phase;
//   double rand_particle;
//
//   FLAG_FILLED = (int*)malloc((MESH_X*MESH_Y)*sizeof(*FLAG_FILLED));
//
//   for (x =0; x < MESH_X; x++) {
//     for (y =0; y < MESH_Y; y++) {
//       gidy = x*MESH_Y + y;
//       FLAG_FILLED[gidy] = 0;
//     }
//   }
//
//   for (i=0; i< NUMBINS; i++) {
//     particlesize_bin[i] = (2*i+3)*0.5*BINSIZE;
//     x                   = particlesize_bin[i]*(0.7/25.0);
//     //Normal distribution................................................
// //     numparticles_bin[i] = floor((double)NUMPARTICLES*(double)(1.0/(sqrt(2.0*M_PI*MAXPARTICLE_SIZE*0.25)))*exp(-(x-MAXPARTICLE_SIZE*0.5)*(x-MAXPARTICLE_SIZE*0.5)/(2.0*MAXPARTICLE_SIZE*0.25)));
//     //...................................................................
//     numparticles_bin[i] = floor((double)NUMPARTICLES*(0.0087+0.32101*x-0.37669*pow(x,2.0)+0.16378*pow(x,3.0)-0.02991*pow(x,4.0)+0.00181*pow(x,5.0)));
//
//     printf("%ld %ld %ld\n",i, x, numparticles_bin[i]);
//   }
//   srand(time(NULL));
//   for (i=0; i< NUMBINS; i++) {
//     fill_count=0;
//     if (numparticles_bin[i] > 0) {
//        while (fill_count <= numparticles_bin[i]) {
//           x = floor((MESH_X-1)*(double)(rand())/(double)RAND_MAX);
//           y = floor((MESH_Y-1)*(double)(rand())/(double)RAND_MAX);
//
//           gidy = x*MESH_Y + y;
//           if ((fill_count <= numparticles_bin[i]) && (x < MESH_X) && (y<MESH_Y)) {
// //               FLAG = check_FLAG(FLAG_FILLED, x, y, particlesize_bin[i]);
// //               if ((FLAG == 0)) {
//                 rand_phase = floor((NUMPHASES-1)*(double)(rand())/(double)RAND_MAX);
//                 b = assign_phase(rand_phase);
// //                 fill_phase_circle(particlesize_bin[i], x,    y,    gridinfo, b);
//                 fill_phase_ellipse(particlesize_bin[i], 0.1, b*fabs(Rtheta) + Rtheta, x,    y,    gridinfo, b);
//                 assign_FLAG(FLAG_FILLED, x, y, particlesize_bin[i]);
//                 fill_count++;
//                 printf("%ld %ld %ld, %ld\n",i , fill_count, numparticles_bin[i], x);
// //               }
//           }
//        }
//     }
//     printf("%ld %ld\n", i, fill_count);
//   }
//   free(FLAG_FILLED);
// }
//
// void find_numcubes(char *pattern, struct filling_type* filling_type_phase) {
//   long x;
//   long a;
//
//   for (a=0; a < NUMPHASES-1; a++) {
//     filling_type_phase[a].NUMCUBES = 0;
//   }
//
//   for (x=0; x < strlen(pattern); x++) {
//     for (a=0; a < NUMPHASES; a++) {
//       if ((pattern[x]) == '0'+ a) {
//         filling_type_phase[a].NUMCUBES++;
//         break;
//       }
//     }
//   }
// }
// void fillpattern(char *pattern, struct filling_type* filling_type_phase, long start_cube, long numperiods, struct variables* gridinfo) {
//   long x, a;
//   long periods =0;
//   for (periods=0; periods < numperiods; periods++) {
//     for (x=0; x < strlen(pattern); x++) {
//       for (a=0; a < NUMPHASES-1; a++) {
//         if (pattern[x] == ('0' + a)) {
//           fill_cube_parameters.x_start = start_cube;
//           fill_cube_parameters.x_end   = start_cube + filling_type_phase[a].length;
//           fill_cube_parameters.y_start = 0;
//           fill_cube_parameters.y_end   = 50;
//           fill_phase_cube(fill_cube_parameters, gridinfo, a);
//           start_cube += filling_type_phase[a].length+1;
//           break;
//         }
//       }
//     }
//   }
//   fill_phase_cube(fill_cube_parameters, gridinfo, NUMPHASES-1);
// }
// void tilting_fillpattern(char *pattern, struct filling_type* filling_type_phase, long start_cube, long numperiods, struct variables* gridinfo) {
//   long x, a;
//   long periods =0;
//   for (periods=0; periods < numperiods; periods++) {
//     for (x=0; x < strlen(pattern); x++) {
//       for (a=0; a < NUMPHASES-1; a++) {
//         if (pattern[x] == ('0' + a)) {
//           fill_cube_parameters.x_start = start_cube;
//           fill_cube_parameters.x_end   = start_cube + filling_type_phase[a].length;
//           fill_cube_parameters.y_start = 0;
//           fill_cube_parameters.y_end   = 50;
//           fill_phase_cube(fill_cube_parameters, gridinfo, a);
//           start_cube += filling_type_phase[a].length+1;
//           break;
//         }
//       }
//     }
//   }
//   fill_phase_cube(fill_cube_parameters, gridinfo, NUMPHASES-1);
// }
// void init_filling () {
//   filling_type_phase[0].volume_fraction = 0.5;
//   filling_type_phase[1].volume_fraction = 0.5;
//   filling_type_phase[2].volume_fraction = 0.5;
// //   filling_type_phase[3].volume_fraction = 0.33333333333;
//
// //   filling_type_phase[0].NUMCUBES = 2;
// //   filling_type_phase[1].NUMCUBES = 2;
// //   filling_type_phase[2].NUMCUBES = 2;
// //   filling_type_phase[3].NUMCUBES = 1;
// }
//
//
