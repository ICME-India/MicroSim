#ifndef FILLING_H_
#define FILLING_H_

void fill_phase_sphere(struct fill_sphere fill_sphere_parameters, struct fields* gridinfo, long b)
{
    long x, y, z, index;
    long a;
    double sum;

    long x_center, y_center, z_center;
    double radius;

    x_center = fill_sphere_parameters.x_center;
    y_center = fill_sphere_parameters.y_center;
    z_center = fill_sphere_parameters.z_center;
    radius   = fill_sphere_parameters.radius;
//
    for (x = 0; x < rows_x; x++)
    {
        for (y = 0; y < rows_y; y++)
        {
            for (z = 0; z < rows_z; z++)
            {
                index = x*layer_size + y*rows_z + z;
                if (b > 0)
                {
                    if (((x-x_center)*(x-x_center) + (y-y_center)*(y-y_center) + (z-z_center)*(z-z_center) <= radius*radius))
                    {
                        gridinfo[index].phia[b].x = 1.00000;
                        for (a = 0; a < NUMPHASES; a++)
                        {
                            if (b != a)
                            {
                                gridinfo[index].phia[a].x = 0.00000;
                            }
                        }
                    }
                    else if (gridinfo[index].phia[b].x != 1.0000)
                    {
                        gridinfo[index].phia[b].x = 0.00000;
                    }
                }
                else
                {
                    sum = 0.0;
                    for (a = 1; a < NUMPHASES; a++)
                        sum += gridinfo[index].phia[a].x;
                    if (sum > 1.0)
                    {
                        printf("Wrong filling operation, will fill it with matrix\n");
                        gridinfo[index].phia[b].x = 0.0;
                        for (a = 1; a < NUMPHASES; a++)
                        {
                            if (a != b)
                            {
                                gridinfo[index].phia[a].x = 0.00000;
                            }
                        }
                    }
                    else
                        gridinfo[index].phia[b].x  = 1.0 - sum;
                }
            }
        }
    }
}

void fill_phase_sphere_occupancy(struct fill_sphere fill_sphere_parameters, struct fields* gridinfo, long b, int* occupancy)
{
    long x, y, z, index;
    long a;
    double sum;

    long x_center, y_center, z_center;
    double radius;

    x_center = fill_sphere_parameters.x_center;
    y_center = fill_sphere_parameters.y_center;
    z_center = fill_sphere_parameters.z_center;
    radius   = fill_sphere_parameters.radius;

    for (x = 0; x < rows_x; x++)
    {
        for (y = 0; y < rows_y; y++)
        {
            for (z = 0; z < rows_z; z++)
            {
                index = x*layer_size + y*rows_z + z;
                if (b > 0)
                {
                    if (((x-x_center)*(x-x_center) + (y-y_center)*(y-y_center) + (z- z_center)*(z-z_center) <= radius*radius))
                    {
                        gridinfo[index].phia[b].x = 1.00000;
                        occupancy[index] = 1;
                        for (a = 0; a < NUMPHASES; a++)
                        {
                            if (b != a)
                            {
                                gridinfo[index].phia[a].x = 0.00000;
                            }
                        }
                    }
                    else if (gridinfo[index].phia[b].x != 1.0000)
                    {
                        gridinfo[index].phia[b].x = 0.00000;
                    }
                }
                else
                {
                    sum = 0.0;
                    for (a = 1; a < NUMPHASES; a++)
                        sum += gridinfo[index].phia[a].x;
                    if (sum > 1.0)
                    {
                        printf("Wrong filling operation, will fill it with matrix\n");
                        gridinfo[index].phia[b].x = 0.0;
                        for (a = 1; a < NUMPHASES; a++)
                        {
                            if (a != b)
                            {
                                gridinfo[index].phia[a].x = 0.00000;
                            }
                        }
                    }
                    else
                        gridinfo[index].phia[b].x   = 1.0 - sum;
                }
            }
        }
    }
}

void fill_phase_cylinder(struct fill_cylinder fill_cylinder_parameters, struct fields* gridinfo, long b)
{
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

    for(x = 0; x < rows_x; x++)
    {
        for (z = 0; z < rows_z; z++)
        {
            for (y = 0; y < rows_y; y++)
            {
                index = x*layer_size + z*rows_y + y;
                if (b > 0)
                {
                    if ((z >= z_start) && (z <= z_end))
                    {
                        if (((x-x_center)*(x-x_center) + (y-y_center)*(y-y_center)) <= (radius*radius))
                        {
                            gridinfo[index].phia[b].x = 1.00000;
                            for (a = 0; a < NUMPHASES; a++)
                            {
                                if (b != a)
                                {
                                    gridinfo[index].phia[a].x = 0.00000;
                                }
                            }
                        }
                        else if (gridinfo[index].phia[b].x != 1.0000)
                            gridinfo[index].phia[b].x = 0.00000;
                    }
                }
                else
                {
                    sum = 0.0;
                    for (a = 1; a < NUMPHASES; a++)
                        sum += gridinfo[index].phia[a].x;

                    if (sum > 1.0)
                    {
                        printf("Wrong filling operation, will fill it with matrix\n");
                        gridinfo[index].phia[b].x = 0.0;
                        for (a = 1; a < NUMPHASES; a++)
                            if (a != b)
                                gridinfo[index].phia[a].x = 0.00000;
                    }
                    else
                        gridinfo[index].phia[b].x = 1.0 - sum;
                }
            }
        }
    }
}

void fill_phase_cylinder_occupancy(struct fill_cylinder fill_cylinder_parameters, struct fields* gridinfo, long b, int* occupancy)
{
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

    for(x = 0; x < rows_x; x++)
    {
        for (y = 0; y < rows_y; y++)
        {
            for (z = 0; z < rows_z; z++)
            {
                index = x*layer_size + y*rows_z + z;
                if (b > 0)
                {
                    if ((z >= z_start) && (z <= z_end))
                    {
                        if (((x-x_center)*(x-x_center) + (y-y_center)*(y-y_center) <= radius*radius))
                        {
                            gridinfo[index].phia[b].x = 1.00000;
                            occupancy[index] = 1;
                            for (a = 0; a < NUMPHASES; a++)
                            {
                                if (b != a)
                                {
                                    gridinfo[index].phia[a].x = 0.00000;
                                }
                            }
                        }
                        else if (gridinfo[index].phia[b].x != 1.0000)
                            gridinfo[index].phia[b].x = 0.00000;
                    }
                }
                else
                {
                    sum = 0.0;
                    for (a = 1; a < NUMPHASES; a++)
                        sum += gridinfo[index].phia[a].x;

                    if (sum > 1.0)
                    {
                        printf("Wrong filling operation, will fill it with matrix\n");
                        gridinfo[index].phia[b].x = 0.0;
                        for (a = 1; a < NUMPHASES; a++)
                            if (a != b)
                                gridinfo[index].phia[a].x = 0.00000;
                    }
                    else
                        gridinfo[index].phia[b].x = 1.0 - sum;
                }
            }
        }
    }
}

// void fill_phase_cylinder_random(struct fields* gridinfo, long b)
// {
//     int overlap = 1;
//     int *occupancy;
//     occupancy = (int*)malloc(sizeof(int)*MESH_X*MESH_Y*MESH_Z);
//
//     for (int i = 0; i < MESH_X*MESH_Y*MESH_Z; i++)
//         occupancy[i] = 0;
//
//     double volume_domain = (double) MESH_X*MESH_Y*MESH_Z;
//     double volume_per_particle = (double)(ppt_radius*ppt_radius)*PI;
//
//     int num_particles = ceil(volume_domain*volume_fraction/volume_per_particle);
//     printf("Domain volume = %lf, ppt_radius = %ld, volume_per_particle = %lf, Number of particles = %d\n", volume_domain, ppt_radius, volume_per_particle, num_particles);
//
//     long particle_index = 1, random_count = 0, random_limit = 1e+4;
//
//     struct fill_cylinder temp_cyl;
//     temp_cyl.z_start  = 0;
//     temp_cyl.z_end    = 0;
//
//     while (particle_index <= num_particles)
//     {
//         long centx = MESH_X*ran(&SEED);
//         long centy = MESH_Y*ran(&SEED);
//         long centz = 0;
//
//         temp_cyl.x_center = centx;
//         temp_cyl.y_center = centy;
//
//         double mdev = 0.5*ppt_radius;
//         double rad  = (double)ppt_radius + (1.0*ran(&SEED) - 0.5)*mdev;
//
//         temp_cyl.radius = rad;
//         printf("particle_index = %d, x_center = %ld, y_center = %ld, z_start = %ld, z_end = %ld, radius = %lf\n", particle_index, centx, centy, temp_cyl.z_start, temp_cyl.z_end, rad);
//
//         random_count++;
//         if (random_count > random_limit)
//         {
//             printf("Random filling attempt count limit = %ld reached. Only %ld particles could be filled.\n", random_limit, particle_index-1);
//             break;
//         }
//
//         if (((centx - rad) < 0.0) || ((centx + rad) > (MESH_X-1)) ||
//             ((centy - rad) < 0.0) || ((centy + rad) > (MESH_Y-1)) )
//             continue;
//
//         overlap = checkOverlap(MESH_X, MESH_Y, MESH_Z, centx, centy, centz, rad, occupancy);
//
//         if(overlap == 0)
//         {
//             fill_phase_cylinder_occupancy(temp_cyl, gridinfo, b, occupancy);
//             fill_phase_cylinder_occupancy(temp_cyl, gridinfo, 0, occupancy);
//             fill_phase_cylinder(temp_cyl, gridinfo, b);
//         }
//         else
//         {
//             overlap = 1;
//             continue;
//         }
//
//         particle_index++;
//     }
//     free(occupancy);
// }

void fill_phase_sphere_random(struct fields* gridinfo, long b)
{
    int overlap = 1;
    int *occupancy;
    occupancy = (int*)malloc(sizeof(int)*MESH_X*MESH_Y*MESH_Z);

    for (int i = 0; i < MESH_X*MESH_Y*MESH_Z; i++)
        occupancy[i] = 0;

    double volume_domain = (double) MESH_X*MESH_Y*MESH_Z;
    double volume_per_particle = (double)(ppt_radius*ppt_radius*ppt_radius)*(4.0/3.0)*PI;

    int num_particles = ceil(volume_domain*volume_fraction/volume_per_particle);
    printf("Domain volume = %lf, ppt_radius = %ld, volume_per_particle = %lf, Number of particles = %d\n", volume_domain, ppt_radius, volume_per_particle, num_particles);

    long particle_index = 1, random_count = 0, random_limit = 1e+5;

    struct fill_sphere temp_sph;

    while (particle_index <= num_particles)
    {
        long centx = MESH_X*ran(&SEED);
        long centy = MESH_Y*ran(&SEED);
        long centz = MESH_Z*ran(&SEED);

        temp_sph.x_center = centx;
        temp_sph.y_center = centy;
        temp_sph.z_center = centz;

        double mdev = 0.5*ppt_radius;
        double rad  = (double)ppt_radius + (1.0*ran(&SEED) - 0.5)*mdev;

        temp_sph.radius = rad;
        //printf("particle_index = %d, x_center = %ld, y_center = %ld, z_start = %ld, z_end = %ld, radius = %lf\n", particle_index, centx, centy, temp_cyl.z_start, temp_cyl.z_end, rad);

        random_count++;
        if (random_count > random_limit)
        {
            printf("Random filling attempt count limit = %ld reached. Only %ld particles could be filled.\n", random_limit, particle_index-1);
            break;
        }

        if (((centx - rad) < 0.0) || ((centx + rad) > (MESH_X-1)) ||
            ((centy - rad) < 0.0) || ((centy + rad) > (MESH_Y-1)) ||
            ((centz - rad) < 0.0) || ((centz + rad) > (MESH_Z-1)) )
            continue;

        overlap = checkOverlap(MESH_X, MESH_Y, MESH_Z, centx, centy, centz, rad, occupancy);

        if(overlap == 0)
        {
            fill_phase_sphere_occupancy(temp_sph, gridinfo, b, occupancy);
            //fill_phase_cylinder_occupancy(temp_cyl, gridinfo, 0, occupancy);
            //fill_phase_cylinder(temp_cyl, gridinfo, b);
        }
        else
        {
            overlap = 1;
            continue;
        }

        particle_index++;
    }
    free(occupancy);
}

void fill_composition_cube(struct fields* gridinfo)
{
    long x, y, z, index;
    long k;
    long b;
    long PHASE_FILLED=0;

    for(x = 0; x < rows_x; x++)
    {
        for(y = 0; y < rows_y; y++)
        {
            for (z = 0; z < rows_z; z++)
            {
                index = x*layer_size + y*rows_z + z;
                PHASE_FILLED = 0;
                for (b = 1; b < NUMPHASES; b++)
                {
                    if (gridinfo[index].phia[b].x == 1.0)
                    {
                        for (k = 0; k < NUMCOMPONENTS-1; k++)
                            gridinfo[index].comp[k].x = ceq[b][b][k];
                        PHASE_FILLED = 1;
                        break;
                    }
                }
                if (!PHASE_FILLED)
                    for (k = 0; k < NUMCOMPONENTS-1; k++)
                        gridinfo[index].comp[k].x = c0;
            }
        }
    }
}
#endif
