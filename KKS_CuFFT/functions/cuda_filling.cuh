__global__ void fill_phase_cylinder_cuda(struct fill_cylinder fill_cylinder_parameters,
                                         cufftDoubleComplex **phi, cufftDoubleComplex **dfdphi, long b, long MESH_X, long MESH_Y,
                                         long MESH_Z, long NUMPHASES, long NUMCOMPONENTS)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    double sum;

    int idx = k + j*MESH_Z + i*MESH_Y*MESH_Z;

    if (i < MESH_X && j < MESH_Y && k < MESH_Z)
    {
        if (b > 0)
        {
            phi[b][idx].y = 0.0;
            dfdphi[b][idx].y = 0.0;

            if (k >= fill_cylinder_parameters.z_start && k <= fill_cylinder_parameters.z_end)
            {
                if ((double)(i - fill_cylinder_parameters.x_center)*(double)(i - fill_cylinder_parameters.x_center) +
                    (double)(j - fill_cylinder_parameters.y_center)*(double)(j - fill_cylinder_parameters.y_center) <=
                    (fill_cylinder_parameters.radius*fill_cylinder_parameters.radius))
                {
                    phi[b][idx].x = 1.00000;
                    dfdphi[b][idx].x = 1.00000;
                    for (int a = 0; a < NUMPHASES; a++)
                        if (a != b)
                        {
                            phi[a][idx].x = 0.00000;
                            dfdphi[a][idx].x = 0.00000;
                        }
                }
                else if(phi[b][idx].x != 1.00000)
                {
                    phi[b][idx].x = 0.00000;
                    dfdphi[b][idx].x = 0.00000;
                }
            }
        }
        else
        {
            sum = 0.0;
            for (int a = 1; a < NUMPHASES; a++)
                sum += phi[b][idx].x;
            if (sum > 1.0)
            {
                phi[b][idx].x = 0.0;
                dfdphi[b][idx].x = 0.0;
                for (int a = 1; a < NUMPHASES; a++)
                {
                    phi[b][idx].x = 0.0;
                    dfdphi[b][idx].x = 0.0;
                }
            }
            else
            {
                phi[b][idx].x = 1.0 - sum;
                dfdphi[b][idx].x = 1.0 - sum;
            }
        }
    }
    __syncthreads();
}

__global__ void fill_phase_cylinder_occupancy_cuda(struct fill_cylinder fill_cylinder_parameters,
                                                   cufftDoubleComplex **phi, cufftDoubleComplex **dfdphi, long b, long MESH_X, long MESH_Y,
                                                   long MESH_Z, long NUMPHASES, long NUMCOMPONENTS, int *occupancy)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    double sum;

    int idx = k + j*MESH_Z + i*MESH_Y*MESH_Z;

    if (i < MESH_X && j < MESH_Y && k < MESH_Z)
    {
        if (b > 0)
        {
            phi[b][idx].y = 0.0;
            dfdphi[b][idx].y = 0.0;

            if (k >= fill_cylinder_parameters.z_start && k <= fill_cylinder_parameters.z_end)
            {
                if ((double)(i - fill_cylinder_parameters.x_center)*(double)(i - fill_cylinder_parameters.x_center) +
                    (double)(j - fill_cylinder_parameters.y_center)*(double)(j - fill_cylinder_parameters.y_center) <=
                    (fill_cylinder_parameters.radius*fill_cylinder_parameters.radius))
                {
                    phi[b][idx].x = 1.00000;
                    dfdphi[b][idx].x = 1.00000;
                    occupancy[idx] = 1;
                    for (int a = 0; a < NUMPHASES; a++)
                        if (a != b)
                        {
                            phi[a][idx].x = 0.00000;
                            dfdphi[a][idx].x = 0.00000;
                        }
                }
                else if(phi[b][idx].x != 1.00000)
                {
                    phi[b][idx].x = 0.00000;
                    dfdphi[b][idx].x = 0.00000;
                }
            }
        }
        else
        {
            sum = 0.0;
            for (int a = 1; a < NUMPHASES; a++)
                sum += phi[b][idx].x;
            if (sum > 1.0)
            {
                phi[b][idx].x = 0.0;
                dfdphi[b][idx].x = 0.0;
                for (int a = 1; a < NUMPHASES; a++)
                {
                    phi[b][idx].x = 0.0;
                    dfdphi[b][idx].x = 0.0;
                }
            }
            else
            {
                phi[b][idx].x = 1.0 - sum;
                dfdphi[b][idx].x = 1.0 - sum;
            }
        }
    }
    __syncthreads();
}

void fill_phase_cylinder_random(long b)
{
    int *overlap, overlap_h = 0;
    cudaMalloc((void**)&overlap, sizeof(int));
    cudaMemset(overlap, 0, sizeof(int));
    int *occupancy;
    cudaMalloc((void**)&occupancy, sizeof(int)*MESH_X*MESH_Y*MESH_Z);
    cudaMemset(occupancy, 0, MESH_X*MESH_Y*MESH_Z*sizeof(int));

    double volume_domain = (double) MESH_X*MESH_Y*MESH_Z;
    double volume_per_particle = (double)(ppt_radius*ppt_radius)*PI;

    int num_particles = ceil(volume_domain*volume_fraction/volume_per_particle);
    printf("Domain volume = %lf, ppt_radius = %ld, volume_per_particle = %lf, Number of particles = %d\n", volume_domain, ppt_radius, volume_per_particle, num_particles);

    long particle_index = 1, random_count = 0, random_limit = 1e+5;

    struct fill_cylinder temp_cyl;
    temp_cyl.z_start  = 0;
    temp_cyl.z_end    = 0;

    while (particle_index <= num_particles)
    {
        long centx = MESH_X*ran(&SEED);
        long centy = MESH_Y*ran(&SEED);

        temp_cyl.x_center = centx;
        temp_cyl.y_center = centy;

        double mdev = 0.2*ppt_radius;
        double rad  = (double)ppt_radius + (1.0*ran(&SEED) - 0.5)*mdev;

        temp_cyl.radius = rad;
        //printf("particle_index = %d, x_center = %ld, y_center = %ld, z_start = %ld, z_end = %ld, radius = %lf\n", particle_index, centx, centy, temp_cyl.z_start, temp_cyl.z_end, rad);

        random_count++;
        if (random_count > random_limit)
        {
            printf("Random filling attempt count limit = %ld reached. Only %ld particles could be filled.\n", random_limit, particle_index-1);
            break;
        }

        if (((centx - rad) < 0.0) || ((centx + rad) > (MESH_X-1)) ||
            ((centy - rad) < 0.0) || ((centy + rad) > (MESH_Y-1)) )
            continue;

        checkOverlap_cuda<<<Gridsize, Blocksize>>>(MESH_X, MESH_Y, MESH_Z, centx, centy, 0, rad, occupancy, overlap, shield_dist);
        cudaMemcpy(&overlap_h, overlap, sizeof(int), cudaMemcpyDeviceToHost);

        if(overlap_h == 0)
        {
            fill_phase_cylinder_occupancy_cuda<<<Gridsize, Blocksize>>>(temp_cyl, PhiBuff, dfdphi, b, MESH_X, MESH_Y, MESH_Z, NUMPHASES, NUMCOMPONENTS, occupancy);
        }
        else
        {
            cudaMemset(overlap, 0, sizeof(int));
            continue;
        }
        particle_index++;
    }
    cudaFree(occupancy);
    cudaFree(overlap);
}

__global__ void fill_phase_sphere_cuda(struct fill_sphere fill_sphere_parameters,
                                       cufftDoubleComplex **phi, cufftDoubleComplex **dfdphi, long b, long MESH_X, long MESH_Y,
                                       long MESH_Z, long NUMPHASES, long NUMCOMPONENTS)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    double sum;

    long idx = k + j*MESH_Z + i*MESH_Y*MESH_Z;

    if (i < MESH_X && j < MESH_Y && k < MESH_Z)
    {
        if (b > 0)
        {
            phi[b][idx].y = 0.0;
            dfdphi[b][idx].y = 0.0;

            if ((double)(i - fill_sphere_parameters.x_center)*(double)(i - fill_sphere_parameters.x_center) +
                (double)(j -  fill_sphere_parameters.y_center)*(double)(j -  fill_sphere_parameters.y_center) +
                (double)(k -  fill_sphere_parameters.z_center)*(double)(k -  fill_sphere_parameters.z_center) <=
                (fill_sphere_parameters.radius*fill_sphere_parameters.radius))
            {
                phi[b][idx].x = 1.00000;
                dfdphi[b][idx].x = 1.00000;
                for (int a = 0; a < NUMPHASES; a++)
                    if (a != b)
                    {
                        phi[a][idx].x = 0.00000;
                        dfdphi[a][idx].x = 0.00000;
                    }
            }
            else if(phi[b][idx].x != 1.00000)
            {
                phi[b][idx].x = 0.00000;
                dfdphi[b][idx].x = 0.00000;
            }
        }
        else
        {
            sum = 0.0;
            for (int a = 1; a < NUMPHASES; a++)
                sum += phi[b][idx].x;
            if (sum > 1.0)
            {
                phi[b][idx].x = 0.0;
                dfdphi[b][idx].x = 0.0;
                for (int a = 1; a < NUMPHASES; a++)
                {
                    phi[b][idx].x = 0.0;
                    dfdphi[b][idx].x = 0.0;
                }
            }
            else
            {
                phi[b][idx].x = 1.0 - sum;
                dfdphi[b][idx].x = 1.0 - sum;
            }
        }
    }
    __syncthreads();
}

__global__ void fill_phase_sphere_occupancy_cuda(struct fill_sphere fill_sphere_parameters,
                                       cufftDoubleComplex **phi, cufftDoubleComplex **dfdphi, long b, long MESH_X, long MESH_Y,
                                       long MESH_Z, long NUMPHASES, long NUMCOMPONENTS, int* occupancy)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    double sum;

    long idx = k + j*MESH_Z + i*MESH_Y*MESH_Z;

    if (i < MESH_X && j < MESH_Y && k < MESH_Z)
    {
        if (b > 0)
        {
            phi[b][idx].y = 0.0;
            dfdphi[b][idx].y = 0.0;

            if ((double)(i - fill_sphere_parameters.x_center)*(double)(i - fill_sphere_parameters.x_center) +
                (double)(j -  fill_sphere_parameters.y_center)*(double)(j -  fill_sphere_parameters.y_center) +
                (double)(k -  fill_sphere_parameters.z_center)*(double)(k -  fill_sphere_parameters.z_center) <=
                (fill_sphere_parameters.radius*fill_sphere_parameters.radius))
            {
                phi[b][idx].x = 1.00000;
                dfdphi[b][idx].x = 1.00000;
                occupancy[idx] = 1;
                for (int a = 0; a < NUMPHASES; a++)
                    if (a != b)
                    {
                        phi[a][idx].x = 0.00000;
                        dfdphi[a][idx].x = 0.00000;
                    }
            }
            else if(phi[b][idx].x != 1.00000)
            {
                phi[b][idx].x = 0.00000;
                dfdphi[b][idx].x = 0.00000;
            }
        }
        else
        {
            sum = 0.0;
            for (int a = 1; a < NUMPHASES; a++)
                sum += phi[b][idx].x;
            if (sum > 1.0)
            {
                phi[b][idx].x = 0.0;
                dfdphi[b][idx].x = 0.0;
                for (int a = 1; a < NUMPHASES; a++)
                {
                    phi[b][idx].x = 0.0;
                    dfdphi[b][idx].x = 0.0;
                }
            }
            else
            {
                phi[b][idx].x = 1.0 - sum;
                dfdphi[b][idx].x = 1.0 - sum;
            }
        }
    }
    __syncthreads();
}

void fill_phase_sphere_random(long b)
{
    int *overlap, overlap_h = 0;
    cudaMalloc((void**)&overlap, sizeof(int));
    cudaMemset(overlap, 0, sizeof(int));
    int *occupancy;
    cudaMalloc((void**)&occupancy, sizeof(int)*MESH_X*MESH_Y*MESH_Z);
    cudaMemset(occupancy, 0, MESH_X*MESH_Y*MESH_Z*sizeof(int));

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

        double mdev = 0.1*ppt_radius;
        double rad  = (double)ppt_radius + (1.0*ran(&SEED) - 0.5)*mdev;

        temp_sph.radius = rad;
        //printf("particle_index = %d, x_center = %ld, y_center = %ld, z_center = %ld, radius = %lf\n", particle_index, centx, centy, centz, rad);

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

        checkOverlap_cuda<<<Gridsize, Blocksize>>>(MESH_X, MESH_Y, MESH_Z, centx, centy, centz, rad, occupancy, overlap, shield_dist);
        cudaMemcpy(&overlap_h, overlap, sizeof(int), cudaMemcpyDeviceToHost);

        if(overlap_h == 0)
        {
            fill_phase_sphere_occupancy_cuda<<<Gridsize, Blocksize>>>(temp_sph, PhiBuff, dfdphi, b, MESH_X, MESH_Y, MESH_Z, NUMPHASES, NUMCOMPONENTS, occupancy);
        }
        else
        {
            cudaMemset(overlap, 0, sizeof(int));
            continue;
        }
        particle_index++;
    }
    cudaFree(occupancy);
    cudaFree(overlap);
}

__global__ void fill_composition_cuda(cufftDoubleComplex **phi, cufftDoubleComplex **comp,
                                      long MESH_X, long MESH_Y, long MESH_Z,
                                      long NUMPHASES, long NUMCOMPONENTS, double c0, double ceq)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long PHASE_FILLED = 0;

    long idx = k + j*MESH_Z + i*MESH_Y*MESH_Z;

    if (i >= 0 && i < MESH_X && j >= 0 && j < MESH_Y && k >= 0 && k < MESH_Z)
    {
        PHASE_FILLED = 0;

        for (int a = 1; a < NUMPHASES; a++)
        {
            if (phi[a][idx].x == 1.0)
            {
                for (int b = 0; b < NUMCOMPONENTS-1; b++)
                {
                    comp[b][idx].x = ceq;
                    comp[b][idx].y = 0.0;
                }
                PHASE_FILLED = 1;
                break;
            }
        }
        if (!PHASE_FILLED)
            for (int b = 0; b < NUMCOMPONENTS-1; b++)
            {
                comp[b][idx].x = c0;
                comp[b][idx].y = 0.0;
            }
    }
}
