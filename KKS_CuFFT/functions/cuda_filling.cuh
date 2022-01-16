#ifndef CUDA_FILLING_CUH_
#define CUDA_FILLING_CUH_

/* __global__ void fill_phase_cylinder_cuda
 *
 * Fills a single cylinder in the domain.
 * 3D Implementation needs to be fixed.
 *
 * Arguments:
 *  **phi                       : Phase-fraction for all phases.
 *  **dfdphi                        : Derivative of FEF wrt phase-fraction for all phases.
 *  fill_cylinder_parameters        : structure with dimensions and position of cylinder.
 *  b                               : Phase to be filled
 *  MESH_X, MESH_Y, MESH_Z          : Dimensions of the domain
 *  NUMPHASES, NUMCOMPONENTS        : Number of phases and components
 */
__global__ void fill_phase_cylinder_cuda(cufftDoubleComplex **phi,
                                         cufftDoubleComplex **dfdphi,
                                         struct fill_cylinder fill_cylinder_parameters,
                                         long b,
                                         long MESH_X, long MESH_Y, long MESH_Z,
                                         long NUMPHASES, long NUMCOMPONENTS)
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

            if (k >= fill_cylinder_parameters.z_start && k <= fill_cylinder_parameters.z_end)
            {
                if ((double)((i - fill_cylinder_parameters.x_center)*(i - fill_cylinder_parameters.x_center))
                    + (double)((j -  fill_cylinder_parameters.y_center)*(j -  fill_cylinder_parameters.y_center))
                    <= (fill_cylinder_parameters.radius*fill_cylinder_parameters.radius))
                {
                    phi[b][idx].x = 1.00000;
                    dfdphi[b][idx].x = 1.00000;
                    for (int a = 1; a < NUMPHASES; a++)
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

/* __global__ void fill_phase_cylinder_occupancy_cuda
 *
 * Fills a single cylinder in the domain, while tracking occupancy.
 * Required for random filling of the domain.
 * 3D Implementation needs to be fixed.
 *
 * Arguments:
 *  **phi                       : Phase-fraction for all phases.
 *  **dfdphi                        : Derivative of FEF wrt phase-fraction for all phases.
 *  fill_cylinder_parameters        : structure with dimensions and position of cylinder.
 *  b                               : Phase to be filled
 *  MESH_X, MESH_Y, MESH_Z          : Dimensions of the domain
 *  NUMPHASES, NUMCOMPONENTS        : Number of phases and components
 *  *occupancy                      : Matrix of dims. MESH_X*MESH_Y*MESH_Z.
 *                                    Contains flags to check if cell has been filled.
 */
__global__ void fill_phase_cylinder_occupancy_cuda(cufftDoubleComplex **phi, cufftDoubleComplex **dfdphi,
                                                   struct fill_cylinder fill_cylinder_parameters, long b,
                                                   long MESH_X, long MESH_Y, long MESH_Z,
                                                   long NUMPHASES, long NUMCOMPONENTS, int *occupancy)
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

            if (k >= fill_cylinder_parameters.z_start && k <= fill_cylinder_parameters.z_end)
            {
                if ((double)((i - fill_cylinder_parameters.x_center)*(i - fill_cylinder_parameters.x_center))
                    + (double)((j -  fill_cylinder_parameters.y_center)*(j -  fill_cylinder_parameters.y_center))
                    <= (fill_cylinder_parameters.radius))
                {
                    phi[b][idx].x = 1.00000;
                    dfdphi[b][idx].x = 1.00000;
                    occupancy[idx] = 1;
                    for (int a = 1; a < NUMPHASES; a++)
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

/* void fill_phase_cylinder_random
 *
 * Generates random coordinates and dimesions for filling cylinders randomly.
 * Fills to specified volume fraction and radius.
 *
 * Arguments:
 *  b                   : Phase to be filled.
 *  ppt_radius          : Average precipitate radius.
 *  volume_fraction     : Volume fraction threshold. May not be attained.
 *                        Change random_limit to increase number of filling trials to
 *                        increase likelihood of reaching threshold.
 *  shield_dist         : Number of cells to be maintained between particles.
 *                        Measured from cell-centre to circumference of other cells.
 */
void fill_phase_cylinder_random(long b, long ppt_radius, double volume_fraction, long shield_dist)
{
    int *overlap, *overlap_h;
    cudaMalloc((void**)&overlap, sizeof(int));
    overlap_h = (int*)malloc(sizeof(int));
    *overlap_h = 0;
    cudaMemset(overlap, 0, sizeof(int));

    int *occupancy;
    cudaMalloc((void**)&occupancy, sizeof(int)*MESH_X*MESH_Y*MESH_Z);

    double volume_domain = (double) MESH_X*MESH_Y*MESH_Z;
    double volume_per_particle = (double)(ppt_radius*ppt_radius)*PI;

    int num_particles = ceil(volume_domain*volume_fraction/volume_per_particle);
    printf("Domain volume = %lf, ppt_radius = %ld, volume_per_particle = %lf, Number of particles = %d\n",
           volume_domain, ppt_radius, volume_per_particle, num_particles);

    long particle_index = 0, random_count = 0, random_limit = 1e+5;

    double sd_rad_sq = 1.0;

    while (particle_index < num_particles)
    {
        long centx = MESH_X*ran(&SEED);
        long centy = MESH_Y*ran(&SEED);

        double mdev = 0.2*ppt_radius;
        double rad  = (double)ppt_radius + (1.0*ran(&SEED) - 0.5)*mdev;

        // Check to see if within domain bounds
        if (((centx - rad) < 0.0) || ((centx + rad) > (MESH_X-1)) ||
            ((centy - rad) < 0.0) || ((centy + rad) > (MESH_Y-1)) )
            continue;

        random_count++;
        if (random_count > random_limit)
        {
            printf("Random filling attempt count limit = %ld reached. Only %ld particles could be filled.\n",
                   random_limit, particle_index);
            break;
        }

        rad *= rad;
        sd_rad_sq = shield_dist*shield_dist*rad;
        checkOverlap_cuda_cyl<<<Gridsize, Blocksize>>>(MESH_X, MESH_Y, MESH_Z, centx, centy, sd_rad_sq, occupancy, overlap);
        cudaMemcpy(overlap_h, overlap, sizeof(int), cudaMemcpyDeviceToHost);

        // No overlap if 0
        if(*overlap_h == 0)
        {
            fill_cylinder_parameters.x_center = centx;
            fill_cylinder_parameters.y_center = centy;
            fill_cylinder_parameters.z_start  = 0;
            fill_cylinder_parameters.z_end    = 0;
            fill_cylinder_parameters.radius   = rad;

            fill_phase_cylinder_occupancy_cuda<<<Gridsize, Blocksize>>>(phiDev, dfdphiDev,
                                                                        fill_cylinder_parameters, b,
                                                                        MESH_X, MESH_Y, MESH_Z,
                                                                        NUMPHASES, NUMCOMPONENTS, occupancy);

            particle_index++;
        }
        else
        {
            cudaMemset(overlap, 0, sizeof(int));
            continue;
        }
    }
    cudaFree(occupancy);
    cudaFree(overlap);
    free(overlap_h);
}

/* __global__ void fill_phase_sphere_cuda
 *
 * Fills a single sphere in the domain, while tracking occupancy.
 *
 * Arguments:
 *  **phi                       : Phase-fraction for all phases.
 *  **dfdphi                        : Derivative of FEF wrt phase-fraction for all phases.
 *  fill_sphere_parameters        : structure with dimensions and position of sphere.
 *  b                               : Phase to be filled
 *  MESH_X, MESH_Y, MESH_Z          : Dimensions of the domain
 *  NUMPHASES, NUMCOMPONENTS        : Number of phases and components
 */
__global__ void fill_phase_sphere_cuda(cufftDoubleComplex **phi, cufftDoubleComplex **dfdphi,
                                       struct fill_sphere fill_sphere_parameters, long b,
                                       long MESH_X, long MESH_Y, long MESH_Z,
                                       long NUMPHASES, long NUMCOMPONENTS)
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
                for (int a = 1; a < NUMPHASES; a++)
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

/* __global__ void fill_phase_sphere_occupancy_cuda
 *
 * Fills a single sphere in the domain, while tracking occupancy.
 * Required for random filling of the domain.
 *
 * Arguments:
 *  **phi                       : Phase-fraction for all phases.
 *  **dfdphi                        : Derivative of FEF wrt phase-fraction for all phases.
 *  fill_sphere_parameters        : structure with dimensions and position of sphere.
 *  b                               : Phase to be filled
 *  MESH_X, MESH_Y, MESH_Z          : Dimensions of the domain
 *  NUMPHASES, NUMCOMPONENTS        : Number of phases and components
 *  *occupancy                      : Matrix of dims. MESH_X*MESH_Y*MESH_Z.
 *                                    Contains flags to check if cell has been filled.
 */
__global__ void fill_phase_sphere_occupancy_cuda(cufftDoubleComplex **phi, cufftDoubleComplex **dfdphi,
                                                 long xc, long yc, long zc, long radsq, long b,
                                                 long MESH_X, long MESH_Y, long MESH_Z,
                                                 long NUMPHASES, long NUMCOMPONENTS, int* occupancy)
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

            if ((double)((i - xc)*(i - xc)) + (double)((j - yc)*(j - yc)) + (double)((k - zc)*(k - zc))
                <= radsq)
            {
                phi[b][idx].x = 1.00000;
                dfdphi[b][idx].x = 1.00000;
                occupancy[idx] = 1;
                for (int a = 1; a < NUMPHASES; a++)
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

/* void fill_phase_sphere_random
 *
 * Generates random coordinates and dimesions for filling spheres randomly.
 * Fills to specified volume fraction and radius.
 *
 * Arguments:
 *  b                   : Phase to be filled.
 *  ppt_radius          : Average precipitate radius.
 *  volume_fraction     : Volume fraction threshold. May not be attained.
 *                        Change random_limit to increase number of filling trials to
 *                        increase likelihood of reaching threshold.
 *  shield_dist         : Number of cells to be maintained between particles.
 *                        Measured from cell-centre to circumference of other cells.
 */
void fill_phase_sphere_random(long b, long ppt_radius, double volume_fraction, long shield_dist)
{
    int *overlap, *overlap_h;
    cudaMalloc((void**)&overlap, sizeof(int));
    overlap_h = (int*)malloc(sizeof(int));
    *overlap_h = 0;
    cudaMemset(overlap, 0, sizeof(int));

    int *occupancy;
    cudaMalloc((void**)&occupancy, sizeof(int)*MESH_X*MESH_Y*MESH_Z);

    double volume_domain = (double) MESH_X*MESH_Y*MESH_Z;
    double volume_per_particle = (double)(ppt_radius*ppt_radius*ppt_radius)*(4.0/3.0)*PI;

    int num_particles = ceil(volume_domain*volume_fraction/volume_per_particle);
    printf("Domain volume = %lf, ppt_radius = %ld, volume_per_particle = %lf, Number of particles = %d\n",
           volume_domain, ppt_radius, volume_per_particle, num_particles);

    long particle_index = 0, random_count = 0, random_limit = 1e+5;

    double sd_rad_sq = 1.0;

    while (particle_index < num_particles)
    {
        long centx = MESH_X*ran(&SEED);
        long centy = MESH_Y*ran(&SEED);
        long centz = MESH_Z*ran(&SEED);

        double mdev = 0.2*ppt_radius;
        double rad  = (double)ppt_radius + (1.0*ran(&SEED) - 0.5)*mdev;

        if (((centx - rad) < 0.0) || ((centx + rad) > (MESH_X-1)) ||
            ((centy - rad) < 0.0) || ((centy + rad) > (MESH_Y-1)) ||
            ((centz - rad) < 0.0) || ((centz + rad) > (MESH_Z-1)) )
            continue;

        random_count++;
        if (random_count > random_limit)
        {
            printf("Random filling attempt count limit = %ld reached. Only %ld particles could be filled.\n",
                   random_limit, particle_index);
            break;
        }

        rad *= rad;
        sd_rad_sq = shield_dist*shield_dist*rad;
        checkOverlap_cuda_sph<<<Gridsize, Blocksize>>>(MESH_X, MESH_Y, MESH_Z, centx, centy, centz, sd_rad_sq, occupancy, overlap);
        cudaMemcpy(overlap_h, overlap, sizeof(int), cudaMemcpyDeviceToHost);

        if(*overlap_h == 0)
        {
            fill_phase_sphere_occupancy_cuda<<<Gridsize, Blocksize>>>(phiDev, dfdphiDev,
                                                                      centx, centy, centz, rad, b,
                                                                      MESH_X, MESH_Y, MESH_Z,
                                                                      NUMPHASES, NUMCOMPONENTS, occupancy);
            particle_index++;
        }
        else
        {
            cudaMemset(overlap, 0, sizeof(int));
            continue;
        }
    }
    cudaFree(occupancy);
    cudaFree(overlap);
    free(overlap_h);
}

/* void fill_composition_cuda
 *
 * Fills composition fields to match the filled phase-fraction fields.
 * Single call fills all components.
 *
 * Arguments:
 *  **comp                      : Composition fields of all components.
 *  **phi                       : Phase-fraction of all phases.
 *  MESH_X, MESH_Y, MESH_Z          : Dimensions of the domain.
 *  NUMPHASES, NUMCOMPONENTS        : Number of phases and components.
 *  c0                              : Composition in the matrix.
 *  ceq                             : Composition in the precipitate.
 */
__global__ void fill_composition_cuda(cufftDoubleComplex **comp, cufftDoubleComplex **phi,
                                      long MESH_X, long MESH_Y, long MESH_Z,
                                      long NUMPHASES, long NUMCOMPONENTS, double c0, double ceq)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long PHASE_FILLED = 0;

    long idx = k + j*MESH_Z + i*MESH_Y*MESH_Z;

    if (i < MESH_X && j < MESH_Y && k < MESH_Z)
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

    __syncthreads();
}

#endif
