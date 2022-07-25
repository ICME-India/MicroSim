#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <endian.h>

#define  IS_BIG_ENDIAN     (1 == htons(1))
#define  IS_LITTLE_ENDIAN  (!IS_BIG_ENDIAN)

int DIMENSION;

long temp_mx = 0, temp_my = 0, temp_mz = 0;

long MESH_X, MESH_Y, MESH_Z;
double DELTA_X, DELTA_Y, DELTA_Z;

int numProcs;

long NUMPHASES, NUMCOMPONENTS;

long start, end, step;

double DELTA_t;

char where[1000];
char name[1000];
char fileName[1000];

char outputName[1000];

char writeformat[10];

double **comp;
double **phi;

double timeOld = 0.0, radiusOld = 0.0;

double sumInner = 0.0, sumOuter = 0.0, minC = 1.5, maxC = -1.5;

double swap_bytes(double value)
{
    double  src_num = value;
    int64_t tmp_num = htobe64(le64toh(*(int64_t*)&src_num));
    double  dst_num = *(double*)&tmp_num;
    return  dst_num;
}

void reading_input_parameters(char *argv[])
{
    FILE *fr;
    if (fr = fopen(argv[1], "rt"))
        printf("\nReading information from %s\n", argv[1]);
    else
        printf("\nFile %s not found\n", argv[1]);

    int i;
    char tempbuff[1000];
    char tmpstr1[100];
    char tmpstr2[100];

    while (fgets(tempbuff, 1000, fr))
    {
        sscanf(tempbuff, "%100s = %100[^;];", tmpstr1, tmpstr2);

        if (tmpstr1[0] != '#')
        {
            if (strcmp(tmpstr1, "DIMENSION") == 0)
            {
                DIMENSION = atoi(tmpstr2);
            }
            else if (strcmp(tmpstr1, "MESH_X") == 0)
            {
                MESH_X = atol(tmpstr2);
            }
            else if (strcmp(tmpstr1, "MESH_Y") == 0)
            {
                MESH_Y = atol(tmpstr2);
            }
            else if (strcmp(tmpstr1, "MESH_Z") == 0)
            {
                MESH_Z = atol(tmpstr2);
            }
            else if (strcmp(tmpstr1, "DELTA_X") == 0)
            {
                DELTA_X = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1, "DELTA_Y") == 0)
            {
                DELTA_Y = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1, "DELTA_Z") == 0)
            {
                DELTA_Z = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1, "DELTA_t") == 0)
            {
                DELTA_t = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1, "numworkers") == 0)
            {
                numProcs = atoi(tmpstr2);
            }
            else if (strcmp(tmpstr1, "NUMPHASES") == 0)
            {
                NUMPHASES = atol(tmpstr2);
            }
            else if (strcmp(tmpstr1, "NUMCOMPONENTS") == 0)
            {
                NUMCOMPONENTS = atol(tmpstr2);
            }
            else if (strcmp(tmpstr1, "WRITEFORMAT") == 0)
            {
                strcpy(writeformat, tmpstr2);
                for (int i = 0; i < strlen(writeformat); i++)
                    if (writeformat[i] == '\n')
                        writeformat[i] = '\0';
            }
            else if (strcmp(tmpstr1, "START") == 0)
            {
                start = atol(tmpstr2);
            }
            else if (strcmp(tmpstr1, "INTERVAL") == 0)
            {
                step = atol(tmpstr2);
            }
            else if (strcmp(tmpstr1, "END") == 0)
            {
                end = atol(tmpstr2);
            }
            else if (strcmp(tmpstr1, "LOCATION") == 0)
            {
                strcpy(where, tmpstr2);
                for (int i = 0; i < strlen(where); i++)
                    if (where[i] == '\n')
                        where[i] = '\0';            }
                        else if (strcmp(tmpstr1, "FILENAME") == 0)
                        {
                            strcpy(name, tmpstr2);
                            for (int i = 0; i < strlen(name); i++)
                                if (name[i] == '\n')
                                    name[i] = '\0';
                        }
                        else if (strcmp(tmpstr1, "OUTPUTNAME") == 0)
                        {
                            strcpy(outputName, tmpstr2);
                            for (int i = 0; i < strlen(outputName); i++)
                                if (outputName[i] == '\n')
                                    outputName[i] = '\0';
                        }
                        else
                        {
                            printf("Unrecognized parameter : \"%s\"\n", tmpstr1);
                        }
        }
    }

    fclose(fr);
}

int read_vtk_ascii(FILE *fp)
{
    char temp[1000];
    long a = 0, k = 0, x, y, z;
    long index;

    while(fscanf(fp, "%s", temp))
    {
        if (strcmp(temp, "DIMENSIONS") == 0)
        {
            fscanf(fp, "%ld", &temp_mx);
            fscanf(fp, "%ld", &temp_my);
            fscanf(fp, "%ld", &temp_mz);

            printf("Read dimensions: %ld, %ld, %ld\n", temp_mx, temp_my, temp_mz);
            printf("Input dimensions: %ld, %ld, %ld\n", MESH_X, MESH_Y, MESH_Z);

            if (temp_mx != MESH_X || temp_my != MESH_Y || temp_mz != MESH_Z)
            {
                printf("Dimensions do not match. Filling using specified filling file\n");
                //return 0;
            }
        }
        else if (strcmp(temp, "SPACING") == 0)
        {
            double temp_dx = 0, temp_dy = 0, temp_dz = 0;

            fscanf(fp, "%le", &temp_dx);
            fscanf(fp, "%le", &temp_dy);
            fscanf(fp, "%le", &temp_dz);

            if (temp_dx != DELTA_X || temp_dy != DELTA_Y || temp_dz != DELTA_Z)
            {
                printf("Dimensions do not match. Filling using specified filling file\n");
                return 0;
            }
        }
        else if (strcmp(temp, "default") == 0 && a < NUMPHASES)
        {
            printf("Reading phase %d\n", a);
            for (x = 0; x < temp_mx; x++) {
                for (y = 0; y < temp_my; y++) {
                    for (z = 0; z < temp_mz; z++) {
                        index = (y + x*temp_my)*temp_mz + z;
                        fscanf(fp, "%le", &phi[a][index]);
                    }
                }
            }
            a++;
        }
        else if (strcmp(temp, "default") == 0 && k < NUMCOMPONENTS-1 && a >= NUMPHASES)
        {
            printf("Reading component %d\n", k);
            for (x = 0; x < temp_mx; x++) {
                for (y = 0; y < temp_my; y++) {
                    for (z = 0; z < temp_mz; z++) {
                        index = (y + x*temp_my)*temp_mz + z;
                        fscanf(fp, "%le", &comp[k][index]);
                    }
                }
            }
            k++;
        }

        else if (k == NUMCOMPONENTS-1)
            break;
    }
    return 1;
}

int read_vtk_binary(FILE *fp)
{
    char temp[1000];
    long a = 0, k = 0, x, y, z;
    long index;
    double value;

    while(fscanf(fp, "%s", temp))
    {
        if (strcmp(temp, "DIMENSIONS") == 0)
        {
            fscanf(fp, "%ld", &temp_mx);
            fscanf(fp, "%ld", &temp_my);
            fscanf(fp, "%ld", &temp_mz);

            printf("Read dimensions: %ld, %ld, %ld\n", temp_mx, temp_my, temp_mz);
            printf("Input dimensions: %ld, %ld, %ld\n", MESH_X, MESH_Y, MESH_Z);

            if (temp_mx != MESH_X || temp_my != MESH_Y || temp_mz != MESH_Z)
            {
                printf("Dimensions do not match. Filling using specified filling file\n");
                //return 0;
            }
        }
        else if (strcmp(temp, "SPACING") == 0)
        {
            double temp_dx = 0, temp_dy = 0, temp_dz = 0;

            fscanf(fp, "%le", &temp_dx);
            fscanf(fp, "%le", &temp_dy);
            fscanf(fp, "%le", &temp_dz);

            if (temp_dx != DELTA_X || temp_dy != DELTA_Y || temp_dz != DELTA_Z)
            {
                printf("Dimensions do not match. Filling using specified filling file\n");
                return 0;
            }
        }
        else if (strcmp(temp, "default") == 0 && k < NUMCOMPONENTS-1)
        {
            printf("Reading component %d\n", k);
            fscanf(fp, "%lf", &value);

            for (x = 0; x < temp_mx; x++) {
                for (y = 0; y < temp_my; y++) {
                    for (z = 0; z < temp_mz; z++) {
                        index = (y + x*temp_my)*temp_mz + z;
                        fread(&value, sizeof(double), 1, fp);
                        if (IS_LITTLE_ENDIAN)
                            comp[k][index] = swap_bytes(value);
                        else
                            comp[k][index] = value;
                    }
                }
            }
            k++;
        }
        else if (strcmp(temp, "default") == 0 && a < NUMPHASES && k >= NUMCOMPONENTS-1)
        {
            printf("Reading phase %d\n", a);
            fscanf(fp, "%lf", &value);

            for (x = 0; x < temp_mx; x++) {
                for (y = 0; y < temp_my; y++) {
                    for (z = 0; z < temp_mz; z++) {
                        index = (y + x*temp_my)*temp_mz + z;
                        fread(&value, sizeof(double), 1, fp);
                        if (IS_LITTLE_ENDIAN)
                            phi[a][index] = swap_bytes(value);
                        else
                            phi[a][index] = value;
                    }
                }
            }
            a++;
        }
        else if (a == NUMPHASES)
            break;
    }
    return 1;
}

void calcStats(long count)
{
    long index;
    for (long x = 0; x < temp_mx; x++)
    {
        for (long y = 0; y < temp_my; y++)
        {
            for (long z = 0; z < temp_mz; z++)
            {
                index = ((y + x*temp_my)*temp_mz + z);

                if (phi[1][index] >= 0.05)
                    sumInner++;
                if (phi[1][index] >= 0.95)
                    sumOuter++;

                if (comp[0][index] < minC)
                    minC = comp[0][index];
                if (comp[0][index] > maxC)
                    maxC = comp[0][index];
            }
        }
    }
}

int main(int argc, char *argv[])
{
    reading_input_parameters(argv);

    int found = 0;

    comp = malloc(sizeof(double*)*(NUMCOMPONENTS-1));
    phi = malloc(sizeof(double*)*NUMPHASES);

    for (int i = 0; i < NUMCOMPONENTS-1; i++)
        comp[i] = malloc(sizeof(double)*MESH_X*MESH_Y*MESH_Z);

    for (int i = 0; i < NUMPHASES; i++)
        phi[i] = malloc(sizeof(double)*MESH_X*MESH_Y*MESH_Z);

    FILE *fp;

    for (long i = start; i <= end; i += step)
    {
        for (int j = 0; j < numProcs; j++)
        {
            sprintf(fileName, "%s/Processor_%d/%s_%07ld.vtk", where, j, name, i);

            if (fp = fopen(fileName, "rb"))
            {
                printf("Reading %s\n", fileName);

                fscanf(fp, "%*[^\n]\n");
                fscanf(fp, "%*[^\n]\n");

                fscanf(fp, "%s", writeformat);

                if (strcmp(writeformat, "ASCII") == 0)
                {
                    if(!(read_vtk_ascii(fp)))
                    {
                        printf("Error reading %s. Verify your input file and the .vtk file\n", fileName);
                        i = end+1;
                    }
                }
                else if (strcmp(writeformat, "BINARY") == 0)
                {
                    if(!(read_vtk_binary(fp)))
                    {
                        printf("Error reading %s. Verify your input file and the .vtk file\n", fileName);
                        i = end+1;
                    }
                }

                calcStats(i);

                found = 1;

                fclose(fp);
            }
            else
            {
                printf("\nCould not find %s.\n", fileName);
                found = 0;
            }
        }

        if (found)
        {
            FILE *fr = fopen(outputName, "ab");

            double radius, velocity;

            if (DIMENSION == 2)
                radius = 0.5*(sqrt(sumInner*DELTA_X*DELTA_Y/M_PI) + sqrt(sumOuter*DELTA_X*DELTA_Y/M_PI));
            else if (DIMENSION == 3)
                radius = 0.5*(pow(sumInner*4.0*DELTA_X*DELTA_Y*DELTA_Z/M_PI/3.0, 0.33333333) + pow(sumOuter*4.0*DELTA_X*DELTA_Y*DELTA_Z/M_PI/3.0, 0.33333333));

            double time = DELTA_t*(double)i;

            if (time == 0)
                velocity = 0.0;
            else
            {
                velocity = (radius - radiusOld)/(time - timeOld);
                radiusOld = radius;
                timeOld = time;
            }

            fprintf(fr, "%lf\t%lf\t%lf\t%lf\t%lf\n", time, radius, velocity, maxC, minC);

            fclose(fr);

            sumInner = 0.0;
            sumOuter = 0.0;
            maxC = -1.5;
            minC = 1.5;
        }
    }

    for (int i = 0; i < NUMCOMPONENTS-1; i++)
        free(comp[i]);
    free(comp);

    for (int i = 0; i < NUMPHASES; i++)
        free(phi[i]);
    free(phi);
}
