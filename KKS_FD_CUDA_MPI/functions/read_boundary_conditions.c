#include "inputReader.h"

void assign_buffer_points_conditions(long j, int BOUNDARY_LEFT, int BOUNDARY_RIGHT, int BOUNDARY_FRONT, int BOUNDARY_BACK, int BOUNDARY_TOP, int BOUNDARY_BOTTOM,
                                     domainInfo *simDomain, controls *simControls, simParameters *simParams)
{
    //     //X-
    //     simControls->boundary[0][j].points[0] = 2;
    //     simControls->boundary[0][j].points[1] = 1;
    //     simControls->boundary[0][j].points[2] = 0;
    //     //X-
    //
    //     //X+
    //     simControls->boundary[1][j].points[0] = rows_x-3;
    //     simControls->boundary[1][j].points[1] = rows_x-2;
    //     simControls->boundary[1][j].points[2] = rows_x-1;
    //     //X+
    //
    //     if (simDomain->DIMENSION >= 2)
    //     {
    //         //Y+
    //         simControls->boundary[2][j].points[0] = rows_y-3;
    //         simControls->boundary[2][j].points[1] = rows_y-2;
    //         simControls->boundary[2][j].points[2] = rows_y-1;
    //         //Y+
    //
    //         //Y-
    //         simControls->boundary[3][j].points[0] = 2;
    //         simControls->boundary[3][j].points[1] = 1;
    //         simControls->boundary[3][j].points[2] = 0;
    //         //Y-
    //     }
    //
    //     if (simDomain->DIMENSION == 3)
    //     {
    //         //Z+
    //         simControls->boundary[4][j].points[0] = rows_z-3;
    //         simControls->boundary[4][j].points[1] = rows_z-2;
    //         simControls->boundary[4][j].points[2] = rows_z-1;
    //         //Z+
    //
    //         //Z-
    //         simControls->boundary[5][j].points[0] = 2;
    //         simControls->boundary[5][j].points[1] = 1;
    //         simControls->boundary[5][j].points[2] = 0;
    //         //Z-
    //     }

    simControls->boundary[0][j].type = BOUNDARY_LEFT;
    simControls->boundary[1][j].type = BOUNDARY_RIGHT;
    simControls->boundary[2][j].type = BOUNDARY_FRONT;
    simControls->boundary[3][j].type = BOUNDARY_BACK;
    simControls->boundary[4][j].type = BOUNDARY_TOP;
    simControls->boundary[5][j].type = BOUNDARY_BOTTOM;

    if ((simControls->boundary[0][j].type == 3) || (simControls->boundary[1][j].type == 3) ) {
        BOUNDARY_LEFT    = 3;
        BOUNDARY_RIGHT   = 3;
        simControls->boundary[0][j].type = 3;
        simControls->boundary[1][j].type = 3;
    }
    if ((simControls->boundary[2][j].type == 3) || (simControls->boundary[3][j].type == 3) ) {
        BOUNDARY_FRONT   = 3;
        BOUNDARY_BACK    = 3;
        simControls->boundary[2][j].type = 3;
        simControls->boundary[3][j].type = 3;
    }
    if ((simControls->boundary[4][j].type == 3) || (simControls->boundary[5][j].type == 3) ) {
        BOUNDARY_TOP     = 3;
        BOUNDARY_BOTTOM  = 3;
        simControls->boundary[4][j].type = 3;
        simControls->boundary[5][j].type = 3;
    }

    if (simDomain->DIMENSION == 2) {
        simControls->boundary[4][j].type = 0;
        simControls->boundary[5][j].type = 0;
    }

    long i;

    for (i = 0; i < 6; i++)
    {
        if ((simControls->boundary[0][j].type == 2) || (simControls->boundary[0][j].type == 0))
        {
            fprintf(stdout, "DIRICHLET and FREE boundary conditions have not yet been implemented. Please wait for next release. Will default to PERIODIC\n");
        }
    }
}


void initialize_boundary_conditions(char *tmpstr, domainInfo *simDomain, controls *simControls, simParameters *simParams)
{
    char **tmp;
    char *str1, *str2, *token;
    char *saveptr1, *saveptr2;
    int i, l, j;

    long len = 7;
    long phase;
    int  DIM;

    tmp = (char**)malloc(sizeof(char*)*len);
    for (i = 0; i < len; ++i) {
        tmp[i] = (char*)malloc(sizeof(char)*10);
    }

    for (i = 0, str1 = tmpstr; ; i++, str1 = NULL) {
        token = strtok_r(str1, "{,}", &saveptr1);
        if (token == NULL)
            break;
        strcpy(tmp[i],token);
    }

    int BOUNDARY_LEFT   = atoi(tmp[2]);
    int BOUNDARY_RIGHT  = atoi(tmp[1]);
    int BOUNDARY_FRONT  = atoi(tmp[3]);
    int BOUNDARY_BACK   = atoi(tmp[4]);
    int BOUNDARY_TOP    = atoi(tmp[5]);
    int BOUNDARY_BOTTOM = atoi(tmp[6]);

    if (strcmp(tmp[0], "phi") == 0) {
        assign_buffer_points_conditions(0, BOUNDARY_LEFT, BOUNDARY_RIGHT, BOUNDARY_FRONT, BOUNDARY_BACK, BOUNDARY_TOP, BOUNDARY_BOTTOM, simDomain, simControls, simParams);
    }
    if (strcmp(tmp[0], "mu") == 0) {
        assign_buffer_points_conditions(1, BOUNDARY_LEFT, BOUNDARY_RIGHT, BOUNDARY_FRONT, BOUNDARY_BACK, BOUNDARY_TOP, BOUNDARY_BOTTOM, simDomain, simControls, simParams);
    }
    if (strcmp(tmp[0], "c") == 0) {
        assign_buffer_points_conditions(2, BOUNDARY_LEFT, BOUNDARY_RIGHT, BOUNDARY_FRONT, BOUNDARY_BACK, BOUNDARY_TOP, BOUNDARY_BOTTOM, simDomain, simControls, simParams);
    }
    if (strcmp(tmp[0], "T") == 0) {
        assign_buffer_points_conditions(3, BOUNDARY_LEFT, BOUNDARY_RIGHT, BOUNDARY_FRONT, BOUNDARY_BACK, BOUNDARY_TOP, BOUNDARY_BOTTOM, simDomain, simControls, simParams);
    }

    for (i = 0; i < len; ++i) {
        free(tmp[i]);
    }
    free(tmp);
    tmp = NULL;
}

void PRINT_BOUNDARY_CONDITIONS(FILE *fp, domainInfo *simDomain, controls *simControls, simParameters *simParams) {
    char **key, var[6];

    int i, j;

    key = (char**)malloc(sizeof(char*)*6);
    for (i = 0; i < 6; ++i) {
        key[i] = (char*)malloc(sizeof(char)*10);
    }

    char Scalars[4][20] = {"phi", "mu", "c", "T"};

    for (j = 0; j < 3; j++) {
        sprintf(key[0], "BOUNDARY_LEFT[%s]",   Scalars[j]);
        sprintf(key[1], "BOUNDARY_RIGHT[%s]",  Scalars[j]);
        sprintf(key[2], "BOUNDARY_FRONT[%s]",  Scalars[j]);
        sprintf(key[3], "BOUNDARY_BACK[%s]",   Scalars[j]);
        sprintf(key[4], "BOUNDARY_TOP[%s]",    Scalars[j]);
        sprintf(key[5], "BOUNDARY_BOTTOM[%s]", Scalars[j]);

        var[0] = simControls->boundary[0][j].type;
        var[1] = simControls->boundary[1][j].type;
        var[2] = simControls->boundary[2][j].type;
        var[3] = simControls->boundary[3][j].type;
        var[4] = simControls->boundary[4][j].type;
        var[5] = simControls->boundary[5][j].type;

        for (i=0 ; i<6; i++) {
            if (var[i] == 0) {
                fprintf(fp, "%s = %s\n", key[i], "FREE");
                fprintf(fp,"\n");
                fprintf(fp, "DIRICHLET and FREE boundary conditions have not yet been implemented. Please wait for next release. Will default to PERIODIC\n");
            }
            if (var[i] == 1) {
                fprintf(fp, "%s = %s\n", key[i], "NEUMANN");
                fprintf(fp,"\n");
            }
            if (var[i] == 2) {
                fprintf(fp, "%s = %s\n", key[i], "DIRICHLET");
                fprintf(fp,"\n");
                fprintf(fp, "DIRICHLET and FREE boundary conditions have not yet been implemented. Please wait for next release. Will default to PERIODIC\n");
            }
            if (var[i] == 3) {
                fprintf(fp, "%s = %s\n", key[i], "PERIODIC");
                fprintf(fp,"\n");

            }
        }
    }
    for (i = 0; i < 6; ++i) {
        free(key[i]);
    }
    free(key);
    key = NULL;
}

void read_boundary_conditions(domainInfo *simDomain, controls *simControls,
                              simParameters *simParams, int rank, char *argv[])
{
    FILE *fr, *fp;
    if (fr = fopen(argv[1], "rt"))
    {
        if (!(rank))
            printf("\nReading boundary conditions from %s\n", argv[1]);
    }
    else
    {
        if (!(rank))
            printf("\nFile %s not found\n", argv[1]);
    }

    char tempbuff[1000];
    char outfile[100];
    char tmpstr1[100];
    char tmpstr2[100];

    while (fgets(tempbuff, 1000, fr))
    {
        sscanf(tempbuff,"%100s = %100[^;];", tmpstr1, tmpstr2);
        if (tmpstr1[0] != '#')
        {
            if ((strcmp(tmpstr1, "BOUNDARY") == 0) && (simDomain->numPhases > 0))
            {
                initialize_boundary_conditions(tmpstr2, simDomain, simControls, simParams);
            }
            else if ((strcmp(tmpstr1, "BOUNDARY_VALUE") == 0) && (simDomain->numPhases > 0))
            {
                //initialize_boundary_points_values(tmpstr2);
            }
        }
    }

    fclose(fr);

    strcpy(tmpstr2, argv[1]);
    strcpy(tmpstr1,strtok(tmpstr2, "."));
    sprintf(outfile, "%s.bd", tmpstr1);

    fp = fopen(outfile, "w");

    PRINT_BOUNDARY_CONDITIONS(fp, simDomain, simControls, simParams);

    fclose(fp);
}
