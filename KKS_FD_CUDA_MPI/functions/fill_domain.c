#include "inputReader.h"

void readFill(fillParameters *simFill, char *argv[], int rank)
{
    FILE *fr;

    if (fr = fopen(argv[2], "rt"))
    {
        if (!(rank))
            printf("\nReading filling parameters from %s\n", argv[2]);
    }
    else
    {
        if (!(rank))
            printf("\n File %s not found\n", argv[2]);
    }

    long i;
    char tempbuff[1000];
    char tmpstr1[100], tmpstr2[100];

    char **tmp;
    char *str1, *token;
    char *saveptr1;

    simFill->countFill = 0;

    while (fgets(tempbuff, 1000, fr))
    {
        sscanf(tempbuff, "%100s = %100[^;];", tmpstr1, tmpstr2);
        if (tmpstr1[0] != '#')
        {
            if ((strcmp(tmpstr1, "FILLCYLINDER") == 0) || (strcmp(tmpstr1, "FILLSPHERE") == 0) || (strcmp(tmpstr1, "FILLCUBE") == 0) || (strcmp(tmpstr1, "FILLCYLINDERRANDOM") == 0) || (strcmp(tmpstr1, "FILLSPHERERANDOM") == 0))
            {
                simFill->countFill++;
            }
        }
    }

    fclose(fr);

    simFill->fillType     = (fill*)malloc(sizeof(fill)*simFill->countFill);
    simFill->xC           = (long*)malloc(sizeof(long)*simFill->countFill);
    simFill->yC           = (long*)malloc(sizeof(long)*simFill->countFill);
    simFill->zC           = (long*)malloc(sizeof(long)*simFill->countFill);
    simFill->xS           = (long*)malloc(sizeof(long)*simFill->countFill);
    simFill->xE           = (long*)malloc(sizeof(long)*simFill->countFill);
    simFill->yS           = (long*)malloc(sizeof(long)*simFill->countFill);
    simFill->yE           = (long*)malloc(sizeof(long)*simFill->countFill);
    simFill->zS           = (long*)malloc(sizeof(long)*simFill->countFill);
    simFill->zE           = (long*)malloc(sizeof(long)*simFill->countFill);
    simFill->radius       = (long*)malloc(sizeof(long)*simFill->countFill);
    simFill->phase        = (long*)malloc(sizeof(long)*simFill->countFill);
    simFill->major_axis   = (double*)malloc(sizeof(double)*simFill->countFill);
    simFill->eccentricity = (double*)malloc(sizeof(double)*simFill->countFill);
    simFill->rot_angle    = (double*)malloc(sizeof(double)*simFill->countFill);
    simFill->seed         = (long*)malloc(sizeof(long)*simFill->countFill);
    simFill->volFrac      = (double*)malloc(sizeof(double)*simFill->countFill);
    simFill->shieldDist   = (long*)malloc(sizeof(long)*simFill->countFill);
    simFill->radVar       = (double*)malloc(sizeof(double)*simFill->countFill);

    long j = 0;

    fr = fopen(argv[2], "rt");

    while (fgets(tempbuff, 1000, fr))
    {
        sscanf(tempbuff, "%100s = %100[^;];", tmpstr1, tmpstr2);

        if (tmpstr1[0] != '#')
        {
            if (j >= simFill->countFill)
                break;

            if (strcmp(tmpstr1, "FILLCYLINDER") == 0)
            {
                tmp = (char**)malloc(sizeof(char*)*6);
                for (i = 0; i < 6; i++)
                    tmp[i] = (char*)malloc(sizeof(char)*10);

                for (i = 0, str1 = tmpstr2; ; i++, str1 = NULL)
                {
                    token = strtok_r(str1, "{,}", &saveptr1);
                    if (token == NULL)
                        break;
                    strcpy(tmp[i],token);
                }

                simFill->fillType[j] = FILLCYLINDER;
                simFill->phase[j]    = atol(tmp[0]);
                simFill->xC[j]       = atol(tmp[1]);
                simFill->yC[j]       = atol(tmp[2]);
                simFill->zS[j]       = atol(tmp[3]);
                simFill->zE[j]       = atol(tmp[4]);
                simFill->radius[j]   = atol(tmp[5]);

                if (!(rank))
                    printf("Read cylinder parameters - (%ld, %ld, %ld, %ld, %ld, %ld)\n", simFill->phase[j], simFill->xC[j], simFill->yC[j], simFill->zS[j], simFill->zE[j], simFill->radius[j]);

                j++;

                for (i = 0; i < 6; i++)
                    free(tmp[i]);
                free(tmp);
            }

            else if (strcmp(tmpstr1, "FILLSPHERE") == 0)
            {
                tmp = (char**)malloc(sizeof(char*)*5);
                for (i = 0; i < 5; i++)
                    tmp[i] = (char*)malloc(sizeof(char)*10);

                for (i = 0, str1 = tmpstr2; ; i++, str1 = NULL)
                {
                    token = strtok_r(str1, "{,}", &saveptr1);
                    if (token == NULL)
                        break;
                    strcpy(tmp[i],token);
                }

                simFill->fillType[j] = FILLSPHERE;
                simFill->phase[j]    = atol(tmp[0]);
                simFill->xC[j]       = atol(tmp[1]);
                simFill->yC[j]       = atol(tmp[2]);
                simFill->zC[j]       = atol(tmp[3]);
                simFill->radius[j]   = atol(tmp[4]);

                j++;

                if (!(rank))
                    printf("Read sphere parameters\n");

                for (i = 0; i < 5; i++)
                    free(tmp[i]);
                free(tmp);
            }

            else if (strcmp(tmpstr1, "FILLCUBE") == 0)
            {
                tmp = (char**)malloc(sizeof(char*)*7);
                for (i = 0; i < 7; i++)
                    tmp[i] = (char*)malloc(sizeof(char)*10);

                for (i = 0, str1 = tmpstr2; ; i++, str1 = NULL)
                {
                    token = strtok_r(str1, "{,}", &saveptr1);
                    if (token == NULL)
                        break;
                    strcpy(tmp[i],token);
                }

                simFill->fillType[j] = FILLCUBE;
                simFill->phase[j]    = atol(tmp[0]);
                simFill->xS[j]       = atol(tmp[1]);
                simFill->yS[j]       = atol(tmp[2]);
                simFill->zS[j]       = atol(tmp[3]);
                simFill->xE[j]       = atol(tmp[4]);
                simFill->yE[j]       = atol(tmp[5]);
                simFill->zE[j]       = atol(tmp[6]);

                j++;

                if (!(rank))
                    printf("Read cube parameters\n");

                for (i = 0; i < 7; i++)
                    free(tmp[i]);
                free(tmp);
            }

            else if (strcmp(tmpstr1, "FILLELLIPSE") == 0)
            {
                tmp = (char**)malloc(sizeof(char*)*8);
                for (i = 0; i < 8; i++)
                    tmp[i] = (char*)malloc(sizeof(char)*10);

                for (i = 0, str1 = tmpstr2; ; i++, str1 = NULL)
                {
                    token = strtok_r(str1, "{,}", &saveptr1);
                    if (token == NULL)
                        break;
                    strcpy(tmp[i],token);
                }

                simFill->fillType[j]     = FILLELLIPSE;
                simFill->phase[j]        = atol(tmp[0]);
                simFill->xC[j]           = atol(tmp[1]);
                simFill->yC[j]           = atol(tmp[2]);
                simFill->zC[j]           = atol(tmp[3]);
                simFill->major_axis[j]   = atol(tmp[5]);
                simFill->eccentricity[j] = atol(tmp[6]);
                simFill->rot_angle[j]    = atol(tmp[7]);

                j++;

                if (!(rank))
                    printf("Read ellipse parameters\n");

                for (i = 0; i < 8; i++)
                    free(tmp[i]);
                free(tmp);
            }

            else if (strcmp(tmpstr1, "FILLCYLINDERRANDOM") == 0)
            {
                tmp = (char**)malloc(sizeof(char*)*6);
                for (i = 0; i < 6; i++)
                    tmp[i] = (char*)malloc(sizeof(char)*10);

                for (i = 0, str1 = tmpstr2; ; i++, str1 = NULL)
                {
                    token = strtok_r(str1, "{,}", &saveptr1);
                    if (token == NULL)
                        break;
                    strcpy(tmp[i],token);
                }

                simFill->fillType[j]   = FILLCYLINDERRANDOM;
                simFill->phase[j]      = atol(tmp[0]);
                simFill->radius[j]     = atol(tmp[1]);
                simFill->volFrac[j]    = atof(tmp[2]);
                simFill->shieldDist[j] = atol(tmp[3]);
                simFill->radVar[j]     = atof(tmp[4]);

                j++;

                if (!(rank))
                    printf("Read random cylinder filling parameters\n");

                for (i = 0; i < 6; i++)
                    free(tmp[i]);
                free(tmp);
            }

            else if (strcmp(tmpstr1, "FILLSPHERERANDOM") == 0)
            {
                tmp = (char**)malloc(sizeof(char*)*6);
                for (i = 0; i < 6; i++)
                    tmp[i] = (char*)malloc(sizeof(char)*10);

                for (i = 0, str1 = tmpstr2; ; i++, str1 = NULL)
                {
                    token = strtok_r(str1, "{,}", &saveptr1);
                    if (token == NULL)
                        break;
                    strcpy(tmp[i],token);
                }

                simFill->fillType[j]   = FILLSPHERERANDOM;
                simFill->phase[j]      = atol(tmp[0]);
                simFill->radius[j]     = atol(tmp[1]);
                simFill->volFrac[j]    = atof(tmp[2]);
                simFill->shieldDist[j] = atol(tmp[3]);
                simFill->radVar[j]     = atof(tmp[4]);

                j++;

                if (!(rank))
                    printf("Read random sphere filling parameters\n");

                for (i = 0; i < 6; i++)
                    free(tmp[i]);
                free(tmp);
            }

            else
                printf("Did not find valid input in %s\n", argv[2]);
        }
    }

    fclose(fr);
}
