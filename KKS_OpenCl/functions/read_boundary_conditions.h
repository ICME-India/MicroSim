#ifndef READ_BOUNDARY_CONDITIONS_H_
#define READ_BOUNDARY_CONDITIONS_H_

void read_boundary_conditions(char *argv[]) {
  FILE *fr;
  int i;
  char tempbuff[1000];
  
  char tmpstr1[100];
  char tmpstr2[100];
  char **tmp;
  
  bool decision;
  
  char *str1, *str2, *token, *subtoken;
  char *saveptr1, *saveptr2;
  
  long k, j;
  long index;
  long length;
  long phase;
  
  fr = fopen(argv[1], "rt");
  
  if(fr == NULL) {
    printf("file %s not found", argv[1]);
  }
  
  while(fgets(tempbuff,1000,fr)) {
    sscanf(tempbuff, "%100s = %100[^;];", tmpstr1, tmpstr2);
//     printf("%s\n",  tmpstr1);
//     printf("%s\n",  tmpstr2);
    if(tmpstr1[0] != '#') {
      if ((strcmp(tmpstr1, "BOUNDARY") == 0) && (NUMPHASES > 0)) {
        initialize_boundary_conditions(tmpstr2);
      }
      //else if ((strcmp(tmpstr1, "BOUNDARY_VALUE") == 0) && (NUMPHASES > 0)) {
        //initialize_boundary_points_values(tmpstr2);
      //}
    }
  }
  fclose(fr);
  
  char outfile[100];
  
  strcpy(tmpstr2, argv[1]);
  
  strcpy(tmpstr1,strtok(tmpstr2, "."));
  
  sprintf(outfile, "%s.bd", tmpstr1);
  
  fr = fopen(outfile, "w");
  
  PRINT_BOUNDARY_CONDITIONS(fr);
  
  fclose(fr);
}
#endif