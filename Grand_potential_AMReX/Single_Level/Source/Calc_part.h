#ifndef _CALC_PART_H_
#define _CALC_PART_H_

#include "Variables.h"

void Calc_partition(int &max_x,int &max_y,int &max_z){
   
   #if(AMREX_SPACEDIM==2)
        for(int j=1; j<=numworkers; j++){
            if(numworkers%j == 0){
                numb.push_back(j);
            }
        }

        if(numb.size()%2 != 0){
            val1 = floor((numb.size())/2);
            val2 = floor((numb.size())/2);
        }
        else{
            val1 = floor((float(numb.size())-1.0)/2.0);
            val2 = ceil((float(numb.size())-1.0)/2.0); 
        }

        Print()<<"Partitions are: "<< numb[val1] <<" and "<<numb[val2]<<"\n";
        max_x = ceil((Real(ncellx))/(Real(numb[val1])));
        max_y = ceil((Real(ncelly))/(Real(numb[val2])));
        max_z = 0;
    #endif

    #if(AMREX_SPACEDIM>2)
        while(numworkers%2==0){
           numb.push_back(2);
            numworkers=numworkers/2;
        }

        for(int p=3; p<=sqrt(numworkers); p=p+2){
            while(numworkers%p==0){
                numb.push_back(p);
                numworkers=numworkers/p;
            }
        }

        if(numworkers>2){
            numb.push_back(numworkers);
        }

        reverse(numb.begin(),numb.end());

        if(numb.size()==2){
            numb.push_back(1);
        }

        if(numb.size()==1){
            numb.push_back(1);
            numb.push_back(1);
        }

        Array1D<int,0,2> ral{numb[0],numb[1],numb[2]};
        int dal{3};
        int p{2};

        while(p>=0){ 
            if(dal==numb.size()){
                break;
            }       
            if(dal<numb.size()){
                ral(p)=ral(p)*numb[dal];
                dal++;
                p--;
            }
            if(p==-1 && dal<numb.size()){
                p=2;
            }
        }

        Print()<<"Partitions are: "<< ral(0) <<" , "<<ral(1)<<" , "<<ral(2)<<"\n";
        max_x = ceil((Real(ncellx))/(Real(ral(0))));
        max_y = ceil((Real(ncelly))/(Real(ral(1))));
        max_z = ceil((Real(ncellz))/(Real(ral(2))));

    #endif
}

#endif