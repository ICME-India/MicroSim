


void FunctionF_4_SplineCPU(struct propmatf4spline *propf4splin, double *temp, int DELTATflag) { 
  
  int ip, is, is1, is2, y;
  double tem[rows_y], DELTAT;
  
  if ( DELTATflag ) { 
    DELTAT = ( pfmvar.deltat ) * ( -pfmdat.TGRADIENT * pfmdat.velocity ); 
    for ( y = 0; y < rows_y; y++ ) { 
      tem[y] = temp[y] + DELTAT; 
    }
  }
  else { 
    
    for ( y = 0; y < rows_y; y++ ) { 
      tem[y] = temp[y]; 
    }
    
  }
    
  
  for ( y = 0; y < rows_y; y++ ) { 
    
    function_A(tem[y], ceq); 
    
    for ( ip = 0; ip < npha; ip++ ) { 
      for ( is1 = 0; is1 < nsol; is1++ ) { 
        for ( is2 = 0; is2 <nsol; is2++ ) {
          
          propf4splin[y].A[ip][is1][is2] = A[ip][is1][is2]; 
          
        }
      }
    }
    
    for ( ip = 0; ip < npha; ip++ ) { 
      for ( is = 0; is < nsol; is++ ) { 
        
        propf4splin[y].B[ip][is] =  function_B(tem[y], is, ip);
        
      }
    }
    
    for ( ip = 0; ip < npha; ip++ ) { 
      
      propf4splin[y].C[a] = function_C(tem[y], ip);
      
    }
    
  }
  
}

    
    
    
    
