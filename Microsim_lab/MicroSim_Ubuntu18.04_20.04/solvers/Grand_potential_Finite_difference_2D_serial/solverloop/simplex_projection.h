void projection_on_simplex(double *divphi) {
  long a, b;
  Deltaphi = 0.0;
  count_phases = 0;
  double Deltaphi_alpha = 0.0;
  //Find the number of phases which are okay.
  double sum_phib = 0.0;
  for (a=0; a < NUMPHASES; a++) {
    if ((fabs(divphi[a]) > 0.0) && (gridinfo[center].phia[a] + grad->deltaphi[a]) > 0.0  && (gridinfo[center].phia[a] + grad->deltaphi[a]) < 1.0) {
      count_phases++;
      sum_phib += gridinfo[center].phia[a];
    }
  }
  
  for (a=0; a < NUMPHASES; a++) {
    if ((gridinfo[center].phia[a] + grad->deltaphi[a]) < 0.0) {
      Deltaphi_alpha = fabs(grad->deltaphi[a] + gridinfo[center].phia[a]);
      grad->deltaphi[a] += Deltaphi_alpha;
      Deltaphi          += Deltaphi_alpha;
    }
  }
  if (Deltaphi > 0.0) {
    for (b=0; b < NUMPHASES; b++) {
      if (((fabs(divphi[b]) > 0.0)) && ((gridinfo[center].phia[b] + grad->deltaphi[b]) > 0.0)  && ((gridinfo[center].phia[b] + grad->deltaphi[b]) < 1.0)) {
        if (fabs(sum_phib) > 0.0) {
          grad->deltaphi[b] -= Deltaphi*(gridinfo[center].phia[b])/(sum_phib);
        } else {
          grad->deltaphi[b] -= Deltaphi/(count_phases);
        }
      }
    }
  }
  for (a=0; a < NUMPHASES; a++) {
    if (fabs(divphi[a]) > 0.0) {
      if ((gridinfo[center].phia[a] + grad->deltaphi[a]) > 1.0) {
        grad->deltaphi[a]  = (1.0-gridinfo[center].phia[a]);
        //Correct all the other phases due to this correction,
        //If you are bulk, all other phases must go to zero.
        for (b=0; b < NUMPHASES; b++) {
          if ((fabs(divphi[b]) > 0.0) && b!=a) {
            grad->deltaphi[b] = -gridinfo[center].phia[b];
          }
        }
        Deltaphi = 0.0;
        break;
      }
    }
  }
}
void projection_on_simplex_without_weights(double *divphi) {
  long a, b;
  Deltaphi = 0.0;
  count_phases = 0;
  
  //Find the number of phases which are okay.
  for (a=0; a < NUMPHASES; a++) {
    if ((fabs(divphi[a]) > 0.0) && (gridinfo[center].phia[a] + grad->deltaphi[a]) > 0.0  && (gridinfo[center].phia[a] + grad->deltaphi[a]) < 1.0) {
      count_phases++;
    }
  }
  
  for (a=0; a < NUMPHASES; a++) {
   if (fabs(divphi[a]) > 0.0) {
    if ((gridinfo[center].phia[a] + grad->deltaphi[a]) > 1.0) {
// 				      Deltaphi =  (grad->deltaphi[a] - (1.0-gridinfo[center].phia[a]));
// 				      grad->deltaphi[a] -= Deltaphi;
      grad->deltaphi[a]  = (1.0-gridinfo[center].phia[a]);
      
      //Correct all the other phases due to this correction,
      //If you are bulk, all other phases must go to zero.
      for (b=0; b < NUMPHASES; b++) {
        if ((fabs(divphi[b]) > 0.0) && b!=a) {
          grad->deltaphi[b] = -gridinfo[center].phia[b];
        }
      }
      break;
    }
    if ((gridinfo[center].phia[a] + grad->deltaphi[a]) < 0.0) {
      Deltaphi = fabs(grad->deltaphi[a] + gridinfo[center].phia[a]);
      grad->deltaphi[a] += Deltaphi;
      
      //Correct all the phases due to the correction in the given phase
        for (b=0; b < NUMPHASES; b++) {
          if (((fabs(divphi[b]) > 0.0) && b!=a) && (gridinfo[center].phia[b] + grad->deltaphi[b]) > 0.0  && (gridinfo[center].phia[b] + grad->deltaphi[b]) < 1.0) {
            grad->deltaphi[b] -= Deltaphi/(count_phases);
          }
        }
      }
    }
  }
}