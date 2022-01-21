#ifndef GRADIENTS_H_
#define GRADIENTS_H_
void swaplayers() {
//Swapping links, phase-field
    tmp             = gradient[-1];
    gradient[-1]    = gradient[0];
    gradient[0]     = gradient[1];
    gradient[1]     = gradient[2];
    gradient[2]     = tmp;
}
#endif