#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct neutron {
  double x,y,z,w;
};

double r2(struct neutron s) {
  return(s.x * s.x + s.y * s.y + s.z*s.z);
}

double rand01() {
  return ((double)rand())/((double)RAND_MAX);
}

void c3() {
  double lmfp = 1;
  int i=0;
  
  double maxbin = 5;
  double binsize = 0.05;
  int nbins = maxbin/binsize;
  int *histogram = malloc(sizeof(int)*nbins);
  for(i=0; i<nbins; i++)
    histogram[nbins] = 0;
  
  int iterations = 1000000;
  for(i=0; i<iterations; i++) {
    double l = -lmfp * log(1.0 - rand01());
    int bin = l/binsize;
    if(bin >= nbins) continue;
    histogram[bin] ++;
  }
  
  for(i=0; i<nbins; i++)
    printf("%f %d\n", i*binsize, histogram[i]);
  
  free(histogram);
}

void gain() {
  int n = 1000000;
  int totalrounds = 30;
  struct neutron *neutrons = malloc(sizeof(struct neutron)*n);
  double *weights = malloc(sizeof(double)*totalrounds);
  int i=0, rnd=0;

  double a = 0.075;
  double alphaprimed = 0.9;
  double a_in = a*alphaprimed;
  double growthfactor = 1.231;
  double lmfp = 0.0268;
  
  
  // Initialize
  for(i=0; i<n; i++) {
    neutrons[i].x = (a+a_in)*0.5;
    neutrons[i].y = 0;
    neutrons[i].z = 0;
    neutrons[i].w = 1;
  }
  
  for(rnd=0; rnd < totalrounds; rnd++) {
    // Update
    for(i=0; i<n; i++) {
      double costheta = rand()%2 == 0 ? rand01() : -rand01();
      double phi = rand01()*M_PI*2.0;
      double l = -lmfp * log(1.0 - rand01());
    
      // check for hollow region
      double hollowa = 1;
      double hollowb = 2.0*(neutrons[i].x * sqrt(1.0 - costheta*costheta) * cos(phi) + neutrons[i].y * sqrt(1.0 - costheta*costheta) * sin(phi) + neutrons[i].z * costheta);
      double hollowc = neutrons[i].x*neutrons[i].x + neutrons[i].y*neutrons[i].y + neutrons[i].z*neutrons[i].z - a_in*a_in;
      if(hollowb*hollowb - 4.0*hollowa*hollowc < 0) {
        
      }
      else {        
        double l1 = (-hollowb + sqrt(hollowb*hollowb - 4.0*hollowa*hollowc))/2.0/hollowa;
        double l2 = (-hollowb - sqrt(hollowb*hollowb - 4.0*hollowa*hollowc))/2.0/hollowa;
        if(l1 > 0 && l2 > 0) {
          double lshort = (l1 < l2) ? l1 : l2;
          double llong = (l1 < l2) ? l2 : l1;
          if(l > lshort) { //goes into hollow
            l += llong-lshort;
          }
        }
      }

      // move it
      neutrons[i].x += l * sqrt(1.0 - costheta*costheta) * cos(phi);
      neutrons[i].y += l * sqrt(1.0 - costheta*costheta) * sin(phi);
      neutrons[i].z += l * costheta;
    
      if(r2(neutrons[i]) > a*a) { // Escaped
        neutrons[i].w = 0;
      }
      else {
        neutrons[i].w *= growthfactor;
      }
    }
  
    // Add up weights
    weights[rnd] = 0;
    int nonzeroweights = 0;
    int zeroweights = 0;
    struct neutron** zeroneutrons = malloc(sizeof(struct neutron*)*n);
    struct neutron** nonzeroneutrons = malloc(sizeof(struct neutron*)*n);
    
    for(i=0; i<n; i++) {
      if(neutrons[i].w != 0) {
        weights[rnd] += neutrons[i].w;
        nonzeroneutrons[nonzeroweights] = &neutrons[i];
        nonzeroweights++;
      }
      else {
        zeroneutrons[zeroweights] = &neutrons[i];
        zeroweights++;
      }
    }
    
    // Splitting
    if(2.0*nonzeroweights < n) {
      for(i=0; i<nonzeroweights; i++) {
        nonzeroneutrons[i]->w /= 2;
        zeroneutrons[i]->x = nonzeroneutrons[i]->x;
        zeroneutrons[i]->y = nonzeroneutrons[i]->y;
        zeroneutrons[i]->z = nonzeroneutrons[i]->z;
        zeroneutrons[i]->w = nonzeroneutrons[i]->w;
      }
    }
    
    free(zeroneutrons);
    free(nonzeroneutrons);
  }
  
  for(rnd=0; rnd<totalrounds; rnd++) {
    double prev;
    if(rnd == 0)
      prev = n;
    else
      prev = weights[rnd-1];

    printf("%f (ratio = %f)\n", weights[rnd], weights[rnd]/prev);
  }
  
  double ratio = 0;
  int ratio_norm = 0;
  for(rnd=totalrounds-10; rnd<totalrounds; rnd++) {
    ratio += weights[rnd] / weights[rnd-1];
    ratio_norm++;
  }
  ratio /= (double)ratio_norm;
  printf("Average over last %d: %f\n", ratio_norm, ratio);
  
  
  free(neutrons);
  free(weights);
  
}

int main() {
  //c3();
  gain();
  
  return 0;
}