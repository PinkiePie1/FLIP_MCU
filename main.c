#include "main.h"
#include "SandSim.h"
#include <stdio.h>

void main(){
    int dummy=0;
    printf("hi.\n");

    InitParticles();



    for (int i=0;i<9000;i++){
 //       if (i%30==0){
 //           visualize_grid();
 //           printf("%d\n",dummy++);  
 //       }
        ParticleIntegrate(0.0f,-9.8f);
        PushParticlesApart(2);
        density_update();
        particles_to_grid();
        compute_grid_forces(20);
        grid_to_particles();
    }

for (int j =0;j<5;j++)
    {
    for (int i=0;i<70;i++){
        if (i%5==0){
            visualize_grid();
            printf("%d\n",dummy++);
           
        }
        ParticleIntegrate(5.0f,-9.8f);
        PushParticlesApart(2);
        density_update();
        particles_to_grid();
        compute_grid_forces(20);
        grid_to_particles();
    }

    for (int i=0;i<70;i++){
        if (i%5==0){
            visualize_grid();
            printf("%d\n",dummy++);         
        }
        ParticleIntegrate(-5.0f,-9.8f);
        PushParticlesApart(2);
        density_update();
        particles_to_grid();
        compute_grid_forces(20);
        grid_to_particles();
    }

	}
    
    return;
}
