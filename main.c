#include "main.h"
#include "SandSim.h"
#include <stdio.h>

void main(){
    int dummy=0;
    printf("hi.\n");

    InitParticles();

    for (int i=0;i<600;i++){
        ParticleIntegrate(SIM_IQ(0.0f), SIM_IQ(9.8f));
        PushParticlesApart(2);
        particles_to_grid();
        density_update();
        compute_grid_forces(20);
        grid_to_particles();
    }

    visualize_grid();
    printf("%d\n",dummy++);  

for (int j =0;j<3;j++)
    {
    for (int i=0;i<70;i++){
        if (i%5==0){
            visualize_grid();
            printf("%d\n",dummy++);
           
        }
        ParticleIntegrate(SIM_IQ(0.5f), SIM_IQ(9.8f));
        PushParticlesApart(2);
        particles_to_grid();
        density_update();
        compute_grid_forces(20);
        grid_to_particles();
    }

    for (int i=0;i<10;i++){
        if (i%2==0){
            visualize_grid();
            printf("%d\n",dummy++);
           
        }
        ParticleIntegrate(SIM_IQ(0.0f), SIM_IQ(9.8f));
        PushParticlesApart(2);
        particles_to_grid();
        density_update();
        compute_grid_forces(20);
        grid_to_particles();
    }

    for (int i=0;i<20;i++){
        if (i%2==0){
            visualize_grid();
            printf("%d\n",dummy++);         
        }
        ParticleIntegrate(SIM_IQ(-20.0f), SIM_IQ(9.8f));
        PushParticlesApart(2);
        particles_to_grid();
        density_update();
        compute_grid_forces(20);
        grid_to_particles();
    }

    for (int i=0;i<70;i++){
        if (i%2==0){
            visualize_grid();
            printf("%d\n",dummy++);
           
        }
        ParticleIntegrate(SIM_IQ(0.0f), SIM_IQ(9.8f));
        PushParticlesApart(2);
        particles_to_grid();
        density_update();
        compute_grid_forces(20);
        grid_to_particles();
    }

	}

    for (int i=0;i<600;i++){
        ParticleIntegrate(SIM_IQ(0.0f), SIM_IQ(9.8f));
        PushParticlesApart(2);
        particles_to_grid();
        density_update();
        compute_grid_forces(20);
        grid_to_particles();
    }

    visualize_grid();
  
    
    return;
}
