#include "SandSim.h"
#include "main.h"

//加快运算速度，一些除法可以先算。
const float invertSpacing = 1.0f/Spacing;
const float initialDensity = 4.0f;

//各种数据的缓存
float particlePos[NumberOfParticles*2]; //粒子的位置，x位置为2*n，y位置为2n+1
float particleVel[NumberOfParticles*2]; //粒子的速度，x速度为2*n，y速度为2n+1

float uVel[(CellNumX + 1) * CellNumY];   //u分量 (x方向)
float vVel[CellNumX * (CellNumY + 1)];   //v分量 (y方向)
float uPrev[(CellNumX + 1) * CellNumY];  //上一帧u
float vPrev[CellNumX * (CellNumY + 1)];  //上一帧v
float uWeights[(CellNumX + 1) * CellNumY];
float vWeights[CellNumX * (CellNumY + 1)];
unsigned int gridType[CellCount]; //网格是气体还是液体

unsigned int Count[CellCount+1U]; //计算粒子碰撞所用的缓存，其中每个entry对应了应该去哪里找粒子。
unsigned int density[CellCount];//框内粒子情况，可用此作为显示。
float particlePosId[NumberOfParticles]; //经过hashgrid排序之后的粒子位置
void printLocation(unsigned int n);

static int clamp_index(int value, int min_value, int max_value) {
    if (value < min_value) {
        return min_value;
    }
    if (value > max_value) {
        return max_value;
    }
    return value;
}

#define U_INDEX(x,y) ((x) * CellNumY + (y))
#define V_INDEX(x,y) ((x) * (CellNumY + 1) + (y))

static void accumulate_u(float x, float y, float value) {
    float fx = x * invertSpacing;
    float fy = y * invertSpacing - 0.5f;
    int x0 = (int)floorf(fx);
    int y0 = (int)floorf(fy);
    x0 = clamp_index(x0, 0, CellNumX - 1);
    y0 = clamp_index(y0, 0, CellNumY - 2);
    float tx = fx - (float)x0;
    float ty = fy - (float)y0;
    tx = fminf(fmaxf(tx, 0.0f), 1.0f);
    ty = fminf(fmaxf(ty, 0.0f), 1.0f);
    int x1 = x0 + 1;
    int y1 = y0 + 1;

    float w00 = (1.0f - tx) * (1.0f - ty);
    float w10 = tx * (1.0f - ty);
    float w01 = (1.0f - tx) * ty;
    float w11 = tx * ty;

    uVel[U_INDEX(x0, y0)] += value * w00;
    uWeights[U_INDEX(x0, y0)] += w00;
    uVel[U_INDEX(x1, y0)] += value * w10;
    uWeights[U_INDEX(x1, y0)] += w10;
    uVel[U_INDEX(x0, y1)] += value * w01;
    uWeights[U_INDEX(x0, y1)] += w01;
    uVel[U_INDEX(x1, y1)] += value * w11;
    uWeights[U_INDEX(x1, y1)] += w11;
}

static void accumulate_v(float x, float y, float value) {
    float fx = x * invertSpacing - 0.5f;
    float fy = y * invertSpacing;
    int x0 = (int)floorf(fx);
    int y0 = (int)floorf(fy);
    x0 = clamp_index(x0, 0, CellNumX - 2);
    y0 = clamp_index(y0, 0, CellNumY - 1);
    float tx = fx - (float)x0;
    float ty = fy - (float)y0;
    tx = fminf(fmaxf(tx, 0.0f), 1.0f);
    ty = fminf(fmaxf(ty, 0.0f), 1.0f);
    int x1 = x0 + 1;
    int y1 = y0 + 1;

    float w00 = (1.0f - tx) * (1.0f - ty);
    float w10 = tx * (1.0f - ty);
    float w01 = (1.0f - tx) * ty;
    float w11 = tx * ty;

    vVel[V_INDEX(x0, y0)] += value * w00;
    vWeights[V_INDEX(x0, y0)] += w00;
    vVel[V_INDEX(x1, y0)] += value * w10;
    vWeights[V_INDEX(x1, y0)] += w10;
    vVel[V_INDEX(x0, y1)] += value * w01;
    vWeights[V_INDEX(x0, y1)] += w01;
    vVel[V_INDEX(x1, y1)] += value * w11;
    vWeights[V_INDEX(x1, y1)] += w11;
}

static float sample_u(float x, float y, const float *grid) {
    float fx = x * invertSpacing;
    float fy = y * invertSpacing - 0.5f;
    int x0 = (int)floorf(fx);
    int y0 = (int)floorf(fy);
    x0 = clamp_index(x0, 0, CellNumX - 1);
    y0 = clamp_index(y0, 0, CellNumY - 2);
    float tx = fx - (float)x0;
    float ty = fy - (float)y0;
    tx = fminf(fmaxf(tx, 0.0f), 1.0f);
    ty = fminf(fmaxf(ty, 0.0f), 1.0f);
    int x1 = x0 + 1;
    int y1 = y0 + 1;

    float w00 = (1.0f - tx) * (1.0f - ty);
    float w10 = tx * (1.0f - ty);
    float w01 = (1.0f - tx) * ty;
    float w11 = tx * ty;

    return grid[U_INDEX(x0, y0)] * w00 +
        grid[U_INDEX(x1, y0)] * w10 +
        grid[U_INDEX(x0, y1)] * w01 +
        grid[U_INDEX(x1, y1)] * w11;
}

static float sample_v(float x, float y, const float *grid) {
    float fx = x * invertSpacing - 0.5f;
    float fy = y * invertSpacing;
    int x0 = (int)floorf(fx);
    int y0 = (int)floorf(fy);
    x0 = clamp_index(x0, 0, CellNumX - 2);
    y0 = clamp_index(y0, 0, CellNumY - 1);
    float tx = fx - (float)x0;
    float ty = fy - (float)y0;
    tx = fminf(fmaxf(tx, 0.0f), 1.0f);
    ty = fminf(fmaxf(ty, 0.0f), 1.0f);
    int x1 = x0 + 1;
    int y1 = y0 + 1;

    float w00 = (1.0f - tx) * (1.0f - ty);
    float w10 = tx * (1.0f - ty);
    float w01 = (1.0f - tx) * ty;
    float w11 = tx * ty;

    return grid[V_INDEX(x0, y0)] * w00 +
        grid[V_INDEX(x1, y0)] * w10 +
        grid[V_INDEX(x0, y1)] * w01 +
        grid[V_INDEX(x1, y1)] * w11;
}

void InitParticles(){
    //在每个格子里生成一个。如果全填满了则不再生成。
    int p_num = 0;

    float half_d = Spacing/2.0f; //半个间距。
    for (int i = 1; i < CellNumY ; i++){
        for(int j = 1; j < CellNumX; j++){

            int particleIndex = p_num;
            particlePos[XID(particleIndex)] = j*Spacing+half_d;
            particlePos[YID(particleIndex)] = i*Spacing+half_d;
            particleVel[XID(particleIndex)] = 0.0f;
            particleVel[YID(particleIndex)] = 0.0f;
            if(p_num++ >= NumberOfParticles){return;} 
            PRINT("generated a particle %d at %f, %f. \n",particleIndex,j*Spacing+half_d,i*Spacing+half_d);
            
        }
    }

    return;
}

/*
 * @brief 根据所提供的加速度仿真粒子的行动
 *
 * @param xAccleration,yAcceleration xy方向的加速度
 * @return none
 *
 * @note 也会仿真墙壁。 
 */ 
void ParticleIntegrate(float xAcceleration, float yAcceleration){
    for(unsigned int i = 0; i < NumberOfParticles ; i++){
        //先算速度
        particleVel[XID(i)] += xAcceleration*dt;
        particleVel[YID(i)] += yAcceleration*dt;
        //再算位置
        particlePos[XID(i)] += particleVel[XID(i)] * dt;
        particlePos[YID(i)] += particleVel[YID(i)] * dt; 

        // 边界。因为网格数据结构的原因，第一行和第一列是没有速度量的，所以得把粒子挤出去，第一行和第一列不放粒子。
        // X方向边界
        if(particlePos[XID(i)] < Spacing){
            particlePos[XID(i)] = Spacing + 0.01f;
            particleVel[XID(i)] *= BOUNCYNESS; // 反弹阻尼
        }
        if(particlePos[XID(i)] >= CellNumX*Spacing-ParticleRadius){
            particlePos[XID(i)] = CellNumX*Spacing-ParticleRadius - 0.01f;
            particleVel[XID(i)] *= BOUNCYNESS;
        }
        
        // Y方向边界
        if(particlePos[YID(i)] < Spacing){
            particlePos[YID(i)] = Spacing + 0.01f;
            particleVel[YID(i)] *= BOUNCYNESS;
        }
        if(particlePos[YID(i)] >= CellNumY*Spacing-ParticleRadius){
            particlePos[YID(i)] = CellNumY*Spacing-ParticleRadius - 0.01f;
            particleVel[YID(i)] *= BOUNCYNESS;
        }

    }
    return;
}

/*
 * @brief 相互推开粒子，模拟粒子之间的碰撞
 *
 * @param nIters 仿真重复的次数
 * @return none
 *
 * @note 也会仿真墙壁。 
 */ 
void PushParticlesApart(unsigned int nIters){

    //清零缓存
    memset(Count,0,sizeof(Count));
    memset(particlePosId,0,sizeof(particlePosId));

    //数数
    for (unsigned int i=0;i<NumberOfParticles;i++){
        float x = particlePos[XID(i)];
        float y = particlePos[YID(i)];

        unsigned int xi = floor(x*invertSpacing); //x和y的网格坐标，即这个粒子在网格中所处的位置
        unsigned int yi = floor(y*invertSpacing); 

        Count[INDEX(xi,yi)] ++;
    }

    unsigned int first = 0;
    
    //partial sum
    for (unsigned int i=0;i<CellCount;i++){
        first += Count[i];
        Count[i] = first;
    } 
    Count[CellCount] = first; //guard
    
    //将排序好的粒子塞入对应的位置缓存中
    for (unsigned int i=0;i<NumberOfParticles;i++){
        float x = particlePos[XID(i)];//取出xy坐标
        float y = particlePos[YID(i)];

        unsigned int xi = floor(x*invertSpacing); //x和y的网格坐标，即这个粒子在网格中所处的位置
        unsigned int yi = floor(y*invertSpacing); 

        int gridindex = --Count[INDEX(xi,yi)]; //对应位置先减1，再放入排序好了的目标缓存
        particlePosId[gridindex] = i; //放入排序的缓存
        PRINT("index is %d %d",xi,yi);
    }
    
//开始把粒子推开

float minDist = 2.0f * ParticleRadius;
float minDist2 = minDist * minDist;
do {
    for (unsigned int i = 0; i < NumberOfParticles; i++) {
        float px = particlePos[XID(i)];
        float py = particlePos[YID(i)];

        unsigned int pxi = floor(px * invertSpacing);
        unsigned int pyi = floor(py * invertSpacing);
        unsigned int x0 = fmax(pxi - 1, 0);
        unsigned int y0 = fmax(pyi - 1, 0);
        unsigned int x1 = fmin(pxi + 1, CellNumX- 1);
        unsigned int y1 = fmin(pyi + 1, CellNumY - 1);

        for (int xi = x0; xi <= x1; xi++) {
            for (int yi = y0; yi <= y1; yi++) {
                unsigned int cellNr = INDEX(xi,yi);
                unsigned int first = Count[cellNr];
                unsigned int last = Count[cellNr + 1];
                for (int j = first; j < last; j++) { //遍历xi，yi格子中的所有粒子，并将其和第i个粒子对比，看是否有碰撞。
                    int id = particlePosId[j];
                    if (id == i)
                        continue;
                    float qx = particlePos[XID(id)];
                    float qy = particlePos[YID(id)];

                    float dx = qx - px;
                    float dy = qy - py;
                    float d2 = dx * dx + dy * dy;
                    if (d2 > minDist2) 
                        continue;
                    else if(d2 == 0.0){
                        PRINT("kiss detected on %d,%d",i,id);
                        continue;
                    }
                    PRINT("bounce detected on %d,%d",i,id);
                    float d = sqrt(d2);
                    float s = 0.5 * (minDist - d) / d;
                    dx *= s;
                    dy *= s;
                    particlePos[XID(i)] -= dx;
                    particlePos[YID(i)] -= dy;
                    particlePos[XID(id)] += dx;
                    particlePos[YID(id)] += dy;
                }

        }
    }
    

    }

    nIters--; 
    } while ( nIters > 0 );
    
    //推完了，再处理一次边界
    for(int i = 0; i < NumberOfParticles; i++){
        if(particlePos[XID(i)] < Spacing){
            particlePos[XID(i)] = Spacing + 0.01f;
        }
        if(particlePos[XID(i)] >= CellNumX*Spacing-ParticleRadius){
            particlePos[XID(i)] = CellNumX*Spacing-ParticleRadius - 0.01f;
        }
        
        // Y方向边界
        if(particlePos[YID(i)] < Spacing){
            particlePos[YID(i)] = Spacing + 0.01f;
        }
        if(particlePos[YID(i)] >= CellNumY*Spacing-ParticleRadius){
            particlePos[YID(i)] = CellNumY*Spacing-ParticleRadius - 0.01f;
        }
    }
    return;
}

/*
 * @brief 寻找有粒子的框。
 *
 * @param none
 * @return none
 *
 * @note 
 */ 
 void density_update(){

    memset(density,0,sizeof(density));
    //数数
    for (unsigned int i=0;i<NumberOfParticles;i++){
        float x = particlePos[XID(i)];
        float y = particlePos[YID(i)];
        unsigned int xi = floor(x*invertSpacing); //x和y的网格坐标，即这个粒子在网格中所处的位置
        unsigned int yi = floor(y*invertSpacing); 
        density[INDEX(xi,yi)] ++;
    }
    return;
 }

 /*
 * @brief 粒子速度换算到网格。
 *
 * @param none
 * @return none
 *
 * @note 
 */ 
void particles_to_grid() {
    //保存上一帧投影后的网格速度供FLIP使用
    memcpy(uPrev, uVel, sizeof(uVel));
    memcpy(vPrev, vVel, sizeof(vVel));

    //清零密度，网格类型和网格速度
    memset(gridType,0,sizeof(gridType));
    memset(uVel,0,sizeof(uVel));
    memset(vVel,0,sizeof(vVel));
    memset(uWeights,0,sizeof(uWeights));
    memset(vWeights,0,sizeof(vWeights));

    for(int p=0; p<NumberOfParticles; p++){

        float x = particlePos[XID(p)];
        float y = particlePos[YID(p)];
        //计算粒子对应的网格位置和delta x y
        int xcell = (int)floorf(x * invertSpacing);
        int ycell = (int)floorf(y * invertSpacing);
        xcell = clamp_index(xcell, 0, CellNumX - 1);
        ycell = clamp_index(ycell, 0, CellNumY - 1);

        PRINT("calculate particle %d at %.2f,%.2f. xcell: %d, ycell: %d.\n",p,x,y,xcell,ycell);
        
        int gridIndex = INDEX(xcell,ycell);//对应的网格在缓存中的位置。
        gridType[gridIndex] = 1;//对应网格是液体网格

        accumulate_u(x, y, particleVel[XID(p)]);
        accumulate_v(x, y, particleVel[YID(p)]);

    }
    //传递完了之后应该除掉总权重。
    for (int i=0;i<CellNumX + 1;i++){
       for (int j=0;j<CellNumY;j++){
        if(uWeights[U_INDEX(i,j)]){ //显然如果这里是空气那就不用算了，不应该处以0
            uVel[U_INDEX(i,j)] /= uWeights[U_INDEX(i,j)];
         }
       }
    }
    for (int i=0;i<CellNumX;i++){
       for (int j=0;j<CellNumY + 1;j++){
        if(vWeights[V_INDEX(i,j)]){ //显然如果这里是空气那就不用算了，不应该处以0
            vVel[V_INDEX(i,j)] /= vWeights[V_INDEX(i,j)];
        }
       }
    }
    return;
}

 /*
 * @brief 计算不可压缩性
 *
 * @param none
 * @return none
 *
 * @note 
 */ 
void compute_grid_forces(unsigned int nIters) {
    //根据网格速度，计算无法压缩的液体,nIters为迭代数量，需要算n次。
    do{
        //对每个grid进行迭代
        for(int xcell = 1; xcell < CellNumX - 1; xcell++){
            for(int ycell = 1; ycell < CellNumY - 1; ycell++){
                if(gridType[INDEX(xcell,ycell)] == 0){
                    continue;//如果这个grid是空气，就跳过。
                } else {
                //如果grid是液体，可以计算divergence。首先求出这个grid的xy坐标。

                //计算divergence
                float d = uVel[U_INDEX(xcell + 1, ycell)]
                    - uVel[U_INDEX(xcell, ycell)]
                    + vVel[V_INDEX(xcell, ycell + 1)]
                    - vVel[V_INDEX(xcell, ycell)];//计算divergence
                d = d * overRelaxiation * dt * invertSpacing;//overrelax并缩放到时间步长

                if(density[INDEX(xcell,ycell)] > initialDensity){
                    float compression = (float)density[INDEX(xcell,ycell)]-initialDensity;
                    compression *= stiffnessCoefficient;
                    d = d - compression;//drift相关,大概是不应该这样写的，密度对每个框框还有个权重。
                }

                if (isnan(d)){
                    ;
                    printf("something wrong with divergence.\n xcell: %d, ycell: %d, gridType: %d\n",xcell,ycell,gridType[INDEX(xcell,ycell)]);
                }

                //处理墙壁与空气邻接
                int s1 = gridType[INDEX(xcell - 1, ycell)] ? 1 : 0;
                int s2 = gridType[INDEX(xcell, ycell - 1)] ? 1 : 0;
                int s3 = gridType[INDEX(xcell + 1, ycell)] ? 1 : 0;
                int s4 = gridType[INDEX(xcell, ycell + 1)] ? 1 : 0;

                //因为只考虑边界，所以可以用判断的形式确定s的值，这意味着我们只有四面墙。
                if(xcell == 1){s1=0;}
                if(ycell == 1){s2=0;}
                if(xcell == CellNumX-2){s3=0;}
                if(ycell == CellNumY-2){s4=0;}

                int s_sum = s1 + s2 + s3 + s4;
                if (s_sum == 0) {
                    continue;
                }
                float s = 1.0f / s_sum;

                uVel[U_INDEX(xcell, ycell)] += d * s1 * s;
                uVel[U_INDEX(xcell + 1, ycell)] -= d * s3 * s;
                vVel[V_INDEX(xcell, ycell)] += d * s2 * s;
                vVel[V_INDEX(xcell, ycell + 1)] -= d * s4 * s;



            }
            }
        }
        nIters--;       
    } while (nIters>0);

    //墙壁边界条件：边界速度为0
    for (int y = 0; y < CellNumY; y++) {
        uVel[U_INDEX(0, y)] = 0.0f;
        uVel[U_INDEX(CellNumX, y)] = 0.0f;
    }
    for (int x = 0; x < CellNumX; x++) {
        vVel[V_INDEX(x, 0)] = 0.0f;
        vVel[V_INDEX(x, CellNumY)] = 0.0f;
    }
    return;
}
 /*
 * @brief 把网格速度还给粒子
 *
 * @param none
 * @return none
 *
 * @note 
 */ 

void grid_to_particles() {
    for(int p=0; p<NumberOfParticles; p++){
        float pic_vx = sample_u(particlePos[XID(p)], particlePos[YID(p)], uVel);
        float pic_vy = sample_v(particlePos[XID(p)], particlePos[YID(p)], vVel);
        float prev_vx = sample_u(particlePos[XID(p)], particlePos[YID(p)], uPrev);
        float prev_vy = sample_v(particlePos[XID(p)], particlePos[YID(p)], vPrev);

        float flip_vx = particleVel[XID(p)] + (pic_vx - prev_vx);
        float flip_vy = particleVel[YID(p)] + (pic_vy - prev_vy);
        const float flipBlend = 0.2f;

        particleVel[XID(p)] = pic_vx * (1.0f - flipBlend) + flip_vx * flipBlend;
        particleVel[YID(p)] = pic_vy * (1.0f - flipBlend) + flip_vy * flipBlend;
    }

    return;
}

void visualize_grid() {
    //显示
    // 初始化缓冲区为全'-'
    char visual_buffer[CellNumY][CellNumX+1];
    memset(visual_buffer, '-', sizeof(visual_buffer));

    printf("\e[1;1H\e[2J");//清空屏幕；。
    
    // 标记粒子位置
    for(int p=0; p<NumberOfParticles; p++){
        int x = floor(particlePos[XID(p)]*invertSpacing);
        int y = floor(particlePos[YID(p)]*invertSpacing);

        visual_buffer[y][x] = 'x';
    }

    // 添加字符串终止符
    for(int j=0; j<CellNumY; j++){
        visual_buffer[j][CellNumX] = '\0';
    }

    
    printf("PIC Simulation (X: %d, Y: %d, Particles: %d)\n", 
           CellNumX, CellNumY, NumberOfParticles);
    for(int j=0; j<CellNumY; j++){
        printf("%s\n", visual_buffer[j]);
    }
    printLocation(0);
    printLocation(1);
    fflush(stdout); 
    
}

void printLocation(unsigned int n){
    printf("Particle %d location:%.2f,%.2f, speed is %.2f,%.2f \n",n,particlePos[XID(n)],particlePos[YID(n)],particleVel[XID(n)],particleVel[YID(n)]);
}
