#include "SandSim.h"
#include "main.h"

const float invertSpacing = 1.0f / Spacing;
const float particleInvSpacing = 1.0f / (2.2f * ParticleRadius);

#define FLUID_CELL 0U
#define AIR_CELL 1U
#define SOLID_CELL 2U

float particlePos[NumberOfParticles * 2];
float particleVel[NumberOfParticles * 2];

static float uVel[CellCount];
static float vVel[CellCount];
static float uPrev[CellCount];
static float vPrev[CellCount];
static float uWeight[CellCount];
static float vWeight[CellCount];

static float pressure[CellCount];
static float particleDensity[CellCount];
static float solidMask[CellCount];
static unsigned int cellType[CellCount];
static float particleRestDensity = 0.0f;

#define PNumX (CellNumX * 2U)
#define PNumY (CellNumY * 2U)
#define PCellCount (PNumX * PNumY)

static unsigned int numCellParticles[PCellCount];
static unsigned int firstCellParticle[PCellCount + 1U];
static unsigned int cellParticleIds[NumberOfParticles];

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

static float clampf_local(float value, float min_value, float max_value) {
    if (value < min_value) {
        return min_value;
    }
    if (value > max_value) {
        return max_value;
    }
    return value;
}

static void setup_solid_mask(void) {
    for (unsigned int x = 0; x < CellNumX; x++) {
        for (unsigned int y = 0; y < CellNumY; y++) {
            float s = 1.0f;
            if (x == 0U || x == CellNumX - 1U || y == 0U) {
                s = 0.0f;
            }
            solidMask[INDEX(x, y)] = s;
        }
    }
}

void InitParticles() {
    memset(particleVel, 0, sizeof(particleVel));
    memset(uVel, 0, sizeof(uVel));
    memset(vVel, 0, sizeof(vVel));
    memset(uPrev, 0, sizeof(uPrev));
    memset(vPrev, 0, sizeof(vPrev));
    memset(pressure, 0, sizeof(pressure));
    memset(particleDensity, 0, sizeof(particleDensity));
    particleRestDensity = 0.0f;
    setup_solid_mask();

    const float h = Spacing;
    const float r = ParticleRadius;
    const float dx = 2.0f * r;
    const float dy = 0.86602540378f * dx; // sqrt(3)/2

    unsigned int p_num = 0;
    for (unsigned int i = 0; i < CellNumX && p_num < NumberOfParticles; i++) {
        for (unsigned int j = 0; j < CellNumY && p_num < NumberOfParticles; j++) {
            float px = h + r + dx * (float)i + ((j % 2U == 0U) ? 0.0f : r);
            float py = h + r + dy * (float)j;
            if (px > (CellNumX - 1U) * h - r || py > (CellNumY - 1U) * h - r) {
                continue;
            }
            particlePos[XID(p_num)] = px;
            particlePos[YID(p_num)] = py;
            p_num++;
        }
    }

    for (; p_num < NumberOfParticles; p_num++) {
        particlePos[XID(p_num)] = h + r;
        particlePos[YID(p_num)] = h + r;
    }
}

void ParticleIntegrate(float xAcceleration, float yAcceleration) {
    const float minX = Spacing + ParticleRadius;
    const float maxX = (CellNumX - 1U) * Spacing - ParticleRadius;
    const float minY = Spacing + ParticleRadius;
    const float maxY = (CellNumY - 1U) * Spacing - ParticleRadius;

    for (unsigned int i = 0; i < NumberOfParticles; i++) {
        particleVel[XID(i)] += xAcceleration * dt;
        particleVel[YID(i)] += yAcceleration * dt;
        particlePos[XID(i)] += particleVel[XID(i)] * dt;
        particlePos[YID(i)] += particleVel[YID(i)] * dt;

        float x = particlePos[XID(i)];
        float y = particlePos[YID(i)];

        if (x < minX) {
            x = minX;
            particleVel[XID(i)] = 0.0f;
        }
        if (x > maxX) {
            x = maxX;
            particleVel[XID(i)] = 0.0f;
        }
        if (y < minY) {
            y = minY;
            particleVel[YID(i)] = 0.0f;
        }
        if (y > maxY) {
            y = maxY;
            particleVel[YID(i)] = 0.0f;
        }

        particlePos[XID(i)] = x;
        particlePos[YID(i)] = y;
    }
}

void PushParticlesApart(unsigned int nIters) {
    memset(numCellParticles, 0, sizeof(numCellParticles));

    for (unsigned int i = 0; i < NumberOfParticles; i++) {
        int xi = clamp_index((int)floorf(particlePos[XID(i)] * particleInvSpacing), 0, (int)PNumX - 1);
        int yi = clamp_index((int)floorf(particlePos[YID(i)] * particleInvSpacing), 0, (int)PNumY - 1);
        numCellParticles[(unsigned int)xi * PNumY + (unsigned int)yi]++;
    }

    unsigned int first = 0;
    for (unsigned int i = 0; i < PCellCount; i++) {
        first += numCellParticles[i];
        firstCellParticle[i] = first;
    }
    firstCellParticle[PCellCount] = first;

    for (unsigned int i = 0; i < NumberOfParticles; i++) {
        int xi = clamp_index((int)floorf(particlePos[XID(i)] * particleInvSpacing), 0, (int)PNumX - 1);
        int yi = clamp_index((int)floorf(particlePos[YID(i)] * particleInvSpacing), 0, (int)PNumY - 1);
        unsigned int cellNr = (unsigned int)xi * PNumY + (unsigned int)yi;
        firstCellParticle[cellNr]--;
        cellParticleIds[firstCellParticle[cellNr]] = i;
    }

    const float minDist = 2.0f * ParticleRadius;
    const float minDist2 = minDist * minDist;

    for (unsigned int iter = 0; iter < nIters; iter++) {
        for (unsigned int i = 0; i < NumberOfParticles; i++) {
            float px = particlePos[XID(i)];
            float py = particlePos[YID(i)];

            int pxi = (int)floorf(px * particleInvSpacing);
            int pyi = (int)floorf(py * particleInvSpacing);
            int x0 = (pxi - 1 > 0) ? pxi - 1 : 0;
            int y0 = (pyi - 1 > 0) ? pyi - 1 : 0;
            int x1 = (pxi + 1 < (int)PNumX - 1) ? pxi + 1 : (int)PNumX - 1;
            int y1 = (pyi + 1 < (int)PNumY - 1) ? pyi + 1 : (int)PNumY - 1;

            for (int xi = x0; xi <= x1; xi++) {
                for (int yi = y0; yi <= y1; yi++) {
                    unsigned int cellNr = (unsigned int)xi * PNumY + (unsigned int)yi;
                    unsigned int firstIdx = firstCellParticle[cellNr];
                    unsigned int lastIdx = firstCellParticle[cellNr + 1U];
                    for (unsigned int j = firstIdx; j < lastIdx; j++) {
                        unsigned int id = cellParticleIds[j];
                        if (id == i) {
                            continue;
                        }

                        float qx = particlePos[XID(id)];
                        float qy = particlePos[YID(id)];
                        float dx = qx - px;
                        float dy = qy - py;
                        float d2 = dx * dx + dy * dy;
                        if (d2 > minDist2 || d2 == 0.0f) {
                            continue;
                        }

                        float d = sqrtf(d2);
                        float s = 0.5f * (minDist - d) / d;
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
    }
}

void density_update(void) {
    memset(particleDensity, 0, sizeof(particleDensity));

    const float h = Spacing;
    const float h1 = invertSpacing;
    const float h2 = 0.5f * h;

    for (unsigned int i = 0; i < NumberOfParticles; i++) {
        float x = clampf_local(particlePos[XID(i)], h, (CellNumX - 1U) * h);
        float y = clampf_local(particlePos[YID(i)], h, (CellNumY - 1U) * h);

        int x0 = (int)floorf((x - h2) * h1);
        int y0 = (int)floorf((y - h2) * h1);
        x0 = clamp_index(x0, 0, (int)CellNumX - 2);
        y0 = clamp_index(y0, 0, (int)CellNumY - 2);
        int x1 = x0 + 1;
        int y1 = y0 + 1;

        float tx = ((x - h2) - x0 * h) * h1;
        float ty = ((y - h2) - y0 * h) * h1;
        float sx = 1.0f - tx;
        float sy = 1.0f - ty;

        particleDensity[INDEX(x0, y0)] += sx * sy;
        particleDensity[INDEX(x1, y0)] += tx * sy;
        particleDensity[INDEX(x1, y1)] += tx * ty;
        particleDensity[INDEX(x0, y1)] += sx * ty;
    }

    if (particleRestDensity == 0.0f) {
        float sum = 0.0f;
        unsigned int numFluid = 0;
        for (unsigned int i = 0; i < CellCount; i++) {
            if (cellType[i] == FLUID_CELL) {
                sum += particleDensity[i];
                numFluid++;
            }
        }
        if (numFluid > 0U) {
            particleRestDensity = sum / (float)numFluid;
        }
    }
}

void particles_to_grid(void) {
    memcpy(uPrev, uVel, sizeof(uVel));
    memcpy(vPrev, vVel, sizeof(vVel));

    memset(uVel, 0, sizeof(uVel));
    memset(vVel, 0, sizeof(vVel));
    memset(uWeight, 0, sizeof(uWeight));
    memset(vWeight, 0, sizeof(vWeight));

    for (unsigned int i = 0; i < CellCount; i++) {
        cellType[i] = (solidMask[i] == 0.0f) ? SOLID_CELL : AIR_CELL;
    }

    for (unsigned int i = 0; i < NumberOfParticles; i++) {
        int xi = clamp_index((int)floorf(particlePos[XID(i)] * invertSpacing), 0, (int)CellNumX - 1);
        int yi = clamp_index((int)floorf(particlePos[YID(i)] * invertSpacing), 0, (int)CellNumY - 1);
        unsigned int cellNr = INDEX((unsigned int)xi, (unsigned int)yi);
        if (cellType[cellNr] == AIR_CELL) {
            cellType[cellNr] = FLUID_CELL;
        }
    }

    const float h = Spacing;
    const float h1 = invertSpacing;
    const float h2 = 0.5f * h;

    for (int component = 0; component < 2; component++) {
        float dx = (component == 0) ? 0.0f : h2;
        float dy = (component == 0) ? h2 : 0.0f;
        float *f = (component == 0) ? uVel : vVel;
        float *w = (component == 0) ? uWeight : vWeight;

        for (unsigned int i = 0; i < NumberOfParticles; i++) {
            float x = clampf_local(particlePos[XID(i)], h, (CellNumX - 1U) * h);
            float y = clampf_local(particlePos[YID(i)], h, (CellNumY - 1U) * h);

            int x0 = (int)floorf((x - dx) * h1);
            int y0 = (int)floorf((y - dy) * h1);
            x0 = clamp_index(x0, 0, (int)CellNumX - 2);
            y0 = clamp_index(y0, 0, (int)CellNumY - 2);
            int x1 = x0 + 1;
            int y1 = y0 + 1;

            float tx = ((x - dx) - x0 * h) * h1;
            float ty = ((y - dy) - y0 * h) * h1;
            float sx = 1.0f - tx;
            float sy = 1.0f - ty;

            float w0 = sx * sy;
            float w1 = tx * sy;
            float w2 = tx * ty;
            float w3 = sx * ty;

            float pv = particleVel[2 * i + (unsigned int)component];
            unsigned int nr0 = INDEX((unsigned int)x0, (unsigned int)y0);
            unsigned int nr1 = INDEX((unsigned int)x1, (unsigned int)y0);
            unsigned int nr2 = INDEX((unsigned int)x1, (unsigned int)y1);
            unsigned int nr3 = INDEX((unsigned int)x0, (unsigned int)y1);

            f[nr0] += pv * w0; w[nr0] += w0;
            f[nr1] += pv * w1; w[nr1] += w1;
            f[nr2] += pv * w2; w[nr2] += w2;
            f[nr3] += pv * w3; w[nr3] += w3;
        }

        for (unsigned int i = 0; i < CellCount; i++) {
            if (w[i] > 0.0f) {
                f[i] /= w[i];
            }
        }

        for (unsigned int x = 0; x < CellNumX; x++) {
            for (unsigned int y = 0; y < CellNumY; y++) {
                unsigned int idx = INDEX(x, y);
                unsigned int solid = (cellType[idx] == SOLID_CELL) ? 1U : 0U;
                if (component == 0) {
                    unsigned int leftSolid = (x > 0U && cellType[INDEX(x - 1U, y)] == SOLID_CELL) ? 1U : 0U;
                    if (solid || leftSolid) {
                        uVel[idx] = uPrev[idx];
                    }
                } else {
                    unsigned int bottomSolid = (y > 0U && cellType[INDEX(x, y - 1U)] == SOLID_CELL) ? 1U : 0U;
                    if (solid || bottomSolid) {
                        vVel[idx] = vPrev[idx];
                    }
                }
            }
        }
    }
}

void compute_grid_forces(unsigned int nIters) {
    memset(pressure, 0, sizeof(pressure));
    memcpy(uPrev, uVel, sizeof(uVel));
    memcpy(vPrev, vVel, sizeof(vVel));

    const float cp = 1000.0f * Spacing / dt;

    for (unsigned int iter = 0; iter < nIters; iter++) {
        for (unsigned int x = 1; x < CellNumX - 1U; x++) {
            for (unsigned int y = 1; y < CellNumY - 1U; y++) {
                unsigned int center = INDEX(x, y);
                if (cellType[center] != FLUID_CELL) {
                    continue;
                }

                unsigned int left = INDEX(x - 1U, y);
                unsigned int right = INDEX(x + 1U, y);
                unsigned int bottom = INDEX(x, y - 1U);
                unsigned int top = INDEX(x, y + 1U);

                float sx0 = solidMask[left];
                float sx1 = solidMask[right];
                float sy0 = solidMask[bottom];
                float sy1 = solidMask[top];
                float s = sx0 + sx1 + sy0 + sy1;
                if (s == 0.0f) {
                    continue;
                }

                float div = uVel[right] - uVel[center] + vVel[top] - vVel[center];

                if (particleRestDensity > 0.0f) {
                    float compression = particleDensity[center] - particleRestDensity;
                    if (compression > 0.0f) {
                        div -= compression;
                    }
                }

                float p = -div / s;
                p *= overRelaxiation;
                pressure[center] += cp * p;

                uVel[center] -= sx0 * p;
                uVel[right] += sx1 * p;
                vVel[center] -= sy0 * p;
                vVel[top] += sy1 * p;
            }
        }
    }
}

void grid_to_particles(void) {
    const float flipRatio = 0.9f;
    const float h = Spacing;
    const float h1 = invertSpacing;
    const float h2 = 0.5f * h;

    for (int component = 0; component < 2; component++) {
        float dx = (component == 0) ? 0.0f : h2;
        float dy = (component == 0) ? h2 : 0.0f;
        float *f = (component == 0) ? uVel : vVel;
        float *prevF = (component == 0) ? uPrev : vPrev;
        int offset = (component == 0) ? (int)CellNumY : 1;

        for (unsigned int i = 0; i < NumberOfParticles; i++) {
            float x = clampf_local(particlePos[XID(i)], h, (CellNumX - 1U) * h);
            float y = clampf_local(particlePos[YID(i)], h, (CellNumY - 1U) * h);

            int x0 = clamp_index((int)floorf((x - dx) * h1), 0, (int)CellNumX - 2);
            int y0 = clamp_index((int)floorf((y - dy) * h1), 0, (int)CellNumY - 2);
            int x1 = x0 + 1;
            int y1 = y0 + 1;

            float tx = ((x - dx) - x0 * h) * h1;
            float ty = ((y - dy) - y0 * h) * h1;
            float sx = 1.0f - tx;
            float sy = 1.0f - ty;

            float d0 = sx * sy;
            float d1 = tx * sy;
            float d2 = tx * ty;
            float d3 = sx * ty;

            unsigned int nr0 = INDEX((unsigned int)x0, (unsigned int)y0);
            unsigned int nr1 = INDEX((unsigned int)x1, (unsigned int)y0);
            unsigned int nr2 = INDEX((unsigned int)x1, (unsigned int)y1);
            unsigned int nr3 = INDEX((unsigned int)x0, (unsigned int)y1);

            float valid0 = (cellType[nr0] != AIR_CELL || (int)nr0 - offset >= 0 && cellType[(unsigned int)((int)nr0 - offset)] != AIR_CELL) ? 1.0f : 0.0f;
            float valid1 = (cellType[nr1] != AIR_CELL || (int)nr1 - offset >= 0 && cellType[(unsigned int)((int)nr1 - offset)] != AIR_CELL) ? 1.0f : 0.0f;
            float valid2 = (cellType[nr2] != AIR_CELL || (int)nr2 - offset >= 0 && cellType[(unsigned int)((int)nr2 - offset)] != AIR_CELL) ? 1.0f : 0.0f;
            float valid3 = (cellType[nr3] != AIR_CELL || (int)nr3 - offset >= 0 && cellType[(unsigned int)((int)nr3 - offset)] != AIR_CELL) ? 1.0f : 0.0f;

            float d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;
            if (d <= 0.0f) {
                continue;
            }

            float picV = (valid0 * d0 * f[nr0] + valid1 * d1 * f[nr1] + valid2 * d2 * f[nr2] + valid3 * d3 * f[nr3]) / d;
            float corr = (valid0 * d0 * (f[nr0] - prevF[nr0]) + valid1 * d1 * (f[nr1] - prevF[nr1])
                        + valid2 * d2 * (f[nr2] - prevF[nr2]) + valid3 * d3 * (f[nr3] - prevF[nr3])) / d;
            float oldV = particleVel[2 * i + (unsigned int)component];
            float flipV = oldV + corr;
            particleVel[2 * i + (unsigned int)component] = (1.0f - flipRatio) * picV + flipRatio * flipV;
        }
    }
}

void visualize_grid() {
    char visual_buffer[CellNumY][CellNumX + 1];
    memset(visual_buffer, '-', sizeof(visual_buffer));

    printf("\e[1;1H\e[2J");

    for (unsigned int p = 0; p < NumberOfParticles; p++) {
        int x = clamp_index((int)floorf(particlePos[XID(p)] * invertSpacing), 0, (int)CellNumX - 1);
        int y = clamp_index((int)floorf(particlePos[YID(p)] * invertSpacing), 0, (int)CellNumY - 1);
        visual_buffer[y][x] = 'x';
    }

    for (unsigned int j = 0; j < CellNumY; j++) {
        visual_buffer[j][CellNumX] = '\0';
    }

    printf("FLIP Simulation (X: %d, Y: %d, Particles: %d)\n", CellNumX, CellNumY, NumberOfParticles);
    for (unsigned int j = 0; j < CellNumY; j++) {
        printf("%s\n", visual_buffer[j]);
    }
    printLocation(0);
    printLocation(1);
    fflush(stdout);
}

void printLocation(unsigned int n) {
    printf("Particle %d location:%.2f,%.2f, speed is %.2f,%.2f \n", n, particlePos[XID(n)], particlePos[YID(n)], particleVel[XID(n)], particleVel[YID(n)]);
}
