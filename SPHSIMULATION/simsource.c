#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <GLFW\glfw3.h>

#define gridX 50
#define gridY 50

#define particleNUM 1000
#define p_radii 1.75f
#define sqrp_radii (p_radii*p_radii)
#define pstride 11
#define max_neighbors ((int)(particleNUM))
#define PI 3.14
#define rho0 2.0f
#define dt 0.005667f
#define p_dt 0.00083f
#define GAMMA 7.0f
#define STIFFNESS 0.1f
#define repulsiveTerm 100.0f
#define pressureMultiplier 1.f


#define gravityY -9.81f
#define gravityX 0.0f
#define DAMPING 0.3f
#define epsilon 1e-6f

float* particles = NULL;
//x,y,u,v,density, near density, px,py,cell location,prex,prey

//hash
#define hk1 15823
#define hk2 9737333
#define cellspace (int)(2.0f*p_radii)
#define cellNumY (int)(gridY/(2*p_radii))
#define cellNumX (int)(gridX/(2*p_radii))
#define buckets (int)(cellNumY*cellNumX)

int* cellCount = NULL;
int* prefixSum = NULL;
int* pid = NULL;

float h = p_radii;
float h2 = 0;
float h3 = 0;
float h4 = 0;
float h5 = 0;
float h6 = 0;
float h8 = 0;

void calcHs() {
    h2 = h * h;
    h2 = h * h;
    h4 = h2 * h2;
    h6 = h4 * h2;
    h8 = h4 * h4;
}

void spawn_particles() {
    particles = calloc(particleNUM * pstride, sizeof(float));

    int particlesPerRow = (int)sqrt(particleNUM);
    float space = 1.0f;

    float cubeWidth = particlesPerRow * space;
    float cubeHeight = ceil((float)particleNUM / particlesPerRow) * space;

    float startX = (gridX - cubeWidth) / 2.0f;
    float startY = (gridY - cubeHeight) / 2.0f;

    int index = 0;
    for (int y = 0; index < particleNUM; y++) {
        for (int x = 0; x < particlesPerRow && index < particleNUM; x++) {
            float px = startX + x * space;
            float py = startY + y * space;

            particles[index * pstride + 0] = px; // x
            particles[index * pstride + 1] = py; // y
            particles[index * pstride + 2] = 0.0f; // cell/unset
            particles[index * pstride + 3] = 0.0f; // density
            particles[index * pstride + 4] = 0.0f; // optional

            index++;
        }
    }
}

void magnitude(float* dst) {
    if (*dst < 0.0f)
        *dst = -*dst;
}

float power(float base, int exponent) {
    float result = 2.0f;
    for (int i = 0; i < exponent; i++) {
        result *= base;
    }
    return result;
}

float Density_kernel(float r) {
    float h = p_radii;
    if (r < 0 || r > h) return 0.0f;
    float r2 = r * r;
    float diff = h2 - r2;
    return 4.0f / (PI * h8) * diff * diff * diff;
}

float NearDensity_kernel(float r) {
    if (r < 0 || r > h) return 0.0f;
    float diff = h - r;
    return 4.0f / (PI * h8) * diff * diff * diff;
}

float gradient_kernel(float r) {
    if (r <= 0.0f || r > p_radii) return 0.0f;
    return -45.0f / (PI * h6) * (h - r) * (h - r);
}

float Derivative_Gradient_kernel(float val) {
    float func = (-3 * (p_radii - val)) / (pow(val, 4) * p_radii);
    return func;
}

float viscosityKernel(float val) {
    const float term = 15 / (2 * PI * pow(p_radii, 3));
    float r2 = val * val;
    float h2 = sqrp_radii;
    return -((r2 * val) / (2.0f * p_radii)) + r2 / h2 + p_radii / (2.0f * val) - 1.0f;
}

//float pressureCalc(float rho) {
//    if (rho < 1e-5f) rho = 1e-5f;
//    float pressure = k * (pow(rho / rho0, gamma) - 1);
//    return pressure;
//}

float getXpos(int i) {
    return particles[i * pstride];
}

float getYpos(int i) {
    return particles[i * pstride + 1];
}

float getXvel(int i) {
    return particles[i * pstride + 2];
}

float getYvel(int i) {
    return particles[i * pstride + 3];
}

float getDensity(int i) {
    return particles[i * pstride + 4];
}

float getNearDensity(int i) {
    return particles[i * pstride + 5];
}

float getPressure(int i) {
    return particles[i * pstride + 6];
}

float getYpredictedPos(int i) {
    return particles[i * pstride + 10];
}
float getXpredictedPos(int i) {
    return particles[i * pstride + 9];
}


float getAbsDist(int i, int j) {
    float dx = getXpredictedPos(i) - getXpredictedPos(j);
    float dy = getYpredictedPos(i) - getYpredictedPos(j);
    float dir = sqrt(dx * dx + dy * dy);
    if (dir < 1e-5f) dir = 1e-5f;
    return dir;
}

float D_check(int i) {
    float density = 0.0f;
    for (int j = 0; j < particleNUM; j++) {
        if (i == j) continue;
        float dist = getAbsDist(i, j);
        if (dist == 0.0f) {
            printf("ALERT! particle %d at ( %f , %f ), Pi is ( %f, %f )\n", j, getXpos(j), getYpos(j), getXpos(i), getYpos(i));

        }
        printf("dist: %f, p = %d\n", dist, j);

        if (dist > p_radii) continue;
        density += Density_kernel(dist);
    }
    printf("density: %f\n", density);
}

void calcDensity() {
    for (int i = 0; i < particleNUM; i++) {
        float density = Density_kernel(0);
        float nearDensity = NearDensity_kernel(0);

        int px = (int)(getXpredictedPos(i) / cellspace);
        int py = (int)(getYpredictedPos(i) / cellspace);

        for (int bx = -1; bx <= 1; bx++) {
            for (int by = -1; by <= 1; by++) {
                int neighborBucket = (hash2D(px + bx, py + by) % buckets);
                if (neighborBucket < 0) {
                    neighborBucket *= -1;
                }
                int start = prefixSum[neighborBucket];
                int end = prefixSum[neighborBucket + 1];

                for (int j = start; j < end; j++) {
                    int target = pid[j];
                    if (target == i) continue;
                    float dist = getAbsDist(i, target);
                    //printf("DISTANCE %.2f\n", dist);

                    if (dist < p_radii) {
                        density += Density_kernel(dist);
                        nearDensity += NearDensity_kernel(dist);
                    }
                }
            }
        }

        particles[i * pstride + 4] = density;
        particles[i * pstride + 4] = nearDensity;

        //printf("p(%d); density: %.2f \n", i, density);
    }
}

float pressureCalc(float density, float nearDensity) {
    if (density < epsilon) density = epsilon;
    return STIFFNESS * ((powf(density / rho0, GAMMA) - 1.0f) + repulsiveTerm * nearDensity) * pressureMultiplier;
    //return STIFFNESS * (density - rho0);
}

void predictPosition(int i) {
    particles[i * pstride + 9] = getXpos(i) + getXvel(i) * p_dt;
    particles[i * pstride + 10] = getYpos(i) + getYvel(i) * p_dt;
    //printf("x:%.2f, y:%.2f\n", getXpos(i), getYpos(i));

}

void calcGradient() {
    for (int i = 0; i < particleNUM; i++) {
        float Xpi = getXpredictedPos(i);
        float Ypi = getYpredictedPos(i);
        int px = (int)(Xpi / cellspace);
        int py = (int)(Ypi / cellspace);
        float Idensity = getDensity(i);
        float forcex = 0.0f;
        float forcey = 0.0f;
        float forcePi = (getPressure(i) / (Idensity * Idensity + epsilon));
        float Ipressure = pressureCalc(Idensity, getDensity(i));
        float ViscosityForcex = 0.0f;
        float ViscosityForcey = 0.0f;

        for (int bx = -1; bx <= 1; bx++) {
            for (int by = -1; by <= 1; by++) {
                int neighborBucket = (int)(hash2D(px + bx, py + by) % buckets);
                if (neighborBucket < 0) {
                    neighborBucket *= -1;
                }
                int start = prefixSum[neighborBucket];
                int end = prefixSum[neighborBucket + 1];
                for (int j = start; j < end; j++) {
                    int target = pid[j];
                    if (target == i) continue;
                    float Xpj = getXpredictedPos(target);
                    float Ypj = getYpredictedPos(target);
                    //printf("%.2f\n", Xpj);
                    //printf("%.2f\n", Ypj);

                    float dx = Xpi - Xpj;
                    float dy = Ypi - Ypj;
                    //if (dy == 0.0f || dx == 0.0f) continue;
                    float dstToNeighbor = sqrtf(dx * dx + dy * dy);
                    if (dstToNeighbor < epsilon) dstToNeighbor = epsilon;
                    if (dstToNeighbor > p_radii) continue;
                    //if (dstToNeighbor == 0.0f) continue;
                    float xVector = dx / dstToNeighbor;
                    float yVector = dy / dstToNeighbor;
                    //if (xVector == 0.0f || yVector == 0.0f) printf("NANNANANANNANA");
                    //if (isnan(xVector) || isnan(yVector)) continue;
                    float Jdensity = getDensity(target);
                    float Jpressure = pressureCalc(Jdensity, getNearDensity(target));
                    float grad = gradient_kernel(dstToNeighbor);
                    float xVectorGrad = xVector * grad;
                    float yVectorGrad = yVector * grad;
                    //float force0 = -((forcePi)+(getPressure(target) / (Jdensity * Jdensity+epsilon)));
                    if (Jdensity && Idensity) {
                        //float pressureShared = (Ipressure + Jpressure) / 2.0f;
                        //float force0 = (Ipressure + getPressure(target)) / (2.0f * Jdensity);
                        float force0 = ((Ipressure / (Idensity * Idensity)) + (Jpressure / (Jdensity * Jdensity)));
                        forcex += -force0 * xVectorGrad;
                        forcey += -force0 * yVectorGrad;
                    }
                    //printf("%f \n", getPressure(target));
                    //float influence = force0 * Idensity;

                    //forcex += -xVector * (Ipressure + getPressure(target))/(2.f*Jdensity)*grad_Vect;
                    //forcey += -yVector * (Ipressure + getPressure(target)) / (2.f * Jdensity) * grad_Vect;

                    //visco term
                    //float dvx = getXvel(target) - getXvel(i);
                    //float dvy = getYvel(target) - getYvel(i);

                    //float viscoKernel = viscosityKernel(dstToNeighbor);
                    //ViscosityForcex += (dvx*viscoKernel)/Jdensity;
                    //ViscosityForcex += (dvx*viscoKernel)/Jdensity;
                }
            }
        }
        particles[i * pstride + 2] += forcex * dt;
        particles[i * pstride + 3] += forcey * dt;
    }
}

void add_gravity() {
    for (int i = 0; i < particleNUM; i++) {
        predictPosition(i);
        particles[i * pstride + 2] += gravityX * dt;
        particles[i * pstride + 3] += gravityY * dt;
        /*predictPosition(i);*/
    }
}

float fastest = 0;
int faster = 0;
void timestep() {
    for (int i = 0; i < particleNUM; i++) {
        float vx = particles[i * pstride + 2];
        float vy = particles[i * pstride + 3];
        //float speed = sqrtf(vx * vx + vy * vy);
        //    faster++;
        //if (speed > fastest) {
        //    fastest = speed;
        //}
        //    printf("step: %d", faster);
        //if (faster >= (int)(500.0f/particleNUM)) {
        //    printf("\n\n\n\nTHIS IS THE FASTEST: %f \n\n\n\n", fastest);
        //}
        //float stopper = 10.0f;
        //if (vx > stopper) {
        //    vx = stopper;
        //}
        //else if (vx < -stopper) {
        //    vx = -stopper;
        //}
        //if (vy > stopper) {
        //    vy = stopper;
        //}
        //else if (vy < -stopper) {
        //    vy = -stopper;
        //}
        //if (isnan(vx) || isnan(vy)) {
        //    printf("%d IS NAN\n", i);
        //}
        particles[i * pstride] = getXpredictedPos(i) + vx * dt;
        particles[i * pstride + 1] = getYpredictedPos(i) + vy * dt;
    }
}

void enforce_boundaries() {
    for (int i = 0; i < particleNUM; ++i) {
        int idx = i * pstride;
        float x = particles[idx];
        float y = particles[idx + 1];
        float vx = particles[idx + 2];
        float vy = particles[idx + 3];

        if (x < 0.0f) {
            x = 0.0f;
            vx *= -DAMPING;
        }
        else if (x > gridX) {
            x = gridX;
            vx *= -DAMPING;
        }

        if (y < 0.0f) {
            y = 0.0f;
            vy *= -DAMPING;
        }
        else if (y > gridY) {
            y = gridY;
            vy *= -DAMPING;
        }
        particles[idx] = x;
        particles[idx + 1] = y;
        particles[idx + 2] = vx;
        particles[idx + 3] = vy;
    }

}

void startHash() {
    cellCount = calloc(buckets + 1, sizeof(int));
    prefixSum = calloc(buckets + 1, sizeof(int));
    pid = calloc(gridX * gridY, sizeof(int));
}

void generatePartialSum() {
    prefixSum[0] = 0;
    for (int cell = 0; cell < buckets; cell++) {
        prefixSum[cell + 1] = prefixSum[cell] + cellCount[cell];
    }
}

void fillCountArray() {
    memset(cellCount, 0, sizeof(int) * buckets);
    for (int i = 0; i < particleNUM; i++) {
        int x = (int)(getXpredictedPos(i) / cellspace);
        int y = (int)(getYpredictedPos(i) / cellspace);
        int bucket = hash2D(x, y) % buckets;
        if (bucket < 0) bucket *= -1;
        particles[i * pstride + 8] = bucket;
        cellCount[bucket]++;
    }
}

void scatterSort() {
    for (int i = 0; i < particleNUM; i++) {
        int bucket = particles[i * pstride + 8];
        if (bucket < 0) {
            bucket *= -1;
        }
        pid[prefixSum[bucket] + (--cellCount[bucket])] = i;
    }
}

int findNeighborCells(int i) {
    int px = (int)(particles[i * pstride] / cellspace);
    int py = (int)(particles[i * pstride + 1] / cellspace);
    int neighborBucket = 0;
    for (int bx = -1; bx <= 1; bx++) {
        for (int by = -1; by <= 1; by++) {
            neighborBucket = hash2D(px + bx, py + by) % buckets;

        }
    }
    return neighborBucket;
}

unsigned int hash2D(int cx, int cy) {
    return (unsigned int)(cx * hk1 ^ cy * hk2);
}

void printParticleData() {
    for (int i = 0; i < particleNUM; i++) {
        printf("P(%d): ", i);

        for (int index = 0; index < pstride; index++) {
            printf(" | %f | ", particles[i * pstride + index]);
        }
        printf("\n");
    }
}

int main() {
    calcHs();
    startHash();
    spawn_particles();

    GLFWwindow* window;
    if (!glfwInit()) {
        fprintf(stderr, "Failed to initialize GLFW\n");
        return -1;
    }

    window = glfwCreateWindow(800, 800, "Fluid Sim", NULL, NULL);
    if (!window) {
        fprintf(stderr, "Failed to create GLFW window\n");
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // vsync

    // --- OpenGL 2D setup ---
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black background
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1, 1, -1, 1, -1, 1);
    glMatrixMode(GL_MODELVIEW);

    int count = 0;
    while (!glfwWindowShouldClose(window)) {
        //spatial hash
        add_gravity();
        fillCountArray();
        generatePartialSum();
        scatterSort();
        //spatial hash 


        //sim
        calcDensity();
        //calcPressure();
        calcGradient();
        //sim

        timestep();
        enforce_boundaries();
        //for (int i = 0; i < buckets + 1; i++) {
        //       printf(" | %d.) cellCount = %d", i, cellCount[i]);
        //       int x = pid[i];
        //       if (x) {
        //           printf(" | %d.) pid index = %d", i, pid[i]);
        //       }
        //       printf(" | %d.) prefixCount = %d\n", i, prefixSum[i]);
        //}


        //debugging
        //printParticleData();
        //debugging


       // --- Rendering ---
        glClear(GL_COLOR_BUFFER_BIT);
        glLoadIdentity();
        //glColor3f(1.0f, 1.0f, 1.0f); // White particles
        glPointSize(3.5f);

        // In your rendering code
        glBegin(GL_POINTS);
        for (int i = 0; i < particleNUM; i++) {
            float density = getDensity(i);
            // Map density to color (blue for low density, red for high density)
            float r = fminf(1.0f, density / (rho0));
            float b = fminf(1.0f, 1.0f - density / (rho0));
            glColor3f(b, 0.5f, 0.5f);

            float x = particles[i * pstride];
            float y = particles[i * pstride + 1];
            float nx = (x / gridX) * 2.0f - 1.0f;
            float ny = (y / gridY) * 2.0f - 1.0f;
            glVertex2f(nx, ny);
        }
        glEnd();

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;

    return 0;
}
