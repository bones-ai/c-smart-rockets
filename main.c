/* 
 * =================================================================
 * 
 * Smart Rockets:
 * A simple genetic algorithm simulation
 *
 * Description:
 * This program simulates a population of "rockets" that evolve
 * over generations to find an optimal path to a target. It uses a
 * genetic algorithm to evolve the rocket's behavior.
 *
 * Libraries Used:
 * Raylib - For graphics rendering - https://www.raylib.com/
 *
 * Links:
 * Repo - https://github.com/bones-ai/c-smart-rockets
 * More Simulations - https://github.com/bones-ai
 *
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "raylib.h"

// =================================================================
// MARK: Configs
// =================================================================

#define WINDOW_W 1400
#define WINDOW_H 900
#define SIM_FPS 120
#define NUM_ROCKETS 5000
#define MAX_FRAMES 250
#define MUTATION_PROBABILITY 0.1f
#define MUTATION_MAGNITUDE 0.1f
#define POP_RANDOM 0.15f
#define BORDERS_SIZE 20
#define ROCKET_SIZE 10.0f
#define TARGET_RADIUS 20.0f

// Spawns
const Vector2 ROCKET_SPAWN_LOC = {BORDERS_SIZE * 2, WINDOW_H / 2};
const Vector2 TARGET_LOC = {WINDOW_W - BORDERS_SIZE * 2, WINDOW_H / 2};

// Colors
const Color COLOR_ROCKET = DARKBLUE;
const Color COLOR_ROCKET_DEAD = {194, 194, 194, 115};
const Color COLOR_ROCKET_COMPLETE = PINK;
const Color COLOR_ROCKET_BORDER = WHITE;
const Color COLOR_OBSTACLE = GRAY;
const Color COLOR_OBSTACLE_BORDER = BLACK;
const Color COLOR_TARGET = DARKGREEN;
const Color COLOR_TARGET_BORDER = BLACK;

// =================================================================
// MARK: Types
// =================================================================

typedef enum
{
    ALIVE,
    DEAD,
    COMPLETE
} RocketState;

typedef enum
{
    SUCCESS,
    ERR_GENE_POOL_INIT_BAD,
    ERR_POP_INIT_BAD,
    ERR_OBSTACLES_INIT_BAD,
    ERR_BAD_RANDOM_POP_CONFIG,
} StatusCode;

typedef struct
{
    Vector2 pos;
    Vector2 vel;
    Vector2 acc;

    RocketState state;
    float fitness;
    Vector2 genes[MAX_FRAMES];
} Rocket;

typedef struct
{
    Vector2 topLeft;
    Vector2 bottomRight;
} Obstacle;

typedef struct
{
    int gen;
    int frameIdx;
    size_t numObstacles;
    bool isSlowMode;
    Rocket *agents;
    Obstacle *obstacles;
} Sim;

typedef struct {
    int* pool;
    int size;
} GenePool;

// =================================================================
// MARK: Prototypes
// =================================================================

float getgetRandFloat(void);
float getDistanceSquared(const Vector2 *first, const Vector2 *second);
bool isInRange(const Vector2 *first, const Vector2 *second, float radius);
Vector2 rotatePoint(Vector2 *point, const Vector2 *center, float angle);

Rocket createRocket(void);
void rocketApplyForce(Rocket *r, Vector2 *force);
float getFitness(const Rocket *r, int frameIdx);
void updateRocket(Rocket *r, Obstacle obstacles[], size_t numObstacles, int frameIdx);
void drawRocket(const Rocket *r);
Rocket crossover(const Rocket *first, const Rocket *second);
void mutate(Rocket *r);

Obstacle createObstacle(float topLeftX, float topLeftY, float bottomRightX, float bottomRightY);

Sim initSim(Obstacle obstacles[], size_t numObstacles);
void updateSim(Sim *sim);
void drawSim(const Sim *sim);
GenePool createGenePool(Sim *sim);
void startNewGeneration(Sim *sim);
void endCurrentGeneration(Sim *sim);
void terminateSim(Sim *sim, StatusCode e);

// =================================================================
// MARK: Utils
// =================================================================

// Return a random float between -1.0 and 1.0
float getRandFloat(void)
{
    float val = (float)rand() / (float)RAND_MAX;
    return (rand() % 2) ? val : -val;
}

// Return the squared distance between two points
float getDistanceSquared(const Vector2 *first, const Vector2 *second)
{
    int dx = first->x - second->x;
    int dy = first->y - second->y;
    return dx * dx + dy * dy;
}

// Return true if the distance between two points is less than a radius
bool isInRange(const Vector2 *first, const Vector2 *second, float radius)
{
    return getDistanceSquared(first, second) <= (radius * radius);
}

// Rotate a point around center by an angle
Vector2 rotatePoint(Vector2 *point, const Vector2 *center, float angle)
{
    Vector2 rotated;
    float s = sinf(angle);
    float c = cosf(angle);

    // Translate point back to origin
    point->x -= center->x;
    point->y -= center->y;

    // Rotate point
    float xnew = point->x * c - point->y * s;
    float ynew = point->x * s + point->y * c;

    // Translate point back
    rotated.x = xnew + center->x;
    rotated.y = ynew + center->y;

    return rotated;
}

// =================================================================
// MARK: Rocket Impl
// =================================================================

// Create a new rocket with random genes
Rocket createRocket(void)
{
    Rocket rocket;
    rocket.pos = ROCKET_SPAWN_LOC;
    rocket.vel = (Vector2){0, 0};
    rocket.acc = (Vector2){0, 0};
    rocket.state = ALIVE;
    rocket.fitness = 0.0f;

    for (int i = 0; i < MAX_FRAMES; i++) {
        rocket.genes[i] = (Vector2){ getRandFloat(), getRandFloat() };
    }

    return rocket;
}

// Apply a force to transform the rocket's position
void rocketApplyForce(Rocket *r, Vector2* force)
{
    r->acc.x += force->x;
    r->acc.y += force->y;

    r->vel.x += r->acc.x;
    r->vel.y += r->acc.y;

    r->pos.x += r->vel.x;
    r->pos.y += r->vel.y;

    r->acc = (Vector2) { 0, 0 };
}

// Compute the fitness of a rocket
// Fitness is a measure of how well the rocket performed
// This isn't a perfect fitness function
// It fails when there's a obstacle just before the target
// A flood fill algorithm would yield better results
float getFitness(const Rocket *r, int frameIdx)
{
    float distToTarget = getDistanceSquared(&r->pos, &TARGET_LOC);
    return 1.0f / fmaxf(distToTarget, 1.0f);
}

// Update a rocket's position and state
void updateRocket(Rocket *r, Obstacle obstacles[], size_t numObstacles, int frameIdx)
{
    // Validate rocket state
    if (r->state != ALIVE || frameIdx >= MAX_FRAMES)
    {
        return;
    }

    // Update fitness
    r->fitness = getFitness(r, frameIdx);

    // Obstacle Collision
    for (int i = 0; i < numObstacles; i ++)
    {
        bool xHit = r->pos.x > obstacles[i].topLeft.x && r->pos.x < obstacles[i].bottomRight.x;
        bool yHit = r->pos.y > obstacles[i].topLeft.y && r->pos.y < obstacles[i].bottomRight.y;
        if (xHit && yHit)
        {
            r->state = DEAD;
            return;
        }
    }

    // Target Collision
    if (isInRange(&r->pos, &TARGET_LOC, TARGET_RADIUS))
    {
        r->state = COMPLETE;
        return;
    }

    // Translate the rocket based on gene
    rocketApplyForce(r, &r->genes[frameIdx]);
}

// Draw a rocket on the screen
// Rockets are drawn as triangles
void drawRocket(const Rocket *r)
{
    // Define the triangle's size
    float height = ROCKET_SIZE * 2.0f;
    float base = ROCKET_SIZE * 1.5f;

    // Calculate the three points of the triangle
    Vector2 top = {r->pos.x, r->pos.y - height/2};
    Vector2 bottomLeft = {r->pos.x - base/2, r->pos.y + height/2};
    Vector2 bottomRight = {r->pos.x + base/2, r->pos.y + height/2};

    // Calculate the rotation angle based on the rocket's velocity
    float rotation = atan2f(r->vel.y, r->vel.x) + PI/2;
    Vector2 rotatedTop = rotatePoint(&top, &r->pos, rotation);
    Vector2 rotatedBottomLeft = rotatePoint(&bottomLeft, &r->pos, rotation);
    Vector2 rotatedBottomRight = rotatePoint(&bottomRight, &r->pos, rotation);

    // Draw the triangle
    Color fillColor;
    switch (r->state) {
        case COMPLETE:
            fillColor = COLOR_ROCKET_COMPLETE;
            break;
        case DEAD:
            fillColor = COLOR_ROCKET_DEAD;
            break;
        default:
            fillColor = COLOR_ROCKET;
            break;
    }
    DrawTriangle(rotatedTop, rotatedBottomLeft, rotatedBottomRight, fillColor);
}

// Create a new rocket by merging the genes of two parent rockets
// The split point is randomly chosen
// The child rocket's genes are a mix of the parent's genes with a chance of mutation
Rocket crossover(const Rocket *first, const Rocket *second)
{
    Rocket child;
    int splitPoint = rand() % MAX_FRAMES;

    child.pos = ROCKET_SPAWN_LOC;
    child.vel = (Vector2){0, 0};
    child.acc = (Vector2){0, 0};
    child.state = ALIVE;
    child.fitness = 0.0f;

    // memcpy(dest, src, size);
    // Copy the first half of the first parent's genes
    memcpy(child.genes, first->genes, splitPoint * sizeof(Vector2));
    // Copy the second half of the second parent's genes
    memcpy(&child.genes[splitPoint], &second->genes[splitPoint], (MAX_FRAMES - splitPoint) * sizeof(Vector2));

    return child;
}

// Mutate a rocket's genes
// Each gene has a chance of being mutated
void mutate(Rocket *r)
{
    // Mutate the rocket genes
    for (int i = 0; i < MAX_FRAMES; i ++)
    {
        float mutataionChance = (float)rand() / (float)RAND_MAX;
        if (mutataionChance <= MUTATION_PROBABILITY)
        {
            // Alter the gene by a random amount
            r->genes[i].x += getRandFloat() * MUTATION_MAGNITUDE;
            r->genes[i].y += getRandFloat() * MUTATION_MAGNITUDE;
        }
    }
}

// =================================================================
// MARK: Obstacle Impl
// =================================================================

// Create a new obstacle
// This doesn't do any validation
Obstacle createObstacle(float topLeftX, float topLeftY, float bottomRightX, float bottomRightY)
{
    return (Obstacle){ 
        (Vector2){ topLeftX, topLeftY }, 
        (Vector2){ bottomRightX, bottomRightY } 
    };
}

// =================================================================
// MARK: Sim Impl
// =================================================================

// Create a gene pool based on the fitness of the rockets
GenePool createGenePool(Sim *sim)
{
    // Normalize fitness values to be between 0.0 to 1.0
    // This ensures that the sum of all fitness values is 1.0
    float totalFitness = 0.0f;
    for (int i = 0; i < NUM_ROCKETS; i++)
    {
        totalFitness += sim->agents[i].fitness;
    }
    for (int i = 0; i < NUM_ROCKETS; i++)
    {
        sim->agents[i].fitness = sim->agents[i].fitness / totalFitness;
    }

    // Calculate total number of slots needed
    size_t totalSlots = 0;
    for (int i = 0; i < NUM_ROCKETS; i++) {
        // Add 0.5f to round up
        totalSlots += (int)(sim->agents[i].fitness * 1000 + 0.5f);
    }

    // Allocate gene pool dynamically
    int *genepool = malloc(totalSlots * sizeof(int));
    if (genepool == NULL) {
        terminateSim(sim, ERR_GENE_POOL_INIT_BAD);
    }

    // Fill the gene pool with the index of the rockets
    // This results in a weighted random selection of the rockets
    // This is also known as Roulette Wheel Selection
    int index = 0;
    for (int i = 0; i < NUM_ROCKETS; i++) {
        int times = (int)(sim->agents[i].fitness * 1000 + 0.5f);
        for (int j = 0; j < times; j++) {
            if (index < totalSlots) {
                genepool[index++] = i;
            } else {
                TraceLog(LOG_INFO, "Wrong genepool total slots calc");
                break;
            }
        }
    }

    return (GenePool) { genepool, totalSlots };
}

// Initialize the simulation with a population of rockets
Sim initSim(Obstacle obstacles[], size_t numObstacles)
{
    Sim sim;
    sim.gen = 0;
    sim.frameIdx = 0;
    sim.isSlowMode = true;
    sim.numObstacles = numObstacles;

    sim.agents = malloc(NUM_ROCKETS * sizeof(Rocket));
    if (sim.agents == NULL) {
        terminateSim(&sim, ERR_POP_INIT_BAD);
    }

    sim.obstacles = malloc(numObstacles * sizeof(Obstacle));
    if (sim.obstacles == NULL) {
        terminateSim(&sim, ERR_OBSTACLES_INIT_BAD);
    }
    for (int i = 0;i < numObstacles; i ++)
    {
        sim.obstacles[i] = obstacles[i];
    }

    for (int i = 0; i < NUM_ROCKETS; i++) {
        sim.agents[i] = createRocket();
    }

    return sim;
}

// Draw the rockets, target, and obstacles
void drawSim(const Sim *sim)
{
    // Target
    DrawCircle(TARGET_LOC.x, TARGET_LOC.y, TARGET_RADIUS, COLOR_TARGET);
    DrawCircleLines(TARGET_LOC.x, TARGET_LOC.y, TARGET_RADIUS, COLOR_TARGET_BORDER);

    // Obstacles
    for (int i = 0; i < sim->numObstacles; i++)
    {
        DrawRectangle(
            sim->obstacles[i].topLeft.x,
            sim->obstacles[i].topLeft.y,
            sim->obstacles[i].bottomRight.x - sim->obstacles[i].topLeft.x,
            sim->obstacles[i].bottomRight.y - sim->obstacles[i].topLeft.y,
            COLOR_OBSTACLE
        );
        DrawRectangleLinesEx(
            (Rectangle){
                sim->obstacles[i].topLeft.x,
                sim->obstacles[i].topLeft.y,
                sim->obstacles[i].bottomRight.x - sim->obstacles[i].topLeft.x,
                sim->obstacles[i].bottomRight.y - sim->obstacles[i].topLeft.y,
            },
            3,
            COLOR_OBSTACLE_BORDER
        );
    }

    // Rockets
    for (int i = 0; i < NUM_ROCKETS; i++)
    {
        drawRocket(&sim->agents[i]);
    }
}

// Start a new generation
// Creates a new population of rockets based on the current population
void startNewGeneration(Sim *sim)
{
    GenePool genepool = createGenePool(sim);
    Rocket *children = malloc(NUM_ROCKETS * sizeof(Rocket));
    int numRandom = NUM_ROCKETS * POP_RANDOM;
    if (numRandom >= NUM_ROCKETS)
    {
        // This is a bad config
        terminateSim(sim, ERR_BAD_RANDOM_POP_CONFIG);
    }

    // Crossover and Mutate
    for (int i = 0; i < (NUM_ROCKETS - numRandom); i ++)
    {
        int parent1 = genepool.pool[rand() % genepool.size];
        int parent2 = genepool.pool[rand() % genepool.size];
        children[i] = crossover(&sim->agents[parent1], &sim->agents[parent2]);
        mutate(&children[i]);
    }

    // Random Pop
    // This helps maintain some diversity in the population
    int i = (NUM_ROCKETS - numRandom);
    while(i < NUM_ROCKETS)
    {
        children[i] = createRocket();
        i ++;
    }

    free(sim->agents);
    sim->agents = children;

    // We don't need this genepool anymore
    free(genepool.pool);
}

// Marks the end of the current generation
void endCurrentGeneration(Sim *sim)
{
    sim->gen += 1;
    sim->frameIdx = 0;
}

// Simulation tick
void updateSim(Sim *sim)
{
    sim->frameIdx += 1;
    if (sim->frameIdx >= MAX_FRAMES)
    {
        endCurrentGeneration(sim);
        startNewGeneration(sim);
    }

    for (int i = 0; i < NUM_ROCKETS; i++)
    {
        updateRocket(&sim->agents[i], sim->obstacles, sim->numObstacles, sim->frameIdx);
    }
}

// End the simulation
void terminateSim(Sim *sim, StatusCode e)
{
    if (sim->agents) { free(sim->agents); }
    if (sim->obstacles) { free(sim->obstacles); }

    TraceLog(LOG_INFO, "Terminating. Status: %d\n", e);
    exit(e == SUCCESS ? 0 : 1);
}

// =================================================================
// MARK: Main
// =================================================================

int main(void)
{
    // Seed rand
    srand((unsigned int)time(NULL));

    InitWindow(WINDOW_W, WINDOW_H, "rockets");
    SetTargetFPS(SIM_FPS);

    // Simulation Init
    const int obstacleWidth = 30;
    Obstacle obstacles[] = {
        // Left Wall
        createObstacle(-100, -100, BORDERS_SIZE / 2, WINDOW_H + 100),
        // Right Wall
        createObstacle(WINDOW_W - BORDERS_SIZE / 2, -100, WINDOW_W + 100, WINDOW_H + 100),
        // Top Wall
        createObstacle(-100, -100, WINDOW_W + 100, BORDERS_SIZE / 2),
        // Bottom Wall
        createObstacle(-100, WINDOW_H - BORDERS_SIZE / 2, WINDOW_W + 100, WINDOW_H + 100),

        // 1st Obstacle - 2 slits, 3 blocks
        createObstacle(WINDOW_W * 0.3 - obstacleWidth, BORDERS_SIZE / 2 - 3, 
            WINDOW_W * 0.3 + obstacleWidth, WINDOW_H * 0.15),
        createObstacle(WINDOW_W * 0.3 - obstacleWidth, WINDOW_H * 0.25, 
            WINDOW_W * 0.3 + obstacleWidth, WINDOW_H * 0.75),
        createObstacle(WINDOW_W * 0.3 - obstacleWidth, WINDOW_H * 0.85, 
            WINDOW_W * 0.3 + obstacleWidth, WINDOW_H - BORDERS_SIZE / 2 + 3),

        // 2nd Obstacle - 1 slit, 2 blocks
        createObstacle(WINDOW_W * 0.55 - obstacleWidth, BORDERS_SIZE / 2 - 3, 
            WINDOW_W * 0.55 + obstacleWidth, WINDOW_H * 0.5 - 50),
        createObstacle(WINDOW_W * 0.55 - obstacleWidth, WINDOW_H * 0.5 + 50, 
            WINDOW_W * 0.55 + obstacleWidth, WINDOW_H - BORDERS_SIZE / 2 + 3),

        // Final obstacle
        createObstacle(WINDOW_W * 0.75 - obstacleWidth, WINDOW_H * 0.3, 
            WINDOW_W * 0.75 + obstacleWidth, WINDOW_H * 0.7)
    };
    Sim sim = initSim(obstacles, sizeof(obstacles)/sizeof(obstacles[0]));

    while (!WindowShouldClose())
    {
        // Update
        updateSim(&sim);

        // Speed Control
        if (IsKeyPressed(KEY_SPACE))
        {
            sim.isSlowMode = !sim.isSlowMode;
            if (sim.isSlowMode)
            {
                SetTargetFPS(SIM_FPS);
            }
            else
            {
                SetTargetFPS(0);
            }
        }

        // Draw
        BeginDrawing();
            ClearBackground(RAYWHITE);

            drawSim(&sim);

            // Stats
            DrawText(
                TextFormat(
                    "GEN: %d\n\n\nFPS: %d\n\n\nSLOW: %s",
                    sim.gen, GetFPS(),
                    (sim.isSlowMode) ? "True" : "False"
                ),
                WINDOW_W * 0.75, 50, 40, DARKGRAY
            );
        EndDrawing();
    }

    terminateSim(&sim, SUCCESS);
    CloseWindow();
    return 0;
}