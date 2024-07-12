# Smart Rockets
A simple genetic algorithm simulation written in C and [Raylib](https://www.raylib.com/).

![screenshot](/screenshot.png)

## Background

### Key Components
1. **Population**: A set of potential solutions to the problem (in this case, rockets).
2. **Genes**: The basic units of information that define an individual (here, movement instructions for rockets).
3. **Fitness Function**: A measure of how well an individual solves the problem (in this simulation, how close a rocket gets to the target and how quickly).
4. **Selection**: The process of choosing individuals for reproduction based on their fitness.
5. **Crossover**: The creation of offspring by combining genetic information from two parents.
6. **Mutation**: Random changes in the genes to maintain genetic diversity.

## Simulation

### Overview
- Each rocket represents an individual in the population.
- The rocket's path is determined by a series of force vectors (genes).
- Fitness is calculated based on how close a rocket gets to the target.
- Rockets that perform better have a higher chance of passing their genes to the next generation.
- The rockets get progressively better through the generations

### Fitness Function Limitations
- The current fitness function, based on the inverse of distance to the target, has limitations. It can lead to suboptimal behavior when obstacles are positioned near the target.
- Rockets that collide with obstacles close to the target may receive high fitness scores, despite failing to reach the actual goal.
- This can result in the population converging on strategies that don't successfully navigate all obstacles.
- A more sophisticated approach, such as using a flood-fill algorithm to calculate distances through valid paths, could significantly improve learning effectiveness.
- Such an enhancement would better reward rockets that successfully navigate around obstacles, even if their final position is technically farther from the target.

### Usage
1. Install [Raylib](https://github.com/raysan5/raylib?tab=readme-ov-file#build-and-installation).
2. Run the simulation with
    ```
    eval cc main.c $(pkg-config --libs --cflags raylib) && ./a.out
    ```
3. Use `Spacebar` to control how fast the simulation runs.

## Resources
Here are some excellent resources to learn more about genetic algorithms:
- **Video Series**: [Genetic Algorithm](https://www.youtube.com/watch?v=9zfeTw-uFCw&list=PLRqwX-V7Uu6bJM3VgzjNV5YxVxUwzALHV&ab_channel=TheCodingTrain) by The Coding Train.
- **Book**: [Nature of Code](http://natureofcode.com/book/) by Daniel Shiffman.
- [Smart Rockets in p5.js on YouTube](https://www.youtube.com/watch?v=bGz7mv2vD6g)

### Comments
- Genetic algorithms (GAs) are search heuristics inspired by Charles Darwin's theory of natural evolution. They reflect the process of natural selection where the fittest individuals are selected for reproduction to produce offspring of the next generation.
- The goal of the simulation is not to get as many rockets as possible to the target but to search a space for the target.

## About Me
- A list of all my projects - [bones-ai.bearblog.dev/projects/](https://bones-ai.bearblog.dev/projects/)
- [Twitter](https://twitter.com/BonesaiDev)
