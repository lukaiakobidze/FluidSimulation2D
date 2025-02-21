using System;
using System.Collections.Generic;
using SFML.Graphics;
using SFML.System;
using SFML.Window;

class Program
{
    static void Main()
    {
        RenderWindow window = new RenderWindow(new VideoMode(800, 600), "Liquid Simulation");
        window.Closed += (sender, e) => window.Close();

        List<Particle> particles = new List<Particle>();
        Dictionary<Vector2i, List<Particle>> grid = new Dictionary<Vector2i, List<Particle>>();
        float cellSize = 2.0f;

        // Initialize particles
        for (int i = 0; i < 40; i++)
        {
            particles.Add(new Particle(i % 20 * 10, 100 + (i / 20) * 10, 5));
            //particles.Add(new Particle(i * 10, 100, 5));
        }
        particles[15].Shape.FillColor = Color.Red;


        Clock clock = new Clock();
        float targetFPS = 60.0f; // Target frames per second
        float timePerFrame = 1.0f / targetFPS; // Time per frame in seconds

        float accumulator = 0.0f;

        while (window.IsOpen)
        {
            window.DispatchEvents();
            window.Clear(Color.Black);

            float deltaTime = clock.Restart().AsSeconds();
            accumulator += deltaTime;

            // Step 1: Run simulation at target FPS
            while (accumulator >= timePerFrame)
            {
                // Update simulation here
                UpdateSimulation(grid, particles, cellSize);
                Console.WriteLine(particles[15].Velocity);
                // Decrease accumulator by time per frame
                accumulator -= timePerFrame;
            }

            // Draw the particles after update
            foreach (var particle in particles)
            {
                window.Draw(particle.Shape);
            }

            window.Display();
        }
    }

    static void UpdateSimulation(Dictionary<Vector2i, List<Particle>> grid, List<Particle> particles, float cellSize)
    {
        grid.Clear();
        foreach (Particle p in particles)
        {
            Vector2i cell = GetGridCell(p.Shape.Position, cellSize);
            if (!grid.ContainsKey(cell))
                grid[cell] = new List<Particle>();

            grid[cell].Add(p);
        }

        // Update particles
        foreach (var particle in particles)
        {
            particle.Update(grid, cellSize);
        }
    }

    static Vector2i GetGridCell(Vector2f position, float cellSize)
    {
        return new Vector2i((int)(position.X / cellSize), (int)(position.Y / cellSize));
    }
}

class Particle
{
    public CircleShape Shape;
    public Vector2f Velocity;
    public float Mass;
    public float Density;
    public float Pressure;
    public float InternalEnergy;
    public Vector2f Acceleration;
    public Vector2f gravityAc = new Vector2f(0.0f, 0.0f);
    public const float alpha = 0.5f;
    public const float beta = 0.5f;
    public const float cs = 1500;
    public Particle(float x, float y, int size)
    {
        Shape = new CircleShape(size) { FillColor = Color.Blue, Position = new Vector2f(x, y) };
        Velocity = new Vector2f(0, 0);
        Density = 0.0f;
        Mass = 1000 * (5 * 5);
        InternalEnergy = 1;
        Acceleration = new Vector2f(0.0f, 0.1f);
    }

    public void Update(Dictionary<Vector2i, List<Particle>> grid, float cellSize)
    {
        ApplyGravity();

        CalculateDensity(grid, cellSize);
        CalculatePressure();
        CalculateInternalEnergy(grid, cellSize);
        //CalculateAccelerationMomentum(grid, cellSize);
        //ApplyArtificialViscosityForNeighbors(grid, cellSize, alpha, beta, cs);
        PreventParticleOverlap(grid, cellSize);
        Velocity += Acceleration;
        Velocity += gravityAc;
        Shape.Position += Velocity;
    }
    private void PreventParticleOverlap(Dictionary<Vector2i, List<Particle>> grid, float h)
    {
        // Get the grid cell for the current particle
        Vector2i myCell = GetGridCell(Shape.Position, h);

        // Get neighboring grid cells
        List<Vector2i> neighborCells = GetNeighborCells(myCell);

        // Loop through each neighboring cell
        foreach (Vector2i cell in neighborCells)
        {
            if (!grid.ContainsKey(cell)) continue;

            // Loop through all particles in the neighboring cell
            foreach (Particle particle in grid[cell])
            {
                // Skip if it's the current particle
                if (particle == this) continue;

                // Calculate distance between the particles
                float distance = MathF.Sqrt(MathF.Pow(this.Shape.Position.X - particle.Shape.Position.X, 2) +
                                            MathF.Pow(this.Shape.Position.Y - particle.Shape.Position.Y, 2));

                // If the distance is less than the sum of their radii, they are overlapping
                float minDistance = this.Shape.Radius + particle.Shape.Radius;

                if (distance < minDistance)
                {
                    // Calculate the overlap distance
                    float overlap = minDistance - distance;

                    // Calculate the direction of the separation vector
                    Vector2f direction = this.Shape.Position - particle.Shape.Position;
                    float length = MathF.Sqrt(direction.X * direction.X + direction.Y * direction.Y);
                    if (length > 0) // Avoid division by zero
                        direction /= length;

                    // Apply the correction to both particles
                    Vector2f correction = direction * overlap / 2.0f;

                    // Apply the position correction
                    this.Shape.Position += correction;
                    particle.Shape.Position -= correction;

                    // Apply velocity adjustment to prevent them from getting stuck together
                    Vector2f relativeVelocity = this.Velocity - particle.Velocity;
                    float velocityDot = DotProduct(relativeVelocity, direction);

                    if (velocityDot < 0) // Only apply if the particles are moving towards each other
                    {
                        // Apply a repulsive force based on relative velocity to separate the particles
                        float restitution = 0.5f; // Restitution coefficient, determines how "bouncy" the collision is
                        Vector2f repulsiveForce = restitution * velocityDot * direction;

                        this.Velocity -= repulsiveForce;
                        particle.Velocity += repulsiveForce;
                    }
                }
            }
        }
    }


    private void ApplyGravity()
    {
        // Gravity acceleration (pointing downwards)
        const float gravity = 0.1f;

        // Only apply gravity if particle is not on the ground
        if (Shape.Position.Y < 590)  // Modify this condition based on your ground Y position
        {
            // Apply gravity to the Y-component of acceleration, only if not touching the ground
            gravityAc.Y = gravity;
        }
        else
        {
            // If the particle is touching the ground, stop gravity and make sure it's not falling below
            if (Shape.Position.Y > 590) // Check if particle is below the ground level
            {
                Shape.Position = new Vector2f(Shape.Position.X, 590);  // Correct particle position to the ground level
                Velocity.Y = -Velocity.Y / 5;  // Stop downward velocity once on the ground
                Acceleration.Y = 0.0f;  // Stop gravity effect
                gravityAc.Y = 0;
            }
        }
    }

    public void ApplyArtificialViscosityForNeighbors(Dictionary<Vector2i, List<Particle>> grid, float h, float alpha, float beta, float cs)
    {
        Vector2f viscosityForce = new Vector2f(0.0f, 0.0f);

        Vector2i myCell = GetGridCell(Shape.Position, h);
        List<Vector2i> neighborCells = GetNeighborCells(myCell);

        foreach (Vector2i cell in neighborCells)
        {
            if (!grid.ContainsKey(cell)) continue;

            // Loop over particles in the neighboring cells
            foreach (Particle i in grid[cell])
            {
                if (i == this) continue;

                float r = MathF.Sqrt(MathF.Pow(this.Shape.Position.X - i.Shape.Position.X, 2) +
                                     MathF.Pow(this.Shape.Position.Y - i.Shape.Position.Y, 2));

                if (r < 0.0001f) continue; // Avoid division by zero

                // Calculate mu (relative velocity and distance factor)
                float mu = (h * DotProduct(this.Velocity - i.Velocity, this.Shape.Position - i.Shape.Position)) /
                           (r * r + 0.01f * h * h);

                // Calculate viscosity term
                float viscosity = (-alpha * cs * mu + beta * mu * mu) / this.Density;

                // Calculate the viscosity force and add it to the total force on this particle
                viscosityForce -= i.Mass * viscosity * KernelGradient(this.Shape.Position, i.Shape.Position, h);
            }
        }

        // Apply the viscosity force to the particle's acceleration
        this.Acceleration += viscosityForce;
    }

    public void CalculateDensity(Dictionary<Vector2i, List<Particle>> grid, float h)
    {
        Density = 0;
        Vector2i myCell = GetGridCell(Shape.Position, h);

        List<Vector2i> neighborCells = GetNeighborCells(myCell);
        foreach (var cell in neighborCells)
        {
            if (!grid.ContainsKey(cell)) continue;

            foreach (Particle i in grid[cell])
            {
                if (i == this) continue;

                float r = MathF.Sqrt(MathF.Pow(i.Shape.Position.X - Shape.Position.X, 2) +
                                     MathF.Pow(i.Shape.Position.Y - Shape.Position.Y, 2));

                Density += i.Mass * KernelFunction(r, h);
            }
        }
    }

    public void CalculatePressure()
    {
        Pressure = 6 * InternalEnergy * Density;
    }
    public void CalculateInternalEnergy(Dictionary<Vector2i, List<Particle>> grid, float h)
    {
        float energyChange = 0;
        Vector2i myCell = GetGridCell(Shape.Position, h);
        List<Vector2i> neighborCells = GetNeighborCells(myCell);
        foreach (var cell in neighborCells)
        {
            if (!grid.ContainsKey(cell)) continue;

            foreach (Particle i in grid[cell])
            {
                if (i == this) { continue; }

                float r = MathF.Sqrt(MathF.Pow(i.Shape.Position.X - Shape.Position.X, 2) +
                                     MathF.Pow(i.Shape.Position.Y - Shape.Position.Y, 2));
                Vector2f gradW = KernelGradient(Shape.Position, i.Shape.Position, h);
                float velocityDotGradW = DotProduct(Velocity - i.Velocity, gradW);

                float pressureTerm = (Pressure / (Density * Density)) + (i.Pressure / (i.Density * i.Density));

                energyChange += i.Mass * pressureTerm * velocityDotGradW;

            }
        }
    }
    public void CalculateAccelerationMomentum(Dictionary<Vector2i, List<Particle>> grid, float h)
    {
        Vector2f pressureForce = new Vector2f(0.0f, 0.0f);
        Vector2i myCell = GetGridCell(Shape.Position, h);
        List<Vector2i> neighborCells = GetNeighborCells(myCell);
        foreach (var cell in neighborCells)
        {
            if (!grid.ContainsKey(cell)) continue;

            foreach (Particle i in grid[cell])
            {
                if (i == this) continue;

                float r = MathF.Sqrt(MathF.Pow(i.Shape.Position.X - Shape.Position.X, 2) +
                                     MathF.Pow(i.Shape.Position.Y - Shape.Position.Y, 2));
                Vector2f gradW = KernelGradient(Shape.Position, i.Shape.Position, h);

                float pressureTerm = (Pressure / (Density * Density)) + (i.Pressure / (i.Density * i.Density));

                pressureForce -= i.Mass * pressureTerm * gradW;
            }
        }

        Acceleration += pressureForce;
    }

    public float DotProduct(Vector2f v1, Vector2f v2)
    {
        return v1.X * v2.X + v1.Y * v2.Y;
    }
    public Vector2f KernelGradient(Vector2f p1, Vector2f p2, float h)
    {
        
        Vector2f direction = p2 - p1;
        float r = MathF.Sqrt(MathF.Pow(p2.X - p1.X, 2) + MathF.Pow(p2.Y - p1.Y, 2));

        if (r >= h || r == 0) return new Vector2f(0, 0); 

        direction /= r;

        float gradientMagnitude = -945.0f / (64.0f * MathF.PI * MathF.Pow(h, 9)) * r * MathF.Pow(h * h - r * r, 2);

        return direction * gradientMagnitude;
    }
    public float KernelFunction(float r, float h)
    {
        if (r >= h) return 0;

        float term = (h * h) - (r * r);
        return (315.0f / (64.0f * MathF.PI * MathF.Pow(h, 9))) * MathF.Pow(term, 3);
    }

    static Vector2i GetGridCell(Vector2f position, float cellSize)
    {
        return new Vector2i((int)(position.X / cellSize), (int)(position.Y / cellSize));
    }

    static List<Vector2i> GetNeighborCells(Vector2i cell)
    {
        return new List<Vector2i>
        {
            cell,
            new Vector2i(cell.X - 1, cell.Y), new Vector2i(cell.X + 1, cell.Y),
            new Vector2i(cell.X, cell.Y - 1), new Vector2i(cell.X, cell.Y + 1),
            new Vector2i(cell.X - 1, cell.Y - 1), new Vector2i(cell.X + 1, cell.Y + 1),
            new Vector2i(cell.X - 1, cell.Y + 1), new Vector2i(cell.X + 1, cell.Y - 1)
        };
    }
}
