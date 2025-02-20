using System;
using System.Collections.Generic;
using SFML.Graphics;
using SFML.System;
using SFML.Window;


class Program
{
    static void Main()
    {
        RenderWindow window = new RenderWindow(new VideoMode(800, 600), "Liquid Simlutaion");
        window.Closed += (sender, e) => window.Close();

        List<Particle> particles = new List<Particle>();

        for (int i = 0; i < 10; i++)
        {
            particles.Add(new Particle(i * 20, 100, 5));
        }

        while (window.IsOpen)
        {
            window.DispatchEvents();
            window.Clear(Color.Black);

            foreach (var particle in particles)
            {
                particle.Update();
                window.Draw(particle.Shape);
            }
            window.Display();
        }
    }
}

class Particle
{
    public CircleShape Shape;
    public Vector2f Velocity;

    public Particle(float x, float y, int size)
    {
        Shape = new CircleShape(size) { FillColor = Color.Blue, Position = new Vector2f(x, y) };

        Velocity = new Vector2f(0, 0);
    }

    public void Update()
    {
        Velocity += new Vector2f(0, 0.0001f);
        Shape.Position += Velocity;

        if (Shape.Position.Y > 580)
        {
            Shape.Position = new Vector2f(Shape.Position.X, 580);
            Velocity = new Vector2f(0, 0);
        }
    }
}