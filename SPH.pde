/*
(Partial) implementation of Clavet et al (2005): "Particle-based Viscoelastic Fluid Simulation"
See https://hmbd.wordpress.com/2019/06/10/sph-fluid-simulation-in-processing/ for further documentation and animations.
Press Mouse to add fluid, change int model in line 11 to see different models.
*/

Fluid fluid;
float kernel_radius = 20; // rayon d'influence des particules
float timestep = 1;
int model = 1; // 1: dam break + lid-driven cavity, 2: obstacles, 3: hélice
float stepsPerFrame = 10;

float barAngle = 0; // Angle initial de l'hélice
float barAngularSpeed = 0.03; // Vitesse angulaire (rad/frame)
float barRadius = width * 0.8; // Longueur de l'hélice

void setup() {
  size(350, 350);
  colorMode(HSB);

  fluid = new Fluid(timestep, kernel_radius);

  // model 1 : dam break + lid-driven cavity
  if (model == 1) {
    stepsPerFrame = 1;
    timestep = 1;
    for (float x = 0; x<width; x+=kernel_radius*0.25) {
      for (float y = height/2; y<height; y+=kernel_radius*0.25) {
        Particle p = new Particle(x, y, false);
        p.vel = new PVector(random(-1, 1), random(-1, 1)); // Vitesse initiale aléatoire
        fluid.particles.add(p);
      }
    }
    createBorders();
  }

  // model 2 : obstacles et fluides opposés
  if (model == 2) { 
    stepsPerFrame = 1;
    timestep = 1;

    // Fluide gauche
    for (float x = 0; x < width*0.3; x += kernel_radius * 0.25) {
      for (float y = height * 0; y < height; y += kernel_radius * 0.25) {
        Particle p = new Particle(x, y, false);
        p.vel = new PVector(1, 0); // Vitesse initiale
        fluid.particles.add(p);
      }
    }

    // Obstacles
    createBorders();
    for (float y = height * 0.70; y < height * 0.80; y += kernel_radius * 0.1) {
      fluid.particles.add(new Particle(width * 0.5, y, true)); // Obstacle central vertical haut
    }
    for (float y = height * 0.85; y < height * 0.95; y += kernel_radius * 0.1) {
      fluid.particles.add(new Particle(width * 0.5, y, true)); // Obstacle central vertical bas
    }
    for (float x = width * 0.1; x < width * 0.2; x += kernel_radius * 0.1) {
      fluid.particles.add(new Particle(x, height * 0.85, true)); // Obstacle gauche horizontal
    }
    for (float x = width * 0.8; x < width * 0.9; x += kernel_radius * 0.1) {
      fluid.particles.add(new Particle(x, height * 0.85, true)); // Obstacle droit horizontal
    }

    // Fluide droit opposé
    for (float x = width * 0.7; x < width; x += kernel_radius * 0.25) {
      for (float y = height * 0.5; y < height; y += kernel_radius * 0.25) {
        Particle p = new Particle(x, y, false);
        p.vel = new PVector(-1, 0); // Vitesse initiale
        fluid.particles.add(p);
      }
    }
  }

  // model 3 : hélice
  if (model == 3) {
    stepsPerFrame = 1;
    timestep = 1;
    for (float x = 0; x<width; x+=kernel_radius*0.25) {
      for (float y = height/2; y<height; y+=kernel_radius*0.25) {
        Particle p = new Particle(x, y, false);
        p.vel = new PVector(random(-1, 1), random(-1, 1)); // Vitesse initiale aléatoire
        fluid.particles.add(p);
      }
    }
    createBorders();

    // Initialisation de l'hélice (2 lignes de particules)
    for (float y = -barRadius; y <= barRadius; y += kernel_radius * 0.1) {
      float x = width * 0.5;
      float posY = height * 0.8 + y;
      Particle barParticle = new Particle(x, posY, true);
      barParticle.isPartOfHelix = true;
      fluid.particles.add(barParticle);
    }
    for (float x = -barRadius; x <= barRadius; x += kernel_radius * 0.1) {
      float posX = width * 0.5 + x;
      float y = height * 0.8;
      Particle barParticle = new Particle(posX, y, true);
      barParticle.isPartOfHelix = true;
      fluid.particles.add(barParticle);
    }
  }

  print("Particles: ");
  println(fluid.particles.size());
}

void draw() {
  background(170);

  if (mousePressed) {
    for (int i=1; i<5; i++) {
      fluid.particles.add(new Particle(mouseX + random(kernel_radius), mouseY + random(kernel_radius), false));
    }
  }

  for (int i=0; i<stepsPerFrame; i++) {
    fluid.step();
  }

  fluid.plot();
}

void createBorders() {
  for (float x=0; x < width; x+= kernel_radius * 0.1) {
    fluid.particles.add(new Particle(x, height, true));
  }
  for (float y=0; y < height; y+= kernel_radius * 0.1) {
    fluid.particles.add(new Particle(0, y, true));
    fluid.particles.add(new Particle(width, y, true));
  }
}
