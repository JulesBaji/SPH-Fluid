/*
(Partial) implementation of Clavet et al (2005): "Particle-based Viscoelastic Fluid Simulation"
See https://hmbd.wordpress.com/2019/06/10/sph-fluid-simulation-in-processing/ for further documentation and animations.
Press Mouse to add fluid, change int model in line 11 to see different models.
*/

Fluid fluid;
float kernel_radius = 20; // rayon d'influence des particules
float timestep = 1;
int model = 1; // 1 - dam break, 2 - tube to fill
float stepsPerFrame = 10;

float barAngle = 0; // Angle initial de la barre
float barAngularSpeed = 0.04; // Vitesse angulaire (rad/frame)
float barRadius = width * 0.8; // Longueur de la demi-barre

void setup() {
  size(500, 500);
  // fullScreen();
  colorMode(HSB);

  fluid = new Fluid(timestep, kernel_radius);

  // model 1: dam break:
  if (model == 1) {
    stepsPerFrame = 1;
    timestep = 1;
    for (float x = 0; x<width; x+=kernel_radius*0.25) {
      for (float y = height/2; y<height; y+=kernel_radius*0.25) {
        Particle p = new Particle(x, y, false);
        p.vel = new PVector(random(-1, 1), random(-1, 1));
        fluid.particles.add(p);
      }
    }
    createBorders();

    // Initialisation de la barre rigide
    for (float y = -barRadius; y <= barRadius; y += kernel_radius * 0.1) {
      float x = width * 0.5; // Position initiale au centre
      float posY = height * 0.8 + y; // Position le long de la barre
      Particle barParticle = new Particle(x, posY, true);
      barParticle.isPartOfBar = true;
      fluid.particles.add(barParticle);
    }
    for (float x = -barRadius; x <= barRadius; x += kernel_radius * 0.1) {
      float posX = width * 0.5 + x; // Position initiale au centre
      float y = height * 0.8; // Position le long de la barre
      Particle barParticle = new Particle(posX, y, true);
      barParticle.isPartOfBar = true;
      fluid.particles.add(barParticle);
    }
  }

  if (model == 2) { // Modèle pour visualiser les tourbillons
    stepsPerFrame = 1;
    timestep = 1;

    // Initialisation du fluide : moitié gauche du canvas
    for (float x = 0; x < width*0.3; x += kernel_radius * 0.25) {
      for (float y = height * 0; y < height; y += kernel_radius * 0.25) {
        Particle p = new Particle(x, y, false);
        p.vel = new PVector(1, 0); // Ajouter un mouvement initial horizontal
        fluid.particles.add(p);
      }
    }

    // Obstacles rigides pour provoquer des perturbations
    createBorders();
    for (float y = height * 0.9; y < height * 0.98; y += kernel_radius * 0.1) {
      fluid.particles.add(new Particle(width * 0.5, y, true)); // Obstacle central vertical
    }
    for (float y = height * 0.75; y < height * 0.85; y += kernel_radius * 0.1) {
      fluid.particles.add(new Particle(width * 0.5, y, true)); // Obstacle central vertical
    }

    // Flux opposé pour créer des cisaillements
    for (float x = width * 0.7; x < width; x += kernel_radius * 0.25) {
      for (float y = height * 0.5; y < height; y += kernel_radius * 0.25) {
        Particle p = new Particle(x, y, false);
        p.vel = new PVector(-1, 0); // Mouvement initial horizontal opposé
        fluid.particles.add(p);
      }
    }
  }

  print("Particles: ");
  println(fluid.particles.size());
}

void draw() {
  background(255);

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
