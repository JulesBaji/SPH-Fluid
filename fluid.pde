class Fluid {
  ArrayList<Particle> particles;
  PVector g;
  float dt;
  float h; // kernel radius for each particle
  float rho_0;
  Grid grid;
  float test;

  Fluid(float dt, float h) {
    particles = new ArrayList<Particle>();
    g = new PVector(0, 0.07);
    g.mult(dt);
    this.dt = dt;
    this.h = h;
    rho_0 = 2e-3;
    grid = new Grid(width, height, h);
    test = 0;
  }

  void step() {
    // Étape 1 : Appliquer la gravité et autres forces externes
    for (Particle p : particles) {
      if (!p.rigid) {
        // Lid-driven cavity : Appliquer une vitesse horizontale dans une zone spécifique
        if (p.pos.x > 0 && p.pos.x < width/2 && p.pos.y > height - 2 * h && p.pos.y < height) {
          p.vel.add(new PVector(2e-2, 0));
        }
        if (p.pos.x > width/2 && p.pos.x < width && p.pos.y > height - 2 * h && p.pos.y < height) {
          p.vel.add(new PVector(-2e-2, 0));
        }
        // Appliquer la gravité
        p.vel.add(g);
      }
    }

    // Étape 2 : Calcul des positions prédites
    for (Particle p : particles) {
      if (p.rigid) continue;
      p.pos_prev = p.pos.copy(); // Sauvegarder la position précédente
      p.pos.add(PVector.mult(p.vel, dt)); // Calcul de la position prédite
    }

    // Étape 3 : Mettre à jour la grille pour les recherches de voisins
    grid.clear_grid();
    for (Particle p : particles) {
      grid.add_particle(p);
      p.find_neighbors(grid);
    }

    // Étape 4 : Calcul des densités
    computeDensity(); // Fonction DFSPH pour calculer la densité locale

    // Étape 5 : Correction de la densité (incompressibilité)
    correctDensity(); // Fonction DFSPH pour ajuster les positions et garantir ρ ≈ ρ₀

    // Étape 6 : Correction de la divergence
    // correctDivergence(); // Fonction DFSPH pour ajuster les vitesses et garantir ∇·v = 0

    // Étape 7 : Résolution des collisions par impulsion
    for (Particle p : particles) {
      if (p.rigid) continue;

      // Collision avec les bordures
      if (p.pos.x < 0) {
        p.pos.x = 0;
        PVector normal = new PVector(1, 0); // Normale à la bordure gauche
        applyImpulse(p, normal);
      } else if (p.pos.x > width) {
        p.pos.x = width;
        PVector normal = new PVector(-1, 0); // Normale à la bordure droite
        applyImpulse(p, normal);
      }

      if (p.pos.y < 0) {
        p.pos.y = 0;
        PVector normal = new PVector(0, 1); // Normale à la bordure inférieure
        applyImpulse(p, normal);
      } else if (p.pos.y > height) {
        p.pos.y = height;
        PVector normal = new PVector(0, -1); // Normale à la bordure supérieure
        applyImpulse(p, normal);
      }
    }

    // Étape 8 : Mise à jour des vitesses (post-correction)
    for (Particle p : particles) {
      if (p.rigid) continue;
      // Calculer la nouvelle vitesse en fonction des positions corrigées
      p.vel = PVector.sub(p.pos, p.pos_prev).div(dt);

      // Couplage avec la microrotation
      coupleVelocityAndMicrorotation(p);
    }
  }

  void plot() {
    for (Particle p : particles) {
      p.plot();
    }
  }

  void coupleVelocityAndMicrorotation(Particle p) {
    float vorticity = 0; // Initialiser la vorticité locale pour la particule

    for (Particle neighbor : p.neighbors) {
      if (neighbor != p && !neighbor.rigid) {
        // Direction et distance entre particules
        PVector r_ij = PVector.sub(neighbor.pos, p.pos);
        float r_ij_mag = r_ij.mag();

        if (r_ij_mag < h && r_ij_mag > 0) { // Voisins dans le rayon d'influence
          r_ij.normalize(); // Normaliser la direction

          // Différence de vitesse entre les particules
          PVector velocityDiff = PVector.sub(neighbor.vel, p.vel);

          // Contribution scalaire à la vorticité (produit vectoriel 2D)
          vorticity += velocityDiff.x * r_ij.y - velocityDiff.y * r_ij.x;
        }
      }
    }

    // Appliquer le torque basé sur la vorticité
    p.torque = vorticity;

    // Mettre à jour la microrotation (Δω = τ · Δt / Θ)
    p.microrotation += (p.torque * dt) / p.momentInertia;
    p.microrotation *= 0.98; // Damping de la microrotation

    // Mise à jour de la microrotation en tenant compte du transfert avec la vorticité
    // Transfert vitesse -> microrotation
    float transferCoefficient = 0.0000001; // Contrôle la force du couplage
    p.microrotation += (vorticity * transferCoefficient) * dt;

    // Transfert microrotation -> vitesse
    PVector rotationalEffect = new PVector(-p.microrotation * transferCoefficient, p.microrotation * transferCoefficient);
    p.vel.add(rotationalEffect);
  }

  void applyImpulse(Particle p, PVector normal) {
    // Calcul de la vitesse normale (projection de la vitesse sur la normale)
    float vn = PVector.dot(p.vel, normal);

    // Si la particule s'éloigne déjà de la surface, pas besoin d'impulsion
    if (vn >= 0) return;

    // Réponse impulsive : inverser et réduire la composante normale de la vitesse
    float restitution = 0.8; // Coefficient de restitution (élasticité)
    float impulseMagnitude = -(1 + restitution) * vn;

    // Appliquer l'impulsion
    PVector impulse = PVector.mult(normal, impulseMagnitude);
    p.vel.add(impulse);
  }

  void computeDensity() {
    for (Particle p : particles) {
      if (p.rigid) continue;
      p.rho = 0; // Réinitialiser la densité
      for (Particle neighbor : p.neighbors) {
        float r_ij_mag = PVector.dist(p.pos, neighbor.pos);
        if (r_ij_mag < h) {
          // Ajout de la contribution au noyau SPH
          float test = 1e5; 
          p.rho += neighbor.mass * kernel(r_ij_mag, h);
        }
      }
    }
  }

  void correctDensity() {
    for (Particle p : particles) {
      if (p.rigid) continue;
      float constraint = (p.rho - rho_0) / rho_0; // Écart de densité relatif

      // Correction de position
      PVector correction = new PVector(0, 0);
      for (Particle neighbor : p.neighbors) {
        float r_ij_mag = PVector.dist(p.pos, neighbor.pos);
        if (r_ij_mag < h) {
          // Calcul de la correction pondérée par le gradient du noyau
          PVector grad = gradientKernel(p.pos, neighbor.pos, h);
          // Multiplier gradient ?
          correction.add(PVector.mult(grad, neighbor.mass * constraint));
        }
      }
      
      // Appliquer la correction inverse
      p.pos.sub(PVector.mult(correction, 1e7));
    }
  }

  void correctDivergence() {
    for (Particle p : particles) {
      if (p.rigid) continue;
      float divergence = 0;

      // Calcul de la divergence de vitesse
      for (Particle neighbor : p.neighbors) {
        float r_ij_mag = PVector.dist(p.pos, neighbor.pos);
        if (r_ij_mag < h) {
          PVector velDiff = PVector.sub(p.vel, neighbor.vel);
          PVector grad = gradientKernel(p.pos, neighbor.pos, h);
          divergence += neighbor.mass * PVector.dot(velDiff, grad);
          // divergence += PVector.dot(velDiff, grad);
        }
      }
      // divergence /= p.rho;

      // Appliquer la correction de divergence
      PVector correction = new PVector(0, 0);
      for (Particle neighbor : p.neighbors) {
        float r_ij_mag = PVector.dist(p.pos, neighbor.pos);
        if (r_ij_mag < h) {
          PVector grad = gradientKernel(p.pos, neighbor.pos, h);
          // grad.mult(1e15);
          correction.add(PVector.mult(grad, neighbor.mass * divergence));
        }
      }
      // p.vel.sub(correction);
      p.vel.sub(PVector.mult(correction, 1));
    }
  }

  // Cubic Spline Kernel
  float kernel(float r, float h) {
    float q = r / h;
    if (q >= 0 && q <= 1) {
      return (315.0f / (64.0f * PI * pow(h, 3))) * pow((1 - q * q), 3);
    }
    return 0;
  }
  
  // Cubic Spline Kernel
  PVector gradientKernel(PVector pos_i, PVector pos_j, float h) {
    PVector r = PVector.sub(pos_i, pos_j);
    float r_mag = r.mag();
    float q = r_mag / h;
    if (q >= 0 && q <= 1) {
      float factor = (-45.0f / (PI * pow(h, 6))) * pow((1 - q), 2);
      return r.normalize().mult(factor);
    }
    return new PVector(0, 0);
  }
}