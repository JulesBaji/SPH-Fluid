class Fluid {
  ArrayList<Particle> particles;
  PVector g;
  float dt;
  float h; // kernel radius for each particle
  float k;
  float k_near;
  float rho_0;
  Grid grid;

  Fluid(float dt, float h) {
    particles = new ArrayList<Particle>();
    g = new PVector(0, 0.01);
    g.mult(dt);
    this.dt = dt;
    this.h = h;
    k = 0.008;
    k_near = 0.01;
    rho_0 = 10;
    grid = new Grid(width, height, h);
  }

  void step() {
    // Étape 1 : Appliquer la gravité et autres forces externes
    for (Particle p : particles) {
      // add gravity
      if (!p.rigid) {
        // lid-driven cavity
        if(p.pos.x > 0 && p.pos.x < width && p.pos.y > height - 5 && p.pos.y < height)
          p.vel.add(new PVector(.1, 0));  
        p.vel.add(g);
      }
    }

    // computeVorticity();
    // applyVorticityConfinement();

    // Étape 3 : Calculer la microrotation
    computeMicrorotation();

    // Étape 4 : Coupler la vitesse avec la microrotation
    coupleVelocityAndMicrorotation();
    amplifyTurbulence();

    // Étape 2 : Mettre à jour la position en fonction de la vitesse
    for (Particle p : particles) {
      // save previous pos
      p.pos_prev = p.pos.copy();
      // advance to predicted position
      p.pos.add(PVector.mult(p.vel, dt));
    }

    // adjust springs
    // not yet implemented

    // apply spring displacements
    // not yet implemented

    // double density relaxation
    // prepare search for neighbors:
    grid.clear_grid();
    for (Particle p : particles) {
      grid.add_particle(p);
    }

    for (Particle p : particles) {
      p.find_neighbors(grid); // much faster than p.find_neighbors(particles, h);
      p.rho = 0;
      p.rho_near = 0;

      // compute density and near-density:
      for (Particle p_neighbor : p.neighbors) {
        if ((p_neighbor != p) ) {
          PVector r_ij_vec = PVector.sub(p_neighbor.pos, p.pos);
          float r_ij = r_ij_vec.mag();
          float q = r_ij / h;
          if (q < 1) {
            // rho = rho + (1-q)^2
            float temp = 1-q;
            float power_of_temp = temp * temp;
            p.rho += power_of_temp;
            // rho_near = rho_near + (1-q)^3
            power_of_temp *= temp;
            p.rho_near += power_of_temp;
          }
        }
      }

      // compute pressure and near-pressure:
      p.pressure = k*(p.rho - rho_0);
      p.pressure_near = k_near * p.rho_near;

      p.dx.mult(0);
      for (Particle p_neighbor : p.neighbors) {
        PVector r_ij_vec = PVector.sub(p_neighbor.pos, p.pos);
        float r_ij = r_ij_vec.mag();
        float q = r_ij / h;
        if (q < 1) {
          r_ij_vec.normalize();
          PVector D = r_ij_vec;
          D.mult( dt * dt * ( p.pressure * (1-q) + p.pressure_near * (1-q) * (1-q) ) );
          D.mult( 0.5 );
          if (!p_neighbor.rigid) {
            p_neighbor.pos.add(D);
          }
          p.dx.sub(D);
        }
      } // for (Particle p_neighbor : p.neighbors) {
      if (!p.rigid) {
        p.pos.add(p.dx);
      }
    } // for (Particle p : particles) {

    // resolve collisions
    // move back from out of border (could also be done with impulses, but I'm lazy here)
    for (Particle p : particles) {
      if (p.pos.x < 0) p.pos.x = 0;
      if (p.pos.x > width) p.pos.x = width;
      if (p.pos.y < 0) p.pos.y = 0;
      if (p.pos.y > height) p.pos.y = height;

      if (p.pos.x < h || p.pos.x > width - h || p.pos.y < h || p.pos.y > height - h) {
        p.vel.mult(0.9); // Réduire légèrement la vitesse
      }
    }

    // compute new velocity
    for (Particle p : particles) {
      p.vel = PVector.sub(p.pos, p.pos_prev);
      p.vel.div(dt);
    }
  }

  void plot() {
    for (Particle p : particles) {
      p.plot();
    }
  }

  void computeVorticity() {
    for (Particle p : particles) {
      float vorticityZ = 0; // La composante z de la vorticité

      for (Particle neighbor : p.neighbors) {
        if (neighbor != p && !neighbor.rigid) {
          // Vecteur directionnel entre les particules
          PVector r_ij = PVector.sub(neighbor.pos, p.pos);
          float r_ij_mag = r_ij.mag();

          if (r_ij_mag < h && r_ij_mag > 0) {
            r_ij.normalize(); // Normaliser le vecteur directionnel

            // Différence de vitesse
            PVector velocityDiff = PVector.sub(neighbor.vel, p.vel);

            // Contribution à la vorticité (seule la composante z du produit vectoriel)
            float contributionZ = velocityDiff.x * r_ij.y - velocityDiff.y * r_ij.x;

            vorticityZ += contributionZ; // Ajouter à la vorticité totale
          }
        }
      }

      // Stocker la composante z de la vorticité dans la particule
      p.vorticity = vorticityZ;
    }
  }

  void applyVorticityConfinement() {
    // Calculer la force de confinement pour chaque particule
    for (Particle p : particles) {
      if (p.rigid) {
        continue; // Ignorer les particules rigides
      }
      PVector eta = new PVector(0, 0); // Gradient de vorticité

      // Calculer le gradient de la vorticité
      for (Particle neighbor : p.neighbors) {
        if (neighbor != p) {
          // Différence de position
          PVector r_ij = PVector.sub(neighbor.pos, p.pos);
          float r_ij_mag = r_ij.mag();

          if (r_ij_mag < h && r_ij_mag > 0) {
            r_ij.normalize();
            // Différence de vorticité scalaire
            float vorticityDiff = neighbor.vorticity - p.vorticity;
            // Ajouter au gradient (pondéré par la distance et normalisé par le noyau SPH)
            eta.add(PVector.mult(r_ij, vorticityDiff));
          }
        }
      }

      // Normaliser le gradient
      if (eta.mag() > 0) {
        eta.normalize();
      }

      // Calculer la force de confinement
      // Produit vectoriel : eta (gradient) × vorticité (z en 2D)
      PVector confinementForce = new PVector(-eta.y * p.vorticity, eta.x * p.vorticity);

      // Facteur d'amplification (epsilon)
      float epsilon = 0.0001f;
      confinementForce.mult(epsilon);

      // Ajouter la force de confinement à la vitesse de la particule
      p.vel.add(confinementForce);
    }
  }

  void computeMicrorotation() {
    for (Particle p : particles) {
      PVector vorticity = new PVector(0, 0); // Calculer la vorticité locale

      for (Particle neighbor : p.neighbors) {
        if (neighbor != p) {
          // Direction et distance entre particules
          PVector r_ij = PVector.sub(neighbor.pos, p.pos);
          float r_ij_mag = r_ij.mag();

          if (r_ij_mag < h && r_ij_mag > 0) {
            r_ij.normalize();

            // Différence de vitesse entre particules
            PVector velocityDiff = PVector.sub(neighbor.vel, p.vel);

            // Contribution au champ de vorticité
            vorticity.x += velocityDiff.y * r_ij.x - velocityDiff.x * r_ij.y;
          }
        }
      }

      // Mettre à jour la microrotation (simple transfert pour fluides inviscides)
      p.torque.set(vorticity); // Stocker comme "torque" temporaire
      p.microrotation.add(PVector.mult(p.torque, dt / p.momentInertia));
    }
  }

  void coupleVelocityAndMicrorotation() {
    float nu_t = 0.0000005; // Coefficient de transfert entre vitesse et rotation
    for (Particle p : particles) {
      if (p.rigid) {
        continue; // Ignorer les particules rigides
      }
      // Transfert rotationnel vers vitesse linéaire
      PVector rotationalEffect = PVector.mult(p.microrotation, nu_t);
      p.vel.add(rotationalEffect);

      // Transfert inverse vers microrotation
      PVector velocityCurl = new PVector(-p.vel.y, p.vel.x).mult(nu_t);
      p.microrotation.add(velocityCurl);
    }
  }

  void amplifyTurbulence() {
    float amplificationFactor = 1e-5; // Facteur d'amplification

    for (Particle p : particles) {
      if (!p.rigid) {
        // Si la microrotation est au-dessus d'un seuil
        if (p.microrotation.mag() > 1000.0) {
          // Ajouter une perturbation à la vitesse
          PVector perturbation = p.microrotation.copy();
          perturbation.mult(amplificationFactor);
          p.vel.add(perturbation);
        }
      }
    }
  }
}
