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
        if(p.pos.x > 0 && p.pos.x < width && p.pos.y > height - 2*h && p.pos.y < height)
          p.vel.add(new PVector(1.2e-3, 0));  
        p.vel.add(g);
      }
    }

    for (Particle p : particles) {
      if (p.rigid) continue;
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
      if (p.rigid) continue;
      grid.add_particle(p);
    }

    for (Particle p : particles) {
      if (p.rigid) continue;
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
    for (Particle p : particles) {
      if (p.rigid) continue;
      if (p.pos.x < 0) {
        p.pos.x = 0; 
        PVector normal = new PVector(1, 0); 
        applyImpulse(p, normal);
      } else if (p.pos.x > width) {
        p.pos.x = width;
        PVector normal = new PVector(-1, 0); 
        applyImpulse(p, normal);
      }

      if (p.pos.y < 0) {
        p.pos.y = 0;
        PVector normal = new PVector(0, 1); 
        applyImpulse(p, normal);
      } else if (p.pos.y > height) {
        p.pos.y = height;
        PVector normal = new PVector(0, -1); 
        applyImpulse(p, normal);
      }
    }

    // compute new velocity
    for (Particle p : particles) {
      if(p.rigid) continue;
      p.vel = PVector.sub(p.pos, p.pos_prev);
      computeMicrorotation(p);
      coupleVelocityAndMicrorotation(p);
      p.vel.div(dt);
    }
  }

  void plot() {
    for (Particle p : particles) {
      p.plot();
    }
  }

  void computeMicrorotation(Particle p) {
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
    p.microrotation *= 0.995; // Damping de la microrotation
  }

  void coupleVelocityAndMicrorotation(Particle p) {
    float nu_t = 0.0000001; // Coefficient de transfert entre rotation et vitesse
    // Transfert microrotation -> vitesse
    PVector rotationalEffect = new PVector(-p.microrotation * nu_t, p.microrotation * nu_t);
    p.vel.add(rotationalEffect);

    // Transfert vitesse -> microrotation
    float velocityCurl = (p.vel.y - p.vel.x) * nu_t; // Approximation simple pour 2D
    p.microrotation += velocityCurl;
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
}