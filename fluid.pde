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

    computeMicrorotation();
    coupleVelocityAndMicrorotation();

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
        p.vel.mult(0.99); // Réduire légèrement la vitesse
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
}
