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

    for (Particle p : particles) {
      // add gravity
      if (!p.rigid) {
        if(p.pos.x > 0 && p.pos.x < 20 && p.pos.y > height - 10 && p.pos.y < height)
          p.vel.add(new PVector(1, 0));  
        p.vel.add(g);
      }
    }

    // apply viscosity
    // not yet implemented

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

    computeVorticity();


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

    for (Particle p : particles) {
      // compute new velocity
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
}
