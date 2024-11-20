class Particle {
  PVector pos;
  PVector vel;
  PVector pos_prev;
  ArrayList<Particle> neighbors;
  float rho;
  float rho_near;
  float pressure;
  float pressure_near;
  PVector dx;
  boolean rigid;

  float vorticity; // Nouvelle propriété pour stocker la vorticité
  PVector microrotation; // Nouvelle propriété : microrotation (ω)
  PVector torque;      // Couple appliqué sur la particule (τ)
  float momentInertia; // Moment d'inertie (Θ)

  Particle(float x, float y, boolean r) {
    pos = new PVector(x, y);
    vel = new PVector();
    neighbors = new ArrayList<Particle>();
    rho = 0;
    rho_near = 0;
    dx = new PVector(0, 0);
    rigid = r;

    this.vorticity = 0; // Initialiser à zéro
    this.microrotation = new PVector(0, 0); // Initialiser la rotation à zéro
    this.torque = new PVector(0, 0);
    this.momentInertia = 2; // Valeur expérimentale pour fluides micropolaires
  }

  void plot() {
    strokeWeight(8);
    if (rigid) {
      stroke(0);
    } else {
      stroke(150 - microrotation.mag()*0.01, 255, 200);
    }
    point(pos.x, pos.y);
  }

  void find_neighbors(ArrayList<Particle> particles, float h) {
    // legacy detection by checking the whole particle Arraylist, works without Grid / Field classes
    // but is SLOW for large particle numbers. Currently not used.
    neighbors.clear();
    for (int j=0; j<particles.size(); j++) {
      Particle p2 = particles.get(j);
      if ( (abs(pos.x - p2.pos.x)<h) && (abs(pos.y - p2.pos.y)<h) ) {
        neighbors.add(p2);
      }
    }
  }

  void find_neighbors(Grid g) {
    Field f = g.get_Field(pos.x, pos.y);
    if (f != null) {
      neighbors = f.particles;
    }
  }
}
