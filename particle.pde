class Particle {
  PVector pos;
  PVector vel;
  PVector pos_prev;
  ArrayList<Particle> neighbors;
  float rho;          // Densité de la particule
  PVector dx;
  boolean rigid;
  float mass;          // Masse de la particule

  float microrotation; // Microrotation scalaire (ω)
  float torque;        // Couple appliqué (τ)
  float momentInertia; // Moment d'inertie (Θ)

  Particle(float x, float y, boolean r) {
    pos = new PVector(x, y);
    vel = new PVector();
    neighbors = new ArrayList<Particle>();
    rho = 0;
    dx = new PVector(0, 0);
    rigid = r;

    this.microrotation = 0; 
    this.torque = 0;        
    this.momentInertia = 1; 
    this.mass = 1;          
  }

  void plot() {
    strokeWeight(8);
    if (rigid) {
      stroke(0);
    } else {
      stroke(150 - microrotation/100, 255, 200); // Coloration en fonction de la microrotation
      // stroke(150 - vel.mag()*10, 255, 200); // Coloration en fonction de la vitesse
      // stroke(150 - rho*1e4, 255, 200); // Coloration en fonction de la densité
    }
    point(pos.x, pos.y);
  }

  void find_neighbors(Grid g) {
    Field f = g.get_Field(pos.x, pos.y);
    if (f != null) {
      neighbors = f.particles;
    }
  }
}
