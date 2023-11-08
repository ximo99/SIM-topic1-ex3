// Authors:  //<>//
// Ximo Casanova

// Problem description:
// Muelle

// Definitions:
enum IntegratorType 
{
  EXPLICIT_EULER, 
  SIMPLECTIC_EULER, 
  HEUN, 
  RK2, 
  RK4
}

// Parameters of the numerical integration:
final boolean REAL_TIME = false;
final float SIM_STEP = 0.1;   // Simulation time-step (s)
IntegratorType _integrator = IntegratorType.EXPLICIT_EULER;   // ODE integration method

// Display values:
final boolean FULL_SCREEN = false;
final int DRAW_FREQ = 50;   // Draw frequency (Hz or Frame-per-second)
int DISPLAY_SIZE_X = 500;   // Display width (pixels)
int DISPLAY_SIZE_Y = 500;   // Display height (pixels)

// Draw values:
final int [] BACKGROUND_COLOR = {200, 200, 255};
final int [] REFERENCE_COLOR = {0, 255, 0};
final int [] OBJECTS_COLOR = {255, 0, 0};
final float OBJECTS_SIZE = 1.0;   // Size of the objects (m)
final float PIXELS_PER_METER = 20.0;   // Display length that corresponds with 1 meter (pixels)
final PVector DISPLAY_CENTER = new PVector(0.0, 0.0);   // World position that corresponds with the center of the display (m)

// Parameters of the problem:
final float M = 1.0;   // Particle mass (kg)
final float Kd = 0.1;  // Constante de rozamiento
final float Ks = 1;    // Constante de elasticidad

final PVector l = new PVector(0.0, 3.0);  // Elongación en reposo
final PVector C = new PVector(DISPLAY_SIZE_X/2, DISPLAY_SIZE_Y/4);

final float Gc = 9.801;   // Gravity constant (m/(s*s))
final PVector G = new PVector(0.0, -Gc);   // Acceleration due to gravity (m/(s*s))

// Time control:
int _lastTimeDraw = 0;   // Last measure of time in draw() function (ms)
float _deltaTimeDraw = 0.0;   // Time between draw() calls (s)
float _simTime = 0.0;   // Simulated time (s)
float _elapsedTime = 0.0;   // Elapsed (real) time (s)

// Output control:
PrintWriter _output;

// Auxiliary variables:
float energy;   // Total energy of the particle (J)

// Variables to be solved:
PVector s = new PVector();   // Position of the particle (m)
PVector v = new PVector();   // Velocity of the particle (m/s)
PVector a = new PVector();   // Accleration of the particle (m/(s*s))
PVector s0 = new PVector(10.0,0.0); // Posición inicial

// Main code:

// Converts distances from world length to pixel length
float worldToPixels(float dist)
{
  return dist*PIXELS_PER_METER;
}

// Converts distances from pixel length to world length
float pixelsToWorld(float dist)
{
  return dist/PIXELS_PER_METER;
}

// Converts a point from world coordinates to screen coordinates
void worldToScreen(PVector worldPos, PVector screenPos)
{
  screenPos.x = 0.5*DISPLAY_SIZE_X + (worldPos.x - DISPLAY_CENTER.x)*PIXELS_PER_METER;
  screenPos.y = 0.5*DISPLAY_SIZE_Y - (worldPos.y - DISPLAY_CENTER.y)*PIXELS_PER_METER;
}

// Converts a point from screen coordinates to world coordinates
void screenToWorld(PVector screenPos, PVector worldPos)
{
  worldPos.x = ((screenPos.x - 0.5*DISPLAY_SIZE_X)/PIXELS_PER_METER) + DISPLAY_CENTER.x;
  worldPos.y = ((0.5*DISPLAY_SIZE_Y - screenPos.y)/PIXELS_PER_METER) + DISPLAY_CENTER.y;
}

void drawStaticEnvironment()
{
  background(BACKGROUND_COLOR[0], BACKGROUND_COLOR[1], BACKGROUND_COLOR[2]);

  fill(REFERENCE_COLOR[0], REFERENCE_COLOR[1], REFERENCE_COLOR[2]);
  strokeWeight(2);
  textSize(16);

  PVector screenPos = new PVector();
  worldToScreen(new PVector(), screenPos);
  
  fill(0,0,0);
  text("Integrador: "+ _integrator, 100, DISPLAY_SIZE_Y - 100);
}

void drawMovingElements()
{
  fill(OBJECTS_COLOR[0], OBJECTS_COLOR[1], OBJECTS_COLOR[2]);
  strokeWeight(2);

  PVector screenPos = new PVector();
  worldToScreen(s, screenPos);
  
  PVector screenPos2 = new PVector();
  worldToScreen(C, screenPos2);
  
  line(screenPos2.x, screenPos2.y, screenPos.x , screenPos.y);
  fill(OBJECTS_COLOR[0], OBJECTS_COLOR[1], OBJECTS_COLOR[2]);
  circle(screenPos.x, screenPos.y, worldToPixels(OBJECTS_SIZE));
}

void PrintInfo()
{
  println("Energy: " + energy + " J");
  println("Elapsed time = " + _elapsedTime + " s");
  println("Simulated time = " + _simTime + " s \n");
}

void initSimulation()
{
  _simTime = 0.0;
  _elapsedTime = 0.0;
  
  s = s0.copy();
  
  v.set(0.0,0.0);
  a.set(0.0,0.0);
}

void updateSimulation()
{
    
    
  switch (_integrator)
  {
  case EXPLICIT_EULER:
    updateSimulationExplicitEuler();
    break;

  case SIMPLECTIC_EULER:
    updateSimulationSimplecticEuler();
    break;

  case HEUN:
    updateSimulationHeun();
    break;

  case RK2:
    updateSimulationRK2();
    break;

  case RK4:
    updateSimulationRK4();
    break;
  }
  
  _simTime += SIM_STEP;
}

void updateSimulationExplicitEuler()
{
  // Calcular la derivada en el principio del intervalo
  a = calculateAcceleration(s, v);
  
  // Calcular la posición siguiente a partir de la velocidad en el principio del intervalo
  s.add(PVector.mult(v, SIM_STEP));
  
  // Calcular la velocidad siguiente a partir de la derivada en el principio del intervalo
  v.add(PVector.mult(a, SIM_STEP));
}

void updateSimulationSimplecticEuler()
{
  // Calcular la derivada en el principio del intervalo
  a = calculateAcceleration(s, v);
  
  // Calcular la velocidad siguiente a partir de la derivada en el principio del intervalo  
  v.add(PVector.mult(a, SIM_STEP));
  
  // Calcular la posición siguiente a partir de la velocidad en el principio del intervalo  
  s.add(PVector.mult(v, SIM_STEP));  
}

void updateSimulationHeun()
{
  // Parte 1, integración numérica de la velocidad
  // Calcular aceleracion a
  a = calculateAcceleration(s, v);
  
  // Paso de Euler, actualizo s2 y v2 (velocidad y posición al final del intervalo)
  PVector s2 = new PVector();
  s2 = s;
  s2.add(PVector.mult(v, SIM_STEP));
  PVector v2 = new PVector();
  v2 = v;
  
  // Cálculo de la velocidad promedio a partir de v y v2
  PVector v_prom = PVector.mult(PVector.add(v, v2), 0.5);
  
  //actualizar s con la v promedio
  s.add(PVector.mult(v_prom, SIM_STEP));
  
  // Parte 2, integración de la acceleración
  // Calcular la aceleración a2 (aceleración al final del intervalo)
  PVector a2 = new PVector();
  a2 = calculateAcceleration(s2, v2);
  
  // Cálculo de la aceleración promedio a partir de a y a2
  PVector a_prom = PVector.mult(PVector.add(a, a2), 0.5);
  
  //Actualizar la velocidad con la aceleración promedio
  v.add(PVector.mult(a_prom, SIM_STEP));
}

void updateSimulationRK2()
{
  //Calcular acceleracion a
  a = calculateAcceleration(s,v);
  
  // k1s = v(t) * h
  PVector k1s = PVector.mult(v, SIM_STEP);
  
  // k1v = a(s(t), v(t)) * h
  PVector k1v = PVector.mult(a, SIM_STEP);
  
  PVector s2 = PVector.add(s, PVector.mult(k1s, 0.5));
  PVector v2 = PVector.add(v, PVector.mult(k1v, 0.5));
  PVector a2 = calculateAcceleration(s2,v2);
  
  // k2v = a(s(t) + k1s / 2, v(t) + k1v / 2) * h
  PVector k2v = PVector.mult(a2,SIM_STEP);
  
  // k2s = (v(t)+k1v/2)*h
  PVector k2s = PVector.mult(PVector.add(v,PVector.mult(k1v,0.5)),SIM_STEP);
  
  v.add(k2v);
  s.add(k2s);
}

void updateSimulationRK4()
{
  // Calcular acceleracion a
  a = calculateAcceleration(s,v);

  // k1v = a(s(t), v(t)) * h
  PVector k1v = PVector.mult(a,SIM_STEP);
  
  // k1s = v(t) * h
  PVector k1s = PVector.mult(v, SIM_STEP);
  
  // k2v = a(s(t) + k1s / 2, v(t) + k1v / 2) * h  
  PVector s2 = PVector.add(s, PVector.mult(k1s, 0.5));
  PVector v2 = PVector.add(v, PVector.mult(k1v, 0.5));
  PVector a2 = calculateAcceleration(s2,v2);
  PVector k2v = PVector.mult(a2,SIM_STEP);
  PVector k2s = PVector.mult(PVector.add(v,PVector.mult(k1v,0.5)),SIM_STEP);
  
  PVector s3 = PVector.add(s, PVector.mult(k2s, 0.5));
  PVector v3 = PVector.add(v, PVector.mult(k2v, 0.5));
  PVector a3 = calculateAcceleration(s3,v3);
  
  // k3v = a(s(t)+k2s/2, v(t)+k2v/2)*h
  PVector k3v = PVector.mult(a3,SIM_STEP);
  
  // k3s = (v(t)+k2v/2)*h
  PVector k3s = PVector.mult(PVector.add(v,PVector.mult(k2v,0.5)),SIM_STEP);
  
  PVector s4 = PVector.add(s, k3s);
  PVector v4 = PVector.add(v, k3v);
  PVector a4 = calculateAcceleration(s4,v4);
  
  // k4v = a(s(t)+k3s, v(t)+k3v)*h
  PVector k4v = PVector.mult(a4,SIM_STEP);
  
  // k4s = (v(t)+k3v)*h
  PVector k4s = PVector.mult(PVector.add(v,k3s),SIM_STEP);
  
  // v(t+h) = v(t) + (1/6)*k1v + (1/3)*k2v + (1/3)*k3v +(1/6)*k4v
  v.add(PVector.mult(k1v,1/6.0));
  v.add(PVector.mult(k2v, 1/3.0));
  v.add(PVector.mult(k3v, 1/3.0));
  v.add(PVector.mult(k4v, 1/6.0));
  
  // s(t+h) = s(t) + (1/6)*k1s + (1/3)*k2s + (1/3)*k3s +(1/6)*k4s
  s.add(PVector.mult(k1s,1/6.0));  
  s.add(PVector.mult(k2s, 1/3.0));
  s.add(PVector.mult(k3s, 1/3.0));
  s.add(PVector.mult(k4s, 1/6.0));
}

PVector calculateAcceleration(PVector s, PVector v)
{
  PVector Fp = new PVector(0,0);  // Fuerza peso
  PVector Fe;  // Fuerza elastica
  PVector Fr;  // Fuerza de rozamiento del plano
  PVector F;   // Fuerza total
  
  Fp = PVector.mult(G, M); 
  Fr = PVector.mult(v, Kd);
  Fr.dot(v);
  
  Fe = PVector.mult(PVector.sub(PVector.sub(s, C), l), -Ks);
  F = PVector.sub(Fe, Fr);
  
  F.add(Fp);
  a = PVector.div(F , M);
  
  return a;
}

void calculateEnergy()
{  
  float Ec = pow(v.mag(), 2)*M/2;
  float Ep_G = M*Gc*s.y;
  float Ep_E = Ks*(pow(PVector.sub(PVector.sub(s, C), l).mag(), 2))/2;
  
  energy = Ec + Ep_G + Ep_E;
}

void settings()
{
  if (FULL_SCREEN)
  {
    fullScreen();
    DISPLAY_SIZE_X = displayWidth;
    DISPLAY_SIZE_Y = displayHeight;
  } 
  else
    size(DISPLAY_SIZE_X, DISPLAY_SIZE_Y);
}

void setup()
{
  frameRate(DRAW_FREQ);
  _lastTimeDraw = millis();
  screenToWorld(C, C);
  
  _output = createWriter("data.txt");

  initSimulation();
}

void draw()
{
  int now = millis();
  _deltaTimeDraw = (now - _lastTimeDraw)/1000.0;
  _elapsedTime += _deltaTimeDraw;
  _lastTimeDraw = now;

  println("\nDraw step = " + _deltaTimeDraw + " s - " + 1.0/_deltaTimeDraw + " Hz");

  if (REAL_TIME)
  {
    float expectedSimulatedTime = 1.0*_deltaTimeDraw;
    float expectedIterations = expectedSimulatedTime/SIM_STEP;
    int iterations = 0; 

    for (; iterations < floor(expectedIterations); iterations++)
      updateSimulation();

    if ((expectedIterations - iterations) > random(0.0, 1.0))
    {
      updateSimulation();
      iterations++;
    }
  } 
  else
    updateSimulation();

  drawStaticEnvironment();
  drawMovingElements();

  calculateEnergy();
  PrintInfo();
  
  if (_simTime%0.5 < 0.1)
    _output.println(_integrator + " Posición: (" + s.x + ", " + s.y + "); Velocidad: (" + v.x + ", " + v.y + "); Energía: " + energy );
}

void mouseClicked() 
{
  PVector aux = new PVector(mouseX, mouseY);
  
  screenToWorld(aux, s0);
  
  
}

void keyPressed()
{
  switch(key)
  {
    case 'e':
    case 'E':
      _integrator = IntegratorType.EXPLICIT_EULER;
    break;
    
    case 's':
    case 'S':
      _integrator = IntegratorType.SIMPLECTIC_EULER;
    break;
    
    case 'h':
    case 'H':
      _integrator = IntegratorType.HEUN;
    break;
    
    case '2':
      _integrator = IntegratorType.RK2;
    break;
    
    case '4':
      _integrator = IntegratorType.RK4;
    break;
    
    case 'r':
    case 'R':
      initSimulation();
    break;
    
    case 'q':
    case 'Q':
      stop();
    break;
  }
}

void stop()
{
  _output.flush();
  _output.close();
  exit();
}
