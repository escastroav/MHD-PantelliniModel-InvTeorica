const double size = 1;
const int N = 50*size;
const double R0 = 1/size;

const double h = N*R0;
const double g = 3*R0;
const double H = h*0.2; // kT/mg=rmax/5

const double tf = sqrt(H/g);
//
class Crandom{
  unsigned long long u,v,w;
public:
  Crandom(unsigned long long j);
  unsigned long long int64();
  double r() {return 5.42101086242752217E-20 * int64();}
  unsigned int int32(){return (unsigned int) int64();};
  double exponencial(float tau);
  double gauss(float mu,float sigma);
};
/*
Crandom::Crandom(unsigned long long j);
unsigned long long Crandom::int64();
double Crandom::exponencial(float tau);
double Crandom:: gauss(float mu,float sigma);
*/
//
class Ball
{
private:
  int i = 0;
  double z=0,
    Vx=0,
    Vy=0,
    Vz=0;
  double az=0;
public:
  void InitBall(int i0, double z0, double Vx0, double Vy0, double Vz0);
  double GetZ(void){return z;};void SetZ(double zi){z=zi;};
  double GetVx(void){return Vx;};void SetVx(double Vxi){Vx=Vxi;};
  double GetVy(void){return Vy;};void SetVy(double Vyi){Vy=Vyi;};
  double GetVz(void){return Vz;};void SetVz(double Vzi){Vz=Vzi;};
  int GetI(void){return i;};
  void MovementEq(void);
  void Integrate(double dt);
  void PrintBall(void);
  friend class Collider;
};

/*
void Ball::InitBall(int i0, double z0, double Vx0, double Vy0, double Vz0);
void Ball::MovementEq(void);
void Ball::Integrate(double dt);
void Ball::PrintBall(void);
*/

class Collider
{
private:
  double t[N];
  double tmin=0.0;
  double I = 0;
  
public:
  Collider(void);
  double GetTmin(void){return tmin;};
  int GetIndex(void){return I;};
  void GetCollisionTime(Ball & ball1, Ball & ball2);
  void CollisionTimeBound(Ball & ball0);
  void CollideBalls(Ball & ball1, Ball & ball2, Crandom rand, double dt);
  void CollideBound(Ball & ball0);
  void InitPositions(Ball * balls, Crandom ran);
  void MeasureDensity(Ball * balls);
  bool AreBallsSorted(Ball * balls);
};

/*
Collider::Collider(void);
void Collider::GetCollisionTime(Ball & ball1, Ball & ball2);
void Collider::CollisionTimeBound(Ball & ball0);
void Collider::CollideBalls(Ball & ball1, Ball & ball2, Crandom ran, double dt);
void Collider::CollideBound(Ball & ball0);
void Collider::InitPositions(Ball * balls, Crandom ran);
void Collider::MeasureDensity(Ball * balls);
bool Collider::AreBallsSorted(Ball * balls);
void InicieAnimacion(void);
void InicieCuadro(void);
void TermineCuadro(void);
*/
