#pragma once
using namespace std;

void ConsoleMenu();
void FileRead();

class VectorLenghtExeption : public exception {
private:
	char* msg;
public:
	VectorLenghtExeption(const char* g);
	~VectorLenghtExeption(){}
	char* getMessage() const { return msg; }	
};
void ThrowVectorExcetpion();



class Element {
public:
	virtual ~Element(){}
};



class Point : public Element {
private:
	double x, y, z;
public:
	friend ostream& operator <<(ostream& out, const Point & m);
	friend istream& operator >>(istream& in, Point & m);

	Point();
	Point(double x, double y, double z);
	Point(const Point & s);
	~Point() {}

	void SetX(double x);
	void SetY(double y);
	void SetZ(double z);
	double GetX() const;
	double GetY() const;
	double GetZ() const;

	Point& operator = (const Point& s);
	bool operator == (const Point& s) const;
};

class EqualPointException : public exception {
private:
	Point A;
	Point B;
public:
	EqualPointException(const Point& s, const Point& t);
	~EqualPointException() {}
	void getMessage() const { 
		cout << "Точките " << A << " и " << B << " съвпадат.";
	}
};
void ThrowPointException(const Point& t, const Point& p);


class Vector : public Point {
public:
	friend ostream& operator <<(ostream& out, const Vector & m);
	friend istream& operator >>(istream& in, Vector & m);

	Vector();
	Vector(double x, double y, double z);
	Vector(Point S, Point T);
	Vector(const Vector& s);
	~Vector(){}

	double VectorLen()const;
	Vector VectorDirection();
	Vector VectorProjection(const Vector& s);
	Vector CrossProduct(const Vector& s);
	bool NullVector();
	bool ParallelVector(const Vector& s);
	bool PerpendicularVector(const Vector& s);

	Vector& operator = (const Vector& s);
	Vector operator +(const Vector & s);
	Vector operator -(const Vector & s);
	Vector operator *(double s);
	Vector operator ^(const Vector & s);
	double operator ()(const Vector & s,const Vector& t);
};

double operator *(const Vector& t, const Vector & s);


class Line : public Vector {
private:
	Point A;
	Point B;
	Vector X;
public:
	friend ostream& operator <<(ostream& out, const Line & m);
	friend istream& operator >>(istream& in, Line & m);

	Line();
	Line(Point &a, Point &b);
	Line(Point &a, Vector &b);
	Line(Line &s);
	~Line() {}

	Vector LineDirection();
	Vector NormalVector();
	double LinesAngle(Line& b);

	void SetA(const Point& A);
	void SetB(const Point& B);
	void SetX(const Vector& X);
	Point GetA() const;
	Point GetB() const;
	Vector GetX() const;

	Line& operator = (const Line& s);
};


class Segment : public Line {
public:
	friend ostream& operator <<(ostream& out, const Segment & m);
	friend istream& operator >>(istream& in, Segment & m);

	Segment();
	Segment(Point &a, Point &b);
	Segment(Segment& s);
	~Segment() {}

	float SegmentLen();
	Point SegmentMid();

	Segment& operator = (const Segment& s);

};

class Triangle : public Point {
private:
	Point A;
	Point B;
	Point C;
public:
	friend ostream& operator <<(ostream& out, const Triangle & m);
	friend istream& operator >>(istream& in, Triangle & m);

	Triangle();
	Triangle(Point A, Point B, Point C);
	Triangle(const Triangle & s);
	~Triangle(){}

	void TriangleType();
	float TriangleArea();
	float TrianglePerimeter();
	Point Medicenter();

	Triangle SetPointA(const Point& A);
	Triangle SetPointB(const Point& B);
	Triangle SetPointC(const Point& C);

	Point GetA() const;
	Point GetB() const;
	Point GetC() const;

	Triangle& operator = (const Triangle& s);
};


class Tetrahedron : public Triangle {
private:
	Point A, B, C, D;
public:
	friend ostream& operator <<(ostream& out, const Tetrahedron & m);
	friend istream& operator >>(istream& in, Tetrahedron & m);

	Tetrahedron();
	Tetrahedron(Point a, Point b, Point c, Point d);
	Tetrahedron(const Tetrahedron& s);
	~Tetrahedron(){}

	bool RegularTetrahedron();
	bool OrthocentricTetrahedron();
	double TetrahedronArea();
	double TetrahedronVolume();

	Tetrahedron SetPointA(const Point& A);
	Tetrahedron SetPointB(const Point& B);
	Tetrahedron SetPointC(const Point& C);
	Tetrahedron SetPointD(const Point& D);
	Point GetA() const;
	Point GetB() const;
	Point GetC() const;
	Point GetD() const;
};


float angle(Point x, Point y, Point z);

bool operator ==(const Vector& s, const Vector& t);


bool operator + (const Point& P, const Line& s);
bool operator ||(const Line& t, const Line& s);
bool operator ==(const Line& t, const Line& s);
bool operator !=(const Line& t, const Line& s);
bool operator &&(const Line& t, const Line& s);
bool operator |(const Line& t, const Line& s);

bool operator ==(const Point& s,const Segment& t);


bool operator < (const Point& t, const Triangle& s);
bool operator > (const Point& t, const Triangle& s);
bool operator == (const Point& t, const Triangle& s);


bool operator < (const Point& t, const Tetrahedron& s);
bool operator > (const Point& t, const Tetrahedron& s);
bool operator == (const Point& t, const Tetrahedron& s);


ostream& operator <<(ostream& out, const Point & m);
ostream& operator <<(ostream& out, const Vector & m);
ostream& operator <<(ostream& out, const Triangle & m);
ostream& operator <<(ostream& out, const Line & m);
ostream& operator <<(ostream& out, const Segment & m);
ostream& operator <<(ostream& out, const Tetrahedron & m);

istream& operator >>(istream& in, Point & m);
istream& operator >>(istream& in, Vector & m);
istream& operator >>(istream& in, Triangle & m);
istream& operator >>(istream& in, Line & m);
istream& operator >>(istream& in, Segment & m);
istream& operator >>(istream& in, Tetrahedron & m);
