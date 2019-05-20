#include "pch.h"
#include <iostream>
#include <cmath>
#include <locale>
#include <fstream>
#include "Vectors.h"

using namespace std;

#pragma warning (disable : 4996)


VectorLenghtExeption::VectorLenghtExeption(const char* n) {
	msg = new char[strlen(n) + 1];
	strcpy(msg, n);
}
void ThrowVectorExcetpion() {
	throw VectorLenghtExeption("Vector Lenght Exception - векторът е с дължина 0.");
}



EqualPointException::EqualPointException(const Point& s, const Point& t){
	A = s;
	B = t;
}
void ThrowPointException(const Point& t,const Point& p) {
	throw EqualPointException(t, p);
}




ostream& operator <<(ostream& out, const Point & m) {
	cout << "(" << m.GetX() << ", " << m.GetY() << ", " << m.GetZ() << ")";
	return out;
}
ostream& operator <<(ostream& out, const Vector & m) {
	cout << "(" << m.GetX() << ", " << m.GetY() << ", " << m.GetZ() << ")";
	return out;
}
ostream& operator <<(ostream& out, const Triangle & m) {
	cout  << "A = " << m.A << ", " << "B = " << m.B << ", " << "C = " << m.C << endl;
	return out;
}
ostream& operator <<(ostream& out, const Line & m) {
	cout << "Точки: " << m.GetA() << m.GetB() << "Посока: " << m.GetX() << endl;
	return out;
}
ostream& operator <<(ostream& out, const Segment & m) {
	cout << "Начало на отсечка: " << m.GetA() << ", край на отсечка: " << m.GetB() << endl;
	return out;
}
ostream& operator <<(ostream& out, const Tetrahedron & m){
	cout << "A = " << m.GetA() << ", " << "B = " << m.GetB() << ", " << "C = " << m.GetC() << ", " << "D = " << m.GetD();
	return out;
}




istream& operator >>(istream& in, Point & m) {
	double x, y, z;
	cout << endl;
	cout << "Моля въведете стойност за x на точката: " << endl;
	cin >> x;
	cout << "Моля въведете стойност за y на точката: " << endl;
	cin >> y;
	cout << "Моля въведете стойност за z на точката: " << endl;
	cin >> z;

	m.SetX(x);
	m.SetY(y);
	m.SetZ(z);

	return in;
}
istream& operator >>(istream& in, Vector & m) {
	double x, y, z;

	cout << endl;
	cout << "Моля въведете стойност за x на вектора: " << endl;
	cin >> x;
	cout << "Моля въведете стойност за y на вектора: " << endl;
	cin >> y;
	cout << "Моля въведете стойност за z на вектора: " << endl;
	cin >> z;

	m.SetX(x);
	m.SetY(y);
	m.SetZ(z);

	return in;
}
istream& operator >>(istream& in, Triangle & m) {
	cout << "Въведете точките на триъгълника: " << endl;
	Point a;
	Point b;
	Point c;
	cin >> a >> b >> c;

	m.SetPointA(a);
	m.SetPointB(b);
	m.SetPointC(c);

	return in;
}
istream& operator >>(istream& in, Line & m) {
	cout<< endl << "Въведете точки за линията: " << endl;
	Point a;
	Point b;
	cin >> a >> b;
	Vector x(a, b);

	m.SetA(a);
	m.SetB(b);
	m.SetX(x);
	
	return in;
}
istream& operator >>(istream& in, Segment & m) {
	cout << endl << "Въведете точки за отсечката: " << endl;
	Point a;
	Point b;
	cin >> a >> b;
	Vector x(a, b);

	m.SetA(a);
	m.SetB(b);
	m.SetX(x);

	return in;
}
istream& operator >>(istream& in, Tetrahedron & m) {
	cout << "Въведете точките на тетраедъра: " << endl;
	Point A, B, C, D;
	cin >> A >> B >> C >> D;

	m.SetPointA(A);
	m.SetPointB(B);
	m.SetPointC(C);
	m.SetPointD(D);

	return in;
}



float angle(Point x, Point y, Point z) {

	Vector AB(x, y);
	Vector BC(y, z);
	Vector AC(x, z);

	double a = BC.VectorLen();
	double b = AC.VectorLen();
	double c = AB.VectorLen();

	double A = (pow(b, 2) + pow(c, 2) - pow(a, 2)) / (2 * c * b);


	float angle = acos(A) * (180.0 / 3.14159);

	return angle;
}





//CLASS POINT
Point::Point() {
	this->x = 0;
	this->y = 0;
	this->z = 0;

}

Point::Point(double x, double y, double z) {
	this->x = x;
	this->y = y;
	this->z = z;
}

Point::Point(const Point & s) {
	this->x = s.x;
	this->y = s.y;
	this->z = s.z;
}

void Point::SetX(double x) {
	this->x = x;
}
void Point::SetY(double y) {
	this->y = y;
}
void Point::SetZ(double z) {
	this->z = z;
}

double Point::GetX() const {
	return this->x;
}
double Point::GetY() const {
	return this->y;
}
double Point::GetZ() const {
	return this->z;
}



Point& Point::operator =(const Point& s) {
	this->x = s.x;
	this->y = s.y;
	this->z = s.z;
	return *this;
}


bool Point::operator == (const Point& s) const {
	if (this->x == s.x && this->y == s.y && this->z == s.z) {
		return true;
	}
	else return false;
}


bool operator < (const Point& t, const Triangle& s) {

	Triangle ABC(s.GetA(), s.GetB(), s.GetC());
	Triangle ABP(s.GetA(), s.GetB(), t);
	Triangle ACP(s.GetA(), s.GetC(), t);
	Triangle BCP(s.GetB(), s.GetC(), t);

	float AreaABC = ABC.TriangleArea();
	float AreaABP = ABP.TriangleArea();
	float AreaACP = ACP.TriangleArea();
	float AreaBCP = BCP.TriangleArea();

	if (AreaABP == 0 || AreaACP == 0 || AreaBCP == 0) {
		return false;
	}

	else if(round((AreaABP + AreaACP + AreaBCP) * 100) / 100 == round((AreaABC) * 100) / 100) {
		return true;
	}

	else return false;
}


bool operator > (const Point& t, const Triangle& s) {
	Triangle ABC(s.GetA(), s.GetB(), s.GetC());
	Triangle ABP(s.GetA(), s.GetB(), t);
	Triangle ACP(s.GetA(), s.GetC(), t);
	Triangle BCP(s.GetB(), s.GetC(), t);

	float AreaABC = ABC.TriangleArea();
	float AreaABP = ABP.TriangleArea();
	float AreaACP = ACP.TriangleArea();
	float AreaBCP = BCP.TriangleArea();


	if (AreaABP == 0 || AreaACP == 0 || AreaBCP == 0) {
		return false;
	}

	else if (round((AreaABP + AreaACP + AreaBCP) * 100) / 100 == round((AreaABC) * 100) / 100) {
		return false;
	}

	else return true;

}


bool operator == (const Point& t, const Triangle& s) {

	Triangle ABC(s.GetA(), s.GetB(), s.GetC());
	Triangle ABP(s.GetA(), s.GetB(), t);
	Triangle ACP(s.GetA(), s.GetC(), t);
	Triangle BCP(s.GetB(), s.GetC(), t);

	float AreaABC = round(ABC.TriangleArea()*10)/10;
	float AreaABP = round(ABP.TriangleArea()*10)/10;
	float AreaACP = round(ACP.TriangleArea()*10)/10;
	float AreaBCP = round(BCP.TriangleArea()*10)/10;

	
	if (AreaABP == 0 || AreaACP == 0 || AreaBCP == 0) {
		return true;
	}

	else return false;
}






//CLASS VECTOR
Vector::Vector() {
	SetX(0);
	SetY(0);
	SetZ(0);
}

Vector::Vector(double x, double y, double z) {
	SetX(x);
	SetY(y);
	SetZ(z);
}

Vector::Vector(Point S, Point T) {
	SetX(T.GetX() - S.GetX());
	SetY(T.GetY() - S.GetY());
	SetZ(T.GetZ() - S.GetZ());
}

Vector::Vector(const Vector& s) {
	SetX(s.GetX());
	SetY(s.GetY());
	SetZ(s.GetZ());
}

Vector& Vector::operator = (const Vector& s) {
	SetX(s.GetX());
	SetY(s.GetY());
	SetZ(s.GetZ());
	return *this;
}


double Vector::VectorLen() const{

	double lenght = sqrt(pow(GetX(), 2) + pow(GetY(), 2) + pow(GetZ(), 2));
	return lenght;
}


Vector Vector::VectorDirection() {
	Vector A(GetX() / VectorLen(), GetY() / VectorLen(), GetZ() / VectorLen());
	return A;
}


Vector Vector::VectorProjection(const Vector& s) {

	double mult = (this->GetX()*s.GetX()) + (this->GetY()*s.GetY()) + (this->GetZ()*s.GetZ());
	double len = pow(s.VectorLen(),2);
	double result = mult / len;

	Vector A(result*s.GetX(), result*s.GetY(), result*s.GetZ());

	return A;
}


Vector Vector::CrossProduct(const Vector& s) {
	Vector A( GetY()*s.GetZ() - GetZ()*s.GetY(), GetZ()*s.GetX() - GetX()*s.GetZ(), GetX()*s.GetY() - GetY()*s.GetX() );
	return A;
}


bool Vector::NullVector() {
	if (GetX() == 0 && GetY() == 0 && GetZ() == 0) {
		return true;
	}
	else return false;
}


bool Vector::ParallelVector(const Vector& s) {

	double check1 = this->GetX() / s.GetX();
	double check2 = this->GetY() / s.GetY();
	double check3 = this->GetZ() / s.GetZ();

	if (check1 == check2 && check1 == check3 && check2 == check3) {
		return true;
	}
	else return false;
}


bool Vector::PerpendicularVector(const Vector& s) {

	double a = GetX()*s.GetX();
	double b = GetY()*s.GetY();
	double c = GetZ()*s.GetZ();

	if ((a + b + c) == 0) {
		return true;
	}
	else return false;
}


bool operator ==(const Vector& s, const Vector& t) {
	if (s.GetX() == t.GetX() && s.GetY() == t.GetY() && s.GetZ() == t.GetZ()) {
		return true;
	}
	else return false;
}

Vector Vector::operator +(const Vector & s) {
	Vector A;
	A.SetX( GetX() + s.GetX() );
	A.SetY( GetY() + s.GetY() );
	A.SetZ( GetZ() + s.GetZ() );

	return A;
}	

Vector Vector::operator -(const Vector & s) {
	Vector A;
	A.SetX( GetX() - s.GetX() );
	A.SetY( GetY() - s.GetY() );
	A.SetZ( GetZ() - s.GetZ() );

	return A;
}

Vector Vector::operator *(double s) {
	Vector A;
	A.SetX( GetX() *s );
	A.SetY( GetY() *s );
	A.SetZ( GetZ() *s );

	return A;
}

double operator *(const Vector& t, const Vector & s) {

	double a = t.GetX() * s.GetX();
	double b = t.GetY() * s.GetY();
	double c = t.GetZ() * s.GetZ();

	return a+b+c;
}

Vector Vector::operator ^ (const Vector & s) {
	Vector A;
	double a = GetX();
	double b = GetY();
	double c = GetZ();

	A.SetX( b*s.GetZ() - c*s.GetY() );
	A.SetY(-a*s.GetZ() + c*s.GetX() );
	A.SetZ( a*s.GetY() - b*s.GetX() );

	return A;
}

double Vector::operator ()(const Vector& s, const Vector& t) {
	double arr[3][3];
	arr[0][0] = GetX();
	arr[0][1] = GetY();
	arr[0][2] = GetZ();

	arr[1][0] = s.GetX();
	arr[1][1] = s.GetY();
	arr[1][2] = s.GetZ();

	arr[2][0] = t.GetX();
	arr[2][1] = t.GetY();
	arr[2][2] = t.GetZ();
	
	double r1, r2, r3;

	r1 = arr[0][0] * ((arr[1][1] * arr[2][2])- (arr[2][1] * arr[1][2]));

	r2 = arr[0][1] * ((arr[1][0] * arr[2][2])- (arr[2][0] * arr[1][2]));

	r3 = arr[0][2] * ((arr[1][0] * arr[2][1])- (arr[2][0] * arr[1][1]));

	double detval = r1 - r2 + r3;

	return -(detval);
}



//CLASS TRIANGLE
Triangle::Triangle() {
	Point a(0, 0, 0);
	this->A = a;
	this->B = a;
	this->C = a;
}

Triangle::Triangle(Point A, Point B, Point C) {
	this->A = A;
	this->B = B;
	this->C = C;
}

Triangle::Triangle(const Triangle & s) {
	this->A = s.A;
	this->B = s.B;
	this->C = s.C;
}

Triangle Triangle::SetPointA(const Point& A) {
	this->A = A;
	return *this;
}
Triangle Triangle::SetPointB(const Point& B) {
	this->B = B;
	return *this;
}
Triangle Triangle::SetPointC(const Point& C) {
	this->C = C;
	return *this;
}


Point Triangle::GetA() const {
	return A;
}
Point Triangle::GetB() const {
	return B;
}
Point Triangle::GetC() const {
	return C;
}


Triangle& Triangle::operator = (const Triangle& s) {
	this->A = s.A;
	this->B = s.B;
	this->C = s.C;
	return *this;
}


void Triangle::TriangleType() {

	Vector AB(this->A, this->B);
	Vector AC(this->A, this->C);
	Vector BC(this->B, this->C);

	bool side1 = (AB.VectorLen() == AC.VectorLen() && AB.VectorLen() != BC.VectorLen());
	bool side2 = (AB.VectorLen() == BC.VectorLen() && AB.VectorLen() != AC.VectorLen());
	bool side3 = (AC.VectorLen() == BC.VectorLen() && AC.VectorLen() != AB.VectorLen());

	double angleA = round(angle(A, B, C)*100)/100;
	double angleB = round(angle(B, A, C)*100)/100;
	double angleC = round(angle(C, A, B)*100)/100;

	
	if ((side1 == true || side2 == true || side3 == true) && (angleA == 90.0 || angleB == 90.0 || angleC == 90.0)) {
		cout << "Триъгълникът е правоъгълен равнобедрен." << endl;
	}
	else if ((side1 == true || side2 == true || side3 == true) && (angleA > 90.0 || angleB > 90.0 || angleC > 90.0)) {
		cout << "Триъгълникът е тъпоъгълен равнобедрен." << endl;
	}
	else if ((side1 == true || side2 == true || side3 == true) && (angleA < 90.0 && angleB < 90.0 && angleC < 90.0)) {
		cout << "Триъгълникът е остроъгълен равнобедрен." << endl;
	}


	else if (AB.VectorLen() == BC.VectorLen() && AB.VectorLen() == AC.VectorLen()) {
		cout << "Триъгълникът е равностранен." << endl;
	}


	else if (angleA == 90.0 || angleB == 90.0 || angleC == 90.0) {
		cout << "Триъгълникът е правоъгълен." << endl;
	}
	else if (angleA > 90.0 || angleB > 90.0 || angleC > 90.0) {
		cout << "Триъгълникът е тъпоъгълен." << endl;
	}
	else if (angleA < 90.0 && angleB < 90.0 && angleC < 90.0) {
		cout << "Триъгълникът е остроъгълен." << endl;
	}

}
float Triangle::TriangleArea() {
	Vector AB(this->A, this->B);
	Vector AC(this->A, this->C);
	Vector BC(this->B, this->C);

	float semiper = (AB.VectorLen() + AC.VectorLen() + BC.VectorLen()) / 2;
	float area = sqrt( semiper * (semiper - AB.VectorLen()) * (semiper - AC.VectorLen()) * (semiper - BC.VectorLen()) );

	return area;
}

float Triangle::TrianglePerimeter() {

	Vector AB(this->A, this->B);
	Vector AC(this->A, this->C);
	Vector BC(this->B, this->C);

	float perimeter = AB.VectorLen() + AC.VectorLen() + BC.VectorLen();

	return perimeter;
}

Point Triangle::Medicenter() {
	Point M( ( (this->A.GetX() + this->B.GetX() + this->C.GetX()) )/ 3, \
		     ( (this->A.GetY() + this->B.GetY() + this->C.GetY()) )/ 3, \
		     ( (this->A.GetZ() + this->B.GetZ() + this->C.GetZ()) )/ 3  );

	return M;
}



//CLASS LINE

Line::Line() {
	A.SetX(0);
	A.SetY(0);
	A.SetZ(0);
	B.SetX(0);
	B.SetY(0);
	B.SetZ(0);
	X.SetX(0);
	X.SetY(0);
	X.SetZ(0);
}
Line::Line(Point &a, Point &b) {
	A.SetX(a.GetX());
	A.SetY(a.GetY());
	A.SetZ(a.GetZ());

	B.SetX(b.GetX());
	B.SetY(b.GetY());
	B.SetZ(b.GetZ());

	X.SetX(b.GetX() - a.GetX());
	X.SetY(b.GetY() - a.GetY());
	X.SetZ(b.GetZ() - a.GetZ());
}
Line::Line(Point &a, Vector &b) {
	A.SetX(a.GetX());
	A.SetY(a.GetY());
	A.SetZ(a.GetZ());

	X.SetX(b.GetX());
	X.SetY(b.GetY());
	X.SetZ(b.GetZ());

	B.SetX(a.GetX() + X.GetX());
	B.SetY(a.GetY() + X.GetY());
	B.SetZ(a.GetZ() + X.GetZ());
}
Line::Line(Line &s) {
	A = s.GetA();
	B = s.GetB();
	X = s.GetX();
}


void Line::SetA(const Point& A) {
	this->A = A;
}
void Line::SetB(const Point& B) {
	this->B = B;
}
void Line::SetX(const Vector& X) {
	this->X = X;
}


Point Line::GetA() const {
	return this->A;
}
Point Line::GetB() const {
	return this->B;
}
Vector Line::GetX() const{
	return this->X;
}


Vector Line::LineDirection() {
	return this->X;
}


Vector Line::NormalVector() {
	double l = X.GetX();
	double m = X.GetY();
	double n = X.GetZ();
	double x1 = A.GetX();
	double y1 = A.GetY();
	double z1 = A.GetZ();

	double k = (l*x1 + m*y1 + n*z1) / ( pow(l, 2) + pow(m, 2) + pow(n, 2) );

	Vector NormalVector(l*k + x1, m*k + y1, n*k + z1);

	return NormalVector;
}


double Line::LinesAngle(Line& b) {
	Vector u = GetX();
	Vector v = b.GetX();
	
	double dotprod = (u.GetX()*v.GetX()) + (u.GetY()*v.GetY()) + (u.GetZ()*v.GetZ());

	double lenprod = u.VectorLen()*v.VectorLen();

	double cosine = dotprod / lenprod;
	
	return acos(cosine);
}


Line& Line::operator = (const Line& s){
	this->A = s.A;
	this->B = s.B;
	this->X = s.X;
	return *this;
}


bool operator + (const Point& P, const Line& s) {
	Vector a(P, s.GetA());
	Vector b(P, s.GetB());

	float lenAP = round(a.VectorLen()*100)/100;
	float lenBP = round(b.VectorLen()*100)/100;
	float lenAB = round(s.GetX().VectorLen()*100)/100;


	if (lenAB == lenAP + lenBP || lenAP == lenBP + lenAB || lenBP == lenAP + lenAB) {
		return true;
	}
	else return false;
}


bool operator ||(const Line& t, const Line& s) {
	bool result = false;

	for (double N = 0.0625; N <= 10; N+=0.0625) {
		if (t.GetX()*(-N) == s.GetX() || s.GetX()*(-N) == t.GetX()) {
			result = true;
			break;
		}
		if (t.GetX()*(N) == s.GetX() || s.GetX()*(N) == t.GetX()) {
			result = true;
			break;
		}
	}
	return result;
}


bool operator ==(const Line& t, const Line& s) {
	
	if ((t.GetA() + s == true || t.GetB() + s == true) && (t||s) == true) {
		return true;
	}
	else return false;	
}


bool operator !=(const Line& t, const Line& s) {
	Vector A = t.GetX().CrossProduct(s.GetX());

	if (A.GetX() == 0 && A.GetY() == 0 && A.GetZ() == 0) {
		return false;
	}
	else return true;
}


bool operator &&(const Line& t, const Line& s) {

	Vector V(t.GetA().GetX() - s.GetA().GetX(), t.GetA().GetY() - s.GetA().GetY(), t.GetA().GetZ() - s.GetA().GetZ());
	Vector S = t.GetX().CrossProduct(s.GetX());

	float dotprod = V.GetX()*S.GetX() + V.GetY()*S.GetY() + V.GetZ()*S.GetZ();

	if (dotprod == 0) {
		return true;
	}
	else return false;
}


bool operator |(const Line& t, const Line& s) {
	if (t&&s) {
		float dotprod = t.GetX().GetX()*s.GetX().GetX() + t.GetX().GetY()*s.GetX().GetY() + t.GetX().GetZ()*s.GetX().GetZ();
		if (dotprod == 0) {
			return true;
		}
	}
	else return false;
}






//CLASS SEGMENT
Segment::Segment() {
	Point a(0, 0, 0);
	Vector b(a, a);
	this->SetA(a);
	this->SetB(a);
	this->SetX(b);
}
Segment::Segment(Point &a, Point &b) {
	this->SetA(a);
	this->SetB(b);
	Vector c(GetA(), GetB());
	SetX(c);
}
Segment::Segment(Segment& s) {
	Point a(s.GetA());
	Point b(s.GetB());
	Vector x(s.GetX());
	this->SetA(a);
	this->SetB(b);
	this->SetX(x);
}


float Segment::SegmentLen() {
	return GetX().VectorLen();
}


Point Segment::SegmentMid() {
	return Point((GetA().GetX() + GetB().GetX()) / 2, (GetA().GetY() + GetB().GetY()) / 2, (GetA().GetZ() + GetB().GetZ()) / 2);
}


bool operator ==(const Point& s,const Segment& t) {
	Vector AP(s, t.GetA());
	Vector BP(s, t.GetB());

	if (round(AP.VectorLen()*1000)/1000 + round(BP.VectorLen()*1000)/1000 == round(t.GetX().VectorLen() * 1000) / 1000) {
		return true;
	}
	else return false;
}	


Segment& Segment::operator = (const Segment& s) {
	SetA(s.GetA());
	SetB(s.GetB());
	SetX(s.GetX());
	return *this;
}




//CLASS TETRAHEDRON
Tetrahedron::Tetrahedron() {
	Point F(0, 0, 0);
	A = F;
	B = F;
	C = F;
	D = F;
}

Tetrahedron::Tetrahedron(Point a, Point b, Point c, Point d) {
	A = a;
	B = b;
	C = c;
	D = d;
}

Tetrahedron::Tetrahedron(const Tetrahedron& s) {
	this->A = s.A;
	this->B = s.B;
	this->C = s.C;
	this->D = s.D;
}


Tetrahedron Tetrahedron::SetPointA(const Point& A) {
	this->A = A;
	return *this;
}
Tetrahedron Tetrahedron::SetPointB(const Point& B) {
	this->B = B;
	return *this;
}
Tetrahedron Tetrahedron::SetPointC(const Point& C) {
	this->C = C;
	return *this;
}
Tetrahedron Tetrahedron::SetPointD(const Point& D) {
	this->D = D;
	return *this;
}


Point Tetrahedron::GetA() const {
	return A;
}
Point Tetrahedron::GetB() const {
	return B;
}
Point Tetrahedron::GetC() const {
	return C;
}
Point Tetrahedron::GetD() const {
	return D;
}


bool Tetrahedron::RegularTetrahedron() {
	Vector AB(A, B);
	Vector AC(A, C);
	Vector BC(B, C);
	Vector AD(A, D);
	Vector BD(B, D);
	Vector CD(C, D);

	if (AB.VectorLen() == AC.VectorLen() && AC.VectorLen() == BC.VectorLen() && AD.VectorLen() == CD.VectorLen() && CD.VectorLen() == BC.VectorLen()) {
		return true;
	}
	else return false;
}


bool Tetrahedron::OrthocentricTetrahedron() {
	Vector AB(A, B);
	Vector AC(A, C);
	Vector BC(B, C);
	Vector AD(A, D);
	Vector BD(B, D);
	Vector CD(C, D);

	if ((pow(AB.VectorLen(), 2) + pow(CD.VectorLen(), 2)) == (pow(AC.VectorLen(), 2) + pow(BD.VectorLen(), 2)) && (pow(AC.VectorLen(), 2) + pow(BD.VectorLen(), 2)) == (pow(AD.VectorLen(), 2) + pow(BC.VectorLen(), 2))) {
		return true;
	}
	else return false;
}	


double Tetrahedron::TetrahedronArea() {
	double area;
	if (this->RegularTetrahedron() == true) {
		Vector AB(A, B);
		area = (AB.VectorLen() / 4)*sqrt(3);
		return area;
	}
	else {
		Triangle ABC(A, B, C);
		Triangle ABD(A, B, D);
		Triangle ACD(A, C, D);
		Triangle BCD(B, C, D);
		area = ABC.TriangleArea() + ABC.TriangleArea() + ACD.TriangleArea() + BCD.TriangleArea();
		return area;
	}
}


double Tetrahedron::TetrahedronVolume() {
	double volume;
	Vector AB(A, B);
	Vector AC(A, C);
	Vector AD(A, D);
	if (this->RegularTetrahedron() == true) {
		volume = AB.VectorLen() / (6 * sqrt(2));
		return volume;
	}

	else {
		double arr[3][3];
		arr[0][0] = AB.GetX();
		arr[0][1] = AB.GetY();
		arr[0][2] = AB.GetZ();

		arr[1][0] = AC.GetX();
		arr[1][1] = AC.GetY();
		arr[1][2] = AC.GetZ();

		arr[2][0] = AD.GetX();
		arr[2][1] = AD.GetY();
		arr[2][2] = AD.GetZ();

		double r1, r2, r3;

		r1 = arr[0][0] * ((arr[1][1] * arr[2][2]) - (arr[2][1] * arr[1][2]));

		r2 = arr[0][1] * ((arr[1][0] * arr[2][2]) - (arr[2][0] * arr[1][2]));

		r3 = arr[0][2] * ((arr[1][0] * arr[2][1]) - (arr[2][0] * arr[1][1]));

		double detval = r1 - r2 + r3;

		cout << detval << ' ' << r1 << ' ' << r2 << ' ' << r3 << endl;
		return abs(detval/6);
	}
}


bool operator < (const Point& t, const Tetrahedron& s) {
	Tetrahedron S(s.GetA(), s.GetB(), s.GetC(), s.GetD());

	Triangle T1(s.GetA(), s.GetB(), s.GetC());
	Triangle T2(s.GetA(), s.GetB(), s.GetD());
	Triangle T3(s.GetA(), s.GetC(), s.GetD());
	Triangle T4(s.GetD(), s.GetC(), s.GetD());

	Tetrahedron S1(s.GetA(), s.GetB(), s.GetC(), t);
	Tetrahedron S2(s.GetA(), s.GetD(), s.GetC(), t);
	Tetrahedron S3(s.GetD(), s.GetB(), s.GetC(), t);
	Tetrahedron S4(s.GetA(), s.GetB(), s.GetD(), t);

	if ( t>T1== true && t>T2 == true && t>T3 == true && t>T4 == true && \
		(S1.TetrahedronVolume() + S2.TetrahedronVolume() + S3.TetrahedronVolume() + S4.TetrahedronVolume()) == S.TetrahedronVolume()) {

		return true;
	}	
	else return false;
}


bool operator > (const Point& t, const Tetrahedron& s) {
	Tetrahedron S(s.GetA(), s.GetB(), s.GetC(), s.GetD());

	Triangle T1(s.GetA(), s.GetB(), s.GetC());
	Triangle T2(s.GetA(), s.GetB(), s.GetD());
	Triangle T3(s.GetA(), s.GetC(), s.GetD());
	Triangle T4(s.GetD(), s.GetC(), s.GetD());

	Tetrahedron S1(s.GetA(), s.GetB(), s.GetC(), t);
	Tetrahedron S2(s.GetA(), s.GetD(), s.GetC(), t);
	Tetrahedron S3(s.GetD(), s.GetB(), s.GetC(), t);
	Tetrahedron S4(s.GetA(), s.GetB(), s.GetD(), t);

	if (t > T1 == true && t > T2 == true && t > T3 == true && t > T4 == true && \
		(S1.TetrahedronVolume() + S2.TetrahedronVolume() + S3.TetrahedronVolume() + S4.TetrahedronVolume()) > S.TetrahedronVolume()) {

		return true;
	}
	else return false;
}


bool operator == (const Point& t, const Tetrahedron& s) {

	Triangle T1(s.GetA(), s.GetB(), s.GetC());
	Triangle T2(s.GetA(), s.GetB(), s.GetD());
	Triangle T3(s.GetA(), s.GetC(), s.GetD());
	Triangle T4(s.GetD(), s.GetC(), s.GetD());

	if ( (t==T1) == true || (t == T2) == true || (t<T3) == true || (t == T4) == true) {
		return true;
	}
	else return false;
}




void ConsoleMenu() {
	char choice{ 'y' };

	while (choice == 'y') {

		cout << endl << "Изберете вид геометричен обект: " << endl << "1 - Точка" << endl << "2 - Вектор" << endl << "3 - Линия" << endl;
		cout << "4 - Отсечка" << endl << "5 - Триъгълник" << endl << "6 - Тетраедър" << endl << endl;

		int object;
		cin >> object;


		if (object == 1) { //ТОЧКА
			Point A;
			cin >> A;

			int PointOperation;
			cout << endl << "Моля изберете операция с точката :" << endl << "1 - Проверка дали две точки съвпадат." << endl;
			cin >> PointOperation;

			if (PointOperation == 1) {
				Point B;
				cin >> B;

				if (A == B) {
					cout << endl << "Точките съвпадат." << endl;
				}
				else cout << endl << "Точките не съвпадат" << endl;
			}
		}

		else if (object == 2) { //ВЕКТОР
			char contchoice{ 'y' };

			cout << endl << "Вектор А: ";
			Vector A;
			cin >> A;
			while (contchoice == 'y') {


				int VectorChoice;
				cout << endl << "Моля изберете операция с вектора :" << endl << "1 - Изчисляване на дължина на вектор" << endl \
					<< "2 - Изчисляване на посока на вектор" << endl << "3 - Проверка за нулев вектор" << endl \
					<< "4 - Проекция на вектор върху друг вектор" << endl << "5 - Проверка за успоредност на два вектора" << endl \
					<< "6 - Проверка за перпендикулярност на два вектора" << endl << "7 - Събиране на два вектора" << endl \
					<< "8 - Умножение на вектор с реално число" << endl << "9 - Скаларно произведение на два вектора" << endl \
					<< "10 - Векторно произведение на два вектора" << endl << "11 - Смесено произведение на три вектора" << endl << endl;
				cin >> VectorChoice;


				if (VectorChoice == 1) {
					cout << endl << "Векторът е с дължина: " << A.VectorLen() << endl;
				}


				else if (VectorChoice == 2) {

					if (A.VectorLen() == 0) {
						try {
							ThrowVectorExcetpion();
						}
						catch (VectorLenghtExeption & e) {
							cout << endl << e.getMessage() << endl;
						}
					}
					else cout << endl << "Посоката на вектора е: " << A.VectorDirection() << endl;
				}


				else if (VectorChoice == 3) {

					if (A.NullVector() == true) {
						cout << endl << "Векторът е нулев." << endl;
					}
					else cout << endl << "Векторът не е нулев." << endl;

				}


				else if (VectorChoice == 4) {
					cout << endl << "Вектор B: ";
					Vector B;
					cin >> B;

					if (A.VectorLen() == 0 || B.VectorLen() == 0) {
						try {
							ThrowVectorExcetpion();
						}
						catch (VectorLenghtExeption & e) {
							cout << endl << e.getMessage() << endl;
						}
					}

					else cout << endl << "Проекцията на вектор А върху вектор В е: " << B.VectorProjection(A) << endl;
				}


				else if (VectorChoice == 5) {
					cout << endl << "Вектор B: ";
					Vector B;
					cin >> B;

					if (A.VectorLen() == 0 || B.VectorLen() == 0) {
						try {
							ThrowVectorExcetpion();
						}
						catch (VectorLenghtExeption & e) {
							cout << endl << e.getMessage() << endl;
						}
					}

					else {
						if (A.ParallelVector(B) == true) {
							cout << endl << "Векторите са успоредни." << endl;
						}
						else cout << endl << "Векторите не са успоредни." << endl;
					}
				}


				else if (VectorChoice == 6) {
					cout << endl << "Вектор B: ";
					Vector B;
					cin >> B;

					if (A.VectorLen() == 0 || B.VectorLen() == 0) {
						try {
							ThrowVectorExcetpion();
						}
						catch (VectorLenghtExeption & e) {
							cout << endl << e.getMessage() << endl;
						}
					}

					else {
						if (A.PerpendicularVector(B) == true) {
							cout << endl << "Векторите са перпендикулярни." << endl;
						}
						else cout << endl << "Векторите не са перпендикулярни." << endl;
					}
				}


				else if (VectorChoice == 7) {
					cout << endl << "Вектор B: ";
					Vector B;
					cin >> B;

					cout << endl << A << " + " << B << " = " << A + B << endl;
				}


				else if (VectorChoice == 8) {
					double b;
					cout << endl << "Въведете реално число: ";
					cin >> b;

					cout << endl << A << " * " << b << " = " << A * b << endl;
				}


				else if (VectorChoice == 9) {
					cout << endl << "Вектор B: ";
					Vector B;
					cin >> B;

					cout << endl << A << " * " << B << " = " << A * B << endl;
				}


				else if (VectorChoice == 10) {
					cout << endl << "Вектор B: ";
					Vector B;
					cin >> B;

					cout << endl << A << " ^ " << B << " = " << (A ^ B) << endl;
				}


				else if (VectorChoice == 11) {
					cout << endl << "Вектор B: ";
					Vector B;
					cin >> B;

					cout << endl << "Вектор C: ";
					Vector C;
					cin >> C;

					cout << endl << "(" << A << " x " << B << ") " << " . " << C << " = " << A(C, B) << endl;
				}


				cout << endl << "Желаете ли да изберете друга операция? - y/n" << endl;
				cin >> contchoice;

			}
		}

		else if (object == 3) { // ЛИНИЯ

			char contchoice{ 'y' };

			Line A;
			cout << endl << "Изберете начин за инициализация на линията: 1 - две точки, 2 - точка и вектор." << endl;

			int initialize;
			cin >> initialize;

			if (initialize == 1) {
				cin >> A;
			}

			else if (initialize == 2) {
				Point B;
				Vector C;
				cin >> B >> C;
				A = Line(B, C);
			}


			while (contchoice == 'y') {
				cout << endl << "Изберете операция с линията: " << endl << "1 - Намиране на посока на правата." << endl \
					<< "2 - Намиране на нормален вектор." << endl << "3 - Намиране на ъгъл между две прави." << endl \
					<< "4 - Проверка дали дадена точка лежи на дадената права." << endl \
					<< "5 - Проверка дали дадена права е успоредна на друга права." << endl \
					<< "6 - Проверка дали права съвпада с друга права." << endl \
					<< "7 - Проверка дали права е кръстосана с друга права." << endl << "8 - Проверка дали две прави се пресичат." << endl \
					<< "9 - Проверка дали две прави са перпендикулярни. " << endl << endl;
					

				int LineChoice;
				cin >> LineChoice;



				if (LineChoice == 1) {
					cout << endl << "Посоката на правата е: " << A.GetX() << endl;
				}

				else if (LineChoice == 2) {
					cout << endl << "Нормален вектор е : " << A.NormalVector() << endl;
				}

				else if (LineChoice == 3) {
					Line B;
					cin >> B;
					cout << endl << "Ъгълът между правите е: " << A.LinesAngle(B) << "rad " << endl;
				}

				else if (LineChoice == 4) {
					cout << "Точка А: " << endl;
					Point B;
					cin >> B;
					if ((B + A) == true) {
						cout << endl << "Точката лежи на правата." << endl;
					}
					else cout << endl << "Точката не лежи на правата" << endl;
				}

				else if (LineChoice == 5) {
					Line B;
					cin >> B;
					if ((A || B) == true) {
						cout << endl << "Правите са успоредни." << endl;
					}
					else cout << endl << "Правите не са успоредни." << endl;
				}

				else if (LineChoice == 6) {
					Line B;
					cin >> B;
					if ((A == B) == true) {
						cout << endl << "Правите съвпадат." << endl;
					}
					else cout << endl << "Правите не съвпадат." << endl;
				}

				else if (LineChoice == 7) {
					Line B;
					cin >> B;
					if ((A != B) == true) {
						cout << endl << "Правите са кръстосани." << endl;
					}
					else cout << endl << "Правите не са кръстосани." << endl;
				}

				else if (LineChoice == 8) {
					Line B;
					cin >> B;
					if ((A && B) == true) {
						cout << endl << "Правите се пресичат." << endl;
					}
					else cout << endl << "Правите не се пресичат." << endl;
				}

				else if (LineChoice == 9) {
					Line B;
					cin >> B;
					if ((A | B) == true) {
						cout << endl << "Правите са перпендикулярни." << endl;
					}
					else cout << endl << "Правите не са перпендикулярни." << endl;
				}

				cout << endl << "Желаете ли да изберете нова операция? - y/n" << endl;
				cin >> contchoice;
			}
		}

		else if (object == 4) {

			Segment A;
			cin >> A;

			char contchoice{ 'y' };
			int SegmentOperation;

			while (contchoice == 'y') {

				cout << endl << "Моля изберете операция с отсечката: " << endl << "1 - Намиране на дължина на отсечка." << endl\
					<< "2 - Намиране на среда на отсечка." << endl << "3 - Проверка дали точка лежи на отсечката" << endl << endl;
				cin >> SegmentOperation;

				if (SegmentOperation == 1) {
					cout << endl << "Дължината на отсечката е : " << A.SegmentLen() << endl;
				}

				else if (SegmentOperation == 2) {
					cout << endl << "Средата на отсечката е : " << A.SegmentMid() << endl;
				}

				else if (SegmentOperation == 3) {
					Point B;
					cin >> B;
					if ((B == A) == true) {
						cout << endl << "Точката лежи на правата. " << endl;
					}
					else cout << endl << "Точката не лежи на правата. " << endl;
				}

				cout << endl << "Желаете ли да изберете друга опрецария? - y/n" << endl;
				cin >> contchoice;
			}
		}

		else if (object == 5) {
			Triangle A;
			cin >> A;

			bool except = false;

			try {
				if (A.GetA() == A.GetB()) {
					throw EqualPointException(A.GetA(), A.GetB());
				}
				else if (A.GetB() == A.GetC()) {
					throw EqualPointException(A.GetB(), A.GetC());
				}
				else if (A.GetA() == A.GetC()) {
					throw EqualPointException(A.GetA(), A.GetC());
				}
			}
			catch (EqualPointException & e) {
				cout << endl;
				e.getMessage();
				cout << endl;
				except = true;
			}


			char contchoice{ 'y' };
			int TriangleOperation;

			while (contchoice == 'y') {

				if (except == true) break;


				cout << endl << "Изберете операция с триъгълника: " << endl << "1 - Определяне на вида на триъгълника." \
					<< endl << "2 - Намиране на лицето на триъгълника." << endl << "3 - Намиране на периметъра на триъгълника." \
					<< endl << "4 - Намиране на медицентъра на триъгълника." << endl \
					<< "5 - Проверка дали точка лежи в триъгълника." << endl << "6 - Проверка дали точка е извън триъгълника." \
					<< endl << "7 - Проверка дали точка лежи на някоя от страните на триъгълника." << endl << endl;

				cin >> TriangleOperation;

				if (TriangleOperation == 1) {
					cout << endl;
					A.TriangleType();
					cout << endl;
				}

				else if (TriangleOperation == 2) {
					cout << endl << "Лицето на триъгълника е: " << A.TriangleArea() << endl;
				}

				else if (TriangleOperation == 3) {
					cout << endl << "Периметъра на триъгълника е: " << A.TrianglePerimeter() << endl;
				}

				else if (TriangleOperation == 4) {
					cout << endl << "Медицентъра на триъгълника е: " << A.Medicenter() << endl;
				}


				else if (TriangleOperation == 5) {
					Point P;
					cin >> P;

					if ((P < A) == true) {
						cout << endl << "Точката лежи в триъгълника. " << endl;
					}
					else cout << endl << "Точката не лежи в триъгълника. " << endl;
				}

				else if (TriangleOperation == 6) {
					Point P;
					cin >> P;

					if ((P > A) == true) {
						cout << endl << "Точката лежи извън триъгълника. " << endl;
					}
					else cout << endl << "Точката не лежи в триъгълника. " << endl;
				}

				else if (TriangleOperation == 7) {
					Point P;
					cin >> P;

					if ((P == A) == true) {
						cout << endl << "Точката лежи на една от страните в триъгълника. " << endl;
					}
					else cout << endl << "Точката не лежи на някоя от страните в триъгълника. " << endl;
				}

				cout << endl << "Желаете ли да изберете нова операция? - y/n." << endl;
				cin >> contchoice;
			}
		}

		else if (object == 6) {
			Tetrahedron A;
			cin >> A;

			bool except = false;

			try {
				if (A.GetA() == A.GetB()) {
					throw EqualPointException(A.GetA(), A.GetB());
				}
				else if (A.GetB() == A.GetC()) {
					throw EqualPointException(A.GetB(), A.GetC());
				}
				else if (A.GetA() == A.GetC()) {
					throw EqualPointException(A.GetA(), A.GetC());
				}
				else if (A.GetA() == A.GetD()) {
					throw EqualPointException(A.GetA(), A.GetD());
				}
				else if (A.GetD() == A.GetC()) {
					throw EqualPointException(A.GetD(), A.GetC());
				}
				else if (A.GetB() == A.GetD()) {
					throw EqualPointException(A.GetB(), A.GetD());
				}
			}
			catch (EqualPointException & e) {
				cout << endl;
				e.getMessage();
				cout << endl;
				except = true;
			}


			char contchoice{ 'y' };
			int TetrahedronOperation;

			while (contchoice == 'y') {

				if (except == true) break;


				cout << endl << "Изберете операция с тетраедъра: " << endl << "1 - Проверка дали тетраедъра е правилен. " \
					<< endl << "2 - Проверка дали тетраедъра е ортогонален." << endl << "3 - Намиране на околна повърхнина на тетраедъра." \
					<< endl << "4 - Намиране на обем на тетраедъра." << endl \
					<< "5 - Проверка дали точка лежи в тетраедъра." << endl << "6 - Проверка дали точка е извън тетраедъра." \
					<< endl << "7 - Проверка дали точка лежи на някоя от страните на тетраедъра." << endl << endl;

				cin >> TetrahedronOperation;

				if (TetrahedronOperation == 1) {

					if (A.RegularTetrahedron() == true) {
						cout << endl << "Тетраедърът е правилен." << endl;
					}
					else cout << endl << "Тетраедърът не е правилен." << endl;

				}

				else if (TetrahedronOperation == 2) {
					if (A.OrthocentricTetrahedron() == true) {
						cout << endl << "Тетраедърът е ортогонален." << endl;
					}
					else cout << endl << "Тетраедърът не е ортогонален." << endl;
				}

				else if (TetrahedronOperation == 3) {
					cout << endl << "Околната повърхнина на тетраедъра е " << A.TetrahedronArea() << endl;
				}

				else if (TetrahedronOperation == 4) {
					cout << endl << "Обемът на тетраедъра е " << A.TetrahedronVolume() << endl;
				}

				else if (TetrahedronOperation == 5) {
					Point B;
					cin >> B;

					if (B < A) {
						cout << endl << "Точката лежи в тетраедъра." << endl;
					}
					else cout << endl << "Точката не лежи в тетраедъра." << endl;
				}

				else if (TetrahedronOperation == 6) {
					Point B;
					cin >> B;

					if (B > A) {
						cout << endl << "Точката лежи в извън тетраедъра." << endl;
					}
					else cout << endl << "Точката не лежи извън тетраедъра." << endl;
				}

				else if (TetrahedronOperation == 7) {
					Point B;
					cin >> B;

					if (B == A) {
						cout << endl << "Точката лежи на една от страните на тетраедъра." << endl;
					}
					else cout << endl << "Точката не лежи на никоя от страните на тетраедъра." << endl;
				}

				cout << endl << "Желаете ли да изберете друга операция? - y/n." << endl;
				cin >> contchoice;
			}
		}

		else {
			cout << endl << "Грешна операция, въведете отново." << endl;
			cin >> object;
		}

		cout << endl << "Желаете ли да изберете друг обект? - y/n" << endl;
		cin >> choice;
	}
}


void FileRead() {

	ifstream input_data;
	input_data.open("D:\\Vectors.txt");

	int object;
	int operation;

	input_data >> object;
	
	if (object == 1) {
		double x, y, z;

		input_data >> x >> y >> z;
		Point A(x, y, z);
		
		input_data >> operation;

		if (operation == 1) {
			double x1, y1, z1;

			input_data >> x1 >> y1 >> z1;
			Point B(x1, y1, z1);

			if (A == B) {
				cout << endl << "Точките " << A << " и " << B << " съвпадат." << endl;
			}
			else cout << endl << "Точките " << A << " и " << B << " не съвпадат" << endl;
		}

	}

	else if (object == 2) {
		double x, y, z;
		double x1, y1, z1;

		input_data >> x >> y >> z;
		Vector A(x, y, z);

		input_data >> operation;

		if (operation == 1) {
			cout << endl << "Векторът " << A << " е с дължина: " << A.VectorLen() << endl;
		}

		else if (operation == 2) {
			if (A.VectorLen() == 0) {
				try {
					ThrowVectorExcetpion();
				}
				catch (VectorLenghtExeption & e) {
					cout << endl << e.getMessage() << endl;
				}
			}
			else cout << endl << "Посоката на вектора " << A << " е: " << A.VectorDirection() << endl;
		}

		else if (operation == 3) {
			if (A.NullVector() == true) {
				cout << endl << "Векторът " << A << " е нулев." << endl;
			}
			else cout << endl << "Векторът " << A << " не е нулев." << endl;
		}

		else if (operation == 4) {

			input_data >> x1 >> y1 >> z1;
			Vector B(x1, y1, z1);

			if (A.VectorLen() == 0 || B.VectorLen() == 0) {
				try {
					ThrowVectorExcetpion();
				}
				catch (VectorLenghtExeption & e) {
					cout << endl << e.getMessage() << endl;
				}
			}
			else cout << endl << "Проекцията на вектор " << A << " върху вектор " << B << " е: " << B.VectorProjection(A) << endl;
		
		}

		else if (operation == 5) {

			input_data >> x1 >> y1 >> z1;
			Vector B(x1, y1, z1);

			if (A.VectorLen() == 0 || B.VectorLen() == 0) {
				try {
					ThrowVectorExcetpion();
				}
				catch (VectorLenghtExeption & e) {
					cout << endl << e.getMessage() << endl;
				}
			}

			else {
				if (A.ParallelVector(B) == true) {
					cout << endl << "Векторите " << A << " и " << B << " са успоредни." << endl;
				}
				else cout << endl << "Векторите " << A << " и " << B << " не са успоредни." << endl;
			}
		}

		else if (operation == 6) {


			input_data >> x1 >> y1 >> z1;
			Vector B(x1, y1, z1);

			if (A.VectorLen() == 0 || B.VectorLen() == 0) {
				try {
					ThrowVectorExcetpion();
				}
				catch (VectorLenghtExeption & e) {
					cout << endl << e.getMessage() << endl;
				}
			}

			else {
				if (A.PerpendicularVector(B) == true) {
					cout << endl << "Векторите " << A << " и " << B << " са перпендикулярни." << endl;
				}
				else cout << endl << "Векторите " << A << " и " << B << " не са перпендикулярни." << endl;
			}
		}

		else if (operation == 7) {


			input_data >> x1 >> y1 >> z1;
			Vector B(x1, y1, z1);

			cout << endl << A << " + " << B << " = " << A + B << endl;
		}

		else if (operation == 8) {
			double num;
			input_data >> num;

			cout << endl << A << " * " << num << " = " << A * num << endl;
		}

		else if (operation == 9) {

			input_data >> x1 >> y1 >> z1;
			Vector B(x1, y1, z1);

			cout << endl << A << " * " << B << " = " << A * B << endl;
		}

		else if (operation == 10) {

			input_data >> x1 >> y1 >> z1;
			Vector B(x1, y1, z1);

			cout << endl << A << " ^ " << B << " = " << (A ^ B) << endl;
		}

		else if (operation == 11) {

		input_data >> x1 >> y1 >> z1;
		Vector B(x1, y1, z1);

		input_data >> x1 >> y1 >> z1;
		Vector C(x1, y1, z1);

		cout << endl << "(" << A << " x " << B << ") " << " . " << C << " = " << A(C, B) << endl;
		}
	}

	else if (object == 3) {
		double x, y, z;
		double x1, y1, z1;

		input_data >> x >> y >> z;
		Point B(x, y, z);

		input_data >> x1 >> y1 >> z1;
		Point C(x1, y1, z1);

		Line A(B, C);
		
		input_data >> operation;

		if (operation == 1) {
			cout << endl << "Посоката на правата е: " << A.GetX() << endl;
		}

		else if (operation == 2) {
			cout << endl << "Нормален вектор е : " << A.NormalVector() << endl;
		}

		else if (operation == 3) {

			input_data >> x >> y >> z;
			Point S(x, y, z);

			input_data >> x1 >> y1 >> z1;
			Point T(x1, y1, z1);

			Line U(S, T);

			cout << endl << "Ъгълът между правите е: " << A.LinesAngle(U) << "rad " << endl;
		}

		else if (operation == 4) {

			input_data >> x >> y >> z;
			Point B(x, y, z);

			if ((B + A) == true) {
				cout << endl << "Точката " << B << " лежи на правата." << endl;
			}
			else cout << endl << "Точката " << B << " не лежи на правата" << endl;
		}

		else if (operation == 5) {

			input_data >> x >> y >> z;
			Point S(x, y, z);

			input_data >> x1 >> y1 >> z1;
			Point T(x1, y1, z1);

			Line U(S, T);

			if ((U || A) == true) {
				
				cout << endl << "Правите са успоредни." << endl;
			}
			else cout << endl << "Правите не са успоредни." << endl;
		}

		else if (operation == 6) {

			input_data >> x >> y >> z;
			Point S(x, y, z);

			input_data >> x1 >> y1 >> z1;
			Point T(x1, y1, z1);

			Line U(S, T);


			if ((A == U) == true) {
				cout << endl << "Правите съвпадат." << endl;
			}
			else cout << endl << "Правите не съвпадат." << endl;
		}

		else if (operation == 7) {

			input_data >> x >> y >> z;
			Point S(x, y, z);

			input_data >> x1 >> y1 >> z1;
			Point T(x1, y1, z1);

			Line U(S, T);

			if ((A != U) == true) {
				cout << endl << "Правите са кръстосани." << endl;
			}
			else cout << endl << "Правите не са кръстосани." << endl;
		}

		else if (operation == 8) {

			input_data >> x >> y >> z;
			Point S(x, y, z);

			input_data >> x1 >> y1 >> z1;
			Point T(x1, y1, z1);

			Line U(S, T);

			if ((A && U) == true) {
				cout << endl << "Правите се пресичат." << endl;
			}
			else cout << endl << "Правите не се пресичат." << endl;
		}

		else if (operation == 9) {

			input_data >> x >> y >> z;
			Point S(x, y, z);

			input_data >> x1 >> y1 >> z1;
			Point T(x1, y1, z1);

			Line U(S, T);

			if ((A | U) == true) {
				cout << endl << "Правите са перпендикулярни." << endl;
			}
			else cout << endl << "Правите не са перпендикулярни." << endl;
		}
	}

	else if (object == 4) {
		double x, y, z;
		double x1, y1, z1;

		input_data >> x >> y >> z;
		Point B(x, y, z);

		input_data >> x1 >> y1 >> z1;
		Point C(x1, y1, z1);

		Segment A(B, C);

		input_data >> operation;

		if (operation == 1) {
			cout << endl << "Дължината на отсечката е : " << A.SegmentLen() << endl;
		}

		else if (operation == 2) {
			cout << endl << "Средата на отсечката е : " << A.SegmentMid() << endl;
		}

		else if (operation == 3) {
			input_data >> x >> y >> z;
			Point B(x, y, z);

			if ((B == A) == true) {
				cout << endl << "Точката лежи на правата. " << endl;
			}
			else cout << endl << "Точката не лежи на правата. " << endl;
		}
	}

	else if (object == 5) {
		double x, y, z;

		input_data >> x >> y >> z;
		Point A(x, y, z);
		input_data >> x >> y >> z;
		Point B(x, y, z);
		input_data >> x >> y >> z;
		Point C(x, y, z);

		Triangle ABC(A, B, C);
		cout << ABC;
		bool except = false;

		try {
			if (ABC.GetA() == ABC.GetB()) {
				throw EqualPointException(ABC.GetA(), ABC.GetB());
			}
			else if (ABC.GetB() == ABC.GetC()) {
				throw EqualPointException(ABC.GetB(), ABC.GetC());
			}
			else if (ABC.GetA() == ABC.GetC()) {
				throw EqualPointException(ABC.GetA(), ABC.GetC());
			}
		}
		catch (EqualPointException & e) {
			cout << endl;
			e.getMessage();
			cout << endl;
			except = true;
		}

		input_data >> operation;

		if (except == false) {
		
			if (operation == 1) {
				ABC.TriangleType();
			}

			else if (operation == 2) {
				cout << endl << "Лицето на триъгълника е: " << ABC.TriangleArea() << endl;
			}

			else if (operation == 3) {
				cout << endl << "Периметъра на триъгълника е: " << ABC.TrianglePerimeter() << endl;
			}

			else if (operation == 4) {
				cout << endl << "Медицентъра на триъгълника е: " << ABC.Medicenter() << endl;
			}

			else if (operation == 5) {

				input_data >> x >> y >> z;
				Point P(x, y, z);

				if ((P < ABC) == true) {
					cout << endl << "Точката лежи в триъгълника. " << endl;
				}
				else cout << endl << "Точката не лежи в триъгълника. " << endl;
			}

			else if (operation == 6) {

				input_data >> x >> y >> z;
				Point P(x, y, z);

				if ((P > ABC) == true) {
					cout << endl << "Точката лежи извън триъгълника. " << endl;
				}
				else cout << endl << "Точката не лежи в триъгълника. " << endl;
			}

			else if (operation == 7) {

				input_data >> x >> y >> z;
				Point P(x, y, z);

				if ((P == ABC) == true) {
					cout << endl << "Точката лежи на една от страните в триъгълника. " << endl;
				}
				else cout << endl << "Точката не лежи на някоя от страните в триъгълника. " << endl;
			}
		}
	}

	else if (object == 6) {
		double x, y, z;

		input_data >> x >> y >> z;
		Point A(x, y, z);
		input_data >> x >> y >> z;
		Point B(x, y, z);
		input_data >> x >> y >> z;
		Point C(x, y, z);
		input_data >> x >> y >> z;
		Point D(x, y, z);

		Tetrahedron ABCD(A, B, C, D);

		bool except = false;

		input_data >> operation;

		try {
			if (ABCD.GetA() == ABCD.GetB()) {
				throw EqualPointException(ABCD.GetA(), ABCD.GetB());
			}
			else if (ABCD.GetB() == ABCD.GetC()) {
				throw EqualPointException(ABCD.GetB(), ABCD.GetC());
			}
			else if (ABCD.GetA() == ABCD.GetC()) {
				throw EqualPointException(ABCD.GetA(), ABCD.GetC());
			}
			else if (ABCD.GetA() == ABCD.GetD()) {
				throw EqualPointException(ABCD.GetA(), ABCD.GetD());
			}
			else if (ABCD.GetD() == ABCD.GetC()) {
				throw EqualPointException(ABCD.GetD(), ABCD.GetC());
			}
			else if (ABCD.GetB() == ABCD.GetD()) {
				throw EqualPointException(ABCD.GetB(), ABCD.GetD());
			}
		}
		catch (EqualPointException & e) {
			cout << endl;
			e.getMessage();
			cout << endl;
			except = true;
		}

		if (except == false) {

			if (operation == 1) {

				if (ABCD.RegularTetrahedron() == true) {
					cout << endl << "Тетраедърът е правилен." << endl;
				}
				else cout << endl << "Тетраедърът не е правилен." << endl;

			}

			else if (operation == 2) {
				if (ABCD.OrthocentricTetrahedron() == true) {
					cout << endl << "Тетраедърът е ортогонален." << endl;
				}
				else cout << endl << "Тетраедърът не е ортогонален." << endl;
			}
			
			else if (operation == 3) {
				cout << endl << "Околната повърхнина на тетраедъра е " << ABCD.TetrahedronArea() << endl;
			}

			else if (operation == 4) {
				cout << endl << "Обемът на тетраедъра е " << ABCD.TetrahedronVolume() << endl;
			}

			else if (operation == 5) {

				input_data >> x >> y >> z;
				Point A(x, y, z);

				if (A < ABCD) {
					cout << endl << "Точката лежи в тетраедъра." << endl;
				}
				else cout << endl << "Точката не лежи в тетраедъра." << endl;
			}

			else if (operation == 6) {

				input_data >> x >> y >> z;
				Point A(x, y, z);

				if (A > ABCD) {
					cout << endl << "Точката лежи в извън тетраедъра." << endl;
				}
				else cout << endl << "Точката не лежи извън тетраедъра." << endl;
			}

			else if (operation == 7) {

				input_data >> x >> y >> z;
				Point A(x, y, z);

				if (A == ABCD) {
					cout << endl << "Точката лежи на една от страните на тетраедъра." << endl;
				}
				else cout << endl << "Точката не лежи на никоя от страните на тетраедъра." << endl;

			}
		}
	}


	input_data.close();
}
