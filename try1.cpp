#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
using namespace std;

class Complex{
    private:
        double real;
        double imag;
    public:
    //CONSTRUCTORS
    Complex(double re, double im) {real = re; imag = im;};
    Complex() {real = 0; imag = 0;}; 
    Complex(double x) {real = x; imag = 0;} ;
    /*Complex(Polar p) {
        real = p.giveMagn() * cos(p.giveAngle());
        imag = p.giveMagn() * sin(p.giveAngle());
    }*/    //questo non funziona non so perchè ancora
    //BASIC FUNCTIONS
    double angle(){
        double num;
        if(real == 0)
             {
                if(imag > 0)
                    {return M_PI_2;}
                return M_PI_2 + M_PI ;
             }
        num = atan(imag / real);
        if(real > 0)
            {return num;}
        return num + M_PI;
    }

    double magnitude(){
        return sqrt(real*real + imag*imag);
    }
    double giveReal() {return real;} ;
    double giveImag() {return imag;} ;
    Complex conjugate() {
        Complex p;
        p.real = this->real;
        p.imag = - (this->imag);
        return p;
    }
    //OPERATORS
    Complex operator +(const Complex& c2){
        Complex ris;
        ris.real = this->real + c2.real;
        ris.imag = this->imag + c2.imag;
        return ris;
    }
    Complex operator -(const Complex& c2){
        Complex ris;
        ris.real = this->real - c2.real;
        ris.imag = this->imag - c2.imag;
        return ris;
    }

    Complex operator *(const Complex& c2){
        Complex ris;
        ris.real = (this->real * c2.real) - (this->imag * c2.imag);
        ris.imag = (this->real * c2.imag) + (this->imag * c2.real);
        return ris;
    }
    Complex operator *(const int& x){
        Complex ris;
        ris.real = this->real * x;
        ris.imag = this->imag * x;
        return ris;
    }
    

    Complex operator /( Complex& c2){
        Complex ris;
        ris.real = ( ( (*this) * (c2.Complex::conjugate()) ).Complex::giveReal() ) / ( c2.Complex::magnitude() * c2.Complex::magnitude() ) ;
        ris.imag = ( ( (*this) * (c2.Complex::conjugate()) ).Complex::giveImag() ) / ( c2.Complex::magnitude() * c2.Complex::magnitude() ) ;
        return ris;
    }


};

class Polar{
    private:
        double angle;
        double magnitude;
    public:
    //Constructors
    Polar(double magn, double theta) {angle = theta; magnitude = magn;} ;
    Polar() {magnitude = 0; angle = 0;} ;
    Polar(double x) {magnitude = x; angle = 0;} ;
    Polar(Complex z) {magnitude = z.Complex::magnitude(); angle = z.Complex::angle() ;} ;
    //Basic functions
    double giveMagn() {return magnitude;} ;
    double giveAngle() {return angle;} ;
};

int main(){
    double re, im;
    cin >> re;
    cin >> im;
    Complex z(re, im);
    cout << "z ----- REAL: " << z.Complex::giveReal() << " IMAG: " << z.Complex::giveImag() << endl ;
    cin >> re;
    cin >> im;
    Complex w(re, im);
    cout << "w ----- REAL: " << w.Complex::giveReal() << " IMAG: " << w.Complex::giveImag() << endl ;
    
    Complex v;
    v = z+w;
    cout << "z+w --- REAL: " << v.Complex::giveReal() << " IMAG: " << v.Complex::giveImag() << endl ;
    v = z-w;
    cout << "z-w --- REAL: " << v.Complex::giveReal() << " IMAG: " << v.Complex::giveImag() << endl ;
    v = z*w;
    cout << "z*w --- REAL: " << v.Complex::giveReal() << " IMAG: " << v.Complex::giveImag() << endl ;
    double x;
    cin >> x;
    v = z*x;
    cout << "z*x --- REAL: " << v.Complex::giveReal() << " IMAG: " << v.Complex::giveImag() << endl ;
    v = z/w;
    cout << "z/w --- REAL: " << v.Complex::giveReal() << " IMAG: " << v.Complex::giveImag() << endl ;
    Polar zp(z);
    
    cout << "THETA: " << zp.Polar::giveAngle() << endl << "MAGNITUDE: " << zp.Polar::giveMagn() << endl;
    return 0;
}