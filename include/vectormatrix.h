#pragma once;

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

/****************************************************************************
** Mathematic Calculation Method
****************************************************************************/
namespace EITS{

#ifndef M_PI
	#define M_PI 3.141592653589793238
#endif

	template<class T> static inline T Rad2Deg(T rad) { return ( (rad)*(180.0/M_PI) ); }
	template<class T> static inline T Deg2Rad(T deg) { return ( (deg)*(M_PI/180.0) ); }

	template<class T>T lu(int n, T *a, int *ip)
	{
		int i, j, k, ii, ik;
		T t, u, det;
		T *weight = new T[n];

		det = 0;
		for (k = 0; k < n; k++) {
			ip[k] = k;
			u = 0;  
			for (j = 0; j < n; j++) {
				t = fabs(a[n*k+j]);  if (t > u) u = t;
			}
			if (u == 0) goto ENDPOINT;
			weight[k] = 1 / u;
		}
		det = 1;
		for (k = 0; k < n; k++) {
			u = -1;
			for (i = k; i < n; i++) {
				ii = ip[i]; 
				t = fabs(a[n*ii+k]) * weight[ii];
				if (t > u) {  u = t;  j = i;  }
			}
			ik = ip[j];
			if (j != k) {
				ip[j] = ip[k];  ip[k] = ik;
				det = -det;
			}
			u = a[n*ik+k];  det *= u;
			if (u == 0) goto ENDPOINT;
			for (i = k + 1; i < n; i++) {
				ii = ip[i];
				t = (a[n*ii+k] /= u);
				for (j = k + 1; j < n; j++)
					a[n*ii+j] -= t * a[n*ik+j];
			}
		}
	ENDPOINT:
		delete []weight;
		return det;
	}
	template<class T>T inverseLU(T *a_inv, T *a,int n)
	{
		int i, j, k, ii;
		T t, det;
		int *ip = new int[n];
		det = lu(n, a, ip);
		if (det != 0)
			for (k = 0; k < n; k++) {
				for (i = 0; i < n; i++) {
					ii = ip[i];  t = (ii == k);
					for (j = 0; j < i; j++)
						t -= a[n*ii+j] * a_inv[n*j+k];
					a_inv[n*i+k] = t;
				}
				for (i = n - 1; i >= 0; i--) {
					t = a_inv[n*i+k];  ii = ip[i];
					for (j = i + 1; j < n; j++)
						t -= a[n*ii+j] * a_inv[n*j+k];
					a_inv[n*i+k] = t / a[n*ii+i];
				}
			}
		delete []ip;
		return det;
	}
}

namespace EITS{
/****************************************************************************
** 2D Vector Class
****************************************************************************/
	template <class T>
	class Vector2{

	public:
		union { struct { T x, y; }; T X[2];};

		Vector2():x(0),y(0){}
		Vector2(T _x, T _y):x(_x),y(_y){}
		Vector2(T _X[2]):x(_X[0]),y(_X[1]){}

		Vector2 &operator=(const Vector2 &v){this->x=v.x;this->y=v.y;return(*this);}
		Vector2 &operator=(const T v){this->x=v;this->y=v;return(*this);}
		Vector2 &operator+=(const Vector2 &u){x+=u.x;y+=u.y;return(*this);}
		Vector2 &operator-=(const Vector2 &u){x-=u.x;y-=u.y;return(*this);}
		Vector2 &operator*=(const Vector2 &u){x*=u.x;y*=u.y;return(*this);}
		Vector2 &operator*=(const T a){x*=a;y*=a;return(*this);}
		Vector2 &operator/=(const T a){x/=a;y/=a;return(*this);}
		bool operator==(const Vector2 &u){return((x==u.x)&&(y==u.y));}
		bool operator!=(const Vector2 &u){return((x!=u.x)||(y!=u.y));}

		Vector2 operator +(){return *this;}
		Vector2 operator -(){return Vector2(-x,-y);}
		T &operator[](const int _elem){return X[_elem];}
		T &operator()(const int _elem){return X[_elem];}
		T abs(){T result=0;for(int i=0;i<2;i++)result+=X[i]*X[i];return sqrt((double)result);}

		friend Vector2 operator *(const Vector2<T> &u, const T a){Vector2<T> result(u);result*=a;return result;}
		friend Vector2 operator *(const T a, const Vector2<T> &u){Vector2<T> result(u);result*=a;return result;}
		friend Vector2 operator +(const Vector2<T> &u, const Vector2<T> &v){Vector2 result(u);result+=v;return result;}
		friend Vector2 operator -(const Vector2<T> &u, const Vector2<T> &v){Vector2 result(u);result-=v;return result;}
		friend Vector2 operator /(const Vector2<T> &u, const T a){Vector2 result(u);result/=a;return result;}
	};
	template <class T>T operator *(const Vector2<T> &u, const Vector2<T> &v){T result=u.x*v.x+u.y*v.y; return result;}
	template <class T>std::ostream& operator<<(std::ostream &stream, Vector2<T> &v){return stream <<v.x<<","<<v.y;}
}

/****************************************************************************
** 3D double Vector Class
****************************************************************************/
namespace EITS{
	template <class T>
	class Vector3{
	public:
		union { struct { T x, y, z; }; T X[3];};
		Vector3():x(0),y(0),z(0){}
		Vector3(T _x, T _y, T _z):x(_x),y(_y),z(_z){}
		Vector3(T _X[3]):x(_X[0]),y(_X[1]),z(_X[2]){}
		Vector3 &operator=(T v){x=v;y=v;z=v;return(*this);}
		Vector3 &operator=(const Vector3 &v){x=v.x;y=v.y;z=v.z;return(*this);}
		Vector3 &operator+=(const Vector3 &u){x+=u.x;y+=u.y;z+=u.z;return(*this);}
		Vector3 &operator-=(const Vector3 &u){x-=u.x;y-=u.y;z-=u.z;return(*this);}
		Vector3 &operator*=(const Vector3 &u){x*=u.x;y*=u.y;z*=u.z;return(*this);}
		Vector3 &operator*=(T a){x*=a;y*=a;z*=a;return(*this);}
		Vector3 &operator/=(T a){x/=a;y/=a;z/=a;return(*this);}
		bool operator==(const Vector3 &u){return((x==u.x)&&(y==u.y)&&(z==u.z));}
		bool operator!=(const Vector3 &u){return((x!=u.x)||(y!=u.y)||(z!=u.z));}
		bool operator!=(const T &u){return((x!=u)||(y!=u)||(z!=u));}
		bool operator<(const Vector3 &u){return((x<u.x)&&(y<u.y)&&(z<u.z));}
		bool operator>(const Vector3 &u){return((x>u.x)&&(y>u.y)&&(z>u.z));}
		bool operator<(const T &u){return(x<u)&&(y<u)&&(z<u);}
		bool operator>(const T &u){return(x>u)&&(y>u)&&(z>u);}
		bool operator<=(const T &u){return(x<=u)&&(y<=u)&&(z<=u);}
		bool operator>=(const T &u){return(x>=u)&&(y>=u)&&(z>=u);}

		Vector3 operator +(){return *this;}
		Vector3 operator -(){return Vector3(-x,-y,-z);}
		T &operator[](int _elem){return X[_elem];}
		T &operator()(int _elem){return X[_elem];}
		T abs(){T result=0;for(int i=0;i<3;i++)result+=X[i]*X[i];return sqrt((double)result);}

		friend Vector3 operator *(const Vector3<T> &u, const T a){Vector3<T> result(u);result*=a;return result;}
		friend Vector3 operator *(const T a, const Vector3<T> &u){Vector3<T> result(u);result*=a;return result;}
		friend Vector3 operator +(const Vector3<T> &u, const Vector3<T> &v){Vector3 result(u);result+=v;return result;}
		friend Vector3 operator -(const Vector3<T> &u, const Vector3<T> &v){Vector3 result(u);result-=v;return result;}
		friend Vector3 operator /(const Vector3<T> &u, const T a){Vector3 result(u);result/=a;return result;}
		friend Vector3 operator %(const Vector3<T> &u, const Vector3<T> &v){return Vector3(u.y*v.z-u.z*v.y,u.z*v.x-u.x*v.z,u.x*v.y-u.y*v.x);}
	};
	template <class T>T operator *(const Vector3<T> &u, const Vector3<T> &v){T result; result=u.x*v.x+u.y*v.y+u.z*v.z; return result;}
	template <class T> std::ostream &operator<<(std::ostream &stream, const Vector3<T> &v){return stream <<v.x<<","<<v.y<<","<<v.z;}
}


/****************************************************************************
** 4D double Vector Class
****************************************************************************/
namespace EITS{
	template <class T>
	class Vector4{
	public:
		union { struct { T x, y, z, w; }; T X[4];};

		Vector4():x(0),y(0),z(0),w(0){}
		Vector4(T _x, T _y, T _z, T _w):x(_x),y(_y),z(_z),w(_w){}
		Vector4(T _X[4]):x(_X[0]),y(_X[1]),z(_X[2]),w(_X[3]){}

		Vector4 &operator=(const Vector4 &v){this->x=v.x;this->y=v.y;this->z=v.z;this->w=v.w;return(*this);}
		Vector4 &operator=(const T v){this->x=v;this->y=v;this->z=v;this->w=v;return(*this);}
		Vector4 &operator+=(const Vector4 &u){x+=u.x;y+=u.y;z+=u.z;w+=u.w;return(*this);}
		Vector4 &operator-=(const Vector4 &u){x-=u.x;y-=u.y;z-=u.z;w-=u.w;return(*this);}
		Vector4 &operator*=(const Vector4 &u){x*=u.x;y*=u.y;z*=u.z;w*=u.w;return(*this);}
		Vector4 &operator*=(T a){x*=a;y*=a;z*=a;w*=a;return(*this);}
		Vector4 &operator/=(T a){x/=a;y/=a;z/=a;w/=a;return(*this);}
		const bool operator==(const Vector4 &u){return((x==u.x)&&(y==u.y)&&(z==u.z)&&(w==u.w));}
		const bool operator!=(const Vector4 &u){return((x!=u.x)||(y!=u.y)||(z!=u.z)||(w!=u.w));}

		Vector4 operator +(){return *this;}
		Vector4 operator -(){return Vector4(-x,-y,-z,-w);}
		T &operator[](int _elem){return X[_elem];}
		T &operator()(int _elem){return X[_elem];}
		T abs(){double result=0;for(int i=0;i<4;i++)result+=X[i]*X[i];return sqrt((double)result);}

		friend Vector4 operator *(Vector4<T> &u, T a){Vector4 result(u);result*=a;return result;}
		friend Vector4 operator *(T a, Vector4<T> &u){Vector4 result(u);result*=a;return result;}
		friend Vector4 operator +(Vector4<T> &u, Vector4<T> &v){Vector4 result(u);result+=v;return result;}
		friend Vector4 operator -(Vector4<T> &u, Vector4<T> &v){Vector4 result(u);result-=v;return result;}
		friend Vector4 operator /(Vector4<T> &u, T a){Vector4 result(u);result/=a;return result;}
	};
	template <class T>T operator *(Vector4<T> &u, Vector4<T> &v){T result; result=u.x*v.x+u.y*v.y+u.z*v.z; return result;}
	template <class T>std::ostream& operator<<(std::ostream& stream, Vector4<T> &v){return stream <<v.x<<","<<v.y<<","<<v.z<<","<<v.w;}
}

/****************************************************************************
** ND double Vector Class
****************************************************************************/
namespace EITS{
	template <class T>
	class VectorN{
	public:
		int n;
		T *X;
		VectorN():n(0),X(0){}
		VectorN(int _n):n(_n){if(_n>1)X=new T[n];for(int i=0;i<n;i++)X[i]=0;}
		VectorN(const VectorN &u){n=u.n;X=new T[n];for(int i=0;i<n;i++)X[i]=u.X[i];}
		~VectorN(){if(n>1)delete []X;n=0;}

		VectorN &operator=(const VectorN &u){
			if(&u==this)return *this;
			if(n==u.n){for(int i=0;i<n;i++)X[i]=u.X[i];return *this;}
			else {if(n>1)delete []X;n=u.n;X=new T[n];
			for(int i=0;i<n;i++)X[i]=u.X[i];}
			return *this;
		}
		VectorN &operator=(const Vector3<T> &u){
			if(n>=3){for(int i=0;i<3;i++)X[i]=u.X[i];return *this;}
			else return *this;
		}
		VectorN &operator=(const T u){
			for(int i=0;i<n;i++)X[i]=u;return *this;
		}
		VectorN &operator+=(const VectorN &u){if(n==u.n){for(int i=0;i<n;i++)X[i]+=u.X[i];}return(*this);}
		VectorN &operator-=(const VectorN &u){if(n==u.n){for(int i=0;i<n;i++)X[i]-=u.X[i];}return(*this);}
		VectorN &operator*=(T a){for(int i=0;i<n;i++)X[i]*=a;return(*this);}
		VectorN &operator/=(T a){for(int i=0;i<n;i++)X[i]/=a;return(*this);}
		bool operator==(const VectorN &u){if(n!=u.n)return false;for(int i=0;i<n;i++){if(X[i]!=u.X[i])return false;}return true;}
		bool operator!=(const VectorN &u){if(n!=u.n)return true;for(int i=0;i<n;i++){if(X[i]!=u.X[i])return true;}return false;}
		VectorN &operator +(){return(*this);}
		VectorN &operator -(){for(int i=0;i<n;i++)X[i]=-X[i];return(*this);}
		T &operator[](int _elem){return X[_elem];}
		T &operator()(int _elem){return X[_elem];}
		T abs(){
			T result=0;
			for(int i=0;i<n;i++)result+=X[i]*X[i];
			return sqrt((double)result);
		}
		void malloc(int _n){
			n=_n; X=new T[n];for(int i=0;i<n;i++)X[i]=0;
		}
		void free(){
			if(n>1)delete []X;
			this->n=0;
		}

		friend VectorN operator *(VectorN<T> &u, T a){VectorN result(u);result*=a;return result;}
		friend VectorN operator *(T a, VectorN<T> &u){VectorN<T> result(u);result*=a;return result;}
		friend VectorN operator +(VectorN<T> &u, VectorN<T> &v){VectorN result(u);result+=v;return result;}
		friend VectorN operator -(VectorN<T> &u, VectorN<T> &v){VectorN result(u);result-=v;return result;}
		friend VectorN operator /(VectorN<T> &u, T a){VectorN result(u);result/=a;return result;}
		friend VectorN operator %(VectorN<T> &u, VectorN<T> &v){
			VectorN result(3);
			if(u.n!=3)return u;
			result.X[0]=u.X[1]*v.X[2]-u.X[2]*v.X[1];
			result.X[1]=u.X[2]*v.X[0]-u.X[0]*v.X[2];
			result.X[2]=u.X[0]*v.X[1]-u.X[1]*v.X[0];
			return result;
		}
	};
	template <class T>T operator *(VectorN<T> &u, VectorN<T> &v){T result=0; for(int i=0;i<u.n;i++)result+=u.X[i]*v.X[i]; return result;}
	template <class T>std::ostream& operator<<(std::ostream& stream, VectorN<T> &v){for(int i=0;i<v.n-1;i++)stream<<v.X[i]<<",";stream<<v.X[v.n-1]; return stream;}
}

namespace EITS{
	typedef	Vector2<int> Vector2i;
	typedef	Vector2<double> Vector2d;
	typedef	Vector2<float> Vector2f;

	typedef	Vector3<int> Vector3i;
	typedef	Vector3<double> Vector3d;
	typedef	Vector3<float> Vector3f;

	typedef	Vector4<int> Vector4i;
	typedef	Vector4<double> Vector4d;
	typedef	Vector4<float> Vector4f;

	typedef	VectorN<int> VectorNi;
	typedef	VectorN<double> VectorNd;
	typedef	VectorN<float> VectorNf;
}

/****************************************************************************
** N*M dobule Matrix Class
****************************************************************************/
namespace EITS{
	template<class T>
	class Matrix{
	public:
		int n;
		int m;
		T *X;

		Matrix():n(0),m(0),X(0){}
		Matrix(int _n,int _m):n(_n),m(_m){if(_n*_m>=1)X=new T[n*m+1];memset(X, 0, sizeof(T)*n*m);}
		Matrix(const Matrix &u){n=u.n;m=u.m;X=new T[n*m+1];memcpy(X, u.X, sizeof(T)*n*m);}
		~Matrix(){if(n*m>=1)delete []X;n=m=0;}
		Matrix &operator=(const Matrix &u){
			if(&u==this)return *this;
			if(n==u.n&&m==u.m){memcpy(X, u.X, sizeof(T)*n*m);return *this;}
			else {if(n*m>=1)delete []X;n=u.n;m=u.m;X=new T[n*m+1];
			memcpy(X, u.X, sizeof(T)*n*m);}return *this;
		}
		Matrix &operator+=(const Matrix &u){if(n==u.n&&m==u.m){for(int i=0;i<n*m;i++)X[i]+=u.X[i];}return(*this);}
		Matrix &operator-=(const Matrix &u){if(n==u.n&&m==u.m){for(int i=0;i<n*m;i++)X[i]-=u.X[i];}return(*this);}
		Matrix &operator*=(const T a){for(int i=0;i<n*m;i++)X[i]*=a;return(*this);}
		Matrix &operator/=(const T a){for(int i=0;i<n*m;i++)X[i]/=a;return(*this);}
		const bool operator==(const Matrix &u){if(m!=u.m||n!=u.n)return false;for(int i=0;i<n*m;i++){if(X[i]!=u.X[i])return false;}return true;}
		const bool operator!=(const Matrix &u){if(m!=u.m||n!=u.n)return true;for(int i=0;i<n*m;i++){if(X[i]!=u.X[i])return true;}return false;}
		Matrix &operator +(){return(*this);}
		Matrix &operator -(){for(int i=0;i<n*m;i++)X[i]=-X[i];return(*this);}
		Matrix<T> inv(void){
			Matrix<T> src=*this;Matrix<T> dst(n,m);
			if(src.n==1&&src.m==1)
				dst.X[0]=1.0/src.X[0];
			else inverseLU(dst.X, src.X,n);
/*
			for(int i=0;i<n;i++){
				for(int j=0;j<n;j++){
					dst.X[n*j+i]=(i==j)?1.0:0.0;
				}
			}
			//‘|‚«o‚µ–@
			double buf;
			for(int i=0;i<n;i++){
				buf=1.0/src.X[n*i+i];
				for(int j=0;j<n;j++){
					src.X[n*j+i]*=buf;
					dst.X[n*j+i]*=buf;
				}
				for(int j=0;j<n;j++){
					if(i!=j){
						buf=src.X[n*i+j];
						for(int k=0;k<n;k++){
							src.X[n*k+j]-=src.X[n*k+i]*buf;
							dst.X[n*k+j]-=dst.X[n*k+i]*buf;
						}
					}
				}
			}
*/
			return dst;
		}
		T det(void){
			int *ip;Matrix src(*this);
			T result;ip = new int[n];
			result=lu(src.n,src.X,ip);
			delete []ip;
			return result;
		}
		Matrix trn(void){
			Matrix result(m,n);
			for(int ii=0;ii<m;ii++)
				for(int jj=0;jj<n;jj++)
					result.X[m*jj+ii]=X[n*ii+jj];
			return result;
		}

		T &operator()(int _row, int _col){return X[m*_row+_col];}

		void free(){if(this->n*this->m)delete []X;this->n=this->m=0;}
		void malloc(int _n ,int _m){n=_n;m=_m;X=new T[n*m];for(int i=0;i<n*m;i++)X[i]=0;}
		void identity(){
			int i;
			memset(this->X, 0, sizeof(T) * this->n*this->m);
			for(i=0;i<this->m;i++)
				if(i<this->n)this->X[m*i+i]=1;
		}
		friend Matrix operator *(const Matrix<T> &u, const T a){Matrix result(u);result*=a;return result;}
		friend Matrix operator *(const T a, const Matrix<T> &u){Matrix result(u);result*=a;return result;}
		friend Matrix operator *(const Matrix<T> &u, const Matrix<T> &v){
			Matrix result(v.n,u.m);
			for(int i=0;i<u.m;i++)
				for(int j=0;j<v.n;j++)
					for(int k=0;k<u.n;k++)
						result.X[v.n*i+j]+=u.X[u.n*i+k]*v.X[v.n*k+j];
			return result;
		}
		friend Matrix operator +(const Matrix<T> &u, const Matrix<T> &v){Matrix result(u);result+=v;return result;}
		friend Matrix operator -(const Matrix<T> &u, const Matrix<T> &v){Matrix result(u);result-=v;return result;}
		friend Matrix operator /(const Matrix<T> &u, const T a){Matrix result=u;result/=a;return result;}

	};
	template<class T>VectorN<T> operator *(const Matrix<T> &u, const VectorN<T> &v){
		VectorN<T> result(u.m);
		if(u.n!=v.n)return v;
		for(int i=0;i<u.m;i++){
			for(int j=0;j<v.n;j++)
				result.X[i]+=u.X[v.n*i+j]*v.X[j];
		}
		return result;
	}
	template<class T> Vector3<T> operator *(const Matrix<T> &u, const Vector3<T> &v){
		Vector3<T> result;
		if(u.n!=3)return v;
		for(int i=0;i<3;i++){
			result.X[i]=0;
			for(int j=0;j<3;j++)
				result.X[i]+=u.X[3*i+j]*v.X[j];
		}
		return result;
	}
	template <class T>std::ostream& operator<<(std::ostream& stream, Matrix<T> &v){
		for(int j=0;j<v.m;j++){
			for(int i=0;i<v.n;i++){
				stream<<v.X[v.n*j+i];
				if(i < v.n - 1)
					stream<<" ";
				else
					stream<<std::endl;
			}
		}
		return stream;
	}

	typedef	Matrix<int> Matrixi;
	typedef	Matrix<double> Matrixd;
	typedef	Matrix<float> Matrixf;
}