#pragma once

#include <LC/math/linalg/Vector/Vector.hpp>

namespace lc
{
	// ================================================================
	//                          QUATERNION
	// ================================================================

	template <lc::math::concepts::ScalarType T = double>
	class Quaternion
	{
	public:

		// ================================
		//          CONSTRUCTORS
		// ================================

		// default constructors: sets the components to 0
		constexpr Quaternion() = default;

		// default copy constructor
		constexpr Quaternion(const Quaternion<T>&) = default;

		// default move constructor
		constexpr Quaternion(Quaternion<T>&&) = default;

		// construct a quaternion from components
		constexpr Quaternion(const T& a, const T& b = T(), const T& c = T(), const T& d = T());

		// copy a quaternion of different type
		// TODO check convertible
		template <lc::math::concepts::ScalarType T2>
		constexpr Quaternion(const Quaternion<T2>& other);

		// construct a quaternion from a 3d vector and a real part
		// TODO check convertible
		template <lc::math::concepts::ScalarType T2>
		constexpr Quaternion(const Vector3<T2>& v, const T& real_part = T(1));

		// construct a quaternion from a 4d vector where the fourth component is the real part:
		// a = v[3]
		// b = v[0];
		// c = v[1];
		// d = v[2];
		template <lc::math::concepts::ScalarType T2>
		constexpr Quaternion(const Vector4<T2>& v);

		// ================================
		//             ACCESS
		// ================================

		// convert a quaternion to a 3d vector ignoring the real part
		template <class T2 = T>
		constexpr Vector3<T2> get_vec3(void) const;

		// convert a quaternion to a 4d vector there the fourth component is the real part
		template <class T2 = T>
		constexpr Vector4<T2> get_vec4(void) const;

		// ================================
		//           OPEARTIONS
		// ================================

		// get the norm of the quaternion
		constexpr T norm(void) const;

		// get the squared norm of the quaternion
		constexpr T norm2(void) const;

		// normalize the quaternion
		constexpr Quaternion<T>& normalize(void);

		// get a normilized copy of the quaternion
		constexpr Quaternion<T> normalized(void) const;

		// invert the quaternion
		constexpr Quaternion<T>& invert(void);

		// returns the inverse of the quaternion
		constexpr Quaternion<T> inverse(void) const;

		// conjugate the quaternion
		constexpr Quaternion<T>& conjugate(void);

		// reurns the conjugated of the quaternion
		constexpr Quaternion<T> conj(void) const;

		// set this to opposite (this = -this)
		constexpr Quaternion<T>& flip(void);

		// get the opposite
		constexpr Quaternion<T> opposite(void) const;

		// ================================
		//           OPERATORS
		// ================================

		// copy assignment
		constexpr Quaternion<T>& operator=(const Quaternion<T>&) = default;

		// move assignment
		constexpr Quaternion<T>& operator=(Quaternion<T>&&) = default;

		// get the opposite
		constexpr Quaternion<T> operator-(void) const;

		// quaternion addition
		constexpr Quaternion<T>& operator+=(const T& right);

		// quaternion addition
		constexpr Quaternion<T>& operator+=(const Quaternion<T>& right);

		// quaternion addition
		constexpr Quaternion<T> operator+(const T& right) const;

		// quaternion addition
		constexpr Quaternion<T> operator+(const Quaternion<T>& right) const;

		// quaternion subtraction
		constexpr Quaternion<T>& operator-=(const T& right);

		// quaternion subtraction
		constexpr Quaternion<T>& operator-=(const Quaternion<T>& right);

		// quaternion subtraction
		constexpr Quaternion<T> operator-(const T& right) const;

		// quaternion subtraction
		constexpr Quaternion<T> operator-(const Quaternion<T>& right) const;

		// scale the quaternion
		constexpr Quaternion<T>& operator*=(const T& right);

		// quaternion multipication
		constexpr Quaternion<T>& operator*=(const Quaternion<T>& right);

		// scale the quaternion
		constexpr Quaternion<T> operator*(const T& right) const;

		// quaternion multipication
		constexpr Quaternion<T> operator*(const Quaternion<T>& right) const;

		// scale the quaternion
		constexpr Quaternion<T>& operator/=(const T& right);

		// quaternion division
		constexpr Quaternion<T>& operator/=(const Quaternion<T>& right);

		// scale the quaternion
		constexpr Quaternion<T> operator/(const T& right) const;

		// quaternion division
		constexpr Quaternion<T> operator/(const Quaternion<T>& right) const;



	public:

		// ================================
		//             DATA
		// ================================

		// the real component (a + bi + cj + dk)
		T a = T();

		// the "i" component (a + bi + cj + dk)
		T b = T();

		// the "j" component (a + bi + cj + dk)
		T c = T();
		
		// the "k" component (a + bi + cj + dk)
		T d = T();
	};

	// ================================================================
	//                       external functions
	// ================================================================

	// ================================
	//             BASIC
	// ================================

	// norm of a quaternion
	template <class T>
	T abs(const lc::Quaternion<T>& q);

	// ================================
	//          EXPONENTIAL
	// ================================

	// quaterion base 'e' exponential
	template <class T>
	constexpr lc::Quaternion<T> exp(const lc::Quaternion<T>& q);

	// quaterion base 'e' logarithm
	template <class T>
	constexpr lc::Quaternion<T> ln(lc::Quaternion<T> q);

	// quaterion base 'e' logarithm
	template <class T>
	constexpr lc::Quaternion<T> log(const lc::Quaternion<T>& q);

	// ================================
	//          TRIGONOMETRIC
	// ================================

	// TODO

	// ================================
	//           HYPERBOLIC
	// ================================

	// TODO

	// ================================
	//              I/O
	// ================================

	// formatted output
	template <class T>
	std::ostream& operator<<(std::ostream& ostream, const lc::Quaternion<T>& q);
}






// ================================================================================================================================
// ================================================================================================================================
//                                                           INL
// ================================================================================================================================
// ================================================================================================================================

namespace lc
{
	// ================================================================
	//                          QUATERNION
	// ================================================================

	// ================================
	//          CONSTRUCTORS
	// ================================

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T>::Quaternion(const T& a, const T& b, const T& c, const T& d) :
		a(a),
		b(b),
		c(c),
		d(d)
	{
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	template <lc::math::concepts::ScalarType T2>
	inline constexpr Quaternion<T>::Quaternion(const Quaternion<T2>& other) :
		a(other.a),
		b(other.b),
		c(other.c),
		d(other.d)
	{
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	template <lc::math::concepts::ScalarType T2>
	inline constexpr Quaternion<T>::Quaternion(const Vector3<T2>& v, const T& real_part) :
		a(real_part),
		b(v[0]),
		c(v[1]),
		d(v[2])
	{
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	template <lc::math::concepts::ScalarType T2>
	inline constexpr Quaternion<T>::Quaternion(const Vector4<T2>& v) :
		a(v[3]),
		b(v[0]),
		c(v[1]),
		d(v[2])
	{
	}

	// ================================
	//             ACCESS
	// ================================

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	template <class T2>
	inline constexpr Vector3<T2> Quaternion<T>::get_vec3(void) const
	{
		return { b, c, d };
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	template <class T2>
	inline constexpr Vector4<T2> Quaternion<T>::get_vec4(void) const
	{
		return { b, c, d, a };
	}

	// ================================
	//           OPEARTIONS
	// ================================

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr T Quaternion<T>::norm(void) const
	{
		return sqrt(this->norm2());
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr T Quaternion<T>::norm2(void) const
	{
		T sum_v = 0;
		sum_v += sqr(abs(a));
		sum_v += sqr(abs(b));
		sum_v += sqr(abs(c));
		sum_v += sqr(abs(d));
		return sum_v;
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T>& Quaternion<T>::normalize(void)
	{
		auto norm = this->norm();
		if (norm == T(0))
			;// nothing to do// *this = T(0);
		else
			*this /= norm;
		return *this;
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T> Quaternion<T>::normalized(void) const
	{
		auto tmp = *this;
		tmp.normalize();
		return tmp;
		return *this;
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T>& Quaternion<T>::invert(void)
	{
		*this = this->inverse();
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T> Quaternion<T>::inverse(void) const
	{
		auto conjugate = this->conj();
		return conjugate / this->norm2();
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T>& Quaternion<T>::conjugate(void)
	{
		b = -b;
		c = -c;
		d = -d;
		return *this;
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T> Quaternion<T>::conj(void) const
	{
		auto tmp = *this;
		tmp.conjugate();
		return tmp;
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T>& Quaternion<T>::flip(void)
	{
		a = -a;
		b = -b;
		c = -c;
		d = -d;
		return *this;
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T> Quaternion<T>::opposite(void) const
	{
		auto tmp = *this;
		tmp.flip();
		return tmp;
	}

	// ================================
	//           OPERATORS
	// ================================

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T> Quaternion<T>::operator-(void) const
	{
		return this->opposite();
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T>& Quaternion<T>::operator+=(const T& right)
	{
		a += right;
		b += right;
		c += right;
		d += right;
		return *this;
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T>& Quaternion<T>::operator+=(const Quaternion<T>& right)
	{
		a += right.a;
		b += right.b;
		c += right.c;
		d += right.d;
		return *this;
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T> Quaternion<T>::operator+(const T& right) const
	{
		auto tmp = *this;
		tmp += right;
		return tmp;
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T> Quaternion<T>::operator+(const Quaternion<T>& right) const
	{
		auto tmp = *this;
		tmp += right;
		return tmp;
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T>& Quaternion<T>::operator-=(const T& right)
	{
		a -= right;
		b -= right;
		c -= right;
		d -= right;
		return *this;
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T>& Quaternion<T>::operator-=(const Quaternion<T>& right)
	{
		a -= right.a;
		b -= right.b;
		c -= right.c;
		d -= right.d;
		return *this;
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T> Quaternion<T>::operator-(const T& right) const
	{
		auto tmp = *this;
		tmp -= right;
		return tmp;
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T> Quaternion<T>::operator-(const Quaternion<T>& right) const
	{
		auto tmp = *this;
		tmp -= right;
		return tmp;
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T>& Quaternion<T>::operator*=(const T& right)
	{
		a *= right;
		b *= right;
		c *= right;
		d *= right;

		return *this;
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T>& Quaternion<T>::operator*=(const Quaternion<T>& right)
	{
		*this = *this + right;

		return *this;
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T> Quaternion<T>::operator*(const T& right) const
	{
		Quaternion<T> tmp = *this;
		tmp *= right;
		return tmp;
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T> Quaternion<T>::operator*(const Quaternion<T>& right) const
	{
		return {
			a * right.a - b * right.b - c * right.c - d * right.d,
			a * right.b + b * right.a + c * right.d - d * right.c,
			a * right.c - b * right.d + c * right.a + d * right.b,
			a * right.d + b * right.c - c * right.b + d * right.a
		};
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T>& Quaternion<T>::operator/=(const T& right)
	{
		a /= right;
		b /= right;
		c /= right;
		d /= right;

		return *this;
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T>& Quaternion<T>::operator/=(const Quaternion<T>& right)
	{
		*this *= right.inverse();
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T> Quaternion<T>::operator/(const T& right) const
	{
		auto tmp = *this;
		tmp /= right;
		return tmp;
	}

	////////////////////////////////////////////////////////////////
	template <lc::math::concepts::ScalarType T>
	inline constexpr Quaternion<T> Quaternion<T>::operator/(const Quaternion<T>& right) const
	{
		auto tmp = *this;
		tmp /= right;
		return tmp;
	}
}

// ================================================================
//                       external functions
// ================================================================

// ================================
//             BASIC
// ================================

////////////////////////////////////////////////////////////////
template <class T>
T abs(const lc::Quaternion<T>& q)
{
	return q.norm();
}

// ================================
//          EXPONENTIAL
// ================================

namespace lc::priv
{
	// TODO sposta
	// nota funziona solo per tipi reali e potrebbe non essere corretta per complessi o ipercomplessi!!!!
	template <class T>
	inline T sinx_over_x(T x, const T& approx_space)
	{
		T real_approx_space = abs(approx_space);
		T real_x = abs(x);
		if (real_x >= real_approx_space)
			return sin(x) / x;
		else
		{
			T t = real_x / real_approx_space;
			return t * sin(real_approx_space) / real_approx_space + (T(1) - t);
		}

		// TODO assolutamente: al momento � un interpolazione lineare da 0 a approx_space, ma in realt�
		// in questo modo se ci sono dei Nan o valori sballati ce li portiamo dietro... � meglio fare
		// sempre 1 in un primo tratto e poi l'interpolazone in un secondo:
		//       v- primo tratto costante, cosi' vanno via i Nan
		// 1->  ___
		//          \
		//           \  <- secondo tratto interpolato
		//            \
		//               ... poi da qui continua sin(x)/x esplicito ... 
		// inoltre 't' � calcolato male perch� potrebbe portarsi dietro informazioni sulla derivata che
		// invece vanno eliminate in quanto deve essere un semplice parametro REALE per l'interpolazione
	}
}

////////////////////////////////////////////////////////////////
template <class T>
inline constexpr lc::Quaternion<T> lc::exp(const lc::Quaternion<T>& q)
{
	using namespace lc;

	auto v = (Vector3<T>)q;
	auto norm_v = v.norm();
	return
		Quaternion<T>(exp(q.a)) *
		(
			Quaternion<T>(cos(norm_v)) +
			Quaternion<T>(v, 0) * priv::sinx_over_x<T>(norm_v, (T)1e-3)// !!!!!
		);
}

////////////////////////////////////////////////////////////////
template <class T>
inline constexpr lc::Quaternion<T> lc::ln(lc::Quaternion<T> q)
{
	// TODO da ricontrollare
	using namespace lc;

	auto norm_q = q.norm();

	q /= norm_q;
	auto v = (Vector3<T>)q;
	
	auto sin_teta = v.norm();

	v *= asin(sin_teta) / sin_teta;

	// TODO clamp i valori prima di arcoseno

	return Quaternion<T>(v, log(norm_q));

	// TODO pensa: i quaternioni non sono associativi, come fa il log ad essere associativo??????????
    // https://en.wikipedia.org/wiki/Quaternion#Exponential,_logarithm,_and_power_functions
}

////////////////////////////////////////////////////////////////
template <class T>
inline constexpr lc::Quaternion<T> lc::log(const lc::Quaternion<T>& q)
{
	return ln(q);
}

// ================================
//          TRIGONOMETRIC
// ================================

// ================================
//           HYPERBOLIC
// ================================

// ================================
//              I/O
// ================================

////////////////////////////////////////////////////////////////
template <class T>
inline std::ostream& lc::operator<<(std::ostream& ostream, const lc::Quaternion<T>& q)
{
	ostream << q.a << " + " << q.b << "i + " << q.c << "j + " << q.d << "k";

	return ostream;
}