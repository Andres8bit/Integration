package integration

type NumVal float64
type Func func(x NumVal) NumVal

// Trapizoid:
// Uses the Trapizoidal method to approximate the integral
// of any polynomial funstion f, where the values of f[i]
// represent the coefs. of the polynomial.
func Trapizoid(f []NumVal, a, b NumVal, n int) NumVal {
	h := (b - a) / NumVal(n)
	size := len(f) - 1
	fa := nested_Mutliplication(f[:], a, size)
	fb := nested_Mutliplication(f[:], b, size)
	sum := 0.5 * (fa + fb)

	for i := 1; i < n; i++ {
		x := a + NumVal(i)*h
		fx := nested_Mutliplication(f[:], x, size)
		sum = sum + fx
	}
	sum = sum * h

	return sum

}

// Trapizoid:
// Uses the Trapizoidal method to approximate the integral
// of any single variable function f. That is conitnous on [a,b] & 1:1
func Trapizoid2(f func(x NumVal) NumVal, a, b NumVal, n int) NumVal {
	h := (b - a) / NumVal(n)
	fa := f(a)
	fb := f(b)
	sum := 0.5 * (fa + fb)

	for i := 1; i < n; i++ {
		x := a + NumVal(i)*h
		fx := f(x)
		sum = sum + fx

	}
	sum = sum * h

	return sum
}

// ======== Helper functions ===========
//returns coef(x) using nested multiplication in O(n) time.
// used for the first trapiziod method which takes an array of ceof.
// representing a polynomial equation of degree len(coef)
func nested_Mutliplication(coef []NumVal, x NumVal, N int) NumVal {
	p := coef[N]
	for i := N - 1; i >= 0; i-- {
		p = coef[i] + x*p
	}
	return p
}
