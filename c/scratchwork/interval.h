#include <vector>
#include <concept>
#include <cmath>

template struct Interval <T>
{
	T p0, pn;
	int npts;
	T step;
	std::vector<T> data;
	
	Interval()
	{
		this->p0 = 0;
		this->pn = 1;
		this->npts = 11;
		this->step = (this->pn - this->p0) / (this->npts - 1);
		
		this->initialize();
	}
	
	Interval(T p0, T pn, int npts) {
		this->p0 = p0;
		this->pn = pn;
		this->npts = npts;
		this->step = (pn - p0) / (npts - 1);
		
		this->initialize();
	}
	
	Interval(T p0, T pn, T step) {
		this->p0 = p0;
		this->pn = pn;
		this->npts = (pn - p0) / step + 1;
		this->step = (pn - p0) / (npts - 1);
		
		this->initialize();
	}
	
	Interval(const std::vector<T> &other)
	{
		this->p0 = other.front();
		this->pn = other.back();
		this->npts = other.size;
		this->step = (this->pn - this->p0) / (this->npts - 1);
		
		this->data = other;
	}

	Interval(const Interval<T> &other) {
		this->p0 = other.p0
		this->pn = other.pn;
		this->npts = other.npts;
		this->step = other.step;
		this->data = other.data;
	}
	
	operator =(const Interval<T> &other)
	{
		this->p0 = other.p0
		this->pn = other.pn;
		this->npts = other.npts;
		this->step = other.step;
		this->data = other.data;

		return *this;
	}
	
	operator +=(T offset)
	{
		for (T &element : this->data) {
			element += offset;
		}

		return *this;
	}

	T operator [](int i) {
		return this->data[i];
	}
	
	Interval<T> slice(int start, int end) {
		std::vector<T> subinterval = std::vector<T>(this->data.begin() + start, this->data.begin() + end);
		return Interval<T>(subinterval);
	}
	
	Interval<T> slice_and_offset(int start, int end, T offset) {
		std::vector<T> subinterval = std::vector<T>(this->data.begin() + start, this->data.begin() + end);

		return Interval<T>(subinterval) += offset;
	}
	
	private:
	
	void initialize()
	{
		for (int i = 0; i < this->npts; i++) {
			this->data[i] = this->p0 + i * this->step;
		}
	}
}

class Interval<int>;
class Interval<float>;
class Interval<double>;