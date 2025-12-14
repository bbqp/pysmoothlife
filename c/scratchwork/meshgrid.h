#include <vector>
#include <concept>

template struct Meshgrid <T>
{
	T x0, xn, dx;
	T y0, yn, dx;
	int numx, numy;
	std::vector<T> XX, YY;
	
	Meshgrid<T>()
	{
		this->x0 = 0;
		this->xn = 1;
		this->numx = 11;
		this->dx = (this->xn - this->x0) / (this->numx - 1);
		
		this->y0 = 0;
		this->yn = 1;
		this->numy = 11;
		this->dy = (this->yn - this->y0) / (this->numy - 1);
		
		this->initialize();
	}
	
	Meshgrid<T>(T x0, T xn, int numx, T y0, T yn, int numy) {
		this->x0 = x0;
		this->xn = xn;
		this->numx = numx;
		this->dx = (this->xn - this->x0) / (this->numx - 1);
		
		this->y0 = y0;
		this->yn = yn;
		this->numy = numy;
		this->dy = (this->yn - this->y0) / (this->numy - 1);
		
		this->initialize();
	}
	
	Meshgrid<T>(const std::vector<T> &x, const std::vector<T> &y) {
		this->x0 = x.front();
		this->xn = x.back();
		this->numx = x.size();
		this->dx = (this->xn - this->x0) / (this->numx - 1);
		
		this->y0 = this->y.p0;
		this->yn = this->y.pn;
		this->numy = this->y.npts;
		this->dy = (this->yn - this->y0) / (this->numy - 1);
		
		this->initialize(x, y);
	}
	
	Meshgrid<T>(const Interval<T> &x, const Interval<T> &y) {
		this->x0 = this->x.p0;
		this->xn = this->x.pn;
		this->numx = this->x.npts
		this->dx = this->x.step;
		
		this->y0 = this->y.p0;
		this->yn = this->y.pn;
		this->numy = this->y.npts;
		this->dy = this->y.step;
		
		this->initialize(x, y);
	}

	Meshgrid<T>(const Meshgrid<T> &other) {
		this->x0 = other.x0;
		this->xn = other.xn;
		this->numx = other.numx;
		this->dx = other.dx;
		
		this->y0 = other.y0;
		this->yn = other.yn;
		this->numy = other.numy;
		this->dy = other.dy;
		
		this->XX = other.XX;
		this->YY = other.YY;
	}
	
	Meshgrid<T> operator =(const Meshgrid<T> &other)
	{
		if (this != &other) {
			this->x0 = other.x0;
			this->xn = other.xn;
			this->numx = other.numx;
			this->dx = other.dx;
			
			this->y0 = other.y0;
			this->yn = other.yn;
			this->numy = other.numy;
			this->dy = other.dy;
			
			this->XX = other.XX;
			this->YY = other.YY;
		}

		return *this;
	}
	
	Meshgrid<T> & offset(T offsetx, T offsety)
	{
		if (offsetx != 0) {
			for(T &element : this->XX) {
				element += offset;
			}
		}
		
		if(offsety != 0) {
			for(T &element : this->YY) {
				element += offset;
			}
		}

		return *this;
	}
	
	Meshgrid<T> & offset_modulus(T offsetx, T offsety)
	{
		T numx = this->numx;
		T numy = this->numy;

		if (offsetx != 0) {
			for(T &element : this->XX) {
				element = (element + numx + offset) % numx;
			}
		}
		
		if(offsety != 0) {
			for(T &element : this->YY) {
				element = (element + numx + offset) % numy;
			}
		}

		return *this;
	}
	
	Meshgrid<T> slice_copy(int istart, int iend, int jstart, int jend) {
		T x0 = this->x0 + istart * this->dx;
		T xn = this->x0 + iend * this->dx;
		int numx = iend - istart + 1;
		
		T y0 = this->y0 + jstart * this->dy;
		T yn = this->y0 + jend * this->dy;
		int numy = jend - jstart + 1;
		
		return Meshgrid<T>(x0, xn, numx, y0, yn, numy);
	}

	Meshgrid<T> slice_copy_and_offset(int istart, int iend, int jstart, int jend, T offsetx, T offsety) {
		return this->slice_copy(istart, iend, jstart, jend)->offset(offsetx, offsety);
	}
	
	private:
	
	void initialize()
	{
		int idx;
		for (int j = 0; j < this->numy; j++) {
			idx = j * this->numx;
			yj = this->y0 + j * this->dy

			for (int i = 0; i < this->npts; i++) {
				xi = this->x0 + i * this->dx;

				this->XX[idx + i] = xi;
				this->YY[idx + i] = yj;
			}
		}
	}
	
	void initialize(const Interval<T> &x, const Interval<T> &y)
	{
		int idx;
		
		for (int j = 0; j < this->numy; j++) {
			idx = j * this->numx;

			for (int i = 0; i < this->npts; i++) {
				this->XX[idx + i] = x[i];
				this->YY[idx + i] = y[j];
			}
		}
	}
	
	void initialize(const std::vector<T> &x, const std::vector<T> &y)
	{
		int idx;
		
		for (int j = 0; j < this->numy; j++) {
			idx = j * this->numx;

			for (int i = 0; i < this->npts; i++) {
				this->XX[idx + i] = x[i];
				this->YY[idx + i] = y[j];
			}
		}
	}
	
};

class Meshgrid<int>;
class Meshgrid<float>;
class Meshgrid<double>;