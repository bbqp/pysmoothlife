#include <vector>
#include <concept>
#include <functional>

template class MeshArray <std::floating_point T>
{
	private:
	
	Meshgrid<T> domain;
	std::vector<T> state;
	
	public:
	
	MeshArray(const Meshgrid &mgrid) {
		this->domain = mgrid;
		this->state = 
	}
	
	void apply(std::function<T, T> f, std::vector<T> &out)
	{
		return return_state;
	}
	
	~MeshArray();
}