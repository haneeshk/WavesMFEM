#include "save.hpp"

void saveData(std::string folder, const int cycle, const Vector& u, const GridFunction &nodes2){
	
	std::ostringstream oss;
oss << std::setw(5) << 						std::setfill('0') << cycle;
	std::string padded = oss.str();  
	std::ofstream u_out(folder+"/u_data_t" + padded + ".dat");
		
		  for (int i = 0; i < u.Size(); i++)
		  {
			  double x1 = nodes2(i);
			  double val = u(i);
			  u_out << x1 << " " << val << "\n";
		  }
		  u_out.close();
	
	return;

}
