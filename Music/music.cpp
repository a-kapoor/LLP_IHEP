#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib> 
#include <ctime> 

using namespace std;

extern"C" {
    void initialize_music_();
    void muon_transport_(double *x0,double *y0, double *z0, double *theta, double *phi, double *energy, double *depth, double *dR, double *time);
}

template <class T>
T GetRandAB (T a, T b) {
    T result;
    result = (a + static_cast <float> (rand()) / ( static_cast <float> (RAND_MAX/(b-a))));
    return (result);
}

int main()
{    
    
    //Important to initialize
    //std::srand(static_cast<unsigned int>(std::time(nullptr))); 
    srand( (unsigned)time( NULL ) );
    initialize_music_();
        
    //cout<<"Initial four vector"<<endl;
    
    //std::cout << std::left << std::setw(20) << "px" << std::setw(20) << "py"<< std::setw(20) << "pz"<< std::setw(20) << "E" << '\n';
    //std::cout << std::left << std::setw(20) << "Richard" << std::setw(20) << 1.0 << '\n';
    cout<<"---------------------------------------------------------------------------------------------"<<endl;
    cout<<"------------------- O U T P U T        R E C O R D-----------------------------------------------"<<endl;
    std::cout << std::left << std::setw(20) << "Event" << std::setw(20) << "Muon" << std::setw(20) << "px" << std::setw(20) << "py"<< std::setw(20) << "pz"<< std::setw(20) << "E" << '\n';
    for(int ll=0;ll<20;ll++){
	//cout<<"|--------------- E V E N T"<<endl;
	for(int lll=0;lll<GetRandAB<int>(1,4);lll++){
	    double x0=GetRandAB<float>(100,200);
	    double y0=GetRandAB<float>(100,200);
	    double z0=GetRandAB<float>(100,1000);
	    double theta=GetRandAB<float>(0,3);
	    double phi=GetRandAB<float>(0,3);
	    double energy=GetRandAB<float>(5,200);
	    double depth=GetRandAB<float>(100,300);
	    double dR=2.6;
	    double time=0;
	    
	    float massmu=0.10566;
	    float p=std::sqrt((energy*energy) - (massmu*massmu));
	    float px=p*std::sin(theta)*std::cos(phi);
	    float py=p*std::sin(theta)*std::sin(phi);
	    float pz=p*std::cos(theta);
	    
	    //cout<<"----------------------- I N P U T      M U O N----------------------------------------------------------"<<endl;
	    //cout<<"Final four vector"<<endl;
	    
	    //std::cout << std::left << std::setw(20) << "px" << std::setw(20) << "py"<< std::setw(20) << "pz"<< std::setw(20) << "E" <<std::setw(20) << "t" << '\n';
	    //std::cout << std::left << std::setw(20) << px << std::setw(20) << py<< std::setw(20) << pz<< std::setw(20) << energy <<std::setw(20) << time << '\n';
	    
	    //cout<<"----------------------------------------------------------------------------------------------------------"<<endl;
	    
	    //cout<<"energy ini = "<<energy<<endl;
	    //cout<<"x0 ini = "<<x0<<endl;
	    //cout<<"time ini = "<<time<<endl;
	    //cout<<"dR ini = "<<dR<<endl;
	    
	    //cout<<"Initial four vector"<<endl;
	    
	    //std::cout << std::left << std::setw(20) << "px" << std::setw(20) << "py"<< std::setw(20) << "pz"<< std::setw(20) << "E" << '\n';
	    //std::cout << std::left << std::setw(20) << px << std::setw(20) << py<< std::setw(20) << pz<< std::setw(20) << E << '\n';
	    
	    /////Call Fortran from C++
	    muon_transport_(&x0,&y0, &z0, &theta, &phi, &energy, &depth, &dR, &time);
	    
	    //cout<<"energy final = "<<energy<<endl;
	    //cout<<"x0 final = "<<x0<<endl;
	    //cout<<"dR final = "<<dR<<endl;
	    
	    
	    p=std::sqrt((energy*energy) - (massmu*massmu));
	    px=p*std::sin(theta)*std::cos(phi);
	    py=p*std::sin(theta)*std::sin(phi);
	    pz=p*std::cos(theta);
	    
	    //cout<<"Final four vector"<<endl;
	    
	    std::cout << std::left << std::setw(20) << ll << std::setw(20) << lll << std::setw(20) << px << std::setw(20) << py<< std::setw(20) << pz<< std::setw(20) << energy << '\n';
	}
	cout<<"----------------------------------------------------------------------------------------------------------"<<endl;
    }

    return 0;
}
