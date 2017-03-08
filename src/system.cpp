#include <stdio.h>
#include <math.h>
#include <string>
#include <map>
#include <vector>
#include <chrono>

using namespace std;

class System {
	public:
		System();
		vector<Molecule> molecules; // added 2-5-17
		Molecule proto;
		Constants constants;
		Pbc pbc;
        Stats stats;
        vector<vector<int>> atommap;

        // defines the "previous checkpoint time object
        std::chrono::time_point<std::chrono::system_clock> previous_checkpoint = std::chrono::system_clock::now();       

        // a function for de-bugging. Prints the current datetime and a string of text supplied in code.
        void checkpoint(string thetext) {
            if (constants.checkpoints_option == "on") {
            std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
            std::time_t thetime = std::chrono::system_clock::to_time_t(now);
            double time_elapsed = (std::chrono::duration_cast<std::chrono::microseconds>(now - previous_checkpoint).count()) /1000.0; // in ms

            printf("\nCheckpoint: %s ---> %.4f ms from last:  %s \n",std::ctime(&thetime), time_elapsed  , thetext.c_str());
            previous_checkpoint = std::chrono::system_clock::now();
            }
        }
};

System::System() {
}
