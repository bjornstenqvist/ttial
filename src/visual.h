#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>

class ProgressBar {
    private:
        unsigned int ticks_total;
        unsigned int tick = 0;
        bool show;
        double bar_width;
        const char done = '=';
        const std::string marker_end = ">";
        const char none = ' ';
        const std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

    public:
        ProgressBar(unsigned int total, bool show_progress=true, int width=70) : ticks_total(total), show(show_progress), bar_width(width-2) {}
        unsigned int operator++() { return tick++; }

        void display() {
            if(show) {
                double where = double(tick) / double(ticks_total);
                int pos = int(bar_width * where);
                std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
                int time = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time).count();

                if(pos > bar_width - 1 ) pos = bar_width - 1;
                std::cout << "[";
                for(int n = 0; n < pos; n++)
                    std::cout << done;
                std::cout << marker_end; //">";
                for(int n = pos+1; n < bar_width; n++)
                    std::cout << none;

                std::cout << "] " << int(where * 100.0) << "% "
                    << double(time) / 1000.0 << "s\r";
                std::cout.flush();
            }
        }

        void end_of_process() {
            if(show)
                std::cout << std::endl;
        }
};