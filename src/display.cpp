#ifndef  __DISPLAY_HPP__
#define  __DISPLAY_HPP__

#include "display.hpp"

/**
 * @brief displays the progress of a process on the console
 * @param percent  the percent of the process that is completed
 *
 * The function takes a single argument which represents the percentage of the process that is completed,
 * and then displays a simple progress bar on the console that indicates the percentage of completion.
 * It uses '=' to represent the completed portion of the process and blanks to represent the remaining portion.
 * The function uses std::cout to output the progress bar and flush the buffer after each update.
 */

void displayProgress(unsigned long percent) {             // Display the progression
    std::cout << " REAL progression        [";
    for (int i = 0; i < 49; i++) {
        if (i < (int) percent / 2) {
            std::cout << "=";
        } else {
            std::cout << " ";
        }
    }
    std::cout << "] " <<  (int) percent +1 << "%\r";
    std::cout.flush();
}

#endif