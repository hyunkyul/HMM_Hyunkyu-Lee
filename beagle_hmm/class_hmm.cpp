/*
 * =====================================================================================
 *
 *       Filename:  class_hmm.cpp
 *
 *    Description:  duplicate beagle hmm using C++ with classes
 *
 *        Version:  1.0
 *        Created:  03/23/2015 12:46:48
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Hyunkyu Lee, Buhm Han*
 *   Organization:  Han lab
 *
 * =====================================================================================
 */


// libraries
#include <iostream>

class level	
{
private:
	int transition;
public:
	int set_trans(int);
	int trans() {return transition;}
} lvl1,lvl2,lvl3,lvl4;

int level::set_trans(int x)
{
	transition=x;
	return 0;
}

int main() 
{
	int i;
		
	std::cout << "Initializing : complete \n";
	lvl1.set_trans(2);
	lvl2.set_trans(3);
	lvl3.set_trans(3);
	lvl4.set_trans(5);

	std::cout << "transition map	: " << lvl1.trans() << lvl2. trans() << lvl3.trans() << lvl4.trans() <<"\n\n";
}


