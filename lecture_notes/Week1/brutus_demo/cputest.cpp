/* Simple cpu code that prints the hostname of the current node */
#include <stdio.h>
#include <unistd.h>
#include "timer.hpp"

int main(int  argc, char *argv[])
{
	char name[256];
	
	timer t;

	t.start();
	gethostname(name, 256);
	sleep(1);
	t.stop();

	printf("Hello, I am running on host: %s\n", name);
	printf("Elapsed time = %lf seconds\n", t.get_timing());
	
	return 0;
}
