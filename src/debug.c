#include <stdio.h>
#include <unistd.h>

void waitForAttach(int myid) {
  volatile int i = 0;
  char host[256];

  if (myid == 0) {
    printf("PID %d is ready for attach\n", getpid());
    fflush(stdout);

    // Set i to something else using GDB.
    while (i == 0) {
      sleep(5);
    }
  }
}