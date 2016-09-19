#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <unistd.h>
#include <sys/mman.h>

#include <errno.h>


void swap_bytes_8( char * x )
{
  char c[4];
  
  c[0] = x[0];
  c[1] = x[1];
  c[2] = x[2];
  c[3] = x[3];

  x[0] = x[7];
  x[1] = x[6];
  x[2] = x[5];
  x[3] = x[4];

  x[4] = c[3];
  x[5] = c[2];
  x[6] = c[1];
  x[7] = c[0];
}

void mmap_open_(int *fd, char *name_in, char *flag)
{
  char name[120];
  int i;

  printf("called mmap_open\n");
 
  i = 0;
  while( i<120 && name_in[i] != '$' ) { name[i] = name_in[i]; i++; }

  if ( i >= 120 ) { perror("mmap_open: name too long\n"); exit(1); }
  name[i] = '\0';

  printf("trying to open: %s %c\n", name, *flag);

  if ( *flag == 'r' )
    *fd = open(name, O_RDONLY);
  else if( *flag == 'w' )
    *fd = open(name, O_RDWR);
  else {
    perror("mmap_open: unknown falg\n");
    exit(1);
  }

  printf("mmap_open ok, %d\n", *fd);

}


void mmap_close_( int *fd )
{
  close( *fd );
}


void mmap_read_(double *buf, int *fd, int *offset, 
	       int *idim, int *jdim, int *j0, int *j1)
{
  int size;
  double * db_buf;
  int i, j;

  size = (*idim)*(*jdim) * sizeof(double);

  db_buf = mmap(NULL, size, PROT_READ, MAP_PRIVATE, *fd, *offset);

  for(j=(*j0)-1; j<*j1; j++)
    for(i=0;i<*idim;i++)
      buf[i+(*idim)*j] = db_buf[i+(*idim)*j];

  munmap(db_buf, size);

}


void mmap_write_(double *buf, int *fd, int *offset, 
	       int *idim, int *jdim, int *j0, int *j1)
{
  int size;
  char * cbuf;
  double * db_buf;
  int i, j;

  printf("called mmap_write with %d %d %d %d %d %d\n",
	 *fd, *offset, *idim, *jdim, *j0, *j1);

  size = (*idim)*(*jdim) * sizeof(double) + 8 ;

  printf("mmap %d %d %d %d %d %d\n",
	 NULL, size, PROT_WRITE, MAP_SHARED, *fd, *offset );

  //getchar();
  cbuf = mmap(NULL, size, PROT_WRITE|PROT_READ, MAP_SHARED, *fd, *offset);
  if ( cbuf == MAP_FAILED ) {
    perror("mmap_write: mmap error\n"); exit(1); }

  msync(cbuf, size, MS_INVALIDATE);

  db_buf = cbuf + 4;
  printf("mmap ok\n");

  for(j=(*j0)-1; j<(*j1); j++)
    for(i=0;i<(*idim);i++) {
      double x = buf[i+(*idim)*j];
      double y = db_buf[i+(*idim)*j];
      swap_bytes_8(&y);
      printf("loop: %d %d %g %g\n", i, j, x, y);
      swap_bytes_8(&x);
      db_buf[i+(*idim)*j] = x;
    }

  printf("%f %f %f %f\n", db_buf[0], db_buf[1], db_buf[2], db_buf[3] );

  //getchar();

  msync(db_buf, size, MS_SYNC);
  munmap(db_buf, size);

}


//#define TEST_MMAP_UTILS
#ifdef TEST_MMAP_UTILS

int main()
{
  int fd;
  char * name="foobar";

  double a[16];
  int offset, idim, jdim, j0, j1;

  a[0] = 1;
  a[1] = 2;
  a[2] = 3;
  a[3] = 4;

  a[4] = 11;
  a[5] = 12;
  a[6] = 13;
  a[7] = 14;

  offset = 0;
  idim = 4;
  jdim = 4;
  j0 = 2;
  j1 = 2;

  fd = open(name, O_RDWR);
  if ( fd < 0 ) { perror("main: can't open\n"); exit(1); }

  mmap_write( a, &fd, &offset, &idim, &jdim, &j0, &j1);

  close( fd );
  

  return 0;
}

#endif
